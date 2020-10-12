/*
 * By Junkang Huang, Dec. 2020
 * 
 * performing parallel DTRSM operation
 * based on:
 *      Sec 3.1 of Matrix Computations  ---by H.Golub & F.Van Loan
 * 
 */

# include "hpl.h"

/* pack A[0:bm, 0:bn] into buffer, column major */
void packblock
(
    double *                        buffer,
    const double *                  A,
    const int                       lda,
    const int                       bm,
    const int                       bn
)
{
    int i, j;
    for (j = 0; j < bn; ++j)
    {
         for (i = 0; i < bm; ++i) 
        {
            *Mptr(buffer, i, j, bm) = *Mptr(A, i, j, lda);
        }
    }
}

/*
 * cal_pre():
 * 
 *      the LU factors of A is stored firmly in fatcors. cal_pre() calculate the inverse of L
 *    and U parallelly based on block DTRSM, then stores them separately in preL and preU.
 */
void cal_pre
(
   HPL_T_grid *                     GRID,           /* (input) process grid information */
   const double *                   factors,        /* (input) LU factors of matrix */
   double *                         preL,           /* (output) preconditioning matrices L-1 */
   double *                         preU,           /* (output) preconditioning matrices U-1 */
   const int                        mp,             /* (input) local number of rows */
   const int                        nq,             /* (input) local number of columns */
   const int                        rmp,            /* (input) local number of rows(expanded) */
   const int                        n,              /* (input) global size of matrix */
   const int                        nb              /* (input) blocking parameter */
)
{
/* local variables */
    int q, i, j, k, ilr, ilc, ipr, ipc, irank, jlr, jlc, jpr, jpc, jrank, klr, klc, kpr, kpc, krank;
    int myrow, mycol, nprow, npcol;
    double *buffer, *buffer2;

    /* some grid parameters */
    myrow = GRID->myrow;
    mycol = GRID->mycol;
    nprow = GRID->nprow;
    npcol = GRID->npcol;

    /* message passing buffer of a block */
    buffer  = (double *)malloc(nb*nb*sizeof(double));
    buffer2 = (double *)malloc(nb*nb*sizeof(double));

    /* number of blocks (including incomplete blocks) along row or column */
    q = (n % nb) ? n / nb + 1 : n / nb;

/* calculate L-1 */

    /* load preL with unit matrix In */
    for (j = 0; j < n; ++j)
    {
        HPL_indxg2lp(&jlr, &jpr, j, nb, nb, 0, nprow);
        HPL_indxg2lp(&jlc, &jpc, j, nb, nb, 0, npcol);
        if (myrow == jpr && mycol == jpc)
        {
            *Mptr(preL, jlr, jlc, rmp) = 1.0;
        }
    }
    for (j = 0; j < q; ++j) 
    {
        /* calculate the local index and process index of global diagonal block index [j,j] */
        HPL_indxg2lp(&jlr, &jpr, j*nb, nb, nb, 0, nprow);
        HPL_indxg2lp(&jlc, &jpc, j*nb, nb, nb, 0, npcol);

        /* pack block data into buffer */
        if (myrow == jpr && mycol == jpc)
        {
            packblock(buffer, Mptr(factors, jlr, jlc, rmp), rmp, nb, nb);
        }

        /* solve for preL[j, 0:j+1] */
        for (i = 0; i < j + 1; ++i)
        {
            /* calculate local index and process index of global block index [j,i] */
            HPL_indxg2lp(&ilr, &ipr, j*nb, nb, nb, 0, nprow);
            HPL_indxg2lp(&ilc, &ipc, i*nb, nb, nb, 0, npcol);

            /* send block factors[j,j] to process which contains preL[j,i] */
            if (myrow == jpr && mycol == jpc) 
            {
                if (ipc != jpc)
                {
                    HPL_dsend(buffer, nb*nb, ipc, 0, GRID->row_comm);
                }
            }
            if (myrow == ipr && mycol == ipc)
            {
                if (ipc != jpc)
                {
                    HPL_drecv(buffer, nb*nb, jpc, 0, GRID->row_comm);
                }
                /* perform serial DTRSM */
                HPL_dtrsm(HplColumnMajor, HplLeft, HplLower, HplNoTrans, HplUnit, nb, 
                    nb, 1, buffer, nb, Mptr(preL, ilr, ilc, rmp), rmp);
            }
        }

        /* update preL[j+1:, :], preL[j+1:, j+1:] should be filled with zeros in where 
         * storing U elements in fact, so just update preL[j+1:,:j+1] will be enough. 
         */
        for (k = j + 1; k < q; ++k)
        {
            /* calculate the local index and process index of global block index [k,j] */
            HPL_indxg2lp(&jlr, &jpr, k*nb, nb, nb, 0, nprow);
            HPL_indxg2lp(&jlc, &jpc, j*nb, nb, nb, 0, npcol);

            /* calculate rank from grid indices */
            if (GRID->order == HPL_ROW_MAJOR)
            {
                jrank = jpr * npcol + jpc;
            }else
            {
                jrank = jpc * nprow + jpr;
            } 

            /* pack block data into buffer */
            if (GRID->iam == jrank)
            {
                packblock(buffer, Mptr(factors, jlr, jlc, rmp), rmp, nb, nb);
            }           

            for (i = 0; i < j + 1; ++i)
            {
                /* calculate local index and process index of global block index [k,i] */
                HPL_indxg2lp(&ilr, &ipr, k*nb, nb, nb, 0, nprow);
                HPL_indxg2lp(&ilc, &ipc, i*nb, nb, nb, 0, npcol);

                /* calculate rank from grid indices */
                if (GRID->order == HPL_ROW_MAJOR)
                {
                    irank = ipr * npcol + ipc;
                }else
                {
                    irank = ipc * nprow + ipr;
                } 

                /* calculate the local index and process index of global block index [j,i] */
                HPL_indxg2lp(&klr, &kpr, j*nb, nb, nb, 0, nprow);
                HPL_indxg2lp(&klc, &kpc, i*nb, nb, nb, 0, npcol);

                /* calculate rank from grid indices */
                if (GRID->order == HPL_ROW_MAJOR)
                {
                    krank = kpr * npcol + kpc;
                }else
                {
                    krank = kpc * nprow + kpr;
                }

                if (GRID->iam == jrank && jrank != irank)
                {/* send factors[k,j] to the process containing preL[k,i] through rank */
                    HPL_dsend(buffer, nb*nb, irank, 1, GRID->all_comm);
                }
                if (GRID->iam == krank)
                {/* send preL[j,i] to the process containing preL[k,i] through rank */
                    packblock(buffer2, Mptr(preL, klr, klc, rmp), rmp, nb, nb);
                    if (GRID->iam != irank)
                        HPL_dsend(buffer2, nb*nb, irank, 2, GRID->all_comm);
                }
                if (GRID->iam == irank)
                {
                    /* collect factors[k,j] and preL[j,i] */
                    if (irank != jrank)
                        HPL_drecv(buffer, nb*nb, jrank, 1, GRID->all_comm);
                    if (irank != krank)
                        HPL_drecv(buffer2, nb*nb, krank, 2, GRID->all_comm);
                    /* perform serial DGEMM to update */
                    HPL_dgemm(HplColumnMajor, HplNoTrans, HplNoTrans, nb, nb, nb, -1, buffer,
                         nb, buffer2, nb, 1, Mptr(preL, ilr, ilc, rmp), rmp);
                }
            }
        }
    }

/* calculate U-1 */

    /* load preU with unit matrix In */
    for (j = 0; j < n; ++j)
    {
        HPL_indxg2lp(&jlr, &jpr, j, nb, nb, 0, nprow);
        HPL_indxg2lp(&jlc, &jpc, j, nb, nb, 0, npcol);
        if (myrow == jpr && mycol == jpc)
        {
            *Mptr(preU, jlr, jlc, rmp) = 1.0;
        }
    }

    for (j = 0; j < q; ++j) 
    {
        /* calculate the local index and process index of global block index [j,j] */
        HPL_indxg2lp(&jlr, &jpr, j*nb, nb, nb, 0, nprow);
        HPL_indxg2lp(&jlc, &jpc, j*nb, nb, nb, 0, npcol);

        /* pack block data into buffer */
        if (myrow == jpr && mycol == jpc)
        {
            packblock(buffer, Mptr(factors, jlr, jlc, rmp), rmp, nb, nb);
        }

        /* solve for preU[0:j+1, j] */
        for (i = 0; i < j + 1; ++i)
        {
            /* calculate local index and process index of global block index [i,j] */
            HPL_indxg2lp(&ilr, &ipr, i*nb, nb, nb, 0, nprow);
            HPL_indxg2lp(&ilc, &ipc, j*nb, nb, nb, 0, npcol);

            /* send block factors[j,j] to process which contains preU[i,j] */
            if (myrow == jpr && mycol == jpc) 
            {
                if (jpr != ipr)
                {
                    HPL_dsend(buffer, nb*nb, ipr, 0, GRID->col_comm);
                }
            }
            if (myrow == ipr && mycol == ipc)
            {
                if (jpr != ipr)
                {
                    HPL_drecv(buffer, nb*nb, jpr, 0, GRID->col_comm);
                }
                /* perform serial DTRSM */
                HPL_dtrsm(HplColumnMajor, HplRight, HplUpper, HplNoTrans, HplNonUnit, nb, 
                    nb, 1, buffer, nb, Mptr(preU, ilr, ilc, rmp), rmp);
            }
        }

        /* update preU[:, j+1:], preU[j+1:, j+1:] should be filled with zeros in where 
         * storing L elements in fact, so just update preU[:j+1,j+1:] will be enough. 
         */
        for (k = j + 1; k < q; ++k)
        {
            /* calculate the local index and process index of global block index [j,k] */
            HPL_indxg2lp(&jlr, &jpr, j*nb, nb, nb, 0, nprow);
            HPL_indxg2lp(&jlc, &jpc, k*nb, nb, nb, 0, npcol);
            /* calculate rank from grid indices */
            if (GRID->order == HPL_ROW_MAJOR)
            {
                jrank = jpr * npcol + jpc;
            }else
            {
                jrank = jpc * nprow + jpr;
            } 

            /* pack block data into buffer */
            if (GRID->iam == jrank)
            {
                packblock(buffer, Mptr(factors, jlr, jlc, rmp), rmp, nb, nb);
            }       

            for (i = 0; i < j + 1; ++i)
            {
                /* calculate local index and process index of global block index [i,k] */
                HPL_indxg2lp(&ilr, &ipr, i*nb, nb, nb, 0, nprow);
                HPL_indxg2lp(&ilc, &ipc, k*nb, nb, nb, 0, npcol);

                /* calculate rank from grid indices */
                if (GRID->order == HPL_ROW_MAJOR)
                {
                    irank = ipr * npcol + ipc;
                }else
                {
                    irank = ipc * nprow + ipr;
                } 

                /* calculate the local index and process index of global block index [i,j] */
                HPL_indxg2lp(&klr, &kpr, i*nb, nb, nb, 0, nprow);
                HPL_indxg2lp(&klc, &kpc, j*nb, nb, nb, 0, npcol);

                /* calculate rank from grid indices */
                if (GRID->order == HPL_ROW_MAJOR)
                {
                    krank = kpr * npcol + kpc;
                }else
                {
                    krank = kpc * nprow + kpr;
                }

                if (GRID->iam == jrank && jrank != irank)
                {/* send factors[j,k] to the process containing preU[i,k] through rank */
                    HPL_dsend(buffer, nb*nb, irank, 1, GRID->all_comm);
                }

                if (GRID->iam == krank)
                {/* send preU[i,j] to the process containing preU[i,k] through rank */
                    packblock(buffer2, Mptr(preU, klr, klc, rmp), rmp, nb, nb);
                    if (GRID->iam != irank)
                        HPL_dsend(buffer2, nb*nb, irank, 2, GRID->all_comm);
                }

                if (GRID->iam == irank)
                {
                    /* collect factors[j,k] and preU[i,j] */
                    if (irank != jrank)
                        HPL_drecv(buffer, nb*nb, jrank, 1, GRID->all_comm);
                    if (irank != krank)
                        HPL_drecv(buffer2, nb*nb, krank, 2, GRID->all_comm);
                    /* perform serial DGEMM to update */
                    HPL_dgemm(HplColumnMajor, HplNoTrans, HplNoTrans, nb, nb, nb, -1, buffer2,
                         nb, buffer, nb, 1, Mptr(preU, ilr, ilc, rmp), rmp);
                }
            }
        }
    }

    if (buffer)  free(buffer);
    if (buffer2) free(buffer2);

/* End of cal_pre() */
}