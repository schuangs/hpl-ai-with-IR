/*
 * By Junkang Huang, Dec. 2020
 * 
 * performing parallel DTRSM operation
 * based on:
 *      Sec 3.1 of Matrix Computations  ---by H.Golub & F.Van Loan
 * 
 */

# include "hpl.h"

void HPL_my_pdtrsm
(
   HPL_T_grid *                     GRID,           /* (input) process grid information */
   const double *                   factors,        /* (input) LU factors of matrix */
   double *                         pre,            /* (output) preconditioning matrices, U-1 and L-1 */
   const int                        mp,             /* (input) local number of rows */
   const int                        nq,             /* (input) local number of columns */
   const int                        n,              /* (input) global size of matrix */
   const int                        nb              /* blocking parameter */
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

    /* message passing buffer */
    buffer  = (double *)malloc(nb*nb*sizeof(double));
    buffer2 = (double *)malloc(nb*nb*sizeof(double));

    /* number of blocks along row or column */
    q = n / nb + 1;

    /* make pre a unit matrix */
    memset(pre, 0, mp*nq*sizeof(double));
    for (j = 0; j < n; ++j)
    {
        HPL_indxg2lp(&jlr, &jpr, j, nb, nb, 0, nprow);
        HPL_indxg2lp(&jlc, &jpc, j, nb, nb, 0, npcol);
        if (myrow == jpr && mycol == jpc)
        {
            *Mptr(pre, jlr, jlc, mp) = 1.0;
        }
    }

/* calculate L-1 */
    for (j = 0; j < q; ++j) 
    {
        /* calculate the local index and process index of global block index [j,j] */
        HPL_indxg2lp(&jlr, &jpr, j*nb, nb, nb, 0, nprow);
        HPL_indxg2lp(&jlc, &jpc, j*nb, nb, nb, 0, npcol);

        /* pack block data into buffer */
        if (myrow == jpr && mycol == jpc)
        {
            packblock(buffer, Mptr(factors, jlr, jlc, mp), nb);
        }

        /* solve for pre[j, 0:j] */
        for (i = 0; i < j; ++i)
        {
            /* calculate local index and process index of global block index [j,i] */
            HPL_indxg2lp(&ilr, &ipr, j*nb, nb, 0, nprow);
            HPL_indxg2lp(&ilc, &ipc, i*nb, nb, 0, npcol);

            /* send block factors[j,j] to process which contains pre[j,i] */
            if (myrow == jpr && mycol == jpc) 
            {
                HPL_send(buffer, nb*nb, ipc, 0, GRID->row_comm);
            }
            else if (myrow == ipr && mycol == ipc)
            {
                HPL_recv(buffer, nb*nb, jpc, 0, GRID->row_comm);
                /* perform serial DTRSM */
                HPL_dtrsm(HplColumnMajor, HplLeft, HplLower, HplNoTrans, HplUnit, nb, 
                    nb, 1, buffer, nb, Mptr(pre, ilr, ilc, mp), mp);
            }
        }

        /* update pre[j+1:, :], pre[j:, j+1:] should be filled with zeros in where 
         * storing U elements in fact, so just update pre[j+1:,:j] will be enough. 
         */
        for (k = j + 1; k < q; ++k)
        {
            /* calculate the local index and process index of global block index [k,j] */
            HPL_indxg2lp(&jlr, &jpr, j*nb, nb, nb, 0, nprow);
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
                packblock(buffer, Mptr(factors, jlr, jlc, mp), nb);
            }           

            for (i = 0; i < j; ++i)
            {
                /* calculate local index and process index of global block index [k,i] */
                HPL_indxg2lp(&ilr, &ipr, k*nb, nb, 0, nprow);
                HPL_indxg2lp(&ilc, &ipc, i*nb, nb, 0, npcol);

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
                {/* send factors[k,j] to the process containing pre[j,i] through rank */
                    HPL_send(buffer, nb*nb, irank, 1, GRID->all_comm);
                }
                else if (GRID->iam == krank)
                {/* send pre[j,i] to the process containing pre[j,i] through rank */
                    packblock(buffer, Mptr(pre, klr, klc, mp), nb);
                    if (GRID->iam != irank)
                        HPL_send(buffer, nb*nb, irank, 2, GRID->all_comm);
                }
                else if (GRID->iam == irank)
                {
                    /* collect factors[k,j] and pre[j,i] */
                    if (irank != jrank)
                        HPL_recv(buffer, nb*nb, jrank, 1, GRID->all_comm);
                    if (irank != krank)
                        HPL_recv(buffer2, nb*nb, krank, 2, GRID->all_comm);
                    /* perform serial DGEMM to update */
                    HPL_dgemm(HplColumnMajor, HplNoTrans, HplNoTrans, nb, nb, nb, -1, buffer,
                         nb, buffer2, nb, 1, Mptr(pre, ilr, ilc, mp), mp);
                }
            }
        }
    }

    if (buffer)  free(buffer);
    if (buffer2) free(buffer2);

    /* End of HPL_my_pdtrsm() */
}