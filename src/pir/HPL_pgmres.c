/*
 * By Junkang Huang, August, 2020
 * 
 */

/*
 * Include files
 */ 
#include "hpl.h"


#ifdef STDC_HEADERS
void HPL_pgmres
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   HPL_T_pdmat *                    A,
   const double *                   FACT,
   double *                         R
)
#else
void HPL_pdgesv0
( GRID, ALGO, A )
   HPL_T_grid *                     GRID;
   HPL_T_palg *                     ALGO;
   HPL_T_pdmat *                    A;
   const double *                   FACT;
   double *                         R;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pir performs iterative refinement procesure to enhance the accur-
 * acy of solved linear system by LU factorization.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * ALGO    (global input)                HPL_T_palg *
 *         On entry,  ALGO  points to  the data structure containing the
 *         algorithmic parameters to be used for this test.
 * 
 * A       (local input)                 HPL_T_pdmat *
 *         On entry, A points to the data structure containing the local
 *         array  information.  The  residual  vector is  stored  in the 
 *         column between [ A | b ] and X as input. The corection vector
 *         will overwrite the residual vector as output.
 * 
 * FACT    (local input)                 const double *
 *         On entry, FA  points  to  the  data structure  containing the 
 *         L and U factors which are used as preconditioning matrices.
 * 
 * R       (local input/output)          double *
 *         On  entry, R  points to  the data  structure  containing  the 
 *         residual vector  which plays a role as the right-hand side of 
 *         this  routine to obtain  correcting vector. Correcting vector 
 *         will overwrite R as output.
 *
 * ---------------------------------------------------------------------
 */ 

/*
 * .. Local Variables ..
 */
   int                i,  n, nloc, m,  lwork, irc[5], icntl[8], info[3], 
                      inv_info[3];
   int                mp, nq, lda;
   int                ione, izero;
   int                diagi, diagj, local_i, local_j;
   double           * work, * rptr, cntl[5], rinfo[2], one, zero;
   double           * invptrL, * invptrU;
   void             * vL, * vU;
   char               ctran;

/* ..
 * .. Executable Statements ..
 */
   n = A->n;  
   m = 8; one = 1.0; zero = 0.0;
   ione = 1; izero = 0;
   mp = A->mp; nq = A->nq-1; lda = A->ld;
   nloc = Mmax(mp, nq);

/*
 * Set integer control pararmeters
 */
   icntl[0] = 6;     /* stream number for error messages */
   icntl[1] = 6;     /* do not display warning messages */
   icntl[2] = 0;     /* do not display convergence history */
   icntl[3] = 1;     /* left preconditioning */
   icntl[4] = 3;     /* orthogonalization scheme: ICGS */
   icntl[5] = 1;     /* initial guess of solution vector is zero */
   icntl[6] = 2000;  /* maximum number of iterations */
   icntl[7] = 0;     /* strategy to compute residual at restart */

/*
 * Set float control pararmeters
 */
   cntl[0] = 1e-10;  /* convergence tolerance for backward error */
   cntl[1] = 0;      /* normalizing factor ALPHA */
   cntl[2] = 0;      /* normalizing factor BETA */
   cntl[3] = 0;      /* normalizing factor ALPHAP */
   cntl[4] = 0;      /* normalizing factor BETAP */

/*
 * set work space size: lwork
 */
   if (icntl[4] == 0 || icntl[4] == 1)
   {
      lwork = m*m + m*(nloc+5) + 5*nloc + 2;
   }
   else if (icntl[4] == 2 || icntl[4] == 3)
   {
      lwork = m*m + m*(nloc+5) + 5*nloc + m + 1;
   }
   if (icntl[7] != 1)
   {
      lwork += nloc;
   }
/*
 * allocate work space and load rhs into work space
 */
   work = (double *)malloc( lwork*sizeof(double) );
   memcpy(work+nloc, R, mp);
/*
 * load initial guess of solution
 */
   memcpy(work, A->X, nq);

/*
 * allocate space for precondition matrix
 */
   vL = (void*)malloc( ( (size_t)(ALGO->align) + 
                           (size_t)(A->ld+1) * (size_t)(nq+1) ) *
                         sizeof(double) );
   vU = (void*)malloc( ( (size_t)(ALGO->align) + 
                           (size_t)(A->ld+1) * (size_t)(nq+1) ) *
                         sizeof(double) );
   inv_info[0] = (vL == NULL || vU == NULL); 
   inv_info[1] = GRID->myrow; 
   inv_info[2] = GRID->mycol;
   (void) HPL_all_reduce( (void *)(inv_info), 3, HPL_INT, HPL_max,
                          GRID->all_comm );
   
   if( inv_info[0] != 0 )
   {
      /* some processes might have succeeded with allocation */
      if (vL) free(vL);
      if (vU) free(vU);
      return;
   }

   invptrL = (double *)HPL_PTR( vL,
                               ((size_t)(ALGO->align) * sizeof(double) ) );
   invptrU = (double *)HPL_PTR( vU,
                               ((size_t)(ALGO->align) * sizeof(double) ) );

/*
 * fill the space with unit matrix for calculating preconditioning matrix
 */
   memset(invptrL, 0, ( (size_t)(ALGO->align) + 
                           (size_t)(A->ld+1) * (size_t)(nq+1) ) *
                         sizeof(double));
   memset(invptrU, 0, ( (size_t)(ALGO->align) + 
                           (size_t)(A->ld+1) * (size_t)(nq+1) ) *
                         sizeof(double));                     
   for (i = 0; i < n; ++i)
   {
      HPL_indxg2lp(&local_i, &diagi, i, A->nb, A->nb, 0, GRID->nprow);
      HPL_indxg2lp(&local_j, &diagj, i, A->nb, A->nb, 0, GRID->npcol);

      if (GRID->myrow == diagi && GRID->mycol == diagj)
      {
         *(Mptr(invptrL, local_i, local_j, A->ld)) = 1.0;
         *(Mptr(invptrU, local_i, local_j, A->ld)) = 1.0;
      }
   }

   // if (GRID->myrow == 0 && GRID->mycol == 0)
   // {
   //    printf("\n======");
   //    for (int i = 0; i < A->mp; ++i)
   //    {
   //       printf("\n|");
   //       for(int j = 0; j < A->nq-1; ++j)
   //       {
   //          if (*Mptr(invptrL, i, j, A->ld) >= 0)
   //             printf(" ");
   //          printf("%f, ", *Mptr(invptrL, i, j, A->ld));
   //       }
   //       printf("|");
   //    }
   //    printf("\n------");
   // }

/*
 * calculate inverse of L and U, stored in invptrL and invptrU
 */
   HPL_dtrsm(HplColumnMajor, HplLeft, HplUpper, HplNoTrans, HplNonUnit, 
            mp, nq, 1.0, FACT, A->ld, invptrU, A->ld);
   HPL_dtrsm(HplColumnMajor, HplLeft, HplLower, HplNoTrans, HplUnit, 
            mp, nq, 1.0, FACT, A->ld, invptrL, A->ld);

/*
 * drive GMRES with reverse communication
 */
   while (1)
   {
      drive_dgmres(&n, &nloc, &m, &lwork, work,
                   irc, icntl, cntl, info, rinfo);

      /* irc[0] is the reverse communication indicator */
      switch (irc[0])
      {
      case 0:
         /* normal exit */
         break;
      case 1:
         /* matrix-vector product */
         ctran = 'N';
         dgemv_(&ctran, &mp, &nq, &one, A->A, lda, work+irc[1],
                &ione, &zero, work+irc[3], &ione, ione);
         break;
      case 2:
         /* first preconditioning */
         ctran = 'N';
         dgemv_(&ctran, &mp, &nq, &one, invptrL, lda, work+irc[1],
                &ione, &zero, work+irc[3], &ione, ione);
         break;
      case 3:
         /* second preconditioning */
         ctran = 'N';
         dgemv_(&ctran, &mp, &nq, &one, invptrU, lda, work+irc[1],
                &ione, &zero, work+irc[3], &ione, ione);
         break;
      case 4:
         /* dot products or in fact another matrix-vector product */
         ctran = 'Y';
         dgemv_(&ctran, &nloc, irc[4], &one, work+irc[1], &nloc, 
               work+irc[2], &ione, &zero, work+irc[3], &ione, ione);
         break;
      default:
         break;
      }
   }

/*
 * copy correcting vector back to overwrite residual vector
 */
   memcpy(R, work, nloc);

   if (work) free(work);
   if (vL) free(vL);
   if (vU) free(vU);


   printf("INFO[0] = %d\n", info[0]);

/*
 * End of HPL_pgmres
 */
}
