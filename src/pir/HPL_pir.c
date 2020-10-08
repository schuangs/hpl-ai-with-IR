/*
 * By Junkang Huang, August, 2020
 * 
 * based on the paper:
 *    A New  Analysis  Of Iterative Refinement And Its  Application To 
 *    Accurate Solution Of Ill Conditioned Sparse Linear Systems
 *    ---by Carson, Erin & Higham, Nicholas J., 2017
 */

/*
 * Include files
 */ 
#include "hpl.h"


#define IR 5

#define TOL 1e-14       /* Tolerance for GMRES residual */
#define MM 5            /* restart size for GMRES */
#define MAXIT 500       /* maximum number of GMRES iteration */

#ifdef STDC_HEADERS
void HPL_pir
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   HPL_T_pdmat *                    A,
   HPL_T_pmat *                     FA
)
#else
void HPL_pir
(GRID, ALGO, A, FA, MYROW, MYCOL, NPCOL)
   HPL_T_grid *                     GRID;
   HPL_T_palg *                     ALGO;
   HPL_T_pdmat *                    A;
   HPL_T_pmat *                     FA;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pir performs iterative refinement procesure to enhance the accur-
 * acy of  the solution  of linear system obtained  by LU factorization. 
 * Parallel  GMRES algorithm  is used  as the inner solver to solve  the 
 * inner correct equation Ad = r.
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
 * A       (local input/output)          HPL_T_pdmat *
 *         On entry, A points to the data structure containing the local
 *         array information. 
 * FA      (local input/output)          HPL_T_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information. It serves as lower precision version if A.
 *
 * ---------------------------------------------------------------------
 */ 

/*
 * .. Local Variables ..
 */
   int                i, j;
   int                mp, nq, n, nb, npcol, myrow, mycol, tarcol, info[3];
   double           * Bptr, *factors, *res, *d, *t, *pre;

/* ..
 * .. Executable Statements ..
 */
   mp = A->mp; nq = A->nq-1; n = A->n; nb = A->nb;
   npcol = GRID->npcol; mycol = GRID->mycol; myrow = GRID->myrow;
   
   Bptr  = Mptr( A->A,  0, nq, A->ld );

/*
 * allocate  space  for  factors, which  contains LU factors of  lower
 * precision and needed to be stored in double precision space. And pre
 * which contains the preconditioning matrices U-1 and L-1. pre in fact
 * contains the inverses of L and U that are stored in factors.
 */
   factors = (double*)malloc( mp * nq * sizeof(double) );
   pre     = (double*)malloc( mp * nq * sizeof(double) );

/*
 * Convert L and U into double precision for further operations
 */
   for (i = 0; i < mp; ++i)
   {
      for (j = 0; j < nq; ++j)
      {
         *Mptr(factors, i, j, A->ld) = (double)*Mptr(FA->A, i, j, FA->ld);
      }
   }
/*
 * Convert initial solution ( which is obtained through lower precision 
 * LU factorization) into higher precision
 */
   for(i = 0; i < nq; ++i)
   {
      *(A->X + i) = (double)(*(FA->X + i));
   }

/*
 * allocate space for residual vector, correction vectors.
 */
   res = (double*)malloc((size_t)mp * sizeof(double));
   d   = (double*)malloc((size_t)mp * sizeof(double));
   t   = (double*)malloc((size_t)nq * sizeof(double));

/*
 * parallelly calculate preconditioning matrix: pre = U-1L-1
 */

/* 
 * calculate U-1 by Gaussian Eliminations, that is to transform:
 *    U-1[U, I] => [I, U-1], in parallel.
 */
   HPL_my_pdtrsm(GRID, factors, pre, mp, nq, n);


/*
 * tarcol is the process column containing b in [ A | b ]
 */
   tarcol = HPL_indxg2p( n, nb, nb, 0, npcol );
/*
 * Iterative Refinement
 */
   for(i = 0; i < IR; ++i)
   {
      /* Calculate residual in double precision */
      if( mycol ==  tarcol)
      {
         memcpy(res, Bptr, mp*sizeof(double));
         HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, -HPL_rone,
                  A->A, A->ld, A->X, 1, HPL_rone, res, 1 );
      }
      else if( nq > 0 )
      {
         HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, -HPL_rone,
                  A->A, A->ld, A->X, 1, HPL_rzero, res, 1 );
      }
      else { for( j = 0; j < mp; j++ ) res[j] = HPL_rzero; }
      
      if (mp > 0)
         HPL_all_reduce( res, A->mp, HPL_DOUBLE, HPL_sum, GRID->row_comm );

   /* 
    * Solve correction  equation using preconditioned  GMRES  method in mix
    * precision.  But d is distributed just like res and replicated in each
    * process of a process row.
    */
      HPL_pgmres(GRID, ALGO, A, pre, res, t, TOL, MM, MAXIT);

   /*
    * in  order  to update solution, first  need to  transform  t into  the 
    * distributing pattern of X, which is stored in d.
    */


   /* update X with d */
      for (j = 0; j < nq; ++j)
      {
         *(A->X + j) += *(d + j);
      }
   }

   /* free dynamic memories */
   if (vF)  free(vF);
   if (d)   free(d);
   if (t)   free(t);
   if (pre) free(pre);
/*
 * End of HPL_pir
 */
}