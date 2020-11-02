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

# include <math.h>

#define IR 2

#define TOL 1e-13       /* Tolerance for GMRES residual */
#define PRE 1e-15       /* solution tolerance */
#define MM 50            /* restart size for GMRES */
#define MAXIT 10        /* maximum number of GMRES iteration */

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
   int                i, j, sizeA;
   int                mp, nq, n, nb, npcol, nprow, myrow, mycol, tarcol, info[3];
   double           * Bptr, *res, *d;
   void             * vptr;
   double             norm;
   HPL_T_pdmat        factors;

/* ..
 * .. Executable Statements ..
 */
   mp = A->mp; nq = A->nq-1; n = A->n; nb = A->nb;
   npcol = GRID->npcol; nprow = GRID->nprow;
   mycol = GRID->mycol; myrow = GRID->myrow;

   /* point to local rhs area */
   Bptr  = Mptr( A->A,  0, nq, A->ld );

   factors.mp = mp; factors.nq = A->nq; factors.n = n;
   factors.ld = A->ld; factors.nb = nb; factors.info = 0;

/*
 * allocate  space  for  factors, which  contains LU factors of  double
 * precision, which is in fact of lower precision, just stored in double. 
 */

   vptr = (void*)malloc( ( (size_t)(ALGO->align) + 
                           (size_t)(A->ld+1) * (size_t)(A->nq) ) *
                         sizeof(double) );

   factors.A  = (double *)HPL_PTR( vptr,((size_t)(ALGO->align) * sizeof(double) ) );
   factors.X  = Mptr( factors.A, 0, factors.nq, factors.ld );

   /* convert low-precision FA into double precision factors */
   sizeA = A->nq*A->ld;
   for(int i=0;i<sizeA;++i){
      *(factors.A+i)=(double)*(FA->A+i);
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
   res = (double*)malloc(mp * sizeof(double));
   d   = (double*)malloc(nq * sizeof(double));

/*
 * tarcol is the process column containing b
 */
   tarcol = HPL_indxg2p( n, nb, nb, 0, npcol ); 
   
/*
 * Iterative Refinement
 */
   for(i = 0; i < IR; ++i)
   {
      memset(res, 0, mp * sizeof(double));
      if (GRID->iam == 0)
         printf("IR Loop %d\n", i);
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
         HPL_all_reduce( res, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm );
      
      norm = 0;
      for (j = 0; j < mp; ++j)
      {
         norm += res[j]*res[j];
      }
      HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm );
      norm = sqrt(norm);

      if (norm < PRE)
         break;      

   /* 
    * Solve correction  equation using preconditioned  GMRES  method in mix
    * precision.  
    */
      memset(d, 0, nq * sizeof(double));
      HPL_pgmres(GRID, A, &factors, res, d, TOL, MM, MAXIT);
   /* 
    * update X with d
    */
      for (j = 0; j < nq; ++j)
      {
         *(A->X + j) += d[j];
      }
   }

   /* free dynamic memories */
   if (d)         free(d);
   if (vptr)      free(vptr);
   if (res)       free(res);

/*
 * End of HPL_pir
 */
}