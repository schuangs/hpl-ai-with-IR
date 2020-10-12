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


#define IR 10

#define TOL 1e-14       /* Tolerance for GMRES residual */
#define MM 50            /* restart size for GMRES */
#define MAXIT 100       /* maximum number of GMRES iteration */

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
   int                mp, nq, n, nb, npcol, nprow, myrow, mycol, tarcol, info[3];
   int                rmp, rnq;
   double           * Bptr, *factors, *res, *d, *preL, *preU;
   double             norm;

/* ..
 * .. Executable Statements ..
 */
   mp = A->mp; nq = A->nq-1; n = A->n; nb = A->nb;
   npcol = GRID->npcol; nprow = GRID->nprow;
   mycol = GRID->mycol; myrow = GRID->myrow;
   
   Bptr  = Mptr( A->A,  0, nq, A->ld );

/*
 * allocate  space  for  factors, which  contains LU factors of  double
 * precision, which is in fact of lower precision, just stored in double. 
 * And preU and preL which contains the preconditioning matrices U-1 and L-1.
 * preL adn preU in fact contains the inverses of L and U that are stored in factors.
 */

   /* IMPORTANT!:: There may be some incomplete blocks at the rightest column and 
    * bottom-est row. In order to treat them the same way as complete blocks, append
    * some columns of 0s to the right and bottom. So local mp and nq of those 
    * processes who contain incomplete blocks will be expanded into rmp and rnq.
    * 
    * First calculate rmp adn and rnq
    */
   rmp = mp; rnq = nq;
   if (n % nb != 0)
   {/* so there are some incomplete blocks */

      /* determine which process column/row contains the incomplete blocks */
      i = HPL_indxg2p((n/nb)*nb, nb, nb, 0, nprow);
      j = HPL_indxg2p((n/nb)*nb, nb, nb, 0, npcol);
      /* and expand mp and nq in these processes */
      if (myrow == i)
      {
         rmp += nb - n % nb;
      }
      if (mycol == j)
      {
         rnq += nb - n % nb;
      }
   }

   /* allocate spaces, and set 0s */
   factors = (double*)malloc( rmp * rnq * sizeof(double) );
   memset(factors, 0, rmp * rnq * sizeof(double));
   preL    = (double*)malloc( rmp * rnq * sizeof(double) );
   memset(preL,    0, rmp * rnq * sizeof(double));
   preU    = (double*)malloc( rmp * rnq * sizeof(double) );
   memset(preU,    0, rmp * rnq * sizeof(double));

/*
 * Convert L and U into double precision for further operations
 */
   for (j = 0; j < nq; ++j)
   {
      for (i = 0; i < mp; ++i)
      {
         *Mptr(factors, i, j, rmp) = (double)*Mptr(FA->A, i, j, FA->ld);
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
   d   = (double*)malloc((size_t)nq * sizeof(double));

/*
 * parallelly calculate preconditioning matrix: pre = U-1L-1
 */

/* 
 * calculate U-1 by Gaussian Eliminations, that is to transform:
 *    U-1[U, I] => [I, U-1], in parallel.
 */
   cal_pre(GRID, factors, preL, preU, mp, nq, rmp, n, nb);
   // if (GRID->iam == 0)
   // {
   //    printf("----------------------------------------------------------------\n");
   //    for (i = 0; i < mp; ++i)
   //    {
   //       for (j = 0; j < nq; ++j)
   //       {
   //          printf("%8f, ", *Mptr(preU, i, j, rmp));
   //       }
   //       printf("\n");
   //    }
   //    printf("================================================================\n");
   //    fflush(stdout);
   // }
/*
 * tarcol is the process column containing b
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
      
      norm = 0;
      for (j = 0; j < nq; ++j)
      {
         norm += res[j]*res[j];
      }
      HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm );
      if (GRID->iam == 0)
         printf("%.16f\n", norm);

   /* 
    * Solve correction  equation using preconditioned  GMRES  method in mix
    * precision.  
    */
      HPL_pgmres(GRID, A, preL, preU, rmp, res, d, TOL, MM, MAXIT);

   /* 
    * update X with d
    */
      for (j = 0; j < nq; ++j)
      {
         *(A->X + j) += d[j];
      }
      // if (GRID->iam == 0)
      // {
      //    for (int k = 0; k < mp; ++k)
      //    {
      //       printf("%.16f, ", (A->X)[k]);
      //    }
      //    printf("\n");
      // }
   }

   /* free dynamic memories */
   if (d)         free(d);
   if (preL)      free(preL);
   if (preU)      free(preU);
   if (factors)   free(factors);
   if (res)       free(res);

/*
 * End of HPL_pir
 */
}