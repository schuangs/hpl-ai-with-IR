/*
 * By Junkang Huang, August, 2020
 * 
 */

/*
 * Include files
 */ 
#include "hpl.h"


#define IR 5

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
   double           * Bptr, *Bptr1, *factors;
   void             * vF;

/* ..
 * .. Executable Statements ..
 */
   mp = A->mp; nq = A->nq-1; n = A->n; nb = A->nb;
   npcol = GRID->npcol; mycol = GRID->mycol; myrow = GRID->myrow;
   
   Bptr  = Mptr( A->A,  0, nq-1, A->ld );
   Bptr1 = Mptr( A->A,  0, nq,   A->ld);

/*
 * allocate  space for factors, which  contains LU factors of  lower 
 * precision and needed to be stored in double precision space.
 */
   vF = (void*)malloc( ( (size_t)(ALGO->align) + 
                           (size_t)(A->ld + 1) * (size_t)(nq+1) ) *
                         sizeof(double) );
   info[0] = (vF == NULL); info[1] = myrow; info[2] = mycol;
   (void) HPL_all_reduce( (void *)(info), 3, HPL_INT, HPL_max,
                          GRID->all_comm );
   if( info[0] != 0 )
   {
      /* some processes might have succeeded with allocation */
      if (vF) free(vF);
      return;
   }
/*
 * align pointer
 */
   factors = (double *)HPL_PTR( vF,
                               ((size_t)(ALGO->align) * sizeof(double) ) );
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
         memcpy(Bptr1,Bptr,A->mp*sizeof(double));
         HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, -HPL_rone,
                  A->A, A->ld, A->X, 1, HPL_rone, Bptr1, 1 );
      }
      else if( nq > 0 )
      {
         HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, -HPL_rone,
                  A->A, A->ld, A->X, 1, HPL_rzero, Bptr1, 1 );
      }
      else { for( j = 0; j < mp; j++ ) Bptr1[j] = HPL_rzero; }
      
      if (mp > 0)
         HPL_all_reduce( Bptr1, A->mp, HPL_DOUBLE, HPL_sum, GRID->row_comm );

      if (myrow == 0 && mycol == 0)
      {
         printf("\n");
         for (j = 0; j < mp; ++j)
         {
            printf("%f, ", *(Bptr1 + j));
         }
         printf("\n");
      }

   /* 
    * Solve correction equation using preconditioned GMRES method in mix
    * precision.  And  the  correction vector d will overwrite the space 
    * pointed by Bptr1 .
    */
   //    HPL_pgmres(GRID, ALGO, A, factors);

   // /*
   //  * update solution
   //  */
   //    for (j = 0; j < nq; ++j)
   //    {
   //       *(A->X + j) += *(Bptr1 + j);
   //    }
   }

   if (vF) free(vF);
/*
 * End of HPL_pir
 */
}