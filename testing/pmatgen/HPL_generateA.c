/*
 * By Junkang Huang, November, 2020
 * 
 * based on the paper:
 *  Matrices with Tunable Infinity-Norm Condition Number
 *  and No Need for Pivoting in LU Factorization
 *    --- by Fasi Massimiliano and Higham Nicholas J.
 *    [MIMS EPrint 2020.17, 
 *     Manchester Institute for Mathematical Sciences,
 *     The University of Manchester, UK, July 2020.]
 */

/*
 * Include files
 */ 
#include "hpl.h"

/*
 * Generate matrix A with given alpha and beta parameters.
 */
void generateA
(
   const HPL_T_grid *               GRID,
   const int                        N,
   const int                        NB,
   double *                         A,
   const int                        LDA,
   const double                     alpha,
   const double                     beta,
   const double                     scale
)
{
   int iloc, jloc, jdisp, i, j, k1, k2, mp, nq;
   int myrow = GRID->myrow; int mycol = GRID->mycol;
   int nprow = GRID->nprow; int npcol = GRID->npcol;
   double ab = alpha * beta;
   double a = - alpha;
   double b = - beta;

   /* local of local A */
   Mnumroc( mp, N, NB, NB, myrow, 0, nprow );
   Mnumroc( nq, N, NB, NB, mycol, 0, npcol );

   /* iloc, jloc are local indices */ 
   jloc = 0;
   /* i, j are the global indices */
   j = mycol * NB + 1;
   while (jloc < nq)
   {
      /* k1 iterates over different columns in a block column */
      for (k1 = 0; k1 < NB && j <= N; ++k1) 
      {
         jdisp = jloc * LDA;
         i = myrow * NB + 1;
         iloc = 0;
         while (iloc < mp) 
         {
            /* k2 iterates over different rows in a block */
            for (k2 = 0; k2 < NB && i <= N; ++k2) 
            {
               if (i > j)
                  A[jdisp + iloc] = (a + (j-1) * ab);
               else if (i == j)
                  A[jdisp + iloc] = (1 + (i-1) * ab);
               else // (j < i)
                  A[jdisp + iloc] = (b + (i - 1) * ab);
               
               A[jdisp + iloc] *= scale;
               iloc++;
               i++;
            }
            i += NB * (nprow - 1);
         }
         jloc++;
         j++;
      }
      j += NB * (npcol - 1);
   }
}


/*
 * Calculate proper ALPHA and BETA values according to size of matrix
 *   Here we set the condition number unchanged, and set ALPHA = BETA/2.
 *   According to the calculation result by Fasi M., we can set BETA 
 *   approximately to (V / N). V = 2.50 when k = 100, V = 5.19 when k = 10000.
 */
void calculate_ab
(
   double *                   alpha,
   double *                   beta,
   const int                  N
)
{
   /* with infinite condition number k = 100 */
   const double V = 2.50;
   *beta = V / N;
   *alpha = 0.5* *beta;
}

#ifdef STDC_HEADERS
void HPL_generateA
(
   const HPL_T_grid *               GRID,
   const int                        N,
   const int                        NB,
   double *                         A,
   const int                        LDA
)
#else
void HPL_generateA
( GRID, N, NB, A, LDA)
   const HPL_T_grid *               GRID;
   const int                        N;
   const int                        NB;
   double *                         A;
   const int                        LDA;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_generateA generates (or regenerates) a parallel N*N matrix A.
 *  Matrix A is generated based on the method proposed by Fasi M..
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 const HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 *
 * N       (global input)                const int
 *         On entry,  N specifies the global size of the matrix A.
 *         N must be at least zero.
 *
 * NB      (global input)                const int
 *         On entry,  NB specifies the blocking factor used to partition
 *         and distribute the matrix A. NB must be larger than one.
 *
 * A       (local output)                double *
 *         On entry,  A  points  to an array of dimension (LDA,LocQ(N)).
 *         On  exit,  this  array  contains  the  coefficients  of  the 
 *         generated matrix.
 *
 * LDA     (local input)                 const int
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least max(1,LocP(N)).
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
double alpha, beta, scale;

/* ..
 * .. Executable Statements ..
 */

/* calculate the alpha and beta parameters */
calculate_ab(&alpha, &beta, N);

/* scale the hole matrix as suggested */
scale = 65504. / 2;
generateA(GRID, N, NB, A, LDA, alpha, beta, scale);
}