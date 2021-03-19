/*
 * By Junkang Huang, November, 2020
 * 
 * based on the paper:
 * 
 *  Matrices with Tunable Infinity-Norm Condition Number
 *  and No Need for Pivoting in LU Factorization
 * 
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
 * Generate matrix A with alpha and beta parameters given.
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
/* Local variables */
   int iloc, jloc, jdisp, i, j, k1, k2, mp, nq;
   int myrow = GRID->myrow; int mycol = GRID->mycol;
   int nprow = GRID->nprow; int npcol = GRID->npcol;
   double ab = alpha * beta;
   double a = - alpha;
   double b = - beta;

   /* shape of local A */
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
   /* we choose infinite condition number k = 100 */
   const double V = 2.50;
   *beta = V / N;
   *alpha = 0.5* *beta;
}

/* 
 * HPL_generateA generates (or regenerates) a parallel N*N matrix A.
 *  Matrix A is generated based on the method proposed by Fasi M..
 */ 
void HPL_generateA
(
   const HPL_T_grid *               GRID,
   const int                        N,
   const int                        NB,
   double *                         A,
   const int                        LDA
) {
/*
 * .. Local Variables ..
 */
double alpha, beta, scale;

/* calculate the alpha and beta parameters */
calculate_ab(&alpha, &beta, N);

/* scale the hole matrix as suggested */
scale = 65504. / 2;
/* generate matrix A with alpha, beta and scale */
generateA(GRID, N, NB, A, LDA, alpha, beta, scale);

/* End of HPL_generateA() */
}


/*
 * Generate random right-hand side parallelly
 * 
 * The  pseudo-random  generator uses the linear congruential algorithm:
 * X(n+1) = (a * X(n) + c) mod m  as  described  in the  Art of Computer
 * Programming, Knuth 1973, Vol. 2.
 * 
 */
void HPL_generateB
(
   const HPL_T_grid *               GRID,
   const int                        N,
   const int                        NB,
   double *                         B,
   const int                        ISEED
) {
/*
 * .. Local Variables ..
 */
   int                        iadd [2], ia1  [2], ia2  [2], ib1  [2],
                              ic1  [2], ic2  [2],
                              iran1[2], iran2[2],
                              itmp1[2], itmp2[2],
                              jseed[2], mult [2];
   int                        ib, iblk, ik, jump1, jump2,
                              jump7, lmb, tarcol,
                              mblks, mp, mycol, myrow,
                              npcol, nprow;
/* ..
 * .. Executable Statements ..
 */
   (void) HPL_grid_info( GRID, &nprow, &npcol, &myrow, &mycol );

/*
 * tarcol is the process column containing b
 */
   tarcol = HPL_indxg2p( N, NB, NB, 0, npcol ); 
   if (mycol != tarcol) return;

   mult [0] = HPL_MULT0; mult [1] = HPL_MULT1;
   iadd [0] = HPL_IADD0; iadd [1] = HPL_IADD1;
   
   jseed[0] = ISEED;     jseed[1] = 0;
/*
 * Generate an M by N matrix starting in process (0,0)
 */
   Mnumroc( mp, N, NB, NB, myrow, 0, nprow );

/*
 * Local number of blocks and size of the last one
 */
   mblks = ( mp + NB - 1 ) / NB; lmb = mp - ( ( mp - 1 ) / NB ) * NB;
/*
 * Compute multiplier/adder for various jumps in random sequence
 */
   jump1 = 1;  jump2 = nprow * NB; jump7 = myrow * NB;

   HPL_xjumpm( jump1, mult, iadd, jseed, iran1, ia1,   ic1   );
   HPL_xjumpm( jump2, mult, iadd, iran1, itmp1, ia2,   ic2   );
   HPL_xjumpm( jump7, mult, iadd, iran1, iran1, itmp1, itmp2 );
   HPL_setran( 0, iran1 ); HPL_setran( 1, ia1 ); HPL_setran( 2, ic1 );
/*
 * Save value of first number in sequence
 */
   ib1[0] = iran1[0]; ib1[1] = iran1[1];

   for( iblk = 0; iblk < mblks; iblk++ )
   {
      ib = ( iblk == mblks - 1 ? lmb : NB );
      for( ik = 0; ik < ib; B++, ik++ ) *B = HPL_rand();
      HPL_jumpit( ia2, ic2, ib1, iran2 );
      ib1[0] = iran2[0]; ib1[1] = iran2[1];
   }

/*
 * End of HPL_generateB()
 */
}