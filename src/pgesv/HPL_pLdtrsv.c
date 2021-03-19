/*
 * By Junkang Huang, Dec. 2020
 * 
 * solving lower triangular system in parallel
 * following the implementation:
 *      HPL_pdtrsv() of HPL-2.3, which solves upper triangular system in parallel
 * 
 */
#include "hpl.h"


void HPL_pLdtrsv 
(
   HPL_T_grid *                GRID,
   HPL_T_pdmat *               AMAT
) {
/* 
 * Purpose
 * =======
 *
 * HPL_pLdtrsv solves an lower triangular system of linear equations.
 *  
 * The rhs is the last column of the N by N+1 matrix A. The solve starts
 * in the process  column owning the  1th  column of A, so the rhs b may
 * need to be moved one process column to the left at the beginning. The
 * routine therefore needs  a column  vector in every process column but
 * the one owning  b. The result is  replicated in all process rows, and
 * returned in XR, i.e. XR is of size nq = LOCq( N ) in all processes.
 *  
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * AMAT    (local input/output)          HPL_T_pdmat *
 *         On entry,  AMAT  points  to the data structure containing the
 *         local array information.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   MPI_Comm                   Ccomm, Rcomm;
   double                      * A=NULL, * Aprev=NULL, * Aptr, * XC=NULL,
                              * XR=NULL, * Xd=NULL, * Xdprev=NULL,
                              * W=NULL;
   int                        Alcol, Alrow, Anpprev, Anp, Anq, Bcol,
                              Cmsgid, GridIsNotPx1, GridIsNot1xQ, Rmsgid,
                              Wfr=0, colprev, kb, kbprev, lda, mycol,
                              myrow, n, n1, n1p, n1pprev=0, nb, npcol,
                              nprow, rowprev, tmp1, tmp2, np, nq, N;
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_DETAILED_TIMING
   HPL_ptimer( HPL_TIMING_PTRSV );
#endif
   if( ( N = n = AMAT->n ) <= 0 ) return;
   nb = AMAT->nb; lda = AMAT->ld; A = AMAT->A; XR = AMAT->X;

   (void) HPL_grid_info( GRID, &nprow, &npcol, &myrow, &mycol );
   Rcomm = GRID->row_comm; Rmsgid = MSGID_BEGIN_PTRSV;
   Ccomm = GRID->col_comm; Cmsgid = MSGID_BEGIN_PTRSV + 1;
   GridIsNot1xQ = ( nprow > 1 ); GridIsNotPx1 = ( npcol > 1 );
/*
 * Move the rhs in the process column owning the first column of A.
 */
   /* np and nq are the local dimensions of of A */
   Mnumroc( np, n, nb, nb, myrow, 0, nprow );
   Mnumroc( nq, n, nb, nb, mycol, 0, npcol );

   Anp = 0; Anq = 0;

   Alrow = 0;
   Alcol = 0;
   kb = Mmin(n, nb);

   Aptr = (double *)(A); XC = Mptr( Aptr, 0, nq, lda );
   Mindxg2p( n, nb, nb, Bcol, 0, npcol );

   if( ( np > 0 ) && ( Alcol != Bcol ) )
   {
      if( mycol == Bcol  )
      { (void) HPL_dsend( XC, np, Alcol, Rmsgid, Rcomm ); }
      else if( mycol == Alcol )
      { (void) HPL_drecv( XC, np, Bcol,  Rmsgid, Rcomm ); }
   }
   Rmsgid = ( Rmsgid + 2 >
              MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV : Rmsgid + 2 );
   
   /* other rhs area set to 0 */
   if( mycol != Alcol )
   { for( tmp1=0; tmp1 < np; tmp1++ ) XC[tmp1] = HPL_rzero; }
/*
 * Set up lookahead, that is to perform iteration on the first block column
 */

   /* n1 is the global distance between two neighboring blocks of the same
       process in a block row */
   n1 = ( npcol - 1 ) * nb; n1 = Mmax( n1, nb );

   /* allocate work space */
   if( np > 0 )
   {
      W = (double*)malloc( (size_t)(Mmin( n1, np )) * sizeof( double ) );
      if( W == NULL )
      { HPL_pabort( __LINE__, "HPL_pLdtrsv", "Memory allocation failed" ); }
      Wfr = 1;
   }

   /* previous version of parameters */
   Anpprev = Anp; Xdprev = Xd = XR; Aprev = Aptr;
   tmp1    = kb; tmp2 = Mmin( n - kb, n1 );
   MnumrocI( n1pprev, tmp2, Mmin( N, tmp1 ), nb, nb, myrow, 0, nprow );

   /* update local variables */
   if( myrow == Alrow ) { Anp += kb; }
   if( mycol == Alcol )
   {
      if( myrow == Alrow )
      {
         HPL_dtrsv( HplColumnMajor, HplLower, HplNoTrans, HplUnit,
                    kb, Aptr, lda, XC, 1 );
         HPL_dcopy( kb, XC, 1, Xd, 1 );
      }

      /* update local variables */
      Xd += kb; Anq += kb; Aptr += lda * kb;
   }

   /* update global variables */
   rowprev = Alrow; Alrow = MModAdd1( Alrow, nprow );
   colprev = Alcol; Alcol = MModAdd1( Alcol, npcol );
   kbprev  = kb; n -= kb;
   tmp1    = N - n + ( kb = Mmin(n, nb) ); tmp2 = Mmin( n - kb, n1 );
   MnumrocI( n1p, tmp2, Mmin( N, tmp1 ), nb, nb, myrow, 0, nprow );
/*
 * Start the operations
 */
   while( n > 0 )
   {
/*
 * Broadcast  (decreasing-ring)  of  previous solution block in previous
 * process column,  compute  partial update of current block and send it
 * to current process column.
 */
      if( mycol == colprev )
      {
/*
 * Send previous solution block in process row below
 */
         if( myrow == rowprev )
         {
            if( GridIsNot1xQ )
               (void) HPL_dsend( Xdprev, kbprev, MModAdd1( myrow, nprow ),
                                Cmsgid, Ccomm );
         }
         else
         {
            (void) HPL_drecv( Xdprev, kbprev, MModSub1( myrow, nprow ),
                             Cmsgid, Ccomm );
         } 
/*
 * Compute partial update of previous solution block and send it to cur-
 * rent column
 */
         if( n1pprev > 0 )
         {
            HPL_dgemv( HplColumnMajor, HplNoTrans, n1pprev, kbprev,
                       -HPL_rone, Aprev+Anp, lda, Xdprev, 1, HPL_rone,
                       XC+Anp, 1 );
            if( GridIsNotPx1 )
               (void) HPL_dsend( XC+Anp, n1pprev, Alcol, Rmsgid, Rcomm );
         }
/*
 * Finish  the (decreasing-ring) broadcast of the solution block in pre-
 * vious process column
 */
         if( ( myrow != rowprev ) &&
             ( myrow != MModSub1( rowprev, nprow ) ) )
            (void) HPL_dsend( Xdprev, kbprev, MModAdd1( myrow, nprow ),
                             Cmsgid, Ccomm );
      }
      else if( mycol == Alcol )
      {
/*
 * Current  column  receives  and accumulates partial update of previous
 * solution block
 */
         if( n1pprev > 0 )
         {
            (void) HPL_drecv( W, n1pprev, colprev, Rmsgid, Rcomm );
            HPL_daxpy( n1pprev, HPL_rone, W, 1, XC+Anp, 1 );
         }
      }
/*
 * Solve current diagonal block 
 */
      if( ( mycol == Alcol ) && ( myrow == Alrow ) )
      {
         HPL_dtrsv( HplColumnMajor, HplLower, HplNoTrans, HplUnit,
                    kb, Aptr+Anp, lda, XC+Anp, 1 );
         HPL_dcopy( kb, XC+Anp, 1, Xd, 1 );
      }
/*
 *  Finish previous update
 */
      if( ( mycol == colprev ) && ( ( tmp1 =  Anp + n1pprev ) < np ) )
         HPL_dgemv( HplColumnMajor, HplNoTrans, np - tmp1, kbprev, -HPL_rone,
                    Aprev + tmp1, lda, Xdprev, 1, HPL_rone, XC + tmp1, 1 );
                    
/*
 *  Save info of current step and update info for the next step
 */
      if( mycol == Alcol ) 
      {
         Aprev = Aptr; Aptr += lda * kb;
         Anq += kb; 
         Xdprev = Xd;  Xd = XR + Anq; 
      }
      if( myrow == Alrow ) { Anpprev = Anp; Anp += kb; }

      rowprev = Alrow; colprev = Alcol;
      n1pprev = n1p;   kbprev  = kb; n -= kb;
      Alrow = MModAdd1( Alrow, nprow ); Alcol = MModAdd1( Alcol, npcol );

      tmp1    = N - n + ( kb = Mmin(n, nb) ); tmp2 = Mmin( n - kb, n1 );
      MnumrocI( n1p, tmp2, Mmin( N, tmp1 ), nb, nb, myrow, 0, nprow );

      Rmsgid = ( Rmsgid+2 > MSGID_END_PTRSV ? 
                 MSGID_BEGIN_PTRSV   : Rmsgid+2 );
      Cmsgid = ( Cmsgid+2 > MSGID_END_PTRSV ?
                 MSGID_BEGIN_PTRSV+1 : Cmsgid+2 );
   }
/*
 * Replicate last solution block
 */
   if( mycol == colprev )
      (void) HPL_broadcast( (void *)Xdprev, kbprev, HPL_DOUBLE, rowprev,
                            Ccomm );

   if( Wfr  ) free( W  );
#ifdef HPL_DETAILED_TIMING
   HPL_ptimer( HPL_TIMING_PTRSV );
#endif
/*
 * End of HPL_pdtrsv
 */
}
