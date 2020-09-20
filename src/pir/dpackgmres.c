/* dPackgmres.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

   http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdio.h>
#include <math.h>
/* Table of constant values */

static int c__1 = 1;
static double c_b276 = 1.;
static double c_b280 = -1.;
static double c_b305 = 0.;

/* * */
/* *  Copyright (C) CERFACS 1998 */
/* * */
/* *  SOFTWARE LICENSE AGREEMENT NOTICE - THIS SOFTWARE IS BEING PROVIDED TO */
/* *  YOU BY CERFACS UNDER THE FOLLOWING LICENSE. BY DOWN-LOADING, INSTALLING */
/* *  AND/OR USING THE SOFTWARE YOU AGREE THAT YOU HAVE READ, UNDERSTOOD AND */
/* *  WILL COMPLY WITH THESE FOLLOWING TERMS AND CONDITIONS. */
/* * */
/* *  1 - This software program provided in source code format ("the " Source */
/* *  Code ") and any associated documentation (the " Documentation ") are */
/* *  licensed, not sold, to you. */
/* * */
/* *  2 - CERFACS grants you a personal, non-exclusive, non-transferable and */
/* *  royalty-free right to use, copy or modify the Source Code and */
/* *  Documentation, provided that you agree to comply with the terms and */
/* *  restrictions of this agreement. You may modify the Source Code and */
/* *  Documentation to make source code derivative works, object code */
/* *  derivative works and/or documentation derivative Works (called " */
/* *  Derivative Works "). The Source Code, Documentation and Derivative */
/* *  Works (called " Licensed Software ") may be used by you for personal */
/* *  and non-commercial use only. " non-commercial use " means uses that are */
/* *  not or will not result in the sale, lease or rental of the Licensed */
/* *  Software and/or the use of the Licensed Software in any commercial */
/* *  product or service. CERFACS reserves all rights not expressly granted */
/* *  to you. No other licenses are granted or implied. */
/* * */
/* *  3 - The Source Code and Documentation are and will remain the sole */
/* *  property of CERFACS. The Source Code and Documentation are copyrighted */
/* *  works. You agree to treat any modification or derivative work of the */
/* *  Licensed Software as if it were part of the Licensed Software itself. */
/* *  In return for this license, you grant CERFACS a non-exclusive perpetual */
/* *  paid-up royalty-free license to make, sell, have made, copy, distribute */
/* *  and make derivative works of any modification or derivative work you */
/* *  make of the Licensed Software. */
/* * */
/* *  4- The licensee shall acknowledge the contribution of the Source Code */
/* *  (using the reference [1]) in any publication of material dependent upon */
/* *  upon the use of the Source Code. The licensee shall use reasonable */
/* *  endeavours to notify the authors of the package of this publication. */
/* * */
/* *  [1] V. Frayss�, L. Giraud, S. Gratton, and J. Langou, A set of GMRES */
/* *    routines for double and complex arithmetics on high performance */
/* *    computers, CERFACS Technical Report TR/PA/03/3, public domain software */
/* *    available on www.cerfacs/algor/Softs, 2003 */
/* * */
/* *  5- CERFACS has no obligation to support the Licensed Software it is */
/* *  providing under this license. */
/* * */
/* *  THE LICENSED SOFTWARE IS PROVIDED " AS IS " AND CERFACS MAKE NO */
/* *  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE, */
/* *  BUT NOT LIMITATION, CERFACS MAKE NO REPRESENTATIONS OR WARRANTIES OF */
/* *  MERCHANTIBILY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF */
/* *  THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY THIRD */
/* *  PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. CERFACS WILL NOT */
/* *  BE LIABLE FOR ANY CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES, OR ANY */
/* *  OTHER RELIEF, OR FOR ANY CLAIM BY ANY THIRD PARTY, ARISING FROM YOUR */
/* *  USE OF THE LICENSED SOFTWARE. */
/* * */
/* *  6- For information regarding a commercial license for the Source Code */
/* *  and Documentation, please contact Mrs Campassens (campasse@cerfacs.fr) */
/* * */
/* *  7- This license is effective until terminated. You may terminate this */
/* *  license at any time by destroying the Licensed Software. */
/* * */
/* *    I agree all the terms and conditions of the above license agreement */
/* * */


int ErrorBlas(char *srname, int *info)
{
  /*  -- LAPACK auxiliary routine (version 2.0) --   
      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
      Courant Institute, Argonne National Lab, and Rice University   
      September 30, 1994   


      Purpose   
      =======   
      XERBLA  is an error handler for the LAPACK routines.   
      It is called by an LAPACK routine if an input parameter has an   
      invalid value.  A message is printed and execution stops.   
      Installers may consider modifying the STOP statement in order to   
      call system-specific exception-handling facilities.   

      Arguments   
      =========   

      SRNAME  (input) CHARACTER*6   
      The name of the routine which called XERBLA.   

      INFO    (input) INTEGER   
      The position of the invalid parameter in the parameter list   

      of the calling routine.   
      ===================================================================== 
  */

  printf("** On entry to %6s, parameter number %d had an illegal value\n",
	 srname, *info);
  return 0;
}

int Compare_char(char *ca, char *cb)
{
  /* System generated locals */
  int ret_val;

  /* Local variables */
  static int inta, intb,zcode;

  /*  -- LAPACK auxiliary routine (version 2.0) --   
      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
      Courant Institute, Argonne National Lab, and Rice University   
      January 31, 1994   
      Purpose   
      =======   
      LSAME returns .TRUE. if CA is the same letter as CB regardless of   
      case.   

      Arguments   
      =========   
      CA      (input) CHARACTER*1   
      CB      (input) CHARACTER*1   
      CA and CB specify the single characters to be compared.   
      ===================================================================== 
  */

  ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
  if (ret_val) {
    return ret_val;
  }

  /*     Now test for equivalence if both characters are alphabetic. */
  zcode = 'Z';
  /*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
	 machines, on which ICHAR returns a value with bit 8 set.   
	 ICHAR('A') on Prime machines returns 193 which is the same as   
	 ICHAR('A') on an EBCDIC machine. */

  inta = *(unsigned char *)ca;
  intb = *(unsigned char *)cb;

  if (zcode == 90 || zcode == 122) {

    /*        ASCII is assumed - ZCODE is the ASCII code of either lower or   
	      upper case 'Z'. */

    if (inta >= 97 && inta <= 122) {
      inta += -32;
    }
    if (intb >= 97 && intb <= 122) {
      intb += -32;
    }

  } else if (zcode == 233 || zcode == 169) {

    /*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
	      or   
	      upper case 'Z'. */

    if ((inta >= 129 && inta <= 137) || (inta >= 145 && inta <= 153) || (inta 
									 >= 162 && inta <= 169)) {
      inta += 64;
    }
    if ((intb >= 129 && intb <= 137) || (intb >= 145 && intb <= 153) || (intb 
									 >= 162 && intb <= 169)) {
      intb += 64;
    }

  } else if (zcode == 218 || zcode == 250) {

    /*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
	      e   
	      plus 128 of either lower or upper case 'Z'. */

    if (inta >= 225 && inta <= 250) {
      inta += -32;
    }
    if (intb >= 225 && intb <= 250) {
      intb += -32;
    }
  }
  ret_val = inta == intb;
  return ret_val;
}



int gcopy(int *n, double *dx, int *incx, 
	  double *dy, int *incy){
  int i__1;
  static int i, m, ix, iy, mp1;

  --dx;
  --dy;

  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }

  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= *n; ++i) {
    dy[iy] = dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

 L20:
  m = *n % 7;
  if (m == 0) {
    goto L40;
  }
  i__1 = m;
  for (i = 1; i <= m; ++i) {
    dy[i] = dx[i];
    /* L30: */
  }
  if (*n < 7) {
    return 0;
  }
 L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i = mp1; i <= *n; i += 7) {
    dy[i] = dx[i];
    dy[i+1] = dx[i+1];
    dy[i+2] = dx[i+2];
    dy[i+3] = dx[i+3];
    dy[i+4] = dx[i+4];
    dy[i+5] = dx[i+5];
    dy[i+6] = dx[i+6];
    /* L50: */
  }
  return 0;

}



int gconstantmut(int *n, double *da, double *dx, 
		 int *incx, double *dy, int *incy)
{
  int i__1;
  static int i, m, ix, iy, mp1;

  --dx;
  --dy;

  if (*n <= 0) {
    return 0;
  }
  if (*da == 0.) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }


  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= *n; ++i) {
    dy[iy] += *da * dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

 L20:
  m = *n % 4;
  if (m == 0) {
    goto L40;
  }
  i__1 = m;
  for (i = 1; i <= m; ++i) {
    dy[i] += *da * dx[i];
    /* L30: */
  }
  if (*n < 4) {
    return 0;
  }
 L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i = mp1; i <= *n; i += 4) {
    dy[i] += *da * dx[i];
    dy[i+1] += *da * dx[i+1];
    dy[i+2] += *da * dx[i+2];
    dy[i+3] += *da * dx[i+3];
    /* L50: */
  }
  return 0;
} 




int MatVectorProduct(char *trans, int *m, int *n, double *
		     alpha, double *a, int *lda, double *x, int *incx, 
		     double *beta, double *y, int *incy)
{

  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  static int info;
  static double temp;
  static int lenx, leny, i, j;
  static int ix, iy, jx, jy, kx, ky;


  /*  Purpose   
      =======   
      DGEMV  performs one of the matrix-vector operations   
      y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   
      where alpha and beta are scalars, x and y are vectors and A is an   
      m by n matrix.   

      Parameters   
      ==========   

      TRANS  - CHARACTER*1.   
      On entry, TRANS specifies the operation to be performed as   
      follows:   

      TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

      TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

      TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

      Unchanged on exit.   

      M      - INT.   
      On entry, M specifies the number of rows of the matrix A.   
      M must be at least zero.   
      Unchanged on exit.   

      N      - INT.   
      On entry, N specifies the number of columns of the matrix A. 
  
      N must be at least zero.   
      Unchanged on exit.   

      ALPHA  - DOUBLE PRECISION.   
      On entry, ALPHA specifies the scalar alpha.   
      Unchanged on exit.   

      A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
      Before entry, the leading m by n part of the array A must   
      contain the matrix of coefficients.   
      Unchanged on exit.   

      LDA    - INT.   
      On entry, LDA specifies the first dimension of A as declared 
  
      in the calling (sub) program. LDA must be at least   
      max( 1, m ).   
      Unchanged on exit.   

      X      - DOUBLE PRECISION array of DIMENSION at least   
      ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
      and at least   
      ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
      Before entry, the incremented array X must contain the   
      vector x.   
      Unchanged on exit.   

      INCX   - INT.   
      On entry, INCX specifies the increment for the elements of   
      X. INCX must not be zero.   
      Unchanged on exit.   

      BETA   - DOUBLE PRECISION.   
      On entry, BETA specifies the scalar beta. When BETA is   
      supplied as zero then Y need not be set on input.   
      Unchanged on exit.   

      Y      - DOUBLE PRECISION array of DIMENSION at least   
      ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
      and at least   
      ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
      Before entry with BETA non-zero, the incremented array Y   
      must contain the vector y. On exit, Y is overwritten by the 
  
      updated vector y.   

      INCY   - INT.   
      On entry, INCY specifies the increment for the elements of   
      Y. INCY must not be zero.   
      Unchanged on exit.   
      Level 2 Blas routine.   

      -- Written on 22-October-1986.   
      Jack Dongarra, Argonne National Lab.   
      Jeremy Du Croz, Nag Central Office.   
      Sven Hammarling, Nag Central Office.   
      Richard Hanson, Sandia National Labs.   

      Function Body */
  --x;
  --y;

  info = 0;
  if (! Compare_char(trans, "N") && ! Compare_char(trans, "T") && ! 
      Compare_char(trans, "C")) {
    info = 1;
  } else if (*m < 0) {
    info = 2;
  } else if (*n < 0) {
    info = 3;
  } else if (*lda < fmax((double)1,(double)(*m))) {
    info = 6;
  } else if (*incx == 0) {
    info = 8;
  } else if (*incy == 0) {
    info = 11;
  }
  if (info != 0) {
    ErrorBlas("DGEMV ", &info);
    return 0;
  }

  if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
    return 0;
  }

  /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
	 up the start points in  X  and  Y. */

  if (Compare_char(trans, "N")) {
    lenx = *n;
    leny = *m;
  } else {
    lenx = *m;
    leny = *n;
  }
  if (*incx > 0) {
    kx = 1;
  } else {
    kx = 1 - (lenx - 1) * *incx;
  }
  if (*incy > 0) {
    ky = 1;
  } else {
    ky = 1 - (leny - 1) * *incy;
  }

  /*     Start the operations. In this version the elements of A are   
	 accessed sequentially with one pass through A.   
	 First form  y := beta*y. */

  if (*beta != 1.) {
    if (*incy == 1) {
      if (*beta == 0.) {
	i__1 = leny;
	for (i = 1; i <= leny; ++i) {
	  y[i]  = 0.;
	  /* L10: */
	}
      } else {
	i__1 = leny;
	for (i = 1; i <= leny; ++i) {
	  y[i]  = *beta * y[i] ;
	  /* L20: */
	}
      }
    } else {
      iy = ky;
      if (*beta == 0.) {
	i__1 = leny;
	for (i = 1; i <= leny; ++i) {
	  y[iy]  = 0.;
	  iy += *incy;
	  /* L30: */
	}
      } else {
	i__1 = leny;
	for (i = 1; i <= leny; ++i) {
	  y[iy]  = *beta * y[iy];
	  iy += *incy;
	  /* L40: */
	}
      }
    }
  }
  if (*alpha == 0.) {
    return 0;
  }
  if (Compare_char(trans, "N")) {

    /*        Form  y := alpha*A*x + y. */
    jx = kx;
    if (*incy == 1) {
      i__1 = *n;
      for (j = 1; j <= *n; ++j) {
	if (x[jx]  != 0.) {
	  temp = *alpha * x[jx];
	  i__2 = *m;
	  for (i = 1; i <= *m; ++i) {
	    y[i] += temp * a[i-1 + (j-1)* ( *lda)];//A(i,j);
	    /* L50: */
	  }
	}
	jx += *incx;
	/* L60: */
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= *n; ++j) {
	if (x[jx] != 0.) {
	  temp = *alpha * x[jx];
	  iy = ky;
	  i__2 = *m;
	  for (i = 1; i <= *m; ++i) {
	    y[iy] += temp * a[i-1 + (j-1)* ( *lda)];
	    iy += *incy;
	    /* L70: */
	  }
	}
	jx += *incx;
	/* L80: */
      }
    }
  } else {

    /*        Form  y := alpha*A'*x + y. */
    jy = ky;
    if (*incx == 1) {
      i__1 = *n;
      for (j = 1; j <= *n; ++j) {
	temp = 0.;
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	  temp += a[i-1 + (j-1)* ( *lda)] * x[i];
	  /* L90: */
	}
	y[iy] += *alpha * temp;
	jy += *incy;
	/* L100: */
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= *n; ++j) {
	temp = 0.;
	ix = kx;
	i__2 = *m;
	for (i = 1; i <= *m; ++i) {
	  temp += a[i-1 + (j-1)* ( *lda)] * x[ix];
	  ix += *incx;
	  /* L110: */
	}
	y[iy] += *alpha * temp;
	jy += *incy;
	/* L120: */
      }
    }
  }
  return 0;
  /*     End of DGEMV . */
} 


double gnorm2(int *n, double *x, int *incx){

#define abs(x) ((x) >= 0 ? (x) : -(x))
  
  /* System generated locals */
  int i__1, i__2;
  double ret_val, d__1;
  double sqrt(double);

  /* Local variables */
  static double norm, scale, absxi;
  static int ix;
  static double ssq;

  //  -- This version written on 25-October-1982.   
  //   Modified on 14-October-1993 to inline the call to DLASSQ.   
  //  Sven Hammarling, Nag Ltd.   
    
  //    Function Body */

  --x;

  if (*n < 1 || *incx < 1) {
    norm = 0.;
  } else if (*n == 1) {
    norm = abs(x[1]);
  } else {
    scale = 0.;
    ssq = 1.;

    i__1 = (*n - 1) * *incx + 1;
    i__2 = *incx;
    for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
      if (x[ix] != 0.) {
	absxi = (d__1 = x[ix], abs(d__1));
	if (scale < absxi) {
	  /* Computing 2nd power */
	  d__1 = scale / absxi;
	  ssq = ssq * (d__1 * d__1) + 1.;
	  scale = absxi;
	} else {
	  /* Computing 2nd power */
	  d__1 = absxi / scale;
	  ssq += d__1 * d__1;
	}
      }
      /* L10: */
    }
    norm = scale * sqrt(ssq);
  }
  ret_val = norm;
  return ret_val;
}

int UpperLowerSolver(char *uplo, char *trans, char *diag, int *n, 
		     double *a, int *lda, double *x, int *incx){
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  static int info;
  static double temp;
  static int i, j;
  static int ix, jx, kx;
  static int nounit;

  /*  Purpose   
      =======   

      DTRSV  solves one of the systems of equations   
      A*x = b,   or   A'*x = b,   
      where b and x are n element vectors and A is an n by n unit, or   
      non-unit, upper or lower triangular matrix.   

      Parameters   
      ==========   

      UPLO   - CHARACTER*1.   
      On entry, UPLO specifies whether the matrix is an upper or   
      lower triangular matrix as follows:   

      UPLO = 'U' or 'u'   A is an upper triangular matrix.   

      UPLO = 'L' or 'l'   A is a lower triangular matrix.   

      Unchanged on exit.   

      TRANS  - CHARACTER*1.   
      On entry, TRANS specifies the equations to be solved as   
      follows:   

      TRANS = 'N' or 'n'   A*x = b.   

      TRANS = 'T' or 't'   A'*x = b.   

      TRANS = 'C' or 'c'   A'*x = b.   

      Unchanged on exit.   

      DIAG   - CHARACTER*1.   
      On entry, DIAG specifies whether or not A is unit   
      triangular as follows:   

      DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

      DIAG = 'N' or 'n'   A is not assumed to be unit   
      triangular.   

      Unchanged on exit.   

      N      - INT.   
      On entry, N specifies the order of the matrix A.   
      N must be at least zero.   
      Unchanged on exit.   

      A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
      Before entry with  UPLO = 'U' or 'u', the leading n by n   
      upper triangular part of the array A must contain the upper 
  
      triangular matrix and the strictly lower triangular part of 
  
      A is not referenced.   
      Before entry with UPLO = 'L' or 'l', the leading n by n   
      lower triangular part of the array A must contain the lower 
  
      triangular matrix and the strictly upper triangular part of 
  
      A is not referenced.   
      Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
      A are not referenced either, but are assumed to be unity.   
      Unchanged on exit.   

      LDA    - INT.   
      On entry, LDA specifies the first dimension of A as declared 
  
      in the calling (sub) program. LDA must be at least   
      max( 1, n ).   
      Unchanged on exit.   

      X      - DOUBLE PRECISION array of dimension at least   
      ( 1 + ( n - 1 )*abs( INCX ) ).   
      Before entry, the incremented array X must contain the n   
      element right-hand side vector b. On exit, X is overwritten 
  
      with the solution vector x.   

      INCX   - INT.   
      On entry, INCX specifies the increment for the elements of   
      X. INCX must not be zero.   
      Unchanged on exit.     

      -- Written on 22-October-1986.   
      Jack Dongarra, Argonne National Lab.   
      Jeremy Du Croz, Nag Central Office.   
      Sven Hammarling, Nag Central Office.   
      Richard Hanson, Sandia National Labs.   
  
      Function Body */

  --x;

  info = 0;
  if (! Compare_char(uplo, "U") && ! Compare_char(uplo, "L")) {
    info = 1;
  } else if (! Compare_char(trans, "N") && ! Compare_char(trans, "T") &&
	     ! Compare_char(trans, "C")) {
    info = 2;
  } else if (! Compare_char(diag, "U") && ! Compare_char(diag, "N")) {
    info = 3;
  } else if (*n < 0) {
    info = 4;
  } else if (*lda < fmax((double)1,(double)(*n))) {
    info = 6;
  } else if (*incx == 0) {
    info = 8;
  }
  if (info != 0) {
    ErrorBlas("DTRSV ", &info);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0) {
    return 0;
  }

  nounit=Compare_char(diag, "N");

  /*     Set up the start point in X if the increment is not unity. This   
	 will be  ( N - 1 )*INCX  too small for descending loops. */

  if (*incx <= 0) {
    kx = 1 - (*n - 1) * *incx;
  } else if (*incx != 1) {
    kx = 1;
  }

  /*     Start the operations. In this version the elements of A are   
	 accessed sequentially with one pass through A. */

  if (Compare_char(trans, "N")) {

    /*        Form  x := inv( A )*x. */

    if (Compare_char(uplo, "U")) {
      if (*incx == 1) {
	for (j = *n; j >= 1; --j) {
	  if (x[j] != 0.) {
	    if (nounit) {
	      x[j] /= a[j-1 + (j-1)* ( *lda)];
	    }
	    temp = x[j];
	    for (i = j - 1; i >= 1; --i) {
	      x[i] -= temp * a[i-1 + (j-1)* ( *lda)];
	      /* L10: */
	    }
	  }
	  /* L20: */
	}
      } else {
	jx = kx + (*n - 1) * *incx;
	for (j = *n; j >= 1; --j) {
	  if (x[jx] != 0.) {
	    if (nounit) {
	      x[jx] /= a[j-1 + (j-1)* ( *lda)];
	    }
	    temp = x[jx];
	    ix = jx;
	    for (i = j - 1; i >= 1; --i) {
	      ix -= *incx;
	      x[ix] -= temp * a[i-1 + (j-1)* ( *lda)];
	      /* L30: */
	    }
	  }
	  jx -= *incx;
	  /* L40: */
	}
      }
    } else {
      if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	  if (x[j] != 0.) {
	    if (nounit) {
	      x[j] /= a[j-1 + (j-1)* ( *lda)];
	    }
	    temp = x[j];
	    i__2 = *n;
	    for (i = j + 1; i <= *n; ++i) {
	      x[i] -= temp * a[i-1 + (j-1)* ( *lda)];
	      /* L50: */
	    }
	  }
	  /* L60: */
	}
      } else {
	jx = kx;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	  if (x[jx] != 0.) {
	    if (nounit) {
	      x[jx] /= a[j-1 + (j-1)* ( *lda)];
	    }
	    temp = x[jx];
	    ix = jx;
	    i__2 = *n;
	    for (i = j + 1; i <= *n; ++i) {
	      ix += *incx;
	      x[ix] -= temp * a[i-1 + (j-1)* ( *lda)];
	      /* L70: */
	    }
	  }
	  jx += *incx;
	  /* L80: */
	}
      }
    }
  } else {

    /*        Form  x := inv( A' )*x. */

    if (Compare_char(uplo, "U")) {
      if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	  temp = x[j];
	  i__2 = j - 1;
	  for (i = 1; i <= j-1; ++i) {
	    temp -= a[i-1 + (j-1)* ( *lda)] * x[i];
	    /* L90: */
	  }
	  if (nounit) {
	    temp /= a[j-1 + (j-1)* ( *lda)];
	  }
	  x[j]  = temp;
	  /* L100: */
	}
      } else {
	jx = kx;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	  temp = x[jx];
	  ix = kx;
	  i__2 = j - 1;
	  for (i = 1; i <= j-1; ++i) {
	    temp -= a[i-1 + (j-1)* ( *lda)] * x[ix];
	    ix += *incx;
	    /* L110: */
	  }
	  if (nounit) {
	    temp /= a[j-1 + (j-1)* ( *lda)];
	  }
	  x[jx] = temp;
	  jx += *incx;
	  /* L120: */
	}
      }
    } else {
      if (*incx == 1) {
	for (j = *n; j >= 1; --j) {
	  temp = x[j];
	  i__1 = j + 1;
	  for (i = *n; i >= j+1; --i) {
	    temp -= a[i-1 + (j-1)* ( *lda)] * x[i];
	    /* L130: */
	  }
	  if (nounit) {
	    temp /= a[j-1 + (j-1)* ( *lda)];
	  }
	  x[j] = temp;
	  /* L140: */
	}
      } else {
	kx += (*n - 1) * *incx;
	jx = kx;
	for (j = *n; j >= 1; --j) {
	  temp = x[jx];
	  ix = kx;
	  i__1 = j + 1;
	  for (i = *n; i >= j+1; --i) {
	    temp -= a[i-1 + (j-1)* ( *lda)] * x[ix];
	    ix -= *incx;
	    /* L150: */
	  }
	  if (nounit) {
	    temp /= a[j-1 + (j-1)* ( *lda)];
	  }
	  x[jx] = temp;
	  jx -= *incx;
	  /* L160: */
	}
      }
    }
  }
  return 0;
}


double Sign(double *a, double *b){
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

int GivensRot(double *da, double *db, double *c, 
	      double *s){
  
  static double c_b4 = 1.;
  /* System generated locals */
  double d__1, d__2;
  /* Builtin functions */
  double sqrt(double), d_sign(double *, double *);

  /* Local variables */
  static double r, scale, z, roe;
    
  /*     construct givens plane rotation.   
	 jack dongarra, linpack, 3/11/78. */

  roe = *db;
  if (abs(*da) > abs(*db)) {
    roe = *da;
  }
  scale = abs(*da) + abs(*db);
  if (scale != 0.) {
    goto L10;
  }
  *c = 1.;
  *s = 0.;
  r = 0.;
  z = 0.;
  goto L20;
 L10:
  /* Computing 2nd power */
  d__1 = *da / scale;
  /* Computing 2nd power */
  d__2 = *db / scale;
  r = scale * sqrt(d__1 * d__1 + d__2 * d__2);
  r = Sign(&c_b4, &roe) * r;
  *c = *da / r;
  *s = *db / r;
  z = 1.;
  if (abs(*da) > abs(*db)) {
    z = *s;
  }
  if (abs(*db) >= abs(*da) && *c != 0.) {
    z = 1. / *c;
  }
 L20:
  *da = r;
  *db = z;
  return 0;
}

int PlaneRotation(int *n, double *dx, int *incx, 
		  double *dy, int *incy, double *c, double *s)
{
  /* System generated locals */
  int i__1;
  /* Local variables */
  static int i;
  static double dtemp;
  static int ix, iy;

  /*     applies a plane rotation.   
	 jack dongarra, linpack, 3/11/78.   
	 modified 12/3/93, array(1) declarations changed to array(*)   */

  --dx;
  --dy;

  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }

  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= *n; ++i) {
    dtemp = *c * dx[ix]  + *s * dy[iy] ;
    dy[iy]  = *c * dy[iy]  - *s * dx[ix] ;
    dx[ix]  = dtemp;
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*       code for both increments equal to 1 */

 L20:
  i__1 = *n;
  for (i = 1; i <= *n; ++i) {
    dtemp = *c * dx[i]  + *s * dy[i];
    dy[i]  = *c * dy[i] - *s * dx[i];
    dx[i]  = dtemp;
    /* L30: */
  }
  return 0;
}

/* Subroutine */ int drive_dgmres(int *n, int *nloc, int *m, 
				  int *lwork, double *work, int *irc, int *icntl, 
				  double *cntl, int *info, double *rinfo){
  /* Initialized data */
  static int icheck = 0;

  /* System generated locals */
  double r__1;
    
  /* Builtin functions */
  double sqrt(double);

  /* Local variables */
  static int xcurrent, ycurrent;
  static double sa, sb;
  static double rc, rn, rx;
  static int newrestart;
  static double spa, spb;
  static int ierr, bptr, hptr, vptr, wptr, xptr, r0ptr, iwarn, ihist;
  static int rotcos, dotptr, rotsin, comprsd, sizewrk;


  /*  Purpose */
  /*  ======= */
  /*    drive_dgmres is the driver routine for solving the linear system */
  /*  Ax = b using the *  Generalized Minimal Residual iterative method */
  /*  with preconditioning. */
  /*  This solver is implemented with a reverse communication scheme: control */
  /*  is returned to the user for computing the (preconditioned) */
  /*  matrix-vector product. */
  /*  See the User's Guide for an example of use. */


  /* Written : June 1996 */
  /* Authors : Luc Giraud, Serge Gratton, V. Fraysse */
  /*             Parallel Algorithms - CERFACS */

  /* Updated : April 1997 */
  /* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton */
  /*             Parallel Algorithms - CERFACS */

  /* Updated : June 1998 */
  /* Authors : Valerie Fraysse, Luc Giraud, Serge Gratton */
  /*             Parallel Algorithms - CERFACS */
  /* Purpose : Make clear that the warning and error messages come from the */
  /*           dgmres modules. */

  /* Updated : December 2002 - L. Giraud, J.Langou */
  /* Purpose : Add the capability to avoid explicit residual calculation at restart */

  /* Updated : March 2005 - L. Giraud */
  /* Purpose : small bug in the computation of the restart if too large vs the */
  /*           size of the workspace */

  /*  Arguments */
  /*  ========= */

  /*  n      (input) INT. */
  /*          On entry, the dimension of the problem. */
  /*          Unchanged on exit. */

  /*  nloc   (input) INT. */
  /*          On entry, the dimension of the local problem. */
  /*          In a parallel distributed envirionment, this corresponds */
  /*          to the size of the subset of entries of the right hand side */
  /*          and solution allocated to the calling process. */
  /*          Unchanged on exit. */


  /*  m      (input) INT */
  /*          Restart parameter, <= N. This parameter controls the amount */
  /*          of memory required for matrix H (see WORK and H). */
  /*          Unchanged on exit. */

  /*  lwork  (input) INT */
  /*          size of the workspace */
  /*          if (icntl(5) = 0 or 1 ) */
  /*            lwork >= m*m + m*(n+5) + 5*n+2, if icntl(8) = 1 */
  /*            lwork >= m*m + m*(n+5) + 6*n+2, if icntl(8) = 0 */
  /*          if (icntl(5) = 2 or 3 ) */
  /*            lwork >= m*m + m*(n+5) + 5*n+m+1, if icntl(8) = 1 */
  /*            lwork >= m*m + m*(n+5) + 6*n+m+1, if icntl(8) = 0 */

  /*  work   (workspace) double precision/double precision array, length lwork */
  /*          work contains the required vector and matrices stored in the */
  /*          following order : */
  /*            x  (n,1)       : computed solution. */
  /*            b  (n,1)       : right hand side. */
  /*            r0 (n,1)       : vector workspace. */
  /*            w  (n,1)       : vector workspace. */
  /*            V  (n,m)       : Krylov basis. */
  /*            H  (m+1,m+1)   : Hessenberg matrix (full storage). */
  /*            yCurrent (m,1) : solution of the current LS */
  /*            xCurrent (n,1) : current iterate */
  /*            rotSin (m,1)   : Sine of the Givens rotation */
  /*            rotCos (m,1)   : Cosine of the Givens rotation */

  /*  irc     (input/output) INT array. length 5 */
  /*            irc(1) : REVCOM   used for reverse communication */
  /*                             (type of external operation) */
  /*            irc(2) : COLX     used for reverse communication */
  /*            irc(3) : COLY     used for reverse communication */
  /*            irc(4) : COLZ     used for reverse communication */
  /*            irc(5) : NBSCAL   used for reverse communication */

  /*  icntl   (input) INT array. length 7 */
  /*            icntl(1) : stdout for error messages */
  /*            icntl(2) : stdout for warnings */
  /*            icntl(3) : stdout for convergence history */
  /*            icntl(4) : 0 - no preconditioning */
  /*                       1 - left preconditioning */
  /*                       2 - right preconditioning */
  /*                       3 - double side preconditioning */
  /*                       4 - error, default set in Init */
  /*            icntl(5) : 0 - modified Gram-Schmidt */
  /*                       1 - iterative modified Gram-Schmidt */
  /*                       2 - classical Gram-Schmidt */
  /*                       3 - iterative classical Gram-Schmidt */
  /*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set) */
  /*                       1 - user supplied initial guess */
  /*            icntl(7) : maximum number of iterations */
  /*            icntl(8) : 0 - use recurence formula at restart */
  /*                       1 - default compute the true residual at each restart */

  /*  cntl    (input) double precision array, length 5 */
  /*            cntl(1) : tolerance for convergence */
  /*            cntl(2) : scaling factor for normwise perturbation on A */
  /*            cntl(3) : scaling factor for normwise perturbation on b */
  /*            cntl(4) : scaling factor for normwise perturbation on the */
  /*                      preconditioned matrix */
  /*            cntl(5) : scaling factor for normwise perturbation on */
  /*                      preconditioned right hand side */

  /*  info    (output) INT array, length 3 */
  /*            info(1) :  0 - normal exit */
  /*                      -1 - n < 1 */
  /*                      -2 - m < 1 */
  /*                      -3 - lwork too small */
  /*                      -4 - convergence not achieved after icntl(7) iterations */
  /*                      -5 - precondition type not set by user */
  /*            info(2) : if info(1)=0 - number of iteration to converge */
  /*                      if info(1)=-3 - minimum workspace size necessary */
  /*            info(3) : optimal size for the workspace */

  /*  rinfo   (output) double precision array, length 2 */
  /*            if info(1)=0 */
  /*              rinfo(1) : backward error for the preconditioned system */
  /*              rinfo(2) : backward error for the unpreconditioned system */

  /* Input variables */
  /* --------------- */
  /* Output variables */
  /* ---------------- */
  /* Input/Output variables */
  /* ---------------------- */
  /* Local variables */
  /* --------------- */

  /* Parameter adjustments */
  --rinfo;
  --info;
  --cntl;
  --icntl;
  --irc;
  --work;

  // printf("index = %d, %.16f, ", irc[4], *(work+irc[4]+1));

  /* Function Body */

  /*       intrinsic ifix, float */

  /*       Executable statements : */

  ierr = icntl[1];
  iwarn = icntl[2];
  ihist = icntl[3];
  comprsd = icntl[8];

  if (ierr < 0) {
    ierr = 6;
  }

  if (comprsd == 1) {
    sizewrk = *m * *m + *m * (*nloc + 5) + *nloc * 5 + 1;
  } else {
    sizewrk = *m * *m + *m * (*nloc + 5) + *nloc * 6 + 1;
  }

  if (icheck == 0) {
   
    /* Check the value of the arguments */
    if (*n < 1 || *nloc < 1) {
      printf(" ERROR GMRES: \n");
      printf("     N < 1 \n");
      info[1] = -1;
      irc[1] = 0;
      return 0;
    }
    if (*m < 1) {
      printf(" ERROR GMRES : \n");
      printf("     M < 1 \n");
      info[1] = -2;
      irc[1] = 0;
      return 0;
    }
    if (icntl[4] != 0 && icntl[4] != 1 && icntl[4] != 2 && icntl[4] != 3) 
      {
	printf(" ERROR GMRES : \n");
	printf("     Undefined preconditioner \n");
	info[1] = -5;
	irc[1] = 0;
	return 0;
      }

    if (icntl[5] < 0 || icntl[5] > 3) {
      icntl[5] = 0;
      if (iwarn != 0) {
	printf(" Warning GMRES : \n");
	printf("     Undefined orthogonalisation \n");
	printf("     Default MGS \n");
      }
    }

    if (icntl[5] == 2 || icntl[5] == 3) {
      /* the workspace should be large enough to store the m dot-products */
      sizewrk += *m;
    } else {
      ++sizewrk;
    }

    if (iwarn != 0) {
      printf(" Warning GMRES : \n");
      printf(" For M = %d: \n",*m);
      printf(" optimal value : \n");
      printf(" for LWORK =  : %d \n",sizewrk);
    }

    if (icntl[6] != 0 && icntl[6] != 1) {
      icntl[6] = 0;
      if (iwarn != 0) {
	printf(" Warning GMRES : \n");
	printf(" Undefined intial guess \n");
	printf(" Default x0 = 0  \n");
      }
    }
    if (icntl[7] <= 0) {
      icntl[7] = *n;
      if (iwarn != 0) {
	printf(" Warning GMRES : \n");
	printf(" Negative max number of iterations \n");
	printf(" Default N  \n");
      }
    }
    if (icntl[8] != 0 && icntl[8] != 1) {
      icntl[8] = 1;
      printf(" Warning GMRES : \n");
      printf("  Undefined strategy for the residua \n");
      printf(" at restart  \n");
      printf(" Default 1  \n");
    }
    /* Check if the restart parameter is correct and if the size of the */
    /*  workspace is big enough for the restart. */
    /* If not try to fix correctly the parameters */

    if (*m > *n || *lwork < sizewrk) {
      if (*m > *n) {
	*m = *n;
	if (iwarn != 0) {
	  printf(" Warning GMRES : \n");
	  printf(" Parameter M bigger than N \n");
	  printf(" New value for M %d : \n",*m);
	}
	if (comprsd == 1) {
	  sizewrk = *m * *m + *m * (*nloc + 5) + *nloc * 5 + 1;
	} else {
	  sizewrk = *m * *m + *m * (*nloc + 5) + *nloc * 6 + 1;
	}
	if (icntl[5] == 2 || icntl[5] == 3) {
	  /* the workspace should be large enough to store the m dot-products */
	  sizewrk += *m;
	} else {
	  ++sizewrk;
	}
      }
      if (*lwork < sizewrk && *n == *nloc) {
	/* Compute the maximum size of the m according to the memory space */
	rn = (double) (*n);
	rx = rn + 5.f;
	rc = rn * 5.f + 1 - (double) (*lwork);

	/* Update the linear part of the second order equation to be solved */
	if (icntl[5] == 2 || icntl[5] == 3) {
	  rx += 1;
	}
	/* Update the constant part of the second order equation to be solved */

	if (icntl[8] == 0) {
	  rc += rn;
	}
	/* Solution of the the second order equation */
	/* Computing 2nd power */
	r__1 = rx;
	newrestart = (int) ((-rx + sqrt(r__1 * r__1 - rc * 4.f)) /
			    2.f);
	if (newrestart > 0) {
	  *m = newrestart;
	  if (iwarn != 0) {
	    printf(" Warning GMRES : \n");
	    printf(" Workspace too small for \n");
	    printf(" New value for M %d : \n",*m);
	  }
	} else {
	  printf(" ERROR GMRES : \n");
	  printf(" Not enough space for the problem \n");
	  printf(" the space does not permit any m: %d \n",*m);
	  info[1] = -3;
	  irc[1] = 0;
	  return 0;
	}
      }
      if (*lwork < sizewrk && *n != *nloc) {
	printf(" ERROR GMRES : \n");
	printf(" Not enough space for the problem \n");
	info[1] = -3;
	irc[1] = 0;
	return 0;
      }
    }

    info[3] = sizewrk;
    icheck = 1;
    /* save the parameters the the history file */

    if (ihist != 0) {
      printf(" \n \n ");
      printf(" CONVERGENCE HISTORY FOR GMRES : \n");
      printf(" \n");
      printf(" Errors are displayed in unit: %d \n",ierr);
      if (iwarn == 0) {
	printf(" Warnings are not displayed \n");
      } else {
	printf(" Warnings are displayed in unit: %d \n",iwarn);
      }
      printf(" matrix size: %d \n",*n);
      printf(" Local matrix size: %d \n",*nloc);
      printf(" Restart: %d \n",*m);
      if (icntl[4] == 0) {
	printf(" No preconditioning  \n");
      } else if (icntl[4] == 1) {
	printf(" Left preconditioning  \n");
      } else if (icntl[4] == 2) {
	printf(" Right preconditioning  \n");
      } else if (icntl[4] == 3) {
	printf(" Left and Right preconditioning \n");
      }
      if (icntl[5] == 0) {
	printf(" Modified Gram-Schmidt  \n");
      } else if (icntl[5] == 1) {
	printf(" Iterative modified Gram-Schmidt  \n");
      } else if (icntl[5] == 2) {
	printf(" Classical Gram-Schmidt  \n");
      } else {
	printf(" Iterative Classical Gram-Schmidt \n");
      }
      if (icntl[6] == 0) {
	printf(" Default initial guess x_0 = 0 \n");
      } else {
	printf(" User supplied initial guess \n");
      }
      if (icntl[8] == 1) {
	printf(" True residual computed at restart \n");
      } else {
	printf(" Recurrence residual at restart \n");
      }
      printf(" Maximum number of iterations: %d \n",icntl[7]);
      printf(" Tolerance for convergence: %.5e \n",cntl[1]);
      printf(" Backward error on the unpreconditioned system Ax= b \n");
      sa = cntl[2];
      sb = cntl[3];
      if (sa == 0. && sb == 0.) {
	printf(" the residual is normalised by ||b|| \n");
      } else {
	printf(" the residual is normalised by  %.2e * ||x|| + %.2e \n",sa,sb);
      }
      spa = cntl[4];
      spb = cntl[5];
      printf(" Backward error on the preconditioned system (P1)A(P2)x= (P1)b \n");
      if (spa == 0. && spb == 0.) {
	printf(" the residual is normalised by ||b|| \n");
      } else {
	printf(" the residual is normalised by  %.2e * ||(P2)y|| + %.2e \n",spa,spb);
      }
      printf(" Optimal size for the local workspace: %d \n",info[3]);
    }

  }
    
  /* setup some pointers on the workspace */
  xptr = 1;
  bptr = xptr + *nloc;
  r0ptr = bptr + *nloc;
  wptr = r0ptr + *nloc;
  vptr = wptr + *nloc;
  if (comprsd == 1) {
    hptr = vptr + *m * *nloc;
  } else {
    hptr = vptr + (*m + 1) * *nloc;
  }
  dotptr = hptr + (*m + 1) * (*m + 1);
  if (icntl[5] == 2 || icntl[5] == 3) {
    ycurrent = dotptr + *m;
  } else {
    ycurrent = dotptr + 1;
  }
  xcurrent = ycurrent + *m;
  rotsin = xcurrent + *nloc;
  rotcos = rotsin + *m;
    
  dgmres_(nloc, m, &work[bptr], &work[xptr], &work[hptr], &work[wptr], &
	  work[r0ptr], &work[vptr], &work[dotptr], &work[ycurrent], &work[
									  xcurrent], &work[rotsin], &work[rotcos], &irc[1], &icntl[1], &
	  cntl[1], &info[1], &rinfo[1]);

  if (irc[1] == 0) {
    icheck = 0;
  }
  return 0;
} 


/* Subroutine */ int dgmres_(int *n, int *m, double *b, 
			     double *x, double *h__, double *w, double *r0, 
			     double *v, double *dot, double *ycurrent, double *
			     xcurrent, double *rotsin, double *rotcos, int *irc, 
			     int *icntl, double *cntl, int *info, double *rinfo)
{
  /* Initialized data */
   
  static int retlbl = 0;

  /* System generated locals */
  int h_dim1, h_offset, v_dim1, v_offset, i__1;
  double d__1, d__2;

  /* Builtin functions */
  double sqrt(double);

  /* Local variables */
  static double auxhjp1j, dnormres;
  static int typeprec, j, initguess;
  static double be, bn;
  static int jh;
  static double sa, sb, bea, dvi, spa, spb, aux, beta, truenormres, 
    dloo, temp;
  static int bptr;
  static int hptr, vptr, wptr, xptr, yptr;
  static int r0ptr;
  static int iwarn, ihist;
  static double auxhjj, dnormw, dnormx;
  static int northo, dotptr, xcuptr, comprsd, itermax, iorthog, iterout;



  /*  Purpose */
  /*  ======= */
  /*  dgmres solves the linear system Ax = b using the */
  /*  Generalized Minimal Residual iterative method */

  /* When preconditioning is used we solve : */
  /*     M_1^{-1} A M_2^{-1} y = M_1^{-1} b */
  /*     x = M_2^{-1} y */

  /*   Convergence test based on the normwise backward error for */
  /*  the preconditioned system */

  /* Written : June 1996 */
  /* Authors : Luc Giraud, Serge Gratton, V. Fraysse */
  /*             Parallel Algorithms - CERFACS */

  /* Updated : April 1997 */
  /* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton */
  /*             Parallel Algorithms - CERFACS */

  /* Updated : March 1998 */
  /* Purpose : Pb with F90 on DEC ws */
  /*           cure : remove "ZDSCAL" when used to initialize vectors to zero */

  /* Updated : May 1998 */
  /* Purpose : r0(1) <-- r0'r0 : pb when used with DGEMV for the dot product */
  /*           cure : w(1) <--  r0'r0 */

  /* Updated : June 1998 */
  /* Purpose : Make clear that the warning and error messages come from the */
  /*           dgmres modules. */

  /* Updated : February 2001 - L. Giraud */
  /* Purpose : In complex version, initializations to zero performed  in complex */
  /*           arithmetic to avoid implicit conversion by the compiler. */

  /* Updated : July 2001 - L. Giraud, J. Langou */
  /* Purpose : Avoid to compute the approximate solution at each step of */
  /*           the Krylov space construction when spA is zero. */

  /* Updated : November 2002 - S. Gratton */
  /* Purpose : Use Givens rotations conform to the classical definition. */
  /*           No impact one the convergence history. */

  /* Updated : November 2002 - L. Giraud */
  /* Purpose : Properly handle the situation when the convergence is obtained */
  /*           exactly at the "IterMax" iteration */

  /* Updated : December 2002 - L. Giraud, J.Langou */
  /* Purpose : Add the capability to avoid explicit residual calculation at restart */

  /* Updated : January  2003 - L. Giraud, S. Gratton */
  /* Purpose : Use Givens rotations from BLAS. */

  /* Updated : March    2003 - L. Giraud */
  /* Purpose : Set back retlbl to zero, if initial guess is solution */
  /*           or right-hand side is zero */

  /* Updated : September 2003 - L. Giraud */
  /* Purpose : Include room in the workspace to store the results of the dot products */
  /*           Fix the bugs that appeared when M > Nloc */

  /*  Arguments */
  /*  ========= */

  /*  n       (input) INT. */
  /*           On entry, the dimension of the problem. */
  /*           Unchanged on exit. */

  /*  m        (input) INT */
  /*           Restart parameter, <= N. This parameter controls the amount */
  /*           of memory required for matrix H (see WORK and H). */
  /*           Unchanged on exit. */

  /*  b        (input) double precision/double precision */
  /*           Right hand side of the linear system. */

  /*  x        (output) double precision/double precision */
  /*           Computed solution of the linear system. */

  /*  H        (workspace)  double precision/double precision */
  /*           Hessenberg matrix built within dgmres */

  /*  w        (workspace)  double precision/double precision */
  /*           Vector used as temporary storage */

  /*  r0       (workspace)  double precision/double precision */
  /*           Vector used as temporary storage */

  /*  V        (workspace)  double precision/double precision */
  /*           Basis computed by the Arnoldi's procedure. */

  /*  dot      (workspace) double precision/double precision */
  /*           Store the results of the dot product calculation */

  /*  yCurrent (workspace) double precision/double precision */
  /*           solution of the current LS */

  /*  xCurrent (workspace) double precision/double precision */
  /*           current iterate */

  /*  rotSin   (workspace) double precision/double precision */
  /*           Sine of the Givens rotation */

  /*  rotCos   (workspace) double precision */
  /*           Cosine of the Givens rotation */

  /*  irc      (input/output) INT array. length 3 */
  /*             irc(1) : REVCOM   used for reverse communication */
  /*                              (type of external operation) */
  /*             irc(2) : COLX     used for reverse communication */
  /*             irc(3) : COLY     used for reverse communication */
  /*             irc(4) : COLZ     used for reverse communication */
  /*             irc(5) : NBSCAL   used for reverse communication */

  /*  icntl    (input) INT array. length 7 */
  /*             icntl(1) : stdout for error messages */
  /*             icntl(2) : stdout for warnings */
  /*             icntl(3) : stdout for convergence history */
  /*             icntl(4) : 0 - no preconditioning */
  /*                        1 - left preconditioning */
  /*                        2 - right preconditioning */
  /*                        3 - double side preconditioning */
  /*                        4 - error, default set in Init */
  /*             icntl(5) : 0 - modified Gram-Schmidt */
  /*                        1 - iterative modified Gram-Schmidt */
  /*                        2 - classical Gram-Schmidt */
  /*                        3 - iterative classical Gram-Schmidt */
  /*             icntl(6) : 0 - default initial guess x_0 = 0 (to be set) */
  /*                        1 - user supplied initial guess */
  /*             icntl(7) : maximum number of iterations */
  /*             icntl(8) : 1 - default compute the true residual at each restart */
  /*                        0 - use recurence formula at restart */

  /*  cntl     (input) double precision array, length 5 */
  /*             cntl(1) : tolerance for convergence */
  /*             cntl(2) : scaling factor for normwise perturbation on A */
  /*             cntl(3) : scaling factor for normwise perturbation on b */
  /*             cntl(4) : scaling factor for normwise perturbation on the */
  /*                       preconditioned matrix */
  /*             cntl(5) : scaling factor for normwise perturbation on */
  /*                       preconditioned right hand side */

  /*  info     (output) INT array, length 2 */
  /*             info(1) :  0 - normal exit */
  /*                       -1 - n < 1 */
  /*                       -2 - m < 1 */
  /*                       -3 - lwork too small */
  /*                       -4 - convergence not achieved after icntl(7) iterations */
  /*                       -5 - precondition type not set by user */
  /*             info(2) : if info(1)=0 - number of iteration to converge */
  /*                       if info(1)=-3 - minimum workspace size necessary */
  /*             info(3) : optimal size for the workspace */

  /* rinfo     (output) double precision array, length 2 */
  /*             if info(1)=0 */
  /*               rinfo(1) : backward error for the preconditioned system */
  /*               rinfo(2) : backward error for the unpreconditioned system */



  /* Reverse communication variables */
  /* ------------------------------- */
  /* Parameter adjustments */
  v_dim1 = *n;
  v_offset = 1 + v_dim1;
  v -= v_offset;
  h_dim1 = *m + 1;
  h_offset = 1 + h_dim1;
  h__ -= h_offset;
  --b;
  --x;
  --w;
  --r0;
  --dot;
  --ycurrent;
  --xcurrent;
  --rotsin;
  --rotcos;
  --irc;
  --icntl;
  --cntl;
  --info;
  --rinfo;

  /* Function Body */

  /*       Executable statements */

  /* setup some pointers on the workspace */
  xptr = 1;
  bptr = xptr + *n;
  r0ptr = bptr + *n;
  wptr = r0ptr + *n;
  vptr = wptr + *n;
  if (icntl[8] == 1) {
    hptr = vptr + *m * *n;
  } else {
    hptr = vptr + (*m + 1) * *n;
  }
  dotptr = hptr + (*m + 1) * (*m + 1);
  if (icntl[5] == 2 || icntl[5] == 3) {
    yptr = dotptr + *m;
  } else {
    yptr = dotptr + 1;
  }
  xcuptr = yptr + *m;

  iwarn = icntl[2];
  ihist = icntl[3];
  typeprec = icntl[4];
  iorthog = icntl[5];
  initguess = icntl[6];
  itermax = icntl[7];

  if (retlbl == 0) {
    comprsd = icntl[8];
  }

  if (retlbl != 0) {
    if (retlbl == 5) {
      goto L5;
    } else if (retlbl == 6) {
      goto L6;
    } else if (retlbl == 8) {
      goto L8;
    } else if (retlbl == 11) {
      goto L11;
    } else if (retlbl == 16) {
      goto L16;
    } else if (retlbl == 18) {
      goto L18;
    } else if (retlbl == 21) {
      goto L21;
    } else if (retlbl == 26) {
      goto L26;
    } else if (retlbl == 31) {
      goto L31;
    } else if (retlbl == 32) {
      goto L32;
    } else if (retlbl == 33) {
      goto L33;
    } else if (retlbl == 34) {
      goto L34;
    } else if (retlbl == 36) {
      goto L36;
    } else if (retlbl == 37) {
      goto L37;
    } else if (retlbl == 38) {
      goto L38;
    } else if (retlbl == 41) {
      goto L41;
    } else if (retlbl == 43) {
      goto L43;
    } else if (retlbl == 46) {
      goto L46;
    } else if (retlbl == 48) {
      goto L48;
    } else if (retlbl == 51) {
      goto L51;
    } else if (retlbl == 52) {
      goto L52;
    } else if (retlbl == 61) {
      goto L61;
    } else if (retlbl == 66) {
      goto L66;
    } else if (retlbl == 68) {
      goto L68;
    }
  }

  
  /* intialization of various variables */

  iterout = 0;
  beta = 0.;

  if (initguess == 0) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      x[j] = 0.;
    }
  }

  /*        bn = dnrm2(n,b,1) */

  irc[1] = 4;
  irc[2] = bptr;
  irc[3] = bptr;
  irc[4] = dotptr;
  irc[5] = 1;
  retlbl = 5;
  return 0;
 L5:
  bn = sqrt(dot[1]);

  if (bn == 0.) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      x[j] = 0.;
    }
    if (iwarn != 0) {
      printf(" WARNING GMRES : \n");
      printf(" Null right hand side : \n");
      printf(" solution set to zero : \n");
    }
    jh = 0;
    bea = 0.;
    be = 0.;
    printf(" %d %.3e \n",jh,bea);
    printf(" %.3e \n",be);
    info[1] = 0;
    info[2] = 0;
    rinfo[1] = 0.;
    rinfo[2] = 0.;
    irc[1] = 0;
    retlbl = 0;
    return 0;
  }

  /* Compute the scaling factor for the backward error on the */
  /*  unpreconditioned sytem */

  sa = cntl[2];
  sb = cntl[3];
  if (sa == 0. && sb == 0.) {
    sb = bn;
  }
  /* Compute the scaling factor for the backward error on the */
  /*  preconditioned sytem */

  spa = cntl[4];
  spb = cntl[5];
  if (spa == 0. && spb == 0.) {
    if (typeprec == 0 || typeprec == 2) {
      spb = bn;
    } else {
      irc[1] = 2;
      irc[2] = bptr;
      irc[4] = r0ptr;
      retlbl = 6;
      return 0;
    }
  }
 L6:
  if (spa == 0. && spb == 0.) {
    if (typeprec == 3 || typeprec == 1) {
      irc[1] = 4;
      irc[2] = r0ptr;
      irc[3] = r0ptr;
      irc[4] = dotptr;
      irc[5] = 1;
      retlbl = 8;
      return 0;
    }
  }
 L8:
  if (spa == 0. && spb == 0.) {
    if (typeprec == 3 || typeprec == 1) {
      spb = sqrt(dot[1]);

    }
  }
    
  /* Compute the first residual */
  /*           Y = AX : r0 <-- A x */

  /* The residual is computed only if the initial guess is not zero */

  if (initguess != 0) {
    irc[1] = 1;
    irc[2] = xptr;
    irc[4] = r0ptr;
    retlbl = 11;
    return 0;
  }
 L11:
  if (initguess != 0) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      r0[j] = b[j] - r0[j];
    }
  } else {
    gcopy(n, &b[1], &c__1, &r0[1], &c__1);
  }

  /* Compute the preconditioned residual if necessary */
  /*      M_1Y = X : w <-- M_1^{-1} r0 */

  if (typeprec == 0 || typeprec == 2) {
    gcopy(n, &r0[1], &c__1, &w[1], &c__1);
  } else {
    irc[1] = 2;
    irc[2] = r0ptr;
    irc[4] = wptr;
    retlbl = 16;
    return 0;
  }
 L16:

  /*       beta = dnrm2(n,w,1) */
  irc[1] = 4;
  irc[2] = wptr;
  irc[3] = wptr;
  irc[4] = dotptr;
  irc[5] = 1;
  retlbl = 18;
  return 0;
 L18:
  beta = sqrt(dot[1]);

  if (beta == 0.) {
    /*  The residual is exactly zero : x is the exact solution */
    info[1] = 0;
    info[2] = 0;
    rinfo[1] = 0.;
    rinfo[2] = 0.;
    irc[1] = 0;
    retlbl = 0;
    jh = 0;
    bea = 0.;
    be = 0.;
    printf(" %d %.3e \n",jh,bea);
    printf(" %.3e \n",be);
    if (iwarn != 0) {
      printf(" WARNING GMRES : \n");
      printf(" Initial residual is zero \n");
      printf(" initial guess is solution\n");
    }
    return 0;
  }

  aux = 1. / beta;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    v[j + v_dim1] = 0.;
  }
  gconstantmut(n, &aux, &w[1], &c__1, &v[v_dim1 + 1], &c__1);

  /*       Most outer loop : dgmres iteration */
  /*       REPEAT */
 L7:

  h__[(*m + 1) * h_dim1 + 1] = beta;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    h__[j + 1 + (*m + 1) * h_dim1] = 0.;
  }
  /*        Construction of the hessenberg matrix WORK and of the orthogonal */
  /*        basis V such that AV=VH */
  jh = 1;
 L10:
  /* Remark : this  do loop has been written with a while do */
  /*          because the */
  /*               " do jH=1,restart " */
  /*         fails with the reverse communication. */
  /*      do  jH=1,restart */

  /* Compute the preconditioned residual if necessary */

  if (typeprec == 2 || typeprec == 3) {
    /*           Y = M_2^{-1}X : w <-- M_2^{-1} V(1,jH) */
    irc[1] = 3;
    irc[2] = vptr + (jh - 1) * *n;
    irc[4] = wptr;
    retlbl = 21;
    return 0;
  } else {
    gcopy(n, &v[jh * v_dim1 + 1], &c__1, &w[1], &c__1);
  }
 L21:
  /*           Y = AX : r0 <-- A w */
  irc[1] = 1;
  irc[2] = wptr;
  irc[4] = r0ptr;
  retlbl = 26;
  return 0;
 L26:
  /*      MY = X : w <-- M_1^{-1} r0 */

  if (typeprec == 0 || typeprec == 2) {
    gcopy(n, &r0[1], &c__1, &w[1], &c__1);
  } else {
    irc[1] = 2;
    irc[2] = r0ptr;
    irc[4] = wptr;
    retlbl = 31;
    return 0;
  }
 L31:

  /* Orthogonalization using either MGS or IMGS */
  /* initialize the Hessenberg matrix to zero in order to be able to use */
  /*     IMGS as orthogonalization procedure. */
  i__1 = jh;
  for (j = 1; j <= i__1; ++j) {
    h__[j + jh * h_dim1] = 0.;
  }
  northo = 0;
 L19:
  ++northo;
  dloo = 0.;

  if (iorthog == 0 || iorthog == 1) {
    /* MGS */
    /*           do j=1,jH */
    j = 1;
    /*           REPEAT */
  }
 L23:
  if (iorthog == 0 || iorthog == 1) {

    /*             dVi     = ddot(n,V(1,j),1,w,1) */

    irc[1] = 4;
    irc[2] = vptr + (j - 1) * *n;
    irc[3] = wptr;
    irc[4] = dotptr;
    irc[5] = 1;
    retlbl = 32;
    return 0;
  }
 L32:
  if (iorthog == 0 || iorthog == 1) {
    dvi = dot[1];
    h__[j + jh * h_dim1] += dvi;
    /* Computing 2nd power */
    d__1 = abs(dvi);
    dloo += d__1 * d__1;
    aux = dvi * -1.;
    gconstantmut(n, &aux, &v[j * v_dim1 + 1], &c__1, &w[1], &c__1);
    ++j;
    if (j <= jh) {
      goto L23;
    }
    /*          enddo_j */
  } else {
    /* CGS */
    /* Gathered dot product calculation */
    /*           call dgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1) */

    irc[1] = 4;
    irc[2] = vptr;
    irc[3] = wptr;
    irc[4] = dotptr;
    irc[5] = jh;
    retlbl = 34;
    return 0;
  }
 L34:
  if (iorthog == 2 || iorthog == 3) {

    gconstantmut(&jh, &c_b276, &dot[1], &c__1, &h__[jh * h_dim1 + 1], &c__1);
    MatVectorProduct("N", n, &jh, &c_b280, &v[v_dim1 + 1], n, &dot[1], &c__1, &
		     c_b276, &w[1], &c__1);
    /* Computing 2nd power */
    d__1 = gnorm2(&jh, &dot[1], &c__1);
    dloo = d__1 * d__1;
  }

  irc[1] = 4;
  irc[2] = wptr;
  irc[3] = wptr;
  irc[4] = dotptr;
  irc[5] = 1;
  retlbl = 33;
  return 0;
 L33:
  dnormw = sqrt(dot[1]);

  if (iorthog == 1 || iorthog == 3) {
    /* IMGS / CGS orthogonalisation */
    dloo = sqrt(dloo);
    /* check the orthogonalization quality */
    if (dnormw * 2.f <= dloo && northo < 3) {
      goto L19;
    }
  }

  h__[jh + 1 + jh * h_dim1] = dnormw;
  if (jh < *m || icntl[8] == 0) {
    aux = 1. / dnormw;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      v[j + (jh + 1) * v_dim1] = 0.;
    }
    gconstantmut(n, &aux, &w[1], &c__1, &v[(jh + 1) * v_dim1 + 1], &c__1);
  }
  /* Apply previous Givens rotations to the new column of H */
  i__1 = jh - 1;
  for (j = 1; j <= i__1; ++j) {
    d__1 = rotcos[j];
    PlaneRotation(&c__1, &h__[j + jh * h_dim1], &c__1, &h__[j + 1 + jh * h_dim1], 
		  &c__1, &d__1, &rotsin[j]);
  }
  auxhjj = h__[jh + jh * h_dim1];
  auxhjp1j = h__[jh + 1 + jh * h_dim1];
  GivensRot(&auxhjj, &auxhjp1j, &temp, &rotsin[jh]);
  rotcos[jh] = temp;
  /* Apply current rotation to the rhs of the least squares problem */
  d__1 = rotcos[jh];
  PlaneRotation(&c__1, &h__[jh + (*m + 1) * h_dim1], &c__1, &h__[jh + 1 + (*m + 1) *
								 h_dim1], &c__1, &d__1, &rotsin[jh]);

  /* zabs(H(jH+1,m+1)) is the residual computed using the least squares */
  /*          solver */
  /* Complete the QR factorisation of the Hessenberg matrix by apply the current */
  /* rotation to the last entry of the collumn */
  d__1 = rotcos[jh];
  PlaneRotation(&c__1, &h__[jh + jh * h_dim1], &c__1, &h__[jh + 1 + jh * h_dim1], &
		c__1, &d__1, &rotsin[jh]);
  h__[jh + 1 + jh * h_dim1] = 0.;

  /* Get the Least square residual */

  dnormres = (d__1 = h__[jh + 1 + (*m + 1) * h_dim1], abs(d__1));
  if (spa != 0.) {

    /* Compute the solution of the current linear least squares problem */

    gcopy(&jh, &h__[(*m + 1) * h_dim1 + 1], &c__1, &ycurrent[1], &c__1);
    i__1 = *m + 1;
    UpperLowerSolver("U", "N", "N", &jh, &h__[h_offset], &i__1, &ycurrent[1], &c__1);

    /* Compute the value of the new iterate */

    MatVectorProduct("N", n, &jh, &c_b276, &v[v_offset], n, &ycurrent[1], &c__1, &
		     c_b305, &xcurrent[1], &c__1);

    if (typeprec == 2 || typeprec == 3) {
      /*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent */
      irc[1] = 3;
      irc[2] = xcuptr;
      irc[4] = r0ptr;
      retlbl = 36;
      return 0;
    } else {
      gcopy(n, &xcurrent[1], &c__1, &r0[1], &c__1);
    }
  }
 L36:

  if (spa != 0.) {
    /* Update the current solution */

    gcopy(n, &x[1], &c__1, &xcurrent[1], &c__1);
    gconstantmut(n, &c_b276, &r0[1], &c__1, &xcurrent[1], &c__1);
    irc[1] = 4;
    irc[2] = xcuptr;
    irc[3] = xcuptr;
    irc[4] = dotptr;
    irc[5] = 1;
    retlbl = 38;
    return 0;
  } else {
    dnormx = 1.;
  }
 L38:
  if (spa != 0.) {
    dnormx = sqrt(dot[1]);
  }

  bea = dnormres / (spa * dnormx + spb);

  /* Check the convergence based on the Arnoldi Backward error for the */
  /* preconditioned system */
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {

    /* The Arnoldi Backward error indicates that dgmres might have converge */
    /* enforce the calculation of the true residual at next restart */
    comprsd = 1;

    /*  If the update of X has not yet been performed */
    if (spa == 0.) {
      /* Compute the solution of the current linear least squares problem */

      gcopy(&jh, &h__[(*m + 1) * h_dim1 + 1], &c__1, &ycurrent[1], &
	    c__1);	  
      i__1 = *m + 1;
      UpperLowerSolver("U", "N", "N", &jh, &h__[h_offset], &i__1, &ycurrent[1],&c__1);

      /* Compute the value of the new iterate */

      MatVectorProduct("N", n, &jh, &c_b276, &v[v_offset], n, &ycurrent[1], &c__1,
		       &c_b305, &xcurrent[1], &c__1);

      if (typeprec == 2 || typeprec == 3) {

	/*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent */
	irc[1] = 3;
	irc[2] = xcuptr;
	irc[4] = r0ptr;
	retlbl = 37;
	return 0;
      } else {
	gcopy(n, &xcurrent[1], &c__1, &r0[1], &c__1);
      }
    }
  }
 L37:
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {
    if (spa == 0.) {
      /* Update the current solution */
      gcopy(n, &x[1], &c__1, &xcurrent[1], &c__1);
      gconstantmut(n, &c_b276, &r0[1], &c__1, &xcurrent[1], &c__1);
    }
    gcopy(n, &xcurrent[1], &c__1, &r0[1], &c__1);
    /* Compute the true residual, the Arnoldi one may be unaccurate */

    /*           Y = AX : w  <-- A r0 */

    irc[1] = 1;
    irc[2] = r0ptr;
    irc[4] = wptr;
    retlbl = 41;
    return 0;
  }
 L41:
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      w[j] = b[j] - w[j];
    }
    /* Compute the norm of the unpreconditioned residual */
    irc[1] = 4;
    irc[2] = wptr;
    irc[3] = wptr;
    irc[4] = dotptr;
    irc[5] = 1;
    retlbl = 43;
    return 0;
  }
 L43:
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {
    truenormres = sqrt(dot[1]);

    if (typeprec == 1 || typeprec == 3) {
      /*      MY = X : r0 <-- M_1^{-1} w */
      irc[1] = 2;
      irc[2] = wptr;
      irc[4] = r0ptr;
      retlbl = 46;
      return 0;
    } else {
      gcopy(n, &w[1], &c__1, &r0[1], &c__1);
    }
  }
 L46:
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {

    irc[1] = 4;
    irc[2] = r0ptr;
    irc[3] = r0ptr;
    irc[4] = dotptr;
    irc[5] = 1;
    retlbl = 48;
    return 0;
  }
 L48:
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {
    dnormres = sqrt(dot[1]);

    be = dnormres / (spa * dnormx + spb);
    /* Save the backward error on a file if convergence history requested */
    if (ihist != 0) {
      i__1 = iterout * *m + jh;
      printf("   %d      %.3e       %.3e \n",i__1,bea,be);
      /* io___168.ciunit = ihist;
	 s_wsfe(&io___168);
	    
	 do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(int));
	 do_fio(&c__1, (char *)&bea, (ftnlen)sizeof(double));
	 do_fio(&c__1, (char *)&be, (ftnlen)sizeof(double));
	 e_wsfe();*/
    }

  }


  /* Check again the convergence */
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {
    if (be <= cntl[1] || iterout * *m + jh >= itermax) {
      /* The convergence has been achieved, we restore the solution in x */
      /* and compute the two backward errors. */
      gcopy(n, &xcurrent[1], &c__1, &x[1], &c__1);

      if (sa != 0.) {

	irc[1] = 4;
	irc[2] = xptr;
	irc[3] = xptr;
	irc[4] = dotptr;
	irc[5] = 1;
	retlbl = 51;
	return 0;
      }
    }
  }
 L51:
  if (bea <= cntl[1] || iterout * *m + jh >= itermax) {
    if (be <= cntl[1] || iterout * *m + jh >= itermax) {
      if (sa != 0.) {
	dnormx = sqrt(dot[1]);
      } else {
	dnormx = 1.;
      }
      /* Return the backward errors */
      rinfo[1] = be;
      rinfo[2] = truenormres / (sa * dnormx + sb);
      if (be <= cntl[1]) {
	info[1] = 0;
	if (ihist != 0) {
	  printf(" Convergence achieved \n");
	}
      } else if (be > cntl[1]) {
	if (iwarn != 0) {
	  printf(" WARNING GMRES \n");
	  i__1 = iterout * *m + jh;
	  printf(" No convergence after : %d iterations\n",i__1);

	}
	if (ihist != 0) {
	  printf(" WARNING GMRES \n");
	  i__1 = iterout * *m + jh;
	  printf(" No convergence after : %d iterations\n",i__1);
	}
	info[1] = -4;
      }
      if (ihist != 0) {
	printf("Residu on the preconditioned system: %.5e \n",rinfo[1]);
	printf("Residu on the unpreconditioned system: %.5e \n",rinfo[2]);
      }
      info[2] = iterout * *m + jh;
      if (ihist != 0) {
	printf("info(1) =  %d\n",info[1]);
	printf("Number of iterations (info(2)) = %d  \n",info[2]);
      }
      irc[1] = 0;
      retlbl = 0;
      return 0;
    }
  } else {
    /* Save the backward error on a file if convergence history requested */
    if (ihist != 0) {	 
      i__1 = iterout * *m + jh;
      printf("   %d      %.3e      --  \n",i__1,bea);
    }
  }

  ++jh;
  if (jh <= *m) {
    goto L10;
  }

  ++iterout;

  /* we have completed the Krylov space construction, we restart if */
  /* we have not yet exceeded the maximum number of iterations allowed. */

  if (spa == 0. && bea > cntl[1]) {

    /* Compute the solution of the current linear least squares problem */

    --jh;
    gcopy(&jh, &h__[(*m + 1) * h_dim1 + 1], &c__1, &ycurrent[1], &c__1);
    i__1 = *m + 1;
    UpperLowerSolver("U", "N", "N", &jh, &h__[h_offset], &i__1, &ycurrent[1], &c__1);

    /* Compute the value of the new iterate */

    MatVectorProduct("N", n, &jh, &c_b276, &v[v_offset], n, &ycurrent[1], &c__1, &
		     c_b305, &xcurrent[1], &c__1);

    if (typeprec == 2 || typeprec == 3) {

      /*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent */
      irc[1] = 3;
      irc[2] = xcuptr;
      irc[4] = r0ptr;
      retlbl = 52;
      return 0;
    } else {
      gcopy(n, &xcurrent[1], &c__1, &r0[1], &c__1);
    }
  }
 L52:
  if (spa == 0. && bea > cntl[1]) {
    /* Update the current solution */

    gcopy(n, &x[1], &c__1, &xcurrent[1], &c__1);
    gconstantmut(n, &c_b276, &r0[1], &c__1, &xcurrent[1], &c__1);
  }

  gcopy(n, &xcurrent[1], &c__1, &x[1], &c__1);

  if (comprsd == 1) {
    /* Compute the true residual */

    gcopy(n, &x[1], &c__1, &w[1], &c__1);
    irc[1] = 1;
    irc[2] = wptr;
    irc[4] = r0ptr;
    retlbl = 61;
    return 0;
  }
 L61:
  if (comprsd == 1) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      r0[j] = b[j] - r0[j];
    }

    /* Precondition the new residual if necessary */

    if (typeprec == 1 || typeprec == 3) {

      /*      MY = X : w <-- M_1^{-1} r0 */
      irc[1] = 2;
      irc[2] = r0ptr;
      irc[4] = wptr;
      retlbl = 66;
      return 0;
    } else {
      gcopy(n, &r0[1], &c__1, &w[1], &c__1);
    }
  }
 L66:

  if (comprsd == 1) {
    irc[1] = 4;
    irc[2] = wptr;
    irc[3] = wptr;
    irc[4] = dotptr;
    irc[5] = 1;
    retlbl = 68;
    return 0;
  }
 L68:
  if (comprsd == 1) {
    beta = sqrt(dot[1]);

  } else {
    /* Use recurrence to approximate the residual at restart */
    beta = (d__1 = h__[*m + 1 + (*m + 1) * h_dim1], abs(d__1));
    /* Apply the Givens rotation is the reverse order */
    for (j = *m; j >= 1; --j) {
      h__[j + (*m + 1) * h_dim1] = 0.;
      d__1 = rotcos[j];
      d__2 = -rotsin[j];
      PlaneRotation(&c__1, &h__[j + (*m + 1) * h_dim1], &c__1, &h__[j + 1 + (*m 
									     + 1) * h_dim1], &c__1, &d__1, &d__2);
    }

    /* On applique les vecteurs V */
    i__1 = *m + 1;
    MatVectorProduct("N", n, &i__1, &c_b276, &v[v_offset], n, &h__[(*m + 1) * 
								   h_dim1 + 1], &c__1, &c_b305, &w[1], &c__1);
  }
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    v[j + v_dim1] = 0.;
  }
  aux = 1. / beta;
  gconstantmut(n, &aux, &w[1], &c__1, &v[v_dim1 + 1], &c__1);

  goto L7;

} /* dgmres_ */



/* Subroutine */ int init_dgmres(int *icntl, double *cntl)
{

  /*  Purpose */
  /*  ======= */
  /*    Set default values for the parameters defining the characteristics */
  /* of the Gmres algorithm. */
  /*  See the User's Guide for an example of use. */


  /* Written : April 1997 */
  /* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton */
  /*             Parallel Algorithms - CERFACS */


  /*  Arguments */
  /*  ========= */

  /* icntl    (input) INT array. length 6 */
  /*            icntl(1) : stdout for error messages */
  /*            icntl(2) : stdout for warnings */
  /*            icntl(3) : stdout for convergence history */
  /*            icntl(4) : 0 - no preconditioning */
  /*                       1 - left preconditioning */
  /*                       2 - right preconditioning */
  /*                       3 - double side preconditioning */
  /*                       4 - error, default set in Init */
  /*            icntl(5) : 0 - modified Gram-Schmidt */
  /*                       1 - iterative modified Gram-Schmidt */
  /*                       2 - classical Gram-Schmidt */
  /*                       3 - iterative classical Gram-Schmidt */
  /*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set) */
  /*                       1 - user supplied initial guess */
  /*            icntl(7) : maximum number of iterations */
  /*            icntl(8) : 1 - default compute the true residual at each restart */
  /*                       0 - use recurence formaula at restart */

  /* cntl     (input) double precision array, length 5 */
  /*            cntl(1) : tolerance for convergence */
  /*            cntl(2) : scaling factor for normwise perturbation on A */
  /*            cntl(3) : scaling factor for normwise perturbation on b */
  /*            cntl(4) : scaling factor for normwise perturbation on the */
  /*                      preconditioned matrix */
  /*            cntl(5) : scaling factor for normwise perturbation on */
  /*                      preconditioned right hand side */

  /* Output variables */
  /* ---------------- */

  /* Parameter adjustments */
  --cntl;
  --icntl;

  /* Function Body */
  icntl[1] = 6;
  icntl[2] = 6;
  icntl[3] = 0;
  icntl[4] = 4;
  icntl[5] = 0;
  icntl[6] = 0;
  icntl[7] = -1;
  icntl[8] = 1;

  cntl[1] = 1e-5;
  cntl[2] = 0.;
  cntl[3] = 0.;
  cntl[4] = 0.;
  cntl[5] = 0.;

  return 0;   
} /* init_dgmres__ */
