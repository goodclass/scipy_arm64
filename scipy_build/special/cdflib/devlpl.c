/* devlpl.f -- translated by f2c (version 20190311).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

doublereal devlpl_(doublereal *a, integer *n, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i__;
    static doublereal term;

/* ********************************************************************** */

/*     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X) */
/*              Double precision EVALuate a PoLynomial at X */


/*                              Function */


/*     returns */
/*          A(1) + A(2)*X + ... + A(N)*X**(N-1) */


/*                              Arguments */


/*     A --> Array of coefficients of the polynomial. */
/*                                        A is DOUBLE PRECISION(N) */

/*     N --> Length of A, also degree of polynomial - 1. */
/*                                        N is INTEGER */

/*     X --> Point at which the polynomial is to be evaluated. */
/*                                        X is DOUBLE PRECISION */

/* ********************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    --a;

    /* Function Body */
    term = a[*n];
    for (i__ = *n - 1; i__ >= 1; --i__) {
	term = a[i__] + term * *x;
/* L10: */
    }
    ret_val = term;
    return ret_val;
} /* devlpl_ */

