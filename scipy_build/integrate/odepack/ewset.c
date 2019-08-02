/* ewset.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int ewset_(integer *n, integer *itol, doublereal *rtol, 
	doublereal *atol, doublereal *ycur, doublereal *ewt)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;

/* lll. optimize */
/* ----------------------------------------------------------------------- */
/* this subroutine sets the error weight vector ewt according to */
/*     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n, */
/* with the subscript on rtol and/or atol possibly replaced by 1 above, */
/* depending on the value of itol. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ewt;
    --ycur;
    --rtol;
    --atol;

    /* Function Body */
    switch (*itol) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }
L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	ewt[i__] = rtol[1] * (d__1 = ycur[i__], abs(d__1)) + atol[1];
    }
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	ewt[i__] = rtol[1] * (d__1 = ycur[i__], abs(d__1)) + atol[i__];
    }
    return 0;
L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L35: */
	ewt[i__] = rtol[i__] * (d__1 = ycur[i__], abs(d__1)) + atol[1];
    }
    return 0;
L40:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L45: */
	ewt[i__] = rtol[i__] * (d__1 = ycur[i__], abs(d__1)) + atol[i__];
    }
    return 0;
/* ----------------------- end of subroutine ewset ----------------------- */
} /* ewset_ */

