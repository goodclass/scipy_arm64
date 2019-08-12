/* dcosqb.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dcosqb_(integer *n, doublereal *x, doublereal *wsave)
{
    /* Initialized data */

    static doublereal tsqrt2 = 2.8284271247461900976;

    static doublereal x1;
    extern /* Subroutine */ int dcosqb1_(integer *, doublereal *, doublereal *
	    , doublereal *);

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if (*n < 2) {
	goto L101;
    }
    if (*n == 2) {
	goto L102;
    }
    goto L103;
L101:
    x[1] *= 4.;
    return 0;
L102:
    x1 = (x[1] + x[2]) * 4.;
    x[2] = tsqrt2 * (x[1] - x[2]);
    x[1] = x1;
    return 0;
L103:
    dcosqb1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* dcosqb_ */

/* Subroutine */ int dcosqb1_(integer *n, doublereal *x, doublereal *w, 
	doublereal *xh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, kc, np2, ns2;
    static doublereal xim1;
    static integer modn;
    extern /* Subroutine */ int dfftb_(integer *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --xh;
    --w;
    --x;

    /* Function Body */
    ns2 = (*n + 1) / 2;
    np2 = *n + 2;
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	xim1 = x[i__ - 1] + x[i__];
	x[i__] -= x[i__ - 1];
	x[i__ - 1] = xim1;
/* L101: */
    }
    x[1] += x[1];
    modn = *n % 2;
    if (modn == 0) {
	x[*n] += x[*n];
    }
    dfftb_(n, &x[1], &xh[1]);
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	xh[k] = w[k - 1] * x[kc] + w[kc - 1] * x[k];
	xh[kc] = w[k - 1] * x[k] - w[kc - 1] * x[kc];
/* L102: */
    }
    if (modn == 0) {
	x[ns2 + 1] = w[ns2] * (x[ns2 + 1] + x[ns2 + 1]);
    }
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	x[k] = xh[k] + xh[kc];
	x[kc] = xh[k] - xh[kc];
/* L103: */
    }
    x[1] += x[1];
    return 0;
} /* dcosqb1_ */

