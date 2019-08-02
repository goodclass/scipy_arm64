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

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal x1;
    extern /* Subroutine */ int cosqb1_(integer *, doublereal *, doublereal *,
	     doublereal *);

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    x[1] *= 4.;
    return 0;
L102:
    x1 = (x[1] + x[2]) * 4.;
    x[2] = tsqrt2 * (x[1] - x[2]);
    x[1] = x1;
    return 0;
L103:
    cosqb1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* dcosqb_ */

