/* dcosqf.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dcosqf_(integer *n, doublereal *x, doublereal *wsave)
{
    /* Initialized data */

    static doublereal sqrt2 = 1.4142135623730950488;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal tsqx;
    extern /* Subroutine */ int cosqf1_(integer *, doublereal *, doublereal *,
	     doublereal *);

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L102;
    } else if (i__1 == 0) {
	goto L101;
    } else {
	goto L103;
    }
L101:
    tsqx = sqrt2 * x[2];
    x[2] = x[1] - tsqx;
    x[1] += tsqx;
L102:
    return 0;
L103:
    cosqf1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* dcosqf_ */

