/* dcosqi.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dcosqi_(integer *n, doublereal *wsave)
{
    /* Initialized data */

    static doublereal pih = 1.57079632679489661923;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer k;
    static doublereal fk, dt;
    extern /* Subroutine */ int dffti_(integer *, doublereal *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    dt = pih / (real) (*n);
    fk = 0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	fk += 1.;
	wsave[k] = cos(fk * dt);
/* L101: */
    }
    dffti_(n, &wsave[*n + 1]);
    return 0;
} /* dcosqi_ */

