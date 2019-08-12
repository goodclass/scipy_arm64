/* cosqi.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cosqi_(integer *n, real *wsave)
{
    /* Initialized data */

    static real pih = 1.57079632679491f;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer k;
    static real fk, dt;
    extern /* Subroutine */ int rffti_(integer *, real *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    dt = pih / (real) (*n);
    fk = 0.f;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	fk += 1.f;
	wsave[k] = cos(fk * dt);
/* L101: */
    }
    rffti_(n, &wsave[*n + 1]);
    return 0;
} /* cosqi_ */

