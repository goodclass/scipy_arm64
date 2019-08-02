/* dsinti.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dsinti_(integer *n, doublereal *wsave)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979323846;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal dt;
    static integer np1, ns2;
    extern /* Subroutine */ int dffti_(integer *, doublereal *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 1) {
	return 0;
    }
    ns2 = *n / 2;
    np1 = *n + 1;
    dt = pi / (real) np1;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	wsave[k] = sin(k * dt) * 2.;
/* L101: */
    }
    dffti_(&np1, &wsave[ns2 + 1]);
    return 0;
} /* dsinti_ */

