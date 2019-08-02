/* dcosti.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dcosti_(integer *n, doublereal *wsave)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979323846;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer k, kc;
    static doublereal fk, dt;
    static integer nm1, np1, ns2;
    extern /* Subroutine */ int dffti_(integer *, doublereal *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 3) {
	return 0;
    }
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    dt = pi / (real) nm1;
    fk = 0.;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	fk += 1.;
	wsave[k] = sin(fk * dt) * 2.;
	wsave[kc] = cos(fk * dt) * 2.;
/* L101: */
    }
    dffti_(&nm1, &wsave[*n + 1]);
    return 0;
} /* dcosti_ */

