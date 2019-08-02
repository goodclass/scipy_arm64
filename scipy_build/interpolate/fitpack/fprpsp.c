/* fprpsp.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fprpsp_(integer *nt, integer *np, doublereal *co, 
	doublereal *si, doublereal *c__, doublereal *f, integer *ncoff)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal c1, c2, c3, cn;
    static integer ii, np4, nt4, npp, ncof;

/*  given the coefficients of a spherical spline function, subroutine */
/*  fprpsp calculates the coefficients in the standard b-spline re- */
/*  presentation of this bicubic spline. */
/*  .. */
/*  ..scalar arguments */
/*  ..array arguments */
/*  ..local scalars */
/*  .. */
    /* Parameter adjustments */
    --si;
    --co;
    --f;
    --c__;

    /* Function Body */
    nt4 = *nt - 4;
    np4 = *np - 4;
    npp = np4 - 3;
    ncof = npp * (nt4 - 4) + 6;
    c1 = c__[1];
    cn = c__[ncof];
    j = *ncoff;
    i__1 = np4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__] = c1;
	f[j] = cn;
	--j;
/* L10: */
    }
    i__ = np4;
    j = 1;
    i__1 = nt4;
    for (l = 3; l <= i__1; ++l) {
	ii = i__;
	if (l == 3 || l == nt4) {
	    goto L30;
	}
	i__2 = npp;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    ++j;
	    f[i__] = c__[j];
/* L20: */
	}
	goto L50;
L30:
	if (l == nt4) {
	    c1 = cn;
	}
	c2 = c__[j + 1];
	c3 = c__[j + 2];
	j += 2;
	i__2 = npp;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    f[i__] = c1 + c2 * co[k] + c3 * si[k];
/* L40: */
	}
L50:
	for (k = 1; k <= 3; ++k) {
	    ++ii;
	    ++i__;
	    f[i__] = f[ii];
/* L60: */
	}
/* L70: */
    }
    i__1 = *ncoff;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = f[i__];
/* L80: */
    }
    return 0;
} /* fprpsp_ */

