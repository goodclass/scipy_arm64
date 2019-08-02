/* fprppo.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fprppo_(integer *nu, integer *nv, integer *if1, integer *
	if2, doublereal *cosi, doublereal *ratio, doublereal *c__, doublereal 
	*f, integer *ncoff)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, ii, nu4, nvv, iopt;

/*  given the coefficients of a constrained bicubic spline, as determined */
/*  in subroutine fppola, subroutine fprppo calculates the coefficients */
/*  in the standard b-spline representation of bicubic splines. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    cosi -= 6;
    --f;
    --c__;

    /* Function Body */
    nu4 = *nu - 4;
    nvv = *nv - 7;
    iopt = *if1 + 1;
    i__1 = *ncoff;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f[i__] = 0.f;
/* L10: */
    }
    i__ = 0;
    i__1 = nu4;
    for (l = 1; l <= i__1; ++l) {
	ii = i__;
	if (l > iopt) {
	    goto L80;
	}
	switch (l) {
	    case 1:  goto L20;
	    case 2:  goto L40;
	    case 3:  goto L60;
	}
L20:
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    f[i__] = c__[1];
/* L30: */
	}
	j = 1;
	goto L100;
L40:
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    f[i__] = c__[1] + c__[2] * cosi[k * 5 + 1] + c__[3] * cosi[k * 5 
		    + 2];
/* L50: */
	}
	j = 3;
	goto L100;
L60:
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    f[i__] = c__[1] + *ratio * (c__[2] * cosi[k * 5 + 1] + c__[3] * 
		    cosi[k * 5 + 2]) + c__[4] * cosi[k * 5 + 3] + c__[5] * 
		    cosi[k * 5 + 4] + c__[6] * cosi[k * 5 + 5];
/* L70: */
	}
	j = 6;
	goto L100;
L80:
	if (l == nu4 && *if2 != 0) {
	    goto L120;
	}
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    ++j;
	    f[i__] = c__[j];
/* L90: */
	}
L100:
	for (k = 1; k <= 3; ++k) {
	    ++ii;
	    ++i__;
	    f[i__] = f[ii];
/* L110: */
	}
L120:
	;
    }
    i__1 = *ncoff;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = f[i__];
/* L130: */
    }
    return 0;
} /* fprppo_ */

