/* zabs.f -- translated by f2c (version 20190311).
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

doublereal azabs_(doublereal *zr, doublereal *zi)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal q, s, u, v;

/* ***BEGIN PROLOGUE  AZABS */
/* ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY */

/*     AZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE */
/*     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI) */

/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  AZABS */
    u = abs(*zr);
    v = abs(*zi);
    s = u + v;
/* ----------------------------------------------------------------------- */
/*     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A */
/*     TRUE FLOATING ZERO */
/* ----------------------------------------------------------------------- */
    s *= 1.;
    if (s == 0.) {
	goto L20;
    }
    if (u > v) {
	goto L10;
    }
    q = u / v;
    ret_val = v * sqrt(q * q + 1.);
    return ret_val;
L10:
    q = v / u;
    ret_val = u * sqrt(q * q + 1.);
    return ret_val;
L20:
    ret_val = 0.;
    return ret_val;
} /* azabs_ */

