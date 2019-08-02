/* zshch.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zshch_(doublereal *zr, doublereal *zi, doublereal *cshr, 
	doublereal *cshi, doublereal *cchr, doublereal *cchi)
{
    /* Builtin functions */
    double sinh(doublereal), cosh(doublereal), sin(doublereal), cos(
	    doublereal);

    /* Local variables */
    static doublereal ch, cn, sh, sn;

/* ***BEGIN PROLOGUE  ZSHCH */
/* ***REFER TO  ZBESK,ZBESH */

/*     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y) */
/*     AND CCH=COSH(X+I*Y), WHERE I**2=-1. */

/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  ZSHCH */

    sh = sinh(*zr);
    ch = cosh(*zr);
    sn = sin(*zi);
    cn = cos(*zi);
    *cshr = sh * cn;
    *cshi = ch * sn;
    *cchr = ch * cn;
    *cchi = sh * sn;
    return 0;
} /* zshch_ */

