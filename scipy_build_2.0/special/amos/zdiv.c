/* zdiv.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zdiv_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci)
{
    static doublereal ca, cb, cc, cd, bm;
    extern doublereal azabs_(doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZDIV */
/* ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY */

/*     DOUBLE PRECISION COMPLEX DIVIDE C=A/B. */

/* ***ROUTINES CALLED  AZABS */
/* ***END PROLOGUE  ZDIV */
    bm = 1. / azabs_(br, bi);
    cc = *br * bm;
    cd = *bi * bm;
    ca = (*ar * cc + *ai * cd) * bm;
    cb = (*ai * cc - *ar * cd) * bm;
    *cr = ca;
    *ci = cb;
    return 0;
} /* zdiv_ */

