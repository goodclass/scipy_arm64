/* zmlt.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zmlt_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci)
{
    static doublereal ca, cb;

/* ***BEGIN PROLOGUE  ZMLT */
/* ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY */

/*     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B. */

/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  ZMLT */
    ca = *ar * *br - *ai * *bi;
    cb = *ar * *bi + *ai * *br;
    *cr = ca;
    *ci = cb;
    return 0;
} /* zmlt_ */

