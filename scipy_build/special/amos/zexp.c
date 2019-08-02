/* zexp.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int azexp_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi)
{
    /* Builtin functions */
    double exp(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal ca, cb, zm;

/* ***BEGIN PROLOGUE  AZEXP */
/* ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY */

/*     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A) */

/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  AZEXP */
    zm = exp(*ar);
    ca = zm * cos(*ai);
    cb = zm * sin(*ai);
    *br = ca;
    *bi = cb;
    return 0;
} /* azexp_ */

