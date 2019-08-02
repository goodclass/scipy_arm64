/* zlog.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int azlog_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, integer *ierr)
{
    /* Initialized data */

    static doublereal dpi = 3.141592653589793238462643383;
    static doublereal dhpi = 1.570796326794896619231321696;

    /* Builtin functions */
    double atan(doublereal), log(doublereal);

    /* Local variables */
    static doublereal zm;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal dtheta;

/* ***BEGIN PROLOGUE  AZLOG */
/* ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY */

/*     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A) */
/*     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0) */
/* ***ROUTINES CALLED  AZABS */
/* ***END PROLOGUE  AZLOG */

    *ierr = 0;
    if (*ar == 0.) {
	goto L10;
    }
    if (*ai == 0.) {
	goto L20;
    }
    dtheta = atan(*ai / *ar);
    if (dtheta <= 0.) {
	goto L40;
    }
    if (*ar < 0.) {
	dtheta -= dpi;
    }
    goto L50;
L10:
    if (*ai == 0.) {
	goto L60;
    }
    *bi = dhpi;
    *br = log((abs(*ai)));
    if (*ai < 0.) {
	*bi = -(*bi);
    }
    return 0;
L20:
    if (*ar > 0.) {
	goto L30;
    }
    *br = log((abs(*ar)));
    *bi = dpi;
    return 0;
L30:
    *br = log(*ar);
    *bi = 0.;
    return 0;
L40:
    if (*ar < 0.) {
	dtheta += dpi;
    }
L50:
    zm = azabs_(ar, ai);
    *br = log(zm);
    *bi = dtheta;
    return 0;
L60:
    *ierr = 1;
    return 0;
} /* azlog_ */

