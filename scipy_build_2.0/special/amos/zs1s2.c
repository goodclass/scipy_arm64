/* zs1s2.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zs1s2_(doublereal *zrr, doublereal *zri, doublereal *s1r,
	 doublereal *s1i, doublereal *s2r, doublereal *s2i, integer *nz, 
	doublereal *ascle, doublereal *alim, integer *iuf)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal aa, c1i, as1, as2, c1r, aln, s1di, s1dr;
    static integer idum;
    extern doublereal azabs_(doublereal *, doublereal *);
    extern /* Subroutine */ int azlog_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *), azexp_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZS1S2 */
/* ***REFER TO  ZBESK,ZAIRY */

/*     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE */
/*     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON- */
/*     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION. */
/*     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF */
/*     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER */
/*     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE */
/*     PRECISION ABOVE THE UNDERFLOW LIMIT. */

/* ***ROUTINES CALLED  AZABS,AZEXP,AZLOG */
/* ***END PROLOGUE  ZS1S2 */
/*     COMPLEX CZERO,C1,S1,S1D,S2,ZR */
    *nz = 0;
    as1 = azabs_(s1r, s1i);
    as2 = azabs_(s2r, s2i);
    if (*s1r == 0. && *s1i == 0.) {
	goto L10;
    }
    if (as1 == 0.) {
	goto L10;
    }
    aln = -(*zrr) - *zrr + log(as1);
    s1dr = *s1r;
    s1di = *s1i;
    *s1r = zeror;
    *s1i = zeroi;
    as1 = zeror;
    if (aln < -(*alim)) {
	goto L10;
    }
    azlog_(&s1dr, &s1di, &c1r, &c1i, &idum);
    c1r = c1r - *zrr - *zrr;
    c1i = c1i - *zri - *zri;
    azexp_(&c1r, &c1i, s1r, s1i);
    as1 = azabs_(s1r, s1i);
    ++(*iuf);
L10:
    aa = max(as1,as2);
    if (aa > *ascle) {
	return 0;
    }
    *s1r = zeror;
    *s1i = zeroi;
    *s2r = zeror;
    *s2i = zeroi;
    *nz = 1;
    *iuf = 0;
    return 0;
} /* zs1s2_ */

