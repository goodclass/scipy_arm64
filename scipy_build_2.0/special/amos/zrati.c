/* zrati.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zrati_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *n, doublereal *cyr, doublereal *cyi, doublereal *tol)
{
    /* Initialized data */

    static doublereal czeror = 0.;
    static doublereal czeroi = 0.;
    static doublereal coner = 1.;
    static doublereal conei = 0.;
    static doublereal rt2 = 1.41421356237309505;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal ak;
    static integer id, kk;
    static doublereal az, ap1, ap2, p1i, p2i, t1i, p1r, p2r, t1r, arg, rak, 
	    rho;
    static integer inu;
    static doublereal pti, tti, rzi, ptr, ttr, rzr, rap1, flam, dfnu, fdnu;
    static integer magz, idnu;
    static doublereal fnup;
    extern /* Subroutine */ int zdiv_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static doublereal test, test1, amagz;
    extern doublereal azabs_(doublereal *, doublereal *);
    static integer itime;
    static doublereal cdfnui, cdfnur;

/* ***BEGIN PROLOGUE  ZRATI */
/* ***REFER TO  ZBESI,ZBESK,ZBESH */

/*     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD */
/*     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD */
/*     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B, */
/*     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973, */
/*     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER, */
/*     BY D. J. SOOKNE. */

/* ***ROUTINES CALLED  AZABS,ZDIV */
/* ***END PROLOGUE  ZRATI */
/*     COMPLEX Z,CY(1),CONE,CZERO,P1,P2,T1,RZ,PT,CDFNU */
    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */
    az = azabs_(zr, zi);
    inu = (integer) ((real) (*fnu));
    idnu = inu + *n - 1;
    magz = (integer) ((real) az);
    amagz = (doublereal) ((real) (magz + 1));
    fdnu = (doublereal) ((real) idnu);
    fnup = max(amagz,fdnu);
    id = idnu - magz - 1;
    itime = 1;
    k = 1;
    ptr = 1. / az;
    rzr = ptr * (*zr + *zr) * ptr;
    rzi = -ptr * (*zi + *zi) * ptr;
    t1r = rzr * fnup;
    t1i = rzi * fnup;
    p2r = -t1r;
    p2i = -t1i;
    p1r = coner;
    p1i = conei;
    t1r += rzr;
    t1i += rzi;
    if (id > 0) {
	id = 0;
    }
    ap2 = azabs_(&p2r, &p2i);
    ap1 = azabs_(&p1r, &p1i);
/* ----------------------------------------------------------------------- */
/*     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU */
/*     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT */
/*     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR */
/*     PREMATURELY. */
/* ----------------------------------------------------------------------- */
    arg = (ap2 + ap2) / (ap1 * *tol);
    test1 = sqrt(arg);
    test = test1;
    rap1 = 1. / ap1;
    p1r *= rap1;
    p1i *= rap1;
    p2r *= rap1;
    p2i *= rap1;
    ap2 *= rap1;
L10:
    ++k;
    ap1 = ap2;
    ptr = p2r;
    pti = p2i;
    p2r = p1r - (t1r * ptr - t1i * pti);
    p2i = p1i - (t1r * pti + t1i * ptr);
    p1r = ptr;
    p1i = pti;
    t1r += rzr;
    t1i += rzi;
    ap2 = azabs_(&p2r, &p2i);
    if (ap1 <= test) {
	goto L10;
    }
    if (itime == 2) {
	goto L20;
    }
    ak = azabs_(&t1r, &t1i) * .5;
    flam = ak + sqrt(ak * ak - 1.);
/* Computing MIN */
    d__1 = ap2 / ap1;
    rho = min(d__1,flam);
    test = test1 * sqrt(rho / (rho * rho - 1.));
    itime = 2;
    goto L10;
L20:
    kk = k + 1 - id;
    ak = (doublereal) ((real) kk);
    t1r = ak;
    t1i = czeroi;
    dfnu = *fnu + (doublereal) ((real) (*n - 1));
    p1r = 1. / ap2;
    p1i = czeroi;
    p2r = czeror;
    p2i = czeroi;
    i__1 = kk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ptr = p1r;
	pti = p1i;
	rap1 = dfnu + t1r;
	ttr = rzr * rap1;
	tti = rzi * rap1;
	p1r = ptr * ttr - pti * tti + p2r;
	p1i = ptr * tti + pti * ttr + p2i;
	p2r = ptr;
	p2i = pti;
	t1r -= coner;
/* L30: */
    }
    if (p1r != czeror || p1i != czeroi) {
	goto L40;
    }
    p1r = *tol;
    p1i = *tol;
L40:
    zdiv_(&p2r, &p2i, &p1r, &p1i, &cyr[*n], &cyi[*n]);
    if (*n == 1) {
	return 0;
    }
    k = *n - 1;
    ak = (doublereal) ((real) k);
    t1r = ak;
    t1i = czeroi;
    cdfnur = *fnu * rzr;
    cdfnui = *fnu * rzi;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[k + 1];
	pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[k + 1];
	ak = azabs_(&ptr, &pti);
	if (ak != czeror) {
	    goto L50;
	}
	ptr = *tol;
	pti = *tol;
	ak = *tol * rt2;
L50:
	rak = coner / ak;
	cyr[k] = rak * ptr * rak;
	cyi[k] = -rak * pti * rak;
	t1r -= coner;
	--k;
/* L60: */
    }
    return 0;
} /* zrati_ */

