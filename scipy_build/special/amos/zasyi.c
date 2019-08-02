/* zasyi.f -- translated by f2c (version 20190311).
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

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int zasyi_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *yr, doublereal *yi, integer *
	nz, doublereal *rl, doublereal *tol, doublereal *elim, doublereal *
	alim)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal rtpi = .159154943091895336;
    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;
    static doublereal coner = 1.;
    static doublereal conei = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal s, aa, bb;
    static integer ib;
    static doublereal ak, bk;
    static integer il, jl;
    static doublereal az;
    static integer nn;
    static doublereal p1i, s2i, p1r, s2r, cki, dki, fdn, arg, aez, arm, ckr, 
	    dkr, czi, ezi, sgn;
    static integer inu;
    static doublereal raz, czr, ezr, sqk, sti, rzi, tzi, str, rzr, tzr, ak1i, 
	    ak1r, cs1i, cs2i, cs1r, cs2r, dnu2, rtr1, dfnu, atol;
    extern /* Subroutine */ int zdiv_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *), zmlt_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer koded;
    extern doublereal azabs_(doublereal *, doublereal *);
    extern /* Subroutine */ int azexp_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int azsqrt_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZASYI */
/* ***REFER TO  ZBESI,ZBESK */

/*     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY */
/*     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE */
/*     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN. */
/*     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1. */

/* ***ROUTINES CALLED  D1MACH,AZABS,ZDIV,AZEXP,ZMLT,AZSQRT */
/* ***END PROLOGUE  ZASYI */
/*     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    az = azabs_(zr, zi);
    arm = d1mach_(&c__1) * 1e3;
    rtr1 = sqrt(arm);
    il = min(2,*n);
    dfnu = *fnu + (doublereal) ((real) (*n - il));
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST */
/* ----------------------------------------------------------------------- */
    raz = 1. / az;
    str = *zr * raz;
    sti = -(*zi) * raz;
    ak1r = rtpi * str * raz;
    ak1i = rtpi * sti * raz;
    azsqrt_(&ak1r, &ak1i, &ak1r, &ak1i);
    czr = *zr;
    czi = *zi;
    if (*kode != 2) {
	goto L10;
    }
    czr = zeror;
    czi = *zi;
L10:
    if (abs(czr) > *elim) {
	goto L100;
    }
    dnu2 = dfnu + dfnu;
    koded = 1;
    if (abs(czr) > *alim && *n > 2) {
	goto L20;
    }
    koded = 0;
    azexp_(&czr, &czi, &str, &sti);
    zmlt_(&ak1r, &ak1i, &str, &sti, &ak1r, &ak1i);
L20:
    fdn = 0.;
    if (dnu2 > rtr1) {
	fdn = dnu2 * dnu2;
    }
    ezr = *zr * 8.;
    ezi = *zi * 8.;
/* ----------------------------------------------------------------------- */
/*     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE */
/*     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE */
/*     EXPANSION FOR THE IMAGINARY PART. */
/* ----------------------------------------------------------------------- */
    aez = az * 8.;
    s = *tol / aez;
    jl = (integer) ((real) (*rl + *rl)) + 2;
    p1r = zeror;
    p1i = zeroi;
    if (*zi == 0.) {
	goto L30;
    }
/* ----------------------------------------------------------------------- */
/*     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF */
/*     SIGNIFICANCE WHEN FNU OR N IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = (integer) ((real) (*fnu));
    arg = (*fnu - (doublereal) ((real) inu)) * pi;
    inu = inu + *n - il;
    ak = -sin(arg);
    bk = cos(arg);
    if (*zi < 0.) {
	bk = -bk;
    }
    p1r = ak;
    p1i = bk;
    if (inu % 2 == 0) {
	goto L30;
    }
    p1r = -p1r;
    p1i = -p1i;
L30:
    i__1 = il;
    for (k = 1; k <= i__1; ++k) {
	sqk = fdn - 1.;
	atol = s * abs(sqk);
	sgn = 1.;
	cs1r = coner;
	cs1i = conei;
	cs2r = coner;
	cs2i = conei;
	ckr = coner;
	cki = conei;
	ak = 0.;
	aa = 1.;
	bb = aez;
	dkr = ezr;
	dki = ezi;
	i__2 = jl;
	for (j = 1; j <= i__2; ++j) {
	    zdiv_(&ckr, &cki, &dkr, &dki, &str, &sti);
	    ckr = str * sqk;
	    cki = sti * sqk;
	    cs2r += ckr;
	    cs2i += cki;
	    sgn = -sgn;
	    cs1r += ckr * sgn;
	    cs1i += cki * sgn;
	    dkr += ezr;
	    dki += ezi;
	    aa = aa * abs(sqk) / bb;
	    bb += aez;
	    ak += 8.;
	    sqk -= ak;
	    if (aa <= atol) {
		goto L50;
	    }
/* L40: */
	}
	goto L110;
L50:
	s2r = cs1r;
	s2i = cs1i;
	if (*zr + *zr >= *elim) {
	    goto L60;
	}
	tzr = *zr + *zr;
	tzi = *zi + *zi;
	d__1 = -tzr;
	d__2 = -tzi;
	azexp_(&d__1, &d__2, &str, &sti);
	zmlt_(&str, &sti, &p1r, &p1i, &str, &sti);
	zmlt_(&str, &sti, &cs2r, &cs2i, &str, &sti);
	s2r += str;
	s2i += sti;
L60:
	fdn = fdn + dfnu * 8. + 4.;
	p1r = -p1r;
	p1i = -p1i;
	m = *n - il + k;
	yr[m] = s2r * ak1r - s2i * ak1i;
	yi[m] = s2r * ak1i + s2i * ak1r;
/* L70: */
    }
    if (*n <= 2) {
	return 0;
    }
    nn = *n;
    k = nn - 2;
    ak = (doublereal) ((real) k);
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    ib = 3;
    i__1 = nn;
    for (i__ = ib; i__ <= i__1; ++i__) {
	yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
	yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
	ak += -1.;
	--k;
/* L80: */
    }
    if (koded == 0) {
	return 0;
    }
    azexp_(&czr, &czi, &ckr, &cki);
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	str = yr[i__] * ckr - yi[i__] * cki;
	yi[i__] = yr[i__] * cki + yi[i__] * ckr;
	yr[i__] = str;
/* L90: */
    }
    return 0;
L100:
    *nz = -1;
    return 0;
L110:
    *nz = -2;
    return 0;
} /* zasyi_ */

