/* zseri.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zseri_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *yr, doublereal *yi, integer *
	nz, doublereal *tol, doublereal *elim, doublereal *alim)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;
    static doublereal coner = 1.;
    static doublereal conei = 0.;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer i__, k, l, m;
    static doublereal s, aa;
    static integer ib;
    static doublereal ak;
    static integer il;
    static doublereal az;
    static integer nn;
    static doublereal wi[2], rs, ss;
    static integer nw;
    static doublereal wr[2], s1i, s2i, s1r, s2r, cki, acz, arm, ckr, czi, hzi,
	     raz, czr, sti, hzr, rzi, str, rzr, ak1i, ak1r, rtr1, dfnu;
    static integer idum;
    static doublereal atol, fnup;
    extern /* Subroutine */ int zdiv_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *), zmlt_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer iflag;
    static doublereal coefi, ascle, coefr;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal crscr;
    extern /* Subroutine */ int azlog_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *), zuchk_(doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *);
    extern doublereal d1mach_(integer *), dgamln_(doublereal *, integer *);

/* ***BEGIN PROLOGUE  ZSERI */
/* ***REFER TO  ZBESI,ZBESK */

/*     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY */
/*     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE */
/*     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN. */
/*     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO */
/*     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE */
/*     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE */
/*     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ). */

/* ***ROUTINES CALLED  DGAMLN,D1MACH,ZUCHK,AZABS,ZDIV,AZLOG,ZMLT */
/* ***END PROLOGUE  ZSERI */
/*     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    az = azabs_(zr, zi);
    if (az == 0.) {
	goto L160;
    }
    arm = d1mach_(&c__1) * 1e3;
    rtr1 = sqrt(arm);
    crscr = 1.;
    iflag = 0;
    if (az < arm) {
	goto L150;
    }
    hzr = *zr * .5;
    hzi = *zi * .5;
    czr = zeror;
    czi = zeroi;
    if (az <= rtr1) {
	goto L10;
    }
    zmlt_(&hzr, &hzi, &hzr, &hzi, &czr, &czi);
L10:
    acz = azabs_(&czr, &czi);
    nn = *n;
    azlog_(&hzr, &hzi, &ckr, &cki, &idum);
L20:
    dfnu = *fnu + (doublereal) ((real) (nn - 1));
    fnup = dfnu + 1.;
/* ----------------------------------------------------------------------- */
/*     UNDERFLOW TEST */
/* ----------------------------------------------------------------------- */
    ak1r = ckr * dfnu;
    ak1i = cki * dfnu;
    ak = dgamln_(&fnup, &idum);
    ak1r -= ak;
    if (*kode == 2) {
	ak1r -= *zr;
    }
    if (ak1r > -(*elim)) {
	goto L40;
    }
L30:
    ++(*nz);
    yr[nn] = zeror;
    yi[nn] = zeroi;
    if (acz > dfnu) {
	goto L190;
    }
    --nn;
    if (nn == 0) {
	return 0;
    }
    goto L20;
L40:
    if (ak1r > -(*alim)) {
	goto L50;
    }
    iflag = 1;
    ss = 1. / *tol;
    crscr = *tol;
    ascle = arm * ss;
L50:
    aa = exp(ak1r);
    if (iflag == 1) {
	aa *= ss;
    }
    coefr = aa * cos(ak1i);
    coefi = aa * sin(ak1i);
    atol = *tol * acz / fnup;
    il = min(2,nn);
    i__1 = il;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dfnu = *fnu + (doublereal) ((real) (nn - i__));
	fnup = dfnu + 1.;
	s1r = coner;
	s1i = conei;
	if (acz < *tol * fnup) {
	    goto L70;
	}
	ak1r = coner;
	ak1i = conei;
	ak = fnup + 2.;
	s = fnup;
	aa = 2.;
L60:
	rs = 1. / s;
	str = ak1r * czr - ak1i * czi;
	sti = ak1r * czi + ak1i * czr;
	ak1r = str * rs;
	ak1i = sti * rs;
	s1r += ak1r;
	s1i += ak1i;
	s += ak;
	ak += 2.;
	aa = aa * acz * rs;
	if (aa > atol) {
	    goto L60;
	}
L70:
	s2r = s1r * coefr - s1i * coefi;
	s2i = s1r * coefi + s1i * coefr;
	wr[i__ - 1] = s2r;
	wi[i__ - 1] = s2i;
	if (iflag == 0) {
	    goto L80;
	}
	zuchk_(&s2r, &s2i, &nw, &ascle, tol);
	if (nw != 0) {
	    goto L30;
	}
L80:
	m = nn - i__ + 1;
	yr[m] = s2r * crscr;
	yi[m] = s2i * crscr;
	if (i__ == il) {
	    goto L90;
	}
	zdiv_(&coefr, &coefi, &hzr, &hzi, &str, &sti);
	coefr = str * dfnu;
	coefi = sti * dfnu;
L90:
	;
    }
    if (nn <= 2) {
	return 0;
    }
    k = nn - 2;
    ak = (doublereal) ((real) k);
    raz = 1. / az;
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    if (iflag == 1) {
	goto L120;
    }
    ib = 3;
L100:
    i__1 = nn;
    for (i__ = ib; i__ <= i__1; ++i__) {
	yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
	yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
	ak += -1.;
	--k;
/* L110: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     RECUR BACKWARD WITH SCALED VALUES */
/* ----------------------------------------------------------------------- */
L120:
/* ----------------------------------------------------------------------- */
/*     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE */
/*     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3 */
/* ----------------------------------------------------------------------- */
    s1r = wr[0];
    s1i = wi[0];
    s2r = wr[1];
    s2i = wi[1];
    i__1 = nn;
    for (l = 3; l <= i__1; ++l) {
	ckr = s2r;
	cki = s2i;
	s2r = s1r + (ak + *fnu) * (rzr * ckr - rzi * cki);
	s2i = s1i + (ak + *fnu) * (rzr * cki + rzi * ckr);
	s1r = ckr;
	s1i = cki;
	ckr = s2r * crscr;
	cki = s2i * crscr;
	yr[k] = ckr;
	yi[k] = cki;
	ak += -1.;
	--k;
	if (azabs_(&ckr, &cki) > ascle) {
	    goto L140;
	}
/* L130: */
    }
    return 0;
L140:
    ib = l + 1;
    if (ib > nn) {
	return 0;
    }
    goto L100;
L150:
    *nz = *n;
    if (*fnu == 0.) {
	--(*nz);
    }
L160:
    yr[1] = zeror;
    yi[1] = zeroi;
    if (*fnu != 0.) {
	goto L170;
    }
    yr[1] = coner;
    yi[1] = conei;
L170:
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	yr[i__] = zeror;
	yi[i__] = zeroi;
/* L180: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE */
/*     THE CALCULATION IN CBINU WITH N=N-IABS(NZ) */
/* ----------------------------------------------------------------------- */
L190:
    *nz = -(*nz);
    return 0;
} /* zseri_ */

