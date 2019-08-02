/* zunk1.f -- translated by f2c (version 20190311).
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
static integer c__2 = 2;
static integer c__0 = 0;

/* Subroutine */ int zunk1_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *mr, integer *n, doublereal *yr, doublereal *
	yi, integer *nz, doublereal *tol, doublereal *elim, doublereal *alim)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;
    static doublereal coner = 1.;
    static doublereal pi = 3.14159265358979324;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), cos(doublereal), sin(doublereal),
	     d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, m, ib, ic;
    static doublereal fn;
    static integer il, kk, nw;
    static doublereal c1i, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, ang, 
	    asc, cki, fnf;
    static integer ifn;
    static doublereal ckr;
    static integer iuf;
    static doublereal cyi[2], fmr, csr, sgn;
    static integer inu;
    static doublereal bry[3], cyr[2], sti, rzi, zri, str, rzr, zrr, aphi, 
	    cscl, phii[2], crsc, phir[2];
    static integer init[2];
    static doublereal csrr[3], cssr[3], rast, sumi[2], razr;
    extern /* Subroutine */ int zs1s2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal sumr[2];
    static integer iflag, kflag;
    static doublereal ascle;
    static integer kdflg;
    static doublereal phidi;
    static integer ipard;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal csgni, phidr;
    static integer initd;
    static doublereal cspni, cwrki[48]	/* was [16][3] */, sumdi;
    extern /* Subroutine */ int zuchk_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal cspnr, cwrkr[48]	/* was [16][3] */, sumdr;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int zunik_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal zeta1i[2], zeta2i[2], zet1di, zet2di, zeta1r[2], zeta2r[
	    2], zet1dr, zet2dr;

/* ***BEGIN PROLOGUE  ZUNK1 */
/* ***REFER TO  ZBESK */

/*     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE */
/*     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE */
/*     UNIFORM ASYMPTOTIC EXPANSION. */
/*     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION. */
/*     NZ=-1 MEANS AN OVERFLOW WILL OCCUR */

/* ***ROUTINES CALLED  ZKSCL,ZS1S2,ZUCHK,ZUNIK,D1MACH,AZABS */
/* ***END PROLOGUE  ZUNK1 */
/*     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO, */
/*    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    kdflg = 1;
    *nz = 0;
/* ----------------------------------------------------------------------- */
/*     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN */
/*     THE UNDERFLOW LIMIT */
/* ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    crsc = *tol;
    cssr[0] = cscl;
    cssr[1] = coner;
    cssr[2] = crsc;
    csrr[0] = crsc;
    csrr[1] = coner;
    csrr[2] = cscl;
    bry[0] = d1mach_(&c__1) * 1e3 / *tol;
    bry[1] = 1. / bry[0];
    bry[2] = d1mach_(&c__2);
    zrr = *zr;
    zri = *zi;
    if (*zr >= 0.) {
	goto L10;
    }
    zrr = -(*zr);
    zri = -(*zi);
L10:
    j = 2;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* ----------------------------------------------------------------------- */
/*     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J */
/* ----------------------------------------------------------------------- */
	j = 3 - j;
	fn = *fnu + (doublereal) ((real) (i__ - 1));
	init[j - 1] = 0;
	zunik_(&zrr, &zri, &fn, &c__2, &c__0, tol, &init[j - 1], &phir[j - 1],
		 &phii[j - 1], &zeta1r[j - 1], &zeta1i[j - 1], &zeta2r[j - 1],
		 &zeta2i[j - 1], &sumr[j - 1], &sumi[j - 1], &cwrkr[(j << 4) 
		- 16], &cwrki[(j << 4) - 16]);
	if (*kode == 1) {
	    goto L20;
	}
	str = zrr + zeta2r[j - 1];
	sti = zri + zeta2i[j - 1];
	rast = fn / azabs_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = zeta1r[j - 1] - str;
	s1i = zeta1i[j - 1] - sti;
	goto L30;
L20:
	s1r = zeta1r[j - 1] - zeta2r[j - 1];
	s1i = zeta1i[j - 1] - zeta2i[j - 1];
L30:
	rs1 = s1r;
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	if (abs(rs1) > *elim) {
	    goto L60;
	}
	if (kdflg == 1) {
	    kflag = 2;
	}
	if (abs(rs1) < *alim) {
	    goto L40;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
	aphi = azabs_(&phir[j - 1], &phii[j - 1]);
	rs1 += log(aphi);
	if (abs(rs1) > *elim) {
	    goto L60;
	}
	if (kdflg == 1) {
	    kflag = 1;
	}
	if (rs1 < 0.) {
	    goto L40;
	}
	if (kdflg == 1) {
	    kflag = 3;
	}
L40:
/* ----------------------------------------------------------------------- */
/*     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR */
/*     EXPONENT EXTREMES */
/* ----------------------------------------------------------------------- */
	s2r = phir[j - 1] * sumr[j - 1] - phii[j - 1] * sumi[j - 1];
	s2i = phir[j - 1] * sumi[j - 1] + phii[j - 1] * sumr[j - 1];
	str = exp(s1r) * cssr[kflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s1r * s2i + s2r * s1i;
	s2r = str;
	if (kflag != 1) {
	    goto L50;
	}
	zuchk_(&s2r, &s2i, &nw, bry, tol);
	if (nw != 0) {
	    goto L60;
	}
L50:
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	yr[i__] = s2r * csrr[kflag - 1];
	yi[i__] = s2i * csrr[kflag - 1];
	if (kdflg == 2) {
	    goto L75;
	}
	kdflg = 2;
	goto L70;
L60:
	if (rs1 > 0.) {
	    goto L300;
	}
/* ----------------------------------------------------------------------- */
/*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
/* ----------------------------------------------------------------------- */
	if (*zr < 0.) {
	    goto L300;
	}
	kdflg = 1;
	yr[i__] = zeror;
	yi[i__] = zeroi;
	++(*nz);
	if (i__ == 1) {
	    goto L70;
	}
	if (yr[i__ - 1] == zeror && yi[i__ - 1] == zeroi) {
	    goto L70;
	}
	yr[i__ - 1] = zeror;
	yi[i__ - 1] = zeroi;
	++(*nz);
L70:
	;
    }
    i__ = *n;
L75:
    razr = 1. / azabs_(&zrr, &zri);
    str = zrr * razr;
    sti = -zri * razr;
    rzr = (str + str) * razr;
    rzi = (sti + sti) * razr;
    ckr = fn * rzr;
    cki = fn * rzi;
    ib = i__ + 1;
    if (*n < ib) {
	goto L160;
    }
/* ----------------------------------------------------------------------- */
/*     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO */
/*     ON UNDERFLOW. */
/* ----------------------------------------------------------------------- */
    fn = *fnu + (doublereal) ((real) (*n - 1));
    ipard = 1;
    if (*mr != 0) {
	ipard = 0;
    }
    initd = 0;
    zunik_(&zrr, &zri, &fn, &c__2, &ipard, tol, &initd, &phidr, &phidi, &
	    zet1dr, &zet1di, &zet2dr, &zet2di, &sumdr, &sumdi, &cwrkr[32], &
	    cwrki[32]);
    if (*kode == 1) {
	goto L80;
    }
    str = zrr + zet2dr;
    sti = zri + zet2di;
    rast = fn / azabs_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = zet1dr - str;
    s1i = zet1di - sti;
    goto L90;
L80:
    s1r = zet1dr - zet2dr;
    s1i = zet1di - zet2di;
L90:
    rs1 = s1r;
    if (abs(rs1) > *elim) {
	goto L95;
    }
    if (abs(rs1) < *alim) {
	goto L100;
    }
/* ---------------------------------------------------------------------------- */
/*     REFINE ESTIMATE AND TEST */
/* ------------------------------------------------------------------------- */
    aphi = azabs_(&phidr, &phidi);
    rs1 += log(aphi);
    if (abs(rs1) < *elim) {
	goto L100;
    }
L95:
    if (abs(rs1) > 0.) {
	goto L300;
    }
/* ----------------------------------------------------------------------- */
/*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
/* ----------------------------------------------------------------------- */
    if (*zr < 0.) {
	goto L300;
    }
    *nz = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yr[i__] = zeror;
	yi[i__] = zeroi;
/* L96: */
    }
    return 0;
/* --------------------------------------------------------------------------- */
/*     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE */
/* ---------------------------------------------------------------------------- */
L100:
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[kflag - 1];
    ascle = bry[kflag - 1];
    i__1 = *n;
    for (i__ = ib; i__ <= i__1; ++i__) {
	c2r = s2r;
	c2i = s2i;
	s2r = ckr * c2r - cki * c2i + s1r;
	s2i = ckr * c2i + cki * c2r + s1i;
	s1r = c2r;
	s1i = c2i;
	ckr += rzr;
	cki += rzi;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[i__] = c2r;
	yi[i__] = c2i;
	if (kflag >= 3) {
	    goto L120;
	}
	str = abs(c2r);
	sti = abs(c2i);
	c2m = max(str,sti);
	if (c2m <= ascle) {
	    goto L120;
	}
	++kflag;
	ascle = bry[kflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[kflag - 1];
	s1i *= cssr[kflag - 1];
	s2r *= cssr[kflag - 1];
	s2i *= cssr[kflag - 1];
	c1r = csrr[kflag - 1];
L120:
	;
    }
L160:
    if (*mr == 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0 */
/* ----------------------------------------------------------------------- */
    *nz = 0;
    fmr = (doublereal) ((real) (*mr));
    sgn = -d_sign(&pi, &fmr);
/* ----------------------------------------------------------------------- */
/*     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP. */
/* ----------------------------------------------------------------------- */
    csgni = sgn;
    inu = (integer) ((real) (*fnu));
    fnf = *fnu - (doublereal) ((real) inu);
    ifn = inu + *n - 1;
    ang = fnf * sgn;
    cspnr = cos(ang);
    cspni = sin(ang);
    if (ifn % 2 == 0) {
	goto L170;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L170:
    asc = bry[0];
    iuf = 0;
    kk = *n;
    kdflg = 1;
    --ib;
    ic = ib - 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	fn = *fnu + (doublereal) ((real) (kk - 1));
/* ----------------------------------------------------------------------- */
/*     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K */
/*     FUNCTION ABOVE */
/* ----------------------------------------------------------------------- */
	m = 3;
	if (*n > 2) {
	    goto L175;
	}
L172:
	initd = init[j - 1];
	phidr = phir[j - 1];
	phidi = phii[j - 1];
	zet1dr = zeta1r[j - 1];
	zet1di = zeta1i[j - 1];
	zet2dr = zeta2r[j - 1];
	zet2di = zeta2i[j - 1];
	sumdr = sumr[j - 1];
	sumdi = sumi[j - 1];
	m = j;
	j = 3 - j;
	goto L180;
L175:
	if (kk == *n && ib < *n) {
	    goto L180;
	}
	if (kk == ib || kk == ic) {
	    goto L172;
	}
	initd = 0;
L180:
	zunik_(&zrr, &zri, &fn, &c__1, &c__0, tol, &initd, &phidr, &phidi, &
		zet1dr, &zet1di, &zet2dr, &zet2di, &sumdr, &sumdi, &cwrkr[(m 
		<< 4) - 16], &cwrki[(m << 4) - 16]);
	if (*kode == 1) {
	    goto L200;
	}
	str = zrr + zet2dr;
	sti = zri + zet2di;
	rast = fn / azabs_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zet1dr + str;
	s1i = -zet1di + sti;
	goto L210;
L200:
	s1r = -zet1dr + zet2dr;
	s1i = -zet1di + zet2di;
L210:
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (abs(rs1) > *elim) {
	    goto L260;
	}
	if (kdflg == 1) {
	    iflag = 2;
	}
	if (abs(rs1) < *alim) {
	    goto L220;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
	aphi = azabs_(&phidr, &phidi);
	rs1 += log(aphi);
	if (abs(rs1) > *elim) {
	    goto L260;
	}
	if (kdflg == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L220;
	}
	if (kdflg == 1) {
	    iflag = 3;
	}
L220:
	str = phidr * sumdr - phidi * sumdi;
	sti = phidr * sumdi + phidi * sumdr;
	s2r = -csgni * sti;
	s2i = csgni * str;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L230;
	}
	zuchk_(&s2r, &s2i, &nw, bry, tol);
	if (nw == 0) {
	    goto L230;
	}
	s2r = zeror;
	s2i = zeroi;
L230:
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	c2r = s2r;
	c2i = s2i;
	s2r *= csrr[iflag - 1];
	s2i *= csrr[iflag - 1];
/* ----------------------------------------------------------------------- */
/*     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N */
/* ----------------------------------------------------------------------- */
	s1r = yr[kk];
	s1i = yi[kk];
	if (*kode == 1) {
	    goto L250;
	}
	zs1s2_(&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim, &iuf);
	*nz += nw;
L250:
	yr[kk] = s1r * cspnr - s1i * cspni + s2r;
	yi[kk] = cspnr * s1i + cspni * s1r + s2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	if (c2r != 0. || c2i != 0.) {
	    goto L255;
	}
	kdflg = 1;
	goto L270;
L255:
	if (kdflg == 2) {
	    goto L275;
	}
	kdflg = 2;
	goto L270;
L260:
	if (rs1 > 0.) {
	    goto L300;
	}
	s2r = zeror;
	s2i = zeroi;
	goto L230;
L270:
	;
    }
    k = *n;
L275:
    il = *n - k;
    if (il == 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE */
/*     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP */
/*     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES. */
/* ----------------------------------------------------------------------- */
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    csr = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    fn = (doublereal) ((real) (inu + il));
    i__1 = il;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	fn += -1.;
	c2r = s2r * csr;
	c2i = s2i * csr;
	ckr = c2r;
	cki = c2i;
	c1r = yr[kk];
	c1i = yi[kk];
	if (*kode == 1) {
	    goto L280;
	}
	zs1s2_(&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim, &iuf);
	*nz += nw;
L280:
	yr[kk] = c1r * cspnr - c1i * cspni + c2r;
	yi[kk] = c1r * cspni + c1i * cspnr + c2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	if (iflag >= 3) {
	    goto L290;
	}
	c2r = abs(ckr);
	c2i = abs(cki);
	c2m = max(c2r,c2i);
	if (c2m <= ascle) {
	    goto L290;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= csr;
	s1i *= csr;
	s2r = ckr;
	s2i = cki;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	csr = csrr[iflag - 1];
L290:
	;
    }
    return 0;
L300:
    *nz = -1;
    return 0;
} /* zunk1_ */

