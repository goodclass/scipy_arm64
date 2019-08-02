/* zunk2.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zunk2_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *mr, integer *n, doublereal *yr, doublereal *
	yi, integer *nz, doublereal *tol, doublereal *elim, doublereal *alim)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal aic = 1.26551212348464539;
    static doublereal cipr[4] = { 1.,0.,-1.,0. };
    static doublereal cipi[4] = { 0.,-1.,0.,1. };
    static doublereal zeroi = 0.;
    static doublereal coner = 1.;
    static doublereal cr1r = 1.;
    static doublereal cr1i = 1.73205080756887729;
    static doublereal cr2r = -.5;
    static doublereal cr2i = -.866025403784438647;
    static doublereal hpi = 1.57079632679489662;
    static doublereal pi = 3.14159265358979324;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), log(doublereal), exp(doublereal),
	     d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, ib, ic;
    static doublereal fn;
    static integer il, kk, in, nw;
    static doublereal yy, c1i, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, 
	    aii, ang, asc, car, cki, fnf;
    static integer nai;
    static doublereal air;
    static integer ifn;
    static doublereal csi, ckr;
    static integer iuf;
    static doublereal cyi[2], fmr, sar, csr, sgn, zbi;
    static integer inu;
    static doublereal bry[3], cyr[2], pti, sti, zbr, zni, rzi, ptr, zri, str, 
	    znr, rzr, zrr, daii, aarg;
    static integer ndai;
    static doublereal dair, aphi, argi[2], cscl, phii[2], crsc, argr[2];
    static integer idum;
    static doublereal phir[2], csrr[3], cssr[3], rast, razr;
    extern /* Subroutine */ int zs1s2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static integer iflag, kflag;
    static doublereal argdi, ascle;
    static integer kdflg;
    static doublereal phidi, argdr;
    extern doublereal azabs_(doublereal *, doublereal *);
    static integer ipard;
    static doublereal csgni, phidr, cspni, asumi[2], bsumi[2];
    extern /* Subroutine */ int zuchk_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal cspnr, asumr[2], bsumr[2];
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int zunhj_(doublereal *, doublereal *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), zairy_(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal zeta1i[2], zeta2i[2], zet1di, zet2di, zeta1r[2], zeta2r[
	    2], zet1dr, zet2dr, asumdi, bsumdi, asumdr, bsumdr;

/* ***BEGIN PROLOGUE  ZUNK2 */
/* ***REFER TO  ZBESK */

/*     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE */
/*     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE */
/*     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN) */
/*     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR */
/*     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT */
/*     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC- */
/*     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION. */
/*     NZ=-1 MEANS AN OVERFLOW WILL OCCUR */

/* ***ROUTINES CALLED  ZAIRY,ZKSCL,ZS1S2,ZUCHK,ZUNHJ,D1MACH,AZABS */
/* ***END PROLOGUE  ZUNK2 */
/*     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC, */
/*    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ, */
/*    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR */
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
    yy = zri;
    znr = zri;
    zni = -zrr;
    zbr = zrr;
    zbi = zri;
    inu = (integer) ((real) (*fnu));
    fnf = *fnu - (doublereal) ((real) inu);
    ang = -hpi * fnf;
    car = cos(ang);
    sar = sin(ang);
    c2r = hpi * sar;
    c2i = -hpi * car;
    kk = inu % 4 + 1;
    str = c2r * cipr[kk - 1] - c2i * cipi[kk - 1];
    sti = c2r * cipi[kk - 1] + c2i * cipr[kk - 1];
    csr = cr1r * str - cr1i * sti;
    csi = cr1r * sti + cr1i * str;
    if (yy > 0.) {
	goto L20;
    }
    znr = -znr;
    zbi = -zbi;
L20:
/* ----------------------------------------------------------------------- */
/*     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST */
/*     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY */
/*     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS */
/* ----------------------------------------------------------------------- */
    j = 2;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* ----------------------------------------------------------------------- */
/*     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J */
/* ----------------------------------------------------------------------- */
	j = 3 - j;
	fn = *fnu + (doublereal) ((real) (i__ - 1));
	zunhj_(&znr, &zni, &fn, &c__0, tol, &phir[j - 1], &phii[j - 1], &argr[
		j - 1], &argi[j - 1], &zeta1r[j - 1], &zeta1i[j - 1], &zeta2r[
		j - 1], &zeta2i[j - 1], &asumr[j - 1], &asumi[j - 1], &bsumr[
		j - 1], &bsumi[j - 1]);
	if (*kode == 1) {
	    goto L30;
	}
	str = zbr + zeta2r[j - 1];
	sti = zbi + zeta2i[j - 1];
	rast = fn / azabs_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = zeta1r[j - 1] - str;
	s1i = zeta1i[j - 1] - sti;
	goto L40;
L30:
	s1r = zeta1r[j - 1] - zeta2r[j - 1];
	s1i = zeta1i[j - 1] - zeta2i[j - 1];
L40:
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (abs(rs1) > *elim) {
	    goto L70;
	}
	if (kdflg == 1) {
	    kflag = 2;
	}
	if (abs(rs1) < *alim) {
	    goto L50;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
	aphi = azabs_(&phir[j - 1], &phii[j - 1]);
	aarg = azabs_(&argr[j - 1], &argi[j - 1]);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (abs(rs1) > *elim) {
	    goto L70;
	}
	if (kdflg == 1) {
	    kflag = 1;
	}
	if (rs1 < 0.) {
	    goto L50;
	}
	if (kdflg == 1) {
	    kflag = 3;
	}
L50:
/* ----------------------------------------------------------------------- */
/*     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR */
/*     EXPONENT EXTREMES */
/* ----------------------------------------------------------------------- */
	c2r = argr[j - 1] * cr2r - argi[j - 1] * cr2i;
	c2i = argr[j - 1] * cr2i + argi[j - 1] * cr2r;
	zairy_(&c2r, &c2i, &c__0, &c__2, &air, &aii, &nai, &idum);
	zairy_(&c2r, &c2i, &c__1, &c__2, &dair, &daii, &ndai, &idum);
	str = dair * bsumr[j - 1] - daii * bsumi[j - 1];
	sti = dair * bsumi[j - 1] + daii * bsumr[j - 1];
	ptr = str * cr2r - sti * cr2i;
	pti = str * cr2i + sti * cr2r;
	str = ptr + (air * asumr[j - 1] - aii * asumi[j - 1]);
	sti = pti + (air * asumi[j - 1] + aii * asumr[j - 1]);
	ptr = str * phir[j - 1] - sti * phii[j - 1];
	pti = str * phii[j - 1] + sti * phir[j - 1];
	s2r = ptr * csr - pti * csi;
	s2i = ptr * csi + pti * csr;
	str = exp(s1r) * cssr[kflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s1r * s2i + s2r * s1i;
	s2r = str;
	if (kflag != 1) {
	    goto L60;
	}
	zuchk_(&s2r, &s2i, &nw, bry, tol);
	if (nw != 0) {
	    goto L70;
	}
L60:
	if (yy <= 0.) {
	    s2i = -s2i;
	}
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	yr[i__] = s2r * csrr[kflag - 1];
	yi[i__] = s2i * csrr[kflag - 1];
	str = csi;
	csi = -csr;
	csr = str;
	if (kdflg == 2) {
	    goto L85;
	}
	kdflg = 2;
	goto L80;
L70:
	if (rs1 > 0.) {
	    goto L320;
	}
/* ----------------------------------------------------------------------- */
/*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
/* ----------------------------------------------------------------------- */
	if (*zr < 0.) {
	    goto L320;
	}
	kdflg = 1;
	yr[i__] = zeror;
	yi[i__] = zeroi;
	++(*nz);
	str = csi;
	csi = -csr;
	csr = str;
	if (i__ == 1) {
	    goto L80;
	}
	if (yr[i__ - 1] == zeror && yi[i__ - 1] == zeroi) {
	    goto L80;
	}
	yr[i__ - 1] = zeror;
	yi[i__ - 1] = zeroi;
	++(*nz);
L80:
	;
    }
    i__ = *n;
L85:
    razr = 1. / azabs_(&zrr, &zri);
    str = zrr * razr;
    sti = -zri * razr;
    rzr = (str + str) * razr;
    rzi = (sti + sti) * razr;
    ckr = fn * rzr;
    cki = fn * rzi;
    ib = i__ + 1;
    if (*n < ib) {
	goto L180;
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
    zunhj_(&znr, &zni, &fn, &ipard, tol, &phidr, &phidi, &argdr, &argdi, &
	    zet1dr, &zet1di, &zet2dr, &zet2di, &asumdr, &asumdi, &bsumdr, &
	    bsumdi);
    if (*kode == 1) {
	goto L90;
    }
    str = zbr + zet2dr;
    sti = zbi + zet2di;
    rast = fn / azabs_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = zet1dr - str;
    s1i = zet1di - sti;
    goto L100;
L90:
    s1r = zet1dr - zet2dr;
    s1i = zet1di - zet2di;
L100:
    rs1 = s1r;
    if (abs(rs1) > *elim) {
	goto L105;
    }
    if (abs(rs1) < *alim) {
	goto L120;
    }
/* ---------------------------------------------------------------------------- */
/*     REFINE ESTIMATE AND TEST */
/* ------------------------------------------------------------------------- */
    aphi = azabs_(&phidr, &phidi);
    rs1 += log(aphi);
    if (abs(rs1) < *elim) {
	goto L120;
    }
L105:
    if (rs1 > 0.) {
	goto L320;
    }
/* ----------------------------------------------------------------------- */
/*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
/* ----------------------------------------------------------------------- */
    if (*zr < 0.) {
	goto L320;
    }
    *nz = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yr[i__] = zeror;
	yi[i__] = zeroi;
/* L106: */
    }
    return 0;
L120:
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
	    goto L130;
	}
	str = abs(c2r);
	sti = abs(c2i);
	c2m = max(str,sti);
	if (c2m <= ascle) {
	    goto L130;
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
L130:
	;
    }
L180:
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
    if (yy <= 0.) {
	csgni = -csgni;
    }
    ifn = inu + *n - 1;
    ang = fnf * sgn;
    cspnr = cos(ang);
    cspni = sin(ang);
    if (ifn % 2 == 0) {
	goto L190;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L190:
/* ----------------------------------------------------------------------- */
/*     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS */
/*     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST */
/*     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY */
/*     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS */
/* ----------------------------------------------------------------------- */
    csr = sar * csgni;
    csi = car * csgni;
    in = ifn % 4 + 1;
    c2r = cipr[in - 1];
    c2i = cipi[in - 1];
    str = csr * c2r + csi * c2i;
    csi = -csr * c2i + csi * c2r;
    csr = str;
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
	if (*n > 2) {
	    goto L175;
	}
L172:
	phidr = phir[j - 1];
	phidi = phii[j - 1];
	argdr = argr[j - 1];
	argdi = argi[j - 1];
	zet1dr = zeta1r[j - 1];
	zet1di = zeta1i[j - 1];
	zet2dr = zeta2r[j - 1];
	zet2di = zeta2i[j - 1];
	asumdr = asumr[j - 1];
	asumdi = asumi[j - 1];
	bsumdr = bsumr[j - 1];
	bsumdi = bsumi[j - 1];
	j = 3 - j;
	goto L210;
L175:
	if (kk == *n && ib < *n) {
	    goto L210;
	}
	if (kk == ib || kk == ic) {
	    goto L172;
	}
	zunhj_(&znr, &zni, &fn, &c__0, tol, &phidr, &phidi, &argdr, &argdi, &
		zet1dr, &zet1di, &zet2dr, &zet2di, &asumdr, &asumdi, &bsumdr, 
		&bsumdi);
L210:
	if (*kode == 1) {
	    goto L220;
	}
	str = zbr + zet2dr;
	sti = zbi + zet2di;
	rast = fn / azabs_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zet1dr + str;
	s1i = -zet1di + sti;
	goto L230;
L220:
	s1r = -zet1dr + zet2dr;
	s1i = -zet1di + zet2di;
L230:
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (abs(rs1) > *elim) {
	    goto L280;
	}
	if (kdflg == 1) {
	    iflag = 2;
	}
	if (abs(rs1) < *alim) {
	    goto L240;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
	aphi = azabs_(&phidr, &phidi);
	aarg = azabs_(&argdr, &argdi);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (abs(rs1) > *elim) {
	    goto L280;
	}
	if (kdflg == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L240;
	}
	if (kdflg == 1) {
	    iflag = 3;
	}
L240:
	zairy_(&argdr, &argdi, &c__0, &c__2, &air, &aii, &nai, &idum);
	zairy_(&argdr, &argdi, &c__1, &c__2, &dair, &daii, &ndai, &idum);
	str = dair * bsumdr - daii * bsumdi;
	sti = dair * bsumdi + daii * bsumdr;
	str += air * asumdr - aii * asumdi;
	sti += air * asumdi + aii * asumdr;
	ptr = str * phidr - sti * phidi;
	pti = str * phidi + sti * phidr;
	s2r = ptr * csr - pti * csi;
	s2i = ptr * csi + pti * csr;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L250;
	}
	zuchk_(&s2r, &s2i, &nw, bry, tol);
	if (nw == 0) {
	    goto L250;
	}
	s2r = zeror;
	s2i = zeroi;
L250:
	if (yy <= 0.) {
	    s2i = -s2i;
	}
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
	    goto L270;
	}
	zs1s2_(&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim, &iuf);
	*nz += nw;
L270:
	yr[kk] = s1r * cspnr - s1i * cspni + s2r;
	yi[kk] = s1r * cspni + s1i * cspnr + s2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	str = csi;
	csi = -csr;
	csr = str;
	if (c2r != 0. || c2i != 0.) {
	    goto L255;
	}
	kdflg = 1;
	goto L290;
L255:
	if (kdflg == 2) {
	    goto L295;
	}
	kdflg = 2;
	goto L290;
L280:
	if (rs1 > 0.) {
	    goto L320;
	}
	s2r = zeror;
	s2i = zeroi;
	goto L250;
L290:
	;
    }
    k = *n;
L295:
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
	    goto L300;
	}
	zs1s2_(&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim, &iuf);
	*nz += nw;
L300:
	yr[kk] = c1r * cspnr - c1i * cspni + c2r;
	yi[kk] = c1r * cspni + c1i * cspnr + c2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	if (iflag >= 3) {
	    goto L310;
	}
	c2r = abs(ckr);
	c2i = abs(cki);
	c2m = max(c2r,c2i);
	if (c2m <= ascle) {
	    goto L310;
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
L310:
	;
    }
    return 0;
L320:
    *nz = -1;
    return 0;
} /* zunk2_ */

