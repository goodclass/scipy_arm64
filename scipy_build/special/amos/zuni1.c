/* zuni1.f -- translated by f2c (version 20190311).
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
static integer c__0 = 0;
static integer c__2 = 2;

/* Subroutine */ int zuni1_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *yr, doublereal *yi, integer *
	nz, integer *nlast, doublereal *fnul, doublereal *tol, doublereal *
	elim, doublereal *alim)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;
    static doublereal coner = 1.;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, m, nd;
    static doublereal fn;
    static integer nn, nw;
    static doublereal c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, cyi[2];
    static integer nuf;
    static doublereal bry[3], cyr[2], sti, rzi, str, rzr, aphi, cscl, phii, 
	    crsc, phir;
    static integer init;
    static doublereal csrr[3], cssr[3], rast, sumi, sumr;
    static integer iflag;
    static doublereal ascle;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal cwrki[16];
    extern /* Subroutine */ int zuchk_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal cwrkr[16];
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int zunik_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), zuoik_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    static doublereal zeta1i, zeta2i, zeta1r, zeta2r;

/* ***BEGIN PROLOGUE  ZUNI1 */
/* ***REFER TO  ZBESI,ZBESK */

/*     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC */
/*     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3. */

/*     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC */
/*     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET. */
/*     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER */
/*     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL. */
/*     Y(I)=CZERO FOR I=NLAST+1,N */

/* ***ROUTINES CALLED  ZUCHK,ZUNIK,ZUOIK,D1MACH,AZABS */
/* ***END PROLOGUE  ZUNI1 */
/*     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1, */
/*    *S2,Y,Z,ZETA1,ZETA2 */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    nd = *n;
    *nlast = 0;
/* ----------------------------------------------------------------------- */
/*     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG- */
/*     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE, */
/*     EXP(ALIM)=EXP(ELIM)*TOL */
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
/* ----------------------------------------------------------------------- */
/*     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER */
/* ----------------------------------------------------------------------- */
    fn = max(*fnu,1.);
    init = 0;
    zunik_(zr, zi, &fn, &c__1, &c__1, tol, &init, &phir, &phii, &zeta1r, &
	    zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
    if (*kode == 1) {
	goto L10;
    }
    str = *zr + zeta2r;
    sti = *zi + zeta2i;
    rast = fn / azabs_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti;
    goto L20;
L10:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
L20:
    rs1 = s1r;
    if (abs(rs1) > *elim) {
	goto L130;
    }
L30:
    nn = min(2,nd);
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fn = *fnu + (doublereal) ((real) (nd - i__));
	init = 0;
	zunik_(zr, zi, &fn, &c__1, &c__0, tol, &init, &phir, &phii, &zeta1r, &
		zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
	if (*kode == 1) {
	    goto L40;
	}
	str = *zr + zeta2r;
	sti = *zi + zeta2i;
	rast = fn / azabs_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zeta1r + str;
	s1i = -zeta1i + sti + *zi;
	goto L50;
L40:
	s1r = -zeta1r + zeta2r;
	s1i = -zeta1i + zeta2i;
L50:
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (abs(rs1) > *elim) {
	    goto L110;
	}
	if (i__ == 1) {
	    iflag = 2;
	}
	if (abs(rs1) < *alim) {
	    goto L60;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
	aphi = azabs_(&phir, &phii);
	rs1 += log(aphi);
	if (abs(rs1) > *elim) {
	    goto L110;
	}
	if (i__ == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L60;
	}
	if (i__ == 1) {
	    iflag = 3;
	}
L60:
/* ----------------------------------------------------------------------- */
/*     SCALE S1 IF CABS(S1).LT.ASCLE */
/* ----------------------------------------------------------------------- */
	s2r = phir * sumr - phii * sumi;
	s2i = phir * sumi + phii * sumr;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L70;
	}
	zuchk_(&s2r, &s2i, &nw, bry, tol);
	if (nw != 0) {
	    goto L110;
	}
L70:
	cyr[i__ - 1] = s2r;
	cyi[i__ - 1] = s2i;
	m = nd - i__ + 1;
	yr[m] = s2r * csrr[iflag - 1];
	yi[m] = s2i * csrr[iflag - 1];
/* L80: */
    }
    if (nd <= 2) {
	goto L100;
    }
    rast = 1. / azabs_(zr, zi);
    str = *zr * rast;
    sti = -(*zi) * rast;
    rzr = (str + str) * rast;
    rzi = (sti + sti) * rast;
    bry[1] = 1. / bry[0];
    bry[2] = d1mach_(&c__2);
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    k = nd - 2;
    fn = (doublereal) ((real) k);
    i__1 = nd;
    for (i__ = 3; i__ <= i__1; ++i__) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (*fnu + fn) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (*fnu + fn) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[k] = c2r;
	yi[k] = c2i;
	--k;
	fn += -1.;
	if (iflag >= 3) {
	    goto L90;
	}
	str = abs(c2r);
	sti = abs(c2i);
	c2m = max(str,sti);
	if (c2m <= ascle) {
	    goto L90;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	c1r = csrr[iflag - 1];
L90:
	;
    }
L100:
    return 0;
/* ----------------------------------------------------------------------- */
/*     SET UNDERFLOW AND UPDATE PARAMETERS */
/* ----------------------------------------------------------------------- */
L110:
    if (rs1 > 0.) {
	goto L120;
    }
    yr[nd] = zeror;
    yi[nd] = zeroi;
    ++(*nz);
    --nd;
    if (nd == 0) {
	goto L100;
    }
    zuoik_(zr, zi, fnu, kode, &c__1, &nd, &yr[1], &yi[1], &nuf, tol, elim, 
	    alim);
    if (nuf < 0) {
	goto L120;
    }
    nd -= nuf;
    *nz += nuf;
    if (nd == 0) {
	goto L100;
    }
    fn = *fnu + (doublereal) ((real) (nd - 1));
    if (fn >= *fnul) {
	goto L30;
    }
    *nlast = nd;
    return 0;
L120:
    *nz = -1;
    return 0;
L130:
    if (rs1 > 0.) {
	goto L120;
    }
    *nz = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yr[i__] = zeror;
	yi[i__] = zeroi;
/* L140: */
    }
    return 0;
} /* zuni1_ */

