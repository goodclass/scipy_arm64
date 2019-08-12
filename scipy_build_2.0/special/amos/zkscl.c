/* zkscl.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zkscl_(doublereal *zrr, doublereal *zri, doublereal *fnu,
	 integer *n, doublereal *yr, doublereal *yi, integer *nz, doublereal *
	rzr, doublereal *rzi, doublereal *ascle, doublereal *tol, doublereal *
	elim)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, ic;
    static doublereal as, fn;
    static integer kk, nn, nw;
    static doublereal s1i, s2i, s1r, s2r, acs, cki, elm, csi, ckr, cyi[2], 
	    zdi, csr, cyr[2], zdr, str, alas;
    static integer idum;
    static doublereal helim;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal celmr;
    extern /* Subroutine */ int azlog_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *), zuchk_(doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZKSCL */
/* ***REFER TO  ZBESK */

/*     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE */
/*     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN */
/*     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL. */

/* ***ROUTINES CALLED  ZUCHK,AZABS,AZLOG */
/* ***END PROLOGUE  ZKSCL */
/*     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    ic = 0;
    nn = min(2,*n);
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1r = yr[i__];
	s1i = yi[i__];
	cyr[i__ - 1] = s1r;
	cyi[i__ - 1] = s1i;
	as = azabs_(&s1r, &s1i);
	acs = -(*zrr) + log(as);
	++(*nz);
	yr[i__] = zeror;
	yi[i__] = zeroi;
	if (acs < -(*elim)) {
	    goto L10;
	}
	azlog_(&s1r, &s1i, &csr, &csi, &idum);
	csr -= *zrr;
	csi -= *zri;
	str = exp(csr) / *tol;
	csr = str * cos(csi);
	csi = str * sin(csi);
	zuchk_(&csr, &csi, &nw, ascle, tol);
	if (nw != 0) {
	    goto L10;
	}
	yr[i__] = csr;
	yi[i__] = csi;
	ic = i__;
	--(*nz);
L10:
	;
    }
    if (*n == 1) {
	return 0;
    }
    if (ic > 1) {
	goto L20;
    }
    yr[1] = zeror;
    yi[1] = zeroi;
    *nz = 2;
L20:
    if (*n == 2) {
	return 0;
    }
    if (*nz == 0) {
	return 0;
    }
    fn = *fnu + 1.;
    ckr = fn * *rzr;
    cki = fn * *rzi;
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    helim = *elim * .5;
    elm = exp(-(*elim));
    celmr = elm;
    zdr = *zrr;
    zdi = *zri;

/*     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF */
/*     S2 GETS LARGER THAN EXP(ELIM/2) */

    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	kk = i__;
	csr = s2r;
	csi = s2i;
	s2r = ckr * csr - cki * csi + s1r;
	s2i = cki * csr + ckr * csi + s1i;
	s1r = csr;
	s1i = csi;
	ckr += *rzr;
	cki += *rzi;
	as = azabs_(&s2r, &s2i);
	alas = log(as);
	acs = -zdr + alas;
	++(*nz);
	yr[i__] = zeror;
	yi[i__] = zeroi;
	if (acs < -(*elim)) {
	    goto L25;
	}
	azlog_(&s2r, &s2i, &csr, &csi, &idum);
	csr -= zdr;
	csi -= zdi;
	str = exp(csr) / *tol;
	csr = str * cos(csi);
	csi = str * sin(csi);
	zuchk_(&csr, &csi, &nw, ascle, tol);
	if (nw != 0) {
	    goto L25;
	}
	yr[i__] = csr;
	yi[i__] = csi;
	--(*nz);
	if (ic == kk - 1) {
	    goto L40;
	}
	ic = kk;
	goto L30;
L25:
	if (alas < helim) {
	    goto L30;
	}
	zdr -= *elim;
	s1r *= celmr;
	s1i *= celmr;
	s2r *= celmr;
	s2i *= celmr;
L30:
	;
    }
    *nz = *n;
    if (ic == *n) {
	*nz = *n - 1;
    }
    goto L45;
L40:
    *nz = kk - 2;
L45:
    i__1 = *nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yr[i__] = zeror;
	yi[i__] = zeroi;
/* L50: */
    }
    return 0;
} /* zkscl_ */

