/* zbinu.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zbinu_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *cyr, doublereal *cyi, integer *
	nz, doublereal *rl, doublereal *fnul, doublereal *tol, doublereal *
	elim, doublereal *alim)
{
    /* Initialized data */

    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal az;
    static integer nn, nw;
    static doublereal cwi[2], cwr[2];
    static integer nui, inw;
    static doublereal dfnu;
    extern doublereal azabs_(doublereal *, doublereal *);
    static integer nlast;
    extern /* Subroutine */ int zbuni_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), zseri_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *), zmlri_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *), zasyi_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *), zuoik_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *), zwrsk_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZBINU */
/* ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY */

/*     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE */

/* ***ROUTINES CALLED  AZABS,ZASYI,ZBUNI,ZMLRI,ZSERI,ZUOIK,ZWRSK */
/* ***END PROLOGUE  ZBINU */
    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */

    *nz = 0;
    az = azabs_(zr, zi);
    nn = *n;
    dfnu = *fnu + (doublereal) ((real) (*n - 1));
    if (az <= 2.) {
	goto L10;
    }
    if (az * az * .25 > dfnu + 1.) {
	goto L20;
    }
L10:
/* ----------------------------------------------------------------------- */
/*     POWER SERIES */
/* ----------------------------------------------------------------------- */
    zseri_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol, elim, alim);
    inw = abs(nw);
    *nz += inw;
    nn -= inw;
    if (nn == 0) {
	return 0;
    }
    if (nw >= 0) {
	goto L120;
    }
    dfnu = *fnu + (doublereal) ((real) (nn - 1));
L20:
    if (az < *rl) {
	goto L40;
    }
    if (dfnu <= 1.) {
	goto L30;
    }
    if (az + az < dfnu * dfnu) {
	goto L50;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z */
/* ----------------------------------------------------------------------- */
L30:
    zasyi_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, rl, tol, elim, alim)
	    ;
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L40:
    if (dfnu <= 1.) {
	goto L70;
    }
L50:
/* ----------------------------------------------------------------------- */
/*     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM */
/* ----------------------------------------------------------------------- */
    zuoik_(zr, zi, fnu, kode, &c__1, &nn, &cyr[1], &cyi[1], &nw, tol, elim, 
	    alim);
    if (nw < 0) {
	goto L130;
    }
    *nz += nw;
    nn -= nw;
    if (nn == 0) {
	return 0;
    }
    dfnu = *fnu + (doublereal) ((real) (nn - 1));
    if (dfnu > *fnul) {
	goto L110;
    }
    if (az > *fnul) {
	goto L110;
    }
L60:
    if (az > *rl) {
	goto L80;
    }
L70:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE SERIES */
/* ----------------------------------------------------------------------- */
    zmlri_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L80:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN */
/* ----------------------------------------------------------------------- */
    zuoik_(zr, zi, fnu, kode, &c__2, &c__2, cwr, cwi, &nw, tol, elim, alim);
    if (nw >= 0) {
	goto L100;
    }
    *nz = nn;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cyr[i__] = zeror;
	cyi[i__] = zeroi;
/* L90: */
    }
    return 0;
L100:
    if (nw > 0) {
	goto L130;
    }
    zwrsk_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, cwr, cwi, tol, elim,
	     alim);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L110:
/* ----------------------------------------------------------------------- */
/*     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD */
/* ----------------------------------------------------------------------- */
    nui = (integer) ((real) (*fnul - dfnu)) + 1;
    nui = max(nui,0);
    zbuni_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, &nui, &nlast, fnul, 
	    tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    *nz += nw;
    if (nlast == 0) {
	goto L120;
    }
    nn = nlast;
    goto L60;
L120:
    return 0;
L130:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zbinu_ */

