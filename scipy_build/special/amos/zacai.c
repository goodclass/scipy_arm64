/* zacai.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zacai_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *mr, integer *n, doublereal *yr, doublereal *
	yi, integer *nz, doublereal *rl, doublereal *tol, doublereal *elim, 
	doublereal *alim)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sin(doublereal), cos(
	    doublereal);

    /* Local variables */
    static doublereal az;
    static integer nn, nw;
    static doublereal yy, c1i, c2i, c1r, c2r, arg;
    static integer iuf;
    static doublereal cyi[2], fmr, sgn;
    static integer inu;
    static doublereal cyr[2], zni, znr, dfnu;
    extern /* Subroutine */ int zs1s2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal ascle;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal csgni, csgnr, cspni, cspnr;
    extern /* Subroutine */ int zbknu_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *), zseri_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    ;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int zmlri_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *), zasyi_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZACAI */
/* ***REFER TO  ZAIRY */

/*     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA */

/*         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN) */
/*                 MP=PI*MR*CMPLX(0.0,1.0) */

/*     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT */
/*     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1. */
/*     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND */
/*     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON */
/*     IS CALLED FROM ZAIRY. */

/* ***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,AZABS */
/* ***END PROLOGUE  ZACAI */
/*     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
    *nz = 0;
    znr = -(*zr);
    zni = -(*zi);
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
/*     POWER SERIES FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    zseri_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol, elim, alim);
    goto L40;
L20:
    if (az < *rl) {
	goto L30;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    zasyi_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, rl, tol, elim, 
	    alim);
    if (nw < 0) {
	goto L80;
    }
    goto L40;
L30:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    zmlri_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol);
    if (nw < 0) {
	goto L80;
    }
L40:
/* ----------------------------------------------------------------------- */
/*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION */
/* ----------------------------------------------------------------------- */
    zbknu_(&znr, &zni, fnu, kode, &c__1, cyr, cyi, &nw, tol, elim, alim);
    if (nw != 0) {
	goto L80;
    }
    fmr = (doublereal) ((real) (*mr));
    sgn = -d_sign(&pi, &fmr);
    csgnr = 0.;
    csgni = sgn;
    if (*kode == 1) {
	goto L50;
    }
    yy = -zni;
    csgnr = -csgni * sin(yy);
    csgni *= cos(yy);
L50:
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = (integer) ((real) (*fnu));
    arg = (*fnu - (doublereal) ((real) inu)) * sgn;
    cspnr = cos(arg);
    cspni = sin(arg);
    if (inu % 2 == 0) {
	goto L60;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L60:
    c1r = cyr[0];
    c1i = cyi[0];
    c2r = yr[1];
    c2i = yi[1];
    if (*kode == 1) {
	goto L70;
    }
    iuf = 0;
    ascle = d1mach_(&c__1) * 1e3 / *tol;
    zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
    *nz += nw;
L70:
    yr[1] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
    yi[1] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
    return 0;
L80:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zacai_ */

