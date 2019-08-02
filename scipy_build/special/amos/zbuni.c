/* zbuni.f -- translated by f2c (version 20190311).
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

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int zbuni_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *yr, doublereal *yi, integer *
	nz, integer *nui, integer *nlast, doublereal *fnul, doublereal *tol, 
	doublereal *elim, doublereal *alim)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    static doublereal ax, ay;
    static integer nl, nw;
    static doublereal c1i, c1m, c1r, s1i, s2i, s1r, s2r, cyi[2], gnu, raz, 
	    cyr[2], sti, bry[3], rzi, str, rzr, dfnu, fnui;
    extern /* Subroutine */ int zuni1_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    , zuni2_(doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer iflag;
    static doublereal ascle;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal csclr, cscrr;
    static integer iform;
    extern doublereal d1mach_(integer *);

/* ***BEGIN PROLOGUE  ZBUNI */
/* ***REFER TO  ZBESI,ZBESK */

/*     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT. */
/*     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM */
/*     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING */
/*     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z) */
/*     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2 */

/* ***ROUTINES CALLED  ZUNI1,ZUNI2,AZABS,D1MACH */
/* ***END PROLOGUE  ZBUNI */
/*     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
    *nz = 0;
    ax = abs(*zr) * 1.7321;
    ay = abs(*zi);
    iform = 1;
    if (ay > ax) {
	iform = 2;
    }
    if (*nui == 0) {
	goto L60;
    }
    fnui = (doublereal) ((real) (*nui));
    dfnu = *fnu + (doublereal) ((real) (*n - 1));
    gnu = dfnu + fnui;
    if (iform == 2) {
	goto L10;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN */
/*     -PI/3.LE.ARG(Z).LE.PI/3 */
/* ----------------------------------------------------------------------- */
    zuni1_(zr, zi, &gnu, kode, &c__2, cyr, cyi, &nw, nlast, fnul, tol, elim, 
	    alim);
    goto L20;
L10:
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU */
/*     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I */
/*     AND HPI=PI/2 */
/* ----------------------------------------------------------------------- */
    zuni2_(zr, zi, &gnu, kode, &c__2, cyr, cyi, &nw, nlast, fnul, tol, elim, 
	    alim);
L20:
    if (nw < 0) {
	goto L50;
    }
    if (nw != 0) {
	goto L90;
    }
    str = azabs_(cyr, cyi);
/* ---------------------------------------------------------------------- */
/*     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED */
/* ---------------------------------------------------------------------- */
    bry[0] = d1mach_(&c__1) * 1e3 / *tol;
    bry[1] = 1. / bry[0];
    bry[2] = bry[1];
    iflag = 2;
    ascle = bry[1];
    csclr = 1.;
    if (str > bry[0]) {
	goto L21;
    }
    iflag = 1;
    ascle = bry[0];
    csclr = 1. / *tol;
    goto L25;
L21:
    if (str < bry[1]) {
	goto L25;
    }
    iflag = 3;
    ascle = bry[2];
    csclr = *tol;
L25:
    cscrr = 1. / csclr;
    s1r = cyr[1] * csclr;
    s1i = cyi[1] * csclr;
    s2r = cyr[0] * csclr;
    s2i = cyi[0] * csclr;
    raz = 1. / azabs_(zr, zi);
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    i__1 = *nui;
    for (i__ = 1; i__ <= i__1; ++i__) {
	str = s2r;
	sti = s2i;
	s2r = (dfnu + fnui) * (rzr * str - rzi * sti) + s1r;
	s2i = (dfnu + fnui) * (rzr * sti + rzi * str) + s1i;
	s1r = str;
	s1i = sti;
	fnui += -1.;
	if (iflag >= 3) {
	    goto L30;
	}
	str = s2r * cscrr;
	sti = s2i * cscrr;
	c1r = abs(str);
	c1i = abs(sti);
	c1m = max(c1r,c1i);
	if (c1m <= ascle) {
	    goto L30;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= cscrr;
	s1i *= cscrr;
	s2r = str;
	s2i = sti;
	csclr *= *tol;
	cscrr = 1. / csclr;
	s1r *= csclr;
	s1i *= csclr;
	s2r *= csclr;
	s2i *= csclr;
L30:
	;
    }
    yr[*n] = s2r * cscrr;
    yi[*n] = s2i * cscrr;
    if (*n == 1) {
	return 0;
    }
    nl = *n - 1;
    fnui = (doublereal) ((real) nl);
    k = nl;
    i__1 = nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
	str = s2r;
	sti = s2i;
	s2r = (*fnu + fnui) * (rzr * str - rzi * sti) + s1r;
	s2i = (*fnu + fnui) * (rzr * sti + rzi * str) + s1i;
	s1r = str;
	s1i = sti;
	str = s2r * cscrr;
	sti = s2i * cscrr;
	yr[k] = str;
	yi[k] = sti;
	fnui += -1.;
	--k;
	if (iflag >= 3) {
	    goto L40;
	}
	c1r = abs(str);
	c1i = abs(sti);
	c1m = max(c1r,c1i);
	if (c1m <= ascle) {
	    goto L40;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= cscrr;
	s1i *= cscrr;
	s2r = str;
	s2i = sti;
	csclr *= *tol;
	cscrr = 1. / csclr;
	s1r *= csclr;
	s1i *= csclr;
	s2r *= csclr;
	s2i *= csclr;
L40:
	;
    }
    return 0;
L50:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
L60:
    if (iform == 2) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN */
/*     -PI/3.LE.ARG(Z).LE.PI/3 */
/* ----------------------------------------------------------------------- */
    zuni1_(zr, zi, fnu, kode, n, &yr[1], &yi[1], &nw, nlast, fnul, tol, elim, 
	    alim);
    goto L80;
L70:
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU */
/*     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I */
/*     AND HPI=PI/2 */
/* ----------------------------------------------------------------------- */
    zuni2_(zr, zi, fnu, kode, n, &yr[1], &yi[1], &nw, nlast, fnul, tol, elim, 
	    alim);
L80:
    if (nw < 0) {
	goto L50;
    }
    *nz = nw;
    return 0;
L90:
    *nlast = *n;
    return 0;
} /* zbuni_ */

