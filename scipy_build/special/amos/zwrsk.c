/* zwrsk.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zwrsk_(doublereal *zrr, doublereal *zri, doublereal *fnu,
	 integer *kode, integer *n, doublereal *yr, doublereal *yi, integer *
	nz, doublereal *cwr, doublereal *cwi, doublereal *tol, doublereal *
	elim, doublereal *alim)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, nw;
    static doublereal c1i, c2i, c1r, c2r, act, acw, cti, ctr, pti, sti, ptr, 
	    str, ract, ascle;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal csclr, cinui, cinur;
    extern /* Subroutine */ int zbknu_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *), zrati_(doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal d1mach_(integer *);

/* ***BEGIN PROLOGUE  ZWRSK */
/* ***REFER TO  ZBESI,ZBESK */

/*     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY */
/*     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN */

/* ***ROUTINES CALLED  D1MACH,ZBKNU,ZRATI,AZABS */
/* ***END PROLOGUE  ZWRSK */
/*     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR */
/* ----------------------------------------------------------------------- */
/*     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS */
/*     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE */
/*     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --yi;
    --yr;
    --cwr;
    --cwi;

    /* Function Body */
    *nz = 0;
    zbknu_(zrr, zri, fnu, kode, &c__2, &cwr[1], &cwi[1], &nw, tol, elim, alim)
	    ;
    if (nw != 0) {
	goto L50;
    }
    zrati_(zrr, zri, fnu, n, &yr[1], &yi[1], tol);
/* ----------------------------------------------------------------------- */
/*     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z), */
/*     R(FNU+J-1,Z)=Y(J),  J=1,...,N */
/* ----------------------------------------------------------------------- */
    cinur = 1.;
    cinui = 0.;
    if (*kode == 1) {
	goto L10;
    }
    cinur = cos(*zri);
    cinui = sin(*zri);
L10:
/* ----------------------------------------------------------------------- */
/*     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH */
/*     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE */
/*     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT */
/*     THE RESULT IS ON SCALE. */
/* ----------------------------------------------------------------------- */
    acw = azabs_(&cwr[2], &cwi[2]);
    ascle = d1mach_(&c__1) * 1e3 / *tol;
    csclr = 1.;
    if (acw > ascle) {
	goto L20;
    }
    csclr = 1. / *tol;
    goto L30;
L20:
    ascle = 1. / ascle;
    if (acw < ascle) {
	goto L30;
    }
    csclr = *tol;
L30:
    c1r = cwr[1] * csclr;
    c1i = cwi[1] * csclr;
    c2r = cwr[2] * csclr;
    c2i = cwi[2] * csclr;
    str = yr[1];
    sti = yi[1];
/* ----------------------------------------------------------------------- */
/*     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0D0/CABS(CT) PREVENTS */
/*     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT) */
/* ----------------------------------------------------------------------- */
    ptr = str * c1r - sti * c1i;
    pti = str * c1i + sti * c1r;
    ptr += c2r;
    pti += c2i;
    ctr = *zrr * ptr - *zri * pti;
    cti = *zrr * pti + *zri * ptr;
    act = azabs_(&ctr, &cti);
    ract = 1. / act;
    ctr *= ract;
    cti = -cti * ract;
    ptr = cinur * ract;
    pti = cinui * ract;
    cinur = ptr * ctr - pti * cti;
    cinui = ptr * cti + pti * ctr;
    yr[1] = cinur * csclr;
    yi[1] = cinui * csclr;
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ptr = str * cinur - sti * cinui;
	cinui = str * cinui + sti * cinur;
	cinur = ptr;
	str = yr[i__];
	sti = yi[i__];
	yr[i__] = cinur * csclr;
	yi[i__] = cinui * csclr;
/* L40: */
    }
    return 0;
L50:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zwrsk_ */

