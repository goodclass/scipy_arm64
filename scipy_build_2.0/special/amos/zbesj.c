/* zbesj.f -- translated by f2c (version 20190311).
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

static integer c__4 = 4;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__5 = 5;
static integer c__14 = 14;
static integer c__9 = 9;
static integer c__1 = 1;

/* Subroutine */ int zbesj_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *cyr, doublereal *cyi, integer *
	nz, integer *ierr)
{
    /* Initialized data */

    static doublereal hpi = 1.57079632679489662;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, k1, k2;
    static doublereal aa, bb, fn;
    static integer nl;
    static doublereal az;
    static integer ir;
    static doublereal rl, dig, cii, arg, r1m5;
    static integer inu;
    static doublereal tol, sti, zni, str, znr, alim, elim, atol;
    static integer inuh;
    static doublereal fnul, rtol, ascle;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal csgni, csgnr;
    extern /* Subroutine */ int zbinu_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  ZBESJ */
/* ***DATE WRITTEN   830501   (YYMMDD) */
/* ***REVISION DATE  890801   (YYMMDD) */
/* ***CATEGORY NO.  B5K */
/* ***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT, */
/*             BESSEL FUNCTION OF FIRST KIND */
/* ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES */
/* ***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT */
/* ***DESCRIPTION */

/*                      ***A DOUBLE PRECISION ROUTINE*** */
/*         ON KODE=1, CBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX */
/*         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE */
/*         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE */
/*         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESJ RETURNS THE SCALED */
/*         FUNCTIONS */

/*         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z) */

/*         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND */
/*         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION */
/*         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS */
/*         (REF. 1). */

/*         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION */
/*           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI */
/*           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0D0 */
/*           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION */
/*                    KODE= 1  RETURNS */
/*                             CY(I)=J(FNU+I-1,Z), I=1,...,N */
/*                        = 2  RETURNS */
/*                             CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N */
/*           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1 */

/*         OUTPUT     CYR,CYI ARE DOUBLE PRECISION */
/*           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS */
/*                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE */
/*                    CY(I)=J(FNU+I-1,Z)  OR */
/*                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y))  I=1,...,N */
/*                    DEPENDING ON KODE, Y=AIMAG(Z). */
/*           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW, */
/*                    NZ= 0   , NORMAL RETURN */
/*                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET  ZERO DUE */
/*                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0), */
/*                              I = N-NZ+1,...,N */
/*           IERR   - ERROR FLAG */
/*                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED */
/*                    IERR=1, INPUT ERROR   - NO COMPUTATION */
/*                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z) */
/*                            TOO LARGE ON KODE=1 */
/*                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE */
/*                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT */
/*                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE */
/*                            ACCURACY */
/*                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA- */
/*                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI- */
/*                            CANCE BY ARGUMENT REDUCTION */
/*                    IERR=5, ERROR              - NO COMPUTATION, */
/*                            ALGORITHM TERMINATION CONDITION NOT MET */

/* ***LONG DESCRIPTION */

/*         THE COMPUTATION IS CARRIED OUT BY THE FORMULA */

/*         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0 */

/*         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0 */

/*         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION. */

/*         FOR NEGATIVE ORDERS,THE FORMULA */

/*              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU) */

/*         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE */
/*         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE */
/*         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A */
/*         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER, */
/*         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF */
/*         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY */
/*         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN */
/*         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE, */
/*         LARGE MEANS FNU.GT.CABS(Z). */

/*         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE- */
/*         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS */
/*         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. */
/*         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN */
/*         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG */
/*         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS */
/*         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION. */
/*         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS */
/*         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS */
/*         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE */
/*         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS */
/*         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3 */
/*         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION */
/*         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION */
/*         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN */
/*         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT */
/*         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS */
/*         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC. */
/*         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES. */

/*         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX */
/*         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT */
/*         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE- */
/*         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE */
/*         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))), */
/*         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF */
/*         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY */
/*         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN */
/*         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY */
/*         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER */
/*         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K, */
/*         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS */
/*         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER */
/*         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY */
/*         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER */
/*         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE */
/*         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES, */
/*         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P, */
/*         OR -PI/2+P. */

/* ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ */
/*                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF */
/*                 COMMERCE, 1955. */

/*               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT */
/*                 BY D. E. AMOS, SAND83-0083, MAY, 1983. */

/*               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT */
/*                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983 */

/*               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX */
/*                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85- */
/*                 1018, MAY, 1985 */

/*               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX */
/*                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS. */
/*                 MATH. SOFTWARE, 1986 */

/* ***ROUTINES CALLED  ZBINU,I1MACH,D1MACH */
/* ***END PROLOGUE  ZBESJ */

/*     COMPLEX CI,CSGN,CY,Z,ZN */
    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  ZBESJ */
    *ierr = 0;
    *nz = 0;
    if (*fnu < 0.) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*n < 1) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/*     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU. */
/* ----------------------------------------------------------------------- */
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    tol = max(d__1,1e-18);
    k1 = i1mach_(&c__15);
    k2 = i1mach_(&c__16);
    r1m5 = d1mach_(&c__5);
/* Computing MIN */
    i__1 = abs(k1), i__2 = abs(k2);
    k = min(i__1,i__2);
    elim = ((doublereal) ((real) k) * r1m5 - 3.) * 2.303;
    k1 = i1mach_(&c__14) - 1;
    aa = r1m5 * (doublereal) ((real) k1);
    dig = min(aa,18.);
    aa *= 2.303;
/* Computing MAX */
    d__1 = -aa;
    alim = elim + max(d__1,-41.45);
    rl = dig * 1.2 + 3.;
    fnul = (dig - 3.) * 6. + 10.;
/* ----------------------------------------------------------------------- */
/*     TEST FOR PROPER RANGE */
/* ----------------------------------------------------------------------- */
    az = azabs_(zr, zi);
    fn = *fnu + (doublereal) ((real) (*n - 1));
    aa = .5 / tol;
    bb = (doublereal) ((real) i1mach_(&c__9)) * .5;
    aa = min(aa,bb);
    if (az > aa) {
	goto L260;
    }
    if (fn > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    if (fn > aa) {
	*ierr = 3;
    }
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    cii = 1.;
    inu = (integer) ((real) (*fnu));
    inuh = inu / 2;
    ir = inu - (inuh << 1);
    arg = (*fnu - (doublereal) ((real) (inu - ir))) * hpi;
    csgnr = cos(arg);
    csgni = sin(arg);
    if (inuh % 2 == 0) {
	goto L40;
    }
    csgnr = -csgnr;
    csgni = -csgni;
L40:
/* ----------------------------------------------------------------------- */
/*     ZN IS IN THE RIGHT HALF PLANE */
/* ----------------------------------------------------------------------- */
    znr = *zi;
    zni = -(*zr);
    if (*zi >= 0.) {
	goto L50;
    }
    znr = -znr;
    zni = -zni;
    csgni = -csgni;
    cii = -cii;
L50:
    zbinu_(&znr, &zni, fnu, kode, n, &cyr[1], &cyi[1], nz, &rl, &fnul, &tol, &
	    elim, &alim);
    if (*nz < 0) {
	goto L130;
    }
    nl = *n - *nz;
    if (nl == 0) {
	return 0;
    }
    rtol = 1. / tol;
    ascle = d1mach_(&c__1) * rtol * 1e3;
    i__1 = nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       STR = CYR(I)*CSGNR - CYI(I)*CSGNI */
/*       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR */
/*       CYR(I) = STR */
	aa = cyr[i__];
	bb = cyi[i__];
	atol = 1.;
/* Computing MAX */
	d__1 = abs(aa), d__2 = abs(bb);
	if (max(d__1,d__2) > ascle) {
	    goto L55;
	}
	aa *= rtol;
	bb *= rtol;
	atol = tol;
L55:
	str = aa * csgnr - bb * csgni;
	sti = aa * csgni + bb * csgnr;
	cyr[i__] = str * atol;
	cyi[i__] = sti * atol;
	str = -csgni * cii;
	csgni = csgnr * cii;
	csgnr = str;
/* L60: */
    }
    return 0;
L130:
    if (*nz == -2) {
	goto L140;
    }
    *nz = 0;
    *ierr = 2;
    return 0;
L140:
    *nz = 0;
    *ierr = 5;
    return 0;
L260:
    *nz = 0;
    *ierr = 4;
    return 0;
} /* zbesj_ */

