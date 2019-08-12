/* zbesy.f -- translated by f2c (version 20190311).
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
static integer c__4 = 4;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__5 = 5;

/* Subroutine */ int zbesy_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *cyr, doublereal *cyi, integer *
	nz, doublereal *cwrkr, doublereal *cwrki, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, k, k1, k2;
    static doublereal aa, bb, ey, c1i, c2i, c1r, c2r;
    static integer nz1, nz2;
    static doublereal exi;
    static real r1m5;
    static doublereal exr, sti, tay, tol, str, hcii, elim, atol, rtol, ascle;
    extern /* Subroutine */ int zbesh_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *);
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  ZBESY */
/* ***DATE WRITTEN   830501   (YYMMDD) */
/* ***REVISION DATE  890801   (YYMMDD) */
/* ***CATEGORY NO.  B5K */
/* ***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT, */
/*             BESSEL FUNCTION OF SECOND KIND */
/* ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES */
/* ***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT */
/* ***DESCRIPTION */

/*                      ***A DOUBLE PRECISION ROUTINE*** */

/*         ON KODE=1, CBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX */
/*         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE */
/*         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE */
/*         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESY RETURNS THE SCALED */
/*         FUNCTIONS */

/*         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z) */

/*         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND */
/*         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION */
/*         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS */
/*         (REF. 1). */

/*         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION */
/*           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0), */
/*                    -PI.LT.ARG(Z).LE.PI */
/*           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0 */
/*           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION */
/*                    KODE= 1  RETURNS */
/*                             CY(I)=Y(FNU+I-1,Z), I=1,...,N */
/*                        = 2  RETURNS */
/*                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N */
/*                             WHERE Y=AIMAG(Z) */
/*           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1 */
/*           CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT */
/*           CWRKI    AT LEAST N */

/*         OUTPUT     CYR,CYI ARE DOUBLE PRECISION */
/*           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS */
/*                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE */
/*                    CY(I)=Y(FNU+I-1,Z)  OR */
/*                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N */
/*                    DEPENDING ON KODE. */
/*           NZ     - NZ=0 , A NORMAL RETURN */
/*                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO */
/*                    UNDERFLOW (GENERALLY ON KODE=2) */
/*           IERR   - ERROR FLAG */
/*                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED */
/*                    IERR=1, INPUT ERROR   - NO COMPUTATION */
/*                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS */
/*                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH */
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

/*         Y(FNU,Z)=0.5*(H(1,FNU,Z)-H(2,FNU,Z))/I */

/*         WHERE I**2 = -1 AND THE HANKEL BESSEL FUNCTIONS H(1,FNU,Z) */
/*         AND H(2,FNU,Z) ARE CALCULATED IN CBESH. */

/*         FOR NEGATIVE ORDERS,THE FORMULA */

/*              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU) */

/*         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD */
/*         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE */
/*         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)* */
/*         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS */
/*         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A */
/*         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM */
/*         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, */
/*         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF */
/*         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z). */

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

/* ***ROUTINES CALLED  ZBESH,I1MACH,D1MACH */
/* ***END PROLOGUE  ZBESY */

/*     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV */
/* ***FIRST EXECUTABLE STATEMENT  ZBESY */
    /* Parameter adjustments */
    --cwrki;
    --cwrkr;
    --cyi;
    --cyr;

    /* Function Body */
    *ierr = 0;
    *nz = 0;
    if (*zr == 0. && *zi == 0.) {
	*ierr = 1;
    }
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
    hcii = .5;
    zbesh_(zr, zi, fnu, kode, &c__1, n, &cyr[1], &cyi[1], &nz1, ierr);
    if (*ierr != 0 && *ierr != 3) {
	goto L170;
    }
    zbesh_(zr, zi, fnu, kode, &c__2, n, &cwrkr[1], &cwrki[1], &nz2, ierr);
    if (*ierr != 0 && *ierr != 3) {
	goto L170;
    }
    *nz = min(nz1,nz2);
    if (*kode == 2) {
	goto L60;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	str = cwrkr[i__] - cyr[i__];
	sti = cwrki[i__] - cyi[i__];
	cyr[i__] = -sti * hcii;
	cyi[i__] = str * hcii;
/* L50: */
    }
    return 0;
L60:
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    tol = max(d__1,1e-18);
    k1 = i1mach_(&c__15);
    k2 = i1mach_(&c__16);
/* Computing MIN */
    i__1 = abs(k1), i__2 = abs(k2);
    k = min(i__1,i__2);
    r1m5 = d1mach_(&c__5);
/* ----------------------------------------------------------------------- */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT */
/* ----------------------------------------------------------------------- */
    elim = ((doublereal) ((real) k) * r1m5 - 3.) * 2.303;
    exr = cos(*zr);
    exi = sin(*zr);
    ey = 0.;
    tay = (d__1 = *zi + *zi, abs(d__1));
    if (tay < elim) {
	ey = exp(-tay);
    }
    if (*zi < 0.) {
	goto L90;
    }
    c1r = exr * ey;
    c1i = exi * ey;
    c2r = exr;
    c2i = -exi;
L70:
    *nz = 0;
    rtol = 1. / tol;
    ascle = d1mach_(&c__1) * rtol * 1e3;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       STR = C1R*CYR(I) - C1I*CYI(I) */
/*       STI = C1R*CYI(I) + C1I*CYR(I) */
/*       STR = -STR + C2R*CWRKR(I) - C2I*CWRKI(I) */
/*       STI = -STI + C2R*CWRKI(I) + C2I*CWRKR(I) */
/*       CYR(I) = -STI*HCII */
/*       CYI(I) = STR*HCII */
	aa = cwrkr[i__];
	bb = cwrki[i__];
	atol = 1.;
/* Computing MAX */
	d__1 = abs(aa), d__2 = abs(bb);
	if (max(d__1,d__2) > ascle) {
	    goto L75;
	}
	aa *= rtol;
	bb *= rtol;
	atol = tol;
L75:
	str = (aa * c2r - bb * c2i) * atol;
	sti = (aa * c2i + bb * c2r) * atol;
	aa = cyr[i__];
	bb = cyi[i__];
	atol = 1.;
/* Computing MAX */
	d__1 = abs(aa), d__2 = abs(bb);
	if (max(d__1,d__2) > ascle) {
	    goto L85;
	}
	aa *= rtol;
	bb *= rtol;
	atol = tol;
L85:
	str -= (aa * c1r - bb * c1i) * atol;
	sti -= (aa * c1i + bb * c1r) * atol;
	cyr[i__] = -sti * hcii;
	cyi[i__] = str * hcii;
	if (str == 0. && sti == 0. && ey == 0.) {
	    ++(*nz);
	}
/* L80: */
    }
    return 0;
L90:
    c1r = exr;
    c1i = exi;
    c2r = exr * ey;
    c2i = -exi * ey;
    goto L70;
L170:
    *nz = 0;
    return 0;
} /* zbesy_ */

