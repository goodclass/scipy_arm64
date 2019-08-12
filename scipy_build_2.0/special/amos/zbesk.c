/* zbesk.f -- translated by f2c (version 20190311).
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
static integer c__2 = 2;

/* Subroutine */ int zbesk_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *cyr, doublereal *cyi, integer *
	nz, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static integer k, k1, k2;
    static doublereal aa, bb, fn, az;
    static integer nn;
    static doublereal rl;
    static integer mr, nw;
    static doublereal dig, arg, aln, r1m5, ufl;
    static integer nuf;
    static doublereal tol, alim, elim, fnul;
    extern doublereal azabs_(doublereal *, doublereal *);
    extern /* Subroutine */ int zacon_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *), zbknu_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *), zbunk_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int zuoik_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  ZBESK */
/* ***DATE WRITTEN   830501   (YYMMDD) */
/* ***REVISION DATE  890801   (YYMMDD) */
/* ***CATEGORY NO.  B5K */
/* ***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION OF THE SECOND KIND, */
/*             BESSEL FUNCTION OF THE THIRD KIND */
/* ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES */
/* ***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT */
/* ***DESCRIPTION */

/*                      ***A DOUBLE PRECISION ROUTINE*** */

/*         ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX */
/*         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE */
/*         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0) */
/*         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESK */
/*         RETURNS THE SCALED K FUNCTIONS, */

/*         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N, */

/*         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND */
/*         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND */
/*         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL */
/*         FUNCTIONS (REF. 1). */

/*         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION */
/*           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0), */
/*                    -PI.LT.ARG(Z).LE.PI */
/*           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0D0 */
/*           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1 */
/*           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION */
/*                    KODE= 1  RETURNS */
/*                             CY(I)=K(FNU+I-1,Z), I=1,...,N */
/*                        = 2  RETURNS */
/*                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N */

/*         OUTPUT     CYR,CYI ARE DOUBLE PRECISION */
/*           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS */
/*                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE */
/*                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR */
/*                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N */
/*                    DEPENDING ON KODE */
/*           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW. */
/*                    NZ= 0   , NORMAL RETURN */
/*                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE */
/*                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0), */
/*                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0 */
/*                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS */
/*                              IN THE SEQUENCE. */

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

/*         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS */
/*         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD */
/*         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT */
/*         HALF PLANE BY THE RELATION */

/*         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z) */
/*         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1 */

/*         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION. */

/*         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED */
/*         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS. */

/*         FOR NEGATIVE ORDERS, THE FORMULA */

/*                       K(-FNU,Z) = K(FNU,Z) */

/*         CAN BE USED. */

/*         CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS */
/*         AVAILABLE. */

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
/*                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983. */

/*               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX */
/*                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85- */
/*                 1018, MAY, 1985 */

/*               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX */
/*                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS. */
/*                 MATH. SOFTWARE, 1986 */

/* ***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,AZABS,I1MACH,D1MACH */
/* ***END PROLOGUE  ZBESK */

/*     COMPLEX CY,Z */
/* ***FIRST EXECUTABLE STATEMENT  ZBESK */
    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */
    *ierr = 0;
    *nz = 0;
    if (*zi == 0.f && *zr == 0.f) {
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
    nn = *n;
/* ----------------------------------------------------------------------- */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/*     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU */
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
    fnul = (dig - 3.) * 6. + 10.;
    rl = dig * 1.2 + 3.;
/* ----------------------------------------------------------------------------- */
/*     TEST FOR PROPER RANGE */
/* ----------------------------------------------------------------------- */
    az = azabs_(zr, zi);
    fn = *fnu + (doublereal) ((real) (nn - 1));
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
/*     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE */
/* ----------------------------------------------------------------------- */
/*     UFL = DEXP(-ELIM) */
    ufl = d1mach_(&c__1) * 1e3;
    if (az < ufl) {
	goto L180;
    }
    if (*fnu > fnul) {
	goto L80;
    }
    if (fn <= 1.) {
	goto L60;
    }
    if (fn > 2.) {
	goto L50;
    }
    if (az > tol) {
	goto L60;
    }
    arg = az * .5;
    aln = -fn * log(arg);
    if (aln > elim) {
	goto L180;
    }
    goto L60;
L50:
    zuoik_(zr, zi, fnu, kode, &c__2, &nn, &cyr[1], &cyi[1], &nuf, &tol, &elim,
	     &alim);
    if (nuf < 0) {
	goto L180;
    }
    *nz += nuf;
    nn -= nuf;
/* ----------------------------------------------------------------------- */
/*     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK */
/*     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I */
/* ----------------------------------------------------------------------- */
    if (nn == 0) {
	goto L100;
    }
L60:
    if (*zr < 0.) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0. */
/* ----------------------------------------------------------------------- */
    zbknu_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, &tol, &elim, &alim);
    if (nw < 0) {
	goto L200;
    }
    *nz = nw;
    return 0;
/* ----------------------------------------------------------------------- */
/*     LEFT HALF PLANE COMPUTATION */
/*     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2. */
/* ----------------------------------------------------------------------- */
L70:
    if (*nz != 0) {
	goto L180;
    }
    mr = 1;
    if (*zi < 0.) {
	mr = -1;
    }
    zacon_(zr, zi, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &rl, &fnul, &
	    tol, &elim, &alim);
    if (nw < 0) {
	goto L200;
    }
    *nz = nw;
    return 0;
/* ----------------------------------------------------------------------- */
/*     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL */
/* ----------------------------------------------------------------------- */
L80:
    mr = 0;
    if (*zr >= 0.) {
	goto L90;
    }
    mr = 1;
    if (*zi < 0.) {
	mr = -1;
    }
L90:
    zbunk_(zr, zi, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &tol, &elim, &
	    alim);
    if (nw < 0) {
	goto L200;
    }
    *nz += nw;
    return 0;
L100:
    if (*zr < 0.) {
	goto L180;
    }
    return 0;
L180:
    *nz = 0;
    *ierr = 2;
    return 0;
L200:
    if (nw == -1) {
	goto L180;
    }
    *nz = 0;
    *ierr = 5;
    return 0;
L260:
    *nz = 0;
    *ierr = 4;
    return 0;
} /* zbesk_ */

