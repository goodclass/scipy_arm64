/* zbesh.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zbesh_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *m, integer *n, doublereal *cyr, doublereal *
	cyi, integer *nz, integer *ierr)
{
    /* Initialized data */

    static doublereal hpi = 1.57079632679489662;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), d_sign(doublereal *, doublereal 
	    *), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, k1, k2;
    static doublereal aa, bb, fn;
    static integer mm;
    static doublereal az;
    static integer ir, nn;
    static doublereal rl;
    static integer mr, nw;
    static doublereal dig, arg, aln, fmm, r1m5, ufl, sgn;
    static integer nuf, inu;
    static doublereal tol, sti, zni, zti, str, znr, alim, elim, atol, rhpi;
    static integer inuh;
    static doublereal fnul, rtol, ascle;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal csgni;
    extern /* Subroutine */ int zacon_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);
    static doublereal csgnr;
    extern /* Subroutine */ int zbknu_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *), zbunk_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int zuoik_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  ZBESH */
/* ***DATE WRITTEN   830501   (YYMMDD) */
/* ***REVISION DATE  890801   (YYMMDD) */
/* ***CATEGORY NO.  B5K */
/* ***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT, */
/*             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS */
/* ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES */
/* ***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT */
/* ***DESCRIPTION */

/*                      ***A DOUBLE PRECISION ROUTINE*** */
/*         ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX */
/*         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1 */
/*         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX */
/*         Z.NE.CMPLX(0.0,0.0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. */
/*         ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS */

/*         CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1. */

/*         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND */
/*         LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE */
/*         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1). */

/*         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION */
/*           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0), */
/*                    -PT.LT.ARG(Z).LE.PI */
/*           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0D0 */
/*           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION */
/*                    KODE= 1  RETURNS */
/*                             CY(J)=H(M,FNU+J-1,Z),   J=1,...,N */
/*                        = 2  RETURNS */
/*                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M)) */
/*                                  J=1,...,N  ,  I**2=-1 */
/*           M      - KIND OF HANKEL FUNCTION, M=1 OR 2 */
/*           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N.GE.1 */

/*         OUTPUT     CYR,CYI ARE DOUBLE PRECISION */
/*           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS */
/*                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE */
/*                    CY(J)=H(M,FNU+J-1,Z)  OR */
/*                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N */
/*                    DEPENDING ON KODE, I**2=-1. */
/*           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW, */
/*                    NZ= 0   , NORMAL RETURN */
/*                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE */
/*                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0) */
/*                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR */
/*                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY */
/*                              HALF PLANES, NZ STATES ONLY THE NUMBER */
/*                              OF UNDERFLOWS. */
/*           IERR   - ERROR FLAG */
/*                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED */
/*                    IERR=1, INPUT ERROR   - NO COMPUTATION */
/*                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU TOO */
/*                            LARGE OR CABS(Z) TOO SMALL OR BOTH */
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

/*         THE COMPUTATION IS CARRIED OUT BY THE RELATION */

/*         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP)) */
/*             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1 */

/*         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE */
/*         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED */
/*         TO THE LEFT HALF PLANE BY THE RELATION */

/*         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z) */
/*         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1 */

/*         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION. */

/*         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z */
/*         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL */
/*         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING */
/*         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE */
/*         WHOLE Z PLANE FOR Z TO INFINITY. */

/*         FOR NEGATIVE ORDERS,THE FORMULAE */

/*               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I) */
/*               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I) */
/*                         I**2=-1 */

/*         CAN BE USED. */

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
/*         ROUNDOFF,1.0D-18) IS THE NOMINAL PRECISION AND 10**S REPRE- */
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

/* ***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,AZABS,I1MACH,D1MACH */
/* ***END PROLOGUE  ZBESH */

/*     COMPLEX CY,Z,ZN,ZT,CSGN */

    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  ZBESH */
    *ierr = 0;
    *nz = 0;
    if (*zr == 0. && *zi == 0.) {
	*ierr = 1;
    }
    if (*fnu < 0.) {
	*ierr = 1;
    }
    if (*m < 1 || *m > 2) {
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
    fn = *fnu + (doublereal) ((real) (nn - 1));
    mm = 3 - *m - *m;
    fmm = (doublereal) ((real) mm);
    znr = fmm * *zi;
    zni = -fmm * *zr;
/* ----------------------------------------------------------------------- */
/*     TEST FOR PROPER RANGE */
/* ----------------------------------------------------------------------- */
    az = azabs_(zr, zi);
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
    ufl = d1mach_(&c__1) * 1e3;
    if (az < ufl) {
	goto L230;
    }
    if (*fnu > fnul) {
	goto L90;
    }
    if (fn <= 1.) {
	goto L70;
    }
    if (fn > 2.) {
	goto L60;
    }
    if (az > tol) {
	goto L70;
    }
    arg = az * .5;
    aln = -fn * log(arg);
    if (aln > elim) {
	goto L230;
    }
    goto L70;
L60:
    zuoik_(&znr, &zni, fnu, kode, &c__2, &nn, &cyr[1], &cyi[1], &nuf, &tol, &
	    elim, &alim);
    if (nuf < 0) {
	goto L230;
    }
    *nz += nuf;
    nn -= nuf;
/* ----------------------------------------------------------------------- */
/*     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK */
/*     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I */
/* ----------------------------------------------------------------------- */
    if (nn == 0) {
	goto L140;
    }
L70:
    if (znr < 0. || znr == 0. && zni < 0. && *m == 2) {
	goto L80;
    }
/* ----------------------------------------------------------------------- */
/*     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR. */
/*     YN.GE.0. .OR. M=1) */
/* ----------------------------------------------------------------------- */
    zbknu_(&znr, &zni, fnu, kode, &nn, &cyr[1], &cyi[1], nz, &tol, &elim, &
	    alim);
    goto L110;
/* ----------------------------------------------------------------------- */
/*     LEFT HALF PLANE COMPUTATION */
/* ----------------------------------------------------------------------- */
L80:
    mr = -mm;
    zacon_(&znr, &zni, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &rl, &fnul,
	     &tol, &elim, &alim);
    if (nw < 0) {
	goto L240;
    }
    *nz = nw;
    goto L110;
L90:
/* ----------------------------------------------------------------------- */
/*     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL */
/* ----------------------------------------------------------------------- */
    mr = 0;
    if (znr >= 0. && (znr != 0. || zni >= 0. || *m != 2)) {
	goto L100;
    }
    mr = -mm;
    if (znr != 0. || zni >= 0.) {
	goto L100;
    }
    znr = -znr;
    zni = -zni;
L100:
    zbunk_(&znr, &zni, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &tol, &
	    elim, &alim);
    if (nw < 0) {
	goto L240;
    }
    *nz += nw;
L110:
/* ----------------------------------------------------------------------- */
/*     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT) */

/*     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2 */
/* ----------------------------------------------------------------------- */
    d__1 = -fmm;
    sgn = d_sign(&hpi, &d__1);
/* ----------------------------------------------------------------------- */
/*     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = (integer) ((real) (*fnu));
    inuh = inu / 2;
    ir = inu - (inuh << 1);
    arg = (*fnu - (doublereal) ((real) (inu - ir))) * sgn;
    rhpi = 1. / sgn;
/*     ZNI = RHPI*DCOS(ARG) */
/*     ZNR = -RHPI*DSIN(ARG) */
    csgni = rhpi * cos(arg);
    csgnr = -rhpi * sin(arg);
    if (inuh % 2 == 0) {
	goto L120;
    }
/*     ZNR = -ZNR */
/*     ZNI = -ZNI */
    csgnr = -csgnr;
    csgni = -csgni;
L120:
    zti = -fmm;
    rtol = 1. / tol;
    ascle = ufl * rtol;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       STR = CYR(I)*ZNR - CYI(I)*ZNI */
/*       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR */
/*       CYR(I) = STR */
/*       STR = -ZNI*ZTI */
/*       ZNI = ZNR*ZTI */
/*       ZNR = STR */
	aa = cyr[i__];
	bb = cyi[i__];
	atol = 1.;
/* Computing MAX */
	d__1 = abs(aa), d__2 = abs(bb);
	if (max(d__1,d__2) > ascle) {
	    goto L135;
	}
	aa *= rtol;
	bb *= rtol;
	atol = tol;
L135:
	str = aa * csgnr - bb * csgni;
	sti = aa * csgni + bb * csgnr;
	cyr[i__] = str * atol;
	cyi[i__] = sti * atol;
	str = -csgni * zti;
	csgni = csgnr * zti;
	csgnr = str;
/* L130: */
    }
    return 0;
L140:
    if (znr < 0.) {
	goto L230;
    }
    return 0;
L230:
    *nz = 0;
    *ierr = 2;
    return 0;
L240:
    if (nw == -1) {
	goto L230;
    }
    *nz = 0;
    *ierr = 5;
    return 0;
L260:
    *nz = 0;
    *ierr = 4;
    return 0;
} /* zbesh_ */

