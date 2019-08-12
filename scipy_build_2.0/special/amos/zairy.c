/* zairy.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int zairy_(doublereal *zr, doublereal *zi, integer *id, 
	integer *kode, doublereal *air, doublereal *aii, integer *nz, integer 
	*ierr)
{
    /* Initialized data */

    static doublereal tth = .666666666666666667;
    static doublereal c1 = .35502805388781724;
    static doublereal c2 = .258819403792806799;
    static doublereal coef = .183776298473930683;
    static doublereal zeror = 0.;
    static doublereal zeroi = 0.;
    static doublereal coner = 1.;
    static doublereal conei = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), pow_dd(doublereal *, doublereal *), sqrt(
	    doublereal);

    /* Local variables */
    static integer k;
    static doublereal d1, d2;
    static integer k1, k2;
    static doublereal aa, bb, ad, cc, ak, bk, ck, dk, az;
    static integer nn;
    static doublereal rl;
    static integer mr;
    static doublereal s1i, az3, s2i, s1r, s2r, z3i, z3r, dig, fid, cyi[1], 
	    r1m5, fnu, cyr[1], tol, sti, ptr, str, sfac, alim, elim, alaz, 
	    csqi, atrm, ztai, csqr, ztar, trm1i, trm2i, trm1r, trm2r;
    static integer iflag;
    extern /* Subroutine */ int zacai_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    extern doublereal azabs_(doublereal *, doublereal *);
    extern /* Subroutine */ int azexp_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), zbknu_(doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    extern /* Subroutine */ int azsqrt_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZAIRY */
/* ***DATE WRITTEN   830501   (YYMMDD) */
/* ***REVISION DATE  890801   (YYMMDD) */
/* ***CATEGORY NO.  B5K */
/* ***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD */
/* ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES */
/* ***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z */
/* ***DESCRIPTION */

/*                      ***A DOUBLE PRECISION ROUTINE*** */
/*         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR */
/*         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON */
/*         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)* */
/*         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN */
/*         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN */
/*         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z). */

/*         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN */
/*         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED */
/*         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS. */
/*         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF */
/*         MATHEMATICAL FUNCTIONS (REF. 1). */

/*         INPUT      ZR,ZI ARE DOUBLE PRECISION */
/*           ZR,ZI  - Z=CMPLX(ZR,ZI) */
/*           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1 */
/*           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION */
/*                    KODE= 1  RETURNS */
/*                             AI=AI(Z)                ON ID=0 OR */
/*                             AI=DAI(Z)/DZ            ON ID=1 */
/*                        = 2  RETURNS */
/*                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR */
/*                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE */
/*                             ZTA=(2/3)*Z*CSQRT(Z) */

/*         OUTPUT     AIR,AII ARE DOUBLE PRECISION */
/*           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND */
/*                    KODE */
/*           NZ     - UNDERFLOW INDICATOR */
/*                    NZ= 0   , NORMAL RETURN */
/*                    NZ= 1   , AI=CMPLX(0.0D0,0.0D0) DUE TO UNDERFLOW IN */
/*                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1 */
/*           IERR   - ERROR FLAG */
/*                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED */
/*                    IERR=1, INPUT ERROR   - NO COMPUTATION */
/*                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA) */
/*                            TOO LARGE ON KODE=1 */
/*                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED */
/*                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION */
/*                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY */
/*                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION */
/*                            COMPLETE LOSS OF ACCURACY BY ARGUMENT */
/*                            REDUCTION */
/*                    IERR=5, ERROR              - NO COMPUTATION, */
/*                            ALGORITHM TERMINATION CONDITION NOT MET */

/* ***LONG DESCRIPTION */

/*         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL */
/*         FUNCTIONS BY */

/*            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA) */
/*                           C=1.0/(PI*SQRT(3.0)) */
/*                            ZTA=(2/3)*Z**(3/2) */

/*         WITH THE POWER SERIES FOR CABS(Z).LE.1.0. */

/*         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE- */
/*         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES */
/*         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF */
/*         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR), */
/*         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR */
/*         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS */
/*         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION. */
/*         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN */
/*         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT */
/*         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE */
/*         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA */
/*         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, */
/*         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE */
/*         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE */
/*         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT- */
/*         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG- */
/*         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN */
/*         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN */
/*         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, */
/*         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE */
/*         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER */
/*         MACHINES. */

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
/*                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983 */

/*               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX */
/*                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85- */
/*                 1018, MAY, 1985 */

/*               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX */
/*                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS. */
/*                 MATH. SOFTWARE, 1986 */

/* ***ROUTINES CALLED  ZACAI,ZBKNU,AZEXP,AZSQRT,I1MACH,D1MACH */
/* ***END PROLOGUE  ZAIRY */
/*     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3 */
/* ***FIRST EXECUTABLE STATEMENT  ZAIRY */
    *ierr = 0;
    *nz = 0;
    if (*id < 0 || *id > 1) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    az = azabs_(zr, zi);
/* Computing MAX */
    d__1 = d1mach_(&c__4);
    tol = max(d__1,1e-18);
    fid = (doublereal) ((real) (*id));
    if (az > 1.) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     POWER SERIES FOR CABS(Z).LE.1. */
/* ----------------------------------------------------------------------- */
    s1r = coner;
    s1i = conei;
    s2r = coner;
    s2i = conei;
    if (az < tol) {
	goto L170;
    }
    aa = az * az;
    if (aa < tol / az) {
	goto L40;
    }
    trm1r = coner;
    trm1i = conei;
    trm2r = coner;
    trm2i = conei;
    atrm = 1.;
    str = *zr * *zr - *zi * *zi;
    sti = *zr * *zi + *zi * *zr;
    z3r = str * *zr - sti * *zi;
    z3i = str * *zi + sti * *zr;
    az3 = az * aa;
    ak = fid + 2.;
    bk = 3. - fid - fid;
    ck = 4. - fid;
    dk = fid + 3. + fid;
    d1 = ak * dk;
    d2 = bk * ck;
    ad = min(d1,d2);
    ak = fid * 9. + 24.;
    bk = 30. - fid * 9.;
    for (k = 1; k <= 25; ++k) {
	str = (trm1r * z3r - trm1i * z3i) / d1;
	trm1i = (trm1r * z3i + trm1i * z3r) / d1;
	trm1r = str;
	s1r += trm1r;
	s1i += trm1i;
	str = (trm2r * z3r - trm2i * z3i) / d2;
	trm2i = (trm2r * z3i + trm2i * z3r) / d2;
	trm2r = str;
	s2r += trm2r;
	s2i += trm2i;
	atrm = atrm * az3 / ad;
	d1 += ak;
	d2 += bk;
	ad = min(d1,d2);
	if (atrm < tol * ad) {
	    goto L40;
	}
	ak += 18.;
	bk += 18.;
/* L30: */
    }
L40:
    if (*id == 1) {
	goto L50;
    }
    *air = s1r * c1 - c2 * (*zr * s2r - *zi * s2i);
    *aii = s1i * c1 - c2 * (*zr * s2i + *zi * s2r);
    if (*kode == 1) {
	return 0;
    }
    azsqrt_(zr, zi, &str, &sti);
    ztar = tth * (*zr * str - *zi * sti);
    ztai = tth * (*zr * sti + *zi * str);
    azexp_(&ztar, &ztai, &str, &sti);
    ptr = *air * str - *aii * sti;
    *aii = *air * sti + *aii * str;
    *air = ptr;
    return 0;
L50:
    *air = -s2r * c2;
    *aii = -s2i * c2;
    if (az <= tol) {
	goto L60;
    }
    str = *zr * s1r - *zi * s1i;
    sti = *zr * s1i + *zi * s1r;
    cc = c1 / (fid + 1.);
    *air += cc * (str * *zr - sti * *zi);
    *aii += cc * (str * *zi + sti * *zr);
L60:
    if (*kode == 1) {
	return 0;
    }
    azsqrt_(zr, zi, &str, &sti);
    ztar = tth * (*zr * str - *zi * sti);
    ztai = tth * (*zr * sti + *zi * str);
    azexp_(&ztar, &ztai, &str, &sti);
    ptr = str * *air - sti * *aii;
    *aii = str * *aii + sti * *air;
    *air = ptr;
    return 0;
/* ----------------------------------------------------------------------- */
/*     CASE FOR CABS(Z).GT.1.0 */
/* ----------------------------------------------------------------------- */
L70:
    fnu = (fid + 1.) / 3.;
/* ----------------------------------------------------------------------- */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/* ----------------------------------------------------------------------- */
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
    alaz = log(az);
/* -------------------------------------------------------------------------- */
/*     TEST FOR PROPER RANGE */
/* ----------------------------------------------------------------------- */
    aa = .5 / tol;
    bb = (doublereal) ((real) i1mach_(&c__9)) * .5;
    aa = min(aa,bb);
    aa = pow_dd(&aa, &tth);
    if (az > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    azsqrt_(zr, zi, &csqr, &csqi);
    ztar = tth * (*zr * csqr - *zi * csqi);
    ztai = tth * (*zr * csqi + *zi * csqr);
/* ----------------------------------------------------------------------- */
/*     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL */
/* ----------------------------------------------------------------------- */
    iflag = 0;
    sfac = 1.;
    ak = ztai;
    if (*zr >= 0.) {
	goto L80;
    }
    bk = ztar;
    ck = -abs(bk);
    ztar = ck;
    ztai = ak;
L80:
    if (*zi != 0.) {
	goto L90;
    }
    if (*zr > 0.) {
	goto L90;
    }
    ztar = 0.;
    ztai = ak;
L90:
    aa = ztar;
    if (aa >= 0. && *zr > 0.) {
	goto L110;
    }
    if (*kode == 2) {
	goto L100;
    }
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST */
/* ----------------------------------------------------------------------- */
    if (aa > -alim) {
	goto L100;
    }
    aa = -aa + alaz * .25;
    iflag = 1;
    sfac = tol;
    if (aa > elim) {
	goto L270;
    }
L100:
/* ----------------------------------------------------------------------- */
/*     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2 */
/* ----------------------------------------------------------------------- */
    mr = 1;
    if (*zi < 0.) {
	mr = -1;
    }
    zacai_(&ztar, &ztai, &fnu, kode, &mr, &c__1, cyr, cyi, &nn, &rl, &tol, &
	    elim, &alim);
    if (nn < 0) {
	goto L280;
    }
    *nz += nn;
    goto L130;
L110:
    if (*kode == 2) {
	goto L120;
    }
/* ----------------------------------------------------------------------- */
/*     UNDERFLOW TEST */
/* ----------------------------------------------------------------------- */
    if (aa < alim) {
	goto L120;
    }
    aa = -aa - alaz * .25;
    iflag = 2;
    sfac = 1. / tol;
    if (aa < -elim) {
	goto L210;
    }
L120:
    zbknu_(&ztar, &ztai, &fnu, kode, &c__1, cyr, cyi, nz, &tol, &elim, &alim);
L130:
    s1r = cyr[0] * coef;
    s1i = cyi[0] * coef;
    if (iflag != 0) {
	goto L150;
    }
    if (*id == 1) {
	goto L140;
    }
    *air = csqr * s1r - csqi * s1i;
    *aii = csqr * s1i + csqi * s1r;
    return 0;
L140:
    *air = -(*zr * s1r - *zi * s1i);
    *aii = -(*zr * s1i + *zi * s1r);
    return 0;
L150:
    s1r *= sfac;
    s1i *= sfac;
    if (*id == 1) {
	goto L160;
    }
    str = s1r * csqr - s1i * csqi;
    s1i = s1r * csqi + s1i * csqr;
    s1r = str;
    *air = s1r / sfac;
    *aii = s1i / sfac;
    return 0;
L160:
    str = -(s1r * *zr - s1i * *zi);
    s1i = -(s1r * *zi + s1i * *zr);
    s1r = str;
    *air = s1r / sfac;
    *aii = s1i / sfac;
    return 0;
L170:
    aa = d1mach_(&c__1) * 1e3;
    s1r = zeror;
    s1i = zeroi;
    if (*id == 1) {
	goto L190;
    }
    if (az <= aa) {
	goto L180;
    }
    s1r = c2 * *zr;
    s1i = c2 * *zi;
L180:
    *air = c1 - s1r;
    *aii = -s1i;
    return 0;
L190:
    *air = -c2;
    *aii = 0.;
    aa = sqrt(aa);
    if (az <= aa) {
	goto L200;
    }
    s1r = (*zr * *zr - *zi * *zi) * .5;
    s1i = *zr * *zi;
L200:
    *air += c1 * s1r;
    *aii += c1 * s1i;
    return 0;
L210:
    *nz = 1;
    *air = zeror;
    *aii = zeroi;
    return 0;
L270:
    *nz = 0;
    *ierr = 2;
    return 0;
L280:
    if (nn == -1) {
	goto L270;
    }
    *nz = 0;
    *ierr = 5;
    return 0;
L260:
    *ierr = 4;
    *nz = 0;
    return 0;
} /* zairy_ */

