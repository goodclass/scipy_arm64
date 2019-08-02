/* dinvr.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dinvr_0_(int n__, integer *status, doublereal *x, 
	doublereal *fx, logical *qleft, logical *qhi, doublereal *zsmall, 
	doublereal *zbig, doublereal *zabsst, doublereal *zrelst, doublereal *
	zstpmu, doublereal *zabsto, doublereal *zrelto)
{
    /* Format strings */
    static char fmt_10[] = "";
    static char fmt_20[] = "";
    static char fmt_90[] = "";
    static char fmt_130[] = "";
    static char fmt_200[] = "";
    static char fmt_270[] = "";

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal yy, zx, zy, zz, big, xlb, xhi;
    static logical qok;
    static doublereal xub, xlo;
    static logical qup;
    static integer i99999;
    static doublereal fbig;
    static logical qbdd, qlim;
    static doublereal step;
    static logical qdum1, qdum2, qcond;
    static doublereal small;
    static logical qincr;
    static doublereal xsave;
    extern /* Subroutine */ int dzror_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, logical *), dstzr_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal fsmall, abstol, absstp, reltol, relstp, stpmul;

    /* Assigned format variables */
    static char *i99999_fmt;

/* ********************************************************************** */

/*     SUBROUTINE DINVR(STATUS, X, FX, QLEFT, QHI) */
/*          Double precision */
/*          bounds the zero of the function and invokes zror */
/*                    Reverse Communication */


/*                              Function */


/*     Bounds the    function  and  invokes  ZROR   to perform the   zero */
/*     finding.  STINVR  must  have   been  called  before this   routine */
/*     in order to set its parameters. */


/*                              Arguments */


/*     STATUS <--> At the beginning of a zero finding problem, STATUS */
/*                 should be set to 0 and INVR invoked.  (The value */
/*                 of parameters other than X will be ignored on this cal */

/*                 When INVR needs the function evaluated, it will set */
/*                 STATUS to 1 and return.  The value of the function */
/*                 should be set in FX and INVR again called without */
/*                 changing any of its other parameters. */

/*                 When INVR has finished without error, it will return */
/*                 with STATUS 0.  In that case X is approximately a root */
/*                 of F(X). */

/*                 If INVR cannot bound the function, it returns status */
/*                 -1 and sets QLEFT and QHI. */
/*                         INTEGER STATUS */

/*     X <-- The value of X at which F(X) is to be evaluated. */
/*                         DOUBLE PRECISION X */

/*     FX --> The value of F(X) calculated when INVR returns with */
/*            STATUS = 1. */
/*                         DOUBLE PRECISION FX */

/*     QLEFT <-- Defined only if QMFINV returns .FALSE.  In that */
/*          case it is .TRUE. If the stepping search terminated */
/*          unsuccessfully at SMALL.  If it is .FALSE. the search */
/*          terminated unsuccessfully at BIG. */
/*                    QLEFT is LOGICAL */

/*     QHI <-- Defined only if QMFINV returns .FALSE.  In that */
/*          case it is .TRUE. if F(X) .GT. Y at the termination */
/*          of the search and .FALSE. if F(X) .LT. Y at the */
/*          termination of the search. */
/*                    QHI is LOGICAL */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */
    switch(n__) {
	case 1: goto L_dstinv;
	}

    if (*status > 0) {
	goto L310;
    }
    qcond = ! (small <= *x && *x <= big);
    if (qcond) {
	s_stop(" SMALL, X, BIG not monotone in INVR", (ftnlen)35);
    }
    xsave = *x;

/*     See that SMALL and BIG bound the zero and set QINCR */

    *x = small;
/*     GET-FUNCTION-VALUE */
    i99999 = 0;
    i99999_fmt = fmt_10;
    goto L300;
L10:
    fsmall = *fx;
    *x = big;
/*     GET-FUNCTION-VALUE */
    i99999 = 1;
    i99999_fmt = fmt_20;
    goto L300;
L20:
    fbig = *fx;
    qincr = fbig > fsmall;
    if (! qincr) {
	goto L50;
    }
    if (fsmall <= 0.) {
	goto L30;
    }
    *status = -1;
    *qleft = TRUE_;
    *qhi = TRUE_;
    return 0;
L30:
    if (fbig >= 0.) {
	goto L40;
    }
    *status = -1;
    *qleft = FALSE_;
    *qhi = FALSE_;
    return 0;
L40:
    goto L80;
L50:
    if (fsmall >= 0.) {
	goto L60;
    }
    *status = -1;
    *qleft = TRUE_;
    *qhi = FALSE_;
    return 0;
L60:
    if (fbig <= 0.) {
	goto L70;
    }
    *status = -1;
    *qleft = FALSE_;
    *qhi = TRUE_;
    return 0;
L70:
L80:
    *x = xsave;
/* Computing MAX */
    d__1 = absstp, d__2 = relstp * abs(*x);
    step = max(d__1,d__2);
/*      YY = F(X) - Y */
/*     GET-FUNCTION-VALUE */
    i99999 = 2;
    i99999_fmt = fmt_90;
    goto L300;
L90:
    yy = *fx;
    if (! (yy == 0.)) {
	goto L100;
    }
    *status = 0;
    qok = TRUE_;
    return 0;
L100:
    qup = qincr && yy < 0. || ! qincr && yy > 0.;
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     HANDLE CASE IN WHICH WE MUST STEP HIGHER */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    if (! qup) {
	goto L170;
    }
    xlb = xsave;
/* Computing MIN */
    d__1 = xlb + step;
    xub = min(d__1,big);
    goto L120;
L110:
    if (qcond) {
	goto L150;
    }
/*      YY = F(XUB) - Y */
L120:
    *x = xub;
/*     GET-FUNCTION-VALUE */
    i99999 = 3;
    i99999_fmt = fmt_130;
    goto L300;
L130:
    yy = *fx;
    qbdd = qincr && yy >= 0. || ! qincr && yy <= 0.;
    qlim = xub >= big;
    qcond = qbdd || qlim;
    if (qcond) {
	goto L140;
    }
    step = stpmul * step;
    xlb = xub;
/* Computing MIN */
    d__1 = xlb + step;
    xub = min(d__1,big);
L140:
    goto L110;
L150:
    if (! (qlim && ! qbdd)) {
	goto L160;
    }
    *status = -1;
    *qleft = FALSE_;
    *qhi = ! qincr;
    *x = big;
    return 0;
L160:
    goto L240;
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     HANDLE CASE IN WHICH WE MUST STEP LOWER */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
L170:
    xub = xsave;
/* Computing MAX */
    d__1 = xub - step;
    xlb = max(d__1,small);
    goto L190;
L180:
    if (qcond) {
	goto L220;
    }
/*      YY = F(XLB) - Y */
L190:
    *x = xlb;
/*     GET-FUNCTION-VALUE */
    i99999 = 4;
    i99999_fmt = fmt_200;
    goto L300;
L200:
    yy = *fx;
    qbdd = qincr && yy <= 0. || ! qincr && yy >= 0.;
    qlim = xlb <= small;
    qcond = qbdd || qlim;
    if (qcond) {
	goto L210;
    }
    step = stpmul * step;
    xub = xlb;
/* Computing MAX */
    d__1 = xub - step;
    xlb = max(d__1,small);
L210:
    goto L180;
L220:
    if (! (qlim && ! qbdd)) {
	goto L230;
    }
    *status = -1;
    *qleft = TRUE_;
    *qhi = qincr;
    *x = small;
    return 0;
L230:
L240:
    dstzr_(&xlb, &xub, &abstol, &reltol);
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F. */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    *status = 0;
    goto L260;
L250:
    if (! (*status == 1)) {
	goto L290;
    }
L260:
    dzror_(status, x, fx, &xlo, &xhi, &qdum1, &qdum2);
    if (! (*status == 1)) {
	goto L280;
    }
/*     GET-FUNCTION-VALUE */
    i99999 = 5;
    i99999_fmt = fmt_270;
    goto L300;
L270:
L280:
    goto L250;
L290:
    *x = xlo;
    *status = 0;
    return 0;

L_dstinv:
/* ********************************************************************** */

/*      SUBROUTINE DSTINV( SMALL, BIG, ABSSTP, RELSTP, STPMUL, */
/*     +                   ABSTOL, RELTOL ) */
/*      Double Precision - SeT INverse finder - Reverse Communication */


/*                              Function */


/*     Concise Description - Given a monotone function F finds X */
/*     such that F(X) = Y.  Uses Reverse communication -- see invr. */
/*     This routine sets quantities needed by INVR. */

/*          More Precise Description of INVR - */

/*     F must be a monotone function, the results of QMFINV are */
/*     otherwise undefined.  QINCR must be .TRUE. if F is non- */
/*     decreasing and .FALSE. if F is non-increasing. */

/*     QMFINV will return .TRUE. if and only if F(SMALL) and */
/*     F(BIG) bracket Y, i. e., */
/*          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or */
/*          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL) */

/*     if QMFINV returns .TRUE., then the X returned satisfies */
/*     the following condition.  let */
/*               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X)) */
/*     then if QINCR is .TRUE., */
/*          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X)) */
/*     and if QINCR is .FALSE. */
/*          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X)) */


/*                              Arguments */


/*     SMALL --> The left endpoint of the interval to be */
/*          searched for a solution. */
/*                    SMALL is DOUBLE PRECISION */

/*     BIG --> The right endpoint of the interval to be */
/*          searched for a solution. */
/*                    BIG is DOUBLE PRECISION */

/*     ABSSTP, RELSTP --> The initial step size in the search */
/*          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm. */
/*                    ABSSTP is DOUBLE PRECISION */
/*                    RELSTP is DOUBLE PRECISION */

/*     STPMUL --> When a step doesn't bound the zero, the step */
/*                size is multiplied by STPMUL and another step */
/*                taken.  A popular value is 2.0 */
/*                    DOUBLE PRECISION STPMUL */

/*     ABSTOL, RELTOL --> Two numbers that determine the accuracy */
/*          of the solution.  See function for a precise definition. */
/*                    ABSTOL is DOUBLE PRECISION */
/*                    RELTOL is DOUBLE PRECISION */


/*                              Method */


/*     Compares F(X) with Y for the input value of X then uses QINCR */
/*     to determine whether to step left or right to bound the */
/*     desired x.  the initial step size is */
/*          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X. */
/*     Iteratively steps right or left until it bounds X. */
/*     At each step which doesn't bound X, the step size is doubled. */
/*     The routine is careful never to step beyond SMALL or BIG.  If */
/*     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE. */
/*     after setting QLEFT and QHI. */

/*     If X is successfully bounded then Algorithm R of the paper */
/*     'Two Efficient Algorithms with Guaranteed Convergence for */
/*     Finding a Zero of a Function' by J. C. P. Bus and */
/*     T. J. Dekker in ACM Transactions on Mathematical */
/*     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed */
/*     to find the zero of the function F(X)-Y. This is routine */
/*     QRZERO. */

/* ********************************************************************** */
/*     Reset all saved variables to known values */
    absstp = 0.;
    abstol = 0.;
    big = 0.;
    fbig = 0.;
    fsmall = 0.;
    relstp = 0.;
    reltol = 0.;
    small = 0.;
    step = 0.;
    stpmul = 0.;
    xhi = 0.;
    xlb = 0.;
    xlo = 0.;
    xsave = 0.;
    xub = 0.;
    yy = 0.;
    zx = 0.;
    zy = 0.;
    zz = 0.;
    i99999 = 0;
    qbdd = FALSE_;
    qcond = FALSE_;
    qdum1 = FALSE_;
    qdum2 = FALSE_;
    qincr = FALSE_;
    qlim = FALSE_;
    qok = FALSE_;
    qup = FALSE_;
/*     Set initial values */
    small = *zsmall;
    big = *zbig;
    absstp = *zabsst;
    relstp = *zrelst;
    stpmul = *zstpmu;
    abstol = *zabsto;
    reltol = *zrelto;
    return 0;
    s_stop("*** EXECUTION FLOWING INTO FLECS PROCEDURES ***", (ftnlen)47);
/*     TO GET-FUNCTION-VALUE */
L300:
    *status = 1;
    return 0;
L310:
    switch (i99999) {
	case 0: goto L10;
	case 1: goto L20;
	case 2: goto L90;
	case 3: goto L130;
	case 4: goto L200;
	case 5: goto L270;
    }
} /* dinvr_ */

/* Subroutine */ int dinvr_(integer *status, doublereal *x, doublereal *fx, 
	logical *qleft, logical *qhi)
{
    return dinvr_0_(0, status, x, fx, qleft, qhi, (doublereal *)0, (
	    doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal *)0,
	     (doublereal *)0, (doublereal *)0);
    }

/* Subroutine */ int dstinv_(doublereal *zsmall, doublereal *zbig, doublereal 
	*zabsst, doublereal *zrelst, doublereal *zstpmu, doublereal *zabsto, 
	doublereal *zrelto)
{
    return dinvr_0_(1, (integer *)0, (doublereal *)0, (doublereal *)0, (
	    logical *)0, (logical *)0, zsmall, zbig, zabsst, zrelst, zstpmu, 
	    zabsto, zrelto);
    }

