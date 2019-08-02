/* zvode.f -- translated by f2c (version 20190311).
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

/* Common Block Declarations */

union {
    struct {
	doublereal acnrm, ccmxj, conp, crate, drc, el[13], eta, etamax, h__, 
		hmin, hmxi, hnew, hrl1, hscal, prl1, rc, rl1, srur, tau[13], 
		tq[5], tn, uround;
	integer icf, init, ipup, jcur, jstart, jsv, kflag, kuth, l, lmax, lyh,
		 lewt, lacor, lsavf, lwm, liwm, locjs, maxord, meth, miter, 
		msbj, mxhnil, mxstep, n, newh, newq, nhnil, nq, nqnyh, nqwait,
		 nslj, nslp, nyh;
    } _1;
    struct {
	doublereal rvod1[50];
	integer ivod1[33];
    } _2;
} zvod01_;

#define zvod01_1 (zvod01_._1)
#define zvod01_2 (zvod01_._2)

union {
    struct {
	doublereal hu;
	integer ncfn, netf, nfe, nje, nlu, nni, nqu, nst;
    } _1;
    struct {
	doublereal rvod2[1];
	integer ivod2[8];
    } _2;
} zvod02_;

#define zvod02_1 (zvod02_._1)
#define zvod02_2 (zvod02_._2)

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__50 = 50;
static integer c__101 = 101;
static integer c__60 = 60;
static integer c__2 = 2;
static integer c__102 = 102;
static integer c__202 = 202;
static integer c__203 = 203;
static integer c__204 = 204;
static integer c__205 = 205;
static integer c__30 = 30;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__40 = 40;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__17 = 17;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__21 = 21;
static integer c__22 = 22;
static integer c__23 = 23;
static integer c__24 = 24;
static integer c__25 = 25;
static integer c__26 = 26;
static integer c__27 = 27;
static integer c__303 = 303;
static integer c__51 = 51;
static integer c__52 = 52;
static integer c_n1 = -1;
static doublereal c_b642 = 1.;

/* DECK ZVODE */
/* Subroutine */ int zvode_(S_fp f, integer *neq, doublecomplex *y, 
	doublereal *t, doublereal *tout, integer *itol, doublereal *rtol, 
	doublereal *atol, integer *itask, integer *istate, integer *iopt, 
	doublecomplex *zwork, integer *lzw, doublereal *rwork, integer *lrw, 
	integer *iwork, integer *liw, U_fp jac, integer *mf, real *rpar, 
	integer *ipar)
{
    /* Initialized data */

    static integer mord[2] = { 12,5 };
    static integer mxstp0 = 500;
    static integer mxhnl0 = 10;
    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal four = 4.;
    static doublereal pt2 = .2;
    static doublereal hun = 100.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer i_sign(integer *, integer *);
    double sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublereal h0;
    static integer ml;
    static doublereal rh;
    static integer mu;
    static doublereal tp;
    static integer lf0;
    static doublereal big;
    static integer mfa, jco, ier, kgo;
    static char msg[80];
    static doublereal hmx;
    static integer lenj;
    static logical ihit;
    static doublereal ewti, hmax, size;
    static integer lenp, mband, iflag;
    static doublereal atoli;
    static integer leniw, niter, lenwm, imxer;
    static doublereal tcrit, tolsf, rtoli;
    static integer lenrw;
    extern /* Subroutine */ int zvhin_(integer *, doublereal *, doublecomplex 
	    *, doublecomplex *, S_fp, real *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublereal *, integer *, 
	    integer *);
    static integer lenzw;
    static doublereal tnext;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dumach_(void);
    extern /* Subroutine */ int dzscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static integer nslast;
    extern /* Subroutine */ int xerrwd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int zvnlsd_();
    extern /* Subroutine */ int zewset_(integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, doublereal *), zvindy_(doublereal *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *);
    extern doublereal zvnorm_(integer *, doublecomplex *, doublereal *);
    extern /* Subroutine */ int zvstep_(doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    S_fp, U_fp, S_fp, U_fp, real *, integer *);

/* ----------------------------------------------------------------------- */
/* ZVODE: Variable-coefficient Ordinary Differential Equation solver, */
/* with fixed-leading-coefficient implementation. */
/* This version is in complex double precision. */

/* ZVODE solves the initial value problem for stiff or nonstiff */
/* systems of first order ODEs, */
/*     dy/dt = f(t,y) ,  or, in component form, */
/*     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ). */
/* Here the y vector is treated as complex. */
/* ZVODE is a package based on the EPISODE and EPISODEB packages, and */
/* on the ODEPACK user interface standard, with minor modifications. */

/* NOTE: When using ZVODE for a stiff system, it should only be used for */
/* the case in which the function f is analytic, that is, when each f(i) */
/* is an analytic function of each y(j).  Analyticity means that the */
/* partial derivative df(i)/dy(j) is a unique complex number, and this */
/* fact is critical in the way ZVODE solves the dense or banded linear */
/* systems that arise in the stiff case.  For a complex stiff ODE system */
/* in which f is not analytic, ZVODE is likely to have convergence */
/* failures, and for this problem one should instead use DVODE on the */
/* equivalent real system (in the real and imaginary parts of y). */
/* ----------------------------------------------------------------------- */
/* Authors: */
/*               Peter N. Brown and Alan C. Hindmarsh */
/*               Center for Applied Scientific Computing */
/*               Lawrence Livermore National Laboratory */
/*               Livermore, CA 94551 */
/* and */
/*               George D. Byrne (Prof. Emeritus) */
/*               Illinois Institute of Technology */
/*               Chicago, IL 60616 */
/* ----------------------------------------------------------------------- */
/* For references, see DVODE. */
/* ----------------------------------------------------------------------- */
/* Summary of usage. */

/* Communication between the user and the ZVODE package, for normal */
/* situations, is summarized here.  This summary describes only a subset */
/* of the full set of options available.  See the full description for */
/* details, including optional communication, nonstandard options, */
/* and instructions for special situations.  See also the example */
/* problem (with program and output) following this summary. */

/* A. First provide a subroutine of the form: */
/*           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR) */
/*           DOUBLE COMPLEX Y(NEQ), YDOT(NEQ) */
/*           DOUBLE PRECISION T */
/* which supplies the vector function f by loading YDOT(i) with f(i). */

/* B. Next determine (or guess) whether or not the problem is stiff. */
/* Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue */
/* whose real part is negative and large in magnitude, compared to the */
/* reciprocal of the t span of interest.  If the problem is nonstiff, */
/* use a method flag MF = 10.  If it is stiff, there are four standard */
/* choices for MF (21, 22, 24, 25), and ZVODE requires the Jacobian */
/* matrix in some form.  In these cases (MF .gt. 0), ZVODE will use a */
/* saved copy of the Jacobian matrix.  If this is undesirable because of */
/* storage limitations, set MF to the corresponding negative value */
/* (-21, -22, -24, -25).  (See full description of MF below.) */
/* The Jacobian matrix is regarded either as full (MF = 21 or 22), */
/* or banded (MF = 24 or 25).  In the banded case, ZVODE requires two */
/* half-bandwidth parameters ML and MU.  These are, respectively, the */
/* widths of the lower and upper parts of the band, excluding the main */
/* diagonal.  Thus the band consists of the locations (i,j) with */
/* i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1. */

/* C. If the problem is stiff, you are encouraged to supply the Jacobian */
/* directly (MF = 21 or 24), but if this is not feasible, ZVODE will */
/* compute it internally by difference quotients (MF = 22 or 25). */
/* If you are supplying the Jacobian, provide a subroutine of the form: */
/*           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR) */
/*           DOUBLE COMPLEX Y(NEQ), PD(NROWPD,NEQ) */
/*           DOUBLE PRECISION T */
/* which supplies df/dy by loading PD as follows: */
/*     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j), */
/* the partial derivative of f(i) with respect to y(j).  (Ignore the */
/* ML and MU arguments in this case.) */
/*     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with */
/* df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of */
/* PD from the top down. */
/*     In either case, only nonzero elements need be loaded. */

/* D. Write a main program which calls subroutine ZVODE once for */
/* each point at which answers are desired.  This should also provide */
/* for possible use of logical unit 6 for output of error messages */
/* by ZVODE.  On the first call to ZVODE, supply arguments as follows: */
/* F      = Name of subroutine for right-hand side vector f. */
/*          This name must be declared external in calling program. */
/* NEQ    = Number of first order ODEs. */
/* Y      = Double complex array of initial values, of length NEQ. */
/* T      = The initial value of the independent variable. */
/* TOUT   = First point where output is desired (.ne. T). */
/* ITOL   = 1 or 2 according as ATOL (below) is a scalar or array. */
/* RTOL   = Relative tolerance parameter (scalar). */
/* ATOL   = Absolute tolerance parameter (scalar or array). */
/*          The estimated local error in Y(i) will be controlled so as */
/*          to be roughly less (in magnitude) than */
/*             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or */
/*             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2. */
/*          Thus the local error test passes if, in each component, */
/*          either the absolute error is less than ATOL (or ATOL(i)), */
/*          or the relative error is less than RTOL. */
/*          Use RTOL = 0.0 for pure absolute error control, and */
/*          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error */
/*          control.  Caution: Actual (global) errors may exceed these */
/*          local tolerances, so choose them conservatively. */
/* ITASK  = 1 for normal computation of output values of Y at t = TOUT. */
/* ISTATE = Integer flag (input and output).  Set ISTATE = 1. */
/* IOPT   = 0 to indicate no optional input used. */
/* ZWORK  = Double precision complex work array of length at least: */
/*             15*NEQ                      for MF = 10, */
/*             8*NEQ + 2*NEQ**2            for MF = 21 or 22, */
/*             10*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25. */
/* LZW    = Declared length of ZWORK (in user's DIMENSION statement). */
/* RWORK  = Real work array of length at least 20 + NEQ. */
/* LRW    = Declared length of RWORK (in user's DIMENSION statement). */
/* IWORK  = Integer work array of length at least: */
/*             30        for MF = 10, */
/*             30 + NEQ  for MF = 21, 22, 24, or 25. */
/*          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower */
/*          and upper half-bandwidths ML,MU. */
/* LIW    = Declared length of IWORK (in user's DIMENSION statement). */
/* JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24). */
/*          If used, this name must be declared external in calling */
/*          program.  If not used, pass a dummy name. */
/* MF     = Method flag.  Standard values are: */
/*          10 for nonstiff (Adams) method, no Jacobian used. */
/*          21 for stiff (BDF) method, user-supplied full Jacobian. */
/*          22 for stiff method, internally generated full Jacobian. */
/*          24 for stiff method, user-supplied banded Jacobian. */
/*          25 for stiff method, internally generated banded Jacobian. */
/* RPAR   = user-defined real or complex array passed to F and JAC. */
/* IPAR   = user-defined integer array passed to F and JAC. */
/* Note that the main program must declare arrays Y, ZWORK, RWORK, IWORK, */
/* and possibly ATOL, RPAR, and IPAR.  RPAR may be declared REAL, DOUBLE, */
/* COMPLEX, or DOUBLE COMPLEX, depending on the user's needs. */

/* E. The output from the first call (or any call) is: */
/*      Y = Array of computed values of y(t) vector. */
/*      T = Corresponding value of independent variable (normally TOUT). */
/* ISTATE = 2  if ZVODE was successful, negative otherwise. */
/*          -1 means excess work done on this call. (Perhaps wrong MF.) */
/*          -2 means excess accuracy requested. (Tolerances too small.) */
/*          -3 means illegal input detected. (See printed message.) */
/*          -4 means repeated error test failures. (Check all input.) */
/*          -5 means repeated convergence failures. (Perhaps bad */
/*             Jacobian supplied or wrong choice of MF or tolerances.) */
/*          -6 means error weight became zero during problem. (Solution */
/*             component i vanished, and ATOL or ATOL(i) = 0.) */

/* F. To continue the integration after a successful return, simply */
/* reset TOUT and call ZVODE again.  No other parameters need be reset. */

/* ----------------------------------------------------------------------- */
/* EXAMPLE PROBLEM */

/* The program below uses ZVODE to solve the following system of 2 ODEs: */
/* dw/dt = -i*w*w*z, dz/dt = i*z; w(0) = 1/2.1, z(0) = 1; t = 0 to 2*pi. */
/* Solution: w = 1/(z + 1.1), z = exp(it).  As z traces the unit circle, */
/* w traces a circle of radius 10/2.1 with center at 11/2.1. */
/* For convenience, Main passes RPAR = (imaginary unit i) to FEX and JEX. */

/*     EXTERNAL FEX, JEX */
/*     DOUBLE COMPLEX Y(2), ZWORK(24), RPAR, WTRU, ERR */
/*     DOUBLE PRECISION ABERR, AEMAX, ATOL, RTOL, RWORK(22), T, TOUT */
/*     DIMENSION IWORK(32) */
/*     NEQ = 2 */
/*     Y(1) = 1.0D0/2.1D0 */
/*     Y(2) = 1.0D0 */
/*     T = 0.0D0 */
/*     DTOUT = 0.1570796326794896D0 */
/*     TOUT = DTOUT */
/*     ITOL = 1 */
/*     RTOL = 1.D-9 */
/*     ATOL = 1.D-8 */
/*     ITASK = 1 */
/*     ISTATE = 1 */
/*     IOPT = 0 */
/*     LZW = 24 */
/*     LRW = 22 */
/*     LIW = 32 */
/*     MF = 21 */
/*     RPAR = DCMPLX(0.0D0,1.0D0) */
/*     AEMAX = 0.0D0 */
/*     WRITE(6,10) */
/* 10  FORMAT('   t',11X,'w',26X,'z') */
/*     DO 40 IOUT = 1,40 */
/*       CALL ZVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT, */
/*    1             ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR) */
/*       WTRU = 1.0D0/DCMPLX(COS(T) + 1.1D0, SIN(T)) */
/*       ERR = Y(1) - WTRU */
/*       ABERR = ABS(DREAL(ERR)) + ABS(DIMAG(ERR)) */
/*       AEMAX = MAX(AEMAX,ABERR) */
/*       WRITE(6,20) T, DREAL(Y(1)),DIMAG(Y(1)), DREAL(Y(2)),DIMAG(Y(2)) */
/* 20    FORMAT(F9.5,2X,2F12.7,3X,2F12.7) */
/*       IF (ISTATE .LT. 0) THEN */
/*         WRITE(6,30) ISTATE */
/* 30      FORMAT(//'***** Error halt.  ISTATE =',I3) */
/*         STOP */
/*         ENDIF */
/* 40    TOUT = TOUT + DTOUT */
/*     WRITE(6,50) IWORK(11), IWORK(12), IWORK(13), IWORK(20), */
/*    1            IWORK(21), IWORK(22), IWORK(23), AEMAX */
/* 50  FORMAT(/' No. steps =',I4,'   No. f-s =',I5, */
/*    1        '   No. J-s =',I4,'   No. LU-s =',I4/ */
/*    2        ' No. nonlinear iterations =',I4/ */
/*    3        ' No. nonlinear convergence failures =',I4/ */
/*    4        ' No. error test failures =',I4/ */
/*    5        ' Max. abs. error in w =',D10.2) */
/*     STOP */
/*     END */

/*     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR) */
/*     DOUBLE COMPLEX Y(NEQ), YDOT(NEQ), RPAR */
/*     DOUBLE PRECISION T */
/*     YDOT(1) = -RPAR*Y(1)*Y(1)*Y(2) */
/*     YDOT(2) = RPAR*Y(2) */
/*     RETURN */
/*     END */

/*     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR) */
/*     DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR */
/*     DOUBLE PRECISION T */
/*     PD(1,1) = -2.0D0*RPAR*Y(1)*Y(2) */
/*     PD(1,2) = -RPAR*Y(1)*Y(1) */
/*     PD(2,2) = RPAR */
/*     RETURN */
/*     END */

/* The output of this example program is as follows: */

/*   t           w                          z */
/*  0.15708     0.4763242  -0.0356919      0.9876884   0.1564345 */
/*  0.31416     0.4767322  -0.0718256      0.9510565   0.3090170 */
/*  0.47124     0.4774351  -0.1088651      0.8910065   0.4539906 */
/*  0.62832     0.4784699  -0.1473206      0.8090170   0.5877853 */
/*  0.78540     0.4798943  -0.1877789      0.7071067   0.7071069 */
/*  0.94248     0.4817938  -0.2309414      0.5877852   0.8090171 */
/*  1.09956     0.4842934  -0.2776778      0.4539904   0.8910066 */
/*  1.25664     0.4875766  -0.3291039      0.3090169   0.9510566 */
/*  1.41372     0.4919177  -0.3866987      0.1564343   0.9876884 */
/*  1.57080     0.4977376  -0.4524889     -0.0000001   1.0000000 */
/*  1.72788     0.5057044  -0.5293524     -0.1564346   0.9876883 */
/*  1.88496     0.5169274  -0.6215400     -0.3090171   0.9510565 */
/*  2.04204     0.5333540  -0.7356275     -0.4539906   0.8910065 */
/*  2.19911     0.5586542  -0.8823669     -0.5877854   0.8090169 */
/*  2.35619     0.6004188  -1.0806013     -0.7071069   0.7071067 */
/*  2.51327     0.6764486  -1.3664281     -0.8090171   0.5877851 */
/*  2.67035     0.8366909  -1.8175245     -0.8910066   0.4539904 */
/*  2.82743     1.2657121  -2.6260146     -0.9510566   0.3090168 */
/*  2.98451     3.0284506  -4.2182180     -0.9876884   0.1564343 */
/*  3.14159    10.0000699   0.0000663     -1.0000000  -0.0000002 */
/*  3.29867     3.0284170   4.2182053     -0.9876883  -0.1564346 */
/*  3.45575     1.2657041   2.6260067     -0.9510565  -0.3090172 */
/*  3.61283     0.8366878   1.8175205     -0.8910064  -0.4539907 */
/*  3.76991     0.6764469   1.3664259     -0.8090169  -0.5877854 */
/*  3.92699     0.6004178   1.0806000     -0.7071066  -0.7071069 */
/*  4.08407     0.5586535   0.8823662     -0.5877851  -0.8090171 */
/*  4.24115     0.5333535   0.7356271     -0.4539903  -0.8910066 */
/*  4.39823     0.5169271   0.6215398     -0.3090168  -0.9510566 */
/*  4.55531     0.5057041   0.5293523     -0.1564343  -0.9876884 */
/*  4.71239     0.4977374   0.4524890      0.0000002  -1.0000000 */
/*  4.86947     0.4919176   0.3866988      0.1564347  -0.9876883 */
/*  5.02655     0.4875765   0.3291040      0.3090172  -0.9510564 */
/*  5.18363     0.4842934   0.2776780      0.4539907  -0.8910064 */
/*  5.34071     0.4817939   0.2309415      0.5877854  -0.8090169 */
/*  5.49779     0.4798944   0.1877791      0.7071069  -0.7071066 */
/*  5.65487     0.4784700   0.1473208      0.8090171  -0.5877850 */
/*  5.81195     0.4774352   0.1088652      0.8910066  -0.4539903 */
/*  5.96903     0.4767324   0.0718257      0.9510566  -0.3090168 */
/*  6.12611     0.4763244   0.0356920      0.9876884  -0.1564342 */
/*  6.28319     0.4761907   0.0000000      1.0000000   0.0000003 */

/* No. steps = 542   No. f-s =  610   No. J-s =  10   No. LU-s =  47 */
/* No. nonlinear iterations = 607 */
/* No. nonlinear convergence failures =   0 */
/* No. error test failures =  13 */
/* Max. abs. error in w =  0.13E-03 */

/* ----------------------------------------------------------------------- */
/* Full description of user interface to ZVODE. */

/* The user interface to ZVODE consists of the following parts. */

/* i.   The call sequence to subroutine ZVODE, which is a driver */
/*      routine for the solver.  This includes descriptions of both */
/*      the call sequence arguments and of user-supplied routines. */
/*      Following these descriptions is */
/*        * a description of optional input available through the */
/*          call sequence, */
/*        * a description of optional output (in the work arrays), and */
/*        * instructions for interrupting and restarting a solution. */

/* ii.  Descriptions of other routines in the ZVODE package that may be */
/*      (optionally) called by the user.  These provide the ability to */
/*      alter error message handling, save and restore the internal */
/*      COMMON, and obtain specified derivatives of the solution y(t). */

/* iii. Descriptions of COMMON blocks to be declared in overlay */
/*      or similar environments. */

/* iv.  Description of two routines in the ZVODE package, either of */
/*      which the user may replace with his own version, if desired. */
/*      these relate to the measurement of errors. */

/* ----------------------------------------------------------------------- */
/* Part i.  Call Sequence. */

/* The call sequence parameters used for input only are */
/*     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF, */
/* and those used for both input and output are */
/*     Y, T, ISTATE. */
/* The work arrays ZWORK, RWORK, and IWORK are also used for conditional */
/* and optional input and optional output.  (The term output here refers */
/* to the return from subroutine ZVODE to the user's calling program.) */

/* The legality of input parameters will be thoroughly checked on the */
/* initial call for the problem, but not checked thereafter unless a */
/* change in input parameters is flagged by ISTATE = 3 in the input. */

/* The descriptions of the call arguments are as follows. */

/* F      = The name of the user-supplied subroutine defining the */
/*          ODE system.  The system must be put in the first-order */
/*          form dy/dt = f(t,y), where f is a vector-valued function */
/*          of the scalar t and the vector y.  Subroutine F is to */
/*          compute the function f.  It is to have the form */
/*               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR) */
/*               DOUBLE COMPLEX Y(NEQ), YDOT(NEQ) */
/*               DOUBLE PRECISION T */
/*          where NEQ, T, and Y are input, and the array YDOT = f(t,y) */
/*          is output.  Y and YDOT are double complex arrays of length */
/*          NEQ.  Subroutine F should not alter Y(1),...,Y(NEQ). */
/*          F must be declared EXTERNAL in the calling program. */

/*          Subroutine F may access user-defined real/complex and */
/*          integer work arrays RPAR and IPAR, which are to be */
/*          dimensioned in the calling program. */

/*          If quantities computed in the F routine are needed */
/*          externally to ZVODE, an extra call to F should be made */
/*          for this purpose, for consistent and accurate results. */
/*          If only the derivative dy/dt is needed, use ZVINDY instead. */

/* NEQ    = The size of the ODE system (number of first order */
/*          ordinary differential equations).  Used only for input. */
/*          NEQ may not be increased during the problem, but */
/*          can be decreased (with ISTATE = 3 in the input). */

/* Y      = A double precision complex array for the vector of dependent */
/*          variables, of length NEQ or more.  Used for both input and */
/*          output on the first call (ISTATE = 1), and only for output */
/*          on other calls.  On the first call, Y must contain the */
/*          vector of initial values.  In the output, Y contains the */
/*          computed solution evaluated at T.  If desired, the Y array */
/*          may be used for other purposes between calls to the solver. */

/*          This array is passed as the Y argument in all calls to */
/*          F and JAC. */

/* T      = The independent variable.  In the input, T is used only on */
/*          the first call, as the initial point of the integration. */
/*          In the output, after each call, T is the value at which a */
/*          computed solution Y is evaluated (usually the same as TOUT). */
/*          On an error return, T is the farthest point reached. */

/* TOUT   = The next value of t at which a computed solution is desired. */
/*          Used only for input. */

/*          When starting the problem (ISTATE = 1), TOUT may be equal */
/*          to T for one call, then should .ne. T for the next call. */
/*          For the initial T, an input value of TOUT .ne. T is used */
/*          in order to determine the direction of the integration */
/*          (i.e. the algebraic sign of the step sizes) and the rough */
/*          scale of the problem.  Integration in either direction */
/*          (forward or backward in t) is permitted. */

/*          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after */
/*          the first call (i.e. the first call with TOUT .ne. T). */
/*          Otherwise, TOUT is required on every call. */

/*          If ITASK = 1, 3, or 4, the values of TOUT need not be */
/*          monotone, but a value of TOUT which backs up is limited */
/*          to the current internal t interval, whose endpoints are */
/*          TCUR - HU and TCUR.  (See optional output, below, for */
/*          TCUR and HU.) */

/* ITOL   = An indicator for the type of error control.  See */
/*          description below under ATOL.  Used only for input. */

/* RTOL   = A relative error tolerance parameter, either a scalar or */
/*          an array of length NEQ.  See description below under ATOL. */
/*          Input only. */

/* ATOL   = An absolute error tolerance parameter, either a scalar or */
/*          an array of length NEQ.  Input only. */

/*          The input parameters ITOL, RTOL, and ATOL determine */
/*          the error control performed by the solver.  The solver will */
/*          control the vector e = (e(i)) of estimated local errors */
/*          in Y, according to an inequality of the form */
/*                      rms-norm of ( e(i)/EWT(i) )   .le.   1, */
/*          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i), */
/*          and the rms-norm (root-mean-square norm) here is */
/*          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i)) */
/*          is a vector of weights which must always be positive, and */
/*          the values of RTOL and ATOL should all be non-negative. */
/*          The following table gives the types (scalar/array) of */
/*          RTOL and ATOL, and the corresponding form of EWT(i). */

/*             ITOL    RTOL       ATOL          EWT(i) */
/*              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL */
/*              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i) */
/*              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL */
/*              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i) */

/*          When either of these parameters is a scalar, it need not */
/*          be dimensioned in the user's calling program. */

/*          If none of the above choices (with ITOL, RTOL, and ATOL */
/*          fixed throughout the problem) is suitable, more general */
/*          error controls can be obtained by substituting */
/*          user-supplied routines for the setting of EWT and/or for */
/*          the norm calculation.  See Part iv below. */

/*          If global errors are to be estimated by making a repeated */
/*          run on the same problem with smaller tolerances, then all */
/*          components of RTOL and ATOL (i.e. of EWT) should be scaled */
/*          down uniformly. */

/* ITASK  = An index specifying the task to be performed. */
/*          Input only.  ITASK has the following values and meanings. */
/*          1  means normal computation of output values of y(t) at */
/*             t = TOUT (by overshooting and interpolating). */
/*          2  means take one step only and return. */
/*          3  means stop at the first internal mesh point at or */
/*             beyond t = TOUT and return. */
/*          4  means normal computation of output values of y(t) at */
/*             t = TOUT but without overshooting t = TCRIT. */
/*             TCRIT must be input as RWORK(1).  TCRIT may be equal to */
/*             or beyond TOUT, but not behind it in the direction of */
/*             integration.  This option is useful if the problem */
/*             has a singularity at or beyond t = TCRIT. */
/*          5  means take one step, without passing TCRIT, and return. */
/*             TCRIT must be input as RWORK(1). */

/*          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT */
/*          (within roundoff), it will return T = TCRIT (exactly) to */
/*          indicate this (unless ITASK = 4 and TOUT comes before TCRIT, */
/*          in which case answers at T = TOUT are returned first). */

/* ISTATE = an index used for input and output to specify the */
/*          the state of the calculation. */

/*          In the input, the values of ISTATE are as follows. */
/*          1  means this is the first call for the problem */
/*             (initializations will be done).  See note below. */
/*          2  means this is not the first call, and the calculation */
/*             is to continue normally, with no change in any input */
/*             parameters except possibly TOUT and ITASK. */
/*             (If ITOL, RTOL, and/or ATOL are changed between calls */
/*             with ISTATE = 2, the new values will be used but not */
/*             tested for legality.) */
/*          3  means this is not the first call, and the */
/*             calculation is to continue normally, but with */
/*             a change in input parameters other than */
/*             TOUT and ITASK.  Changes are allowed in */
/*             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU, */
/*             and any of the optional input except H0. */
/*             (See IWORK description for ML and MU.) */
/*          Note:  A preliminary call with TOUT = T is not counted */
/*          as a first call here, as no initialization or checking of */
/*          input is done.  (Such a call is sometimes useful to include */
/*          the initial conditions in the output.) */
/*          Thus the first call for which TOUT .ne. T requires */
/*          ISTATE = 1 in the input. */

/*          In the output, ISTATE has the following values and meanings. */
/*           1  means nothing was done, as TOUT was equal to T with */
/*              ISTATE = 1 in the input. */
/*           2  means the integration was performed successfully. */
/*          -1  means an excessive amount of work (more than MXSTEP */
/*              steps) was done on this call, before completing the */
/*              requested task, but the integration was otherwise */
/*              successful as far as T.  (MXSTEP is an optional input */
/*              and is normally 500.)  To continue, the user may */
/*              simply reset ISTATE to a value .gt. 1 and call again. */
/*              (The excess work step counter will be reset to 0.) */
/*              In addition, the user may increase MXSTEP to avoid */
/*              this error return.  (See optional input below.) */
/*          -2  means too much accuracy was requested for the precision */
/*              of the machine being used.  This was detected before */
/*              completing the requested task, but the integration */
/*              was successful as far as T.  To continue, the tolerance */
/*              parameters must be reset, and ISTATE must be set */
/*              to 3.  The optional output TOLSF may be used for this */
/*              purpose.  (Note: If this condition is detected before */
/*              taking any steps, then an illegal input return */
/*              (ISTATE = -3) occurs instead.) */
/*          -3  means illegal input was detected, before taking any */
/*              integration steps.  See written message for details. */
/*              Note:  If the solver detects an infinite loop of calls */
/*              to the solver with illegal input, it will cause */
/*              the run to stop. */
/*          -4  means there were repeated error test failures on */
/*              one attempted step, before completing the requested */
/*              task, but the integration was successful as far as T. */
/*              The problem may have a singularity, or the input */
/*              may be inappropriate. */
/*          -5  means there were repeated convergence test failures on */
/*              one attempted step, before completing the requested */
/*              task, but the integration was successful as far as T. */
/*              This may be caused by an inaccurate Jacobian matrix, */
/*              if one is being used. */
/*          -6  means EWT(i) became zero for some i during the */
/*              integration.  Pure relative error control (ATOL(i)=0.0) */
/*              was requested on a variable which has now vanished. */
/*              The integration was successful as far as T. */

/*          Note:  Since the normal output value of ISTATE is 2, */
/*          it does not need to be reset for normal continuation. */
/*          Also, since a negative input value of ISTATE will be */
/*          regarded as illegal, a negative output value requires the */
/*          user to change it, and possibly other input, before */
/*          calling the solver again. */

/* IOPT   = An integer flag to specify whether or not any optional */
/*          input is being used on this call.  Input only. */
/*          The optional input is listed separately below. */
/*          IOPT = 0 means no optional input is being used. */
/*                   Default values will be used in all cases. */
/*          IOPT = 1 means optional input is being used. */

/* ZWORK  = A double precision complex working array. */
/*          The length of ZWORK must be at least */
/*             NYH*(MAXORD + 1) + 2*NEQ + LWM    where */
/*          NYH    = the initial value of NEQ, */
/*          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a */
/*                   smaller value is given as an optional input), */
/*          LWM = length of work space for matrix-related data: */
/*          LWM = 0             if MITER = 0, */
/*          LWM = 2*NEQ**2      if MITER = 1 or 2, and MF.gt.0, */
/*          LWM = NEQ**2        if MITER = 1 or 2, and MF.lt.0, */
/*          LWM = NEQ           if MITER = 3, */
/*          LWM = (3*ML+2*MU+2)*NEQ     if MITER = 4 or 5, and MF.gt.0, */
/*          LWM = (2*ML+MU+1)*NEQ       if MITER = 4 or 5, and MF.lt.0. */
/*          (See the MF description for METH and MITER.) */
/*          Thus if MAXORD has its default value and NEQ is constant, */
/*          this length is: */
/*             15*NEQ                    for MF = 10, */
/*             15*NEQ + 2*NEQ**2         for MF = 11 or 12, */
/*             15*NEQ + NEQ**2           for MF = -11 or -12, */
/*             16*NEQ                    for MF = 13, */
/*             17*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15, */
/*             16*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15, */
/*              8*NEQ                    for MF = 20, */
/*              8*NEQ + 2*NEQ**2         for MF = 21 or 22, */
/*              8*NEQ + NEQ**2           for MF = -21 or -22, */
/*              9*NEQ                    for MF = 23, */
/*             10*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25. */
/*              9*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25. */

/* LZW    = The length of the array ZWORK, as declared by the user. */
/*          (This will be checked by the solver.) */

/* RWORK  = A real working array (double precision). */
/*          The length of RWORK must be at least 20 + NEQ. */
/*          The first 20 words of RWORK are reserved for conditional */
/*          and optional input and optional output. */

/*          The following word in RWORK is a conditional input: */
/*            RWORK(1) = TCRIT = critical value of t which the solver */
/*                       is not to overshoot.  Required if ITASK is */
/*                       4 or 5, and ignored otherwise.  (See ITASK.) */

/* LRW    = The length of the array RWORK, as declared by the user. */
/*          (This will be checked by the solver.) */

/* IWORK  = An integer work array.  The length of IWORK must be at least */
/*             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or */
/*             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25). */
/*          The first 30 words of IWORK are reserved for conditional and */
/*          optional input and optional output. */

/*          The following 2 words in IWORK are conditional input: */
/*            IWORK(1) = ML     These are the lower and upper */
/*            IWORK(2) = MU     half-bandwidths, respectively, of the */
/*                       banded Jacobian, excluding the main diagonal. */
/*                       The band is defined by the matrix locations */
/*                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU */
/*                       must satisfy  0 .le.  ML,MU  .le. NEQ-1. */
/*                       These are required if MITER is 4 or 5, and */
/*                       ignored otherwise.  ML and MU may in fact be */
/*                       the band parameters for a matrix to which */
/*                       df/dy is only approximately equal. */

/* LIW    = the length of the array IWORK, as declared by the user. */
/*          (This will be checked by the solver.) */

/* Note:  The work arrays must not be altered between calls to ZVODE */
/* for the same problem, except possibly for the conditional and */
/* optional input, and except for the last 2*NEQ words of ZWORK and */
/* the last NEQ words of RWORK.  The latter space is used for internal */
/* scratch space, and so is available for use by the user outside ZVODE */
/* between calls, if desired (but not for use by F or JAC). */

/* JAC    = The name of the user-supplied routine (MITER = 1 or 4) to */
/*          compute the Jacobian matrix, df/dy, as a function of */
/*          the scalar t and the vector y.  It is to have the form */
/*               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, */
/*                               RPAR, IPAR) */
/*               DOUBLE COMPLEX Y(NEQ), PD(NROWPD,NEQ) */
/*               DOUBLE PRECISION T */
/*          where NEQ, T, Y, ML, MU, and NROWPD are input and the array */
/*          PD is to be loaded with partial derivatives (elements of the */
/*          Jacobian matrix) in the output.  PD must be given a first */
/*          dimension of NROWPD.  T and Y have the same meaning as in */
/*          Subroutine F. */
/*               In the full matrix case (MITER = 1), ML and MU are */
/*          ignored, and the Jacobian is to be loaded into PD in */
/*          columnwise manner, with df(i)/dy(j) loaded into PD(i,j). */
/*               In the band matrix case (MITER = 4), the elements */
/*          within the band are to be loaded into PD in columnwise */
/*          manner, with diagonal lines of df/dy loaded into the rows */
/*          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j). */
/*          ML and MU are the half-bandwidth parameters. (See IWORK). */
/*          The locations in PD in the two triangular areas which */
/*          correspond to nonexistent matrix elements can be ignored */
/*          or loaded arbitrarily, as they are overwritten by ZVODE. */
/*               JAC need not provide df/dy exactly.  A crude */
/*          approximation (possibly with a smaller bandwidth) will do. */
/*               In either case, PD is preset to zero by the solver, */
/*          so that only the nonzero elements need be loaded by JAC. */
/*          Each call to JAC is preceded by a call to F with the same */
/*          arguments NEQ, T, and Y.  Thus to gain some efficiency, */
/*          intermediate quantities shared by both calculations may be */
/*          saved in a user COMMON block by F and not recomputed by JAC, */
/*          if desired.  Also, JAC may alter the Y array, if desired. */
/*          JAC must be declared external in the calling program. */
/*               Subroutine JAC may access user-defined real/complex and */
/*          integer work arrays, RPAR and IPAR, whose dimensions are set */
/*          by the user in the calling program. */

/* MF     = The method flag.  Used only for input.  The legal values of */
/*          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25, */
/*          -11, -12, -14, -15, -21, -22, -24, -25. */
/*          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER). */
/*          JSV = SIGN(MF) indicates the Jacobian-saving strategy: */
/*            JSV =  1 means a copy of the Jacobian is saved for reuse */
/*                     in the corrector iteration algorithm. */
/*            JSV = -1 means a copy of the Jacobian is not saved */
/*                     (valid only for MITER = 1, 2, 4, or 5). */
/*          METH indicates the basic linear multistep method: */
/*            METH = 1 means the implicit Adams method. */
/*            METH = 2 means the method based on backward */
/*                     differentiation formulas (BDF-s). */
/*          MITER indicates the corrector iteration method: */
/*            MITER = 0 means functional iteration (no Jacobian matrix */
/*                      is involved). */
/*            MITER = 1 means chord iteration with a user-supplied */
/*                      full (NEQ by NEQ) Jacobian. */
/*            MITER = 2 means chord iteration with an internally */
/*                      generated (difference quotient) full Jacobian */
/*                      (using NEQ extra calls to F per df/dy value). */
/*            MITER = 3 means chord iteration with an internally */
/*                      generated diagonal Jacobian approximation */
/*                      (using 1 extra call to F per df/dy evaluation). */
/*            MITER = 4 means chord iteration with a user-supplied */
/*                      banded Jacobian. */
/*            MITER = 5 means chord iteration with an internally */
/*                      generated banded Jacobian (using ML+MU+1 extra */
/*                      calls to F per df/dy evaluation). */
/*          If MITER = 1 or 4, the user must supply a subroutine JAC */
/*          (the name is arbitrary) as described above under JAC. */
/*          For other values of MITER, a dummy argument can be used. */

/* RPAR     User-specified array used to communicate real or complex */
/*          parameters to user-supplied subroutines.  If RPAR is an */
/*          array, it must be dimensioned in the user's calling program; */
/*          if it is unused or it is a scalar, then it need not be */
/*          dimensioned.  The type of RPAR may be REAL, DOUBLE, COMPLEX, */
/*          or DOUBLE COMPLEX, depending on the user program's needs. */
/*          RPAR is not type-declared within ZVODE, but simply passed */
/*          (by address) to the user's F and JAC routines. */

/* IPAR     User-specified array used to communicate integer parameter */
/*          to user-supplied subroutines.  If IPAR is an array, it must */
/*          be dimensioned in the user's calling program. */
/* ----------------------------------------------------------------------- */
/* Optional Input. */

/* The following is a list of the optional input provided for in the */
/* call sequence.  (See also Part ii.)  For each such input variable, */
/* this table lists its name as used in this documentation, its */
/* location in the call sequence, its meaning, and the default value. */
/* The use of any of this input requires IOPT = 1, and in that */
/* case all of this input is examined.  A value of zero for any */
/* of these optional input variables will cause the default value to be */
/* used.  Thus to use a subset of the optional input, simply preload */
/* locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and */
/* then set those of interest to nonzero values. */

/* NAME    LOCATION      MEANING AND DEFAULT VALUE */

/* H0      RWORK(5)  The step size to be attempted on the first step. */
/*                   The default value is determined by the solver. */

/* HMAX    RWORK(6)  The maximum absolute step size allowed. */
/*                   The default value is infinite. */

/* HMIN    RWORK(7)  The minimum absolute step size allowed. */
/*                   The default value is 0.  (This lower bound is not */
/*                   enforced on the final step before reaching TCRIT */
/*                   when ITASK = 4 or 5.) */

/* MAXORD  IWORK(5)  The maximum order to be allowed.  The default */
/*                   value is 12 if METH = 1, and 5 if METH = 2. */
/*                   If MAXORD exceeds the default value, it will */
/*                   be reduced to the default value. */
/*                   If MAXORD is changed during the problem, it may */
/*                   cause the current order to be reduced. */

/* MXSTEP  IWORK(6)  Maximum number of (internally defined) steps */
/*                   allowed during one call to the solver. */
/*                   The default value is 500. */

/* MXHNIL  IWORK(7)  Maximum number of messages printed (per problem) */
/*                   warning that T + H = T on a step (H = step size). */
/*                   This must be positive to result in a non-default */
/*                   value.  The default value is 10. */

/* ----------------------------------------------------------------------- */
/* Optional Output. */

/* As optional additional output from ZVODE, the variables listed */
/* below are quantities related to the performance of ZVODE */
/* which are available to the user.  These are communicated by way of */
/* the work arrays, but also have internal mnemonic names as shown. */
/* Except where stated otherwise, all of this output is defined */
/* on any successful return from ZVODE, and on any return with */
/* ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return */
/* (ISTATE = -3), they will be unchanged from their existing values */
/* (if any), except possibly for TOLSF, LENZW, LENRW, and LENIW. */
/* On any error return, output relevant to the error will be defined, */
/* as noted below. */

/* NAME    LOCATION      MEANING */

/* HU      RWORK(11) The step size in t last used (successfully). */

/* HCUR    RWORK(12) The step size to be attempted on the next step. */

/* TCUR    RWORK(13) The current value of the independent variable */
/*                   which the solver has actually reached, i.e. the */
/*                   current internal mesh point in t.  In the output, */
/*                   TCUR will always be at least as far from the */
/*                   initial value of t as the current argument T, */
/*                   but may be farther (if interpolation was done). */

/* TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0, */
/*                   computed when a request for too much accuracy was */
/*                   detected (ISTATE = -3 if detected at the start of */
/*                   the problem, ISTATE = -2 otherwise).  If ITOL is */
/*                   left unaltered but RTOL and ATOL are uniformly */
/*                   scaled up by a factor of TOLSF for the next call, */
/*                   then the solver is deemed likely to succeed. */
/*                   (The user may also ignore TOLSF and alter the */
/*                   tolerance parameters in any other way appropriate.) */

/* NST     IWORK(11) The number of steps taken for the problem so far. */

/* NFE     IWORK(12) The number of f evaluations for the problem so far. */

/* NJE     IWORK(13) The number of Jacobian evaluations so far. */

/* NQU     IWORK(14) The method order last used (successfully). */

/* NQCUR   IWORK(15) The order to be attempted on the next step. */

/* IMXER   IWORK(16) The index of the component of largest magnitude in */
/*                   the weighted local error vector ( e(i)/EWT(i) ), */
/*                   on an error return with ISTATE = -4 or -5. */

/* LENZW   IWORK(17) The length of ZWORK actually required. */
/*                   This is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* LENRW   IWORK(18) The length of RWORK actually required. */
/*                   This is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* LENIW   IWORK(19) The length of IWORK actually required. */
/*                   This is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* NLU     IWORK(20) The number of matrix LU decompositions so far. */

/* NNI     IWORK(21) The number of nonlinear (Newton) iterations so far. */

/* NCFN    IWORK(22) The number of convergence failures of the nonlinear */
/*                   solver so far. */

/* NETF    IWORK(23) The number of error test failures of the integrator */
/*                   so far. */

/* The following two arrays are segments of the ZWORK array which */
/* may also be of interest to the user as optional output. */
/* For each array, the table below gives its internal name, */
/* its base address in ZWORK, and its description. */

/* NAME    BASE ADDRESS      DESCRIPTION */

/* YH       1             The Nordsieck history array, of size NYH by */
/*                        (NQCUR + 1), where NYH is the initial value */
/*                        of NEQ.  For j = 0,1,...,NQCUR, column j+1 */
/*                        of YH contains HCUR**j/factorial(j) times */
/*                        the j-th derivative of the interpolating */
/*                        polynomial currently representing the */
/*                        solution, evaluated at t = TCUR. */

/* ACOR     LENZW-NEQ+1   Array of size NEQ used for the accumulated */
/*                        corrections on each step, scaled in the output */
/*                        to represent the estimated local error in Y */
/*                        on the last step.  This is the vector e in */
/*                        the description of the error control.  It is */
/*                        defined only on a successful return from ZVODE. */

/* ----------------------------------------------------------------------- */
/* Interrupting and Restarting */

/* If the integration of a given problem by ZVODE is to be */
/* interrupted and then later continued, such as when restarting */
/* an interrupted run or alternating between two or more ODE problems, */
/* the user should save, following the return from the last ZVODE call */
/* prior to the interruption, the contents of the call sequence */
/* variables and internal COMMON blocks, and later restore these */
/* values before the next ZVODE call for that problem.  To save */
/* and restore the COMMON blocks, use subroutine ZVSRCO, as */
/* described below in part ii. */

/* In addition, if non-default values for either LUN or MFLAG are */
/* desired, an extra call to XSETUN and/or XSETF should be made just */
/* before continuing the integration.  See Part ii below for details. */

/* ----------------------------------------------------------------------- */
/* Part ii.  Other Routines Callable. */

/* The following are optional calls which the user may make to */
/* gain additional capabilities in conjunction with ZVODE. */
/* (The routines XSETUN and XSETF are designed to conform to the */
/* SLATEC error handling package.) */

/*     FORM OF CALL                  FUNCTION */
/*  CALL XSETUN(LUN)           Set the logical unit number, LUN, for */
/*                             output of messages from ZVODE, if */
/*                             the default is not desired. */
/*                             The default value of LUN is 6. */

/*  CALL XSETF(MFLAG)          Set a flag to control the printing of */
/*                             messages by ZVODE. */
/*                             MFLAG = 0 means do not print. (Danger: */
/*                             This risks losing valuable information.) */
/*                             MFLAG = 1 means print (the default). */

/*                             Either of the above calls may be made at */
/*                             any time and will take effect immediately. */

/*  CALL ZVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of */
/*                             the internal COMMON blocks used by */
/*                             ZVODE. (See Part iii below.) */
/*                             RSAV must be a real array of length 51 */
/*                             or more, and ISAV must be an integer */
/*                             array of length 40 or more. */
/*                             JOB=1 means save COMMON into RSAV/ISAV. */
/*                             JOB=2 means restore COMMON from RSAV/ISAV. */
/*                                ZVSRCO is useful if one is */
/*                             interrupting a run and restarting */
/*                             later, or alternating between two or */
/*                             more problems solved with ZVODE. */

/*  CALL ZVINDY(,,,,,)         Provide derivatives of y, of various */
/*        (See below.)         orders, at a specified point T, if */
/*                             desired.  It may be called only after */
/*                             a successful return from ZVODE. */

/* The detailed instructions for using ZVINDY are as follows. */
/* The form of the call is: */

/*  CALL ZVINDY (T, K, ZWORK, NYH, DKY, IFLAG) */

/* The input parameters are: */

/* T         = Value of independent variable where answers are desired */
/*             (normally the same as the T last returned by ZVODE). */
/*             For valid results, T must lie between TCUR - HU and TCUR. */
/*             (See optional output for TCUR and HU.) */
/* K         = Integer order of the derivative desired.  K must satisfy */
/*             0 .le. K .le. NQCUR, where NQCUR is the current order */
/*             (see optional output).  The capability corresponding */
/*             to K = 0, i.e. computing y(T), is already provided */
/*             by ZVODE directly.  Since NQCUR .ge. 1, the first */
/*             derivative dy/dt is always available with ZVINDY. */
/* ZWORK     = The history array YH. */
/* NYH       = Column length of YH, equal to the initial value of NEQ. */

/* The output parameters are: */

/* DKY       = A double complex array of length NEQ containing the */
/*             computed value of the K-th derivative of y(t). */
/* IFLAG     = Integer flag, returned as 0 if K and T were legal, */
/*             -1 if K was illegal, and -2 if T was illegal. */
/*             On an error return, a message is also written. */
/* ----------------------------------------------------------------------- */
/* Part iii.  COMMON Blocks. */
/* If ZVODE is to be used in an overlay situation, the user */
/* must declare, in the primary overlay, the variables in: */
/*   (1) the call sequence to ZVODE, */
/*   (2) the two internal COMMON blocks */
/*         /ZVOD01/  of length  83  (50 double precision words */
/*                         followed by 33 integer words), */
/*         /ZVOD02/  of length  9  (1 double precision word */
/*                         followed by 8 integer words), */

/* If ZVODE is used on a system in which the contents of internal */
/* COMMON blocks are not preserved between calls, the user should */
/* declare the above two COMMON blocks in his calling program to insure */
/* that their contents are preserved. */

/* ----------------------------------------------------------------------- */
/* Part iv.  Optionally Replaceable Solver Routines. */

/* Below are descriptions of two routines in the ZVODE package which */
/* relate to the measurement of errors.  Either routine can be */
/* replaced by a user-supplied version, if desired.  However, since such */
/* a replacement may have a major impact on performance, it should be */
/* done only when absolutely necessary, and only with great caution. */
/* (Note: The means by which the package version of a routine is */
/* superseded by the user's version may be system-dependent.) */

/* (a) ZEWSET. */
/* The following subroutine is called just before each internal */
/* integration step, and sets the array of error weights, EWT, as */
/* described under ITOL/RTOL/ATOL above: */
/*     SUBROUTINE ZEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT) */
/* where NEQ, ITOL, RTOL, and ATOL are as in the ZVODE call sequence, */
/* YCUR contains the current (double complex) dependent variable vector, */
/* and EWT is the array of weights set by ZEWSET. */

/* If the user supplies this subroutine, it must return in EWT(i) */
/* (i = 1,...,NEQ) a positive quantity suitable for comparison with */
/* errors in Y(i).  The EWT array returned by ZEWSET is passed to the */
/* ZVNORM routine (See below.), and also used by ZVODE in the computation */
/* of the optional output IMXER, the diagonal Jacobian approximation, */
/* and the increments for difference quotient Jacobians. */

/* In the user-supplied version of ZEWSET, it may be desirable to use */
/* the current values of derivatives of y.  Derivatives up to order NQ */
/* are available from the history array YH, described above under */
/* Optional Output.  In ZEWSET, YH is identical to the YCUR array, */
/* extended to NQ + 1 columns with a column length of NYH and scale */
/* factors of h**j/factorial(j).  On the first call for the problem, */
/* given by NST = 0, NQ is 1 and H is temporarily set to 1.0. */
/* NYH is the initial value of NEQ.  The quantities NQ, H, and NST */
/* can be obtained by including in ZEWSET the statements: */
/*     DOUBLE PRECISION RVOD, H, HU */
/*     COMMON /ZVOD01/ RVOD(50), IVOD(33) */
/*     COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST */
/*     NQ = IVOD(28) */
/*     H = RVOD(21) */
/* Thus, for example, the current value of dy/dt can be obtained as */
/* YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is */
/* unnecessary when NST = 0). */

/* (b) ZVNORM. */
/* The following is a real function routine which computes the weighted */
/* root-mean-square norm of a vector v: */
/*     D = ZVNORM (N, V, W) */
/* where: */
/*   N = the length of the vector, */
/*   V = double complex array of length N containing the vector, */
/*   W = real array of length N containing weights, */
/*   D = sqrt( (1/N) * sum(abs(V(i))*W(i))**2 ). */
/* ZVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where */
/* EWT is as set by subroutine ZEWSET. */

/* If the user supplies this function, it should return a non-negative */
/* value of ZVNORM suitable for use in the error control in ZVODE. */
/* None of the arguments should be altered by ZVNORM. */
/* For example, a user-supplied ZVNORM routine might: */
/*   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or */
/*   -ignore some components of V in the norm, with the effect of */
/*    suppressing the error control on those components of Y. */
/* ----------------------------------------------------------------------- */
/* REVISION HISTORY (YYYYMMDD) */
/*  20060517  DATE WRITTEN, modified from DVODE of 20020430. */
/*  20061227  Added note on use for analytic f. */
/* ----------------------------------------------------------------------- */
/* Other Routines in the ZVODE Package. */

/* In addition to Subroutine ZVODE, the ZVODE package includes the */
/* following subroutines and function routines: */
/*  ZVHIN     computes an approximate step size for the initial step. */
/*  ZVINDY    computes an interpolated value of the y vector at t = TOUT. */
/*  ZVSTEP    is the core integrator, which does one step of the */
/*            integration and the associated error control. */
/*  ZVSET     sets all method coefficients and test constants. */
/*  ZVNLSD    solves the underlying nonlinear system -- the corrector. */
/*  ZVJAC     computes and preprocesses the Jacobian matrix J = df/dy */
/*            and the Newton iteration matrix P = I - (h/l1)*J. */
/*  ZVSOL     manages solution of linear system in chord iteration. */
/*  ZVJUST    adjusts the history array on a change of order. */
/*  ZEWSET    sets the error weight vector EWT before each step. */
/*  ZVNORM    computes the weighted r.m.s. norm of a vector. */
/*  ZABSSQ    computes the squared absolute value of a double complex z. */
/*  ZVSRCO    is a user-callable routine to save and restore */
/*            the contents of the internal COMMON blocks. */
/*  ZACOPY    is a routine to copy one two-dimensional array to another. */
/*  ZGETRF and ZGETRS   are routines from LAPACK for solving full */
/*            systems of linear algebraic equations. */
/*  ZGBTRF and ZGBTRS   are routines from LAPACK for solving banded */
/*            linear systems. */
/*  DZSCAL    scales a double complex array by a double prec. scalar. */
/*  DZAXPY    adds a D.P. scalar times one complex vector to another. */
/*  ZCOPY     is a basic linear algebra module from the BLAS. */
/*  DUMACH    sets the unit roundoff of the machine. */
/*  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all */
/*            error messages and warnings.  XERRWD is machine-dependent. */
/* Note: ZVNORM, ZABSSQ, DUMACH, IXSAV, and IUMACH are function routines. */
/* All the others are subroutines. */
/* The intrinsic functions called with double precision complex arguments */
/* are: ABS, DREAL, and DIMAG.  All of these are expected to return */
/* double precision real values. */

/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for labeled COMMON block ZVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */


/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to ZVODE. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* The following internal COMMON blocks contain variables which are */
/* communicated between subroutines in the ZVODE package, or which are */
/* to be saved between calls to ZVODE. */
/* In each block, real variables precede integers. */
/* The block /ZVOD01/ appears in subroutines ZVODE, ZVINDY, ZVSTEP, */
/* ZVSET, ZVNLSD, ZVJAC, ZVSOL, ZVJUST and ZVSRCO. */
/* The block /ZVOD02/ appears in subroutines ZVODE, ZVINDY, ZVSTEP, */
/* ZVNLSD, ZVJAC, and ZVSRCO. */

/* The variables stored in the internal COMMON blocks are as follows: */

/* ACNRM  = Weighted r.m.s. norm of accumulated correction vectors. */
/* CCMXJ  = Threshold on DRC for updating the Jacobian. (See DRC.) */
/* CONP   = The saved value of TQ(5). */
/* CRATE  = Estimated corrector convergence rate constant. */
/* DRC    = Relative change in H*RL1 since last ZVJAC call. */
/* EL     = Real array of integration coefficients.  See ZVSET. */
/* ETA    = Saved tentative ratio of new to old H. */
/* ETAMAX = Saved maximum value of ETA to be allowed. */
/* H      = The step size. */
/* HMIN   = The minimum absolute value of the step size H to be used. */
/* HMXI   = Inverse of the maximum absolute value of H to be used. */
/*          HMXI = 0.0 is allowed and corresponds to an infinite HMAX. */
/* HNEW   = The step size to be attempted on the next step. */
/* HRL1   = Saved value of H*RL1. */
/* HSCAL  = Stepsize in scaling of YH array. */
/* PRL1   = The saved value of RL1. */
/* RC     = Ratio of current H*RL1 to value on last ZVJAC call. */
/* RL1    = The reciprocal of the coefficient EL(1). */
/* SRUR   = Sqrt(UROUND), used in difference quotient algorithms. */
/* TAU    = Real vector of past NQ step sizes, length 13. */
/* TQ     = A real vector of length 5 in which ZVSET stores constants */
/*          used for the convergence test, the error test, and the */
/*          selection of H at a new order. */
/* TN     = The independent variable, updated on each step taken. */
/* UROUND = The machine unit roundoff.  The smallest positive real number */
/*          such that  1.0 + UROUND .ne. 1.0 */
/* ICF    = Integer flag for convergence failure in ZVNLSD: */
/*            0 means no failures. */
/*            1 means convergence failure with out of date Jacobian */
/*                   (recoverable error). */
/*            2 means convergence failure with current Jacobian or */
/*                   singular matrix (unrecoverable error). */
/* INIT   = Saved integer flag indicating whether initialization of the */
/*          problem has been done (INIT = 1) or not. */
/* IPUP   = Saved flag to signal updating of Newton matrix. */
/* JCUR   = Output flag from ZVJAC showing Jacobian status: */
/*            JCUR = 0 means J is not current. */
/*            JCUR = 1 means J is current. */
/* JSTART = Integer flag used as input to ZVSTEP: */
/*            0  means perform the first step. */
/*            1  means take a new step continuing from the last. */
/*            -1 means take the next step with a new value of MAXORD, */
/*                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters. */
/*          On return, ZVSTEP sets JSTART = 1. */
/* JSV    = Integer flag for Jacobian saving, = sign(MF). */
/* KFLAG  = A completion code from ZVSTEP with the following meanings: */
/*               0      the step was successful. */
/*              -1      the requested error could not be achieved. */
/*              -2      corrector convergence could not be achieved. */
/*              -3, -4  fatal error in VNLS (can not occur here). */
/* KUTH   = Input flag to ZVSTEP showing whether H was reduced by the */
/*          driver.  KUTH = 1 if H was reduced, = 0 otherwise. */
/* L      = Integer variable, NQ + 1, current order plus one. */
/* LMAX   = MAXORD + 1 (used for dimensioning). */
/* LOCJS  = A pointer to the saved Jacobian, whose storage starts at */
/*          WM(LOCJS), if JSV = 1. */
/* LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers */
/*          to segments of ZWORK, RWORK, and IWORK. */
/* MAXORD = The maximum order of integration method to be allowed. */
/* METH/MITER = The method flags.  See MF. */
/* MSBJ   = The maximum number of steps between J evaluations, = 50. */
/* MXHNIL = Saved value of optional input MXHNIL. */
/* MXSTEP = Saved value of optional input MXSTEP. */
/* N      = The number of first-order ODEs, = NEQ. */
/* NEWH   = Saved integer to flag change of H. */
/* NEWQ   = The method order to be used on the next step. */
/* NHNIL  = Saved counter for occurrences of T + H = T. */
/* NQ     = Integer variable, the current integration method order. */
/* NQNYH  = Saved value of NQ*NYH. */
/* NQWAIT = A counter controlling the frequency of order changes. */
/*          An order change is about to be considered if NQWAIT = 1. */
/* NSLJ   = The number of steps taken as of the last Jacobian update. */
/* NSLP   = Saved value of NST as of last Newton matrix update. */
/* NYH    = Saved value of the initial value of NEQ. */
/* HU     = The step size in t last used. */
/* NCFN   = Number of nonlinear convergence failures so far. */
/* NETF   = The number of error test failures of the integrator so far. */
/* NFE    = The number of f evaluations for the problem so far. */
/* NJE    = The number of Jacobian evaluations so far. */
/* NLU    = The number of matrix LU decompositions so far. */
/* NNI    = Number of nonlinear iterations so far. */
/* NQU    = The method order last used. */
/* NST    = The number of steps taken for the problem so far. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --zwork;
    --rwork;
    --iwork;
    --rpar;
    --ipar;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Block A. */
/* This code block is executed on every call. */
/* It tests ISTATE and ITASK for legality and branches appropriately. */
/* If ISTATE .gt. 1 but the flag INIT shows that initialization has */
/* not yet been done, an error return occurs. */
/* If ISTATE = 1 and TOUT = T, return immediately. */
/* ----------------------------------------------------------------------- */
    if (*istate < 1 || *istate > 3) {
	goto L601;
    }
    if (*itask < 1 || *itask > 5) {
	goto L602;
    }
    if (*istate == 1) {
	goto L10;
    }
    if (zvod01_1.init != 1) {
	goto L603;
    }
    if (*istate == 2) {
	goto L200;
    }
    goto L20;
L10:
    zvod01_1.init = 0;
    if (*tout == *t) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/* Block B. */
/* The next code block is executed for the initial call (ISTATE = 1), */
/* or for a continuation call with parameter changes (ISTATE = 3). */
/* It contains checking of all input and various initializations. */

/* First check legality of the non-optional input NEQ, ITOL, IOPT, */
/* MF, ML, and MU. */
/* ----------------------------------------------------------------------- */
L20:
    if (*neq <= 0) {
	goto L604;
    }
    if (*istate == 1) {
	goto L25;
    }
    if (*neq > zvod01_1.n) {
	goto L605;
    }
L25:
    zvod01_1.n = *neq;
    if (*itol < 1 || *itol > 4) {
	goto L606;
    }
    if (*iopt < 0 || *iopt > 1) {
	goto L607;
    }
    zvod01_1.jsv = i_sign(&c__1, mf);
    mfa = abs(*mf);
    zvod01_1.meth = mfa / 10;
    zvod01_1.miter = mfa - zvod01_1.meth * 10;
    if (zvod01_1.meth < 1 || zvod01_1.meth > 2) {
	goto L608;
    }
    if (zvod01_1.miter < 0 || zvod01_1.miter > 5) {
	goto L608;
    }
    if (zvod01_1.miter <= 3) {
	goto L30;
    }
    ml = iwork[1];
    mu = iwork[2];
    if (ml < 0 || ml >= zvod01_1.n) {
	goto L609;
    }
    if (mu < 0 || mu >= zvod01_1.n) {
	goto L610;
    }
L30:
/* Next process and check the optional input. --------------------------- */
    if (*iopt == 1) {
	goto L40;
    }
    zvod01_1.maxord = mord[zvod01_1.meth - 1];
    zvod01_1.mxstep = mxstp0;
    zvod01_1.mxhnil = mxhnl0;
    if (*istate == 1) {
	h0 = zero;
    }
    zvod01_1.hmxi = zero;
    zvod01_1.hmin = zero;
    goto L60;
L40:
    zvod01_1.maxord = iwork[5];
    if (zvod01_1.maxord < 0) {
	goto L611;
    }
    if (zvod01_1.maxord == 0) {
	zvod01_1.maxord = 100;
    }
/* Computing MIN */
    i__1 = zvod01_1.maxord, i__2 = mord[zvod01_1.meth - 1];
    zvod01_1.maxord = min(i__1,i__2);
    zvod01_1.mxstep = iwork[6];
    if (zvod01_1.mxstep < 0) {
	goto L612;
    }
    if (zvod01_1.mxstep == 0) {
	zvod01_1.mxstep = mxstp0;
    }
    zvod01_1.mxhnil = iwork[7];
    if (zvod01_1.mxhnil < 0) {
	goto L613;
    }
    if (zvod01_1.mxhnil == 0) {
	zvod01_1.mxhnil = mxhnl0;
    }
    if (*istate != 1) {
	goto L50;
    }
    h0 = rwork[5];
    if ((*tout - *t) * h0 < zero) {
	goto L614;
    }
L50:
    hmax = rwork[6];
    if (hmax < zero) {
	goto L615;
    }
    zvod01_1.hmxi = zero;
    if (hmax > zero) {
	zvod01_1.hmxi = one / hmax;
    }
    zvod01_1.hmin = rwork[7];
    if (zvod01_1.hmin < zero) {
	goto L616;
    }
/* ----------------------------------------------------------------------- */
/* Set work array pointers and check lengths LZW, LRW, and LIW. */
/* Pointers to segments of ZWORK, RWORK, and IWORK are named by prefixing */
/* L to the name of the segment.  E.g., segment YH starts at ZWORK(LYH). */
/* Segments of ZWORK (in order) are denoted  YH, WM, SAVF, ACOR. */
/* Besides optional inputs/outputs, RWORK has only the segment EWT. */
/* Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0). */
/* ----------------------------------------------------------------------- */
L60:
    zvod01_1.lyh = 1;
    if (*istate == 1) {
	zvod01_1.nyh = zvod01_1.n;
    }
    zvod01_1.lwm = zvod01_1.lyh + (zvod01_1.maxord + 1) * zvod01_1.nyh;
    jco = max(0,zvod01_1.jsv);
    if (zvod01_1.miter == 0) {
	lenwm = 0;
    }
    if (zvod01_1.miter == 1 || zvod01_1.miter == 2) {
	lenwm = (jco + 1) * zvod01_1.n * zvod01_1.n;
	zvod01_1.locjs = zvod01_1.n * zvod01_1.n + 1;
    }
    if (zvod01_1.miter == 3) {
	lenwm = zvod01_1.n;
    }
    if (zvod01_1.miter == 4 || zvod01_1.miter == 5) {
	mband = ml + mu + 1;
	lenp = (mband + ml) * zvod01_1.n;
	lenj = mband * zvod01_1.n;
	lenwm = lenp + jco * lenj;
	zvod01_1.locjs = lenp + 1;
    }
    zvod01_1.lsavf = zvod01_1.lwm + lenwm;
    zvod01_1.lacor = zvod01_1.lsavf + zvod01_1.n;
    lenzw = zvod01_1.lacor + zvod01_1.n - 1;
    iwork[17] = lenzw;
    zvod01_1.lewt = 21;
    lenrw = zvod01_1.n + 20;
    iwork[18] = lenrw;
    zvod01_1.liwm = 1;
    leniw = zvod01_1.n + 30;
    if (zvod01_1.miter == 0 || zvod01_1.miter == 3) {
	leniw = 30;
    }
    iwork[19] = leniw;
    if (lenzw > *lzw) {
	goto L628;
    }
    if (lenrw > *lrw) {
	goto L617;
    }
    if (leniw > *liw) {
	goto L618;
    }
/* Check RTOL and ATOL for legality. ------------------------------------ */
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol >= 3) {
	    rtoli = rtol[i__];
	}
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	if (rtoli < zero) {
	    goto L619;
	}
	if (atoli < zero) {
	    goto L620;
	}
/* L70: */
    }
    if (*istate == 1) {
	goto L100;
    }
/* If ISTATE = 3, set flag to signal parameter changes to ZVSTEP. ------- */
    zvod01_1.jstart = -1;
    if (zvod01_1.nq <= zvod01_1.maxord) {
	goto L200;
    }
/* MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. --------- */
    zcopy_(&zvod01_1.n, &zwork[zvod01_1.lwm], &c__1, &zwork[zvod01_1.lsavf], &
	    c__1);
    goto L200;
/* ----------------------------------------------------------------------- */
/* Block C. */
/* The next block is for the initial call only (ISTATE = 1). */
/* It contains all remaining initializations, the initial call to F, */
/* and the calculation of the initial step size. */
/* The error weights in EWT are inverted after being loaded. */
/* ----------------------------------------------------------------------- */
L100:
    zvod01_1.uround = dumach_();
    zvod01_1.tn = *t;
    if (*itask != 4 && *itask != 5) {
	goto L110;
    }
    tcrit = rwork[1];
    if ((tcrit - *tout) * (*tout - *t) < zero) {
	goto L625;
    }
    if (h0 != zero && (*t + h0 - tcrit) * h0 > zero) {
	h0 = tcrit - *t;
    }
L110:
    zvod01_1.jstart = 0;
    if (zvod01_1.miter > 0) {
	zvod01_1.srur = sqrt(zvod01_1.uround);
    }
    zvod01_1.ccmxj = pt2;
    zvod01_1.msbj = 50;
    zvod01_1.nhnil = 0;
    zvod02_1.nst = 0;
    zvod02_1.nje = 0;
    zvod02_1.nni = 0;
    zvod02_1.ncfn = 0;
    zvod02_1.netf = 0;
    zvod02_1.nlu = 0;
    zvod01_1.nslj = 0;
    nslast = 0;
    zvod02_1.hu = zero;
    zvod02_1.nqu = 0;
/* Initial call to F.  (LF0 points to YH(*,2).) ------------------------- */
    lf0 = zvod01_1.lyh + zvod01_1.nyh;
    (*f)(&zvod01_1.n, t, &y[1], &zwork[lf0], &rpar[1], &ipar[1]);
    zvod02_1.nfe = 1;
/* Load the initial value vector in YH. --------------------------------- */
    zcopy_(&zvod01_1.n, &y[1], &c__1, &zwork[zvod01_1.lyh], &c__1);
/* Load and invert the EWT array.  (H is temporarily set to 1.0.) ------- */
    zvod01_1.nq = 1;
    zvod01_1.h__ = one;
    zewset_(&zvod01_1.n, itol, &rtol[1], &atol[1], &zwork[zvod01_1.lyh], &
	    rwork[zvod01_1.lewt]);
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + zvod01_1.lewt - 1] <= zero) {
	    goto L621;
	}
/* L120: */
	rwork[i__ + zvod01_1.lewt - 1] = one / rwork[i__ + zvod01_1.lewt - 1];
    }
    if (h0 != zero) {
	goto L180;
    }
/* Call ZVHIN to set initial step size H0 to be attempted. -------------- */
    zvhin_(&zvod01_1.n, t, &zwork[zvod01_1.lyh], &zwork[lf0], (S_fp)f, &rpar[
	    1], &ipar[1], tout, &zvod01_1.uround, &rwork[zvod01_1.lewt], itol,
	     &atol[1], &y[1], &zwork[zvod01_1.lacor], &h0, &niter, &ier);
    zvod02_1.nfe += niter;
    if (ier != 0) {
	goto L622;
    }
/* Adjust H0 if necessary to meet HMAX bound. --------------------------- */
L180:
    rh = abs(h0) * zvod01_1.hmxi;
    if (rh > one) {
	h0 /= rh;
    }
/* Load H with H0 and scale YH(*,2) by H0. ------------------------------ */
    zvod01_1.h__ = h0;
    dzscal_(&zvod01_1.n, &h0, &zwork[lf0], &c__1);
    goto L270;
/* ----------------------------------------------------------------------- */
/* Block D. */
/* The next code block is for continuation calls only (ISTATE = 2 or 3) */
/* and is to check stop conditions before taking a step. */
/* ----------------------------------------------------------------------- */
L200:
    nslast = zvod02_1.nst;
    zvod01_1.kuth = 0;
    switch (*itask) {
	case 1:  goto L210;
	case 2:  goto L250;
	case 3:  goto L220;
	case 4:  goto L230;
	case 5:  goto L240;
    }
L210:
    if ((zvod01_1.tn - *tout) * zvod01_1.h__ < zero) {
	goto L250;
    }
    zvindy_(tout, &c__0, &zwork[zvod01_1.lyh], &zvod01_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L220:
    tp = zvod01_1.tn - zvod02_1.hu * (one + hun * zvod01_1.uround);
    if ((tp - *tout) * zvod01_1.h__ > zero) {
	goto L623;
    }
    if ((zvod01_1.tn - *tout) * zvod01_1.h__ < zero) {
	goto L250;
    }
    goto L400;
L230:
    tcrit = rwork[1];
    if ((zvod01_1.tn - tcrit) * zvod01_1.h__ > zero) {
	goto L624;
    }
    if ((tcrit - *tout) * zvod01_1.h__ < zero) {
	goto L625;
    }
    if ((zvod01_1.tn - *tout) * zvod01_1.h__ < zero) {
	goto L245;
    }
    zvindy_(tout, &c__0, &zwork[zvod01_1.lyh], &zvod01_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L240:
    tcrit = rwork[1];
    if ((zvod01_1.tn - tcrit) * zvod01_1.h__ > zero) {
	goto L624;
    }
L245:
    hmx = abs(zvod01_1.tn) + abs(zvod01_1.h__);
    ihit = (d__1 = zvod01_1.tn - tcrit, abs(d__1)) <= hun * zvod01_1.uround * 
	    hmx;
    if (ihit) {
	goto L400;
    }
    tnext = zvod01_1.tn + zvod01_1.hnew * (one + four * zvod01_1.uround);
    if ((tnext - tcrit) * zvod01_1.h__ <= zero) {
	goto L250;
    }
    zvod01_1.h__ = (tcrit - zvod01_1.tn) * (one - four * zvod01_1.uround);
    zvod01_1.kuth = 1;
/* ----------------------------------------------------------------------- */
/* Block E. */
/* The next block is normally executed for all calls and contains */
/* the call to the one-step core integrator ZVSTEP. */

/* This is a looping point for the integration steps. */

/* First check for too many steps being taken, update EWT (if not at */
/* start of problem), check for too much accuracy being requested, and */
/* check for H below the roundoff level in T. */
/* ----------------------------------------------------------------------- */
L250:
    if (zvod02_1.nst - nslast >= zvod01_1.mxstep) {
	goto L500;
    }
    zewset_(&zvod01_1.n, itol, &rtol[1], &atol[1], &zwork[zvod01_1.lyh], &
	    rwork[zvod01_1.lewt]);
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + zvod01_1.lewt - 1] <= zero) {
	    goto L510;
	}
/* L260: */
	rwork[i__ + zvod01_1.lewt - 1] = one / rwork[i__ + zvod01_1.lewt - 1];
    }
L270:
    tolsf = zvod01_1.uround * zvnorm_(&zvod01_1.n, &zwork[zvod01_1.lyh], &
	    rwork[zvod01_1.lewt]);
    if (tolsf <= one) {
	goto L280;
    }
    tolsf *= two;
    if (zvod02_1.nst == 0) {
	goto L626;
    }
    goto L520;
L280:
    if (zvod01_1.tn + zvod01_1.h__ != zvod01_1.tn) {
	goto L290;
    }
    ++zvod01_1.nhnil;
    if (zvod01_1.nhnil > zvod01_1.mxhnil) {
	goto L290;
    }
    s_copy(msg, "ZVODE--  Warning: internal T (=R1) and H (=R2) are", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__101, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      such that in the machine, T + H = T on the next step  "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__101, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      (H = step size). solver will continue anyway", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__101, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    zvod01_1.tn, &zvod01_1.h__, (ftnlen)80);
    if (zvod01_1.nhnil < zvod01_1.mxhnil) {
	goto L290;
    }
    s_copy(msg, "ZVODE--  Above warning has been issued I1 times.  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__102, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      it will not be issued again for this problem", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__102, &c__1, &c__1, &zvod01_1.mxhnil, &c__0, &
	    c__0, &zero, &zero, (ftnlen)80);
L290:
/* ----------------------------------------------------------------------- */
/* CALL ZVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR, */
/*              WM, IWM, F, JAC, F, ZVNLSD, RPAR, IPAR) */
/* ----------------------------------------------------------------------- */
    zvstep_(&y[1], &zwork[zvod01_1.lyh], &zvod01_1.nyh, &zwork[zvod01_1.lyh], 
	    &rwork[zvod01_1.lewt], &zwork[zvod01_1.lsavf], &y[1], &zwork[
	    zvod01_1.lacor], &zwork[zvod01_1.lwm], &iwork[zvod01_1.liwm], (
	    S_fp)f, (U_fp)jac, (S_fp)f, (U_fp)zvnlsd_, &rpar[1], &ipar[1]);
    kgo = 1 - zvod01_1.kflag;
/* Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3. */
/*  KFLAG .eq. 0,   -1,  -2 */
    switch (kgo) {
	case 1:  goto L300;
	case 2:  goto L530;
	case 3:  goto L540;
    }
/* ----------------------------------------------------------------------- */
/* Block F. */
/* The following block handles the case of a successful return from the */
/* core integrator (KFLAG = 0).  Test for stop conditions. */
/* ----------------------------------------------------------------------- */
L300:
    zvod01_1.init = 1;
    zvod01_1.kuth = 0;
    switch (*itask) {
	case 1:  goto L310;
	case 2:  goto L400;
	case 3:  goto L330;
	case 4:  goto L340;
	case 5:  goto L350;
    }
/* ITASK = 1.  If TOUT has been reached, interpolate. ------------------- */
L310:
    if ((zvod01_1.tn - *tout) * zvod01_1.h__ < zero) {
	goto L250;
    }
    zvindy_(tout, &c__0, &zwork[zvod01_1.lyh], &zvod01_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
/* ITASK = 3.  Jump to exit if TOUT was reached. ------------------------ */
L330:
    if ((zvod01_1.tn - *tout) * zvod01_1.h__ >= zero) {
	goto L400;
    }
    goto L250;
/* ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary. */
L340:
    if ((zvod01_1.tn - *tout) * zvod01_1.h__ < zero) {
	goto L345;
    }
    zvindy_(tout, &c__0, &zwork[zvod01_1.lyh], &zvod01_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
L345:
    hmx = abs(zvod01_1.tn) + abs(zvod01_1.h__);
    ihit = (d__1 = zvod01_1.tn - tcrit, abs(d__1)) <= hun * zvod01_1.uround * 
	    hmx;
    if (ihit) {
	goto L400;
    }
    tnext = zvod01_1.tn + zvod01_1.hnew * (one + four * zvod01_1.uround);
    if ((tnext - tcrit) * zvod01_1.h__ <= zero) {
	goto L250;
    }
    zvod01_1.h__ = (tcrit - zvod01_1.tn) * (one - four * zvod01_1.uround);
    zvod01_1.kuth = 1;
    goto L250;
/* ITASK = 5.  See if TCRIT was reached and jump to exit. --------------- */
L350:
    hmx = abs(zvod01_1.tn) + abs(zvod01_1.h__);
    ihit = (d__1 = zvod01_1.tn - tcrit, abs(d__1)) <= hun * zvod01_1.uround * 
	    hmx;
/* ----------------------------------------------------------------------- */
/* Block G. */
/* The following block handles all successful returns from ZVODE. */
/* If ITASK .ne. 1, Y is loaded from YH and T is set accordingly. */
/* ISTATE is set to 2, and the optional output is loaded into the work */
/* arrays before returning. */
/* ----------------------------------------------------------------------- */
L400:
    zcopy_(&zvod01_1.n, &zwork[zvod01_1.lyh], &c__1, &y[1], &c__1);
    *t = zvod01_1.tn;
    if (*itask != 4 && *itask != 5) {
	goto L420;
    }
    if (ihit) {
	*t = tcrit;
    }
L420:
    *istate = 2;
    rwork[11] = zvod02_1.hu;
    rwork[12] = zvod01_1.hnew;
    rwork[13] = zvod01_1.tn;
    iwork[11] = zvod02_1.nst;
    iwork[12] = zvod02_1.nfe;
    iwork[13] = zvod02_1.nje;
    iwork[14] = zvod02_1.nqu;
    iwork[15] = zvod01_1.newq;
    iwork[20] = zvod02_1.nlu;
    iwork[21] = zvod02_1.nni;
    iwork[22] = zvod02_1.ncfn;
    iwork[23] = zvod02_1.netf;
    return 0;
/* ----------------------------------------------------------------------- */
/* Block H. */
/* The following block handles all unsuccessful returns other than */
/* those for illegal input.  First the error message routine is called. */
/* if there was an error test or convergence test failure, IMXER is set. */
/* Then Y is loaded from YH, and T is set to TN. */
/* The optional output is loaded into the work arrays before returning. */
/* ----------------------------------------------------------------------- */
/* The maximum number of steps was taken before reaching TOUT. ---------- */
/* Error message removed, see gh-7888 */
L500:
    *istate = -1;
    goto L580;
/* EWT(i) .le. 0.0 for some i (not at start of problem). ---------------- */
L510:
    ewti = rwork[zvod01_1.lewt + i__ - 1];
    s_copy(msg, "ZVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__202, &c__1, &c__1, &i__, &c__0, &c__2, &
	    zvod01_1.tn, &ewti, (ftnlen)80);
    *istate = -6;
    goto L580;
/* Too much accuracy requested for machine precision. ------------------- */
L520:
    s_copy(msg, "ZVODE--  At T (=R1), too much accuracy requested  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__203, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      for precision of machine:   see TOLSF (=R2) ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__203, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    zvod01_1.tn, &tolsf, (ftnlen)80);
    rwork[14] = tolsf;
    *istate = -2;
    goto L580;
/* KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. ----- */
L530:
    s_copy(msg, "ZVODE--  At T(=R1) and step size H(=R2), the error", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__204, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      test failed repeatedly or with abs(H) = HMIN", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__204, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    zvod01_1.tn, &zvod01_1.h__, (ftnlen)80);
    *istate = -4;
    goto L560;
/* KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ---- */
L540:
    s_copy(msg, "ZVODE--  At T (=R1) and step size H (=R2), the    ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__205, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      corrector convergence failed repeatedly     ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__205, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      or with abs(H) = HMIN   ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__205, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    zvod01_1.tn, &zvod01_1.h__, (ftnlen)80);
    *istate = -5;
/* Compute IMXER if relevant. ------------------------------------------- */
L560:
    big = zero;
    imxer = 1;
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	size = z_abs(&zwork[i__ + zvod01_1.lacor - 1]) * rwork[i__ + 
		zvod01_1.lewt - 1];
	if (big >= size) {
	    goto L570;
	}
	big = size;
	imxer = i__;
L570:
	;
    }
    iwork[16] = imxer;
/* Set Y vector, T, and optional output. -------------------------------- */
L580:
    zcopy_(&zvod01_1.n, &zwork[zvod01_1.lyh], &c__1, &y[1], &c__1);
    *t = zvod01_1.tn;
    rwork[11] = zvod02_1.hu;
    rwork[12] = zvod01_1.h__;
    rwork[13] = zvod01_1.tn;
    iwork[11] = zvod02_1.nst;
    iwork[12] = zvod02_1.nfe;
    iwork[13] = zvod02_1.nje;
    iwork[14] = zvod02_1.nqu;
    iwork[15] = zvod01_1.nq;
    iwork[20] = zvod02_1.nlu;
    iwork[21] = zvod02_1.nni;
    iwork[22] = zvod02_1.ncfn;
    iwork[23] = zvod02_1.netf;
    return 0;
/* ----------------------------------------------------------------------- */
/* Block I. */
/* The following block handles all error returns due to illegal input */
/* (ISTATE = -3), as detected before calling the core integrator. */
/* First the error message routine is called.   If the illegal input */
/* is a negative ISTATE, the run is aborted (apparent infinite loop). */
/* ----------------------------------------------------------------------- */
L601:
    s_copy(msg, "ZVODE--  ISTATE (=I1) illegal ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__1, &c__1, &c__1, istate, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    if (*istate < 0) {
	goto L800;
    }
    goto L700;
L602:
    s_copy(msg, "ZVODE--  ITASK (=I1) illegal  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__2, &c__1, &c__1, itask, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L603:
    s_copy(msg, "ZVODE--  ISTATE (=I1) .gt. 1 but ZVODE not initialized      "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__3, &c__1, &c__1, istate, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L604:
    s_copy(msg, "ZVODE--  NEQ (=I1) .lt. 1     ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__4, &c__1, &c__1, neq, &c__0, &c__0, &zero, &zero,
	     (ftnlen)80);
    goto L700;
L605:
    s_copy(msg, "ZVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__5, &c__1, &c__2, &zvod01_1.n, neq, &c__0, &zero, 
	    &zero, (ftnlen)80);
    goto L700;
L606:
    s_copy(msg, "ZVODE--  ITOL (=I1) illegal   ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__6, &c__1, &c__1, itol, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L607:
    s_copy(msg, "ZVODE--  IOPT (=I1) illegal   ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__7, &c__1, &c__1, iopt, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L608:
    s_copy(msg, "ZVODE--  MF (=I1) illegal     ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__8, &c__1, &c__1, mf, &c__0, &c__0, &zero, &zero, 
	    (ftnlen)80);
    goto L700;
L609:
    s_copy(msg, "ZVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__9, &c__1, &c__2, &ml, neq, &c__0, &zero, &zero, (
	    ftnlen)80);
    goto L700;
L610:
    s_copy(msg, "ZVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__10, &c__1, &c__2, &mu, neq, &c__0, &zero, &zero, 
	    (ftnlen)80);
    goto L700;
L611:
    s_copy(msg, "ZVODE--  MAXORD (=I1) .lt. 0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__11, &c__1, &c__1, &zvod01_1.maxord, &c__0, &c__0,
	     &zero, &zero, (ftnlen)80);
    goto L700;
L612:
    s_copy(msg, "ZVODE--  MXSTEP (=I1) .lt. 0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__12, &c__1, &c__1, &zvod01_1.mxstep, &c__0, &c__0,
	     &zero, &zero, (ftnlen)80);
    goto L700;
L613:
    s_copy(msg, "ZVODE--  MXHNIL (=I1) .lt. 0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__13, &c__1, &c__1, &zvod01_1.mxhnil, &c__0, &c__0,
	     &zero, &zero, (ftnlen)80);
    goto L700;
L614:
    s_copy(msg, "ZVODE--  TOUT (=R1) behind T (=R2)      ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__14, &c__1, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    ftnlen)80);
    s_copy(msg, "      integration direction is given by H0 (=R1)  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__14, &c__1, &c__0, &c__0, &c__0, &c__1, &h0, &
	    zero, (ftnlen)80);
    goto L700;
L615:
    s_copy(msg, "ZVODE--  HMAX (=R1) .lt. 0.0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__15, &c__1, &c__0, &c__0, &c__0, &c__1, &hmax, &
	    zero, (ftnlen)80);
    goto L700;
L616:
    s_copy(msg, "ZVODE--  HMIN (=R1) .lt. 0.0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__16, &c__1, &c__0, &c__0, &c__0, &c__1, &
	    zvod01_1.hmin, &zero, (ftnlen)80);
    goto L700;
L617:
    s_copy(msg, "ZVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__17, &c__1, &c__2, &lenrw, lrw, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L618:
    s_copy(msg, "ZVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__18, &c__1, &c__2, &leniw, liw, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L619:
    s_copy(msg, "ZVODE--  RTOL(I1) is R1 .lt. 0.0        ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__19, &c__1, &c__1, &i__, &c__0, &c__1, &rtoli, &
	    zero, (ftnlen)80);
    goto L700;
L620:
    s_copy(msg, "ZVODE--  ATOL(I1) is R1 .lt. 0.0        ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__20, &c__1, &c__1, &i__, &c__0, &c__1, &atoli, &
	    zero, (ftnlen)80);
    goto L700;
L621:
    ewti = rwork[zvod01_1.lewt + i__ - 1];
    s_copy(msg, "ZVODE--  EWT(I1) is R1 .le. 0.0         ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__21, &c__1, &c__1, &i__, &c__0, &c__1, &ewti, &
	    zero, (ftnlen)80);
    goto L700;
L622:
    s_copy(msg, "ZVODE--  TOUT (=R1) too close to T(=R2) to start integration"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__22, &c__1, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    ftnlen)80);
    goto L700;
L623:
    s_copy(msg, "ZVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__23, &c__1, &c__1, itask, &c__0, &c__2, tout, &tp,
	     (ftnlen)80);
    goto L700;
L624:
    s_copy(msg, "ZVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__24, &c__1, &c__0, &c__0, &c__0, &c__2, &tcrit, &
	    zvod01_1.tn, (ftnlen)80);
    goto L700;
L625:
    s_copy(msg, "ZVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__25, &c__1, &c__0, &c__0, &c__0, &c__2, &tcrit, 
	    tout, (ftnlen)80);
    goto L700;
L626:
    s_copy(msg, "ZVODE--  At start of problem, too much accuracy   ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__26, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      requested for precision of machine:   see TOLSF (=R1) "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__26, &c__1, &c__0, &c__0, &c__0, &c__1, &tolsf, &
	    zero, (ftnlen)80);
    rwork[14] = tolsf;
    goto L700;
L627:
    s_copy(msg, "ZVODE--  Trouble from ZVINDY.  ITASK = I1, TOUT = R1.       "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__27, &c__1, &c__1, itask, &c__0, &c__1, tout, &
	    zero, (ftnlen)80);
    goto L700;
L628:
    s_copy(msg, "ZVODE--  ZWORK length needed, LENZW (=I1), exceeds LZW (=I2)"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__17, &c__1, &c__2, &lenzw, lzw, &c__0, &zero, &
	    zero, (ftnlen)80);

L700:
    *istate = -3;
    return 0;

L800:
    s_copy(msg, "ZVODE--  Run aborted:  apparent infinite loop     ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__303, &c__2, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    return 0;
/* ----------------------- End of Subroutine ZVODE ----------------------- */
} /* zvode_ */

/* DECK ZVHIN */
/* Subroutine */ int zvhin_(integer *n, doublereal *t0, doublecomplex *y0, 
	doublecomplex *ydot, S_fp f, real *rpar, integer *ipar, doublereal *
	tout, doublereal *uround, doublereal *ewt, integer *itol, doublereal *
	atol, doublecomplex *y, doublecomplex *temp, doublereal *h0, integer *
	niter, integer *ier)
{
    /* Initialized data */

    static doublereal half = .5;
    static doublereal hun = 100.;
    static doublereal pt1 = .1;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal), d_sign(doublereal *, 
	    doublereal *);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal t1, hg, afi, hlb, hub, hrat, hnew;
    static integer iter;
    static doublereal delyi, atoli, tdist, yddnrm, tround;
    extern doublereal zvnorm_(integer *, doublecomplex *, doublereal *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND, */
/*                        EWT, ITOL, ATOL, Y, TEMP */
/* Call sequence output -- H0, NITER, IER */
/* COMMON block variables accessed -- None */

/* Subroutines called by ZVHIN:  F */
/* Function routines called by ZVHIN: ZVNORM */
/* ----------------------------------------------------------------------- */
/* This routine computes the step size, H0, to be attempted on the */
/* first step, when the user has not supplied a value for this. */

/* First we check that TOUT - T0 differs significantly from zero.  Then */
/* an iteration is done to approximate the initial second derivative */
/* and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1. */
/* A bias factor of 1/2 is applied to the resulting h. */
/* The sign of H0 is inferred from the initial values of TOUT and T0. */

/* Communication with ZVHIN is done with the following variables: */

/* N      = Size of ODE system, input. */
/* T0     = Initial value of independent variable, input. */
/* Y0     = Vector of initial conditions, input. */
/* YDOT   = Vector of initial first derivatives, input. */
/* F      = Name of subroutine for right-hand side f(t,y), input. */
/* RPAR, IPAR = User's real/complex and integer work arrays. */
/* TOUT   = First output value of independent variable */
/* UROUND = Machine unit roundoff */
/* EWT, ITOL, ATOL = Error weights and tolerance parameters */
/*                   as described in the driver routine, input. */
/* Y, TEMP = Work arrays of length N. */
/* H0     = Step size to be attempted, output. */
/* NITER  = Number of iterations (and of f evaluations) to compute H0, */
/*          output. */
/* IER    = The error flag, returned with the value */
/*          IER = 0  if no trouble occurred, or */
/*          IER = -1 if TOUT and T0 are considered too close to proceed. */
/* ----------------------------------------------------------------------- */

/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --temp;
    --y;
    --atol;
    --ewt;
    --ipar;
    --rpar;
    --ydot;
    --y0;

    /* Function Body */

    *niter = 0;
    tdist = (d__1 = *tout - *t0, abs(d__1));
/* Computing MAX */
    d__1 = abs(*t0), d__2 = abs(*tout);
    tround = *uround * max(d__1,d__2);
    if (tdist < two * tround) {
	goto L100;
    }

/* Set a lower bound on h based on the roundoff level in T0 and TOUT. --- */
    hlb = hun * tround;
/* Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. - */
    hub = pt1 * tdist;
    atoli = atol[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	delyi = pt1 * z_abs(&y0[i__]) + atoli;
	afi = z_abs(&ydot[i__]);
	if (afi * hub > delyi) {
	    hub = delyi / afi;
	}
/* L10: */
    }

/* Set initial guess for h as geometric mean of upper and lower bounds. - */
    iter = 0;
    hg = sqrt(hlb * hub);
/* If the bounds have crossed, exit with the mean value. ---------------- */
    if (hub < hlb) {
	*h0 = hg;
	goto L90;
    }

/* Looping point for iteration. ----------------------------------------- */
L50:
/* Estimate the second derivative as a difference quotient in f. -------- */
    d__1 = *tout - *t0;
    h__ = d_sign(&hg, &d__1);
    t1 = *t0 + h__;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	i__2 = i__;
	i__3 = i__;
	i__4 = i__;
	z__2.r = h__ * ydot[i__4].r, z__2.i = h__ * ydot[i__4].i;
	z__1.r = y0[i__3].r + z__2.r, z__1.i = y0[i__3].i + z__2.i;
	y[i__2].r = z__1.r, y[i__2].i = z__1.i;
    }
    (*f)(n, &t1, &y[1], &temp[1], &rpar[1], &ipar[1]);
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L70: */
	i__3 = i__;
	i__4 = i__;
	i__1 = i__;
	z__2.r = temp[i__4].r - ydot[i__1].r, z__2.i = temp[i__4].i - ydot[
		i__1].i;
	z__1.r = z__2.r / h__, z__1.i = z__2.i / h__;
	temp[i__3].r = z__1.r, temp[i__3].i = z__1.i;
    }
    yddnrm = zvnorm_(n, &temp[1], &ewt[1]);
/* Get the corresponding new value of h. -------------------------------- */
    if (yddnrm * hub * hub > two) {
	hnew = sqrt(two / yddnrm);
    } else {
	hnew = sqrt(hg * hub);
    }
    ++iter;
/* ----------------------------------------------------------------------- */
/* Test the stopping conditions. */
/* Stop if the new and previous h values differ by a factor of .lt. 2. */
/* Stop if four iterations have been done.  Also, stop with previous h */
/* if HNEW/HG .gt. 2 after first iteration, as this probably means that */
/* the second derivative value is bad because of cancellation error. */
/* ----------------------------------------------------------------------- */
    if (iter >= 4) {
	goto L80;
    }
    hrat = hnew / hg;
    if (hrat > half && hrat < two) {
	goto L80;
    }
    if (iter >= 2 && hnew > two * hg) {
	hnew = hg;
	goto L80;
    }
    hg = hnew;
    goto L50;

/* Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ---- */
L80:
    *h0 = hnew * half;
    if (*h0 < hlb) {
	*h0 = hlb;
    }
    if (*h0 > hub) {
	*h0 = hub;
    }
L90:
    d__1 = *tout - *t0;
    *h0 = d_sign(h0, &d__1);
    *niter = iter;
    *ier = 0;
    return 0;
/* Error return for TOUT - T0 too small. -------------------------------- */
L100:
    *ier = -1;
    return 0;
/* ----------------------- End of Subroutine ZVHIN ----------------------- */
} /* zvhin_ */

/* DECK ZVINDY */
/* Subroutine */ int zvindy_(doublereal *t, integer *k, doublecomplex *yh, 
	integer *ldyh, doublecomplex *dky, integer *iflag)
{
    /* Initialized data */

    static doublereal hun = 100.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal r__, s;
    static integer ic, jb, jj;
    static doublereal tp;
    static integer jb2, jj1, jp1;
    static doublereal tn1;
    static char msg[80];
    static doublereal tfuzz;
    extern /* Subroutine */ int dzscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerrwd_(char *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, ftnlen);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- T, K, YH, LDYH */
/* Call sequence output -- DKY, IFLAG */
/* COMMON block variables accessed: */
/*     /ZVOD01/ --  H, TN, UROUND, L, N, NQ */
/*     /ZVOD02/ --  HU */

/* Subroutines called by ZVINDY: DZSCAL, XERRWD */
/* Function routines called by ZVINDY: None */
/* ----------------------------------------------------------------------- */
/* ZVINDY computes interpolated values of the K-th derivative of the */
/* dependent variable vector y, and stores it in DKY.  This routine */
/* is called within the package with K = 0 and T = TOUT, but may */
/* also be called by the user for any K up to the current order. */
/* (See detailed instructions in the usage documentation.) */
/* ----------------------------------------------------------------------- */
/* The computed values in DKY are gotten by interpolation using the */
/* Nordsieck history array YH.  This array corresponds uniquely to a */
/* vector-valued polynomial of degree NQCUR or less, and DKY is set */
/* to the K-th derivative of this polynomial at T. */
/* The formula for DKY is: */
/*              q */
/*  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1) */
/*             j=K */
/* where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR. */
/* The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are */
/* communicated by COMMON.  The above sum is done in reverse order. */
/* IFLAG is returned negative if either K or T is out of bounds. */

/* Discussion above and comments in driver explain all variables. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for labeled COMMON block ZVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --dky;

    /* Function Body */

    *iflag = 0;
    if (*k < 0 || *k > zvod01_1.nq) {
	goto L80;
    }
    d__1 = abs(zvod01_1.tn) + abs(zvod02_1.hu);
    tfuzz = hun * zvod01_1.uround * d_sign(&d__1, &zvod02_1.hu);
    tp = zvod01_1.tn - zvod02_1.hu - tfuzz;
    tn1 = zvod01_1.tn + tfuzz;
    if ((*t - tp) * (*t - tn1) > zero) {
	goto L90;
    }

    s = (*t - zvod01_1.tn) / zvod01_1.h__;
    ic = 1;
    if (*k == 0) {
	goto L15;
    }
    jj1 = zvod01_1.l - *k;
    i__1 = zvod01_1.nq;
    for (jj = jj1; jj <= i__1; ++jj) {
/* L10: */
	ic *= jj;
    }
L15:
    c__ = (real) ic;
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	i__2 = i__;
	i__3 = i__ + zvod01_1.l * yh_dim1;
	z__1.r = c__ * yh[i__3].r, z__1.i = c__ * yh[i__3].i;
	dky[i__2].r = z__1.r, dky[i__2].i = z__1.i;
    }
    if (*k == zvod01_1.nq) {
	goto L55;
    }
    jb2 = zvod01_1.nq - *k;
    i__2 = jb2;
    for (jb = 1; jb <= i__2; ++jb) {
	j = zvod01_1.nq - jb;
	jp1 = j + 1;
	ic = 1;
	if (*k == 0) {
	    goto L35;
	}
	jj1 = jp1 - *k;
	i__3 = j;
	for (jj = jj1; jj <= i__3; ++jj) {
/* L30: */
	    ic *= jj;
	}
L35:
	c__ = (real) ic;
	i__3 = zvod01_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L40: */
	    i__1 = i__;
	    i__4 = i__ + jp1 * yh_dim1;
	    z__2.r = c__ * yh[i__4].r, z__2.i = c__ * yh[i__4].i;
	    i__5 = i__;
	    z__3.r = s * dky[i__5].r, z__3.i = s * dky[i__5].i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    dky[i__1].r = z__1.r, dky[i__1].i = z__1.i;
	}
/* L50: */
    }
    if (*k == 0) {
	return 0;
    }
L55:
    i__2 = -(*k);
    r__ = pow_di(&zvod01_1.h__, &i__2);
    dzscal_(&zvod01_1.n, &r__, &dky[1], &c__1);
    return 0;

L80:
    s_copy(msg, "ZVINDY-- K (=I1) illegal      ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__51, &c__1, &c__1, k, &c__0, &c__0, &zero, &zero, 
	    (ftnlen)80);
    *iflag = -1;
    return 0;
L90:
    s_copy(msg, "ZVINDY-- T (=R1) illegal      ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__52, &c__1, &c__0, &c__0, &c__0, &c__1, t, &zero, 
	    (ftnlen)80);
    s_copy(msg, "      T not in interval TCUR - HU (= R1) to TCUR (=R2)      "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__52, &c__1, &c__0, &c__0, &c__0, &c__2, &tp, &
	    zvod01_1.tn, (ftnlen)80);
    *iflag = -2;
    return 0;
/* ----------------------- End of Subroutine ZVINDY ---------------------- */
} /* zvindy_ */

/* DECK ZVSTEP */
/* Subroutine */ int zvstep_(doublecomplex *y, doublecomplex *yh, integer *
	ldyh, doublecomplex *yh1, doublereal *ewt, doublecomplex *savf, 
	doublecomplex *vsav, doublecomplex *acor, doublecomplex *wm, integer *
	iwm, S_fp f, U_fp jac, S_fp psol, S_fp vnls, real *rpar, integer *
	ipar)
{
    /* Initialized data */

    static integer kfc = -3;
    static integer kfh = -7;
    static integer mxncf = 10;
    static doublereal addon = 1e-6;
    static doublereal bias1 = 6.;
    static doublereal bias2 = 6.;
    static doublereal bias3 = 10.;
    static doublereal etacf = .25;
    static doublereal etamin = .1;
    static doublereal etamxf = .2;
    static doublereal etamx1 = 1e4;
    static doublereal etamx2 = 10.;
    static doublereal etamx3 = 10.;
    static doublereal onepsm = 1.00001;
    static doublereal thresh = 1.5;
    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;

    /* Local variables */
    static integer i__, j;
    static doublereal r__;
    static integer i1, i2, jb;
    static doublereal ddn;
    static integer ncf;
    static doublereal dsm, dup, etaq, told;
    static integer iback, nflag;
    static doublereal flotl;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zvset_(void);
    static doublereal etaqm1, etaqp1;
    extern /* Subroutine */ int dzscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal cnquot;
    extern /* Subroutine */ int dzaxpy_(integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal zvnorm_(integer *, doublecomplex *, doublereal *);
    extern /* Subroutine */ int zvjust_(doublecomplex *, integer *, integer *)
	    ;

/* ----------------------------------------------------------------------- */
/* Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV, */
/*                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR */
/* Call sequence output -- YH, ACOR, WM, IWM */
/* COMMON block variables accessed: */
/*     /ZVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13), */
/*               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH, */
/*               L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT */
/*     /ZVOD02/  HU, NCFN, NETF, NFE, NQU, NST */

/* Subroutines called by ZVSTEP: F, DZAXPY, ZCOPY, DZSCAL, */
/*                               ZVJUST, VNLS, ZVSET */
/* Function routines called by ZVSTEP: ZVNORM */
/* ----------------------------------------------------------------------- */
/* ZVSTEP performs one step of the integration of an initial value */
/* problem for a system of ordinary differential equations. */
/* ZVSTEP calls subroutine VNLS for the solution of the nonlinear system */
/* arising in the time step.  Thus it is independent of the problem */
/* Jacobian structure and the type of nonlinear system solution method. */
/* ZVSTEP returns a completion flag KFLAG (in COMMON). */
/* A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10 */
/* consecutive failures occurred.  On a return with KFLAG negative, */
/* the values of TN and the YH array are as of the beginning of the last */
/* step, and H is the last step size attempted. */

/* Communication with ZVSTEP is done with the following variables: */

/* Y      = An array of length N used for the dependent variable vector. */
/* YH     = An LDYH by LMAX array containing the dependent variables */
/*          and their approximate scaled derivatives, where */
/*          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate */
/*          j-th derivative of y(i), scaled by H**j/factorial(j) */
/*          (j = 0,1,...,NQ).  On entry for the first step, the first */
/*          two columns of YH must be set from the initial values. */
/* LDYH   = A constant integer .ge. N, the first dimension of YH. */
/*          N is the number of ODEs in the system. */
/* YH1    = A one-dimensional array occupying the same space as YH. */
/* EWT    = An array of length N containing multiplicative weights */
/*          for local error measurements.  Local errors in y(i) are */
/*          compared to 1.0/EWT(i) in various error tests. */
/* SAVF   = An array of working storage, of length N. */
/*          also used for input of YH(*,MAXORD+2) when JSTART = -1 */
/*          and MAXORD .lt. the current order NQ. */
/* VSAV   = A work array of length N passed to subroutine VNLS. */
/* ACOR   = A work array of length N, used for the accumulated */
/*          corrections.  On a successful return, ACOR(i) contains */
/*          the estimated one-step local error in y(i). */
/* WM,IWM = Complex and integer work arrays associated with matrix */
/*          operations in VNLS. */
/* F      = Dummy name for the user-supplied subroutine for f. */
/* JAC    = Dummy name for the user-supplied Jacobian subroutine. */
/* PSOL   = Dummy name for the subroutine passed to VNLS, for */
/*          possible use there. */
/* VNLS   = Dummy name for the nonlinear system solving subroutine, */
/*          whose real name is dependent on the method used. */
/* RPAR, IPAR = User's real/complex and integer work arrays. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for labeled COMMON block ZVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --yh1;
    --ewt;
    --savf;
    --vsav;
    --acor;
    --wm;
    --iwm;
    --rpar;
    --ipar;

    /* Function Body */

    zvod01_1.kflag = 0;
    told = zvod01_1.tn;
    ncf = 0;
    zvod01_1.jcur = 0;
    nflag = 0;
    if (zvod01_1.jstart > 0) {
	goto L20;
    }
    if (zvod01_1.jstart == -1) {
	goto L100;
    }
/* ----------------------------------------------------------------------- */
/* On the first call, the order is set to 1, and other variables are */
/* initialized.  ETAMAX is the maximum ratio by which H can be increased */
/* in a single step.  It is normally 10, but is larger during the */
/* first step to compensate for the small initial H.  If a failure */
/* occurs (in corrector convergence or error test), ETAMAX is set to 1 */
/* for the next increase. */
/* ----------------------------------------------------------------------- */
    zvod01_1.lmax = zvod01_1.maxord + 1;
    zvod01_1.nq = 1;
    zvod01_1.l = 2;
    zvod01_1.nqnyh = zvod01_1.nq * *ldyh;
    zvod01_1.tau[0] = zvod01_1.h__;
    zvod01_1.prl1 = one;
    zvod01_1.rc = zero;
    zvod01_1.etamax = etamx1;
    zvod01_1.nqwait = 2;
    zvod01_1.hscal = zvod01_1.h__;
    goto L200;
/* ----------------------------------------------------------------------- */
/* Take preliminary actions on a normal continuation step (JSTART.GT.0). */
/* If the driver changed H, then ETA must be reset and NEWH set to 1. */
/* If a change of order was dictated on the previous step, then */
/* it is done here and appropriate adjustments in the history are made. */
/* On an order decrease, the history array is adjusted by ZVJUST. */
/* On an order increase, the history array is augmented by a column. */
/* On a change of step size H, the history array YH is rescaled. */
/* ----------------------------------------------------------------------- */
L20:
    if (zvod01_1.kuth == 1) {
/* Computing MIN */
	d__1 = zvod01_1.eta, d__2 = zvod01_1.h__ / zvod01_1.hscal;
	zvod01_1.eta = min(d__1,d__2);
	zvod01_1.newh = 1;
    }
L50:
    if (zvod01_1.newh == 0) {
	goto L200;
    }
    if (zvod01_1.newq == zvod01_1.nq) {
	goto L150;
    }
    if (zvod01_1.newq < zvod01_1.nq) {
	zvjust_(&yh[yh_offset], ldyh, &c_n1);
	zvod01_1.nq = zvod01_1.newq;
	zvod01_1.l = zvod01_1.nq + 1;
	zvod01_1.nqwait = zvod01_1.l;
	goto L150;
    }
    if (zvod01_1.newq > zvod01_1.nq) {
	zvjust_(&yh[yh_offset], ldyh, &c__1);
	zvod01_1.nq = zvod01_1.newq;
	zvod01_1.l = zvod01_1.nq + 1;
	zvod01_1.nqwait = zvod01_1.l;
	goto L150;
    }
/* ----------------------------------------------------------------------- */
/* The following block handles preliminaries needed when JSTART = -1. */
/* If N was reduced, zero out part of YH to avoid undefined references. */
/* If MAXORD was reduced to a value less than the tentative order NEWQ, */
/* then NQ is set to MAXORD, and a new H ratio ETA is chosen. */
/* Otherwise, we take the same preliminary actions as for JSTART .gt. 0. */
/* In any case, NQWAIT is reset to L = NQ + 1 to prevent further */
/* changes in order for that many steps. */
/* The new H ratio ETA is limited by the input H if KUTH = 1, */
/* by HMIN if KUTH = 0, and by HMXI in any case. */
/* Finally, the history array YH is rescaled. */
/* ----------------------------------------------------------------------- */
L100:
    zvod01_1.lmax = zvod01_1.maxord + 1;
    if (zvod01_1.n == *ldyh) {
	goto L120;
    }
    i1 = (zvod01_1.newq + 1) * *ldyh + 1;
    i2 = (zvod01_1.maxord + 1) * *ldyh;
    if (i1 > i2) {
	goto L120;
    }
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
/* L110: */
	i__2 = i__;
	yh1[i__2].r = zero, yh1[i__2].i = 0.;
    }
L120:
    if (zvod01_1.newq <= zvod01_1.maxord) {
	goto L140;
    }
    flotl = (real) zvod01_1.lmax;
    if (zvod01_1.maxord < zvod01_1.nq - 1) {
	ddn = zvnorm_(&zvod01_1.n, &savf[1], &ewt[1]) / zvod01_1.tq[0];
	d__1 = bias1 * ddn;
	d__2 = one / flotl;
	zvod01_1.eta = one / (pow_dd(&d__1, &d__2) + addon);
    }
    if (zvod01_1.maxord == zvod01_1.nq && zvod01_1.newq == zvod01_1.nq + 1) {
	zvod01_1.eta = etaq;
    }
    if (zvod01_1.maxord == zvod01_1.nq - 1 && zvod01_1.newq == zvod01_1.nq + 
	    1) {
	zvod01_1.eta = etaqm1;
	zvjust_(&yh[yh_offset], ldyh, &c_n1);
    }
    if (zvod01_1.maxord == zvod01_1.nq - 1 && zvod01_1.newq == zvod01_1.nq) {
	ddn = zvnorm_(&zvod01_1.n, &savf[1], &ewt[1]) / zvod01_1.tq[0];
	d__1 = bias1 * ddn;
	d__2 = one / flotl;
	zvod01_1.eta = one / (pow_dd(&d__1, &d__2) + addon);
	zvjust_(&yh[yh_offset], ldyh, &c_n1);
    }
    zvod01_1.eta = min(zvod01_1.eta,one);
    zvod01_1.nq = zvod01_1.maxord;
    zvod01_1.l = zvod01_1.lmax;
L140:
    if (zvod01_1.kuth == 1) {
/* Computing MIN */
	d__2 = zvod01_1.eta, d__3 = (d__1 = zvod01_1.h__ / zvod01_1.hscal, 
		abs(d__1));
	zvod01_1.eta = min(d__2,d__3);
    }
    if (zvod01_1.kuth == 0) {
/* Computing MAX */
	d__1 = zvod01_1.eta, d__2 = zvod01_1.hmin / abs(zvod01_1.hscal);
	zvod01_1.eta = max(d__1,d__2);
    }
/* Computing MAX */
    d__1 = one, d__2 = abs(zvod01_1.hscal) * zvod01_1.hmxi * zvod01_1.eta;
    zvod01_1.eta /= max(d__1,d__2);
    zvod01_1.newh = 1;
    zvod01_1.nqwait = zvod01_1.l;
    if (zvod01_1.newq <= zvod01_1.maxord) {
	goto L50;
    }
/* Rescale the history array for a change in H by a factor of ETA. ------ */
L150:
    r__ = one;
    i__2 = zvod01_1.l;
    for (j = 2; j <= i__2; ++j) {
	r__ *= zvod01_1.eta;
	dzscal_(&zvod01_1.n, &r__, &yh[j * yh_dim1 + 1], &c__1);
/* L180: */
    }
    zvod01_1.h__ = zvod01_1.hscal * zvod01_1.eta;
    zvod01_1.hscal = zvod01_1.h__;
    zvod01_1.rc *= zvod01_1.eta;
    zvod01_1.nqnyh = zvod01_1.nq * *ldyh;
/* ----------------------------------------------------------------------- */
/* This section computes the predicted values by effectively */
/* multiplying the YH array by the Pascal triangle matrix. */
/* ZVSET is called to calculate all integration coefficients. */
/* RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1. */
/* ----------------------------------------------------------------------- */
L200:
    zvod01_1.tn += zvod01_1.h__;
    i1 = zvod01_1.nqnyh + 1;
    i__2 = zvod01_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *ldyh;
	i__1 = zvod01_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L210: */
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__ + *ldyh;
	    z__1.r = yh1[i__4].r + yh1[i__5].r, z__1.i = yh1[i__4].i + yh1[
		    i__5].i;
	    yh1[i__3].r = z__1.r, yh1[i__3].i = z__1.i;
	}
/* L220: */
    }
    zvset_();
    zvod01_1.rl1 = one / zvod01_1.el[1];
    zvod01_1.rc *= zvod01_1.rl1 / zvod01_1.prl1;
    zvod01_1.prl1 = zvod01_1.rl1;

/* Call the nonlinear system solver. ------------------------------------ */

    (*vnls)(&y[1], &yh[yh_offset], ldyh, &vsav[1], &savf[1], &ewt[1], &acor[1]
	    , &iwm[1], &wm[1], (S_fp)f, (U_fp)jac, (S_fp)psol, &nflag, &rpar[
	    1], &ipar[1]);

    if (nflag == 0) {
	goto L450;
    }
/* ----------------------------------------------------------------------- */
/* The VNLS routine failed to achieve convergence (NFLAG .NE. 0). */
/* The YH array is retracted to its values before prediction. */
/* The step size H is reduced and the step is retried, if possible. */
/* Otherwise, an error exit is taken. */
/* ----------------------------------------------------------------------- */
    ++ncf;
    ++zvod02_1.ncfn;
    zvod01_1.etamax = one;
    zvod01_1.tn = told;
    i1 = zvod01_1.nqnyh + 1;
    i__2 = zvod01_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *ldyh;
	i__3 = zvod01_1.nqnyh;
	for (i__ = i1; i__ <= i__3; ++i__) {
/* L420: */
	    i__4 = i__;
	    i__5 = i__;
	    i__1 = i__ + *ldyh;
	    z__1.r = yh1[i__5].r - yh1[i__1].r, z__1.i = yh1[i__5].i - yh1[
		    i__1].i;
	    yh1[i__4].r = z__1.r, yh1[i__4].i = z__1.i;
	}
/* L430: */
    }
    if (nflag < -1) {
	goto L680;
    }
    if (abs(zvod01_1.h__) <= zvod01_1.hmin * onepsm) {
	goto L670;
    }
    if (ncf == mxncf) {
	goto L670;
    }
    zvod01_1.eta = etacf;
/* Computing MAX */
    d__1 = zvod01_1.eta, d__2 = zvod01_1.hmin / abs(zvod01_1.h__);
    zvod01_1.eta = max(d__1,d__2);
    nflag = -1;
    goto L150;
/* ----------------------------------------------------------------------- */
/* The corrector has converged (NFLAG = 0).  The local error test is */
/* made and control passes to statement 500 if it fails. */
/* ----------------------------------------------------------------------- */
L450:
    dsm = zvod01_1.acnrm / zvod01_1.tq[1];
    if (dsm > one) {
	goto L500;
    }
/* ----------------------------------------------------------------------- */
/* After a successful step, update the YH and TAU arrays and decrement */
/* NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved */
/* for use in a possible order increase on the next step. */
/* If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2. */
/* ----------------------------------------------------------------------- */
    zvod01_1.kflag = 0;
    ++zvod02_1.nst;
    zvod02_1.hu = zvod01_1.h__;
    zvod02_1.nqu = zvod01_1.nq;
    i__2 = zvod01_1.nq;
    for (iback = 1; iback <= i__2; ++iback) {
	i__ = zvod01_1.l - iback;
/* L470: */
	zvod01_1.tau[i__] = zvod01_1.tau[i__ - 1];
    }
    zvod01_1.tau[0] = zvod01_1.h__;
    i__2 = zvod01_1.l;
    for (j = 1; j <= i__2; ++j) {
	dzaxpy_(&zvod01_1.n, &zvod01_1.el[j - 1], &acor[1], &c__1, &yh[j * 
		yh_dim1 + 1], &c__1);
/* L480: */
    }
    --zvod01_1.nqwait;
    if (zvod01_1.l == zvod01_1.lmax || zvod01_1.nqwait != 1) {
	goto L490;
    }
    zcopy_(&zvod01_1.n, &acor[1], &c__1, &yh[zvod01_1.lmax * yh_dim1 + 1], &
	    c__1);
    zvod01_1.conp = zvod01_1.tq[4];
L490:
    if (zvod01_1.etamax != one) {
	goto L560;
    }
    if (zvod01_1.nqwait < 2) {
	zvod01_1.nqwait = 2;
    }
    zvod01_1.newq = zvod01_1.nq;
    zvod01_1.newh = 0;
    zvod01_1.eta = one;
    zvod01_1.hnew = zvod01_1.h__;
    goto L690;
/* ----------------------------------------------------------------------- */
/* The error test failed.  KFLAG keeps track of multiple failures. */
/* Restore TN and the YH array to their previous values, and prepare */
/* to try the step again.  Compute the optimum step size for the */
/* same order.  After repeated failures, H is forced to decrease */
/* more rapidly. */
/* ----------------------------------------------------------------------- */
L500:
    --zvod01_1.kflag;
    ++zvod02_1.netf;
    nflag = -2;
    zvod01_1.tn = told;
    i1 = zvod01_1.nqnyh + 1;
    i__2 = zvod01_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *ldyh;
	i__4 = zvod01_1.nqnyh;
	for (i__ = i1; i__ <= i__4; ++i__) {
/* L510: */
	    i__5 = i__;
	    i__1 = i__;
	    i__3 = i__ + *ldyh;
	    z__1.r = yh1[i__1].r - yh1[i__3].r, z__1.i = yh1[i__1].i - yh1[
		    i__3].i;
	    yh1[i__5].r = z__1.r, yh1[i__5].i = z__1.i;
	}
/* L520: */
    }
    if (abs(zvod01_1.h__) <= zvod01_1.hmin * onepsm) {
	goto L660;
    }
    zvod01_1.etamax = one;
    if (zvod01_1.kflag <= kfc) {
	goto L530;
    }
/* Compute ratio of new H to current H at the current order. ------------ */
    flotl = (real) zvod01_1.l;
    d__1 = bias2 * dsm;
    d__2 = one / flotl;
    zvod01_1.eta = one / (pow_dd(&d__1, &d__2) + addon);
/* Computing MAX */
    d__1 = zvod01_1.eta, d__2 = zvod01_1.hmin / abs(zvod01_1.h__), d__1 = max(
	    d__1,d__2);
    zvod01_1.eta = max(d__1,etamin);
    if (zvod01_1.kflag <= -2 && zvod01_1.eta > etamxf) {
	zvod01_1.eta = etamxf;
    }
    goto L150;
/* ----------------------------------------------------------------------- */
/* Control reaches this section if 3 or more consecutive failures */
/* have occurred.  It is assumed that the elements of the YH array */
/* have accumulated errors of the wrong order.  The order is reduced */
/* by one, if possible.  Then H is reduced by a factor of 0.1 and */
/* the step is retried.  After a total of 7 consecutive failures, */
/* an exit is taken with KFLAG = -1. */
/* ----------------------------------------------------------------------- */
L530:
    if (zvod01_1.kflag == kfh) {
	goto L660;
    }
    if (zvod01_1.nq == 1) {
	goto L540;
    }
/* Computing MAX */
    d__1 = etamin, d__2 = zvod01_1.hmin / abs(zvod01_1.h__);
    zvod01_1.eta = max(d__1,d__2);
    zvjust_(&yh[yh_offset], ldyh, &c_n1);
    zvod01_1.l = zvod01_1.nq;
    --zvod01_1.nq;
    zvod01_1.nqwait = zvod01_1.l;
    goto L150;
L540:
/* Computing MAX */
    d__1 = etamin, d__2 = zvod01_1.hmin / abs(zvod01_1.h__);
    zvod01_1.eta = max(d__1,d__2);
    zvod01_1.h__ *= zvod01_1.eta;
    zvod01_1.hscal = zvod01_1.h__;
    zvod01_1.tau[0] = zvod01_1.h__;
    (*f)(&zvod01_1.n, &zvod01_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++zvod02_1.nfe;
    i__2 = zvod01_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L550: */
	i__5 = i__ + (yh_dim1 << 1);
	i__1 = i__;
	z__1.r = zvod01_1.h__ * savf[i__1].r, z__1.i = zvod01_1.h__ * savf[
		i__1].i;
	yh[i__5].r = z__1.r, yh[i__5].i = z__1.i;
    }
    zvod01_1.nqwait = 10;
    goto L200;
/* ----------------------------------------------------------------------- */
/* If NQWAIT = 0, an increase or decrease in order by one is considered. */
/* Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could */
/* be multiplied at order q, q-1, or q+1, respectively. */
/* The largest of these is determined, and the new order and */
/* step size set accordingly. */
/* A change of H or NQ is made only if H increases by at least a */
/* factor of THRESH.  If an order change is considered and rejected, */
/* then NQWAIT is set to 2 (reconsider it after 2 steps). */
/* ----------------------------------------------------------------------- */
/* Compute ratio of new H to current H at the current order. ------------ */
L560:
    flotl = (real) zvod01_1.l;
    d__1 = bias2 * dsm;
    d__2 = one / flotl;
    etaq = one / (pow_dd(&d__1, &d__2) + addon);
    if (zvod01_1.nqwait != 0) {
	goto L600;
    }
    zvod01_1.nqwait = 2;
    etaqm1 = zero;
    if (zvod01_1.nq == 1) {
	goto L570;
    }
/* Compute ratio of new H to current H at the current order less one. --- */
    ddn = zvnorm_(&zvod01_1.n, &yh[zvod01_1.l * yh_dim1 + 1], &ewt[1]) / 
	    zvod01_1.tq[0];
    d__1 = bias1 * ddn;
    d__2 = one / (flotl - one);
    etaqm1 = one / (pow_dd(&d__1, &d__2) + addon);
L570:
    etaqp1 = zero;
    if (zvod01_1.l == zvod01_1.lmax) {
	goto L580;
    }
/* Compute ratio of new H to current H at current order plus one. ------- */
    d__1 = zvod01_1.h__ / zvod01_1.tau[1];
    cnquot = zvod01_1.tq[4] / zvod01_1.conp * pow_di(&d__1, &zvod01_1.l);
    i__5 = zvod01_1.n;
    for (i__ = 1; i__ <= i__5; ++i__) {
/* L575: */
	i__1 = i__;
	i__2 = i__;
	i__3 = i__ + zvod01_1.lmax * yh_dim1;
	z__2.r = cnquot * yh[i__3].r, z__2.i = cnquot * yh[i__3].i;
	z__1.r = acor[i__2].r - z__2.r, z__1.i = acor[i__2].i - z__2.i;
	savf[i__1].r = z__1.r, savf[i__1].i = z__1.i;
    }
    dup = zvnorm_(&zvod01_1.n, &savf[1], &ewt[1]) / zvod01_1.tq[2];
    d__1 = bias3 * dup;
    d__2 = one / (flotl + one);
    etaqp1 = one / (pow_dd(&d__1, &d__2) + addon);
L580:
    if (etaq >= etaqp1) {
	goto L590;
    }
    if (etaqp1 > etaqm1) {
	goto L620;
    }
    goto L610;
L590:
    if (etaq < etaqm1) {
	goto L610;
    }
L600:
    zvod01_1.eta = etaq;
    zvod01_1.newq = zvod01_1.nq;
    goto L630;
L610:
    zvod01_1.eta = etaqm1;
    zvod01_1.newq = zvod01_1.nq - 1;
    goto L630;
L620:
    zvod01_1.eta = etaqp1;
    zvod01_1.newq = zvod01_1.nq + 1;
    zcopy_(&zvod01_1.n, &acor[1], &c__1, &yh[zvod01_1.lmax * yh_dim1 + 1], &
	    c__1);
/* Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ---- */
L630:
    if (zvod01_1.eta < thresh || zvod01_1.etamax == one) {
	goto L640;
    }
    zvod01_1.eta = min(zvod01_1.eta,zvod01_1.etamax);
/* Computing MAX */
    d__1 = one, d__2 = abs(zvod01_1.h__) * zvod01_1.hmxi * zvod01_1.eta;
    zvod01_1.eta /= max(d__1,d__2);
    zvod01_1.newh = 1;
    zvod01_1.hnew = zvod01_1.h__ * zvod01_1.eta;
    goto L690;
L640:
    zvod01_1.newq = zvod01_1.nq;
    zvod01_1.newh = 0;
    zvod01_1.eta = one;
    zvod01_1.hnew = zvod01_1.h__;
    goto L690;
/* ----------------------------------------------------------------------- */
/* All returns are made through this section. */
/* On a successful return, ETAMAX is reset and ACOR is scaled. */
/* ----------------------------------------------------------------------- */
L660:
    zvod01_1.kflag = -1;
    goto L720;
L670:
    zvod01_1.kflag = -2;
    goto L720;
L680:
    if (nflag == -2) {
	zvod01_1.kflag = -3;
    }
    if (nflag == -3) {
	zvod01_1.kflag = -4;
    }
    goto L720;
L690:
    zvod01_1.etamax = etamx3;
    if (zvod02_1.nst <= 10) {
	zvod01_1.etamax = etamx2;
    }
/* L700: */
    r__ = one / zvod01_1.tq[1];
    dzscal_(&zvod01_1.n, &r__, &acor[1], &c__1);
L720:
    zvod01_1.jstart = 1;
    return 0;
/* ----------------------- End of Subroutine ZVSTEP ---------------------- */
} /* zvstep_ */

/* DECK ZVSET */
/* Subroutine */ int zvset_(void)
{
    /* Initialized data */

    static doublereal cortes = .1;
    static doublereal one = 1.;
    static doublereal six = 6.;
    static doublereal two = 2.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal s, t1, t2, t3, t4, t5, t6, em[13], xi, em0;
    static integer jp1;
    static doublereal elp, rxi;
    static integer nqm1, nqm2;
    static doublereal csum, hsum, rxis, alph0, cnqm1;
    static integer iback;
    static doublereal floti, flotl, ahatn0, flotnq;

/* ----------------------------------------------------------------------- */
/* Call sequence communication: None */
/* COMMON block variables accessed: */
/*     /ZVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1), */
/*                 METH, NQ, NQWAIT */

/* Subroutines called by ZVSET: None */
/* Function routines called by ZVSET: None */
/* ----------------------------------------------------------------------- */
/* ZVSET is called by ZVSTEP and sets coefficients for use there. */

/* For each order NQ, the coefficients in EL are calculated by use of */
/*  the generating polynomial lambda(x), with coefficients EL(i). */
/*      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ). */
/* For the backward differentiation formulas, */
/*                                     NQ-1 */
/*      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) . */
/*                                     i = 1 */
/* For the Adams formulas, */
/*                              NQ-1 */
/*      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) , */
/*                              i = 1 */
/*      lambda(-1) = 0,    lambda(0) = 1, */
/* where c is a normalization constant. */
/* In both cases, xi(i) is defined by */
/*      H*xi(i) = t sub n  -  t sub (n-i) */
/*              = H + TAU(1) + TAU(2) + ... TAU(i-1). */


/* In addition to variables described previously, communication */
/* with ZVSET uses the following: */
/*   TAU    = A vector of length 13 containing the past NQ values */
/*            of H. */
/*   EL     = A vector of length 13 in which vset stores the */
/*            coefficients for the corrector formula. */
/*   TQ     = A vector of length 5 in which vset stores constants */
/*            used for the convergence test, the error test, and the */
/*            selection of H at a new order. */
/*   METH   = The basic method indicator. */
/*   NQ     = The current order. */
/*   L      = NQ + 1, the length of the vector stored in EL, and */
/*            the number of columns of the YH array being used. */
/*   NQWAIT = A counter controlling the frequency of order changes. */
/*            An order change is about to be considered if NQWAIT = 1. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */



    flotl = (real) zvod01_1.l;
    nqm1 = zvod01_1.nq - 1;
    nqm2 = zvod01_1.nq - 2;
    switch (zvod01_1.meth) {
	case 1:  goto L100;
	case 2:  goto L200;
    }

/* Set coefficients for Adams methods. ---------------------------------- */
L100:
    if (zvod01_1.nq != 1) {
	goto L110;
    }
    zvod01_1.el[0] = one;
    zvod01_1.el[1] = one;
    zvod01_1.tq[0] = one;
    zvod01_1.tq[1] = two;
    zvod01_1.tq[2] = six * zvod01_1.tq[1];
    zvod01_1.tq[4] = one;
    goto L300;
L110:
    hsum = zvod01_1.h__;
    em[0] = one;
    flotnq = flotl - one;
    i__1 = zvod01_1.l;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L115: */
	em[i__ - 1] = zero;
    }
    i__1 = nqm1;
    for (j = 1; j <= i__1; ++j) {
	if (j != nqm1 || zvod01_1.nqwait != 1) {
	    goto L130;
	}
	s = one;
	csum = zero;
	i__2 = nqm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    csum += s * em[i__ - 1] / (real) (i__ + 1);
/* L120: */
	    s = -s;
	}
	zvod01_1.tq[0] = em[nqm1 - 1] / (flotnq * csum);
L130:
	rxi = zvod01_1.h__ / hsum;
	i__2 = j;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 2 - iback;
/* L140: */
	    em[i__ - 1] += em[i__ - 2] * rxi;
	}
	hsum += zvod01_1.tau[j - 1];
/* L150: */
    }
/* Compute integral from -1 to 0 of polynomial and of x times it. ------- */
    s = one;
    em0 = zero;
    csum = zero;
    i__1 = zvod01_1.nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	floti = (real) i__;
	em0 += s * em[i__ - 1] / floti;
	csum += s * em[i__ - 1] / (floti + one);
/* L160: */
	s = -s;
    }
/* In EL, form coefficients of normalized integrated polynomial. -------- */
    s = one / em0;
    zvod01_1.el[0] = one;
    i__1 = zvod01_1.nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
	zvod01_1.el[i__] = s * em[i__ - 1] / (real) i__;
    }
    xi = hsum / zvod01_1.h__;
    zvod01_1.tq[1] = xi * em0 / csum;
    zvod01_1.tq[4] = xi / zvod01_1.el[zvod01_1.l - 1];
    if (zvod01_1.nqwait != 1) {
	goto L300;
    }
/* For higher order control constant, multiply polynomial by 1+x/xi(q). - */
    rxi = one / xi;
    i__1 = zvod01_1.nq;
    for (iback = 1; iback <= i__1; ++iback) {
	i__ = zvod01_1.l + 1 - iback;
/* L180: */
	em[i__ - 1] += em[i__ - 2] * rxi;
    }
/* Compute integral of polynomial. -------------------------------------- */
    s = one;
    csum = zero;
    i__1 = zvod01_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	csum += s * em[i__ - 1] / (real) (i__ + 1);
/* L190: */
	s = -s;
    }
    zvod01_1.tq[2] = flotl * em0 / csum;
    goto L300;

/* Set coefficients for BDF methods. ------------------------------------ */
L200:
    i__1 = zvod01_1.l;
    for (i__ = 3; i__ <= i__1; ++i__) {
/* L210: */
	zvod01_1.el[i__ - 1] = zero;
    }
    zvod01_1.el[0] = one;
    zvod01_1.el[1] = one;
    alph0 = -one;
    ahatn0 = -one;
    hsum = zvod01_1.h__;
    rxi = one;
    rxis = one;
    if (zvod01_1.nq == 1) {
	goto L240;
    }
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------ */
	hsum += zvod01_1.tau[j - 1];
	rxi = zvod01_1.h__ / hsum;
	jp1 = j + 1;
	alph0 -= one / (real) jp1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 3 - iback;
/* L220: */
	    zvod01_1.el[i__ - 1] += zvod01_1.el[i__ - 2] * rxi;
	}
/* L230: */
    }
    alph0 -= one / (real) zvod01_1.nq;
    rxis = -zvod01_1.el[1] - alph0;
    hsum += zvod01_1.tau[nqm1 - 1];
    rxi = zvod01_1.h__ / hsum;
    ahatn0 = -zvod01_1.el[1] - rxi;
    i__1 = zvod01_1.nq;
    for (iback = 1; iback <= i__1; ++iback) {
	i__ = zvod01_1.nq + 2 - iback;
/* L235: */
	zvod01_1.el[i__ - 1] += zvod01_1.el[i__ - 2] * rxis;
    }
L240:
    t1 = one - ahatn0 + alph0;
    t2 = one + (real) zvod01_1.nq * t1;
    zvod01_1.tq[1] = (d__1 = alph0 * t2 / t1, abs(d__1));
    zvod01_1.tq[4] = (d__1 = t2 / (zvod01_1.el[zvod01_1.l - 1] * rxi / rxis), 
	    abs(d__1));
    if (zvod01_1.nqwait != 1) {
	goto L300;
    }
    cnqm1 = rxis / zvod01_1.el[zvod01_1.l - 1];
    t3 = alph0 + one / (real) zvod01_1.nq;
    t4 = ahatn0 + rxi;
    elp = t3 / (one - t4 + t3);
    zvod01_1.tq[0] = (d__1 = elp / cnqm1, abs(d__1));
    hsum += zvod01_1.tau[zvod01_1.nq - 1];
    rxi = zvod01_1.h__ / hsum;
    t5 = alph0 - one / (real) (zvod01_1.nq + 1);
    t6 = ahatn0 - rxi;
    elp = t2 / (one - t6 + t5);
    zvod01_1.tq[2] = (d__1 = elp * rxi * (flotl + one) * t5, abs(d__1));
L300:
    zvod01_1.tq[3] = cortes * zvod01_1.tq[1];
    return 0;
/* ----------------------- End of Subroutine ZVSET ----------------------- */
} /* zvset_ */

/* DECK ZVJUST */
/* Subroutine */ int zvjust_(doublecomplex *yh, integer *ldyh, integer *iord)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j;
    static doublereal t1, xi;
    static integer jp1, lp1, nqm1, nqm2, nqp1;
    static doublereal prod, hsum, alph0, alph1;
    static integer iback;
    static doublereal xiold;
    extern /* Subroutine */ int dzaxpy_(integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- YH, LDYH, IORD */
/* Call sequence output -- YH */
/* COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N */
/* COMMON block variables accessed: */
/*     /ZVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ, */

/* Subroutines called by ZVJUST: DZAXPY */
/* Function routines called by ZVJUST: None */
/* ----------------------------------------------------------------------- */
/* This subroutine adjusts the YH array on reduction of order, */
/* and also when the order is increased for the stiff option (METH = 2). */
/* Communication with ZVJUST uses the following: */
/* IORD  = An integer flag used when METH = 2 to indicate an order */
/*         increase (IORD = +1) or an order decrease (IORD = -1). */
/* HSCAL = Step size H used in scaling of Nordsieck array YH. */
/*         (If IORD = +1, ZVJUST assumes that HSCAL = TAU(1).) */
/* See References 1 and 2 for details. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;

    /* Function Body */

    if (zvod01_1.nq == 2 && *iord != 1) {
	return 0;
    }
    nqm1 = zvod01_1.nq - 1;
    nqm2 = zvod01_1.nq - 2;
    switch (zvod01_1.meth) {
	case 1:  goto L100;
	case 2:  goto L200;
    }
/* ----------------------------------------------------------------------- */
/* Nonstiff option... */
/* Check to see if the order is being increased or decreased. */
/* ----------------------------------------------------------------------- */
L100:
    if (*iord == 1) {
	goto L180;
    }
/* Order decrease. ------------------------------------------------------ */
    i__1 = zvod01_1.lmax;
    for (j = 1; j <= i__1; ++j) {
/* L110: */
	zvod01_1.el[j - 1] = zero;
    }
    zvod01_1.el[1] = one;
    hsum = zero;
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). ----------------- */
	hsum += zvod01_1.tau[j - 1];
	xi = hsum / zvod01_1.hscal;
	jp1 = j + 1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 3 - iback;
/* L120: */
	    zvod01_1.el[i__ - 1] = zvod01_1.el[i__ - 1] * xi + zvod01_1.el[
		    i__ - 2];
	}
/* L130: */
    }
/* Construct coefficients of integrated polynomial. --------------------- */
    i__1 = nqm1;
    for (j = 2; j <= i__1; ++j) {
/* L140: */
	zvod01_1.el[j] = (real) zvod01_1.nq * zvod01_1.el[j - 1] / (real) j;
    }
/* Subtract correction terms from YH array. ----------------------------- */
    i__1 = zvod01_1.nq;
    for (j = 3; j <= i__1; ++j) {
	i__2 = zvod01_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L160: */
	    i__3 = i__ + j * yh_dim1;
	    i__4 = i__ + j * yh_dim1;
	    i__5 = i__ + zvod01_1.l * yh_dim1;
	    i__6 = j - 1;
	    z__2.r = zvod01_1.el[i__6] * yh[i__5].r, z__2.i = zvod01_1.el[
		    i__6] * yh[i__5].i;
	    z__1.r = yh[i__4].r - z__2.r, z__1.i = yh[i__4].i - z__2.i;
	    yh[i__3].r = z__1.r, yh[i__3].i = z__1.i;
	}
/* L170: */
    }
    return 0;
/* Order increase. ------------------------------------------------------ */
/* Zero out next column in YH array. ------------------------------------ */
L180:
    lp1 = zvod01_1.l + 1;
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L190: */
	i__3 = i__ + lp1 * yh_dim1;
	yh[i__3].r = zero, yh[i__3].i = 0.;
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* Stiff option... */
/* Check to see if the order is being increased or decreased. */
/* ----------------------------------------------------------------------- */
L200:
    if (*iord == 1) {
	goto L300;
    }
/* Order decrease. ------------------------------------------------------ */
    i__3 = zvod01_1.lmax;
    for (j = 1; j <= i__3; ++j) {
/* L210: */
	zvod01_1.el[j - 1] = zero;
    }
    zvod01_1.el[2] = one;
    hsum = zero;
    i__3 = nqm2;
    for (j = 1; j <= i__3; ++j) {
/* Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). --------------- */
	hsum += zvod01_1.tau[j - 1];
	xi = hsum / zvod01_1.hscal;
	jp1 = j + 1;
	i__1 = jp1;
	for (iback = 1; iback <= i__1; ++iback) {
	    i__ = j + 4 - iback;
/* L220: */
	    zvod01_1.el[i__ - 1] = zvod01_1.el[i__ - 1] * xi + zvod01_1.el[
		    i__ - 2];
	}
/* L230: */
    }
/* Subtract correction terms from YH array. ----------------------------- */
    i__3 = zvod01_1.nq;
    for (j = 3; j <= i__3; ++j) {
	i__1 = zvod01_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L240: */
	    i__4 = i__ + j * yh_dim1;
	    i__5 = i__ + j * yh_dim1;
	    i__6 = i__ + zvod01_1.l * yh_dim1;
	    i__2 = j - 1;
	    z__2.r = zvod01_1.el[i__2] * yh[i__6].r, z__2.i = zvod01_1.el[
		    i__2] * yh[i__6].i;
	    z__1.r = yh[i__5].r - z__2.r, z__1.i = yh[i__5].i - z__2.i;
	    yh[i__4].r = z__1.r, yh[i__4].i = z__1.i;
	}
/* L250: */
    }
    return 0;
/* Order increase. ------------------------------------------------------ */
L300:
    i__3 = zvod01_1.lmax;
    for (j = 1; j <= i__3; ++j) {
/* L310: */
	zvod01_1.el[j - 1] = zero;
    }
    zvod01_1.el[2] = one;
    alph0 = -one;
    alph1 = one;
    prod = one;
    xiold = one;
    hsum = zvod01_1.hscal;
    if (zvod01_1.nq == 1) {
	goto L340;
    }
    i__3 = nqm1;
    for (j = 1; j <= i__3; ++j) {
/* Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). --------------- */
	jp1 = j + 1;
	hsum += zvod01_1.tau[jp1 - 1];
	xi = hsum / zvod01_1.hscal;
	prod *= xi;
	alph0 -= one / (real) jp1;
	alph1 += one / xi;
	i__4 = jp1;
	for (iback = 1; iback <= i__4; ++iback) {
	    i__ = j + 4 - iback;
/* L320: */
	    zvod01_1.el[i__ - 1] = zvod01_1.el[i__ - 1] * xiold + zvod01_1.el[
		    i__ - 2];
	}
	xiold = xi;
/* L330: */
    }
L340:
    t1 = (-alph0 - alph1) / prod;
/* Load column L + 1 in YH array. --------------------------------------- */
    lp1 = zvod01_1.l + 1;
    i__3 = zvod01_1.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L350: */
	i__4 = i__ + lp1 * yh_dim1;
	i__5 = i__ + zvod01_1.lmax * yh_dim1;
	z__1.r = t1 * yh[i__5].r, z__1.i = t1 * yh[i__5].i;
	yh[i__4].r = z__1.r, yh[i__4].i = z__1.i;
    }
/* Add correction terms to YH array. ------------------------------------ */
    nqp1 = zvod01_1.nq + 1;
    i__4 = nqp1;
    for (j = 3; j <= i__4; ++j) {
	dzaxpy_(&zvod01_1.n, &zvod01_1.el[j - 1], &yh[lp1 * yh_dim1 + 1], &
		c__1, &yh[j * yh_dim1 + 1], &c__1);
/* L370: */
    }
    return 0;
/* ----------------------- End of Subroutine ZVJUST ---------------------- */
} /* zvjust_ */

/* DECK ZVNLSD */
/* Subroutine */ int zvnlsd_(doublecomplex *y, doublecomplex *yh, integer *
	ldyh, doublecomplex *vsav, doublecomplex *savf, doublereal *ewt, 
	doublecomplex *acor, integer *iwm, doublecomplex *wm, S_fp f, U_fp 
	jac, U_fp pdum, integer *nflag, real *rpar, integer *ipar)
{
    /* Initialized data */

    static doublereal ccmax = .3;
    static doublereal crdown = .3;
    static integer maxcor = 3;
    static integer msbp = 20;
    static doublereal rdiv = 2.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    static integer i__, m;
    static doublereal del, dcon, delp;
    static integer ierpj;
    extern /* Subroutine */ int zvjac_(doublecomplex *, doublecomplex *, 
	    integer *, doublereal *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, S_fp, U_fp, integer *, real *, 
	    integer *);
    static integer iersl;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zvsol_(doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal cscale;
    extern /* Subroutine */ int dzscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), dzaxpy_(integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal zvnorm_(integer *, doublecomplex *, doublereal *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM, */
/*                        F, JAC, NFLAG, RPAR, IPAR */
/* Call sequence output -- YH, ACOR, WM, IWM, NFLAG */
/* COMMON block variables accessed: */
/*     /ZVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF, */
/*                JCUR, METH, MITER, N, NSLP */
/*     /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST */

/* Subroutines called by ZVNLSD: F, DZAXPY, ZCOPY, DZSCAL, ZVJAC, ZVSOL */
/* Function routines called by ZVNLSD: ZVNORM */
/* ----------------------------------------------------------------------- */
/* Subroutine ZVNLSD is a nonlinear system solver, which uses functional */
/* iteration or a chord (modified Newton) method.  For the chord method */
/* direct linear algebraic system solvers are used.  Subroutine ZVNLSD */
/* then handles the corrector phase of this integration package. */

/* Communication with ZVNLSD is done with the following variables. (For */
/* more details, please see the comments in the driver subroutine.) */

/* Y          = The dependent variable, a vector of length N, input. */
/* YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input */
/*              and output.  On input, it contains predicted values. */
/* LDYH       = A constant .ge. N, the first dimension of YH, input. */
/* VSAV       = Unused work array. */
/* SAVF       = A work array of length N. */
/* EWT        = An error weight vector of length N, input. */
/* ACOR       = A work array of length N, used for the accumulated */
/*              corrections to the predicted y vector. */
/* WM,IWM     = Complex and integer work arrays associated with matrix */
/*              operations in chord iteration (MITER .ne. 0). */
/* F          = Dummy name for user-supplied routine for f. */
/* JAC        = Dummy name for user-supplied Jacobian routine. */
/* PDUM       = Unused dummy subroutine name.  Included for uniformity */
/*              over collection of integrators. */
/* NFLAG      = Input/output flag, with values and meanings as follows: */
/*              INPUT */
/*                  0 first call for this time step. */
/*                 -1 convergence failure in previous call to ZVNLSD. */
/*                 -2 error test failure in ZVSTEP. */
/*              OUTPUT */
/*                  0 successful completion of nonlinear solver. */
/*                 -1 convergence failure or singular matrix. */
/*                 -2 unrecoverable error in matrix preprocessing */
/*                    (cannot occur here). */
/*                 -3 unrecoverable error in solution (cannot occur */
/*                    here). */
/* RPAR, IPAR = User's real/complex and integer work arrays. */

/* IPUP       = Own variable flag with values and meanings as follows: */
/*              0,            do not update the Newton matrix. */
/*              MITER .ne. 0, update Newton matrix, because it is the */
/*                            initial step, order was changed, the error */
/*                            test failed, or an update is indicated by */
/*                            the scalar RC or step counter NST. */

/* For more details, see comments in driver subroutine. */
/* ----------------------------------------------------------------------- */
/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for labeled COMMON block ZVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --y;
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --vsav;
    --savf;
    --ewt;
    --acor;
    --iwm;
    --wm;
    --rpar;
    --ipar;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* On the first step, on a change of method order, or after a */
/* nonlinear convergence failure with NFLAG = -2, set IPUP = MITER */
/* to force a Jacobian update when MITER .ne. 0. */
/* ----------------------------------------------------------------------- */
    if (zvod01_1.jstart == 0) {
	zvod01_1.nslp = 0;
    }
    if (*nflag == 0) {
	zvod01_1.icf = 0;
    }
    if (*nflag == -2) {
	zvod01_1.ipup = zvod01_1.miter;
    }
    if (zvod01_1.jstart == 0 || zvod01_1.jstart == -1) {
	zvod01_1.ipup = zvod01_1.miter;
    }
/* If this is functional iteration, set CRATE .eq. 1 and drop to 220 */
    if (zvod01_1.miter == 0) {
	zvod01_1.crate = one;
	goto L220;
    }
/* ----------------------------------------------------------------------- */
/* RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1. */
/* When RC differs from 1 by more than CCMAX, IPUP is set to MITER */
/* to force ZVJAC to be called, if a Jacobian is involved. */
/* In any case, ZVJAC is called at least every MSBP steps. */
/* ----------------------------------------------------------------------- */
    zvod01_1.drc = (d__1 = zvod01_1.rc - one, abs(d__1));
    if (zvod01_1.drc > ccmax || zvod02_1.nst >= zvod01_1.nslp + msbp) {
	zvod01_1.ipup = zvod01_1.miter;
    }
/* ----------------------------------------------------------------------- */
/* Up to MAXCOR corrector iterations are taken.  A convergence test is */
/* made on the r.m.s. norm of each correction, weighted by the error */
/* weight vector EWT.  The sum of the corrections is accumulated in the */
/* vector ACOR(i).  The YH array is not altered in the corrector loop. */
/* ----------------------------------------------------------------------- */
L220:
    m = 0;
    delp = zero;
    zcopy_(&zvod01_1.n, &yh[yh_dim1 + 1], &c__1, &y[1], &c__1);
    (*f)(&zvod01_1.n, &zvod01_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++zvod02_1.nfe;
    if (zvod01_1.ipup <= 0) {
	goto L250;
    }
/* ----------------------------------------------------------------------- */
/* If indicated, the matrix P = I - h*rl1*J is reevaluated and */
/* preprocessed before starting the corrector iteration.  IPUP is set */
/* to 0 as an indicator that this has been done. */
/* ----------------------------------------------------------------------- */
    zvjac_(&y[1], &yh[yh_offset], ldyh, &ewt[1], &acor[1], &savf[1], &wm[1], &
	    iwm[1], (S_fp)f, (U_fp)jac, &ierpj, &rpar[1], &ipar[1]);
    zvod01_1.ipup = 0;
    zvod01_1.rc = one;
    zvod01_1.drc = zero;
    zvod01_1.crate = one;
    zvod01_1.nslp = zvod02_1.nst;
/* If matrix is singular, take error return to force cut in step size. -- */
    if (ierpj != 0) {
	goto L430;
    }
L250:
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L260: */
	i__2 = i__;
	acor[i__2].r = zero, acor[i__2].i = 0.;
    }
/* This is a looping point for the corrector iteration. ----------------- */
L270:
    if (zvod01_1.miter != 0) {
	goto L350;
    }
/* ----------------------------------------------------------------------- */
/* In the case of functional iteration, update Y directly from */
/* the result of the last function evaluation. */
/* ----------------------------------------------------------------------- */
    i__2 = zvod01_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L280: */
	i__1 = i__;
	i__3 = i__;
	z__3.r = zvod01_1.h__ * savf[i__3].r, z__3.i = zvod01_1.h__ * savf[
		i__3].i;
	i__4 = i__ + (yh_dim1 << 1);
	z__2.r = z__3.r - yh[i__4].r, z__2.i = z__3.i - yh[i__4].i;
	z__1.r = zvod01_1.rl1 * z__2.r, z__1.i = zvod01_1.rl1 * z__2.i;
	savf[i__1].r = z__1.r, savf[i__1].i = z__1.i;
    }
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L290: */
	i__3 = i__;
	i__4 = i__;
	i__2 = i__;
	z__1.r = savf[i__4].r - acor[i__2].r, z__1.i = savf[i__4].i - acor[
		i__2].i;
	y[i__3].r = z__1.r, y[i__3].i = z__1.i;
    }
    del = zvnorm_(&zvod01_1.n, &y[1], &ewt[1]);
    i__3 = zvod01_1.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L300: */
	i__4 = i__;
	i__2 = i__ + yh_dim1;
	i__1 = i__;
	z__1.r = yh[i__2].r + savf[i__1].r, z__1.i = yh[i__2].i + savf[i__1]
		.i;
	y[i__4].r = z__1.r, y[i__4].i = z__1.i;
    }
    zcopy_(&zvod01_1.n, &savf[1], &c__1, &acor[1], &c__1);
    goto L400;
/* ----------------------------------------------------------------------- */
/* In the case of the chord method, compute the corrector error, */
/* and solve the linear system with that as right-hand side and */
/* P as coefficient matrix.  The correction is scaled by the factor */
/* 2/(1+RC) to account for changes in h*rl1 since the last ZVJAC call. */
/* ----------------------------------------------------------------------- */
L350:
    i__4 = zvod01_1.n;
    for (i__ = 1; i__ <= i__4; ++i__) {
/* L360: */
	i__2 = i__;
	d__1 = zvod01_1.rl1 * zvod01_1.h__;
	i__1 = i__;
	z__2.r = d__1 * savf[i__1].r, z__2.i = d__1 * savf[i__1].i;
	i__3 = i__ + (yh_dim1 << 1);
	z__4.r = zvod01_1.rl1 * yh[i__3].r, z__4.i = zvod01_1.rl1 * yh[i__3]
		.i;
	i__5 = i__;
	z__3.r = z__4.r + acor[i__5].r, z__3.i = z__4.i + acor[i__5].i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	y[i__2].r = z__1.r, y[i__2].i = z__1.i;
    }
    zvsol_(&wm[1], &iwm[1], &y[1], &iersl);
    ++zvod02_1.nni;
    if (iersl > 0) {
	goto L410;
    }
    if (zvod01_1.meth == 2 && zvod01_1.rc != one) {
	cscale = two / (one + zvod01_1.rc);
	dzscal_(&zvod01_1.n, &cscale, &y[1], &c__1);
    }
    del = zvnorm_(&zvod01_1.n, &y[1], &ewt[1]);
    dzaxpy_(&zvod01_1.n, &one, &y[1], &c__1, &acor[1], &c__1);
    i__2 = zvod01_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L380: */
	i__1 = i__;
	i__3 = i__ + yh_dim1;
	i__5 = i__;
	z__1.r = yh[i__3].r + acor[i__5].r, z__1.i = yh[i__3].i + acor[i__5]
		.i;
	y[i__1].r = z__1.r, y[i__1].i = z__1.i;
    }
/* ----------------------------------------------------------------------- */
/* Test for convergence.  If M .gt. 0, an estimate of the convergence */
/* rate constant is stored in CRATE, and this is used in the test. */
/* ----------------------------------------------------------------------- */
L400:
    if (m != 0) {
/* Computing MAX */
	d__1 = crdown * zvod01_1.crate, d__2 = del / delp;
	zvod01_1.crate = max(d__1,d__2);
    }
    dcon = del * min(one,zvod01_1.crate) / zvod01_1.tq[3];
    if (dcon <= one) {
	goto L450;
    }
    ++m;
    if (m == maxcor) {
	goto L410;
    }
    if (m >= 2 && del > rdiv * delp) {
	goto L410;
    }
    delp = del;
    (*f)(&zvod01_1.n, &zvod01_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++zvod02_1.nfe;
    goto L270;

L410:
    if (zvod01_1.miter == 0 || zvod01_1.jcur == 1) {
	goto L430;
    }
    zvod01_1.icf = 1;
    zvod01_1.ipup = zvod01_1.miter;
    goto L220;

L430:
    *nflag = -1;
    zvod01_1.icf = 2;
    zvod01_1.ipup = zvod01_1.miter;
    return 0;

/* Return for successful step. ------------------------------------------ */
L450:
    *nflag = 0;
    zvod01_1.jcur = 0;
    zvod01_1.icf = 0;
    if (m == 0) {
	zvod01_1.acnrm = del;
    }
    if (m > 0) {
	zvod01_1.acnrm = zvnorm_(&zvod01_1.n, &acor[1], &ewt[1]);
    }
    return 0;
/* ----------------------- End of Subroutine ZVNLSD ---------------------- */
} /* zvnlsd_ */

/* DECK ZVJAC */
/* Subroutine */ int zvjac_(doublecomplex *y, doublecomplex *yh, integer *
	ldyh, doublereal *ewt, doublecomplex *ftem, doublecomplex *savf, 
	doublecomplex *wm, integer *iwm, S_fp f, S_fp jac, integer *ierpj, 
	real *rpar, integer *ipar)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal thou = 1e3;
    static doublereal zero = 0.;
    static doublereal pt1 = .1;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal r__;
    static integer i1, i2, j1;
    static doublereal r0;
    static doublecomplex r1, di;
    static integer ii, jj, ml;
    static doublecomplex yi, yj;
    static integer mu, ml1, np1;
    static doublereal fac;
    static integer mba;
    static doublereal con;
    static integer ier, jok;
    static doublecomplex yjj;
    static integer meb1, lenp, mband;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer meband;
    extern /* Subroutine */ int dzscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), zgbtrf_(integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    integer *), zgetrf_(integer *, integer *, doublecomplex *, 
	    integer *, integer *, integer *), zacopy_(integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal zvnorm_(integer *, doublecomplex *, doublereal *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, */
/*                        F, JAC, RPAR, IPAR */
/* Call sequence output -- WM, IWM, IERPJ */
/* COMMON block variables accessed: */
/*     /ZVOD01/  CCMXJ, DRC, H, HRL1, RL1, SRUR, TN, UROUND, ICF, JCUR, */
/*               LOCJS, MITER, MSBJ, N, NSLJ */
/*     /ZVOD02/  NFE, NST, NJE, NLU */

/* Subroutines called by ZVJAC: F, JAC, ZACOPY, ZCOPY, ZGBTRF, ZGETRF, */
/*                              DZSCAL */
/* Function routines called by ZVJAC: ZVNORM */
/* ----------------------------------------------------------------------- */
/* ZVJAC is called by ZVNLSD to compute and process the matrix */
/* P = I - h*rl1*J , where J is an approximation to the Jacobian. */
/* Here J is computed by the user-supplied routine JAC if */
/* MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5. */
/* If MITER = 3, a diagonal approximation to J is used. */
/* If JSV = -1, J is computed from scratch in all cases. */
/* If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is */
/* considered acceptable, then P is constructed from the saved J. */
/* J is stored in wm and replaced by P.  If MITER .ne. 3, P is then */
/* subjected to LU decomposition in preparation for later solution */
/* of linear systems with P as coefficient matrix. This is done */
/* by ZGETRF if MITER = 1 or 2, and by ZGBTRF if MITER = 4 or 5. */

/* Communication with ZVJAC is done with the following variables.  (For */
/* more details, please see the comments in the driver subroutine.) */
/* Y          = Vector containing predicted values on entry. */
/* YH         = The Nordsieck array, an LDYH by LMAX array, input. */
/* LDYH       = A constant .ge. N, the first dimension of YH, input. */
/* EWT        = An error weight vector of length N. */
/* SAVF       = Array containing f evaluated at predicted y, input. */
/* WM         = Complex work space for matrices.  In the output, it */
/*              contains the inverse diagonal matrix if MITER = 3 and */
/*              the LU decomposition of P if MITER is 1, 2 , 4, or 5. */
/*              Storage of the saved Jacobian starts at WM(LOCJS). */
/* IWM        = Integer work space containing pivot information, */
/*              starting at IWM(31), if MITER is 1, 2, 4, or 5. */
/*              IWM also contains band parameters ML = IWM(1) and */
/*              MU = IWM(2) if MITER is 4 or 5. */
/* F          = Dummy name for the user-supplied subroutine for f. */
/* JAC        = Dummy name for the user-supplied Jacobian subroutine. */
/* RPAR, IPAR = User's real/complex and integer work arrays. */
/* RL1        = 1/EL(2) (input). */
/* IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P */
/*              matrix is found to be singular. */
/* JCUR       = Output flag to indicate whether the Jacobian matrix */
/*              (or approximation) is now current. */
/*              JCUR = 0 means J is not current. */
/*              JCUR = 1 means J is current. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for labeled COMMON block ZVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this subroutine. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --ewt;
    --ftem;
    --savf;
    --wm;
    --iwm;
    --rpar;
    --ipar;

    /* Function Body */

    *ierpj = 0;
    zvod01_1.hrl1 = zvod01_1.h__ * zvod01_1.rl1;
/* See whether J should be evaluated (JOK = -1) or not (JOK = 1). ------- */
    jok = zvod01_1.jsv;
    if (zvod01_1.jsv == 1) {
	if (zvod02_1.nst == 0 || zvod02_1.nst > zvod01_1.nslj + zvod01_1.msbj)
		 {
	    jok = -1;
	}
	if (zvod01_1.icf == 1 && zvod01_1.drc < zvod01_1.ccmxj) {
	    jok = -1;
	}
	if (zvod01_1.icf == 2) {
	    jok = -1;
	}
    }
/* End of setting JOK. -------------------------------------------------- */

    if (jok == -1 && zvod01_1.miter == 1) {
/* If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------ */
	++zvod02_1.nje;
	zvod01_1.nslj = zvod02_1.nst;
	zvod01_1.jcur = 1;
	lenp = zvod01_1.n * zvod01_1.n;
	i__1 = lenp;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	    i__2 = i__;
	    wm[i__2].r = zero, wm[i__2].i = 0.;
	}
	(*jac)(&zvod01_1.n, &zvod01_1.tn, &y[1], &c__0, &c__0, &wm[1], &
		zvod01_1.n, &rpar[1], &ipar[1]);
	if (zvod01_1.jsv == 1) {
	    zcopy_(&lenp, &wm[1], &c__1, &wm[zvod01_1.locjs], &c__1);
	}
    }

    if (jok == -1 && zvod01_1.miter == 2) {
/* If MITER = 2, make N calls to F to approximate the Jacobian. --------- */
	++zvod02_1.nje;
	zvod01_1.nslj = zvod02_1.nst;
	zvod01_1.jcur = 1;
	fac = zvnorm_(&zvod01_1.n, &savf[1], &ewt[1]);
	r0 = thou * abs(zvod01_1.h__) * zvod01_1.uround * (real) zvod01_1.n * 
		fac;
	if (r0 == zero) {
	    r0 = one;
	}
	j1 = 0;
	i__2 = zvod01_1.n;
	for (j = 1; j <= i__2; ++j) {
	    i__1 = j;
	    yj.r = y[i__1].r, yj.i = y[i__1].i;
/* Computing MAX */
	    d__1 = zvod01_1.srur * z_abs(&yj), d__2 = r0 / ewt[j];
	    r__ = max(d__1,d__2);
	    i__1 = j;
	    i__3 = j;
	    z__1.r = y[i__3].r + r__, z__1.i = y[i__3].i;
	    y[i__1].r = z__1.r, y[i__1].i = z__1.i;
	    fac = one / r__;
	    (*f)(&zvod01_1.n, &zvod01_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[
		    1]);
	    i__1 = zvod01_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L220: */
		i__3 = i__ + j1;
		i__4 = i__;
		i__5 = i__;
		z__2.r = ftem[i__4].r - savf[i__5].r, z__2.i = ftem[i__4].i - 
			savf[i__5].i;
		z__1.r = fac * z__2.r, z__1.i = fac * z__2.i;
		wm[i__3].r = z__1.r, wm[i__3].i = z__1.i;
	    }
	    i__3 = j;
	    y[i__3].r = yj.r, y[i__3].i = yj.i;
	    j1 += zvod01_1.n;
/* L230: */
	}
	zvod02_1.nfe += zvod01_1.n;
	lenp = zvod01_1.n * zvod01_1.n;
	if (zvod01_1.jsv == 1) {
	    zcopy_(&lenp, &wm[1], &c__1, &wm[zvod01_1.locjs], &c__1);
	}
    }

    if (jok == 1 && (zvod01_1.miter == 1 || zvod01_1.miter == 2)) {
	zvod01_1.jcur = 0;
	lenp = zvod01_1.n * zvod01_1.n;
	zcopy_(&lenp, &wm[zvod01_1.locjs], &c__1, &wm[1], &c__1);
    }

    if (zvod01_1.miter == 1 || zvod01_1.miter == 2) {
/* Multiply Jacobian by scalar, add identity, and do LU decomposition. -- */
	con = -zvod01_1.hrl1;
	dzscal_(&lenp, &con, &wm[1], &c__1);
	j = 1;
	np1 = zvod01_1.n + 1;
	i__2 = zvod01_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = j;
	    i__4 = j;
	    z__1.r = wm[i__4].r + one, z__1.i = wm[i__4].i;
	    wm[i__3].r = z__1.r, wm[i__3].i = z__1.i;
/* L250: */
	    j += np1;
	}
	++zvod02_1.nlu;
/*     Replaced LINPACK zgefa with LAPACK zgetrf */
/*      CALL ZGEFA (WM, N, N, IWM(31), IER) */
	zgetrf_(&zvod01_1.n, &zvod01_1.n, &wm[1], &zvod01_1.n, &iwm[31], &ier)
		;
	if (ier != 0) {
	    *ierpj = 1;
	}
	return 0;
    }
/* End of code block for MITER = 1 or 2. -------------------------------- */

    if (zvod01_1.miter == 3) {
/* If MITER = 3, construct a diagonal approximation to J and P. --------- */
	++zvod02_1.nje;
	zvod01_1.jcur = 1;
	r__ = zvod01_1.rl1 * pt1;
	i__2 = zvod01_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L310: */
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    z__4.r = zvod01_1.h__ * savf[i__5].r, z__4.i = zvod01_1.h__ * 
		    savf[i__5].i;
	    i__1 = i__ + (yh_dim1 << 1);
	    z__3.r = z__4.r - yh[i__1].r, z__3.i = z__4.i - yh[i__1].i;
	    z__2.r = r__ * z__3.r, z__2.i = r__ * z__3.i;
	    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
	    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
	}
	(*f)(&zvod01_1.n, &zvod01_1.tn, &y[1], &wm[1], &rpar[1], &ipar[1]);
	++zvod02_1.nfe;
	i__3 = zvod01_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = i__;
	    z__2.r = zvod01_1.h__ * savf[i__4].r, z__2.i = zvod01_1.h__ * 
		    savf[i__4].i;
	    i__5 = i__ + (yh_dim1 << 1);
	    z__1.r = z__2.r - yh[i__5].r, z__1.i = z__2.i - yh[i__5].i;
	    r1.r = z__1.r, r1.i = z__1.i;
	    z__2.r = pt1 * r1.r, z__2.i = pt1 * r1.i;
	    i__4 = i__;
	    i__5 = i__;
	    z__4.r = wm[i__4].r - savf[i__5].r, z__4.i = wm[i__4].i - savf[
		    i__5].i;
	    z__3.r = zvod01_1.h__ * z__4.r, z__3.i = zvod01_1.h__ * z__4.i;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    di.r = z__1.r, di.i = z__1.i;
	    i__4 = i__;
	    wm[i__4].r = one, wm[i__4].i = 0.;
	    if (z_abs(&r1) < zvod01_1.uround / ewt[i__]) {
		goto L320;
	    }
	    if (z_abs(&di) == zero) {
		goto L330;
	    }
	    i__4 = i__;
	    z__2.r = pt1 * r1.r, z__2.i = pt1 * r1.i;
	    z_div(&z__1, &z__2, &di);
	    wm[i__4].r = z__1.r, wm[i__4].i = z__1.i;
L320:
	    ;
	}
	return 0;
L330:
	*ierpj = 1;
	return 0;
    }
/* End of code block for MITER = 3. ------------------------------------- */

/* Set constants for MITER = 4 or 5. ------------------------------------ */
    ml = iwm[1];
    mu = iwm[2];
    ml1 = ml + 1;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * zvod01_1.n;

    if (jok == -1 && zvod01_1.miter == 4) {
/* If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------ */
	++zvod02_1.nje;
	zvod01_1.nslj = zvod02_1.nst;
	zvod01_1.jcur = 1;
	i__3 = lenp;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L410: */
	    i__4 = i__;
	    wm[i__4].r = zero, wm[i__4].i = 0.;
	}
	(*jac)(&zvod01_1.n, &zvod01_1.tn, &y[1], &ml, &mu, &wm[ml1], &meband, 
		&rpar[1], &ipar[1]);
	if (zvod01_1.jsv == 1) {
	    zacopy_(&mband, &zvod01_1.n, &wm[ml1], &meband, &wm[
		    zvod01_1.locjs], &mband);
	}
    }

    if (jok == -1 && zvod01_1.miter == 5) {
/* If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. --- */
	++zvod02_1.nje;
	zvod01_1.nslj = zvod02_1.nst;
	zvod01_1.jcur = 1;
	mba = min(mband,zvod01_1.n);
	meb1 = meband - 1;
	fac = zvnorm_(&zvod01_1.n, &savf[1], &ewt[1]);
	r0 = thou * abs(zvod01_1.h__) * zvod01_1.uround * (real) zvod01_1.n * 
		fac;
	if (r0 == zero) {
	    r0 = one;
	}
	i__4 = mba;
	for (j = 1; j <= i__4; ++j) {
	    i__3 = zvod01_1.n;
	    i__5 = mband;
	    for (i__ = j; i__5 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__5) {
		i__1 = i__;
		yi.r = y[i__1].r, yi.i = y[i__1].i;
/* Computing MAX */
		d__1 = zvod01_1.srur * z_abs(&yi), d__2 = r0 / ewt[i__];
		r__ = max(d__1,d__2);
/* L530: */
		i__1 = i__;
		i__2 = i__;
		z__1.r = y[i__2].r + r__, z__1.i = y[i__2].i;
		y[i__1].r = z__1.r, y[i__1].i = z__1.i;
	    }
	    (*f)(&zvod01_1.n, &zvod01_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[
		    1]);
	    i__1 = zvod01_1.n;
	    i__2 = mband;
	    for (jj = j; i__2 < 0 ? jj >= i__1 : jj <= i__1; jj += i__2) {
		i__5 = jj;
		i__3 = jj + yh_dim1;
		y[i__5].r = yh[i__3].r, y[i__5].i = yh[i__3].i;
		i__5 = jj;
		yjj.r = y[i__5].r, yjj.i = y[i__5].i;
/* Computing MAX */
		d__1 = zvod01_1.srur * z_abs(&yjj), d__2 = r0 / ewt[jj];
		r__ = max(d__1,d__2);
		fac = one / r__;
/* Computing MAX */
		i__5 = jj - mu;
		i1 = max(i__5,1);
/* Computing MIN */
		i__5 = jj + ml;
		i2 = min(i__5,zvod01_1.n);
		ii = jj * meb1 - ml;
		i__5 = i2;
		for (i__ = i1; i__ <= i__5; ++i__) {
/* L540: */
		    i__3 = ii + i__;
		    i__6 = i__;
		    i__7 = i__;
		    z__2.r = ftem[i__6].r - savf[i__7].r, z__2.i = ftem[i__6]
			    .i - savf[i__7].i;
		    z__1.r = fac * z__2.r, z__1.i = fac * z__2.i;
		    wm[i__3].r = z__1.r, wm[i__3].i = z__1.i;
		}
/* L550: */
	    }
/* L560: */
	}
	zvod02_1.nfe += mba;
	if (zvod01_1.jsv == 1) {
	    zacopy_(&mband, &zvod01_1.n, &wm[ml1], &meband, &wm[
		    zvod01_1.locjs], &mband);
	}
    }

    if (jok == 1) {
	zvod01_1.jcur = 0;
	zacopy_(&mband, &zvod01_1.n, &wm[zvod01_1.locjs], &mband, &wm[ml1], &
		meband);
    }

/* Multiply Jacobian by scalar, add identity, and do LU decomposition. */
    con = -zvod01_1.hrl1;
    dzscal_(&lenp, &con, &wm[1], &c__1);
    ii = mband;
    i__4 = zvod01_1.n;
    for (i__ = 1; i__ <= i__4; ++i__) {
	i__2 = ii;
	i__1 = ii;
	z__1.r = wm[i__1].r + one, z__1.i = wm[i__1].i;
	wm[i__2].r = z__1.r, wm[i__2].i = z__1.i;
/* L580: */
	ii += meband;
    }
    ++zvod02_1.nlu;
/*     Replaced LINPACK zgbfa with LAPACK zgbtrf */
/*      CALL ZGBFA (WM, MEBAND, N, ML, MU, IWM(31), IER) */
    zgbtrf_(&zvod01_1.n, &zvod01_1.n, &ml, &mu, &wm[1], &meband, &iwm[31], &
	    ier);
    if (ier != 0) {
	*ierpj = 1;
    }
    return 0;
/* End of code block for MITER = 4 or 5. -------------------------------- */

/* ----------------------- End of Subroutine ZVJAC ----------------------- */
} /* zvjac_ */

/* DECK ZACOPY */
/* Subroutine */ int zacopy_(integer *nrow, integer *ncol, doublecomplex *a, 
	integer *nrowa, doublecomplex *b, integer *nrowb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static integer ic;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- NROW, NCOL, A, NROWA, NROWB */
/* Call sequence output -- B */
/* COMMON block variables accessed -- None */

/* Subroutines called by ZACOPY: ZCOPY */
/* Function routines called by ZACOPY: None */
/* ----------------------------------------------------------------------- */
/* This routine copies one rectangular array, A, to another, B, */
/* where A and B may have different row dimensions, NROWA and NROWB. */
/* The data copied consists of NROW rows and NCOL columns. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    a_dim1 = *nrowa;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *nrowb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    i__1 = *ncol;
    for (ic = 1; ic <= i__1; ++ic) {
	zcopy_(nrow, &a[ic * a_dim1 + 1], &c__1, &b[ic * b_dim1 + 1], &c__1);
/* L20: */
    }

    return 0;
/* ----------------------- End of Subroutine ZACOPY ---------------------- */
} /* zacopy_ */

/* DECK ZVSOL */
/* Subroutine */ int zvsol_(doublecomplex *wm, integer *iwm, doublecomplex *x,
	 integer *iersl)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublereal r__;
    static doublecomplex di;
    static integer ml, mu, ier;
    static doublereal phrl1;
    static integer meband;
    extern /* Subroutine */ int zgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen), zgetrs_(char *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- WM, IWM, X */
/* Call sequence output -- X, IERSL */
/* COMMON block variables accessed: */
/*     /ZVOD01/ -- H, HRL1, RL1, MITER, N */

/* Subroutines called by ZVSOL: ZGETRS, ZGBTRS */
/* Function routines called by ZVSOL: None */
/* ----------------------------------------------------------------------- */
/* This routine manages the solution of the linear system arising from */
/* a chord iteration.  It is called if MITER .ne. 0. */
/* If MITER is 1 or 2, it calls ZGETRS to accomplish this. */
/* If MITER = 3 it updates the coefficient H*RL1 in the diagonal */
/* matrix, and then computes the solution. */
/* If MITER is 4 or 5, it calls ZGBTRS. */
/* Communication with ZVSOL uses the following variables: */
/* WM    = Real work space containing the inverse diagonal matrix if */
/*         MITER = 3 and the LU decomposition of the matrix otherwise. */
/* IWM   = Integer work space containing pivot information, starting at */
/*         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band */
/*         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5. */
/* X     = The right-hand side vector on input, and the solution vector */
/*         on output, of length N. */
/* IERSL = Output flag.  IERSL = 0 if no trouble occurred. */
/*         IERSL = 1 if a singular matrix arose with MITER = 3. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block ZVOD01 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --x;
    --iwm;
    --wm;

    /* Function Body */

    *iersl = 0;
    switch (zvod01_1.miter) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L400;
    }
/*     Replaced LINPACK zgesl with LAPACK zgetrs */
/* 100  CALL ZGESL (WM, N, N, IWM(31), X, 0) */
L100:
    zgetrs_("N", &zvod01_1.n, &c__1, &wm[1], &zvod01_1.n, &iwm[31], &x[1], &
	    zvod01_1.n, &ier, (ftnlen)1);
    return 0;

L300:
    phrl1 = zvod01_1.hrl1;
    zvod01_1.hrl1 = zvod01_1.h__ * zvod01_1.rl1;
    if (zvod01_1.hrl1 == phrl1) {
	goto L330;
    }
    r__ = zvod01_1.hrl1 / phrl1;
    i__1 = zvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__5.r = one, z__5.i = 0.;
	z_div(&z__4, &z__5, &wm[i__]);
	z__3.r = one - z__4.r, z__3.i = -z__4.i;
	z__2.r = r__ * z__3.r, z__2.i = r__ * z__3.i;
	z__1.r = one - z__2.r, z__1.i = -z__2.i;
	di.r = z__1.r, di.i = z__1.i;
	if (z_abs(&di) == zero) {
	    goto L390;
	}
/* L320: */
	i__2 = i__;
	z__2.r = one, z__2.i = 0.;
	z_div(&z__1, &z__2, &di);
	wm[i__2].r = z__1.r, wm[i__2].i = z__1.i;
    }

L330:
    i__2 = zvod01_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L340: */
	i__1 = i__;
	i__3 = i__;
	i__4 = i__;
	z__1.r = wm[i__3].r * x[i__4].r - wm[i__3].i * x[i__4].i, z__1.i = wm[
		i__3].r * x[i__4].i + wm[i__3].i * x[i__4].r;
	x[i__1].r = z__1.r, x[i__1].i = z__1.i;
    }
    return 0;
L390:
    *iersl = 1;
    return 0;

L400:
    ml = iwm[1];
    mu = iwm[2];
    meband = (ml << 1) + mu + 1;
/*     Replaced LINPACK zgbsl with LAPACK zgbtrs */
/*      CALL ZGBSL (WM, MEBAND, N, ML, MU, IWM(31), X, 0) */
    zgbtrs_("N", &zvod01_1.n, &ml, &mu, &c__1, &wm[1], &meband, &iwm[31], &x[
	    1], &zvod01_1.n, &ier, (ftnlen)1);
    return 0;
/* ----------------------- End of Subroutine ZVSOL ----------------------- */
} /* zvsol_ */

/* DECK ZVSRCO */
/* Subroutine */ int zvsrco_(doublereal *rsav, integer *isav, integer *job)
{
    /* Initialized data */

    static integer lenrv1 = 50;
    static integer leniv1 = 33;
    static integer lenrv2 = 1;
    static integer leniv2 = 8;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ----------------------------------------------------------------------- */
/* Call sequence input -- RSAV, ISAV, JOB */
/* Call sequence output -- RSAV, ISAV */
/* COMMON block variables accessed -- All of /ZVOD01/ and /ZVOD02/ */

/* Subroutines/functions called by ZVSRCO: None */
/* ----------------------------------------------------------------------- */
/* This routine saves or restores (depending on JOB) the contents of the */
/* COMMON blocks ZVOD01 and ZVOD02, which are used internally by ZVODE. */

/* RSAV = real array of length 51 or more. */
/* ISAV = integer array of length 41 or more. */
/* JOB  = flag indicating to save or restore the COMMON blocks: */
/*        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV). */
/*        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV). */
/*        A call with JOB = 2 presumes a prior call with JOB = 1. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

    if (*job == 2) {
	goto L100;
    }
    i__1 = lenrv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	rsav[i__] = zvod01_2.rvod1[i__ - 1];
    }
    i__1 = lenrv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	rsav[lenrv1 + i__] = zvod02_2.rvod2[i__ - 1];
    }

    i__1 = leniv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	isav[i__] = zvod01_2.ivod1[i__ - 1];
    }
    i__1 = leniv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	isav[leniv1 + i__] = zvod02_2.ivod2[i__ - 1];
    }

    return 0;

L100:
    i__1 = lenrv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	zvod01_2.rvod1[i__ - 1] = rsav[i__];
    }
    i__1 = lenrv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L115: */
	zvod02_2.rvod2[i__ - 1] = rsav[lenrv1 + i__];
    }

    i__1 = leniv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	zvod01_2.ivod1[i__ - 1] = isav[i__];
    }
    i__1 = leniv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	zvod02_2.ivod2[i__ - 1] = isav[leniv1 + i__];
    }

    return 0;
/* ----------------------- End of Subroutine ZVSRCO ---------------------- */
} /* zvsrco_ */

/* DECK ZEWSET */
/* Subroutine */ int zewset_(integer *n, integer *itol, doublereal *rtol, 
	doublereal *atol, doublecomplex *ycur, doublereal *ewt)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  ZEWSET */
/* ***SUBSIDIARY */
/* ***PURPOSE  Set error weight vector. */
/* ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D, ZEWSET-Z) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This subroutine sets the error weight vector EWT according to */
/*      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N, */
/*  with the subscript on RTOL and/or ATOL possibly replaced by 1 above, */
/*  depending on the value of ITOL. */

/* ***SEE ALSO  DLSODE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   060502  DATE WRITTEN, modified from DEWSET of 930809. */
/* ***END PROLOGUE  ZEWSET */

/* ***FIRST EXECUTABLE STATEMENT  ZEWSET */
    /* Parameter adjustments */
    --ewt;
    --ycur;
    --rtol;
    --atol;

    /* Function Body */
    switch (*itol) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }
L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	ewt[i__] = rtol[1] * z_abs(&ycur[i__]) + atol[1];
    }
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	ewt[i__] = rtol[1] * z_abs(&ycur[i__]) + atol[i__];
    }
    return 0;
L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L35: */
	ewt[i__] = rtol[i__] * z_abs(&ycur[i__]) + atol[1];
    }
    return 0;
L40:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L45: */
	ewt[i__] = rtol[i__] * z_abs(&ycur[i__]) + atol[i__];
    }
    return 0;
/* ----------------------- END OF SUBROUTINE ZEWSET ---------------------- */
} /* zewset_ */

/* DECK ZVNORM */
doublereal zvnorm_(integer *n, doublecomplex *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sum;
    extern doublereal zabssq_(doublecomplex *);

/* ***BEGIN PROLOGUE  ZVNORM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Weighted root-mean-square vector norm. */
/* ***TYPE      DOUBLE COMPLEX (SVNORM-S, DVNORM-D, ZVNORM-Z) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This function routine computes the weighted root-mean-square norm */
/*  of the vector of length N contained in the double complex array V, */
/*  with weights contained in the array W of length N: */
/*    ZVNORM = SQRT( (1/N) * SUM( abs(V(i))**2 * W(i)**2 ) */
/*  The squared absolute value abs(v)**2 is computed by ZABSSQ. */

/* ***SEE ALSO  DLSODE */
/* ***ROUTINES CALLED  ZABSSQ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   060502  DATE WRITTEN, modified from DVNORM of 930809. */
/* ***END PROLOGUE  ZVNORM */

/* ***FIRST EXECUTABLE STATEMENT  ZVNORM */
    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    sum = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing 2nd power */
	d__1 = w[i__];
	sum += zabssq_(&v[i__]) * (d__1 * d__1);
    }
    ret_val = sqrt(sum / *n);
    return ret_val;
/* ----------------------- END OF FUNCTION ZVNORM ------------------------ */
} /* zvnorm_ */

/* DECK ZABSSQ */
doublereal zabssq_(doublecomplex *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

/* ***BEGIN PROLOGUE  ZABSSQ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Squared absolute value of a double complex number. */
/* ***TYPE      DOUBLE PRECISION (ZABSSQ-Z) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This function routine computes the square of the absolute value of */
/*  a double precision complex number Z, */
/*    ZABSSQ = DREAL(Z)**2 * DIMAG(Z)**2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   060502  DATE WRITTEN. */
/* ***END PROLOGUE  ZABSSQ */
/* Computing 2nd power */
    d__1 = z__->r;
/* Computing 2nd power */
    d__2 = d_imag(z__);
    ret_val = d__1 * d__1 + d__2 * d__2;
    return ret_val;
/* ----------------------- END OF FUNCTION ZABSSQ ------------------------ */
} /* zabssq_ */

/* DECK DZSCAL */
/* Subroutine */ int dzscal_(integer *n, doublereal *da, doublecomplex *zx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, ix;

/* ***BEGIN PROLOGUE  DZSCAL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Scale a double complex vector by a double prec. constant. */
/* ***TYPE      DOUBLE PRECISION (DZSCAL-Z) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */
/*  Scales a double complex vector by a double precision constant. */
/*  Minor modification of BLAS routine ZSCAL. */
/* ***REVISION HISTORY  (YYMMDD) */
/*   060530  DATE WRITTEN. */
/* ***END PROLOGUE  DZSCAL */

    /* Parameter adjustments */
    --zx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }
/* Code for increment not equal to 1 */
    ix = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	i__3 = ix;
	z__1.r = *da * zx[i__3].r, z__1.i = *da * zx[i__3].i;
	zx[i__2].r = z__1.r, zx[i__2].i = z__1.i;
	ix += *incx;
/* L10: */
    }
    return 0;
/* Code for increment equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	z__1.r = *da * zx[i__3].r, z__1.i = *da * zx[i__3].i;
	zx[i__2].r = z__1.r, zx[i__2].i = z__1.i;
/* L30: */
    }
    return 0;
} /* dzscal_ */

/* DECK DZAXPY */
/* Subroutine */ int dzaxpy_(integer *n, doublereal *da, doublecomplex *zx, 
	integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, ix, iy;

/* ***BEGIN PROLOGUE  DZAXPY */
/* ***PURPOSE  Real constant times a complex vector plus a complex vector. */
/* ***TYPE      DOUBLE PRECISION (DZAXPY-Z) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */
/*  Add a D.P. real constant times a complex vector to a complex vector. */
/*  Minor modification of BLAS routine ZAXPY. */
/* ***REVISION HISTORY  (YYMMDD) */
/*   060530  DATE WRITTEN. */
/* ***END PROLOGUE  DZAXPY */
    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (abs(*da) == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/* Code for unequal increments or equal increments not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iy;
	i__3 = iy;
	i__4 = ix;
	z__2.r = *da * zx[i__4].r, z__2.i = *da * zx[i__4].i;
	z__1.r = zy[i__3].r + z__2.r, z__1.i = zy[i__3].i + z__2.i;
	zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	i__4 = i__;
	z__2.r = *da * zx[i__4].r, z__2.i = *da * zx[i__4].i;
	z__1.r = zy[i__3].r + z__2.r, z__1.i = zy[i__3].i + z__2.i;
	zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
/* L30: */
    }
    return 0;
} /* dzaxpy_ */

/* DECK DUMACH */
doublereal dumach_(void)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal u, comp;
    extern /* Subroutine */ int dumsum_(doublereal *, doublereal *, 
	    doublereal *);

/* ***BEGIN PROLOGUE  DUMACH */
/* ***PURPOSE  Compute the unit roundoff of the machine. */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */
/* *Usage: */
/*        DOUBLE PRECISION  A, DUMACH */
/*        A = DUMACH() */

/* *Function Return Values: */
/*     A : the unit roundoff of the machine. */

/* *Description: */
/*     The unit roundoff is defined as the smallest positive machine */
/*     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH */
/*     in a machine-independent manner. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DUMSUM */
/* ***REVISION HISTORY  (YYYYMMDD) */
/*   19930216  DATE WRITTEN */
/*   19930818  Added SLATEC-format prologue.  (FNF) */
/*   20030707  Added DUMSUM to force normal storage of COMP.  (ACH) */
/* ***END PROLOGUE  DUMACH */

/* ***FIRST EXECUTABLE STATEMENT  DUMACH */
    u = 1.;
L10:
    u *= .5;
    dumsum_(&c_b642, &u, &comp);
    if (comp != 1.) {
	goto L10;
    }
    ret_val = u * 2.;
    return ret_val;
/* ----------------------- End of Function DUMACH ------------------------ */
} /* dumach_ */

/* Subroutine */ int dumsum_(doublereal *a, doublereal *b, doublereal *c__)
{
/*     Routine to force normal storing of A + B, for DUMACH. */
    *c__ = *a + *b;
    return 0;
} /* dumsum_ */

