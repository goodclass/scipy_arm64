/* lsoda.f -- translated by f2c (version 20190311).
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

struct {
    doublereal rowns[209], ccmax, el0, h__, hmin, hmxi, hu, rc, tn, uround;
    integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, 
	    nhnil, ntrep, nslast, nyh, iowns[6], icf, ierpj, iersl, jcur, 
	    jstart, kflag, l, meth, miter, maxord, maxcor, msbp, mxncf, n, nq,
	     nst, nfe, nje, nqu;
} ls0001_;

#define ls0001_1 ls0001_

struct {
    doublereal tsw, rowns2[20], pdnorm;
    integer insufr, insufi, ixpr, iowns2[2], jtyp, mused, mxordn, mxords;
} lsa001_;

#define lsa001_1 lsa001_

/* Table of constant values */

static integer c__60 = 60;
static integer c__103 = 103;
static integer c__0 = 0;
static doublereal c_b42 = 0.;
static integer c__50 = 50;
static integer c__2 = 2;
static integer c__104 = 104;
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__101 = 101;
static integer c__102 = 102;
static integer c__105 = 105;
static integer c__106 = 106;
static integer c__107 = 107;
static integer c__301 = 301;
static integer c__202 = 202;
static integer c__203 = 203;
static integer c__204 = 204;
static integer c__205 = 205;
static integer c__30 = 30;
static integer c__206 = 206;
static integer c__207 = 207;
static integer c__3 = 3;
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
static integer c__28 = 28;
static integer c__29 = 29;
static integer c__302 = 302;
static integer c__303 = 303;

/* Subroutine */ int lsoda_(S_fp f, integer *neq, doublereal *y, doublereal *
	t, doublereal *tout, integer *itol, doublereal *rtol, doublereal *
	atol, integer *itask, integer *istate, integer *iopt, doublereal *
	rwork, integer *lrw, integer *iwork, integer *liw, U_fp jac, integer *
	jt)
{
    /* Initialized data */

    static integer mord[2] = { 12,5 };
    static integer mxstp0 = 500;
    static integer mxhnl0 = 10;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal h0;
    static integer i1, i2;
    static doublereal w0;
    static integer ml;
    static doublereal rh;
    static integer mu;
    static doublereal tp;
    static integer lf0;
    static doublereal big;
    static integer kgo;
    static doublereal ayi, hmx, tol, sum;
    static integer len1, len2;
    extern /* Subroutine */ int prja_();
    static doublereal hmax;
    static logical ihit;
    static integer isav[50];
    static doublereal ewti, size, rsav[240];
    static integer len1c, len1n, len1s, iflag;
    extern /* Subroutine */ int srcma_(doublereal *, integer *, integer *);
    static doublereal atoli;
    static integer leniw, lenwm;
    extern /* Subroutine */ int stoda_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, S_fp, U_fp, U_fp, U_fp);
    static integer imxer;
    static doublereal tcrit;
    static integer lenrw;
    static doublereal rtoli, tdist, tolsf;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int ewset_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), intdy_(doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal tnext;
    extern /* Subroutine */ int solsy_();
    static integer leniwc, lenrwc, lenrwn, lenrws;
    extern doublereal vmnorm_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int xerrwv_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, ftnlen);

/* ----------------------------------------------------------------------- */
/* this is the 24 feb 1997 version of */
/* lsoda.. livermore solver for ordinary differential equations, with */
/*         automatic method switching for stiff and nonstiff problems. */

/* this version is in double precision. */

/* lsoda solves the initial value problem for stiff or nonstiff */
/* systems of first order ode-s, */
/*     dy/dt = f(t,y) ,  or, in component form, */
/*     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq). */

/* this a variant version of the lsode package. */
/* it switches automatically between stiff and nonstiff methods. */
/* this means that the user does not have to determine whether the */
/* problem is stiff or not, and the solver will automatically choose the */
/* appropriate method.  it always starts with the nonstiff method. */

/* authors.. */
/*                linda r. petzold  and  alan c. hindmarsh, */
/*                computing and mathematics research division, l-316 */
/*                lawrence livermore national laboratory */
/*                livermore, ca 94550. */

/* references.. */
/* 1.  alan c. hindmarsh,  odepack, a systematized collection of ode */
/*     solvers, in scientific computing, r. s. stepleman et al. (eds.), */
/*     north-holland, amsterdam, 1983, pp. 55-64. */
/* 2.  linda r. petzold, automatic selection of methods for solving */
/*     stiff and nonstiff systems of ordinary differential equations, */
/*     siam j. sci. stat. comput. 4 (1983), pp. 136-148. */
/* ----------------------------------------------------------------------- */
/* summary of usage. */

/* communication between the user and the lsoda package, for normal */
/* situations, is summarized here.  this summary describes only a subset */
/* of the full set of options available.  see the full description for */
/* details, including alternative treatment of the jacobian matrix, */
/* optional inputs and outputs, nonstandard options, and */
/* instructions for special situations.  see also the example */
/* problem (with program and output) following this summary. */

/* a. first provide a subroutine of the form.. */
/*               subroutine f (neq, t, y, ydot) */
/*               dimension y(neq), ydot(neq) */
/* which supplies the vector function f by loading ydot(i) with f(i). */

/* b. write a main program which calls subroutine lsoda once for */
/* each point at which answers are desired.  this should also provide */
/* for possible use of logical unit 6 for output of error messages */
/* by lsoda.  on the first call to lsoda, supply arguments as follows.. */
/* f      = name of subroutine for right-hand side vector f. */
/*          this name must be declared external in calling program. */
/* neq    = number of first order ode-s. */
/* y      = array of initial values, of length neq. */
/* t      = the initial value of the independent variable. */
/* tout   = first point where output is desired (.ne. t). */
/* itol   = 1 or 2 according as atol (below) is a scalar or array. */
/* rtol   = relative tolerance parameter (scalar). */
/* atol   = absolute tolerance parameter (scalar or array). */
/*          the estimated local error in y(i) will be controlled so as */
/*          to be less than */
/*             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or */
/*             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2. */
/*          thus the local error test passes if, in each component, */
/*          either the absolute error is less than atol (or atol(i)), */
/*          or the relative error is less than rtol. */
/*          use rtol = 0.0 for pure absolute error control, and */
/*          use atol = 0.0 (or atol(i) = 0.0) for pure relative error */
/*          control.  caution.. actual (global) errors may exceed these */
/*          local tolerances, so choose them conservatively. */
/* itask  = 1 for normal computation of output values of y at t = tout. */
/* istate = integer flag (input and output).  set istate = 1. */
/* iopt   = 0 to indicate no optional inputs used. */
/* rwork  = real work array of length at least.. */
/*             22 + neq * max(16, neq + 9). */
/*          see also paragraph e below. */
/* lrw    = declared length of rwork (in user-s dimension). */
/* iwork  = integer work array of length at least  20 + neq. */
/* liw    = declared length of iwork (in user-s dimension). */
/* jac    = name of subroutine for jacobian matrix. */
/*          use a dummy name.  see also paragraph e below. */
/* jt     = jacobian type indicator.  set jt = 2. */
/*          see also paragraph e below. */
/* note that the main program must declare arrays y, rwork, iwork, */
/* and possibly atol. */

/* c. the output from the first call (or any call) is.. */
/*      y = array of computed values of y(t) vector. */
/*      t = corresponding value of independent variable (normally tout). */
/* istate = 2  if lsoda was successful, negative otherwise. */
/*          -1 means excess work done on this call (perhaps wrong jt). */
/*          -2 means excess accuracy requested (tolerances too small). */
/*          -3 means illegal input detected (see printed message). */
/*          -4 means repeated error test failures (check all inputs). */
/*          -5 means repeated convergence failures (perhaps bad jacobian */
/*             supplied or wrong choice of jt or tolerances). */
/*          -6 means error weight became zero during problem. (solution */
/*             component i vanished, and atol or atol(i) = 0.) */
/*          -7 means work space insufficient to finish (see messages). */

/* d. to continue the integration after a successful return, simply */
/* reset tout and call lsoda again.  no other parameters need be reset. */

/* e. note.. if and when lsoda regards the problem as stiff, and */
/* switches methods accordingly, it must make use of the neq by neq */
/* jacobian matrix, j = df/dy.  for the sake of simplicity, the */
/* inputs to lsoda recommended in paragraph b above cause lsoda to */
/* treat j as a full matrix, and to approximate it internally by */
/* difference quotients.  alternatively, j can be treated as a band */
/* matrix (with great potential reduction in the size of the rwork */
/* array).  also, in either the full or banded case, the user can supply */
/* j in closed form, with a routine whose name is passed as the jac */
/* argument.  these alternatives are described in the paragraphs on */
/* rwork, jac, and jt in the full description of the call sequence below. */

/* ----------------------------------------------------------------------- */
/* example problem. */

/* the following is a simple example problem, with the coding */
/* needed for its solution by lsoda.  the problem is from chemical */
/* kinetics, and consists of the following three rate equations.. */
/*     dy1/dt = -.04*y1 + 1.e4*y2*y3 */
/*     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2 */
/*     dy3/dt = 3.e7*y2**2 */
/* on the interval from t = 0.0 to t = 4.e10, with initial conditions */
/* y1 = 1.0, y2 = y3 = 0.  the problem is stiff. */

/* the following coding solves this problem with lsoda, */
/* printing results at t = .4, 4., ..., 4.e10.  it uses */
/* itol = 2 and atol much smaller for y2 than y1 or y3 because */
/* y2 has much smaller values. */
/* at the end of the run, statistical quantities of interest are */
/* printed (see optional outputs in the full description below). */

/*     external fex */
/*     double precision atol, rtol, rwork, t, tout, y */
/*     dimension y(3), atol(3), rwork(70), iwork(23) */
/*     neq = 3 */
/*     y(1) = 1.0d0 */
/*     y(2) = 0.0d0 */
/*     y(3) = 0.0d0 */
/*     t = 0.0d0 */
/*     tout = 0.4d0 */
/*     itol = 2 */
/*     rtol = 1.0d-4 */
/*     atol(1) = 1.0d-6 */
/*     atol(2) = 1.0d-10 */
/*     atol(3) = 1.0d-6 */
/*     itask = 1 */
/*     istate = 1 */
/*     iopt = 0 */
/*     lrw = 70 */
/*     liw = 23 */
/*     jt = 2 */
/*     do 40 iout = 1,12 */
/*       call lsoda(fex,neq,y,t,tout,itol,rtol,atol,itask,istate, */
/*    1     iopt,rwork,lrw,iwork,liw,jdum,jt) */
/*       write(6,20)t,y(1),y(2),y(3) */
/* 20    format(' at t =',e12.4,'   y =',3e14.6) */
/*       if (istate .lt. 0) go to 80 */
/* 40    tout = tout*10.0d0 */
/*     write(6,60)iwork(11),iwork(12),iwork(13),iwork(19),rwork(15) */
/* 60  format(/' no. steps =',i4,'  no. f-s =',i4,'  no. j-s =',i4/ */
/*    1   ' method last used =',i2,'   last switch was at t =',e12.4) */
/*     stop */
/* 80  write(6,90)istate */
/* 90  format(///' error halt.. istate =',i3) */
/*     stop */
/*     end */

/*     subroutine fex (neq, t, y, ydot) */
/*     double precision t, y, ydot */
/*     dimension y(3), ydot(3) */
/*     ydot(1) = -.04d0*y(1) + 1.0d4*y(2)*y(3) */
/*     ydot(3) = 3.0d7*y(2)*y(2) */
/*     ydot(2) = -ydot(1) - ydot(3) */
/*     return */
/*     end */

/* the output of this program (on a cdc-7600 in single precision) */
/* is as follows.. */

/*   at t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02 */
/*   at t =  4.0000e+00   y =  9.055333e-01  2.240655e-05  9.444430e-02 */
/*   at t =  4.0000e+01   y =  7.158403e-01  9.186334e-06  2.841505e-01 */
/*   at t =  4.0000e+02   y =  4.505250e-01  3.222964e-06  5.494717e-01 */
/*   at t =  4.0000e+03   y =  1.831975e-01  8.941774e-07  8.168016e-01 */
/*   at t =  4.0000e+04   y =  3.898730e-02  1.621940e-07  9.610125e-01 */
/*   at t =  4.0000e+05   y =  4.936363e-03  1.984221e-08  9.950636e-01 */
/*   at t =  4.0000e+06   y =  5.161831e-04  2.065786e-09  9.994838e-01 */
/*   at t =  4.0000e+07   y =  5.179817e-05  2.072032e-10  9.999482e-01 */
/*   at t =  4.0000e+08   y =  5.283401e-06  2.113371e-11  9.999947e-01 */
/*   at t =  4.0000e+09   y =  4.659031e-07  1.863613e-12  9.999995e-01 */
/*   at t =  4.0000e+10   y =  1.404280e-08  5.617126e-14  1.000000e+00 */

/*   no. steps = 361  no. f-s = 693  no. j-s =  64 */
/*   method last used = 2   last switch was at t =  6.0092e-03 */
/* ----------------------------------------------------------------------- */
/* full description of user interface to lsoda. */

/* the user interface to lsoda consists of the following parts. */

/* i.   the call sequence to subroutine lsoda, which is a driver */
/*      routine for the solver.  this includes descriptions of both */
/*      the call sequence arguments and of user-supplied routines. */
/*      following these descriptions is a description of */
/*      optional inputs available through the call sequence, and then */
/*      a description of optional outputs (in the work arrays). */

/* ii.  descriptions of other routines in the lsoda package that may be */
/*      (optionally) called by the user.  these provide the ability to */
/*      alter error message handling, save and restore the internal */
/*      common, and obtain specified derivatives of the solution y(t). */

/* iii. descriptions of common blocks to be declared in overlay */
/*      or similar environments, or to be saved when doing an interrupt */
/*      of the problem and continued solution later. */

/* iv.  description of a subroutine in the lsoda package, */
/*      which the user may replace with his own version, if desired. */
/*      this relates to the measurement of errors. */

/* ----------------------------------------------------------------------- */
/* part i.  call sequence. */

/* the call sequence parameters used for input only are */
/*     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, jt, */
/* and those used for both input and output are */
/*     y, t, istate. */
/* the work arrays rwork and iwork are also used for conditional and */
/* optional inputs and optional outputs.  (the term output here refers */
/* to the return from subroutine lsoda to the user-s calling program.) */

/* the legality of input parameters will be thoroughly checked on the */
/* initial call for the problem, but not checked thereafter unless a */
/* change in input parameters is flagged by istate = 3 on input. */

/* the descriptions of the call arguments are as follows. */

/* f      = the name of the user-supplied subroutine defining the */
/*          ode system.  the system must be put in the first-order */
/*          form dy/dt = f(t,y), where f is a vector-valued function */
/*          of the scalar t and the vector y.  subroutine f is to */
/*          compute the function f.  it is to have the form */
/*               subroutine f (neq, t, y, ydot) */
/*               dimension y(1), ydot(1) */
/*          where neq, t, and y are input, and the array ydot = f(t,y) */
/*          is output.  y and ydot are arrays of length neq. */
/*          (in the dimension statement above, 1 is a dummy */
/*          dimension.. it can be replaced by any value.) */
/*          subroutine f should not alter y(1),...,y(neq). */
/*          f must be declared external in the calling program. */

/*          subroutine f may access user-defined quantities in */
/*          neq(2),... and/or in y(neq(1)+1),... if neq is an array */
/*          (dimensioned in f) and/or y has length exceeding neq(1). */
/*          see the descriptions of neq and y below. */

/*          if quantities computed in the f routine are needed */
/*          externally to lsoda, an extra call to f should be made */
/*          for this purpose, for consistent and accurate results. */
/*          if only the derivative dy/dt is needed, use intdy instead. */

/* neq    = the size of the ode system (number of first order */
/*          ordinary differential equations).  used only for input. */
/*          neq may be decreased, but not increased, during the problem. */
/*          if neq is decreased (with istate = 3 on input), the */
/*          remaining components of y should be left undisturbed, if */
/*          these are to be accessed in f and/or jac. */

/*          normally, neq is a scalar, and it is generally referred to */
/*          as a scalar in this user interface description.  however, */
/*          neq may be an array, with neq(1) set to the system size. */
/*          (the lsoda package accesses only neq(1).)  in either case, */
/*          this parameter is passed as the neq argument in all calls */
/*          to f and jac.  hence, if it is an array, locations */
/*          neq(2),... may be used to store other integer data and pass */
/*          it to f and/or jac.  subroutines f and/or jac must include */
/*          neq in a dimension statement in that case. */

/* y      = a real array for the vector of dependent variables, of */
/*          length neq or more.  used for both input and output on the */
/*          first call (istate = 1), and only for output on other calls. */
/*          on the first call, y must contain the vector of initial */
/*          values.  on output, y contains the computed solution vector, */
/*          evaluated at t.  if desired, the y array may be used */
/*          for other purposes between calls to the solver. */

/*          this array is passed as the y argument in all calls to */
/*          f and jac.  hence its length may exceed neq, and locations */
/*          y(neq+1),... may be used to store other real data and */
/*          pass it to f and/or jac.  (the lsoda package accesses only */
/*          y(1),...,y(neq).) */

/* t      = the independent variable.  on input, t is used only on the */
/*          first call, as the initial point of the integration. */
/*          on output, after each call, t is the value at which a */
/*          computed solution y is evaluated (usually the same as tout). */
/*          on an error return, t is the farthest point reached. */

/* tout   = the next value of t at which a computed solution is desired. */
/*          used only for input. */

/*          when starting the problem (istate = 1), tout may be equal */
/*          to t for one call, then should .ne. t for the next call. */
/*          for the initial t, an input value of tout .ne. t is used */
/*          in order to determine the direction of the integration */
/*          (i.e. the algebraic sign of the step sizes) and the rough */
/*          scale of the problem.  integration in either direction */
/*          (forward or backward in t) is permitted. */

/*          if itask = 2 or 5 (one-step modes), tout is ignored after */
/*          the first call (i.e. the first call with tout .ne. t). */
/*          otherwise, tout is required on every call. */

/*          if itask = 1, 3, or 4, the values of tout need not be */
/*          monotone, but a value of tout which backs up is limited */
/*          to the current internal t interval, whose endpoints are */
/*          tcur - hu and tcur (see optional outputs, below, for */
/*          tcur and hu). */

/* itol   = an indicator for the type of error control.  see */
/*          description below under atol.  used only for input. */

/* rtol   = a relative error tolerance parameter, either a scalar or */
/*          an array of length neq.  see description below under atol. */
/*          input only. */

/* atol   = an absolute error tolerance parameter, either a scalar or */
/*          an array of length neq.  input only. */

/*             the input parameters itol, rtol, and atol determine */
/*          the error control performed by the solver.  the solver will */
/*          control the vector e = (e(i)) of estimated local errors */
/*          in y, according to an inequality of the form */
/*                      max-norm of ( e(i)/ewt(i) )   .le.   1, */
/*          where ewt = (ewt(i)) is a vector of positive error weights. */
/*          the values of rtol and atol should all be non-negative. */
/*          the following table gives the types (scalar/array) of */
/*          rtol and atol, and the corresponding form of ewt(i). */

/*             itol    rtol       atol          ewt(i) */
/*              1     scalar     scalar     rtol*abs(y(i)) + atol */
/*              2     scalar     array      rtol*abs(y(i)) + atol(i) */
/*              3     array      scalar     rtol(i)*abs(y(i)) + atol */
/*              4     array      array      rtol(i)*abs(y(i)) + atol(i) */

/*          when either of these parameters is a scalar, it need not */
/*          be dimensioned in the user-s calling program. */

/*          if none of the above choices (with itol, rtol, and atol */
/*          fixed throughout the problem) is suitable, more general */
/*          error controls can be obtained by substituting a */
/*          user-supplied routine for the setting of ewt. */
/*          see part iv below. */

/*          if global errors are to be estimated by making a repeated */
/*          run on the same problem with smaller tolerances, then all */
/*          components of rtol and atol (i.e. of ewt) should be scaled */
/*          down uniformly. */

/* itask  = an index specifying the task to be performed. */
/*          input only.  itask has the following values and meanings. */
/*          1  means normal computation of output values of y(t) at */
/*             t = tout (by overshooting and interpolating). */
/*          2  means take one step only and return. */
/*          3  means stop at the first internal mesh point at or */
/*             beyond t = tout and return. */
/*          4  means normal computation of output values of y(t) at */
/*             t = tout but without overshooting t = tcrit. */
/*             tcrit must be input as rwork(1).  tcrit may be equal to */
/*             or beyond tout, but not behind it in the direction of */
/*             integration.  this option is useful if the problem */
/*             has a singularity at or beyond t = tcrit. */
/*          5  means take one step, without passing tcrit, and return. */
/*             tcrit must be input as rwork(1). */

/*          note..  if itask = 4 or 5 and the solver reaches tcrit */
/*          (within roundoff), it will return t = tcrit (exactly) to */
/*          indicate this (unless itask = 4 and tout comes before tcrit, */
/*          in which case answers at t = tout are returned first). */

/* istate = an index used for input and output to specify the */
/*          the state of the calculation. */

/*          on input, the values of istate are as follows. */
/*          1  means this is the first call for the problem */
/*             (initializations will be done).  see note below. */
/*          2  means this is not the first call, and the calculation */
/*             is to continue normally, with no change in any input */
/*             parameters except possibly tout and itask. */
/*             (if itol, rtol, and/or atol are changed between calls */
/*             with istate = 2, the new values will be used but not */
/*             tested for legality.) */
/*          3  means this is not the first call, and the */
/*             calculation is to continue normally, but with */
/*             a change in input parameters other than */
/*             tout and itask.  changes are allowed in */
/*             neq, itol, rtol, atol, iopt, lrw, liw, jt, ml, mu, */
/*             and any optional inputs except h0, mxordn, and mxords. */
/*             (see iwork description for ml and mu.) */
/*          note..  a preliminary call with tout = t is not counted */
/*          as a first call here, as no initialization or checking of */
/*          input is done.  (such a call is sometimes useful for the */
/*          purpose of outputting the initial conditions.) */
/*          thus the first call for which tout .ne. t requires */
/*          istate = 1 on input. */

/*          on output, istate has the following values and meanings. */
/*           1  means nothing was done, as tout was equal to t with */
/*              istate = 1 on input.  (however, an internal counter was */
/*              set to detect and prevent repeated calls of this type.) */
/*           2  means the integration was performed successfully. */
/*          -1  means an excessive amount of work (more than mxstep */
/*              steps) was done on this call, before completing the */
/*              requested task, but the integration was otherwise */
/*              successful as far as t.  (mxstep is an optional input */
/*              and is normally 500.)  to continue, the user may */
/*              simply reset istate to a value .gt. 1 and call again */
/*              (the excess work step counter will be reset to 0). */
/*              in addition, the user may increase mxstep to avoid */
/*              this error return (see below on optional inputs). */
/*          -2  means too much accuracy was requested for the precision */
/*              of the machine being used.  this was detected before */
/*              completing the requested task, but the integration */
/*              was successful as far as t.  to continue, the tolerance */
/*              parameters must be reset, and istate must be set */
/*              to 3.  the optional output tolsf may be used for this */
/*              purpose.  (note.. if this condition is detected before */
/*              taking any steps, then an illegal input return */
/*              (istate = -3) occurs instead.) */
/*          -3  means illegal input was detected, before taking any */
/*              integration steps.  see written message for details. */
/*              note..  if the solver detects an infinite loop of calls */
/*              to the solver with illegal input, it will cause */
/*              the run to stop. */
/*          -4  means there were repeated error test failures on */
/*              one attempted step, before completing the requested */
/*              task, but the integration was successful as far as t. */
/*              the problem may have a singularity, or the input */
/*              may be inappropriate. */
/*          -5  means there were repeated convergence test failures on */
/*              one attempted step, before completing the requested */
/*              task, but the integration was successful as far as t. */
/*              this may be caused by an inaccurate jacobian matrix, */
/*              if one is being used. */
/*          -6  means ewt(i) became zero for some i during the */
/*              integration.  pure relative error control (atol(i)=0.0) */
/*              was requested on a variable which has now vanished. */
/*              the integration was successful as far as t. */
/*          -7  means the length of rwork and/or iwork was too small to */
/*              proceed, but the integration was successful as far as t. */
/*              this happens when lsoda chooses to switch methods */
/*              but lrw and/or liw is too small for the new method. */

/*          note..  since the normal output value of istate is 2, */
/*          it does not need to be reset for normal continuation. */
/*          also, since a negative input value of istate will be */
/*          regarded as illegal, a negative output value requires the */
/*          user to change it, and possibly other inputs, before */
/*          calling the solver again. */

/* iopt   = an integer flag to specify whether or not any optional */
/*          inputs are being used on this call.  input only. */
/*          the optional inputs are listed separately below. */
/*          iopt = 0 means no optional inputs are being used. */
/*                   default values will be used in all cases. */
/*          iopt = 1 means one or more optional inputs are being used. */

/* rwork  = a real array (double precision) for work space, and (in the */
/*          first 20 words) for conditional and optional inputs and */
/*          optional outputs. */
/*          as lsoda switches automatically between stiff and nonstiff */
/*          methods, the required length of rwork can change during the */
/*          problem.  thus the rwork array passed to lsoda can either */
/*          have a static (fixed) length large enough for both methods, */
/*          or have a dynamic (changing) length altered by the calling */
/*          program in response to output from lsoda. */

/*                       --- fixed length case --- */
/*          if the rwork length is to be fixed, it should be at least */
/*               max (lrn, lrs), */
/*          where lrn and lrs are the rwork lengths required when the */
/*          current method is nonstiff or stiff, respectively. */

/*          the separate rwork length requirements lrn and lrs are */
/*          as follows.. */
/*          if neq is constant and the maximum method orders have */
/*          their default values, then */
/*             lrn = 20 + 16*neq, */
/*             lrs = 22 + 9*neq + neq**2           if jt = 1 or 2, */
/*             lrs = 22 + 10*neq + (2*ml+mu)*neq   if jt = 4 or 5. */
/*          under any other conditions, lrn and lrs are given by.. */
/*             lrn = 20 + nyh*(mxordn+1) + 3*neq, */
/*             lrs = 20 + nyh*(mxords+1) + 3*neq + lmat, */
/*          where */
/*             nyh    = the initial value of neq, */
/*             mxordn = 12, unless a smaller  value is given as an */
/*                      optional input, */
/*             mxords = 5, unless a smaller value is given as an */
/*                      optional input, */
/*             lmat   = length of matrix work space.. */
/*             lmat   = neq**2 + 2              if jt = 1 or 2, */
/*             lmat   = (2*ml + mu + 1)*neq + 2 if jt = 4 or 5. */

/*                       --- dynamic length case --- */
/*          if the length of rwork is to be dynamic, then it should */
/*          be at least lrn or lrs, as defined above, depending on the */
/*          current method.  initially, it must be at least lrn (since */
/*          lsoda starts with the nonstiff method).  on any return */
/*          from lsoda, the optional output mcur indicates the current */
/*          method.  if mcur differs from the value it had on the */
/*          previous return, or if there has only been one call to */
/*          lsoda and mcur is now 2, then lsoda has switched */
/*          methods during the last call, and the length of rwork */
/*          should be reset (to lrn if mcur = 1, or to lrs if */
/*          mcur = 2).  (an increase in the rwork length is required */
/*          if lsoda returned istate = -7, but not otherwise.) */
/*          after resetting the length, call lsoda with istate = 3 */
/*          to signal that change. */

/* lrw    = the length of the array rwork, as declared by the user. */
/*          (this will be checked by the solver.) */

/* iwork  = an integer array for work space. */
/*          as lsoda switches automatically between stiff and nonstiff */
/*          methods, the required length of iwork can change during */
/*          problem, between */
/*             lis = 20 + neq   and   lin = 20, */
/*          respectively.  thus the iwork array passed to lsoda can */
/*          either have a fixed length of at least 20 + neq, or have a */
/*          dynamic length of at least lin or lis, depending on the */
/*          current method.  the comments on dynamic length under */
/*          rwork above apply here.  initially, this length need */
/*          only be at least lin = 20. */

/*          the first few words of iwork are used for conditional and */
/*          optional inputs and optional outputs. */

/*          the following 2 words in iwork are conditional inputs.. */
/*            iwork(1) = ml     these are the lower and upper */
/*            iwork(2) = mu     half-bandwidths, respectively, of the */
/*                       banded jacobian, excluding the main diagonal. */
/*                       the band is defined by the matrix locations */
/*                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu */
/*                       must satisfy  0 .le.  ml,mu  .le. neq-1. */
/*                       these are required if jt is 4 or 5, and */
/*                       ignored otherwise.  ml and mu may in fact be */
/*                       the band parameters for a matrix to which */
/*                       df/dy is only approximately equal. */

/* liw    = the length of the array iwork, as declared by the user. */
/*          (this will be checked by the solver.) */

/* note.. the base addresses of the work arrays must not be */
/* altered between calls to lsoda for the same problem. */
/* the contents of the work arrays must not be altered */
/* between calls, except possibly for the conditional and */
/* optional inputs, and except for the last 3*neq words of rwork. */
/* the latter space is used for internal scratch space, and so is */
/* available for use by the user outside lsoda between calls, if */
/* desired (but not for use by f or jac). */

/* jac    = the name of the user-supplied routine to compute the */
/*          jacobian matrix, df/dy, if jt = 1 or 4.  the jac routine */
/*          is optional, but if the problem is expected to be stiff much */
/*          of the time, you are encouraged to supply jac, for the sake */
/*          of efficiency.  (alternatively, set jt = 2 or 5 to have */
/*          lsoda compute df/dy internally by difference quotients.) */
/*          if and when lsoda uses df/dy, if treats this neq by neq */
/*          matrix either as full (jt = 1 or 2), or as banded (jt = */
/*          4 or 5) with half-bandwidths ml and mu (discussed under */
/*          iwork above).  in either case, if jt = 1 or 4, the jac */
/*          routine must compute df/dy as a function of the scalar t */
/*          and the vector y.  it is to have the form */
/*               subroutine jac (neq, t, y, ml, mu, pd, nrowpd) */
/*               dimension y(1), pd(nrowpd,1) */
/*          where neq, t, y, ml, mu, and nrowpd are input and the array */
/*          pd is to be loaded with partial derivatives (elements of */
/*          the jacobian matrix) on output.  pd must be given a first */
/*          dimension of nrowpd.  t and y have the same meaning as in */
/*          subroutine f.  (in the dimension statement above, 1 is a */
/*          dummy dimension.. it can be replaced by any value.) */
/*               in the full matrix case (jt = 1), ml and mu are */
/*          ignored, and the jacobian is to be loaded into pd in */
/*          columnwise manner, with df(i)/dy(j) loaded into pd(i,j). */
/*               in the band matrix case (jt = 4), the elements */
/*          within the band are to be loaded into pd in columnwise */
/*          manner, with diagonal lines of df/dy loaded into the rows */
/*          of pd.  thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j). */
/*          ml and mu are the half-bandwidth parameters (see iwork). */
/*          the locations in pd in the two triangular areas which */
/*          correspond to nonexistent matrix elements can be ignored */
/*          or loaded arbitrarily, as they are overwritten by lsoda. */
/*               jac need not provide df/dy exactly.  a crude */
/*          approximation (possibly with a smaller bandwidth) will do. */
/*               in either case, pd is preset to zero by the solver, */
/*          so that only the nonzero elements need be loaded by jac. */
/*          each call to jac is preceded by a call to f with the same */
/*          arguments neq, t, and y.  thus to gain some efficiency, */
/*          intermediate quantities shared by both calculations may be */
/*          saved in a user common block by f and not recomputed by jac, */
/*          if desired.  also, jac may alter the y array, if desired. */
/*          jac must be declared external in the calling program. */
/*               subroutine jac may access user-defined quantities in */
/*          neq(2),... and/or in y(neq(1)+1),... if neq is an array */
/*          (dimensioned in jac) and/or y has length exceeding neq(1). */
/*          see the descriptions of neq and y above. */

/* jt     = jacobian type indicator.  used only for input. */
/*          jt specifies how the jacobian matrix df/dy will be */
/*          treated, if and when lsoda requires this matrix. */
/*          jt has the following values and meanings.. */
/*           1 means a user-supplied full (neq by neq) jacobian. */
/*           2 means an internally generated (difference quotient) full */
/*             jacobian (using neq extra calls to f per df/dy value). */
/*           4 means a user-supplied banded jacobian. */
/*           5 means an internally generated banded jacobian (using */
/*             ml+mu+1 extra calls to f per df/dy evaluation). */
/*          if jt = 1 or 4, the user must supply a subroutine jac */
/*          (the name is arbitrary) as described above under jac. */
/*          if jt = 2 or 5, a dummy argument can be used. */
/* ----------------------------------------------------------------------- */
/* optional inputs. */

/* the following is a list of the optional inputs provided for in the */
/* call sequence.  (see also part ii.)  for each such input variable, */
/* this table lists its name as used in this documentation, its */
/* location in the call sequence, its meaning, and the default value. */
/* the use of any of these inputs requires iopt = 1, and in that */
/* case all of these inputs are examined.  a value of zero for any */
/* of these optional inputs will cause the default value to be used. */
/* thus to use a subset of the optional inputs, simply preload */
/* locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and */
/* then set those of interest to nonzero values. */

/* name    location      meaning and default value */

/* h0      rwork(5)  the step size to be attempted on the first step. */
/*                   the default value is determined by the solver. */

/* hmax    rwork(6)  the maximum absolute step size allowed. */
/*                   the default value is infinite. */

/* hmin    rwork(7)  the minimum absolute step size allowed. */
/*                   the default value is 0.  (this lower bound is not */
/*                   enforced on the final step before reaching tcrit */
/*                   when itask = 4 or 5.) */

/* ixpr    iwork(5)  flag to generate extra printing at method switches. */
/*                   ixpr = 0 means no extra printing (the default). */
/*                   ixpr = 1 means print data on each switch. */
/*                   t, h, and nst will be printed on the same logical */
/*                   unit as used for error messages. */

/* mxstep  iwork(6)  maximum number of (internally defined) steps */
/*                   allowed during one call to the solver. */
/*                   the default value is 500. */

/* mxhnil  iwork(7)  maximum number of messages printed (per problem) */
/*                   warning that t + h = t on a step (h = step size). */
/*                   this must be positive to result in a non-default */
/*                   value.  the default value is 10. */

/* mxordn  iwork(8)  the maximum order to be allowed for the nonstiff */
/*                   (adams) method.  the default value is 12. */
/*                   if mxordn exceeds the default value, it will */
/*                   be reduced to the default value. */
/*                   mxordn is held constant during the problem. */

/* mxords  iwork(9)  the maximum order to be allowed for the stiff */
/*                   (bdf) method.  the default value is 5. */
/*                   if mxords exceeds the default value, it will */
/*                   be reduced to the default value. */
/*                   mxords is held constant during the problem. */
/* ----------------------------------------------------------------------- */
/* optional outputs. */

/* as optional additional output from lsoda, the variables listed */
/* below are quantities related to the performance of lsoda */
/* which are available to the user.  these are communicated by way of */
/* the work arrays, but also have internal mnemonic names as shown. */
/* except where stated otherwise, all of these outputs are defined */
/* on any successful return from lsoda, and on any return with */
/* istate = -1, -2, -4, -5, or -6.  on an illegal input return */
/* (istate = -3), they will be unchanged from their existing values */
/* (if any), except possibly for tolsf, lenrw, and leniw. */
/* on any error return, outputs relevant to the error will be defined, */
/* as noted below. */

/* name    location      meaning */

/* hu      rwork(11) the step size in t last used (successfully). */

/* hcur    rwork(12) the step size to be attempted on the next step. */

/* tcur    rwork(13) the current value of the independent variable */
/*                   which the solver has actually reached, i.e. the */
/*                   current internal mesh point in t.  on output, tcur */
/*                   will always be at least as far as the argument */
/*                   t, but may be farther (if interpolation was done). */

/* tolsf   rwork(14) a tolerance scale factor, greater than 1.0, */
/*                   computed when a request for too much accuracy was */
/*                   detected (istate = -3 if detected at the start of */
/*                   the problem, istate = -2 otherwise).  if itol is */
/*                   left unaltered but rtol and atol are uniformly */
/*                   scaled up by a factor of tolsf for the next call, */
/*                   then the solver is deemed likely to succeed. */
/*                   (the user may also ignore tolsf and alter the */
/*                   tolerance parameters in any other way appropriate.) */

/* tsw     rwork(15) the value of t at the time of the last method */
/*                   switch, if any. */

/* nst     iwork(11) the number of steps taken for the problem so far. */

/* nfe     iwork(12) the number of f evaluations for the problem so far. */

/* nje     iwork(13) the number of jacobian evaluations (and of matrix */
/*                   lu decompositions) for the problem so far. */

/* nqu     iwork(14) the method order last used (successfully). */

/* nqcur   iwork(15) the order to be attempted on the next step. */

/* imxer   iwork(16) the index of the component of largest magnitude in */
/*                   the weighted local error vector ( e(i)/ewt(i) ), */
/*                   on an error return with istate = -4 or -5. */

/* lenrw   iwork(17) the length of rwork actually required, assuming */
/*                   that the length of rwork is to be fixed for the */
/*                   rest of the problem, and that switching may occur. */
/*                   this is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* leniw   iwork(18) the length of iwork actually required, assuming */
/*                   that the length of iwork is to be fixed for the */
/*                   rest of the problem, and that switching may occur. */
/*                   this is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* mused   iwork(19) the method indicator for the last successful step.. */
/*                   1 means adams (nonstiff), 2 means bdf (stiff). */

/* mcur    iwork(20) the current method indicator.. */
/*                   1 means adams (nonstiff), 2 means bdf (stiff). */
/*                   this is the method to be attempted */
/*                   on the next step.  thus it differs from mused */
/*                   only if a method switch has just been made. */

/* the following two arrays are segments of the rwork array which */
/* may also be of interest to the user as optional outputs. */
/* for each array, the table below gives its internal name, */
/* its base address in rwork, and its description. */

/* name    base address      description */

/* yh      21             the nordsieck history array, of size nyh by */
/*                        (nqcur + 1), where nyh is the initial value */
/*                        of neq.  for j = 0,1,...,nqcur, column j+1 */
/*                        of yh contains hcur**j/factorial(j) times */
/*                        the j-th derivative of the interpolating */
/*                        polynomial currently representing the solution, */
/*                        evaluated at t = tcur. */

/* acor     lacor         array of size neq used for the accumulated */
/*         (from common   corrections on each step, scaled on output */
/*           as noted)    to represent the estimated local error in y */
/*                        on the last step.  this is the vector e in */
/*                        the description of the error control.  it is */
/*                        defined only on a successful return from lsoda. */
/*                        the base address lacor is obtained by */
/*                        including in the user-s program the */
/*                        following 3 lines.. */
/*                           double precision rls */
/*                           common /ls0001/ rls(218), ils(39) */
/*                           lacor = ils(5) */

/* ----------------------------------------------------------------------- */
/* part ii.  other routines callable. */

/* the following are optional calls which the user may make to */
/* gain additional capabilities in conjunction with lsoda. */
/* (the routines xsetun and xsetf are designed to conform to the */
/* slatec error handling package.) */

/*     form of call                  function */
/*   call xsetun(lun)          set the logical unit number, lun, for */
/*                             output of messages from lsoda, if */
/*                             the default is not desired. */
/*                             the default value of lun is 6. */

/*   call xsetf(mflag)         set a flag to control the printing of */
/*                             messages by lsoda. */
/*                             mflag = 0 means do not print. (danger.. */
/*                             this risks losing valuable information.) */
/*                             mflag = 1 means print (the default). */

/*                             either of the above calls may be made at */
/*                             any time and will take effect immediately. */

/*   call srcma(rsav,isav,job) saves and restores the contents of */
/*                             the internal common blocks used by */
/*                             lsoda (see part iii below). */
/*                             rsav must be a real array of length 240 */
/*                             or more, and isav must be an integer */
/*                             array of length 50 or more. */
/*                             job=1 means save common into rsav/isav. */
/*                             job=2 means restore common from rsav/isav. */
/*                                srcma is useful if one is */
/*                             interrupting a run and restarting */
/*                             later, or alternating between two or */
/*                             more problems solved with lsoda. */

/*   call intdy(,,,,,)         provide derivatives of y, of various */
/*        (see below)          orders, at a specified point t, if */
/*                             desired.  it may be called only after */
/*                             a successful return from lsoda. */

/* the detailed instructions for using intdy are as follows. */
/* the form of the call is.. */

/*   call intdy (t, k, rwork(21), nyh, dky, iflag) */

/* the input parameters are.. */

/* t         = value of independent variable where answers are desired */
/*             (normally the same as the t last returned by lsoda). */
/*             for valid results, t must lie between tcur - hu and tcur. */
/*             (see optional outputs for tcur and hu.) */
/* k         = integer order of the derivative desired.  k must satisfy */
/*             0 .le. k .le. nqcur, where nqcur is the current order */
/*             (see optional outputs).  the capability corresponding */
/*             to k = 0, i.e. computing y(t), is already provided */
/*             by lsoda directly.  since nqcur .ge. 1, the first */
/*             derivative dy/dt is always available with intdy. */
/* rwork(21) = the base address of the history array yh. */
/* nyh       = column length of yh, equal to the initial value of neq. */

/* the output parameters are.. */

/* dky       = a real array of length neq containing the computed value */
/*             of the k-th derivative of y(t). */
/* iflag     = integer flag, returned as 0 if k and t were legal, */
/*             -1 if k was illegal, and -2 if t was illegal. */
/*             on an error return, a message is also written. */
/* ----------------------------------------------------------------------- */
/* part iii.  common blocks. */

/* if lsoda is to be used in an overlay situation, the user */
/* must declare, in the primary overlay, the variables in.. */
/*   (1) the call sequence to lsoda, */
/*   (2) the three internal common blocks */
/*         /ls0001/  of length  257  (218 double precision words */
/*                         followed by 39 integer words), */
/*         /lsa001/  of length  31    (22 double precision words */
/*                         followed by  9 integer words), */
/*         /eh0001/  of length  2 (integer words). */

/* if lsoda is used on a system in which the contents of internal */
/* common blocks are not preserved between calls, the user should */
/* declare the above common blocks in his main program to insure */
/* that their contents are preserved. */

/* if the solution of a given problem by lsoda is to be interrupted */
/* and then later continued, such as when restarting an interrupted run */
/* or alternating between two or more problems, the user should save, */
/* following the return from the last lsoda call prior to the */
/* interruption, the contents of the call sequence variables and the */
/* internal common blocks, and later restore these values before the */
/* next lsoda call for that problem.  to save and restore the common */
/* blocks, use subroutine srcma (see part ii above). */

/* ----------------------------------------------------------------------- */
/* part iv.  optionally replaceable solver routines. */

/* below is a description of a routine in the lsoda package which */
/* relates to the measurement of errors, and can be */
/* replaced by a user-supplied version, if desired.  however, since such */
/* a replacement may have a major impact on performance, it should be */
/* done only when absolutely necessary, and only with great caution. */
/* (note.. the means by which the package version of a routine is */
/* superseded by the user-s version may be system-dependent.) */

/* (a) ewset. */
/* the following subroutine is called just before each internal */
/* integration step, and sets the array of error weights, ewt, as */
/* described under itol/rtol/atol above.. */
/*     subroutine ewset (neq, itol, rtol, atol, ycur, ewt) */
/* where neq, itol, rtol, and atol are as in the lsoda call sequence, */
/* ycur contains the current dependent variable vector, and */
/* ewt is the array of weights set by ewset. */

/* if the user supplies this subroutine, it must return in ewt(i) */
/* (i = 1,...,neq) a positive quantity suitable for comparing errors */
/* in y(i) to.  the ewt array returned by ewset is passed to the */
/* vmnorm routine, and also used by lsoda in the computation */
/* of the optional output imxer, and the increments for difference */
/* quotient jacobians. */

/* in the user-supplied version of ewset, it may be desirable to use */
/* the current values of derivatives of y.  derivatives up to order nq */
/* are available from the history array yh, described above under */
/* optional outputs.  in ewset, yh is identical to the ycur array, */
/* extended to nq + 1 columns with a column length of nyh and scale */
/* factors of h**j/factorial(j).  on the first call for the problem, */
/* given by nst = 0, nq is 1 and h is temporarily set to 1.0. */
/* the quantities nq, nyh, h, and nst can be obtained by including */
/* in ewset the statements.. */
/*     double precision h, rls */
/*     common /ls0001/ rls(218),ils(39) */
/*     nq = ils(35) */
/*     nyh = ils(14) */
/*     nst = ils(36) */
/*     h = rls(212) */
/* thus, for example, the current value of dy/dt can be obtained as */
/* ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is */
/* unnecessary when nst = 0). */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* other routines in the lsoda package. */

/* in addition to subroutine lsoda, the lsoda package includes the */
/* following subroutines and function routines.. */
/*  intdy    computes an interpolated value of the y vector at t = tout. */
/*  stoda    is the core integrator, which does one step of the */
/*           integration and the associated error control. */
/*  cfode    sets all method coefficients and test constants. */
/*  prja     computes and preprocesses the jacobian matrix j = df/dy */
/*           and the newton iteration matrix p = i - h*l0*j. */
/*  solsy    manages solution of linear system in chord iteration. */
/*  ewset    sets the error weight vector ewt before each step. */
/*  vmnorm   computes the weighted max-norm of a vector. */
/*  fnorm    computes the norm of a full matrix consistent with the */
/*           weighted max-norm on vectors. */
/*  bnorm    computes the norm of a band matrix consistent with the */
/*           weighted max-norm on vectors. */
/*  srcma    is a user-callable routine to save and restore */
/*           the contents of the internal common blocks. */
/*  dgetrf and dgetrs   are routines from lapack for solving full */
/*           systems of linear algebraic equations. */
/*  dgbtrf and dgbtrs   are routines from lapack for solving banded */
/*           linear systems. */
/*  daxpy, dscal, idamax, and ddot   are basic linear algebra modules */
/*           (blas) used by the above linpack routines. */
/*  d1mach   computes the unit roundoff in a machine-independent manner. */
/*  xerrwv, xsetun, and xsetf   handle the printing of all error */
/*           messages and warnings.  xerrwv is machine-dependent. */
/* note..  vmnorm, fnorm, bnorm, idamax, ddot, and d1mach are function */
/* routines.  all the others are subroutines. */

/* the intrinsic and external routines used by lsoda are.. */
/* dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write. */

/* a block data subprogram is also included with the package, */
/* for loading some of the variables in internal common. */

/* ----------------------------------------------------------------------- */
/* the following card is for optimized compilation on lll compilers. */
/* lll. optimize */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* the following two internal common blocks contain */
/* (a) variables which are local to any subroutine but whose values must */
/*     be preserved between calls to the routine (own variables), and */
/* (b) variables which are communicated between subroutines. */
/* the structure of each block is as follows..  all real variables are */
/* listed first, followed by all integers.  within each type, the */
/* variables are grouped with those local to subroutine lsoda first, */
/* then those local to subroutine stoda, and finally those used */
/* for communication.  the block ls0001 is declared in subroutines */
/* lsoda, intdy, stoda, prja, and solsy.  the block lsa001 is declared */
/* in subroutines lsoda, stoda, and prja.  groups of variables are */
/* replaced by dummy arrays in the common declarations in routines */
/* where those variables are not used. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --neq;
    --y;
    --rtol;
    --atol;
    --rwork;
    --iwork;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* block a. */
/* this code block is executed on every call. */
/* it tests istate and itask for legality and branches appropriately. */
/* if istate .gt. 1 but the flag init shows that initialization has */
/* not yet been done, an error return occurs. */
/* if istate = 1 and tout = t, jump to block g and return immediately. */
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
    if (ls0001_1.init == 0) {
	goto L603;
    }
    if (*istate == 2) {
	goto L200;
    }
    goto L20;
L10:
    ls0001_1.init = 0;
    if (*tout == *t) {
	goto L430;
    }
L20:
    ls0001_1.ntrep = 0;
/* ----------------------------------------------------------------------- */
/* block b. */
/* the next code block is executed for the initial call (istate = 1), */
/* or for a continuation call with parameter changes (istate = 3). */
/* it contains checking of all inputs and various initializations. */

/* first check legality of the non-optional inputs neq, itol, iopt, */
/* jt, ml, and mu. */
/* ----------------------------------------------------------------------- */
    if (neq[1] <= 0) {
	goto L604;
    }
    if (*istate == 1) {
	goto L25;
    }
    if (neq[1] > ls0001_1.n) {
	goto L605;
    }
L25:
    ls0001_1.n = neq[1];
    if (*itol < 1 || *itol > 4) {
	goto L606;
    }
    if (*iopt < 0 || *iopt > 1) {
	goto L607;
    }
    if (*jt == 3 || *jt < 1 || *jt > 5) {
	goto L608;
    }
    lsa001_1.jtyp = *jt;
    if (*jt <= 2) {
	goto L30;
    }
    ml = iwork[1];
    mu = iwork[2];
    if (ml < 0 || ml >= ls0001_1.n) {
	goto L609;
    }
    if (mu < 0 || mu >= ls0001_1.n) {
	goto L610;
    }
L30:
/* next process and check the optional inputs. -------------------------- */
    if (*iopt == 1) {
	goto L40;
    }
    lsa001_1.ixpr = 0;
    ls0001_1.mxstep = mxstp0;
    ls0001_1.mxhnil = mxhnl0;
    ls0001_1.hmxi = 0.;
    ls0001_1.hmin = 0.;
    if (*istate != 1) {
	goto L60;
    }
    h0 = 0.;
    lsa001_1.mxordn = mord[0];
    lsa001_1.mxords = mord[1];
    goto L60;
L40:
    lsa001_1.ixpr = iwork[5];
    if (lsa001_1.ixpr < 0 || lsa001_1.ixpr > 1) {
	goto L611;
    }
    ls0001_1.mxstep = iwork[6];
    if (ls0001_1.mxstep < 0) {
	goto L612;
    }
    if (ls0001_1.mxstep == 0) {
	ls0001_1.mxstep = mxstp0;
    }
    ls0001_1.mxhnil = iwork[7];
    if (ls0001_1.mxhnil < 0) {
	goto L613;
    }
    if (ls0001_1.mxhnil == 0) {
	ls0001_1.mxhnil = mxhnl0;
    }
    if (*istate != 1) {
	goto L50;
    }
    h0 = rwork[5];
    lsa001_1.mxordn = iwork[8];
    if (lsa001_1.mxordn < 0) {
	goto L628;
    }
    if (lsa001_1.mxordn == 0) {
	lsa001_1.mxordn = 100;
    }
    lsa001_1.mxordn = min(lsa001_1.mxordn,mord[0]);
    lsa001_1.mxords = iwork[9];
    if (lsa001_1.mxords < 0) {
	goto L629;
    }
    if (lsa001_1.mxords == 0) {
	lsa001_1.mxords = 100;
    }
    lsa001_1.mxords = min(lsa001_1.mxords,mord[1]);
    if ((*tout - *t) * h0 < 0.) {
	goto L614;
    }
L50:
    hmax = rwork[6];
    if (hmax < 0.) {
	goto L615;
    }
    ls0001_1.hmxi = 0.;
    if (hmax > 0.) {
	ls0001_1.hmxi = 1. / hmax;
    }
    ls0001_1.hmin = rwork[7];
    if (ls0001_1.hmin < 0.) {
	goto L616;
    }
/* ----------------------------------------------------------------------- */
/* set work array pointers and check lengths lrw and liw. */
/* if istate = 1, meth is initialized to 1 here to facilitate the */
/* checking of work space lengths. */
/* pointers to segments of rwork and iwork are named by prefixing l to */
/* the name of the segment.  e.g., the segment yh starts at rwork(lyh). */
/* segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor. */
/* if the lengths provided are insufficient for the current method, */
/* an error return occurs.  this is treated as illegal input on the */
/* first call, but as a problem interruption with istate = -7 on a */
/* continuation call.  if the lengths are sufficient for the current */
/* method but not for both methods, a warning message is sent. */
/* ----------------------------------------------------------------------- */
L60:
    if (*istate == 1) {
	ls0001_1.meth = 1;
    }
    if (*istate == 1) {
	ls0001_1.nyh = ls0001_1.n;
    }
    ls0001_1.lyh = 21;
    len1n = (lsa001_1.mxordn + 1) * ls0001_1.nyh + 20;
    len1s = (lsa001_1.mxords + 1) * ls0001_1.nyh + 20;
    ls0001_1.lwm = len1s + 1;
    if (*jt <= 2) {
	lenwm = ls0001_1.n * ls0001_1.n + 2;
    }
    if (*jt >= 4) {
	lenwm = ((ml << 1) + mu + 1) * ls0001_1.n + 2;
    }
    len1s += lenwm;
    len1c = len1n;
    if (ls0001_1.meth == 2) {
	len1c = len1s;
    }
    len1 = max(len1n,len1s);
    len2 = ls0001_1.n * 3;
    lenrw = len1 + len2;
    lenrwn = len1n + len2;
    lenrws = len1s + len2;
    lenrwc = len1c + len2;
    iwork[17] = lenrw;
    ls0001_1.liwm = 1;
    leniw = ls0001_1.n + 20;
    leniwc = 20;
    if (ls0001_1.meth == 2) {
	leniwc = leniw;
    }
    iwork[18] = leniw;
    if (*istate == 1 && *lrw < lenrwc) {
	goto L617;
    }
    if (*istate == 1 && *liw < leniwc) {
	goto L618;
    }
    if (*istate == 3 && *lrw < lenrwc) {
	goto L550;
    }
    if (*istate == 3 && *liw < leniwc) {
	goto L555;
    }
    ls0001_1.lewt = len1 + 1;
    lsa001_1.insufr = 0;
    if (*lrw >= lenrw) {
	goto L65;
    }
    lsa001_1.insufr = 2;
    ls0001_1.lewt = len1c + 1;
    xerrwv_("lsoda--  warning.. rwork length is sufficient for now, but  ", &
	    c__60, &c__103, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42,
	     (ftnlen)60);
    xerrwv_("      may not be later.  integration will proceed anyway.   ", &
	    c__60, &c__103, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42,
	     (ftnlen)60);
    xerrwv_("      length needed is lenrw = i1, while lrw = i2.", &c__50, &
	    c__103, &c__0, &c__2, &lenrw, lrw, &c__0, &c_b42, &c_b42, (ftnlen)
	    50);
L65:
    ls0001_1.lsavf = ls0001_1.lewt + ls0001_1.n;
    ls0001_1.lacor = ls0001_1.lsavf + ls0001_1.n;
    lsa001_1.insufi = 0;
    if (*liw >= leniw) {
	goto L70;
    }
    lsa001_1.insufi = 2;
    xerrwv_("lsoda--  warning.. iwork length is sufficient for now, but  ", &
	    c__60, &c__104, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42,
	     (ftnlen)60);
    xerrwv_("      may not be later.  integration will proceed anyway.   ", &
	    c__60, &c__104, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42,
	     (ftnlen)60);
    xerrwv_("      length needed is leniw = i1, while liw = i2.", &c__50, &
	    c__104, &c__0, &c__2, &leniw, liw, &c__0, &c_b42, &c_b42, (ftnlen)
	    50);
L70:
/* check rtol and atol for legality. ------------------------------------ */
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol >= 3) {
	    rtoli = rtol[i__];
	}
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	if (rtoli < 0.) {
	    goto L619;
	}
	if (atoli < 0.) {
	    goto L620;
	}
/* L75: */
    }
    if (*istate == 1) {
	goto L100;
    }
/* if istate = 3, set flag to signal parameter changes to stoda. -------- */
    ls0001_1.jstart = -1;
    if (ls0001_1.n == ls0001_1.nyh) {
	goto L200;
    }
/* neq was reduced.  zero part of yh to avoid undefined references. ----- */
    i1 = ls0001_1.lyh + ls0001_1.l * ls0001_1.nyh;
    i2 = ls0001_1.lyh + (ls0001_1.maxord + 1) * ls0001_1.nyh - 1;
    if (i1 > i2) {
	goto L200;
    }
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
/* L95: */
	rwork[i__] = 0.;
    }
    goto L200;
/* ----------------------------------------------------------------------- */
/* block c. */
/* the next block is for the initial call only (istate = 1). */
/* it contains all remaining initializations, the initial call to f, */
/* and the calculation of the initial step size. */
/* the error weights in ewt are inverted after being loaded. */
/* ----------------------------------------------------------------------- */
L100:
    ls0001_1.uround = d1mach_(&c__4);
    ls0001_1.tn = *t;
    lsa001_1.tsw = *t;
    ls0001_1.maxord = lsa001_1.mxordn;
    if (*itask != 4 && *itask != 5) {
	goto L110;
    }
    tcrit = rwork[1];
    if ((tcrit - *tout) * (*tout - *t) < 0.) {
	goto L625;
    }
    if (h0 != 0. && (*t + h0 - tcrit) * h0 > 0.) {
	h0 = tcrit - *t;
    }
L110:
    ls0001_1.jstart = 0;
    ls0001_1.nhnil = 0;
    ls0001_1.nst = 0;
    ls0001_1.nje = 0;
    ls0001_1.nslast = 0;
    ls0001_1.hu = 0.;
    ls0001_1.nqu = 0;
    lsa001_1.mused = 0;
    ls0001_1.miter = 0;
    ls0001_1.ccmax = .3;
    ls0001_1.maxcor = 3;
    ls0001_1.msbp = 20;
    ls0001_1.mxncf = 10;
/* initial call to f.  (lf0 points to yh(*,2).) ------------------------- */
    lf0 = ls0001_1.lyh + ls0001_1.nyh;
    srcma_(rsav, isav, &c__1);
    (*f)(&neq[1], t, &y[1], &rwork[lf0]);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    srcma_(rsav, isav, &c__2);
    ls0001_1.nfe = 1;
/* load the initial value vector in yh. --------------------------------- */
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L115: */
	rwork[i__ + ls0001_1.lyh - 1] = y[i__];
    }
/* load and invert the ewt array.  (h is temporarily set to 1.0.) ------- */
    ls0001_1.nq = 1;
    ls0001_1.h__ = 1.;
    ewset_(&ls0001_1.n, itol, &rtol[1], &atol[1], &rwork[ls0001_1.lyh], &
	    rwork[ls0001_1.lewt]);
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + ls0001_1.lewt - 1] <= 0.) {
	    goto L621;
	}
/* L120: */
	rwork[i__ + ls0001_1.lewt - 1] = 1. / rwork[i__ + ls0001_1.lewt - 1];
    }
/* ----------------------------------------------------------------------- */
/* the coding below computes the step size, h0, to be attempted on the */
/* first step, unless the user has supplied a value for this. */
/* first check that tout - t differs significantly from zero. */
/* a scalar tolerance quantity tol is computed, as max(rtol(i)) */
/* if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted */
/* so as to be between 100*uround and 1.0e-3. */
/* then the computed value h0 is given by.. */

/*   h0**(-2)  =  1./(tol * w0**2)  +  tol * (norm(f))**2 */

/* where   w0     = max ( abs(t), abs(tout) ), */
/*         f      = the initial value of the vector f(t,y), and */
/*         norm() = the weighted vector norm used throughout, given by */
/*                  the vmnorm function routine, and weighted by the */
/*                  tolerances initially loaded into the ewt array. */
/* the sign of h0 is inferred from the initial values of tout and t. */
/* abs(h0) is made .le. abs(tout-t) in any case. */
/* ----------------------------------------------------------------------- */
    if (h0 != 0.) {
	goto L180;
    }
    tdist = (d__1 = *tout - *t, abs(d__1));
/* Computing MAX */
    d__1 = abs(*t), d__2 = abs(*tout);
    w0 = max(d__1,d__2);
    if (tdist < ls0001_1.uround * 2. * w0) {
	goto L622;
    }
    tol = rtol[1];
    if (*itol <= 2) {
	goto L140;
    }
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
/* Computing MAX */
	d__1 = tol, d__2 = rtol[i__];
	tol = max(d__1,d__2);
    }
L140:
    if (tol > 0.) {
	goto L160;
    }
    atoli = atol[1];
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	ayi = (d__1 = y[i__], abs(d__1));
	if (ayi != 0.) {
/* Computing MAX */
	    d__1 = tol, d__2 = atoli / ayi;
	    tol = max(d__1,d__2);
	}
/* L150: */
    }
L160:
/* Computing MAX */
    d__1 = tol, d__2 = ls0001_1.uround * 100.;
    tol = max(d__1,d__2);
    tol = min(tol,.001);
    sum = vmnorm_(&ls0001_1.n, &rwork[lf0], &rwork[ls0001_1.lewt]);
/* Computing 2nd power */
    d__1 = sum;
    sum = 1. / (tol * w0 * w0) + tol * (d__1 * d__1);
    h0 = 1. / sqrt(sum);
    h0 = min(h0,tdist);
    d__1 = *tout - *t;
    h0 = d_sign(&h0, &d__1);
/* adjust h0 if necessary to meet hmax bound. --------------------------- */
L180:
    rh = abs(h0) * ls0001_1.hmxi;
    if (rh > 1.) {
	h0 /= rh;
    }
/* load h with h0 and scale yh(*,2) by h0. ------------------------------ */
    ls0001_1.h__ = h0;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L190: */
	rwork[i__ + lf0 - 1] = h0 * rwork[i__ + lf0 - 1];
    }
    goto L270;
/* ----------------------------------------------------------------------- */
/* block d. */
/* the next code block is for continuation calls only (istate = 2 or 3) */
/* and is to check stop conditions before taking a step. */
/* ----------------------------------------------------------------------- */
L200:
    ls0001_1.nslast = ls0001_1.nst;
    switch (*itask) {
	case 1:  goto L210;
	case 2:  goto L250;
	case 3:  goto L220;
	case 4:  goto L230;
	case 5:  goto L240;
    }
L210:
    if ((ls0001_1.tn - *tout) * ls0001_1.h__ < 0.) {
	goto L250;
    }
    intdy_(tout, &c__0, &rwork[ls0001_1.lyh], &ls0001_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L220:
    tp = ls0001_1.tn - ls0001_1.hu * (ls0001_1.uround * 100. + 1.);
    if ((tp - *tout) * ls0001_1.h__ > 0.) {
	goto L623;
    }
    if ((ls0001_1.tn - *tout) * ls0001_1.h__ < 0.) {
	goto L250;
    }
    *t = ls0001_1.tn;
    goto L400;
L230:
    tcrit = rwork[1];
    if ((ls0001_1.tn - tcrit) * ls0001_1.h__ > 0.) {
	goto L624;
    }
    if ((tcrit - *tout) * ls0001_1.h__ < 0.) {
	goto L625;
    }
    if ((ls0001_1.tn - *tout) * ls0001_1.h__ < 0.) {
	goto L245;
    }
    intdy_(tout, &c__0, &rwork[ls0001_1.lyh], &ls0001_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L240:
    tcrit = rwork[1];
    if ((ls0001_1.tn - tcrit) * ls0001_1.h__ > 0.) {
	goto L624;
    }
L245:
    hmx = abs(ls0001_1.tn) + abs(ls0001_1.h__);
    ihit = (d__1 = ls0001_1.tn - tcrit, abs(d__1)) <= ls0001_1.uround * 100. *
	     hmx;
    if (ihit) {
	*t = tcrit;
    }
    if (ihit) {
	goto L400;
    }
    tnext = ls0001_1.tn + ls0001_1.h__ * (ls0001_1.uround * 4. + 1.);
    if ((tnext - tcrit) * ls0001_1.h__ <= 0.) {
	goto L250;
    }
    ls0001_1.h__ = (tcrit - ls0001_1.tn) * (1. - ls0001_1.uround * 4.);
    if (*istate == 2 && ls0001_1.jstart >= 0) {
	ls0001_1.jstart = -2;
    }
/* ----------------------------------------------------------------------- */
/* block e. */
/* the next block is normally executed for all calls and contains */
/* the call to the one-step core integrator stoda. */

/* this is a looping point for the integration steps. */

/* first check for too many steps being taken, update ewt (if not at */
/* start of problem), check for too much accuracy being requested, and */
/* check for h below the roundoff level in t. */
/* ----------------------------------------------------------------------- */
L250:
    if (ls0001_1.meth == lsa001_1.mused) {
	goto L255;
    }
    if (lsa001_1.insufr == 1) {
	goto L550;
    }
    if (lsa001_1.insufi == 1) {
	goto L555;
    }
L255:
    if (ls0001_1.nst - ls0001_1.nslast >= ls0001_1.mxstep) {
	goto L500;
    }
    ewset_(&ls0001_1.n, itol, &rtol[1], &atol[1], &rwork[ls0001_1.lyh], &
	    rwork[ls0001_1.lewt]);
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + ls0001_1.lewt - 1] <= 0.) {
	    goto L510;
	}
/* L260: */
	rwork[i__ + ls0001_1.lewt - 1] = 1. / rwork[i__ + ls0001_1.lewt - 1];
    }
L270:
    tolsf = ls0001_1.uround * vmnorm_(&ls0001_1.n, &rwork[ls0001_1.lyh], &
	    rwork[ls0001_1.lewt]);
    if (tolsf <= .01) {
	goto L280;
    }
    tolsf *= 200.;
    if (ls0001_1.nst == 0) {
	goto L626;
    }
    goto L520;
L280:
    if (ls0001_1.tn + ls0001_1.h__ != ls0001_1.tn) {
	goto L290;
    }
    ++ls0001_1.nhnil;
    if (ls0001_1.nhnil > ls0001_1.mxhnil) {
	goto L290;
    }
    xerrwv_("lsoda--  warning..internal t (=r1) and h (=r2) are", &c__50, &
	    c__101, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      such that in the machine, t + h = t on the next step  ", &
	    c__60, &c__101, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42,
	     (ftnlen)60);
    xerrwv_("      (h = step size). solver will continue anyway", &c__50, &
	    c__101, &c__0, &c__0, &c__0, &c__0, &c__2, &ls0001_1.tn, &
	    ls0001_1.h__, (ftnlen)50);
    if (ls0001_1.nhnil < ls0001_1.mxhnil) {
	goto L290;
    }
    xerrwv_("lsoda--  above warning has been issued i1 times.  ", &c__50, &
	    c__102, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      it will not be issued again for this problem", &c__50, &
	    c__102, &c__0, &c__1, &ls0001_1.mxhnil, &c__0, &c__0, &c_b42, &
	    c_b42, (ftnlen)50);
L290:
/* ----------------------------------------------------------------------- */
/*     call stoda(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prja,solsy) */
/* ----------------------------------------------------------------------- */
    stoda_(&neq[1], &y[1], &rwork[ls0001_1.lyh], &ls0001_1.nyh, &rwork[
	    ls0001_1.lyh], &rwork[ls0001_1.lewt], &rwork[ls0001_1.lsavf], &
	    rwork[ls0001_1.lacor], &rwork[ls0001_1.lwm], &iwork[ls0001_1.liwm]
	    , (S_fp)f, (U_fp)jac, (U_fp)prja_, (U_fp)solsy_);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    kgo = 1 - ls0001_1.kflag;
    switch (kgo) {
	case 1:  goto L300;
	case 2:  goto L530;
	case 3:  goto L540;
    }
/* ----------------------------------------------------------------------- */
/* block f. */
/* the following block handles the case of a successful return from the */
/* core integrator (kflag = 0). */
/* if a method switch was just made, record tsw, reset maxord, */
/* set jstart to -1 to signal stoda to complete the switch, */
/* and do extra printing of data if ixpr = 1. */
/* then, in any case, check for stop conditions. */
/* ----------------------------------------------------------------------- */
L300:
    ls0001_1.init = 1;
    if (ls0001_1.meth == lsa001_1.mused) {
	goto L310;
    }
    lsa001_1.tsw = ls0001_1.tn;
    ls0001_1.maxord = lsa001_1.mxordn;
    if (ls0001_1.meth == 2) {
	ls0001_1.maxord = lsa001_1.mxords;
    }
    if (ls0001_1.meth == 2) {
	rwork[ls0001_1.lwm] = sqrt(ls0001_1.uround);
    }
    lsa001_1.insufr = min(lsa001_1.insufr,1);
    lsa001_1.insufi = min(lsa001_1.insufi,1);
    ls0001_1.jstart = -1;
    if (lsa001_1.ixpr == 0) {
	goto L310;
    }
    if (ls0001_1.meth == 2) {
	xerrwv_("lsoda-- a switch to the bdf (stiff) method has occurred     "
		, &c__60, &c__105, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, 
		&c_b42, (ftnlen)60);
    }
    if (ls0001_1.meth == 1) {
	xerrwv_("lsoda-- a switch to the adams (nonstiff) method has occurred"
		, &c__60, &c__106, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, 
		&c_b42, (ftnlen)60);
    }
    xerrwv_("     at t = r1,  tentative step size h = r2,  step nst = i1 ", &
	    c__60, &c__107, &c__0, &c__1, &ls0001_1.nst, &c__0, &c__2, &
	    ls0001_1.tn, &ls0001_1.h__, (ftnlen)60);
L310:
    switch (*itask) {
	case 1:  goto L320;
	case 2:  goto L400;
	case 3:  goto L330;
	case 4:  goto L340;
	case 5:  goto L350;
    }
/* itask = 1.  if tout has been reached, interpolate. ------------------- */
L320:
    if ((ls0001_1.tn - *tout) * ls0001_1.h__ < 0.) {
	goto L250;
    }
    intdy_(tout, &c__0, &rwork[ls0001_1.lyh], &ls0001_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
/* itask = 3.  jump to exit if tout was reached. ------------------------ */
L330:
    if ((ls0001_1.tn - *tout) * ls0001_1.h__ >= 0.) {
	goto L400;
    }
    goto L250;
/* itask = 4.  see if tout or tcrit was reached.  adjust h if necessary. */
L340:
    if ((ls0001_1.tn - *tout) * ls0001_1.h__ < 0.) {
	goto L345;
    }
    intdy_(tout, &c__0, &rwork[ls0001_1.lyh], &ls0001_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
L345:
    hmx = abs(ls0001_1.tn) + abs(ls0001_1.h__);
    ihit = (d__1 = ls0001_1.tn - tcrit, abs(d__1)) <= ls0001_1.uround * 100. *
	     hmx;
    if (ihit) {
	goto L400;
    }
    tnext = ls0001_1.tn + ls0001_1.h__ * (ls0001_1.uround * 4. + 1.);
    if ((tnext - tcrit) * ls0001_1.h__ <= 0.) {
	goto L250;
    }
    ls0001_1.h__ = (tcrit - ls0001_1.tn) * (1. - ls0001_1.uround * 4.);
    if (ls0001_1.jstart >= 0) {
	ls0001_1.jstart = -2;
    }
    goto L250;
/* itask = 5.  see if tcrit was reached and jump to exit. --------------- */
L350:
    hmx = abs(ls0001_1.tn) + abs(ls0001_1.h__);
    ihit = (d__1 = ls0001_1.tn - tcrit, abs(d__1)) <= ls0001_1.uround * 100. *
	     hmx;
/* ----------------------------------------------------------------------- */
/* block g. */
/* the following block handles all successful returns from lsoda. */
/* if itask .ne. 1, y is loaded from yh and t is set accordingly. */
/* istate is set to 2, the illegal input counter is zeroed, and the */
/* optional outputs are loaded into the work arrays before returning. */
/* if istate = 1 and tout = t, there is a return with no action taken, */
/* except that if this has happened repeatedly, the run is terminated. */
/* ----------------------------------------------------------------------- */
L400:
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	y[i__] = rwork[i__ + ls0001_1.lyh - 1];
    }
    *t = ls0001_1.tn;
    if (*itask != 4 && *itask != 5) {
	goto L420;
    }
    if (ihit) {
	*t = tcrit;
    }
L420:
    *istate = 2;
    ls0001_1.illin = 0;
    rwork[11] = ls0001_1.hu;
    rwork[12] = ls0001_1.h__;
    rwork[13] = ls0001_1.tn;
    rwork[15] = lsa001_1.tsw;
    iwork[11] = ls0001_1.nst;
    iwork[12] = ls0001_1.nfe;
    iwork[13] = ls0001_1.nje;
    iwork[14] = ls0001_1.nqu;
    iwork[15] = ls0001_1.nq;
    iwork[19] = lsa001_1.mused;
    iwork[20] = ls0001_1.meth;
    return 0;

L430:
    ++ls0001_1.ntrep;
    if (ls0001_1.ntrep < 5) {
	return 0;
    }
    xerrwv_("lsoda--  repeated calls with istate = 1 and tout = t (=r1)  ", &
	    c__60, &c__301, &c__0, &c__0, &c__0, &c__0, &c__1, t, &c_b42, (
	    ftnlen)60);
    goto L800;
/* ----------------------------------------------------------------------- */
/* block h. */
/* the following block handles all unsuccessful returns other than */
/* those for illegal input.  first the error message routine is called. */
/* if there was an error test or convergence test failure, imxer is set. */
/* then y is loaded from yh, t is set to tn, and the illegal input */
/* counter illin is set to 0.  the optional outputs are loaded into */
/* the work arrays before returning. */
/* ----------------------------------------------------------------------- */
/* the maximum number of steps was taken before reaching tout. ---------- */
/* Error message removed, see gh-7888 */
L500:
    *istate = -1;
    goto L580;
/* ewt(i) .le. 0.0 for some i (not at start of problem). ---------------- */
L510:
    ewti = rwork[ls0001_1.lewt + i__ - 1];
    xerrwv_("lsoda--  at t (=r1), ewt(i1) has become r2 .le. 0.", &c__50, &
	    c__202, &c__0, &c__1, &i__, &c__0, &c__2, &ls0001_1.tn, &ewti, (
	    ftnlen)50);
    *istate = -6;
    goto L580;
/* too much accuracy requested for machine precision. ------------------- */
L520:
    xerrwv_("lsoda--  at t (=r1), too much accuracy requested  ", &c__50, &
	    c__203, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      for precision of machine..  see tolsf (=r2) ", &c__50, &
	    c__203, &c__0, &c__0, &c__0, &c__0, &c__2, &ls0001_1.tn, &tolsf, (
	    ftnlen)50);
    rwork[14] = tolsf;
    *istate = -2;
    goto L580;
/* kflag = -1.  error test failed repeatedly or with abs(h) = hmin. ----- */
L530:
    xerrwv_("lsoda--  at t(=r1) and step size h(=r2), the error", &c__50, &
	    c__204, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      test failed repeatedly or with abs(h) = hmin", &c__50, &
	    c__204, &c__0, &c__0, &c__0, &c__0, &c__2, &ls0001_1.tn, &
	    ls0001_1.h__, (ftnlen)50);
    *istate = -4;
    goto L560;
/* kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ---- */
L540:
    xerrwv_("lsoda--  at t (=r1) and step size h (=r2), the    ", &c__50, &
	    c__205, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      corrector convergence failed repeatedly     ", &c__50, &
	    c__205, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      or with abs(h) = hmin   ", &c__30, &c__205, &c__0, &c__0, &
	    c__0, &c__0, &c__2, &ls0001_1.tn, &ls0001_1.h__, (ftnlen)30);
    *istate = -5;
    goto L560;
/* rwork length too small to proceed. ----------------------------------- */
L550:
    xerrwv_("lsoda--  at current t(=r1), rwork length too small", &c__50, &
	    c__206, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      to proceed.  the integration was otherwise successful.", &
	    c__60, &c__206, &c__0, &c__0, &c__0, &c__0, &c__1, &ls0001_1.tn, &
	    c_b42, (ftnlen)60);
    *istate = -7;
    goto L580;
/* iwork length too small to proceed. ----------------------------------- */
L555:
    xerrwv_("lsoda--  at current t(=r1), iwork length too small", &c__50, &
	    c__207, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    xerrwv_("      to proceed.  the integration was otherwise successful.", &
	    c__60, &c__207, &c__0, &c__0, &c__0, &c__0, &c__1, &ls0001_1.tn, &
	    c_b42, (ftnlen)60);
    *istate = -7;
    goto L580;
/* compute imxer if relevant. ------------------------------------------- */
L560:
    big = 0.;
    imxer = 1;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	size = (d__1 = rwork[i__ + ls0001_1.lacor - 1] * rwork[i__ + 
		ls0001_1.lewt - 1], abs(d__1));
	if (big >= size) {
	    goto L570;
	}
	big = size;
	imxer = i__;
L570:
	;
    }
    iwork[16] = imxer;
/* set y vector, t, illin, and optional outputs. ------------------------ */
L580:
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L590: */
	y[i__] = rwork[i__ + ls0001_1.lyh - 1];
    }
    *t = ls0001_1.tn;
    ls0001_1.illin = 0;
    rwork[11] = ls0001_1.hu;
    rwork[12] = ls0001_1.h__;
    rwork[13] = ls0001_1.tn;
    rwork[15] = lsa001_1.tsw;
    iwork[11] = ls0001_1.nst;
    iwork[12] = ls0001_1.nfe;
    iwork[13] = ls0001_1.nje;
    iwork[14] = ls0001_1.nqu;
    iwork[15] = ls0001_1.nq;
    iwork[19] = lsa001_1.mused;
    iwork[20] = ls0001_1.meth;
    return 0;
/* ----------------------------------------------------------------------- */
/* block i. */
/* the following block handles all error returns due to illegal input */
/* (istate = -3), as detected before calling the core integrator. */
/* first the error message routine is called.  then if there have been */
/* 5 consecutive such returns just before this call to the solver, */
/* the run is halted. */
/* ----------------------------------------------------------------------- */
L601:
    xerrwv_("lsoda--  istate (=i1) illegal ", &c__30, &c__1, &c__0, &c__1, 
	    istate, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L602:
    xerrwv_("lsoda--  itask (=i1) illegal  ", &c__30, &c__2, &c__0, &c__1, 
	    itask, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L603:
    xerrwv_("lsoda--  istate .gt. 1 but lsoda not initialized  ", &c__50, &
	    c__3, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)
	    50);
    goto L700;
L604:
    xerrwv_("lsoda--  neq (=i1) .lt. 1     ", &c__30, &c__4, &c__0, &c__1, &
	    neq[1], &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L605:
    xerrwv_("lsoda--  istate = 3 and neq increased (i1 to i2)  ", &c__50, &
	    c__5, &c__0, &c__2, &ls0001_1.n, &neq[1], &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    goto L700;
L606:
    xerrwv_("lsoda--  itol (=i1) illegal   ", &c__30, &c__6, &c__0, &c__1, 
	    itol, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L607:
    xerrwv_("lsoda--  iopt (=i1) illegal   ", &c__30, &c__7, &c__0, &c__1, 
	    iopt, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L608:
    xerrwv_("lsoda--  jt (=i1) illegal     ", &c__30, &c__8, &c__0, &c__1, jt,
	     &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L609:
    xerrwv_("lsoda--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2)", &c__50, &
	    c__9, &c__0, &c__2, &ml, &neq[1], &c__0, &c_b42, &c_b42, (ftnlen)
	    50);
    goto L700;
L610:
    xerrwv_("lsoda--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2)", &c__50, &
	    c__10, &c__0, &c__2, &mu, &neq[1], &c__0, &c_b42, &c_b42, (ftnlen)
	    50);
    goto L700;
L611:
    xerrwv_("lsoda--  ixpr (=i1) illegal   ", &c__30, &c__11, &c__0, &c__1, &
	    lsa001_1.ixpr, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L612:
    xerrwv_("lsoda--  mxstep (=i1) .lt. 0  ", &c__30, &c__12, &c__0, &c__1, &
	    ls0001_1.mxstep, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L613:
    xerrwv_("lsoda--  mxhnil (=i1) .lt. 0  ", &c__30, &c__13, &c__0, &c__1, &
	    ls0001_1.mxhnil, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L614:
    xerrwv_("lsoda--  tout (=r1) behind t (=r2)      ", &c__40, &c__14, &c__0,
	     &c__0, &c__0, &c__0, &c__2, tout, t, (ftnlen)40);
    xerrwv_("      integration direction is given by h0 (=r1)  ", &c__50, &
	    c__14, &c__0, &c__0, &c__0, &c__0, &c__1, &h0, &c_b42, (ftnlen)50)
	    ;
    goto L700;
L615:
    xerrwv_("lsoda--  hmax (=r1) .lt. 0.0  ", &c__30, &c__15, &c__0, &c__0, &
	    c__0, &c__0, &c__1, &hmax, &c_b42, (ftnlen)30);
    goto L700;
L616:
    xerrwv_("lsoda--  hmin (=r1) .lt. 0.0  ", &c__30, &c__16, &c__0, &c__0, &
	    c__0, &c__0, &c__1, &ls0001_1.hmin, &c_b42, (ftnlen)30);
    goto L700;
L617:
    xerrwv_("lsoda--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)", &
	    c__60, &c__17, &c__0, &c__2, &lenrw, lrw, &c__0, &c_b42, &c_b42, (
	    ftnlen)60);
    goto L700;
L618:
    xerrwv_("lsoda--  iwork length needed, leniw (=i1), exceeds liw (=i2)", &
	    c__60, &c__18, &c__0, &c__2, &leniw, liw, &c__0, &c_b42, &c_b42, (
	    ftnlen)60);
    goto L700;
L619:
    xerrwv_("lsoda--  rtol(i1) is r1 .lt. 0.0        ", &c__40, &c__19, &c__0,
	     &c__1, &i__, &c__0, &c__1, &rtoli, &c_b42, (ftnlen)40);
    goto L700;
L620:
    xerrwv_("lsoda--  atol(i1) is r1 .lt. 0.0        ", &c__40, &c__20, &c__0,
	     &c__1, &i__, &c__0, &c__1, &atoli, &c_b42, (ftnlen)40);
    goto L700;
L621:
    ewti = rwork[ls0001_1.lewt + i__ - 1];
    xerrwv_("lsoda--  ewt(i1) is r1 .le. 0.0         ", &c__40, &c__21, &c__0,
	     &c__1, &i__, &c__0, &c__1, &ewti, &c_b42, (ftnlen)40);
    goto L700;
L622:
    xerrwv_("lsoda--  tout (=r1) too close to t(=r2) to start integration", &
	    c__60, &c__22, &c__0, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    ftnlen)60);
    goto L700;
L623:
    xerrwv_("lsoda--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ", &
	    c__60, &c__23, &c__0, &c__1, itask, &c__0, &c__2, tout, &tp, (
	    ftnlen)60);
    goto L700;
L624:
    xerrwv_("lsoda--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ", &
	    c__60, &c__24, &c__0, &c__0, &c__0, &c__0, &c__2, &tcrit, &
	    ls0001_1.tn, (ftnlen)60);
    goto L700;
L625:
    xerrwv_("lsoda--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ", &
	    c__60, &c__25, &c__0, &c__0, &c__0, &c__0, &c__2, &tcrit, tout, (
	    ftnlen)60);
    goto L700;
L626:
    xerrwv_("lsoda--  at start of problem, too much accuracy   ", &c__50, &
	    c__26, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)
	    50);
    xerrwv_("      requested for precision of machine..  see tolsf (=r1) ", &
	    c__60, &c__26, &c__0, &c__0, &c__0, &c__0, &c__1, &tolsf, &c_b42, 
	    (ftnlen)60);
    rwork[14] = tolsf;
    goto L700;
L627:
    xerrwv_("lsoda--  trouble from intdy. itask = i1, tout = r1", &c__50, &
	    c__27, &c__0, &c__1, itask, &c__0, &c__1, tout, &c_b42, (ftnlen)
	    50);
    goto L700;
L628:
    xerrwv_("lsoda--  mxordn (=i1) .lt. 0  ", &c__30, &c__28, &c__0, &c__1, &
	    lsa001_1.mxordn, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);
    goto L700;
L629:
    xerrwv_("lsoda--  mxords (=i1) .lt. 0  ", &c__30, &c__29, &c__0, &c__1, &
	    lsa001_1.mxords, &c__0, &c__0, &c_b42, &c_b42, (ftnlen)30);

L700:
    if (ls0001_1.illin == 5) {
	goto L710;
    }
    ++ls0001_1.illin;
    *istate = -3;
    return 0;
L710:
    xerrwv_("lsoda--  repeated occurrences of illegal input    ", &c__50, &
	    c__302, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);

L800:
    xerrwv_("lsoda--  run aborted.. apparent infinite loop     ", &c__50, &
	    c__303, &c__2, &c__0, &c__0, &c__0, &c__0, &c_b42, &c_b42, (
	    ftnlen)50);
    return 0;
/* ----------------------- end of subroutine lsoda ----------------------- */
} /* lsoda_ */

