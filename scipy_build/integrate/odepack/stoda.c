/* stoda.f -- translated by f2c (version 20190311).
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
    doublereal conit, crate, el[13], elco[156]	/* was [13][12] */, hold, 
	    rmax, tesco[36]	/* was [3][12] */, ccmax, el0, h__, hmin, 
	    hmxi, hu, rc, tn, uround;
    integer iownd[14], ialth, ipup, lmax, meo, nqnyh, nslp, icf, ierpj, iersl,
	     jcur, jstart, kflag, l, meth, miter, maxord, maxcor, msbp, mxncf,
	     n, nq, nst, nfe, nje, nqu;
} ls0001_;

#define ls0001_1 ls0001_

struct {
    doublereal rownd2, pdest, pdlast, ratio, cm1[12], cm2[5], pdnorm;
    integer iownd2[3], icount, irflag, jtyp, mused, mxordn, mxords;
} lsa001_;

#define lsa001_1 lsa001_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int stoda_(integer *neq, doublereal *y, doublereal *yh, 
	integer *nyh, doublereal *yh1, doublereal *ewt, doublereal *savf, 
	doublereal *acor, doublereal *wm, integer *iwm, S_fp f, U_fp jac, 
	S_fp pjac, S_fp slvs)
{
    /* Initialized data */

    static doublereal sm1[12] = { .5,.575,.55,.45,.35,.25,.2,.15,.1,.075,.05,
	    .025 };

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, m;
    static doublereal r__;
    static integer i1, jb;
    static doublereal rh, rm, dm1, dm2;
    static integer lm1, lm2;
    static doublereal rh1, rh2, del, ddn;
    static integer ncf;
    static doublereal pdh, dsm, dup, exm1, exm2;
    static integer nqm1, nqm2;
    static doublereal dcon, delp;
    static integer lm1p1, lm2p1;
    static doublereal exdn, rhdn;
    static integer iret, isav[50];
    static doublereal told, rhsm;
    static integer newq;
    static doublereal rsav[240], exsm, rhup, rate, exup, rh1it;
    extern /* Subroutine */ int cfode_(integer *, doublereal *, doublereal *);
    static doublereal alpha;
    static integer iredo;
    extern /* Subroutine */ int srcma_(doublereal *, integer *, integer *);
    static doublereal pnorm;
    extern doublereal vmnorm_(integer *, doublereal *, doublereal *);

/* lll. optimize */
    /* Parameter adjustments */
    --neq;
    --y;
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --yh1;
    --ewt;
    --savf;
    --acor;
    --wm;
    --iwm;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* stoda performs one step of the integration of an initial value */
/* problem for a system of ordinary differential equations. */
/* note.. stoda is independent of the value of the iteration method */
/* indicator miter, when this is .ne. 0, and hence is independent */
/* of the type of chord method used, or the jacobian structure. */
/* communication with stoda is done with the following variables.. */

/* y      = an array of length .ge. n used as the y argument in */
/*          all calls to f and jac. */
/* neq    = integer array containing problem size in neq(1), and */
/*          passed as the neq argument in all calls to f and jac. */
/* yh     = an nyh by lmax array containing the dependent variables */
/*          and their approximate scaled derivatives, where */
/*          lmax = maxord + 1.  yh(i,j+1) contains the approximate */
/*          j-th derivative of y(i), scaled by h**j/factorial(j) */
/*          (j = 0,1,...,nq).  on entry for the first step, the first */
/*          two columns of yh must be set from the initial values. */
/* nyh    = a constant integer .ge. n, the first dimension of yh. */
/* yh1    = a one-dimensional array occupying the same space as yh. */
/* ewt    = an array of length n containing multiplicative weights */
/*          for local error measurements.  local errors in y(i) are */
/*          compared to 1.0/ewt(i) in various error tests. */
/* savf   = an array of working storage, of length n. */
/* acor   = a work array of length n, used for the accumulated */
/*          corrections.  on a successful return, acor(i) contains */
/*          the estimated one-step local error in y(i). */
/* wm,iwm = real and integer work arrays associated with matrix */
/*          operations in chord iteration (miter .ne. 0). */
/* pjac   = name of routine to evaluate and preprocess jacobian matrix */
/*          and p = i - h*el0*jac, if a chord method is being used. */
/*          it also returns an estimate of norm(jac) in pdnorm. */
/* slvs   = name of routine to solve linear system in chord iteration. */
/* ccmax  = maximum relative change in h*el0 before pjac is called. */
/* h      = the step size to be attempted on the next step. */
/*          h is altered by the error control algorithm during the */
/*          problem.  h can be either positive or negative, but its */
/*          sign must remain constant throughout the problem. */
/* hmin   = the minimum absolute value of the step size h to be used. */
/* hmxi   = inverse of the maximum absolute value of h to be used. */
/*          hmxi = 0.0 is allowed and corresponds to an infinite hmax. */
/*          hmin and hmxi may be changed at any time, but will not */
/*          take effect until the next change of h is considered. */
/* tn     = the independent variable. tn is updated on each step taken. */
/* jstart = an integer used for input only, with the following */
/*          values and meanings.. */
/*               0  perform the first step. */
/*           .gt.0  take a new step continuing from the last. */
/*              -1  take the next step with a new value of h, */
/*                    n, meth, miter, and/or matrix parameters. */
/*              -2  take the next step with a new value of h, */
/*                    but with other inputs unchanged. */
/*          on return, jstart is set to 1 to facilitate continuation. */
/* kflag  = a completion code with the following meanings.. */
/*               0  the step was successful. */
/*              -1  the requested error could not be achieved. */
/*              -2  corrector convergence could not be achieved. */
/*              -3  fatal error in pjac or slvs. */
/*          a return with kflag = -1 or -2 means either */
/*          abs(h) = hmin or 10 consecutive failures occurred. */
/*          on a return with kflag negative, the values of tn and */
/*          the yh array are as of the beginning of the last */
/*          step, and h is the last step size attempted. */
/* maxord = the maximum order of integration method to be allowed. */
/* maxcor = the maximum number of corrector iterations allowed. */
/* msbp   = maximum number of steps between pjac calls (miter .gt. 0). */
/* mxncf  = maximum number of convergence failures allowed. */
/* meth   = current method. */
/*          meth = 1 means adams method (nonstiff) */
/*          meth = 2 means bdf method (stiff) */
/*          meth may be reset by stoda. */
/* miter  = corrector iteration method. */
/*          miter = 0 means functional iteration. */
/*          miter = jt .gt. 0 means a chord iteration corresponding */
/*          to jacobian type jt.  (the lsoda argument jt is */
/*          communicated here as jtyp, but is not used in stoda */
/*          except to load miter following a method switch.) */
/*          miter may be reset by stoda. */
/* n      = the number of first-order differential equations. */
/* ----------------------------------------------------------------------- */
    ls0001_1.kflag = 0;
    told = ls0001_1.tn;
    ncf = 0;
    ls0001_1.ierpj = 0;
    ls0001_1.iersl = 0;
    ls0001_1.jcur = 0;
    ls0001_1.icf = 0;
    delp = 0.;
    if (ls0001_1.jstart > 0) {
	goto L200;
    }
    if (ls0001_1.jstart == -1) {
	goto L100;
    }
    if (ls0001_1.jstart == -2) {
	goto L160;
    }
/* ----------------------------------------------------------------------- */
/* on the first call, the order is set to 1, and other variables are */
/* initialized.  rmax is the maximum ratio by which h can be increased */
/* in a single step.  it is initially 1.e4 to compensate for the small */
/* initial h, but then is normally equal to 10.  if a failure */
/* occurs (in corrector convergence or error test), rmax is set at 2 */
/* for the next increase. */
/* cfode is called to get the needed coefficients for both methods. */
/* ----------------------------------------------------------------------- */
    ls0001_1.lmax = ls0001_1.maxord + 1;
    ls0001_1.nq = 1;
    ls0001_1.l = 2;
    ls0001_1.ialth = 2;
    ls0001_1.rmax = 1e4;
    ls0001_1.rc = 0.;
    ls0001_1.el0 = 1.;
    ls0001_1.crate = .7;
    ls0001_1.hold = ls0001_1.h__;
    ls0001_1.nslp = 0;
    ls0001_1.ipup = ls0001_1.miter;
    iret = 3;
/* initialize switching parameters.  meth = 1 is assumed initially. ----- */
    lsa001_1.icount = 20;
    lsa001_1.irflag = 0;
    lsa001_1.pdest = 0.;
    lsa001_1.pdlast = 0.;
    lsa001_1.ratio = 5.;
    cfode_(&c__2, ls0001_1.elco, ls0001_1.tesco);
    for (i__ = 1; i__ <= 5; ++i__) {
/* L10: */
	lsa001_1.cm2[i__ - 1] = ls0001_1.tesco[i__ * 3 - 2] * ls0001_1.elco[
		i__ + 1 + i__ * 13 - 14];
    }
    cfode_(&c__1, ls0001_1.elco, ls0001_1.tesco);
    for (i__ = 1; i__ <= 12; ++i__) {
/* L20: */
	lsa001_1.cm1[i__ - 1] = ls0001_1.tesco[i__ * 3 - 2] * ls0001_1.elco[
		i__ + 1 + i__ * 13 - 14];
    }
    goto L150;
/* ----------------------------------------------------------------------- */
/* the following block handles preliminaries needed when jstart = -1. */
/* ipup is set to miter to force a matrix update. */
/* if an order increase is about to be considered (ialth = 1), */
/* ialth is reset to 2 to postpone consideration one more step. */
/* if the caller has changed meth, cfode is called to reset */
/* the coefficients of the method. */
/* if h is to be changed, yh must be rescaled. */
/* if h or meth is being changed, ialth is reset to l = nq + 1 */
/* to prevent further changes in h for that many steps. */
/* ----------------------------------------------------------------------- */
L100:
    ls0001_1.ipup = ls0001_1.miter;
    ls0001_1.lmax = ls0001_1.maxord + 1;
    if (ls0001_1.ialth == 1) {
	ls0001_1.ialth = 2;
    }
    if (ls0001_1.meth == lsa001_1.mused) {
	goto L160;
    }
    cfode_(&ls0001_1.meth, ls0001_1.elco, ls0001_1.tesco);
    ls0001_1.ialth = ls0001_1.l;
    iret = 1;
/* ----------------------------------------------------------------------- */
/* the el vector and related constants are reset */
/* whenever the order nq is changed, or at the start of the problem. */
/* ----------------------------------------------------------------------- */
L150:
    i__1 = ls0001_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L155: */
	ls0001_1.el[i__ - 1] = ls0001_1.elco[i__ + ls0001_1.nq * 13 - 14];
    }
    ls0001_1.nqnyh = ls0001_1.nq * *nyh;
    ls0001_1.rc = ls0001_1.rc * ls0001_1.el[0] / ls0001_1.el0;
    ls0001_1.el0 = ls0001_1.el[0];
    ls0001_1.conit = .5 / (doublereal) (ls0001_1.nq + 2);
    switch (iret) {
	case 1:  goto L160;
	case 2:  goto L170;
	case 3:  goto L200;
    }
/* ----------------------------------------------------------------------- */
/* if h is being changed, the h ratio rh is checked against */
/* rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to */
/* l = nq + 1 to prevent a change of h for that many steps, unless */
/* forced by a convergence or error test failure. */
/* ----------------------------------------------------------------------- */
L160:
    if (ls0001_1.h__ == ls0001_1.hold) {
	goto L200;
    }
    rh = ls0001_1.h__ / ls0001_1.hold;
    ls0001_1.h__ = ls0001_1.hold;
    iredo = 3;
    goto L175;
L170:
/* Computing MAX */
    d__1 = rh, d__2 = ls0001_1.hmin / abs(ls0001_1.h__);
    rh = max(d__1,d__2);
L175:
    rh = min(rh,ls0001_1.rmax);
/* Computing MAX */
    d__1 = 1., d__2 = abs(ls0001_1.h__) * ls0001_1.hmxi * rh;
    rh /= max(d__1,d__2);
/* ----------------------------------------------------------------------- */
/* if meth = 1, also restrict the new step size by the stability region. */
/* if this reduces h, set irflag to 1 so that if there are roundoff */
/* problems later, we can assume that is the cause of the trouble. */
/* ----------------------------------------------------------------------- */
    if (ls0001_1.meth == 2) {
	goto L178;
    }
    lsa001_1.irflag = 0;
/* Computing MAX */
    d__1 = abs(ls0001_1.h__) * lsa001_1.pdlast;
    pdh = max(d__1,1e-6);
    if (rh * pdh * 1.00001 < sm1[ls0001_1.nq - 1]) {
	goto L178;
    }
    rh = sm1[ls0001_1.nq - 1] / pdh;
    lsa001_1.irflag = 1;
L178:
    r__ = 1.;
    i__1 = ls0001_1.l;
    for (j = 2; j <= i__1; ++j) {
	r__ *= rh;
	i__2 = ls0001_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L180: */
	    yh[i__ + j * yh_dim1] *= r__;
	}
    }
    ls0001_1.h__ *= rh;
    ls0001_1.rc *= rh;
    ls0001_1.ialth = ls0001_1.l;
    if (iredo == 0) {
	goto L690;
    }
/* ----------------------------------------------------------------------- */
/* this section computes the predicted values by effectively */
/* multiplying the yh array by the pascal triangle matrix. */
/* rc is the ratio of new to old values of the coefficient  h*el(1). */
/* when rc differs from 1 by more than ccmax, ipup is set to miter */
/* to force pjac to be called, if a jacobian is involved. */
/* in any case, pjac is called at least every msbp steps. */
/* ----------------------------------------------------------------------- */
L200:
    if ((d__1 = ls0001_1.rc - 1., abs(d__1)) > ls0001_1.ccmax) {
	ls0001_1.ipup = ls0001_1.miter;
    }
    if (ls0001_1.nst >= ls0001_1.nslp + ls0001_1.msbp) {
	ls0001_1.ipup = ls0001_1.miter;
    }
    ls0001_1.tn += ls0001_1.h__;
    i1 = ls0001_1.nqnyh + 1;
    i__2 = ls0001_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *nyh;
/* dir$ ivdep */
	i__1 = ls0001_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L210: */
	    yh1[i__] += yh1[i__ + *nyh];
	}
/* L215: */
    }
    pnorm = vmnorm_(&ls0001_1.n, &yh1[1], &ewt[1]);
/* ----------------------------------------------------------------------- */
/* up to maxcor corrector iterations are taken.  a convergence test is */
/* made on the r.m.s. norm of each correction, weighted by the error */
/* weight vector ewt.  the sum of the corrections is accumulated in the */
/* vector acor(i).  the yh array is not altered in the corrector loop. */
/* ----------------------------------------------------------------------- */
L220:
    m = 0;
    rate = 0.;
    del = 0.;
    i__2 = ls0001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L230: */
	y[i__] = yh[i__ + yh_dim1];
    }
    srcma_(rsav, isav, &c__1);
    (*f)(&neq[1], &ls0001_1.tn, &y[1], &savf[1]);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    srcma_(rsav, isav, &c__2);
    ++ls0001_1.nfe;
    if (ls0001_1.ipup <= 0) {
	goto L250;
    }
/* ----------------------------------------------------------------------- */
/* if indicated, the matrix p = i - h*el(1)*j is reevaluated and */
/* preprocessed before starting the corrector iteration.  ipup is set */
/* to 0 as an indicator that this has been done. */
/* ----------------------------------------------------------------------- */
    (*pjac)(&neq[1], &y[1], &yh[yh_offset], nyh, &ewt[1], &acor[1], &savf[1], 
	    &wm[1], &iwm[1], (S_fp)f, (U_fp)jac);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    ls0001_1.ipup = 0;
    ls0001_1.rc = 1.;
    ls0001_1.nslp = ls0001_1.nst;
    ls0001_1.crate = .7;
    if (ls0001_1.ierpj != 0) {
	goto L430;
    }
L250:
    i__2 = ls0001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L260: */
	acor[i__] = 0.;
    }
L270:
    if (ls0001_1.miter != 0) {
	goto L350;
    }
/* ----------------------------------------------------------------------- */
/* in the case of functional iteration, update y directly from */
/* the result of the last function evaluation. */
/* ----------------------------------------------------------------------- */
    i__2 = ls0001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	savf[i__] = ls0001_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
/* L290: */
	y[i__] = savf[i__] - acor[i__];
    }
    del = vmnorm_(&ls0001_1.n, &y[1], &ewt[1]);
    i__2 = ls0001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__] = yh[i__ + yh_dim1] + ls0001_1.el[0] * savf[i__];
/* L300: */
	acor[i__] = savf[i__];
    }
    goto L400;
/* ----------------------------------------------------------------------- */
/* in the case of the chord method, compute the corrector error, */
/* and solve the linear system with that as right-hand side and */
/* p as coefficient matrix. */
/* ----------------------------------------------------------------------- */
L350:
    i__2 = ls0001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L360: */
	y[i__] = ls0001_1.h__ * savf[i__] - (yh[i__ + (yh_dim1 << 1)] + acor[
		i__]);
    }
    (*slvs)(&wm[1], &iwm[1], &y[1], &savf[1]);
    if (ls0001_1.iersl < 0) {
	goto L430;
    }
    if (ls0001_1.iersl > 0) {
	goto L410;
    }
    del = vmnorm_(&ls0001_1.n, &y[1], &ewt[1]);
    i__2 = ls0001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	acor[i__] += y[i__];
/* L380: */
	y[i__] = yh[i__ + yh_dim1] + ls0001_1.el[0] * acor[i__];
    }
/* ----------------------------------------------------------------------- */
/* test for convergence.  if m.gt.0, an estimate of the convergence */
/* rate constant is stored in crate, and this is used in the test. */

/* we first check for a change of iterates that is the size of */
/* roundoff error.  if this occurs, the iteration has converged, and a */
/* new rate estimate is not formed. */
/* in all other cases, force at least two iterations to estimate a */
/* local lipschitz constant estimate for adams methods. */
/* on convergence, form pdest = local maximum lipschitz constant */
/* estimate.  pdlast is the most recent nonzero estimate. */
/* ----------------------------------------------------------------------- */
L400:
    if (del <= pnorm * 100. * ls0001_1.uround) {
	goto L450;
    }
    if (m == 0 && ls0001_1.meth == 1) {
	goto L405;
    }
    if (m == 0) {
	goto L402;
    }
    rm = 1024.;
    if (del <= delp * 1024.) {
	rm = del / delp;
    }
    rate = max(rate,rm);
/* Computing MAX */
    d__1 = ls0001_1.crate * .2;
    ls0001_1.crate = max(d__1,rm);
L402:
/* Computing MIN */
    d__1 = 1., d__2 = ls0001_1.crate * 1.5;
    dcon = del * min(d__1,d__2) / (ls0001_1.tesco[ls0001_1.nq * 3 - 2] * 
	    ls0001_1.conit);
    if (dcon > 1.) {
	goto L405;
    }
/* Computing MAX */
    d__2 = lsa001_1.pdest, d__3 = rate / (d__1 = ls0001_1.h__ * ls0001_1.el[0]
	    , abs(d__1));
    lsa001_1.pdest = max(d__2,d__3);
    if (lsa001_1.pdest != 0.) {
	lsa001_1.pdlast = lsa001_1.pdest;
    }
    goto L450;
L405:
    ++m;
    if (m == ls0001_1.maxcor) {
	goto L410;
    }
    if (m >= 2 && del > delp * 2.) {
	goto L410;
    }
    delp = del;
    srcma_(rsav, isav, &c__1);
    (*f)(&neq[1], &ls0001_1.tn, &y[1], &savf[1]);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    srcma_(rsav, isav, &c__2);
    ++ls0001_1.nfe;
    goto L270;
/* ----------------------------------------------------------------------- */
/* the corrector iteration failed to converge. */
/* if miter .ne. 0 and the jacobian is out of date, pjac is called for */
/* the next try.  otherwise the yh array is retracted to its values */
/* before prediction, and h is reduced, if possible.  if h cannot be */
/* reduced or mxncf failures have occurred, exit with kflag = -2. */
/* ----------------------------------------------------------------------- */
L410:
    if (ls0001_1.miter == 0 || ls0001_1.jcur == 1) {
	goto L430;
    }
    ls0001_1.icf = 1;
    ls0001_1.ipup = ls0001_1.miter;
    goto L220;
L430:
    ls0001_1.icf = 2;
    ++ncf;
    ls0001_1.rmax = 2.;
    ls0001_1.tn = told;
    i1 = ls0001_1.nqnyh + 1;
    i__2 = ls0001_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *nyh;
/* dir$ ivdep */
	i__1 = ls0001_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L440: */
	    yh1[i__] -= yh1[i__ + *nyh];
	}
/* L445: */
    }
    if (ls0001_1.ierpj < 0 || ls0001_1.iersl < 0) {
	goto L680;
    }
    if (abs(ls0001_1.h__) <= ls0001_1.hmin * 1.00001) {
	goto L670;
    }
    if (ncf == ls0001_1.mxncf) {
	goto L670;
    }
    rh = .25;
    ls0001_1.ipup = ls0001_1.miter;
    iredo = 1;
    goto L170;
/* ----------------------------------------------------------------------- */
/* the corrector has converged.  jcur is set to 0 */
/* to signal that the jacobian involved may need updating later. */
/* the local error test is made and control passes to statement 500 */
/* if it fails. */
/* ----------------------------------------------------------------------- */
L450:
    ls0001_1.jcur = 0;
    if (m == 0) {
	dsm = del / ls0001_1.tesco[ls0001_1.nq * 3 - 2];
    }
    if (m > 0) {
	dsm = vmnorm_(&ls0001_1.n, &acor[1], &ewt[1]) / ls0001_1.tesco[
		ls0001_1.nq * 3 - 2];
    }
    if (dsm > 1.) {
	goto L500;
    }
/* ----------------------------------------------------------------------- */
/* after a successful step, update the yh array. */
/* decrease icount by 1, and if it is -1, consider switching methods. */
/* if a method switch is made, reset various parameters, */
/* rescale the yh array, and exit.  if there is no switch, */
/* consider changing h if ialth = 1.  otherwise decrease ialth by 1. */
/* if ialth is then 1 and nq .lt. maxord, then acor is saved for */
/* use in a possible order increase on the next step. */
/* if a change in h is considered, an increase or decrease in order */
/* by one is considered also.  a change in h is made only if it is by a */
/* factor of at least 1.1.  if not, ialth is set to 3 to prevent */
/* testing for that many steps. */
/* ----------------------------------------------------------------------- */
    ls0001_1.kflag = 0;
    iredo = 0;
    ++ls0001_1.nst;
    ls0001_1.hu = ls0001_1.h__;
    ls0001_1.nqu = ls0001_1.nq;
    lsa001_1.mused = ls0001_1.meth;
    i__2 = ls0001_1.l;
    for (j = 1; j <= i__2; ++j) {
	i__1 = ls0001_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L460: */
	    yh[i__ + j * yh_dim1] += ls0001_1.el[j - 1] * acor[i__];
	}
    }
    --lsa001_1.icount;
    if (lsa001_1.icount >= 0) {
	goto L488;
    }
    if (ls0001_1.meth == 2) {
	goto L480;
    }
/* ----------------------------------------------------------------------- */
/* we are currently using an adams method.  consider switching to bdf. */
/* if the current order is greater than 5, assume the problem is */
/* not stiff, and skip this section. */
/* if the lipschitz constant and error estimate are not polluted */
/* by roundoff, go to 470 and perform the usual test. */
/* otherwise, switch to the bdf methods if the last step was */
/* restricted to insure stability (irflag = 1), and stay with adams */
/* method if not.  when switching to bdf with polluted error estimates, */
/* in the absence of other information, double the step size. */

/* when the estimates are ok, we make the usual test by computing */
/* the step size we could have (ideally) used on this step, */
/* with the current (adams) method, and also that for the bdf. */
/* if nq .gt. mxords, we consider changing to order mxords on switching. */
/* compare the two step sizes to decide whether to switch. */
/* the step size advantage must be at least ratio = 5 to switch. */
/* ----------------------------------------------------------------------- */
    if (ls0001_1.nq > 5) {
	goto L488;
    }
    if (dsm > pnorm * 100. * ls0001_1.uround && lsa001_1.pdest != 0.) {
	goto L470;
    }
    if (lsa001_1.irflag == 0) {
	goto L488;
    }
    rh2 = 2.;
    nqm2 = min(ls0001_1.nq,lsa001_1.mxords);
    goto L478;
L470:
    exsm = 1. / (doublereal) ls0001_1.l;
    rh1 = 1. / (pow_dd(&dsm, &exsm) * 1.2 + 1.2e-6);
    rh1it = rh1 * 2.;
    pdh = lsa001_1.pdlast * abs(ls0001_1.h__);
    if (pdh * rh1 > 1e-5) {
	rh1it = sm1[ls0001_1.nq - 1] / pdh;
    }
    rh1 = min(rh1,rh1it);
    if (ls0001_1.nq <= lsa001_1.mxords) {
	goto L474;
    }
    nqm2 = lsa001_1.mxords;
    lm2 = lsa001_1.mxords + 1;
    exm2 = 1. / (doublereal) lm2;
    lm2p1 = lm2 + 1;
    dm2 = vmnorm_(&ls0001_1.n, &yh[lm2p1 * yh_dim1 + 1], &ewt[1]) / 
	    lsa001_1.cm2[lsa001_1.mxords - 1];
    rh2 = 1. / (pow_dd(&dm2, &exm2) * 1.2 + 1.2e-6);
    goto L476;
L474:
    dm2 = dsm * (lsa001_1.cm1[ls0001_1.nq - 1] / lsa001_1.cm2[ls0001_1.nq - 1]
	    );
    rh2 = 1. / (pow_dd(&dm2, &exsm) * 1.2 + 1.2e-6);
    nqm2 = ls0001_1.nq;
L476:
    if (rh2 < lsa001_1.ratio * rh1) {
	goto L488;
    }
/* the switch test passed.  reset relevant quantities for bdf. ---------- */
L478:
    rh = rh2;
    lsa001_1.icount = 20;
    ls0001_1.meth = 2;
    ls0001_1.miter = lsa001_1.jtyp;
    lsa001_1.pdlast = 0.;
    ls0001_1.nq = nqm2;
    ls0001_1.l = ls0001_1.nq + 1;
    goto L170;
/* ----------------------------------------------------------------------- */
/* we are currently using a bdf method.  consider switching to adams. */
/* compute the step size we could have (ideally) used on this step, */
/* with the current (bdf) method, and also that for the adams. */
/* if nq .gt. mxordn, we consider changing to order mxordn on switching. */
/* compare the two step sizes to decide whether to switch. */
/* the step size advantage must be at least 5/ratio = 1 to switch. */
/* if the step size for adams would be so small as to cause */
/* roundoff pollution, we stay with bdf. */
/* ----------------------------------------------------------------------- */
L480:
    exsm = 1. / (doublereal) ls0001_1.l;
    if (lsa001_1.mxordn >= ls0001_1.nq) {
	goto L484;
    }
    nqm1 = lsa001_1.mxordn;
    lm1 = lsa001_1.mxordn + 1;
    exm1 = 1. / (doublereal) lm1;
    lm1p1 = lm1 + 1;
    dm1 = vmnorm_(&ls0001_1.n, &yh[lm1p1 * yh_dim1 + 1], &ewt[1]) / 
	    lsa001_1.cm1[lsa001_1.mxordn - 1];
    rh1 = 1. / (pow_dd(&dm1, &exm1) * 1.2 + 1.2e-6);
    goto L486;
L484:
    dm1 = dsm * (lsa001_1.cm2[ls0001_1.nq - 1] / lsa001_1.cm1[ls0001_1.nq - 1]
	    );
    rh1 = 1. / (pow_dd(&dm1, &exsm) * 1.2 + 1.2e-6);
    nqm1 = ls0001_1.nq;
    exm1 = exsm;
L486:
    rh1it = rh1 * 2.;
    pdh = lsa001_1.pdnorm * abs(ls0001_1.h__);
    if (pdh * rh1 > 1e-5) {
	rh1it = sm1[nqm1 - 1] / pdh;
    }
    rh1 = min(rh1,rh1it);
    rh2 = 1. / (pow_dd(&dsm, &exsm) * 1.2 + 1.2e-6);
    if (rh1 * lsa001_1.ratio < rh2 * 5.) {
	goto L488;
    }
    alpha = max(.001,rh1);
    dm1 = pow_dd(&alpha, &exm1) * dm1;
    if (dm1 <= ls0001_1.uround * 1e3 * pnorm) {
	goto L488;
    }
/* the switch test passed.  reset relevant quantities for adams. -------- */
    rh = rh1;
    lsa001_1.icount = 20;
    ls0001_1.meth = 1;
    ls0001_1.miter = 0;
    lsa001_1.pdlast = 0.;
    ls0001_1.nq = nqm1;
    ls0001_1.l = ls0001_1.nq + 1;
    goto L170;

/* no method switch is being made.  do the usual step/order selection. -- */
L488:
    --ls0001_1.ialth;
    if (ls0001_1.ialth == 0) {
	goto L520;
    }
    if (ls0001_1.ialth > 1) {
	goto L700;
    }
    if (ls0001_1.l == ls0001_1.lmax) {
	goto L700;
    }
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L490: */
	yh[i__ + ls0001_1.lmax * yh_dim1] = acor[i__];
    }
    goto L700;
/* ----------------------------------------------------------------------- */
/* the error test failed.  kflag keeps track of multiple failures. */
/* restore tn and the yh array to their previous values, and prepare */
/* to try the step again.  compute the optimum step size for this or */
/* one lower order.  after 2 or more failures, h is forced to decrease */
/* by a factor of 0.2 or less. */
/* ----------------------------------------------------------------------- */
L500:
    --ls0001_1.kflag;
    ls0001_1.tn = told;
    i1 = ls0001_1.nqnyh + 1;
    i__1 = ls0001_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *nyh;
/* dir$ ivdep */
	i__2 = ls0001_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L510: */
	    yh1[i__] -= yh1[i__ + *nyh];
	}
/* L515: */
    }
    ls0001_1.rmax = 2.;
    if (abs(ls0001_1.h__) <= ls0001_1.hmin * 1.00001) {
	goto L660;
    }
    if (ls0001_1.kflag <= -3) {
	goto L640;
    }
    iredo = 2;
    rhup = 0.;
    goto L540;
/* ----------------------------------------------------------------------- */
/* regardless of the success or failure of the step, factors */
/* rhdn, rhsm, and rhup are computed, by which h could be multiplied */
/* at order nq - 1, order nq, or order nq + 1, respectively. */
/* in the case of failure, rhup = 0.0 to avoid an order increase. */
/* the largest of these is determined and the new order chosen */
/* accordingly.  if the order is to be increased, we compute one */
/* additional scaled derivative. */
/* ----------------------------------------------------------------------- */
L520:
    rhup = 0.;
    if (ls0001_1.l == ls0001_1.lmax) {
	goto L540;
    }
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L530: */
	savf[i__] = acor[i__] - yh[i__ + ls0001_1.lmax * yh_dim1];
    }
    dup = vmnorm_(&ls0001_1.n, &savf[1], &ewt[1]) / ls0001_1.tesco[
	    ls0001_1.nq * 3 - 1];
    exup = 1. / (doublereal) (ls0001_1.l + 1);
    rhup = 1. / (pow_dd(&dup, &exup) * 1.4 + 1.4e-6);
L540:
    exsm = 1. / (doublereal) ls0001_1.l;
    rhsm = 1. / (pow_dd(&dsm, &exsm) * 1.2 + 1.2e-6);
    rhdn = 0.;
    if (ls0001_1.nq == 1) {
	goto L550;
    }
    ddn = vmnorm_(&ls0001_1.n, &yh[ls0001_1.l * yh_dim1 + 1], &ewt[1]) / 
	    ls0001_1.tesco[ls0001_1.nq * 3 - 3];
    exdn = 1. / (doublereal) ls0001_1.nq;
    rhdn = 1. / (pow_dd(&ddn, &exdn) * 1.3 + 1.3e-6);
/* if meth = 1, limit rh according to the stability region also. -------- */
L550:
    if (ls0001_1.meth == 2) {
	goto L560;
    }
/* Computing MAX */
    d__1 = abs(ls0001_1.h__) * lsa001_1.pdlast;
    pdh = max(d__1,1e-6);
    if (ls0001_1.l < ls0001_1.lmax) {
/* Computing MIN */
	d__1 = rhup, d__2 = sm1[ls0001_1.l - 1] / pdh;
	rhup = min(d__1,d__2);
    }
/* Computing MIN */
    d__1 = rhsm, d__2 = sm1[ls0001_1.nq - 1] / pdh;
    rhsm = min(d__1,d__2);
    if (ls0001_1.nq > 1) {
/* Computing MIN */
	d__1 = rhdn, d__2 = sm1[ls0001_1.nq - 2] / pdh;
	rhdn = min(d__1,d__2);
    }
    lsa001_1.pdest = 0.;
L560:
    if (rhsm >= rhup) {
	goto L570;
    }
    if (rhup > rhdn) {
	goto L590;
    }
    goto L580;
L570:
    if (rhsm < rhdn) {
	goto L580;
    }
    newq = ls0001_1.nq;
    rh = rhsm;
    goto L620;
L580:
    newq = ls0001_1.nq - 1;
    rh = rhdn;
    if (ls0001_1.kflag < 0 && rh > 1.) {
	rh = 1.;
    }
    goto L620;
L590:
    newq = ls0001_1.l;
    rh = rhup;
    if (rh < 1.1) {
	goto L610;
    }
    r__ = ls0001_1.el[ls0001_1.l - 1] / (doublereal) ls0001_1.l;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L600: */
	yh[i__ + (newq + 1) * yh_dim1] = acor[i__] * r__;
    }
    goto L630;
L610:
    ls0001_1.ialth = 3;
    goto L700;
/* if meth = 1 and h is restricted by stability, bypass 10 percent test. */
L620:
    if (ls0001_1.meth == 2) {
	goto L622;
    }
    if (rh * pdh * 1.00001 >= sm1[newq - 1]) {
	goto L625;
    }
L622:
    if (ls0001_1.kflag == 0 && rh < 1.1) {
	goto L610;
    }
L625:
    if (ls0001_1.kflag <= -2) {
	rh = min(rh,.2);
    }
/* ----------------------------------------------------------------------- */
/* if there is a change of order, reset nq, l, and the coefficients. */
/* in any case h is reset according to rh and the yh array is rescaled. */
/* then exit from 690 if the step was ok, or redo the step otherwise. */
/* ----------------------------------------------------------------------- */
    if (newq == ls0001_1.nq) {
	goto L170;
    }
L630:
    ls0001_1.nq = newq;
    ls0001_1.l = ls0001_1.nq + 1;
    iret = 2;
    goto L150;
/* ----------------------------------------------------------------------- */
/* control reaches this section if 3 or more failures have occurred. */
/* if 10 failures have occurred, exit with kflag = -1. */
/* it is assumed that the derivatives that have accumulated in the */
/* yh array have errors of the wrong order.  hence the first */
/* derivative is recomputed, and the order is set to 1.  then */
/* h is reduced by a factor of 10, and the step is retried, */
/* until it succeeds or h reaches hmin. */
/* ----------------------------------------------------------------------- */
L640:
    if (ls0001_1.kflag == -10) {
	goto L660;
    }
    rh = .1;
/* Computing MAX */
    d__1 = ls0001_1.hmin / abs(ls0001_1.h__);
    rh = max(d__1,rh);
    ls0001_1.h__ *= rh;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L645: */
	y[i__] = yh[i__ + yh_dim1];
    }
    srcma_(rsav, isav, &c__1);
    (*f)(&neq[1], &ls0001_1.tn, &y[1], &savf[1]);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    srcma_(rsav, isav, &c__2);
    ++ls0001_1.nfe;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L650: */
	yh[i__ + (yh_dim1 << 1)] = ls0001_1.h__ * savf[i__];
    }
    ls0001_1.ipup = ls0001_1.miter;
    ls0001_1.ialth = 5;
    if (ls0001_1.nq == 1) {
	goto L200;
    }
    ls0001_1.nq = 1;
    ls0001_1.l = 2;
    iret = 3;
    goto L150;
/* ----------------------------------------------------------------------- */
/* all returns are made through this section.  h is saved in hold */
/* to allow the caller to change h on the next step. */
/* ----------------------------------------------------------------------- */
L660:
    ls0001_1.kflag = -1;
    goto L720;
L670:
    ls0001_1.kflag = -2;
    goto L720;
L680:
    ls0001_1.kflag = -3;
    goto L720;
L690:
    ls0001_1.rmax = 10.;
L700:
    r__ = 1. / ls0001_1.tesco[ls0001_1.nqu * 3 - 2];
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L710: */
	acor[i__] *= r__;
    }
L720:
    ls0001_1.hold = ls0001_1.h__;
    ls0001_1.jstart = 1;
    return 0;
/* ----------------------- end of subroutine stoda ----------------------- */
} /* stoda_ */

