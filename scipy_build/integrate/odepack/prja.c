/* prja.f -- translated by f2c (version 20190311).
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
    integer iownd[14], iowns[6], icf, ierpj, iersl, jcur, jstart, kflag, l, 
	    meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, 
	    nqu;
} ls0001_;

#define ls0001_1 ls0001_

struct {
    doublereal rownd2, rowns2[20], pdnorm;
    integer iownd2[3], iowns2[2], jtyp, mused, mxordn, mxords;
} lsa001_;

#define lsa001_1 lsa001_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/* Subroutine */ int prja_(integer *neq, doublereal *y, doublereal *yh, 
	integer *nyh, doublereal *ewt, doublereal *ftem, doublereal *savf, 
	doublereal *wm, integer *iwm, S_fp f, S_fp jac)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal r__;
    static integer i1, i2, j1;
    static doublereal r0;
    static integer ii, jj, ml, mu;
    static doublereal yi, yj, hl0;
    static integer ml3, np1;
    static doublereal fac;
    static integer mba, ier;
    static doublereal con, yjj;
    static integer meb1, lenp, isav[50];
    static doublereal rsav[240], srur;
    static integer mband;
    extern /* Subroutine */ int srcma_(doublereal *, integer *, integer *);
    extern doublereal bnorm_(integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *), fnorm_(integer *, doublereal *, 
	    doublereal *);
    static integer meband;
    extern /* Subroutine */ int dgbtrf_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *);
    extern doublereal vmnorm_(integer *, doublereal *, doublereal *);

/* lll. optimize */
/* ----------------------------------------------------------------------- */
/* prja is called by stoda to compute and process the matrix */
/* p = i - h*el(1)*j , where j is an approximation to the jacobian. */
/* here j is computed by the user-supplied routine jac if */
/* miter = 1 or 4 or by finite differencing if miter = 2 or 5. */
/* j, scaled by -h*el(1), is stored in wm.  then the norm of j (the */
/* matrix norm consistent with the weighted max-norm on vectors given */
/* by vmnorm) is computed, and j is overwritten by p.  p is then */
/* subjected to lu decomposition in preparation for later solution */
/* of linear systems with p as coefficient matrix. this is done */
/* by dgetrf if miter = 1 or 2, and by dgbtrf if miter = 4 or 5. */

/* in addition to variables described previously, communication */
/* with prja uses the following.. */
/* y     = array containing predicted values on entry. */
/* ftem  = work array of length n (acor in stoda). */
/* savf  = array containing f evaluated at predicted y. */
/* wm    = real work space for matrices.  on output it contains the */
/*         lu decomposition of p. */
/*         storage of matrix elements starts at wm(3). */
/*         wm also contains the following matrix-related data.. */
/*         wm(1) = sqrt(uround), used in numerical jacobian increments. */
/* iwm   = integer work space containing pivot information, starting at */
/*         iwm(21).   iwm also contains the band parameters */
/*         ml = iwm(1) and mu = iwm(2) if miter is 4 or 5. */
/* el0   = el(1) (input). */
/* pdnorm= norm of jacobian matrix. (output). */
/* ierpj = output error flag,  = 0 if no trouble, .gt. 0 if */
/*         p matrix found to be singular. */
/* jcur  = output flag = 1 to indicate that the jacobian matrix */
/*         (or approximation) is now current. */
/* this routine also uses the common variables el0, h, tn, uround, */
/* miter, n, nfe, and nje. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --neq;
    --y;
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --ewt;
    --ftem;
    --savf;
    --wm;
    --iwm;

    /* Function Body */
    ++ls0001_1.nje;
    ls0001_1.ierpj = 0;
    ls0001_1.jcur = 1;
    hl0 = ls0001_1.h__ * ls0001_1.el0;
    switch (ls0001_1.miter) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
    }
/* if miter = 1, call jac and multiply by scalar. ----------------------- */
L100:
    lenp = ls0001_1.n * ls0001_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	wm[i__ + 2] = 0.;
    }
    srcma_(rsav, isav, &c__1);
    (*jac)(&neq[1], &ls0001_1.tn, &y[1], &c__0, &c__0, &wm[3], &ls0001_1.n);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    srcma_(rsav, isav, &c__2);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	wm[i__ + 2] *= con;
    }
    goto L240;
/* if miter = 2, make n calls to f to approximate j. -------------------- */
L200:
    fac = vmnorm_(&ls0001_1.n, &savf[1], &ewt[1]);
    r0 = abs(ls0001_1.h__) * 1e3 * ls0001_1.uround * (doublereal) ls0001_1.n *
	     fac;
    if (r0 == 0.) {
	r0 = 1.;
    }
    srur = wm[1];
    j1 = 2;
    i__1 = ls0001_1.n;
    for (j = 1; j <= i__1; ++j) {
	yj = y[j];
/* Computing MAX */
	d__1 = srur * abs(yj), d__2 = r0 / ewt[j];
	r__ = max(d__1,d__2);
	y[j] += r__;
	fac = -hl0 / r__;
	srcma_(rsav, isav, &c__1);
	(*f)(&neq[1], &ls0001_1.tn, &y[1], &ftem[1]);
/*       SCIPY error check: */
	if (neq[1] == -1) {
	    return 0;
	}
	srcma_(rsav, isav, &c__2);
	i__2 = ls0001_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L220: */
	    wm[i__ + j1] = (ftem[i__] - savf[i__]) * fac;
	}
	y[j] = yj;
	j1 += ls0001_1.n;
/* L230: */
    }
    ls0001_1.nfe += ls0001_1.n;
L240:
/* compute norm of jacobian. -------------------------------------------- */
    lsa001_1.pdnorm = fnorm_(&ls0001_1.n, &wm[3], &ewt[1]) / abs(hl0);
/* add identity matrix. ------------------------------------------------- */
    np1 = ls0001_1.n + 1;
    j = 3;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[j] += 1.;
/* L250: */
	j += np1;
    }
/* do lu decomposition on p. -------------------------------------------- */
/*     Replaced LINPACK dgefa with LAPACK dgetrf */
/*      call dgefa (wm(3), n, n, iwm(21), ier) */
    dgetrf_(&ls0001_1.n, &ls0001_1.n, &wm[3], &ls0001_1.n, &iwm[21], &ier);
    if (ier != 0) {
	ls0001_1.ierpj = 1;
    }
    return 0;
/* dummy block only, since miter is never 3 in this routine. ------------ */
L300:
    return 0;
/* if miter = 4, call jac and multiply by scalar. ----------------------- */
L400:
    ml = iwm[1];
    mu = iwm[2];
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * ls0001_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	wm[i__ + 2] = 0.;
    }
    srcma_(rsav, isav, &c__1);
    (*jac)(&neq[1], &ls0001_1.tn, &y[1], &ml, &mu, &wm[ml3], &meband);
/*     SCIPY error check: */
    if (neq[1] == -1) {
	return 0;
    }
    srcma_(rsav, isav, &c__2);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L420: */
	wm[i__ + 2] *= con;
    }
    goto L570;
/* if miter = 5, make mband calls to f to approximate j. ---------------- */
L500:
    ml = iwm[1];
    mu = iwm[2];
    mband = ml + mu + 1;
    mba = min(mband,ls0001_1.n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = wm[1];
    fac = vmnorm_(&ls0001_1.n, &savf[1], &ewt[1]);
    r0 = abs(ls0001_1.h__) * 1e3 * ls0001_1.uround * (doublereal) ls0001_1.n *
	     fac;
    if (r0 == 0.) {
	r0 = 1.;
    }
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ls0001_1.n;
	i__3 = mband;
	for (i__ = j; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    yi = y[i__];
/* Computing MAX */
	    d__1 = srur * abs(yi), d__2 = r0 / ewt[i__];
	    r__ = max(d__1,d__2);
/* L530: */
	    y[i__] += r__;
	}
	srcma_(rsav, isav, &c__1);
	(*f)(&neq[1], &ls0001_1.tn, &y[1], &ftem[1]);
/*       SCIPY error check: */
	if (neq[1] == -1) {
	    return 0;
	}
	srcma_(rsav, isav, &c__2);
	i__3 = ls0001_1.n;
	i__2 = mband;
	for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    y[jj] = yh[jj + yh_dim1];
	    yjj = y[jj];
/* Computing MAX */
	    d__1 = srur * abs(yjj), d__2 = r0 / ewt[jj];
	    r__ = max(d__1,d__2);
	    fac = -hl0 / r__;
/* Computing MAX */
	    i__4 = jj - mu;
	    i1 = max(i__4,1);
/* Computing MIN */
	    i__4 = jj + ml;
	    i2 = min(i__4,ls0001_1.n);
	    ii = jj * meb1 - ml + 2;
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* L540: */
		wm[ii + i__] = (ftem[i__] - savf[i__]) * fac;
	    }
/* L550: */
	}
/* L560: */
    }
    ls0001_1.nfe += mba;
L570:
/* compute norm of jacobian. -------------------------------------------- */
    lsa001_1.pdnorm = bnorm_(&ls0001_1.n, &wm[3], &meband, &ml, &mu, &ewt[1]) 
	    / abs(hl0);
/* add identity matrix. ------------------------------------------------- */
    ii = mband + 2;
    i__1 = ls0001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[ii] += 1.;
/* L580: */
	ii += meband;
    }
/* do lu decomposition of p. -------------------------------------------- */
/*     Replaced LINPACK dgefa with LAPACK dgetrf */
/*      call dgbfa (wm(3), meband, n, ml, mu, iwm(21), ier) */
    dgbtrf_(&ls0001_1.n, &ls0001_1.n, &ml, &mu, &wm[3], &meband, &iwm[21], &
	    ier);
    if (ier != 0) {
	ls0001_1.ierpj = 1;
    }
    return 0;
/* ----------------------- end of subroutine prja ------------------------ */
} /* prja_ */

