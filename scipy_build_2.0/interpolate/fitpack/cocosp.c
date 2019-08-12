/* cocosp.f -- translated by f2c (version 20190311).
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

static integer c__3 = 3;

/* Subroutine */ int cocosp_(integer *m, doublereal *x, doublereal *y, 
	doublereal *w, integer *n, doublereal *t, doublereal *e, integer *
	maxtr, integer *maxbin, doublereal *c__, doublereal *sq, doublereal *
	sx, logical *bind, doublereal *wrk, integer *lwrk, integer *iwrk, 
	integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n6, ia, ib, ic, mb, ji, jl, iq, nm, jr, iu, ju, iz, 
	    jib, jjb;
    static doublereal one;
    static integer izz, kwest, lwest;
    extern /* Subroutine */ int fpchec_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *), fpcosp_(integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     logical *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *);

/*  given the set of data points (x(i),y(i)) and the set of positive */
/*  numbers w(i),i=1,2,...,m, subroutine cocosp determines the weighted */
/*  least-squares cubic spline s(x) with given knots t(j),j=1,2,...,n */
/*  which satisfies the following concavity/convexity conditions */
/*      s''(t(j+3))*e(j) <= 0, j=1,2,...n-6 */
/*  the fit is given in the b-spline representation( b-spline coef- */
/*  ficients c(j),j=1,2,...n-4) and can be evaluated by means of */
/*  subroutine splev. */

/*  calling sequence: */
/*     call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk, */
/*    * lwrk,iwrk,kwrk,ier) */

/*  parameters: */
/*    m   : integer. on entry m must specify the number of data points. */
/*          m > 3. unchanged on exit. */
/*    x   : real array of dimension at least (m). before entry, x(i) */
/*          must be set to the i-th value of the independent variable x, */
/*          for i=1,2,...,m. these values must be supplied in strictly */
/*          ascending order. unchanged on exit. */
/*    y   : real array of dimension at least (m). before entry, y(i) */
/*          must be set to the i-th value of the dependent variable y, */
/*          for i=1,2,...,m. unchanged on exit. */
/*    w   : real array of dimension at least (m). before entry, w(i) */
/*          must be set to the i-th value in the set of weights. the */
/*          w(i) must be strictly positive. unchanged on exit. */
/*    n   : integer. on entry n must contain the total number of knots */
/*          of the cubic spline. m+4>=n>=8. unchanged on exit. */
/*    t   : real array of dimension at least (n). before entry, this */
/*          array must contain the knots of the spline, i.e. the position */
/*          of the interior knots t(5),t(6),...,t(n-4) as well as the */
/*          position of the boundary knots t(1),t(2),t(3),t(4) and t(n-3) */
/*          t(n-2),t(n-1),t(n) needed for the b-spline representation. */
/*          unchanged on exit. see also the restrictions (ier=10). */
/*    e   : real array of dimension at least (n). before entry, e(j) */
/*          must be set to 1 if s(x) must be locally concave at t(j+3), */
/*          to (-1) if s(x) must be locally convex at t(j+3) and to 0 */
/*          if no convexity constraint is imposed at t(j+3),j=1,2,..,n-6. */
/*          e(n-5),...,e(n) are not used. unchanged on exit. */
/*  maxtr : integer. on entry maxtr must contain an over-estimate of the */
/*          total number of records in the used tree structure, to indic- */
/*          ate the storage space available to the routine. maxtr >=1 */
/*          in most practical situation maxtr=100 will be sufficient. */
/*          always large enough is */
/*                         n-5       n-6 */
/*              maxtr =  (     ) + (     )  with l the greatest */
/*                          l        l+1 */
/*          integer <= (n-6)/2 . unchanged on exit. */
/*  maxbin: integer. on entry maxbin must contain an over-estimate of the */
/*          number of knots where s(x) will have a zero second derivative */
/*          maxbin >=1. in most practical situation maxbin = 10 will be */
/*          sufficient. always large enough is maxbin=n-6. */
/*          unchanged on exit. */
/*    c   : real array of dimension at least (n). */
/*          on successful exit, this array will contain the coefficients */
/*          c(1),c(2),..,c(n-4) in the b-spline representation of s(x) */
/*    sq  : real. on successful exit, sq contains the weighted sum of */
/*          squared residuals of the spline approximation returned. */
/*    sx  : real array of dimension at least m. on successful exit */
/*          this array will contain the spline values s(x(i)),i=1,...,m */
/*   bind : logical array of dimension at least (n). on successful exit */
/*          this array will indicate the knots where s''(x)=0, i.e. */
/*                s''(t(j+3)) .eq. 0 if  bind(j) = .true. */
/*                s''(t(j+3)) .ne. 0 if  bind(j) = .false., j=1,2,...,n-6 */
/*   wrk  : real array of dimension at least  m*4+n*7+maxbin*(maxbin+n+1) */
/*          used as working space. */
/*   lwrk : integer. on entry,lwrk must specify the actual dimension of */
/*          the array wrk as declared in the calling (sub)program.lwrk */
/*          must not be too small (see wrk). unchanged on exit. */
/*   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1)) */
/*          used as working space. */
/*   kwrk : integer. on entry,kwrk must specify the actual dimension of */
/*          the array iwrk as declared in the calling (sub)program. kwrk */
/*          must not be too small (see iwrk). unchanged on exit. */
/*   ier   : integer. error flag */
/*      ier=0 : successful exit. */
/*      ier>0 : abnormal termination: no approximation is returned */
/*        ier=1  : the number of knots where s''(x)=0 exceeds maxbin. */
/*                 probably causes : maxbin too small. */
/*        ier=2  : the number of records in the tree structure exceeds */
/*                 maxtr. */
/*                 probably causes : maxtr too small. */
/*        ier=3  : the algorithm finds no solution to the posed quadratic */
/*                 programming problem. */
/*                 probably causes : rounding errors. */
/*        ier=10 : on entry, the input data are controlled on validity. */
/*                 the following restrictions must be satisfied */
/*                   m>3, maxtr>=1, maxbin>=1, 8<=n<=m+4,w(i) > 0, */
/*                   x(1)<x(2)<...<x(m), t(1)<=t(2)<=t(3)<=t(4)<=x(1), */
/*                   x(1)<t(5)<t(6)<...<t(n-4)<x(m)<=t(n-3)<=...<=t(n), */
/*                   kwrk>=maxtr*4+2*(maxbin+1), */
/*                   lwrk>=m*4+n*7+maxbin*(maxbin+n+1), */
/*                   the schoenberg-whitney conditions, i.e. there must */
/*                   be a subset of data points xx(j) such that */
/*                     t(j) < xx(j) < t(j+4), j=1,2,...,n-4 */
/*                 if one of these restrictions is found to be violated */
/*                 control is immediately repassed to the calling program */


/*  other subroutines required: */
/*    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno,fpchec */

/*  references: */
/*   dierckx p. : an algorithm for cubic spline fitting with convexity */
/*                constraints, computing 24 (1980) 349-371. */
/*   dierckx p. : an algorithm for least-squares cubic spline fitting */
/*                with convexity and concavity constraints, report tw39, */
/*                dept. computer science, k.u.leuven, 1978. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author: */
/*   p. dierckx */
/*   dept. computer science, k.u.leuven */
/*   celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*   e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  creation date : march 1978 */
/*  latest update : march 1987. */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  set constant */
    /* Parameter adjustments */
    --sx;
    --w;
    --y;
    --x;
    --bind;
    --c__;
    --e;
    --t;
    --wrk;
    --iwrk;

    /* Function Body */
    one = 1.f;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (*m < 4 || *n < 8) {
	goto L40;
    }
    if (*maxtr < 1 || *maxbin < 1) {
	goto L40;
    }
    lwest = *n * 7 + (*m << 2) + *maxbin * (*n + 1 + *maxbin);
    kwest = (*maxtr << 2) + (*maxbin + 1 << 1);
    if (*lwrk < lwest || *kwrk < kwest) {
	goto L40;
    }
    if (w[1] <= 0.f) {
	goto L40;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__ - 1] >= x[i__] || w[i__] <= 0.f) {
	    goto L40;
	}
/* L10: */
    }
    fpchec_(&x[1], m, &t[1], n, &c__3, ier);
    if (*ier == 0) {
	goto L20;
    }
    goto L40;
/*  set numbers e(i) */
L20:
    n6 = *n - 6;
    i__1 = n6;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] > 0.f) {
	    e[i__] = one;
	}
	if (e[i__] < 0.f) {
	    e[i__] = -one;
	}
/* L30: */
    }
/*  we partition the working space and determine the spline approximation */
    nm = *n + *maxbin;
    mb = *maxbin + 1;
    ia = 1;
    ib = ia + (*n << 2);
    ic = ib + nm * *maxbin;
    iz = ic + *n;
    izz = iz + *n;
    iu = izz + *n;
    iq = iu + *maxbin;
    ji = 1;
    ju = ji + *maxtr;
    jl = ju + *maxtr;
    jr = jl + *maxtr;
    jjb = jr + *maxtr;
    jib = jjb + mb;
    fpcosp_(m, &x[1], &y[1], &w[1], n, &t[1], &e[1], maxtr, maxbin, &c__[1], 
	    sq, &sx[1], &bind[1], &nm, &mb, &wrk[ia], &wrk[ib], &wrk[ic], &
	    wrk[iz], &wrk[izz], &wrk[iu], &wrk[iq], &iwrk[ji], &iwrk[ju], &
	    iwrk[jl], &iwrk[jr], &iwrk[jjb], &iwrk[jib], ier);
L40:
    return 0;
} /* cocosp_ */

