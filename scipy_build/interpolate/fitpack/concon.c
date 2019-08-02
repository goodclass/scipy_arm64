/* concon.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int concon_(integer *iopt, integer *m, doublereal *x, 
	doublereal *y, doublereal *w, doublereal *v, doublereal *s, integer *
	nest, integer *maxtr, integer *maxbin, integer *n, doublereal *t, 
	doublereal *c__, doublereal *sq, doublereal *sx, logical *bind, 
	doublereal *wrk, integer *lwrk, integer *iwrk, integer *kwrk, integer 
	*ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ie, iw;
    static doublereal one;
    static integer lww, kwest, lwest;
    extern /* Subroutine */ int fpcoco_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *);

/*  given the set of data points (x(i),y(i)) and the set of positive */
/*  numbers w(i), i=1,2,...,m,subroutine concon determines a cubic spline */
/*  approximation s(x) which satisfies the following local convexity */
/*  constraints  s''(x(i))*v(i) <= 0, i=1,2,...,m. */
/*  the number of knots n and the position t(j),j=1,2,...n is chosen */
/*  automatically by the routine in a way that */
/*       sq = sum((w(i)*(y(i)-s(x(i))))**2) be <= s. */
/*  the fit is given in the b-spline representation (b-spline coef- */
/*  ficients c(j),j=1,2,...n-4) and can be evaluated by means of */
/*  subroutine splev. */

/*  calling sequence: */

/*     call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq, */
/*    * sx,bind,wrk,lwrk,iwrk,kwrk,ier) */

/*  parameters: */
/*    iopt: integer flag. */
/*          if iopt=0, the routine will start with the minimal number of */
/*          knots to guarantee that the convexity conditions will be */
/*          satisfied. if iopt=1, the routine will continue with the set */
/*          of knots found at the last call of the routine. */
/*          attention: a call with iopt=1 must always be immediately */
/*          preceded by another call with iopt=1 or iopt=0. */
/*          unchanged on exit. */
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
/*    v   : real array of dimension at least (m). before entry, v(i) */
/*          must be set to 1 if s(x) must be locally concave at x(i), */
/*          to (-1) if s(x) must be locally convex at x(i) and to 0 */
/*          if no convexity constraint is imposed at x(i). */
/*    s   : real. on entry s must specify an over-estimate for the */
/*          the weighted sum of squared residuals sq of the requested */
/*          spline. s >=0. unchanged on exit. */
/*   nest : integer. on entry nest must contain an over-estimate of the */
/*          total number of knots of the spline returned, to indicate */
/*          the storage space available to the routine. nest >=8. */
/*          in most practical situation nest=m/2 will be sufficient. */
/*          always large enough is  nest=m+4. unchanged on exit. */
/*  maxtr : integer. on entry maxtr must contain an over-estimate of the */
/*          total number of records in the used tree structure, to indic- */
/*          ate the storage space available to the routine. maxtr >=1 */
/*          in most practical situation maxtr=100 will be sufficient. */
/*          always large enough is */
/*                         nest-5      nest-6 */
/*              maxtr =  (       ) + (        )  with l the greatest */
/*                           l          l+1 */
/*          integer <= (nest-6)/2 . unchanged on exit. */
/*  maxbin: integer. on entry maxbin must contain an over-estimate of the */
/*          number of knots where s(x) will have a zero second derivative */
/*          maxbin >=1. in most practical situation maxbin = 10 will be */
/*          sufficient. always large enough is maxbin=nest-6. */
/*          unchanged on exit. */
/*    n   : integer. */
/*          on exit with ier <=0, n will contain the total number of */
/*          knots of the spline approximation returned. if the comput- */
/*          ation mode iopt=1 is used this value of n should be left */
/*          unchanged between subsequent calls. */
/*    t   : real array of dimension at least (nest). */
/*          on exit with ier<=0, this array will contain the knots of the */
/*          spline,i.e. the position of the interior knots t(5),t(6),..., */
/*          t(n-4) as well as the position of the additional knots */
/*          t(1)=t(2)=t(3)=t(4)=x(1) and t(n-3)=t(n-2)=t(n-1)=t(n)=x(m) */
/*          needed for the the b-spline representation. */
/*          if the computation mode iopt=1 is used, the values of t(1), */
/*          t(2),...,t(n) should be left unchanged between subsequent */
/*          calls. */
/*    c   : real array of dimension at least (nest). */
/*          on successful exit, this array will contain the coefficients */
/*          c(1),c(2),..,c(n-4) in the b-spline representation of s(x) */
/*    sq  : real. unless ier>0 , sq contains the weighted sum of */
/*          squared residuals of the spline approximation returned. */
/*    sx  : real array of dimension at least m. on exit with ier<=0 */
/*          this array will contain the spline values s(x(i)),i=1,...,m */
/*          if the computation mode iopt=1 is used, the values of sx(1), */
/*          sx(2),...,sx(m) should be left unchanged between subsequent */
/*          calls. */
/*    bind: logical array of dimension at least nest. on exit with ier<=0 */
/*          this array will indicate the knots where s''(x)=0, i.e. */
/*                s''(t(j+3)) .eq. 0 if  bind(j) = .true. */
/*                s''(t(j+3)) .ne. 0 if  bind(j) = .false., j=1,2,...,n-6 */
/*          if the computation mode iopt=1 is used, the values of bind(1) */
/*          ,...,bind(n-6) should be left unchanged between subsequent */
/*          calls. */
/*   wrk  : real array of dimension at least (m*4+nest*8+maxbin*(maxbin+ */
/*          nest+1)). used as working space. */
/*   lwrk : integer. on entry,lwrk must specify the actual dimension of */
/*          the array wrk as declared in the calling (sub)program.lwrk */
/*          must not be too small (see wrk). unchanged on exit. */
/*   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1)) */
/*          used as working space. */
/*   kwrk : integer. on entry,kwrk must specify the actual dimension of */
/*          the array iwrk as declared in the calling (sub)program. kwrk */
/*          must not be too small (see iwrk). unchanged on exit. */
/*   ier   : integer. error flag */
/*      ier=0 : normal return, s(x) satisfies the concavity/convexity */
/*              constraints and sq <= s. */
/*      ier<0 : abnormal termination: s(x) satisfies the concavity/ */
/*              convexity constraints but sq > s. */
/*        ier=-3 : the requested storage space exceeds the available */
/*                 storage space as specified by the parameter nest. */
/*                 probably causes: nest too small. if nest is already */
/*                 large (say nest > m/2), it may also indicate that s */
/*                 is too small. */
/*                 the approximation returned is the least-squares cubic */
/*                 spline according to the knots t(1),...,t(n) (n=nest) */
/*                 which satisfies the convexity constraints. */
/*        ier=-2 : the maximal number of knots n=m+4 has been reached. */
/*                 probably causes: s too small. */
/*        ier=-1 : the number of knots n is less than the maximal number */
/*                 m+4 but concon finds that adding one or more knots */
/*                 will not further reduce the value of sq. */
/*                 probably causes : s too small. */
/*      ier>0 : abnormal termination: no approximation is returned */
/*        ier=1  : the number of knots where s''(x)=0 exceeds maxbin. */
/*                 probably causes : maxbin too small. */
/*        ier=2  : the number of records in the tree structure exceeds */
/*                 maxtr. */
/*                 probably causes : maxtr too small. */
/*        ier=3  : the algorithm finds no solution to the posed quadratic */
/*                 programming problem. */
/*                 probably causes : rounding errors. */
/*        ier=4  : the minimum number of knots (given by n) to guarantee */
/*                 that the concavity/convexity conditions will be */
/*                 satisfied is greater than nest. */
/*                 probably causes: nest too small. */
/*        ier=5  : the minimum number of knots (given by n) to guarantee */
/*                 that the concavity/convexity conditions will be */
/*                 satisfied is greater than m+4. */
/*                 probably causes: strongly alternating convexity and */
/*                 concavity conditions. normally the situation can be */
/*                 coped with by adding n-m-4 extra data points (found */
/*                 by linear interpolation e.g.) with a small weight w(i) */
/*                 and a v(i) number equal to zero. */
/*        ier=10 : on entry, the input data are controlled on validity. */
/*                 the following restrictions must be satisfied */
/*                   0<=iopt<=1, m>3, nest>=8, s>=0, maxtr>=1, maxbin>=1, */
/*                   kwrk>=maxtr*4+2*(maxbin+1), w(i)>0, x(i) < x(i+1), */
/*                   lwrk>=m*4+nest*8+maxbin*(maxbin+nest+1) */
/*                 if one of these restrictions is found to be violated */
/*                 control is immediately repassed to the calling program */

/*  further comments: */
/*    as an example of the use of the computation mode iopt=1, the */
/*    following program segment will cause concon to return control */
/*    each time a spline with a new set of knots has been computed. */
/*     ............. */
/*     iopt = 0 */
/*     s = 0.1e+60  (s very large) */
/*     do 10 i=1,m */
/*       call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx, */
/*    *  bind,wrk,lwrk,iwrk,kwrk,ier) */
/*       ...... */
/*       s = sq */
/*       iopt=1 */
/* 10  continue */
/*     ............. */

/*  other subroutines required: */
/*    fpcoco,fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno */

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
    --v;
    --w;
    --y;
    --x;
    --bind;
    --c__;
    --t;
    --wrk;
    --iwrk;

    /* Function Body */
    one = 1.f;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (*iopt < 0 || *iopt > 1) {
	goto L30;
    }
    if (*m < 4 || *nest < 8) {
	goto L30;
    }
    if (*s < 0.f) {
	goto L30;
    }
    if (*maxtr < 1 || *maxbin < 1) {
	goto L30;
    }
    lwest = (*nest << 3) + (*m << 2) + *maxbin * (*nest + 1 + *maxbin);
    kwest = (*maxtr << 2) + (*maxbin + 1 << 1);
    if (*lwrk < lwest || *kwrk < kwest) {
	goto L30;
    }
    if (*iopt > 0) {
	goto L20;
    }
    if (w[1] <= 0.f) {
	goto L30;
    }
    if (v[1] > 0.f) {
	v[1] = one;
    }
    if (v[1] < 0.f) {
	v[1] = -one;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__ - 1] >= x[i__] || w[i__] <= 0.f) {
	    goto L30;
	}
	if (v[i__] > 0.f) {
	    v[i__] = one;
	}
	if (v[i__] < 0.f) {
	    v[i__] = -one;
	}
/* L10: */
    }
L20:
    *ier = 0;
/*  we partition the working space and determine the spline approximation */
    ie = 1;
    iw = ie + *nest;
    lww = *lwrk - *nest;
    fpcoco_(iopt, m, &x[1], &y[1], &w[1], &v[1], s, nest, maxtr, maxbin, n, &
	    t[1], &c__[1], sq, &sx[1], &bind[1], &wrk[ie], &wrk[iw], &lww, &
	    iwrk[1], kwrk, ier);
L30:
    return 0;
} /* concon_ */

