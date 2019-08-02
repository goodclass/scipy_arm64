/* curfit.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int curfit_(integer *iopt, integer *m, doublereal *x, 
	doublereal *y, doublereal *w, doublereal *xb, doublereal *xe, integer 
	*k, doublereal *s, integer *nest, integer *n, doublereal *t, 
	doublereal *c__, doublereal *fp, doublereal *wrk, integer *lwrk, 
	integer *iwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k1, k2, ia, ib, ig, iq, iz, ifp;
    static doublereal tol;
    static integer nmin, maxit, lwest;
    extern /* Subroutine */ int fpchec_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *), fpcurf_(integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);

/*  given the set of data points (x(i),y(i)) and the set of positive */
/*  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline */
/*  approximation of degree k on the interval xb <= x <= xe. */
/*  if iopt=-1 curfit calculates the weighted least-squares spline */
/*  according to a given set of knots. */
/*  if iopt>=0 the number of knots of the spline s(x) and the position */
/*  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth- */
/*  ness of s(x) is then achieved by minimalizing the discontinuity */
/*  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,..., */
/*  n-k-1. the amount of smoothness is determined by the condition that */
/*  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non- */
/*  negative constant, called the smoothing factor. */
/*  the fit s(x) is given in the b-spline representation (b-spline coef- */
/*  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of */
/*  subroutine splev. */

/*  calling sequence: */
/*     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk, */
/*    * lwrk,iwrk,ier) */

/*  parameters: */
/*   iopt  : integer flag. on entry iopt must specify whether a weighted */
/*           least-squares spline (iopt=-1) or a smoothing spline (iopt= */
/*           0 or 1) must be determined. if iopt=0 the routine will start */
/*           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,... */
/*           k+1. if iopt=1 the routine will continue with the knots */
/*           found at the last call of the routine. */
/*           attention: a call with iopt=1 must always be immediately */
/*           preceded by another call with iopt=1 or iopt=0. */
/*           unchanged on exit. */
/*   m     : integer. on entry m must specify the number of data points. */
/*           m > k. unchanged on exit. */
/*   x     : real array of dimension at least (m). before entry, x(i) */
/*           must be set to the i-th value of the independent variable x, */
/*           for i=1,2,...,m. these values must be supplied in strictly */
/*           ascending order. unchanged on exit. */
/*   y     : real array of dimension at least (m). before entry, y(i) */
/*           must be set to the i-th value of the dependent variable y, */
/*           for i=1,2,...,m. unchanged on exit. */
/*   w     : real array of dimension at least (m). before entry, w(i) */
/*           must be set to the i-th value in the set of weights. the */
/*           w(i) must be strictly positive. unchanged on exit. */
/*           see also further comments. */
/*   xb,xe : real values. on entry xb and xe must specify the boundaries */
/*           of the approximation interval. xb<=x(1), xe>=x(m). */
/*           unchanged on exit. */
/*   k     : integer. on entry k must specify the degree of the spline. */
/*           1<=k<=5. it is recommended to use cubic splines (k=3). */
/*           the user is strongly dissuaded from choosing k even,together */
/*           with a small s-value. unchanged on exit. */
/*   s     : real.on entry (in case iopt>=0) s must specify the smoothing */
/*           factor. s >=0. unchanged on exit. */
/*           for advice on the choice of s see further comments. */
/*   nest  : integer. on entry nest must contain an over-estimate of the */
/*           total number of knots of the spline returned, to indicate */
/*           the storage space available to the routine. nest >=2*k+2. */
/*           in most practical situation nest=m/2 will be sufficient. */
/*           always large enough is  nest=m+k+1, the number of knots */
/*           needed for interpolation (s=0). unchanged on exit. */
/*   n     : integer. */
/*           unless ier =10 (in case iopt >=0), n will contain the */
/*           total number of knots of the spline approximation returned. */
/*           if the computation mode iopt=1 is used this value of n */
/*           should be left unchanged between subsequent calls. */
/*           in case iopt=-1, the value of n must be specified on entry. */
/*   t     : real array of dimension at least (nest). */
/*           on successful exit, this array will contain the knots of the */
/*           spline,i.e. the position of the interior knots t(k+2),t(k+3) */
/*           ...,t(n-k-1) as well as the position of the additional knots */
/*           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for */
/*           the b-spline representation. */
/*           if the computation mode iopt=1 is used, the values of t(1), */
/*           t(2),...,t(n) should be left unchanged between subsequent */
/*           calls. if the computation mode iopt=-1 is used, the values */
/*           t(k+2),...,t(n-k-1) must be supplied by the user, before */
/*           entry. see also the restrictions (ier=10). */
/*   c     : real array of dimension at least (nest). */
/*           on successful exit, this array will contain the coefficients */
/*           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x) */
/*   fp    : real. unless ier=10, fp contains the weighted sum of */
/*           squared residuals of the spline approximation returned. */
/*   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)). */
/*           used as working space. if the computation mode iopt=1 is */
/*           used, the values wrk(1),...,wrk(n) should be left unchanged */
/*           between subsequent calls. */
/*   lwrk  : integer. on entry,lwrk must specify the actual dimension of */
/*           the array wrk as declared in the calling (sub)program.lwrk */
/*           must not be too small (see wrk). unchanged on exit. */
/*   iwrk  : integer array of dimension at least (nest). */
/*           used as working space. if the computation mode iopt=1 is */
/*           used,the values iwrk(1),...,iwrk(n) should be left unchanged */
/*           between subsequent calls. */
/*   ier   : integer. unless the routine detects an error, ier contains a */
/*           non-positive value on exit, i.e. */
/*    ier=0  : normal return. the spline returned has a residual sum of */
/*             squares fp such that abs(fp-s)/s <= tol with tol a relat- */
/*             ive tolerance set to 0.001 by the program. */
/*    ier=-1 : normal return. the spline returned is an interpolating */
/*             spline (fp=0). */
/*    ier=-2 : normal return. the spline returned is the weighted least- */
/*             squares polynomial of degree k. in this extreme case fp */
/*             gives the upper bound fp0 for the smoothing factor s. */
/*    ier=1  : error. the required storage space exceeds the available */
/*             storage space, as specified by the parameter nest. */
/*             probably causes : nest too small. if nest is already */
/*             large (say nest > m/2), it may also indicate that s is */
/*             too small */
/*             the approximation returned is the weighted least-squares */
/*             spline according to the knots t(1),t(2),...,t(n). (n=nest) */
/*             the parameter fp gives the corresponding weighted sum of */
/*             squared residuals (fp>s). */
/*    ier=2  : error. a theoretically impossible result was found during */
/*             the iteration process for finding a smoothing spline with */
/*             fp = s. probably causes : s too small. */
/*             there is an approximation returned but the corresponding */
/*             weighted sum of squared residuals does not satisfy the */
/*             condition abs(fp-s)/s < tol. */
/*    ier=3  : error. the maximal number of iterations maxit (set to 20 */
/*             by the program) allowed for finding a smoothing spline */
/*             with fp=s has been reached. probably causes : s too small */
/*             there is an approximation returned but the corresponding */
/*             weighted sum of squared residuals does not satisfy the */
/*             condition abs(fp-s)/s < tol. */
/*    ier=10 : error. on entry, the input data are controlled on validity */
/*             the following restrictions must be satisfied. */
/*             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m */
/*             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k) */
/*             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1) */
/*                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe */
/*                       the schoenberg-whitney conditions, i.e. there */
/*                       must be a subset of data points xx(j) such that */
/*                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1 */
/*             if iopt>=0: s>=0 */
/*                         if s=0 : nest >= m+k+1 */
/*             if one of these conditions is found to be violated,control */
/*             is immediately repassed to the calling program. in that */
/*             case there is no approximation returned. */

/*  further comments: */
/*   by means of the parameter s, the user can control the tradeoff */
/*   between closeness of fit and smoothness of fit of the approximation. */
/*   if s is too large, the spline will be too smooth and signal will be */
/*   lost ; if s is too small the spline will pick up too much noise. in */
/*   the extreme cases the program will return an interpolating spline if */
/*   s=0 and the weighted least-squares polynomial of degree k if s is */
/*   very large. between these extremes, a properly chosen s will result */
/*   in a good compromise between closeness of fit and smoothness of fit. */
/*   to decide whether an approximation, corresponding to a certain s is */
/*   satisfactory the user is highly recommended to inspect the fits */
/*   graphically. */
/*   recommended values for s depend on the weights w(i). if these are */
/*   taken as 1/d(i) with d(i) an estimate of the standard deviation of */
/*   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+ */
/*   sqrt(2*m)). if nothing is known about the statistical error in y(i) */
/*   each w(i) can be set equal to one and s determined by trial and */
/*   error, taking account of the comments above. the best is then to */
/*   start with a very large value of s ( to determine the least-squares */
/*   polynomial and the corresponding upper bound fp0 for s) and then to */
/*   progressively decrease the value of s ( say by a factor 10 in the */
/*   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the */
/*   approximation shows more detail) to obtain closer fits. */
/*   to economize the search for a good s-value the program provides with */
/*   different modes of computation. at the first call of the routine, or */
/*   whenever he wants to restart with the initial set of knots the user */
/*   must set iopt=0. */
/*   if iopt=1 the program will continue with the set of knots found at */
/*   the last call of the routine. this will save a lot of computation */
/*   time if curfit is called repeatedly for different values of s. */
/*   the number of knots of the spline returned and their location will */
/*   depend on the value of s and on the complexity of the shape of the */
/*   function underlying the data. but, if the computation mode iopt=1 */
/*   is used, the knots returned may also depend on the s-values at */
/*   previous calls (if these were smaller). therefore, if after a number */
/*   of trials with different s-values and iopt=1, the user can finally */
/*   accept a fit as satisfactory, it may be worthwhile for him to call */
/*   curfit once more with the selected value for s but now with iopt=0. */
/*   indeed, curfit may then return an approximation of the same quality */
/*   of fit but with fewer knots and therefore better if data reduction */
/*   is also an important objective for the user. */

/*  other subroutines required: */
/*    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota */

/*  references: */
/*   dierckx p. : an algorithm for smoothing, differentiation and integ- */
/*                ration of experimental data using spline functions, */
/*                j.comp.appl.maths 1 (1975) 165-184. */
/*   dierckx p. : a fast algorithm for smoothing data on a rectangular */
/*                grid while using spline functions, siam j.numer.anal. */
/*                19 (1982) 1286-1304. */
/*   dierckx p. : an improved algorithm for curve fitting with spline */
/*                functions, report tw54, dept. computer science,k.u. */
/*                leuven, 1981. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author: */
/*    p.dierckx */
/*    dept. computer science, k.u. leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  creation date : may 1979 */
/*  latest update : march 1987 */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  we set up the parameters tol and maxit */
    /* Parameter adjustments */
    --w;
    --y;
    --x;
    --iwrk;
    --c__;
    --t;
    --wrk;

    /* Function Body */
    maxit = 20;
    tol = .001;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (*k <= 0 || *k > 5) {
	goto L50;
    }
    k1 = *k + 1;
    k2 = k1 + 1;
    if (*iopt < -1 || *iopt > 1) {
	goto L50;
    }
    nmin = k1 << 1;
    if (*m < k1 || *nest < nmin) {
	goto L50;
    }
    lwest = *m * k1 + *nest * (*k * 3 + 7);
    if (*lwrk < lwest) {
	goto L50;
    }
    if (*xb > x[1] || *xe < x[*m]) {
	goto L50;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__ - 1] > x[i__]) {
	    goto L50;
	}
/* L10: */
    }
    if (*iopt >= 0) {
	goto L30;
    }
    if (*n < nmin || *n > *nest) {
	goto L50;
    }
    j = *n;
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__] = *xb;
	t[j] = *xe;
	--j;
/* L20: */
    }
    fpchec_(&x[1], m, &t[1], n, k, ier);
    if (*ier == 0) {
	goto L40;
    }
    goto L50;
L30:
    if (*s < 0.f) {
	goto L50;
    }
    if (*s == 0.f && *nest < *m + k1) {
	goto L50;
    }
/* we partition the working space and determine the spline approximation. */
L40:
    ifp = 1;
    iz = ifp + *nest;
    ia = iz + *nest;
    ib = ia + *nest * k1;
    ig = ib + *nest * k2;
    iq = ig + *nest * k2;
    fpcurf_(iopt, &x[1], &y[1], &w[1], m, xb, xe, k, s, nest, &tol, &maxit, &
	    k1, &k2, n, &t[1], &c__[1], fp, &wrk[ifp], &wrk[iz], &wrk[ia], &
	    wrk[ib], &wrk[ig], &wrk[iq], &iwrk[1], ier);
L50:
    return 0;
} /* curfit_ */

