/* clocur.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int clocur_(integer *iopt, integer *ipar, integer *idim, 
	integer *m, doublereal *u, integer *mx, doublereal *x, doublereal *w, 
	integer *k, doublereal *s, integer *nest, integer *n, doublereal *t, 
	integer *nc, doublereal *c__, doublereal *fp, doublereal *wrk, 
	integer *lwrk, integer *iwrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, i1, i2, j1, j2, k1, k2, m1, ib, iq, iz, ia1, ia2, 
	    ig1, ig2, ncc, ifp;
    static doublereal per, tol;
    static integer nmin;
    static doublereal dist;
    static integer maxit, lwest;
    extern /* Subroutine */ int fpchep_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *), fpclos_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

/*  given the ordered set of m points x(i) in the idim-dimensional space */
/*  with x(1)=x(m), and given also a corresponding set of strictly in- */
/*  creasing values u(i) and the set of positive numbers w(i),i=1,2,...,m */
/*  subroutine clocur determines a smooth approximating closed spline */
/*  curve s(u), i.e. */
/*      x1 = s1(u) */
/*      x2 = s2(u)       u(1) <= u <= u(m) */
/*      ......... */
/*      xidim = sidim(u) */
/*  with sj(u),j=1,2,...,idim periodic spline functions of degree k with */
/*  common knots t(j),j=1,2,...,n. */
/*  if ipar=1 the values u(i),i=1,2,...,m must be supplied by the user. */
/*  if ipar=0 these values are chosen automatically by clocur as */
/*      v(1) = 0 */
/*      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m */
/*      u(i) = v(i)/v(m) ,i=1,2,...,m */
/*  if iopt=-1 clocur calculates the weighted least-squares closed spline */
/*  curve according to a given set of knots. */
/*  if iopt>=0 the number of knots of the splines sj(u) and the position */
/*  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth- */
/*  ness of s(u) is then achieved by minimalizing the discontinuity */
/*  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,..., */
/*  n-k-1. the amount of smoothness is determined by the condition that */
/*  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non- */
/*  negative constant, called the smoothing factor. */
/*  the fit s(u) is given in the b-spline representation and can be */
/*  evaluated by means of subroutine curev. */

/*  calling sequence: */
/*     call clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c, */
/*    * fp,wrk,lwrk,iwrk,ier) */

/*  parameters: */
/*   iopt  : integer flag. on entry iopt must specify whether a weighted */
/*           least-squares closed spline curve (iopt=-1) or a smoothing */
/*           closed spline curve (iopt=0 or 1) must be determined. if */
/*           iopt=0 the routine will start with an initial set of knots */
/*           t(i)=u(1)+(u(m)-u(1))*(i-k-1),i=1,2,...,2*k+2. if iopt=1 the */
/*           routine will continue with the knots found at the last call. */
/*           attention: a call with iopt=1 must always be immediately */
/*           preceded by another call with iopt=1 or iopt=0. */
/*           unchanged on exit. */
/*   ipar  : integer flag. on entry ipar must specify whether (ipar=1) */
/*           the user will supply the parameter values u(i),or whether */
/*           (ipar=0) these values are to be calculated by clocur. */
/*           unchanged on exit. */
/*   idim  : integer. on entry idim must specify the dimension of the */
/*           curve. 0 < idim < 11. */
/*           unchanged on exit. */
/*   m     : integer. on entry m must specify the number of data points. */
/*           m > 1. unchanged on exit. */
/*   u     : real array of dimension at least (m). in case ipar=1,before */
/*           entry, u(i) must be set to the i-th value of the parameter */
/*           variable u for i=1,2,...,m. these values must then be */
/*           supplied in strictly ascending order and will be unchanged */
/*           on exit. in case ipar=0, on exit,the array will contain the */
/*           values u(i) as determined by clocur. */
/*   mx    : integer. on entry mx must specify the actual dimension of */
/*           the array x as declared in the calling (sub)program. mx must */
/*           not be too small (see x). unchanged on exit. */
/*   x     : real array of dimension at least idim*m. */
/*           before entry, x(idim*(i-1)+j) must contain the j-th coord- */
/*           inate of the i-th data point for i=1,2,...,m and j=1,2,..., */
/*           idim. since first and last data point must coincide it */
/*           means that x(j)=x(idim*(m-1)+j),j=1,2,...,idim. */
/*           unchanged on exit. */
/*   w     : real array of dimension at least (m). before entry, w(i) */
/*           must be set to the i-th value in the set of weights. the */
/*           w(i) must be strictly positive. w(m) is not used. */
/*           unchanged on exit. see also further comments. */
/*   k     : integer. on entry k must specify the degree of the splines. */
/*           1<=k<=5. it is recommended to use cubic splines (k=3). */
/*           the user is strongly dissuaded from choosing k even,together */
/*           with a small s-value. unchanged on exit. */
/*   s     : real.on entry (in case iopt>=0) s must specify the smoothing */
/*           factor. s >=0. unchanged on exit. */
/*           for advice on the choice of s see further comments. */
/*   nest  : integer. on entry nest must contain an over-estimate of the */
/*           total number of knots of the splines returned, to indicate */
/*           the storage space available to the routine. nest >=2*k+2. */
/*           in most practical situation nest=m/2 will be sufficient. */
/*           always large enough is nest=m+2*k, the number of knots */
/*           needed for interpolation (s=0). unchanged on exit. */
/*   n     : integer. */
/*           unless ier = 10 (in case iopt >=0), n will contain the */
/*           total number of knots of the smoothing spline curve returned */
/*           if the computation mode iopt=1 is used this value of n */
/*           should be left unchanged between subsequent calls. */
/*           in case iopt=-1, the value of n must be specified on entry. */
/*   t     : real array of dimension at least (nest). */
/*           on successful exit, this array will contain the knots of the */
/*           spline curve,i.e. the position of the interior knots t(k+2), */
/*           t(k+3),..,t(n-k-1) as well as the position of the additional */
/*           t(1),t(2),..,t(k+1)=u(1) and u(m)=t(n-k),...,t(n) needed for */
/*           the b-spline representation. */
/*           if the computation mode iopt=1 is used, the values of t(1), */
/*           t(2),...,t(n) should be left unchanged between subsequent */
/*           calls. if the computation mode iopt=-1 is used, the values */
/*           t(k+2),...,t(n-k-1) must be supplied by the user, before */
/*           entry. see also the restrictions (ier=10). */
/*   nc    : integer. on entry nc must specify the actual dimension of */
/*           the array c as declared in the calling (sub)program. nc */
/*           must not be too small (see c). unchanged on exit. */
/*   c     : real array of dimension at least (nest*idim). */
/*           on successful exit, this array will contain the coefficients */
/*           in the b-spline representation of the spline curve s(u),i.e. */
/*           the b-spline coefficients of the spline sj(u) will be given */
/*           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim. */
/*   fp    : real. unless ier = 10, fp contains the weighted sum of */
/*           squared residuals of the spline curve returned. */
/*   wrk   : real array of dimension at least m*(k+1)+nest*(7+idim+5*k). */
/*           used as working space. if the computation mode iopt=1 is */
/*           used, the values wrk(1),...,wrk(n) should be left unchanged */
/*           between subsequent calls. */
/*   lwrk  : integer. on entry,lwrk must specify the actual dimension of */
/*           the array wrk as declared in the calling (sub)program. lwrk */
/*           must not be too small (see wrk). unchanged on exit. */
/*   iwrk  : integer array of dimension at least (nest). */
/*           used as working space. if the computation mode iopt=1 is */
/*           used,the values iwrk(1),...,iwrk(n) should be left unchanged */
/*           between subsequent calls. */
/*   ier   : integer. unless the routine detects an error, ier contains a */
/*           non-positive value on exit, i.e. */
/*    ier=0  : normal return. the close curve returned has a residual */
/*             sum of squares fp such that abs(fp-s)/s <= tol with tol a */
/*             relative tolerance set to 0.001 by the program. */
/*    ier=-1 : normal return. the curve returned is an interpolating */
/*             spline curve (fp=0). */
/*    ier=-2 : normal return. the curve returned is the weighted least- */
/*             squares point,i.e. each spline sj(u) is a constant. in */
/*             this extreme case fp gives the upper bound fp0 for the */
/*             smoothing factor s. */
/*    ier=1  : error. the required storage space exceeds the available */
/*             storage space, as specified by the parameter nest. */
/*             probably causes : nest too small. if nest is already */
/*             large (say nest > m/2), it may also indicate that s is */
/*             too small */
/*             the approximation returned is the least-squares closed */
/*             curve according to the knots t(1),t(2),...,t(n). (n=nest) */
/*             the parameter fp gives the corresponding weighted sum of */
/*             squared residuals (fp>s). */
/*    ier=2  : error. a theoretically impossible result was found during */
/*             the iteration process for finding a smoothing curve with */
/*             fp = s. probably causes : s too small. */
/*             there is an approximation returned but the corresponding */
/*             weighted sum of squared residuals does not satisfy the */
/*             condition abs(fp-s)/s < tol. */
/*    ier=3  : error. the maximal number of iterations maxit (set to 20 */
/*             by the program) allowed for finding a smoothing curve */
/*             with fp=s has been reached. probably causes : s too small */
/*             there is an approximation returned but the corresponding */
/*             weighted sum of squared residuals does not satisfy the */
/*             condition abs(fp-s)/s < tol. */
/*    ier=10 : error. on entry, the input data are controlled on validity */
/*             the following restrictions must be satisfied. */
/*             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,2,...,m */
/*             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(7+idim+5*k), */
/*             nc>=nest*idim, x(j)=x(idim*(m-1)+j), j=1,2,...,idim */
/*             if ipar=0: sum j=1,idim (x(i*idim+j)-x((i-1)*idim+j))**2>0 */
/*                        i=1,2,...,m-1. */
/*             if ipar=1: u(1)<u(2)<...<u(m) */
/*             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k) */
/*                         u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m) */
/*                            (u(1)=0 and u(m)=1 in case ipar=0) */
/*                       the schoenberg-whitney conditions, i.e. there */
/*                       must be a subset of data points uu(j) with */
/*                       uu(j) = u(i) or u(i)+(u(m)-u(1)) such that */
/*                         t(j) < uu(j) < t(j+k+1), j=k+1,...,n-k-1 */
/*             if iopt>=0: s>=0 */
/*                         if s=0 : nest >= m+2*k */
/*             if one of these conditions is found to be violated,control */
/*             is immediately repassed to the calling program. in that */
/*             case there is no approximation returned. */

/*  further comments: */
/*   by means of the parameter s, the user can control the tradeoff */
/*   between closeness of fit and smoothness of fit of the approximation. */
/*   if s is too large, the curve will be too smooth and signal will be */
/*   lost ; if s is too small the curve will pick up too much noise. in */
/*   the extreme cases the program will return an interpolating curve if */
/*   s=0 and the weighted least-squares point if s is very large. */
/*   between these extremes, a properly chosen s will result in a good */
/*   compromise between closeness of fit and smoothness of fit. */
/*   to decide whether an approximation, corresponding to a certain s is */
/*   satisfactory the user is highly recommended to inspect the fits */
/*   graphically. */
/*   recommended values for s depend on the weights w(i). if these are */
/*   taken as 1/d(i) with d(i) an estimate of the standard deviation of */
/*   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+ */
/*   sqrt(2*m)). if nothing is known about the statistical error in x(i) */
/*   each w(i) can be set equal to one and s determined by trial and */
/*   error, taking account of the comments above. the best is then to */
/*   start with a very large value of s ( to determine the weighted */
/*   least-squares point and the upper bound fp0 for s) and then to */
/*   progressively decrease the value of s ( say by a factor 10 in the */
/*   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the */
/*   approximating curve shows more detail) to obtain closer fits. */
/*   to economize the search for a good s-value the program provides with */
/*   different modes of computation. at the first call of the routine, or */
/*   whenever he wants to restart with the initial set of knots the user */
/*   must set iopt=0. */
/*   if iopt=1 the program will continue with the set of knots found at */
/*   the last call of the routine. this will save a lot of computation */
/*   time if clocur is called repeatedly for different values of s. */
/*   the number of knots of the spline returned and their location will */
/*   depend on the value of s and on the complexity of the shape of the */
/*   curve underlying the data. but, if the computation mode iopt=1 is */
/*   used, the knots returned may also depend on the s-values at previous */
/*   calls (if these were smaller). therefore, if after a number of */
/*   trials with different s-values and iopt=1, the user can finally */
/*   accept a fit as satisfactory, it may be worthwhile for him to call */
/*   clocur once more with the selected value for s but now with iopt=0. */
/*   indeed, clocur may then return an approximation of the same quality */
/*   of fit but with fewer knots and therefore better if data reduction */
/*   is also an important objective for the user. */

/*   the form of the approximating curve can strongly be affected  by */
/*   the choice of the parameter values u(i). if there is no physical */
/*   reason for choosing a particular parameter u, often good results */
/*   will be obtained with the choice of clocur(in case ipar=0), i.e. */
/*        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m */
/*   where */
/*        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 ) */
/*   other possibilities for q(i) are */
/*        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2 */
/*        q(i)= sum j=1,idim abs(xj(i)-xj(i-1)) */
/*        q(i)= max j=1,idim abs(xj(i)-xj(i-1)) */
/*        q(i)= 1 */


/*  other subroutines required: */
/*    fpbacp,fpbspl,fpchep,fpclos,fpdisc,fpgivs,fpknot,fprati,fprota */

/*  references: */
/*   dierckx p. : algorithms for smoothing data with periodic and */
/*                parametric splines, computer graphics and image */
/*                processing 20 (1982) 171-184. */
/*   dierckx p. : algorithms for smoothing data with periodic and param- */
/*                etric splines, report tw55, dept. computer science, */
/*                k.u.leuven, 1981. */
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
/*  ..function references.. */
/*  we set up the parameters tol and maxit */
    /* Parameter adjustments */
    --w;
    --u;
    --x;
    --iwrk;
    --t;
    --c__;
    --wrk;

    /* Function Body */
    maxit = 20;
    tol = .001f;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (*iopt < -1 || *iopt > 1) {
	goto L90;
    }
    if (*ipar < 0 || *ipar > 1) {
	goto L90;
    }
    if (*idim <= 0 || *idim > 10) {
	goto L90;
    }
    if (*k <= 0 || *k > 5) {
	goto L90;
    }
    k1 = *k + 1;
    k2 = k1 + 1;
    nmin = k1 << 1;
    if (*m < 2 || *nest < nmin) {
	goto L90;
    }
    ncc = *nest * *idim;
    if (*mx < *m * *idim || *nc < ncc) {
	goto L90;
    }
    lwest = *m * k1 + *nest * (*idim + 7 + *k * 5);
    if (*lwrk < lwest) {
	goto L90;
    }
    i1 = *idim;
    i2 = *m * *idim;
    i__1 = *idim;
    for (j = 1; j <= i__1; ++j) {
	if (x[i1] != x[i2]) {
	    goto L90;
	}
	--i1;
	--i2;
/* L5: */
    }
    if (*ipar != 0 || *iopt > 0) {
	goto L40;
    }
    i1 = 0;
    i2 = *idim;
    u[1] = 0.f;
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dist = 0.f;
	i__2 = *idim;
	for (j1 = 1; j1 <= i__2; ++j1) {
	    ++i1;
	    ++i2;
/* Computing 2nd power */
	    d__1 = x[i2] - x[i1];
	    dist += d__1 * d__1;
/* L10: */
	}
	u[i__] = u[i__ - 1] + sqrt(dist);
/* L20: */
    }
    if (u[*m] <= 0.f) {
	goto L90;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	u[i__] /= u[*m];
/* L30: */
    }
    u[*m] = 1.f;
L40:
    if (w[1] <= 0.f) {
	goto L90;
    }
    m1 = *m - 1;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (u[i__] >= u[i__ + 1] || w[i__] <= 0.f) {
	    goto L90;
	}
/* L50: */
    }
    if (*iopt >= 0) {
	goto L70;
    }
    if (*n <= nmin || *n > *nest) {
	goto L90;
    }
    per = u[*m] - u[1];
    j1 = k1;
    t[j1] = u[1];
    i1 = *n - *k;
    t[i1] = u[*m];
    j2 = j1;
    i2 = i1;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++i1;
	--i2;
	++j1;
	--j2;
	t[j2] = t[i2] - per;
	t[i1] = t[j1] + per;
/* L60: */
    }
    fpchep_(&u[1], m, &t[1], n, k, ier);
    if (*ier == 0) {
	goto L80;
    }
    goto L90;
L70:
    if (*s < 0.f) {
	goto L90;
    }
    if (*s == 0.f && *nest < *m + (*k << 1)) {
	goto L90;
    }
    *ier = 0;
/* we partition the working space and determine the spline approximation. */
L80:
    ifp = 1;
    iz = ifp + *nest;
    ia1 = iz + ncc;
    ia2 = ia1 + *nest * k1;
    ib = ia2 + *nest * *k;
    ig1 = ib + *nest * k2;
    ig2 = ig1 + *nest * k2;
    iq = ig2 + *nest * k1;
    fpclos_(iopt, idim, m, &u[1], mx, &x[1], &w[1], k, s, nest, &tol, &maxit, 
	    &k1, &k2, n, &t[1], &ncc, &c__[1], fp, &wrk[ifp], &wrk[iz], &wrk[
	    ia1], &wrk[ia2], &wrk[ib], &wrk[ig1], &wrk[ig2], &wrk[iq], &iwrk[
	    1], ier);
L90:
    return 0;
} /* clocur_ */

