/* sphere.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int sphere_(integer *iopt, integer *m, doublereal *teta, 
	doublereal *phi, doublereal *r__, doublereal *w, doublereal *s, 
	integer *ntest, integer *npest, doublereal *eps, integer *nt, 
	doublereal *tt, integer *np, doublereal *tp, doublereal *c__, 
	doublereal *fp, doublereal *wrk1, integer *lwrk1, doublereal *wrk2, 
	integer *lwrk2, integer *iwrk, integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, j, la, lf, ki, lh;
    static doublereal pi;
    static integer kn, lq, ib1, ib3;
    static doublereal pi2;
    static integer np4, nt4, lcc, ncc, lff, lbp, lco;
    static doublereal one;
    static integer lbt, lcs, lfp, lro, npp;
    static doublereal tol;
    static integer lsp, lst, ntt, ncof, nreg, ncest, maxit, nrint, kwest, 
	    lwest;
    extern /* Subroutine */ int fpsphe_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *);

/*  subroutine sphere determines a smooth bicubic spherical spline */
/*  approximation s(teta,phi), 0 <= teta <= pi ; 0 <= phi <= 2*pi */
/*  to a given set of data points (teta(i),phi(i),r(i)),i=1,2,...,m. */
/*  such a spline has the following specific properties */

/*    (1) s(0,phi)  = constant   0 <=phi<= 2*pi. */

/*    (2) s(pi,phi) = constant   0 <=phi<= 2*pi */

/*         j             j */
/*        d s(teta,0)   d s(teta,2*pi) */
/*    (3) ----------- = ------------   0 <=teta<=pi, j=0,1,2 */
/*             j             j */
/*        d phi         d phi */

/*        d s(0,phi)    d s(0,0)             d s(0,pi/2) */
/*    (4) ----------  = -------- *cos(phi) + ----------- *sin(phi) */
/*        d teta        d teta               d teta */

/*        d s(pi,phi)   d s(pi,0)            d s(pi,pi/2) */
/*    (5) ----------- = ---------*cos(phi) + ------------*sin(phi) */
/*        d teta        d teta               d teta */

/*  if iopt =-1 sphere calculates a weighted least-squares spherical */
/*  spline according to a given set of knots in teta- and phi- direction. */
/*  if iopt >=0, the number of knots in each direction and their position */
/*  tt(j),j=1,2,...,nt ; tp(j),j=1,2,...,np are chosen automatically by */
/*  the routine. the smoothness of s(teta,phi) is then achieved by mini- */
/*  malizing the discontinuity jumps of the derivatives of the spline */
/*  at the knots. the amount of smoothness of s(teta,phi) is determined */
/*  by the condition that fp = sum((w(i)*(r(i)-s(teta(i),phi(i))))**2) */
/*  be <= s, with s a given non-negative constant. */
/*  the spherical spline is given in the standard b-spline representation */
/*  of bicubic splines and can be evaluated by means of subroutine bispev */

/* calling sequence: */
/*     call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps, */
/*    *  nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier) */

/* parameters: */
/*  iopt  : integer flag. on entry iopt must specify whether a weighted */
/*          least-squares spherical spline (iopt=-1) or a smoothing */
/*          spherical spline (iopt=0 or 1) must be determined. */
/*          if iopt=0 the routine will start with an initial set of knots */
/*          tt(i)=0,tt(i+4)=pi,i=1,...,4;tp(i)=0,tp(i+4)=2*pi,i=1,...,4. */
/*          if iopt=1 the routine will continue with the set of knots */
/*          found at the last call of the routine. */
/*          attention: a call with iopt=1 must always be immediately pre- */
/*                     ceded by another call with iopt=1 or iopt=0. */
/*          unchanged on exit. */
/*  m     : integer. on entry m must specify the number of data points. */
/*          m >= 2. unchanged on exit. */
/*  teta  : real array of dimension at least (m). */
/*  phi   : real array of dimension at least (m). */
/*  r     : real array of dimension at least (m). */
/*          before entry,teta(i),phi(i),r(i) must be set to the spherical */
/*          co-ordinates of the i-th data point, for i=1,...,m.the order */
/*          of the data points is immaterial. unchanged on exit. */
/*  w     : real array of dimension at least (m). before entry, w(i) must */
/*          be set to the i-th value in the set of weights. the w(i) must */
/*          be strictly positive. unchanged on exit. */
/*  s     : real. on entry (in case iopt>=0) s must specify the smoothing */
/*          factor. s >=0. unchanged on exit. */
/*          for advice on the choice of s see further comments */
/*  ntest : integer. unchanged on exit. */
/*  npest : integer. unchanged on exit. */
/*          on entry, ntest and npest must specify an upper bound for the */
/*          number of knots required in the teta- and phi-directions. */
/*          these numbers will also determine the storage space needed by */
/*          the routine. ntest >= 8, npest >= 8. */
/*          in most practical situation ntest = npest = 8+sqrt(m/2) will */
/*          be sufficient. see also further comments. */
/*  eps   : real. */
/*          on entry, eps must specify a threshold for determining the */
/*          effective rank of an over-determined linear system of equat- */
/*          ions. 0 < eps < 1.  if the number of decimal digits in the */
/*          computer representation of a real number is q, then 10**(-q) */
/*          is a suitable value for eps in most practical applications. */
/*          unchanged on exit. */
/*  nt    : integer. */
/*          unless ier=10 (in case iopt >=0), nt will contain the total */
/*          number of knots with respect to the teta-variable, of the */
/*          spline approximation returned. if the computation mode iopt=1 */
/*          is used, the value of nt should be left unchanged between */
/*          subsequent calls. */
/*          in case iopt=-1, the value of nt should be specified on entry */
/*  tt    : real array of dimension at least ntest. */
/*          on successful exit, this array will contain the knots of the */
/*          spline with respect to the teta-variable, i.e. the position */
/*          of the interior knots tt(5),...,tt(nt-4) as well as the */
/*          position of the additional knots tt(1)=...=tt(4)=0 and */
/*          tt(nt-3)=...=tt(nt)=pi needed for the b-spline representation */
/*          if the computation mode iopt=1 is used, the values of tt(1), */
/*          ...,tt(nt) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt=-1 is used, the values tt(5), */
/*          ...tt(nt-4) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  np    : integer. */
/*          unless ier=10 (in case iopt >=0), np will contain the total */
/*          number of knots with respect to the phi-variable, of the */
/*          spline approximation returned. if the computation mode iopt=1 */
/*          is used, the value of np should be left unchanged between */
/*          subsequent calls. */
/*          in case iopt=-1, the value of np (>=9) should be specified */
/*          on entry. */
/*  tp    : real array of dimension at least npest. */
/*          on successful exit, this array will contain the knots of the */
/*          spline with respect to the phi-variable, i.e. the position of */
/*          the interior knots tp(5),...,tp(np-4) as well as the position */
/*          of the additional knots tp(1),...,tp(4) and tp(np-3),..., */
/*          tp(np) needed for the b-spline representation. */
/*          if the computation mode iopt=1 is used, the values of tp(1), */
/*          ...,tp(np) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt=-1 is used, the values tp(5), */
/*          ...tp(np-4) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  c     : real array of dimension at least (ntest-4)*(npest-4). */
/*          on successful exit, c contains the coefficients of the spline */
/*          approximation s(teta,phi). */
/*  fp    : real. unless ier=10, fp contains the weighted sum of */
/*          squared residuals of the spline approximation returned. */
/*  wrk1  : real array of dimension (lwrk1). used as workspace. */
/*          if the computation mode iopt=1 is used the value of wrk1(1) */
/*          should be left unchanged between subsequent calls. */
/*          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the */
/*          values d(i)/max(d(i)),i=1,...,ncof=6+(np-7)*(nt-8) */
/*          with d(i) the i-th diagonal element of the reduced triangular */
/*          matrix for calculating the b-spline coefficients. it includes */
/*          those elements whose square is less than eps,which are treat- */
/*          ed as 0 in the case of presumed rank deficiency (ier<-2). */
/*  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of */
/*          the array wrk1 as declared in the calling (sub)program. */
/*          lwrk1 must not be too small. let */
/*            u = ntest-7, v = npest-7, then */
/*          lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m */
/*  wrk2  : real array of dimension (lwrk2). used as workspace, but */
/*          only in the case a rank deficient system is encountered. */
/*  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of */
/*          the array wrk2 as declared in the calling (sub)program. */
/*          lwrk2 > 0 . a save upper bound  for lwrk2 = 48+21*v+7*u*v+ */
/*          4*(u-1)*v**2 where u,v are as above. if there are enough data */
/*          points, scattered uniformly over the approximation domain */
/*          and if the smoothing factor s is not too small, there is a */
/*          good chance that this extra workspace is not needed. a lot */
/*          of memory might therefore be saved by setting lwrk2=1. */
/*          (see also ier > 10) */
/*  iwrk  : integer array of dimension (kwrk). used as workspace. */
/*  kwrk  : integer. on entry kwrk must specify the actual dimension of */
/*          the array iwrk as declared in the calling (sub)program. */
/*          kwrk >= m+(ntest-7)*(npest-7). */
/*  ier   : integer. unless the routine detects an error, ier contains a */
/*          non-positive value on exit, i.e. */
/*   ier=0  : normal return. the spline returned has a residual sum of */
/*            squares fp such that abs(fp-s)/s <= tol with tol a relat- */
/*            ive tolerance set to 0.001 by the program. */
/*   ier=-1 : normal return. the spline returned is a spherical */
/*            interpolating spline (fp=0). */
/*   ier=-2 : normal return. the spline returned is the weighted least- */
/*            squares constrained polynomial . in this extreme case */
/*            fp gives the upper bound for the smoothing factor s. */
/*   ier<-2 : warning. the coefficients of the spline returned have been */
/*            computed as the minimal norm least-squares solution of a */
/*            (numerically) rank deficient system. (-ier) gives the rank. */
/*            especially if the rank deficiency which can be computed as */
/*            6+(nt-8)*(np-7)+ier, is large the results may be inaccurate */
/*            they could also seriously depend on the value of eps. */
/*   ier=1  : error. the required storage space exceeds the available */
/*            storage space, as specified by the parameters ntest and */
/*            npest. */
/*            probably causes : ntest or npest too small. if these param- */
/*            eters are already large, it may also indicate that s is */
/*            too small */
/*            the approximation returned is the weighted least-squares */
/*            spherical spline according to the current set of knots. */
/*            the parameter fp gives the corresponding weighted sum of */
/*            squared residuals (fp>s). */
/*   ier=2  : error. a theoretically impossible result was found during */
/*            the iteration process for finding a smoothing spline with */
/*            fp = s. probably causes : s too small or badly chosen eps. */
/*            there is an approximation returned but the corresponding */
/*            weighted sum of squared residuals does not satisfy the */
/*            condition abs(fp-s)/s < tol. */
/*   ier=3  : error. the maximal number of iterations maxit (set to 20 */
/*            by the program) allowed for finding a smoothing spline */
/*            with fp=s has been reached. probably causes : s too small */
/*            there is an approximation returned but the corresponding */
/*            weighted sum of squared residuals does not satisfy the */
/*            condition abs(fp-s)/s < tol. */
/*   ier=4  : error. no more knots can be added because the dimension */
/*            of the spherical spline 6+(nt-8)*(np-7) already exceeds */
/*            the number of data points m. */
/*            probably causes : either s or m too small. */
/*            the approximation returned is the weighted least-squares */
/*            spherical spline according to the current set of knots. */
/*            the parameter fp gives the corresponding weighted sum of */
/*            squared residuals (fp>s). */
/*   ier=5  : error. no more knots can be added because the additional */
/*            knot would (quasi) coincide with an old one. */
/*            probably causes : s too small or too large a weight to an */
/*            inaccurate data point. */
/*            the approximation returned is the weighted least-squares */
/*            spherical spline according to the current set of knots. */
/*            the parameter fp gives the corresponding weighted sum of */
/*            squared residuals (fp>s). */
/*   ier=10 : error. on entry, the input data are controlled on validity */
/*            the following restrictions must be satisfied. */
/*            -1<=iopt<=1,  m>=2, ntest>=8 ,npest >=8, 0<eps<1, */
/*            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m */
/*            lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m */
/*            kwrk >= m+(ntest-7)*(npest-7) */
/*            if iopt=-1: 8<=nt<=ntest , 9<=np<=npest */
/*                        0<tt(5)<tt(6)<...<tt(nt-4)<pi */
/*                        0<tp(5)<tp(6)<...<tp(np-4)<2*pi */
/*            if iopt>=0: s>=0 */
/*            if one of these conditions is found to be violated,control */
/*            is immediately repassed to the calling program. in that */
/*            case there is no approximation returned. */
/*   ier>10 : error. lwrk2 is too small, i.e. there is not enough work- */
/*            space for computing the minimal least-squares solution of */
/*            a rank deficient system of linear equations. ier gives the */
/*            requested value for lwrk2. there is no approximation re- */
/*            turned but, having saved the information contained in nt, */
/*            np,tt,tp,wrk1, and having adjusted the value of lwrk2 and */
/*            the dimension of the array wrk2 accordingly, the user can */
/*            continue at the point the program was left, by calling */
/*            sphere with iopt=1. */

/* further comments: */
/*  by means of the parameter s, the user can control the tradeoff */
/*   between closeness of fit and smoothness of fit of the approximation. */
/*   if s is too large, the spline will be too smooth and signal will be */
/*   lost ; if s is too small the spline will pick up too much noise. in */
/*   the extreme cases the program will return an interpolating spline if */
/*   s=0 and the constrained weighted least-squares polynomial if s is */
/*   very large. between these extremes, a properly chosen s will result */
/*   in a good compromise between closeness of fit and smoothness of fit. */
/*   to decide whether an approximation, corresponding to a certain s is */
/*   satisfactory the user is highly recommended to inspect the fits */
/*   graphically. */
/*   recommended values for s depend on the weights w(i). if these are */
/*   taken as 1/d(i) with d(i) an estimate of the standard deviation of */
/*   r(i), a good s-value should be found in the range (m-sqrt(2*m),m+ */
/*   sqrt(2*m)). if nothing is known about the statistical error in r(i) */
/*   each w(i) can be set equal to one and s determined by trial and */
/*   error, taking account of the comments above. the best is then to */
/*   start with a very large value of s ( to determine the least-squares */
/*   polynomial and the corresponding upper bound fp0 for s) and then to */
/*   progressively decrease the value of s ( say by a factor 10 in the */
/*   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the */
/*   approximation shows more detail) to obtain closer fits. */
/*   to choose s very small is strongly discouraged. this considerably */
/*   increases computation time and memory requirements. it may also */
/*   cause rank-deficiency (ier<-2) and endager numerical stability. */
/*   to economize the search for a good s-value the program provides with */
/*   different modes of computation. at the first call of the routine, or */
/*   whenever he wants to restart with the initial set of knots the user */
/*   must set iopt=0. */
/*   if iopt=1 the program will continue with the set of knots found at */
/*   the last call of the routine. this will save a lot of computation */
/*   time if sphere is called repeatedly for different values of s. */
/*   the number of knots of the spline returned and their location will */
/*   depend on the value of s and on the complexity of the shape of the */
/*   function underlying the data. if the computation mode iopt=1 */
/*   is used, the knots returned may also depend on the s-values at */
/*   previous calls (if these were smaller). therefore, if after a number */
/*   of trials with different s-values and iopt=1, the user can finally */
/*   accept a fit as satisfactory, it may be worthwhile for him to call */
/*   sphere once more with the selected value for s but now with iopt=0. */
/*   indeed, sphere may then return an approximation of the same quality */
/*   of fit but with fewer knots and therefore better if data reduction */
/*   is also an important objective for the user. */
/*   the number of knots may also depend on the upper bounds ntest and */
/*   npest. indeed, if at a certain stage in sphere the number of knots */
/*   in one direction (say nt) has reached the value of its upper bound */
/*   (ntest), then from that moment on all subsequent knots are added */
/*   in the other (phi) direction. this may indicate that the value of */
/*   ntest is too small. on the other hand, it gives the user the option */
/*   of limiting the number of knots the routine locates in any direction */
/*   for example, by setting ntest=8 (the lowest allowable value for */
/*   ntest), the user can indicate that he wants an approximation which */
/*   is a cubic polynomial in the variable teta. */

/*  other subroutines required: */
/*    fpback,fpbspl,fpsphe,fpdisc,fpgivs,fprank,fprati,fprota,fporde, */
/*    fprpsp */

/*  references: */
/*   dierckx p. : algorithms for smoothing data on the sphere with tensor */
/*                product splines, computing 32 (1984) 319-342. */
/*   dierckx p. : algorithms for smoothing data on the sphere with tensor */
/*                product splines, report tw62, dept. computer science, */
/*                k.u.leuven, 1983. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author: */
/*    p.dierckx */
/*    dept. computer science, k.u. leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  creation date : july 1983 */
/*  latest update : march 1989 */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpsphe */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --w;
    --r__;
    --phi;
    --teta;
    --tt;
    --c__;
    --tp;
    --wrk1;
    --wrk2;
    --iwrk;

    /* Function Body */
    one = 1.f;
/*  we set up the parameters tol and maxit. */
    maxit = 20;
    tol = .001f;
/*  before starting computations a data check is made. if the input data */
/*  are invalid,control is immediately repassed to the calling program. */
    *ier = 10;
    if (*eps <= 0.f || *eps >= 1.f) {
	goto L80;
    }
    if (*iopt < -1 || *iopt > 1) {
	goto L80;
    }
    if (*m < 2) {
	goto L80;
    }
    if (*ntest < 8 || *npest < 8) {
	goto L80;
    }
    nt4 = *ntest - 4;
    np4 = *npest - 4;
    ncest = nt4 * np4;
    ntt = *ntest - 7;
    npp = *npest - 7;
    ncc = npp * (ntt - 1) + 6;
    nrint = ntt + npp;
    nreg = ntt * npp;
    ncof = npp * 3 + 6;
    ib1 = npp << 2;
    ib3 = ib1 + 3;
    if (ncof > ib1) {
	ib1 = ncof;
    }
    if (ncof > ib3) {
	ib3 = ncof;
    }
/* Computing 2nd power */
    i__1 = npp;
    lwest = npp * 52 + 185 + ntt * 10 + ntt * 14 * npp + (*m + (ntt - 1) * (
	    i__1 * i__1) << 3);
    kwest = *m + nreg;
    if (*lwrk1 < lwest || *kwrk < kwest) {
	goto L80;
    }
    if (*iopt > 0) {
	goto L60;
    }
    pi = atan(one) * 4;
    pi2 = pi + pi;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (w[i__] <= 0.f) {
	    goto L80;
	}
	if (teta[i__] < 0.f || teta[i__] > pi) {
	    goto L80;
	}
	if (phi[i__] < 0.f || phi[i__] > pi2) {
	    goto L80;
	}
/* L20: */
    }
    if (*iopt == 0) {
	goto L60;
    }
    ntt = *nt - 8;
    if (ntt < 0 || *nt > *ntest) {
	goto L80;
    }
    if (ntt == 0) {
	goto L40;
    }
    tt[4] = 0.f;
    i__1 = ntt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__ + 4;
	if (tt[j] <= tt[j - 1] || tt[j] >= pi) {
	    goto L80;
	}
/* L30: */
    }
L40:
    npp = *np - 8;
    if (npp < 1 || *np > *npest) {
	goto L80;
    }
    tp[4] = 0.f;
    i__1 = npp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__ + 4;
	if (tp[j] <= tp[j - 1] || tp[j] >= pi2) {
	    goto L80;
	}
/* L50: */
    }
    goto L70;
L60:
    if (*s < 0.f) {
	goto L80;
    }
L70:
    *ier = 0;
/*  we partition the working space and determine the spline approximation */
    kn = 1;
    ki = kn + *m;
    lq = 2;
    la = lq + ncc * ib3;
    lf = la + ncc * ib1;
    lff = lf + ncc;
    lfp = lff + ncest;
    lco = lfp + nrint;
    lh = lco + nrint;
    lbt = lh + ib3;
    lbp = lbt + *ntest * 5;
    lro = lbp + *npest * 5;
    lcc = lro + *npest;
    lcs = lcc + *npest;
    lst = lcs + *npest;
    lsp = lst + (*m << 2);
    fpsphe_(iopt, m, &teta[1], &phi[1], &r__[1], &w[1], s, ntest, npest, eps, 
	    &tol, &maxit, &ib1, &ib3, &ncest, &ncc, &nrint, &nreg, nt, &tt[1],
	     np, &tp[1], &c__[1], fp, &wrk1[1], &wrk1[lfp], &wrk1[lco], &wrk1[
	    lf], &wrk1[lff], &wrk1[lro], &wrk1[lcc], &wrk1[lcs], &wrk1[la], &
	    wrk1[lq], &wrk1[lbt], &wrk1[lbp], &wrk1[lst], &wrk1[lsp], &wrk1[
	    lh], &iwrk[ki], &iwrk[kn], &wrk2[1], lwrk2, ier);
L80:
    return 0;
} /* sphere_ */

