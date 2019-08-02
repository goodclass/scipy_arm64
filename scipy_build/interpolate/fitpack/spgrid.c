/* spgrid.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int spgrid_(integer *iopt, integer *ider, integer *mu, 
	doublereal *u, integer *mv, doublereal *v, doublereal *r__, 
	doublereal *r0, doublereal *r1, doublereal *s, integer *nuest, 
	integer *nvest, integer *nu, doublereal *tu, integer *nv, doublereal *
	tv, doublereal *c__, doublereal *fp, doublereal *wrk, integer *lwrk, 
	integer *iwrk, integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double atan2(doublereal, doublereal);

    /* Local variables */
    static integer i__, j, l, m, i1, i2, j1, j2, nc;
    static doublereal rb, re, pi, ve, rn, uu, one;
    static integer ldr;
    static doublereal per, tol;
    static integer muu, lww;
    static doublereal half;
    static integer kndu, kndv;
    static doublereal rmin;
    static integer lfpu;
    static doublereal rmax;
    static integer lfpv, jwrk, knru, knrv, maxit, mumin, kwest, lwest;
    extern /* Subroutine */ int fpchec_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *), fpchep_(doublereal *, integer *
	    , doublereal *, integer *, integer *, integer *), fpspgr_(integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *);

/*  given the function values r(i,j) on the latitude-longitude grid */
/*  (u(i),v(j)), i=1,...,mu ; j=1,...,mv , spgrid determines a smooth */
/*  bicubic spline approximation on the rectangular domain 0<=u<=pi, */
/*  vb<=v<=ve (vb = v(1), ve=vb+2*pi). */
/*  this approximation s(u,v) will satisfy the properties */

/*    (1) s(0,v) = s(0,0) = dr(1) */

/*        d s(0,v)           d s(0,0)           d s(0,pi/2) */
/*    (2) -------- = cos(v)* -------- + sin(v)* ----------- */
/*        d u                d u                d u */

/*                 = cos(v)*dr(2)+sin(v)*dr(3) */
/*                                                     vb <= v <= ve */
/*    (3) s(pi,v) = s(pi,0) = dr(4) */

/*        d s(pi,v)           d s(pi,0)           d s(pi,pi/2) */
/*    (4) -------- = cos(v)*  --------- + sin(v)* ------------ */
/*        d u                 d u                 d u */

/*                 = cos(v)*dr(5)+sin(v)*dr(6) */

/*  and will be periodic in the variable v, i.e. */

/*         j           j */
/*        d s(u,vb)   d s(u,ve) */
/*    (5) --------- = ---------   0 <=u<= pi , j=0,1,2 */
/*           j           j */
/*        d v         d v */

/*  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu; */
/*  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the */
/*  smoothness of s(u,v) is achieved by minimalizing the discontinuity */
/*  jumps of the derivatives of the spline at the knots. the amount of */
/*  smoothness of s(u,v) is determined by the condition that */
/*  fp=sumi=1,mu(sumj=1,mv((r(i,j)-s(u(i),v(j)))**2))+(r0-s(0,v))**2 */
/*  + (r1-s(pi,v))**2 <= s, with s a given non-negative constant. */
/*  the fit s(u,v) is given in its b-spline representation and can be */
/*  evaluated by means of routine bispev */

/* calling sequence: */
/*     call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu, */
/*    *  ,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier) */

/* parameters: */
/*  iopt  : integer array of dimension 3, specifying different options. */
/*          unchanged on exit. */
/*  iopt(1):on entry iopt(1) must specify whether a least-squares spline */
/*          (iopt(1)=-1) or a smoothing spline (iopt(1)=0 or 1) must be */
/*          determined. */
/*          if iopt(1)=0 the routine will start with an initial set of */
/*          knots tu(i)=0,tu(i+4)=pi,i=1,...,4;tv(i)=v(1)+(i-4)*2*pi, */
/*          i=1,...,8. */
/*          if iopt(1)=1 the routine will continue with the set of knots */
/*          found at the last call of the routine. */
/*          attention: a call with iopt(1)=1 must always be immediately */
/*          preceded by another call with iopt(1) = 1 or iopt(1) = 0. */
/*  iopt(2):on entry iopt(2) must specify the requested order of conti- */
/*          nuity at the pole u=0. */
/*          if iopt(2)=0 only condition (1) must be fulfilled and */
/*          if iopt(2)=1 conditions (1)+(2) must be fulfilled. */
/*  iopt(3):on entry iopt(3) must specify the requested order of conti- */
/*          nuity at the pole u=pi. */
/*          if iopt(3)=0 only condition (3) must be fulfilled and */
/*          if iopt(3)=1 conditions (3)+(4) must be fulfilled. */
/*  ider  : integer array of dimension 4, specifying different options. */
/*          unchanged on exit. */
/*  ider(1):on entry ider(1) must specify whether (ider(1)=0 or 1) or not */
/*          (ider(1)=-1) there is a data value r0 at the pole u=0. */
/*          if ider(1)=1, r0 will be considered to be the right function */
/*          value, and it will be fitted exactly (s(0,v)=r0). */
/*          if ider(1)=0, r0 will be considered to be a data value just */
/*          like the other data values r(i,j). */
/*  ider(2):on entry ider(2) must specify whether (ider(2)=1) or not */
/*          (ider(2)=0) the approximation has vanishing derivatives */
/*          dr(2) and dr(3) at the pole u=0  (in case iopt(2)=1) */
/*  ider(3):on entry ider(3) must specify whether (ider(3)=0 or 1) or not */
/*          (ider(3)=-1) there is a data value r1 at the pole u=pi. */
/*          if ider(3)=1, r1 will be considered to be the right function */
/*          value, and it will be fitted exactly (s(pi,v)=r1). */
/*          if ider(3)=0, r1 will be considered to be a data value just */
/*          like the other data values r(i,j). */
/*  ider(4):on entry ider(4) must specify whether (ider(4)=1) or not */
/*          (ider(4)=0) the approximation has vanishing derivatives */
/*          dr(5) and dr(6) at the pole u=pi (in case iopt(3)=1) */
/*  mu    : integer. on entry mu must specify the number of grid points */
/*          along the u-axis. unchanged on exit. */
/*          mu >= 1, mu >=mumin=4-i0-i1-ider(2)-ider(4) with */
/*            i0=min(1,ider(1)+1), i1=min(1,ider(3)+1) */
/*  u     : real array of dimension at least (mu). before entry, u(i) */
/*          must be set to the u-co-ordinate of the i-th grid point */
/*          along the u-axis, for i=1,2,...,mu. these values must be */
/*          supplied in strictly ascending order. unchanged on exit. */
/*          0 < u(i) < pi. */
/*  mv    : integer. on entry mv must specify the number of grid points */
/*          along the v-axis. mv > 3 . unchanged on exit. */
/*  v     : real array of dimension at least (mv). before entry, v(j) */
/*          must be set to the v-co-ordinate of the j-th grid point */
/*          along the v-axis, for j=1,2,...,mv. these values must be */
/*          supplied in strictly ascending order. unchanged on exit. */
/*          -pi <= v(1) < pi , v(mv) < v(1)+2*pi. */
/*  r     : real array of dimension at least (mu*mv). */
/*          before entry, r(mv*(i-1)+j) must be set to the data value at */
/*          the grid point (u(i),v(j)) for i=1,...,mu and j=1,...,mv. */
/*          unchanged on exit. */
/*  r0    : real value. on entry (if ider(1) >=0 ) r0 must specify the */
/*          data value at the pole u=0. unchanged on exit. */
/*  r1    : real value. on entry (if ider(1) >=0 ) r1 must specify the */
/*          data value at the pole u=pi. unchanged on exit. */
/*  s     : real. on entry (if iopt(1)>=0) s must specify the smoothing */
/*          factor. s >=0. unchanged on exit. */
/*          for advice on the choice of s see further comments */
/*  nuest : integer. unchanged on exit. */
/*  nvest : integer. unchanged on exit. */
/*          on entry, nuest and nvest must specify an upper bound for the */
/*          number of knots required in the u- and v-directions respect. */
/*          these numbers will also determine the storage space needed by */
/*          the routine. nuest >= 8, nvest >= 8. */
/*          in most practical situation nuest = mu/2, nvest=mv/2, will */
/*          be sufficient. always large enough are nuest=mu+6+iopt(2)+ */
/*          iopt(3), nvest = mv+7, the number of knots needed for */
/*          interpolation (s=0). see also further comments. */
/*  nu    : integer. */
/*          unless ier=10 (in case iopt(1)>=0), nu will contain the total */
/*          number of knots with respect to the u-variable, of the spline */
/*          approximation returned. if the computation mode iopt(1)=1 is */
/*          used, the value of nu should be left unchanged between sub- */
/*          sequent calls. in case iopt(1)=-1, the value of nu should be */
/*          specified on entry. */
/*  tu    : real array of dimension at least (nuest). */
/*          on successful exit, this array will contain the knots of the */
/*          spline with respect to the u-variable, i.e. the position of */
/*          the interior knots tu(5),...,tu(nu-4) as well as the position */
/*          of the additional knots tu(1)=...=tu(4)=0 and tu(nu-3)=...= */
/*          tu(nu)=pi needed for the b-spline representation. */
/*          if the computation mode iopt(1)=1 is used,the values of tu(1) */
/*          ...,tu(nu) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt(1)=-1 is used, the values tu(5), */
/*          ...tu(nu-4) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  nv    : integer. */
/*          unless ier=10 (in case iopt(1)>=0), nv will contain the total */
/*          number of knots with respect to the v-variable, of the spline */
/*          approximation returned. if the computation mode iopt(1)=1 is */
/*          used, the value of nv should be left unchanged between sub- */
/*          sequent calls. in case iopt(1) = -1, the value of nv should */
/*          be specified on entry. */
/*  tv    : real array of dimension at least (nvest). */
/*          on successful exit, this array will contain the knots of the */
/*          spline with respect to the v-variable, i.e. the position of */
/*          the interior knots tv(5),...,tv(nv-4) as well as the position */
/*          of the additional knots tv(1),...,tv(4) and tv(nv-3),..., */
/*          tv(nv) needed for the b-spline representation. */
/*          if the computation mode iopt(1)=1 is used,the values of tv(1) */
/*          ...,tv(nv) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt(1)=-1 is used, the values tv(5), */
/*          ...tv(nv-4) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  c     : real array of dimension at least (nuest-4)*(nvest-4). */
/*          on successful exit, c contains the coefficients of the spline */
/*          approximation s(u,v) */
/*  fp    : real. unless ier=10, fp contains the sum of squared */
/*          residuals of the spline approximation returned. */
/*  wrk   : real array of dimension (lwrk). used as workspace. */
/*          if the computation mode iopt(1)=1 is used the values of */
/*          wrk(1),..,wrk(12) should be left unchanged between subsequent */
/*          calls. */
/*  lwrk  : integer. on entry lwrk must specify the actual dimension of */
/*          the array wrk as declared in the calling (sub)program. */
/*          lwrk must not be too small. */
/*           lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+q */
/*           where q is the larger of (mv+nvest) and nuest. */
/*  iwrk  : integer array of dimension (kwrk). used as workspace. */
/*          if the computation mode iopt(1)=1 is used the values of */
/*          iwrk(1),.,iwrk(5) should be left unchanged between subsequent */
/*          calls. */
/*  kwrk  : integer. on entry kwrk must specify the actual dimension of */
/*          the array iwrk as declared in the calling (sub)program. */
/*          kwrk >= 5+mu+mv+nuest+nvest. */
/*  ier   : integer. unless the routine detects an error, ier contains a */
/*          non-positive value on exit, i.e. */
/*   ier=0  : normal return. the spline returned has a residual sum of */
/*            squares fp such that abs(fp-s)/s <= tol with tol a relat- */
/*            ive tolerance set to 0.001 by the program. */
/*   ier=-1 : normal return. the spline returned is an interpolating */
/*            spline (fp=0). */
/*   ier=-2 : normal return. the spline returned is the least-squares */
/*            constrained polynomial. in this extreme case fp gives the */
/*            upper bound for the smoothing factor s. */
/*   ier=1  : error. the required storage space exceeds the available */
/*            storage space, as specified by the parameters nuest and */
/*            nvest. */
/*            probably causes : nuest or nvest too small. if these param- */
/*            eters are already large, it may also indicate that s is */
/*            too small */
/*            the approximation returned is the least-squares spline */
/*            according to the current set of knots. the parameter fp */
/*            gives the corresponding sum of squared residuals (fp>s). */
/*   ier=2  : error. a theoretically impossible result was found during */
/*            the iteration process for finding a smoothing spline with */
/*            fp = s. probably causes : s too small. */
/*            there is an approximation returned but the corresponding */
/*            sum of squared residuals does not satisfy the condition */
/*            abs(fp-s)/s < tol. */
/*   ier=3  : error. the maximal number of iterations maxit (set to 20 */
/*            by the program) allowed for finding a smoothing spline */
/*            with fp=s has been reached. probably causes : s too small */
/*            there is an approximation returned but the corresponding */
/*            sum of squared residuals does not satisfy the condition */
/*            abs(fp-s)/s < tol. */
/*   ier=10 : error. on entry, the input data are controlled on validity */
/*            the following restrictions must be satisfied. */
/*            -1<=iopt(1)<=1, 0<=iopt(2)<=1, 0<=iopt(3)<=1, */
/*            -1<=ider(1)<=1, 0<=ider(2)<=1, ider(2)=0 if iopt(2)=0. */
/*            -1<=ider(3)<=1, 0<=ider(4)<=1, ider(4)=0 if iopt(3)=0. */
/*            mu >= mumin (see above), mv >= 4, nuest >=8, nvest >= 8, */
/*            kwrk>=5+mu+mv+nuest+nvest, */
/*            lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+ */
/*             max(nuest,mv+nvest) */
/*            0< u(i-1)<u(i)< pi,i=2,..,mu, */
/*            -pi<=v(1)< pi, v(1)<v(i-1)<v(i)<v(1)+2*pi, i=3,...,mv */
/*            if iopt(1)=-1: 8<=nu<=min(nuest,mu+6+iopt(2)+iopt(3)) */
/*                           0<tu(5)<tu(6)<...<tu(nu-4)< pi */
/*                           8<=nv<=min(nvest,mv+7) */
/*                           v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(1)+2*pi */
/*                    the schoenberg-whitney conditions, i.e. there must */
/*                    be subset of grid co-ordinates uu(p) and vv(q) such */
/*                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4 */
/*                     (iopt(2)=1 and iopt(3)=1 also count for a uu-value */
/*                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4 */
/*                     (vv(q) is either a value v(j) or v(j)+2*pi) */
/*            if iopt(1)>=0: s>=0 */
/*                       if s=0: nuest>=mu+6+iopt(2)+iopt(3), nvest>=mv+7 */
/*            if one of these conditions is found to be violated,control */
/*            is immediately repassed to the calling program. in that */
/*            case there is no approximation returned. */

/* further comments: */
/*   spgrid does not allow individual weighting of the data-values. */
/*   so, if these were determined to widely different accuracies, then */
/*   perhaps the general data set routine sphere should rather be used */
/*   in spite of efficiency. */
/*   by means of the parameter s, the user can control the tradeoff */
/*   between closeness of fit and smoothness of fit of the approximation. */
/*   if s is too large, the spline will be too smooth and signal will be */
/*   lost ; if s is too small the spline will pick up too much noise. in */
/*   the extreme cases the program will return an interpolating spline if */
/*   s=0 and the constrained least-squares polynomial(degrees 3,0)if s is */
/*   very large. between these extremes, a properly chosen s will result */
/*   in a good compromise between closeness of fit and smoothness of fit. */
/*   to decide whether an approximation, corresponding to a certain s is */
/*   satisfactory the user is highly recommended to inspect the fits */
/*   graphically. */
/*   recommended values for s depend on the accuracy of the data values. */
/*   if the user has an idea of the statistical errors on the data, he */
/*   can also find a proper estimate for s. for, by assuming that, if he */
/*   specifies the right s, spgrid will return a spline s(u,v) which */
/*   exactly reproduces the function underlying the data he can evaluate */
/*   the sum((r(i,j)-s(u(i),v(j)))**2) to find a good estimate for this s */
/*   for example, if he knows that the statistical errors on his r(i,j)- */
/*   values is not greater than 0.1, he may expect that a good s should */
/*   have a value not larger than mu*mv*(0.1)**2. */
/*   if nothing is known about the statistical error in r(i,j), s must */
/*   be determined by trial and error, taking account of the comments */
/*   above. the best is then to start with a very large value of s (to */
/*   determine the least-squares polynomial and the corresponding upper */
/*   bound fp0 for s) and then to progressively decrease the value of s */
/*   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,... */
/*   and more carefully as the approximation shows more detail) to */
/*   obtain closer fits. */
/*   to economize the search for a good s-value the program provides with */
/*   different modes of computation. at the first call of the routine, or */
/*   whenever he wants to restart with the initial set of knots the user */
/*   must set iopt(1)=0. */
/*   if iopt(1) = 1 the program will continue with the knots found at */
/*   the last call of the routine. this will save a lot of computation */
/*   time if spgrid is called repeatedly for different values of s. */
/*   the number of knots of the spline returned and their location will */
/*   depend on the value of s and on the complexity of the shape of the */
/*   function underlying the data. if the computation mode iopt(1) = 1 */
/*   is used, the knots returned may also depend on the s-values at */
/*   previous calls (if these were smaller). therefore, if after a number */
/*   of trials with different s-values and iopt(1)=1,the user can finally */
/*   accept a fit as satisfactory, it may be worthwhile for him to call */
/*   spgrid once more with the chosen value for s but now with iopt(1)=0. */
/*   indeed, spgrid may then return an approximation of the same quality */
/*   of fit but with fewer knots and therefore better if data reduction */
/*   is also an important objective for the user. */
/*   the number of knots may also depend on the upper bounds nuest and */
/*   nvest. indeed, if at a certain stage in spgrid the number of knots */
/*   in one direction (say nu) has reached the value of its upper bound */
/*   (nuest), then from that moment on all subsequent knots are added */
/*   in the other (v) direction. this may indicate that the value of */
/*   nuest is too small. on the other hand, it gives the user the option */
/*   of limiting the number of knots the routine locates in any direction */
/*   for example, by setting nuest=8 (the lowest allowable value for */
/*   nuest), the user can indicate that he wants an approximation which */
/*   is a simple cubic polynomial in the variable u. */

/*  other subroutines required: */
/*    fpspgr,fpchec,fpchep,fpknot,fpopsp,fprati,fpgrsp,fpsysy,fpback, */
/*    fpbacp,fpbspl,fpcyt1,fpcyt2,fpdisc,fpgivs,fprota */

/*  references: */
/*   dierckx p. : fast algorithms for smoothing data over a disc or a */
/*                sphere using tensor product splines, in "algorithms */
/*                for approximation", ed. j.c.mason and m.g.cox, */
/*                clarendon press oxford, 1987, pp. 51-65 */
/*   dierckx p. : fast algorithms for smoothing data over a disc or a */
/*                sphere using tensor product splines, report tw73, dept. */
/*                computer science,k.u.leuven, 1985. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author: */
/*    p.dierckx */
/*    dept. computer science, k.u. leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  creation date : july 1985 */
/*  latest update : march 1989 */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpchec,fpchep,fpspgr */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --iopt;
    --ider;
    --u;
    --r__;
    --v;
    --tu;
    --c__;
    --tv;
    --wrk;
    --iwrk;

    /* Function Body */
    one = 1.;
    half = .5f;
    pi = atan2(0., -one);
    per = pi + pi;
    ve = v[1] + per;
/*  we set up the parameters tol and maxit. */
    maxit = 20;
    tol = .001f;
/*  before starting computations, a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (iopt[1] < -1 || iopt[1] > 1) {
	goto L200;
    }
    if (iopt[2] < 0 || iopt[2] > 1) {
	goto L200;
    }
    if (iopt[3] < 0 || iopt[3] > 1) {
	goto L200;
    }
    if (ider[1] < -1 || ider[1] > 1) {
	goto L200;
    }
    if (ider[2] < 0 || ider[2] > 1) {
	goto L200;
    }
    if (ider[2] == 1 && iopt[2] == 0) {
	goto L200;
    }
    if (ider[3] < -1 || ider[3] > 1) {
	goto L200;
    }
    if (ider[4] < 0 || ider[4] > 1) {
	goto L200;
    }
    if (ider[4] == 1 && iopt[3] == 0) {
	goto L200;
    }
    mumin = 4;
    if (ider[1] >= 0) {
	--mumin;
    }
    if (iopt[2] == 1 && ider[2] == 1) {
	--mumin;
    }
    if (ider[3] >= 0) {
	--mumin;
    }
    if (iopt[3] == 1 && ider[4] == 1) {
	--mumin;
    }
    if (mumin == 0) {
	mumin = 1;
    }
    if (*mu < mumin || *mv < 4) {
	goto L200;
    }
    if (*nuest < 8 || *nvest < 8) {
	goto L200;
    }
    m = *mu * *mv;
    nc = (*nuest - 4) * (*nvest - 4);
/* Computing MAX */
    i__1 = *nuest, i__2 = *mv + *nvest;
    lwest = *nuest * (*mv + *nvest + 3) + 12 + *nvest * 24 + (*mu << 2) + (*
	    mv << 3) + max(i__1,i__2);
    kwest = *mu + 5 + *mv + *nuest + *nvest;
    if (*lwrk < lwest || *kwrk < kwest) {
	goto L200;
    }
    if (u[1] <= 0.f || u[*mu] >= pi) {
	goto L200;
    }
    if (*mu == 1) {
	goto L30;
    }
    i__1 = *mu;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (u[i__ - 1] >= u[i__]) {
	    goto L200;
	}
/* L20: */
    }
L30:
    if (v[1] < -pi || v[1] >= pi) {
	goto L200;
    }
    if (v[*mv] >= v[1] + per) {
	goto L200;
    }
    i__1 = *mv;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (v[i__ - 1] >= v[i__]) {
	    goto L200;
	}
/* L40: */
    }
    if (iopt[1] > 0) {
	goto L140;
    }
/*  if not given, we compute an estimate for r0. */
    rn = (doublereal) (*mv);
    if (ider[1] < 0) {
	goto L45;
    }
    rb = *r0;
    goto L55;
L45:
    rb = 0.f;
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rb += r__[i__];
/* L50: */
    }
    rb /= rn;
/*  if not given, we compute an estimate for r1. */
L55:
    if (ider[3] < 0) {
	goto L60;
    }
    re = *r1;
    goto L70;
L60:
    re = 0.f;
    j = m;
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	re += r__[j];
	--j;
/* L65: */
    }
    re /= rn;
/*  we determine the range of r-values. */
L70:
    rmin = rb;
    rmax = re;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (r__[i__] < rmin) {
	    rmin = r__[i__];
	}
	if (r__[i__] > rmax) {
	    rmax = r__[i__];
	}
/* L80: */
    }
    wrk[5] = rb;
    wrk[6] = 0.f;
    wrk[7] = 0.f;
    wrk[8] = re;
    wrk[9] = 0.f;
    wrk[10] = 0.f;
    wrk[11] = rmax - rmin;
    wrk[12] = wrk[11];
    iwrk[4] = *mu;
    iwrk[5] = *mu;
    if (iopt[1] == 0) {
	goto L140;
    }
    if (*nu < 8 || *nu > *nuest) {
	goto L200;
    }
    if (*nv < 11 || *nv > *nvest) {
	goto L200;
    }
    j = *nu;
    for (i__ = 1; i__ <= 4; ++i__) {
	tu[i__] = 0.f;
	tu[j] = pi;
	--j;
/* L90: */
    }
    l = 13;
    wrk[l] = 0.f;
    if (iopt[2] == 0) {
	goto L100;
    }
    ++l;
    uu = u[1];
    if (uu > tu[5]) {
	uu = tu[5];
    }
    wrk[l] = uu * half;
L100:
    i__1 = *mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++l;
	wrk[l] = u[i__];
/* L110: */
    }
    if (iopt[3] == 0) {
	goto L120;
    }
    ++l;
    uu = u[*mu];
    if (uu < tu[*nu - 4]) {
	uu = tu[*nu - 4];
    }
    wrk[l] = uu + (pi - uu) * half;
L120:
    ++l;
    wrk[l] = pi;
    muu = l - 12;
    fpchec_(&wrk[13], &muu, &tu[1], nu, &c__3, ier);
    if (*ier != 0) {
	goto L200;
    }
    j1 = 4;
    tv[j1] = v[1];
    i1 = *nv - 3;
    tv[i1] = ve;
    j2 = j1;
    i2 = i1;
    for (i__ = 1; i__ <= 3; ++i__) {
	++i1;
	--i2;
	++j1;
	--j2;
	tv[j2] = tv[i2] - per;
	tv[i1] = tv[j1] + per;
/* L130: */
    }
    l = 13;
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wrk[l] = v[i__];
	++l;
/* L135: */
    }
    wrk[l] = ve;
    i__1 = *mv + 1;
    fpchep_(&wrk[13], &i__1, &tv[1], nv, &c__3, ier);
    if (*ier == 0) {
	goto L150;
    }
    goto L200;
L140:
    if (*s < 0.f) {
	goto L200;
    }
    if (*s == 0.f && (*nuest < *mu + 6 + iopt[2] + iopt[3] || *nvest < *mv + 
	    7)) {
	goto L200;
    }
/*  we partition the working space and determine the spline approximation */
L150:
    ldr = 5;
    lfpu = 13;
    lfpv = lfpu + *nuest;
    lww = lfpv + *nvest;
    jwrk = *lwrk - 12 - *nuest - *nvest;
    knru = 6;
    knrv = knru + *mu;
    kndu = knrv + *mv;
    kndv = kndu + *nuest;
    fpspgr_(&iopt[1], &ider[1], &u[1], mu, &v[1], mv, &r__[1], &m, &rb, &re, 
	    s, nuest, nvest, &tol, &maxit, &nc, nu, &tu[1], nv, &tv[1], &c__[
	    1], fp, &wrk[1], &wrk[2], &wrk[3], &wrk[4], &wrk[lfpu], &wrk[lfpv]
	    , &wrk[ldr], &wrk[11], &iwrk[1], &iwrk[2], &iwrk[3], &iwrk[4], &
	    iwrk[5], &iwrk[knru], &iwrk[knrv], &iwrk[kndu], &iwrk[kndv], &wrk[
	    lww], &jwrk, ier);
L200:
    return 0;
} /* spgrid_ */

