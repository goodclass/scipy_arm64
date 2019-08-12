/* parsur.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int parsur_(integer *iopt, integer *ipar, integer *idim, 
	integer *mu, doublereal *u, integer *mv, doublereal *v, doublereal *f,
	 doublereal *s, integer *nuest, integer *nvest, integer *nu, 
	doublereal *tu, integer *nv, doublereal *tv, doublereal *c__, 
	doublereal *fp, doublereal *wrk, integer *lwrk, integer *iwrk, 
	integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, l1, l2, l3, l4, nc, mf;
    static doublereal ub, vb, ue, ve, tol;
    static integer lww, kndu, kndv, lfpu, lfpv;
    static doublereal peru, perv;
    static integer jwrk, knru, knrv, maxit, mumin, mvmin, kwest, lwest;
    extern /* Subroutine */ int fpchec_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *), fpchep_(doublereal *, integer *
	    , doublereal *, integer *, integer *, integer *), fppasu_(integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *);

/*  given the set of ordered points f(i,j) in the idim-dimensional space, */
/*  corresponding to grid values (u(i),v(j)) ,i=1,...,mu ; j=1,...,mv, */
/*  parsur determines a smooth approximating spline surface s(u,v) , i.e. */
/*    f1 = s1(u,v) */
/*      ...                u(1) <= u <= u(mu) ; v(1) <= v <= v(mv) */
/*    fidim = sidim(u,v) */
/*  with sl(u,v), l=1,2,...,idim bicubic spline functions with common */
/*  knots tu(i),i=1,...,nu in the u-variable and tv(j),j=1,...,nv in the */
/*  v-variable. */
/*  in addition, these splines will be periodic in the variable u if */
/*  ipar(1) = 1 and periodic in the variable v if ipar(2) = 1. */
/*  if iopt=-1, parsur determines the least-squares bicubic spline */
/*  surface according to a given set of knots. */
/*  if iopt>=0, the number of knots of s(u,v) and their position */
/*  is chosen automatically by the routine. the smoothness of s(u,v) is */
/*  achieved by minimalizing the discontinuity jumps of the derivatives */
/*  of the splines at the knots. the amount of smoothness of s(u,v) is */
/*  determined by the condition that */
/*  fp=sumi=1,mu(sumj=1,mv(dist(f(i,j)-s(u(i),v(j)))**2))<=s, */
/*  with s a given non-negative constant. */
/*  the fit s(u,v) is given in its b-spline representation and can be */
/*  evaluated by means of routine surev. */

/* calling sequence: */
/*     call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,nu,tu, */
/*    *  nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier) */

/* parameters: */
/*  iopt  : integer flag. unchanged on exit. */
/*          on entry iopt must specify whether a least-squares surface */
/*          (iopt=-1) or a smoothing surface (iopt=0 or 1)must be */
/*          determined. */
/*          if iopt=0 the routine will start with the initial set of */
/*          knots needed for determining the least-squares polynomial */
/*          surface. */
/*          if iopt=1 the routine will continue with the set of knots */
/*          found at the last call of the routine. */
/*          attention: a call with iopt=1 must always be immediately */
/*          preceded by another call with iopt = 1 or iopt = 0. */
/*  ipar  : integer array of dimension 2. unchanged on exit. */
/*          on entry ipar(1) must specify whether (ipar(1)=1) or not */
/*          (ipar(1)=0) the splines must be periodic in the variable u. */
/*          on entry ipar(2) must specify whether (ipar(2)=1) or not */
/*          (ipar(2)=0) the splines must be periodic in the variable v. */
/*  idim  : integer. on entry idim must specify the dimension of the */
/*          surface. 1 <= idim <= 3. unchanged on exit. */
/*  mu    : integer. on entry mu must specify the number of grid points */
/*          along the u-axis. unchanged on exit. */
/*          mu >= mumin where mumin=4-2*ipar(1) */
/*  u     : real array of dimension at least (mu). before entry, u(i) */
/*          must be set to the u-co-ordinate of the i-th grid point */
/*          along the u-axis, for i=1,2,...,mu. these values must be */
/*          supplied in strictly ascending order. unchanged on exit. */
/*  mv    : integer. on entry mv must specify the number of grid points */
/*          along the v-axis. unchanged on exit. */
/*          mv >= mvmin where mvmin=4-2*ipar(2) */
/*  v     : real array of dimension at least (mv). before entry, v(j) */
/*          must be set to the v-co-ordinate of the j-th grid point */
/*          along the v-axis, for j=1,2,...,mv. these values must be */
/*          supplied in strictly ascending order. unchanged on exit. */
/*  f     : real array of dimension at least (mu*mv*idim). */
/*          before entry, f(mu*mv*(l-1)+mv*(i-1)+j) must be set to the */
/*          l-th co-ordinate of the data point corresponding to the */
/*          the grid point (u(i),v(j)) for l=1,...,idim ,i=1,...,mu */
/*          and j=1,...,mv. unchanged on exit. */
/*          if ipar(1)=1 it is expected that f(mu*mv*(l-1)+mv*(mu-1)+j) */
/*          = f(mu*mv*(l-1)+j), l=1,...,idim ; j=1,...,mv */
/*          if ipar(2)=1 it is expected that f(mu*mv*(l-1)+mv*(i-1)+mv) */
/*          = f(mu*mv*(l-1)+mv*(i-1)+1), l=1,...,idim ; i=1,...,mu */
/*  s     : real. on entry (if iopt>=0) s must specify the smoothing */
/*          factor. s >=0. unchanged on exit. */
/*          for advice on the choice of s see further comments */
/*  nuest : integer. unchanged on exit. */
/*  nvest : integer. unchanged on exit. */
/*          on entry, nuest and nvest must specify an upper bound for the */
/*          number of knots required in the u- and v-directions respect. */
/*          these numbers will also determine the storage space needed by */
/*          the routine. nuest >= 8, nvest >= 8. */
/*          in most practical situation nuest = mu/2, nvest=mv/2, will */
/*          be sufficient. always large enough are nuest=mu+4+2*ipar(1), */
/*          nvest = mv+4+2*ipar(2), the number of knots needed for */
/*          interpolation (s=0). see also further comments. */
/*  nu    : integer. */
/*          unless ier=10 (in case iopt>=0), nu will contain the total */
/*          number of knots with respect to the u-variable, of the spline */
/*          surface returned. if the computation mode iopt=1 is used, */
/*          the value of nu should be left unchanged between subsequent */
/*          calls. in case iopt=-1, the value of nu should be specified */
/*          on entry. */
/*  tu    : real array of dimension at least (nuest). */
/*          on successful exit, this array will contain the knots of the */
/*          splines with respect to the u-variable, i.e. the position of */
/*          the interior knots tu(5),...,tu(nu-4) as well as the position */
/*          of the additional knots tu(1),...,tu(4) and tu(nu-3),..., */
/*          tu(nu) needed for the b-spline representation. */
/*          if the computation mode iopt=1 is used,the values of tu(1) */
/*          ...,tu(nu) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt=-1 is used, the values tu(5), */
/*          ...tu(nu-4) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  nv    : integer. */
/*          unless ier=10 (in case iopt>=0), nv will contain the total */
/*          number of knots with respect to the v-variable, of the spline */
/*          surface returned. if the computation mode iopt=1 is used, */
/*          the value of nv should be left unchanged between subsequent */
/*          calls. in case iopt=-1, the value of nv should be specified */
/*          on entry. */
/*  tv    : real array of dimension at least (nvest). */
/*          on successful exit, this array will contain the knots of the */
/*          splines with respect to the v-variable, i.e. the position of */
/*          the interior knots tv(5),...,tv(nv-4) as well as the position */
/*          of the additional knots tv(1),...,tv(4) and tv(nv-3),..., */
/*          tv(nv) needed for the b-spline representation. */
/*          if the computation mode iopt=1 is used,the values of tv(1) */
/*          ...,tv(nv) should be left unchanged between subsequent calls. */
/*          if the computation mode iopt=-1 is used, the values tv(5), */
/*          ...tv(nv-4) must be supplied by the user, before entry. */
/*          see also the restrictions (ier=10). */
/*  c     : real array of dimension at least (nuest-4)*(nvest-4)*idim. */
/*          on successful exit, c contains the coefficients of the spline */
/*          approximation s(u,v) */
/*  fp    : real. unless ier=10, fp contains the sum of squared */
/*          residuals of the spline surface returned. */
/*  wrk   : real array of dimension (lwrk). used as workspace. */
/*          if the computation mode iopt=1 is used the values of */
/*          wrk(1),...,wrk(4) should be left unchanged between subsequent */
/*          calls. */
/*  lwrk  : integer. on entry lwrk must specify the actual dimension of */
/*          the array wrk as declared in the calling (sub)program. */
/*          lwrk must not be too small. */
/*           lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+ */
/*           4*(mu+mv)+q*idim where q is the larger of mv and nuest. */
/*  iwrk  : integer array of dimension (kwrk). used as workspace. */
/*          if the computation mode iopt=1 is used the values of */
/*          iwrk(1),.,iwrk(3) should be left unchanged between subsequent */
/*          calls. */
/*  kwrk  : integer. on entry kwrk must specify the actual dimension of */
/*          the array iwrk as declared in the calling (sub)program. */
/*          kwrk >= 3+mu+mv+nuest+nvest. */
/*  ier   : integer. unless the routine detects an error, ier contains a */
/*          non-positive value on exit, i.e. */
/*   ier=0  : normal return. the surface returned has a residual sum of */
/*            squares fp such that abs(fp-s)/s <= tol with tol a relat- */
/*            ive tolerance set to 0.001 by the program. */
/*   ier=-1 : normal return. the spline surface returned is an */
/*            interpolating surface (fp=0). */
/*   ier=-2 : normal return. the surface returned is the least-squares */
/*            polynomial surface. in this extreme case fp gives the */
/*            upper bound for the smoothing factor s. */
/*   ier=1  : error. the required storage space exceeds the available */
/*            storage space, as specified by the parameters nuest and */
/*            nvest. */
/*            probably causes : nuest or nvest too small. if these param- */
/*            eters are already large, it may also indicate that s is */
/*            too small */
/*            the approximation returned is the least-squares surface */
/*            according to the current set of knots. the parameter fp */
/*            gives the corresponding sum of squared residuals (fp>s). */
/*   ier=2  : error. a theoretically impossible result was found during */
/*            the iteration process for finding a smoothing surface with */
/*            fp = s. probably causes : s too small. */
/*            there is an approximation returned but the corresponding */
/*            sum of squared residuals does not satisfy the condition */
/*            abs(fp-s)/s < tol. */
/*   ier=3  : error. the maximal number of iterations maxit (set to 20 */
/*            by the program) allowed for finding a smoothing surface */
/*            with fp=s has been reached. probably causes : s too small */
/*            there is an approximation returned but the corresponding */
/*            sum of squared residuals does not satisfy the condition */
/*            abs(fp-s)/s < tol. */
/*   ier=10 : error. on entry, the input data are controlled on validity */
/*            the following restrictions must be satisfied. */
/*            -1<=iopt<=1, 0<=ipar(1)<=1, 0<=ipar(2)<=1, 1 <=idim<=3 */
/*            mu >= 4-2*ipar(1),mv >= 4-2*ipar(2), nuest >=8, nvest >= 8, */
/*            kwrk>=3+mu+mv+nuest+nvest, */
/*            lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2)) */
/*             +4*(mu+mv)+max(nuest,mv)*idim */
/*            u(i-1)<u(i),i=2,..,mu, v(i-1)<v(i),i=2,...,mv */
/*            if iopt=-1: 8<=nu<=min(nuest,mu+4+2*ipar(1)) */
/*                        u(1)<tu(5)<tu(6)<...<tu(nu-4)<u(mu) */
/*                        8<=nv<=min(nvest,mv+4+2*ipar(2)) */
/*                        v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(mv) */
/*                    the schoenberg-whitney conditions, i.e. there must */
/*                    be subset of grid co-ordinates uu(p) and vv(q) such */
/*                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4 */
/*                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4 */
/*                     (see fpchec or fpchep) */
/*            if iopt>=0: s>=0 */
/*                       if s=0: nuest>=mu+4+2*ipar(1) */
/*                               nvest>=mv+4+2*ipar(2) */
/*            if one of these conditions is found to be violated,control */
/*            is immediately repassed to the calling program. in that */
/*            case there is no approximation returned. */

/* further comments: */
/*   by means of the parameter s, the user can control the tradeoff */
/*   between closeness of fit and smoothness of fit of the approximation. */
/*   if s is too large, the surface will be too smooth and signal will be */
/*   lost ; if s is too small the surface will pick up too much noise. in */
/*   the extreme cases the program will return an interpolating surface */
/*   if s=0 and the constrained least-squares polynomial surface if s is */
/*   very large. between these extremes, a properly chosen s will result */
/*   in a good compromise between closeness of fit and smoothness of fit. */
/*   to decide whether an approximation, corresponding to a certain s is */
/*   satisfactory the user is highly recommended to inspect the fits */
/*   graphically. */
/*   recommended values for s depend on the accuracy of the data values. */
/*   if the user has an idea of the statistical errors on the data, he */
/*   can also find a proper estimate for s. for, by assuming that, if he */
/*   specifies the right s, parsur will return a surface s(u,v) which */
/*   exactly reproduces the surface underlying the data he can evaluate */
/*   the sum(dist(f(i,j)-s(u(i),v(j)))**2) to find a good estimate for s. */
/*   for example, if he knows that the statistical errors on his f(i,j)- */
/*   values is not greater than 0.1, he may expect that a good s should */
/*   have a value not larger than mu*mv*(0.1)**2. */
/*   if nothing is known about the statistical error in f(i,j), s must */
/*   be determined by trial and error, taking account of the comments */
/*   above. the best is then to start with a very large value of s (to */
/*   determine the le-sq polynomial surface and the corresponding upper */
/*   bound fp0 for s) and then to progressively decrease the value of s */
/*   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,... */
/*   and more carefully as the approximation shows more detail) to */
/*   obtain closer fits. */
/*   to economize the search for a good s-value the program provides with */
/*   different modes of computation. at the first call of the routine, or */
/*   whenever he wants to restart with the initial set of knots the user */
/*   must set iopt=0. */
/*   if iopt = 1 the program will continue with the knots found at */
/*   the last call of the routine. this will save a lot of computation */
/*   time if parsur is called repeatedly for different values of s. */
/*   the number of knots of the surface returned and their location will */
/*   depend on the value of s and on the complexity of the shape of the */
/*   surface underlying the data. if the computation mode iopt = 1 */
/*   is used, the knots returned may also depend on the s-values at */
/*   previous calls (if these were smaller). therefore, if after a number */
/*   of trials with different s-values and iopt=1,the user can finally */
/*   accept a fit as satisfactory, it may be worthwhile for him to call */
/*   parsur once more with the chosen value for s but now with iopt=0. */
/*   indeed, parsur may then return an approximation of the same quality */
/*   of fit but with fewer knots and therefore better if data reduction */
/*   is also an important objective for the user. */
/*   the number of knots may also depend on the upper bounds nuest and */
/*   nvest. indeed, if at a certain stage in parsur the number of knots */
/*   in one direction (say nu) has reached the value of its upper bound */
/*   (nuest), then from that moment on all subsequent knots are added */
/*   in the other (v) direction. this may indicate that the value of */
/*   nuest is too small. on the other hand, it gives the user the option */
/*   of limiting the number of knots the routine locates in any direction */
/*   for example, by setting nuest=8 (the lowest allowable value for */
/*   nuest), the user can indicate that he wants an approximation with */
/*   splines which are simple cubic polynomials in the variable u. */

/*  other subroutines required: */
/*    fppasu,fpchec,fpchep,fpknot,fprati,fpgrpa,fptrnp,fpback, */
/*    fpbacp,fpbspl,fptrpe,fpdisc,fpgivs,fprota */

/*  author: */
/*    p.dierckx */
/*    dept. computer science, k.u. leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1989 */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fppasu,fpchec,fpchep */
/*  .. */
/*  we set up the parameters tol and maxit. */
    /* Parameter adjustments */
    --ipar;
    --u;
    --f;
    --v;
    --tu;
    --c__;
    --tv;
    --wrk;
    --iwrk;

    /* Function Body */
    maxit = 20;
    tol = .001f;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (*iopt < -1 || *iopt > 1) {
	goto L200;
    }
    if (ipar[1] < 0 || ipar[1] > 1) {
	goto L200;
    }
    if (ipar[2] < 0 || ipar[2] > 1) {
	goto L200;
    }
    if (*idim <= 0 || *idim > 3) {
	goto L200;
    }
    mumin = 4 - (ipar[1] << 1);
    if (*mu < mumin || *nuest < 8) {
	goto L200;
    }
    mvmin = 4 - (ipar[2] << 1);
    if (*mv < mvmin || *nvest < 8) {
	goto L200;
    }
    mf = *mu * *mv;
    nc = (*nuest - 4) * (*nvest - 4);
    lwest = *nuest * (*mv * *idim + 11 + (ipar[1] << 2)) + 4 + *nvest * ((
	    ipar[2] << 2) + 11) + (*mu + *mv << 2) + max(*nuest,*mv) * *idim;
    kwest = *mu + 3 + *mv + *nuest + *nvest;
    if (*lwrk < lwest || *kwrk < kwest) {
	goto L200;
    }
    i__1 = *mu;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (u[i__ - 1] >= u[i__]) {
	    goto L200;
	}
/* L10: */
    }
    i__1 = *mv;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (v[i__ - 1] >= v[i__]) {
	    goto L200;
	}
/* L20: */
    }
    if (*iopt >= 0) {
	goto L100;
    }
    if (*nu < 8 || *nu > *nuest) {
	goto L200;
    }
    ub = u[1];
    ue = u[*mu];
    if (ipar[1] != 0) {
	goto L40;
    }
    j = *nu;
    for (i__ = 1; i__ <= 4; ++i__) {
	tu[i__] = ub;
	tu[j] = ue;
	--j;
/* L30: */
    }
    fpchec_(&u[1], mu, &tu[1], nu, &c__3, ier);
    if (*ier != 0) {
	goto L200;
    }
    goto L60;
L40:
    l1 = 4;
    l2 = l1;
    l3 = *nu - 3;
    l4 = l3;
    peru = ue - ub;
    tu[l2] = ub;
    tu[l3] = ue;
    for (j = 1; j <= 3; ++j) {
	++l1;
	--l2;
	++l3;
	--l4;
	tu[l2] = tu[l4] - peru;
	tu[l3] = tu[l1] + peru;
/* L50: */
    }
    fpchep_(&u[1], mu, &tu[1], nu, &c__3, ier);
    if (*ier != 0) {
	goto L200;
    }
L60:
    if (*nv < 8 || *nv > *nvest) {
	goto L200;
    }
    vb = v[1];
    ve = v[*mv];
    if (ipar[2] != 0) {
	goto L80;
    }
    j = *nv;
    for (i__ = 1; i__ <= 4; ++i__) {
	tv[i__] = vb;
	tv[j] = ve;
	--j;
/* L70: */
    }
    fpchec_(&v[1], mv, &tv[1], nv, &c__3, ier);
    if (*ier != 0) {
	goto L200;
    }
    goto L150;
L80:
    l1 = 4;
    l2 = l1;
    l3 = *nv - 3;
    l4 = l3;
    perv = ve - vb;
    tv[l2] = vb;
    tv[l3] = ve;
    for (j = 1; j <= 3; ++j) {
	++l1;
	--l2;
	++l3;
	--l4;
	tv[l2] = tv[l4] - perv;
	tv[l3] = tv[l1] + perv;
/* L90: */
    }
    fpchep_(&v[1], mv, &tv[1], nv, &c__3, ier);
    if (*ier == 0) {
	goto L150;
    }
    goto L200;
L100:
    if (*s < 0.f) {
	goto L200;
    }
    if (*s == 0.f && (*nuest < *mu + 4 + (ipar[1] << 1) || *nvest < *mv + 4 + 
	    (ipar[2] << 1))) {
	goto L200;
    }
    *ier = 0;
/*  we partition the working space and determine the spline approximation */
L150:
    lfpu = 5;
    lfpv = lfpu + *nuest;
    lww = lfpv + *nvest;
    jwrk = *lwrk - 4 - *nuest - *nvest;
    knru = 4;
    knrv = knru + *mu;
    kndu = knrv + *mv;
    kndv = kndu + *nuest;
    fppasu_(iopt, &ipar[1], idim, &u[1], mu, &v[1], mv, &f[1], &mf, s, nuest, 
	    nvest, &tol, &maxit, &nc, nu, &tu[1], nv, &tv[1], &c__[1], fp, &
	    wrk[1], &wrk[2], &wrk[3], &wrk[4], &wrk[lfpu], &wrk[lfpv], &iwrk[
	    1], &iwrk[2], &iwrk[3], &iwrk[knru], &iwrk[knrv], &iwrk[kndu], &
	    iwrk[kndv], &wrk[lww], &jwrk, ier);
L200:
    return 0;
} /* parsur_ */

