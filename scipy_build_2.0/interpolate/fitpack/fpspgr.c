/* fpspgr.f -- translated by f2c (version 20190311).
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

static integer c__1 = 1;

/* Subroutine */ int fpspgr_(integer *iopt, integer *ider, doublereal *u, 
	integer *mu, doublereal *v, integer *mv, doublereal *r__, integer *mr,
	 doublereal *r0, doublereal *r1, doublereal *s, integer *nuest, 
	integer *nvest, doublereal *tol, integer *maxit, integer *nc, integer 
	*nu, doublereal *tu, integer *nv, doublereal *tv, doublereal *c__, 
	doublereal *fp, doublereal *fp0, doublereal *fpold, doublereal *
	reducu, doublereal *reducv, doublereal *fpintu, doublereal *fpintv, 
	doublereal *dr, doublereal *step, integer *lastdi, integer *nplusu, 
	integer *nplusv, integer *lastu0, integer *lastu1, integer *nru, 
	integer *nrv, integer *nrdatu, integer *nrdatv, doublereal *wrk, 
	integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    double atan2(doublereal, doublereal);

    /* Local variables */
    static integer i__, j, l;
    static doublereal p, f1, f2, f3;
    static integer i1, i2, l1, l2, l3, l4;
    static doublereal p1, p2, p3, vb, pi, ve;
    static integer nn;
    static doublereal rn;
    static integer ju, mu0, mu1;
    static doublereal acc;
    static integer idd[4];
    static doublereal one, per;
    static integer nue, nve, mpm;
    static doublereal drr[6];
    static integer ktu, ich1, ich3;
    static doublereal con1, con4, con9;
    static integer npl1, ifbu, ifbv, iter;
    static doublereal fpms, rmin, rmax;
    static integer ifsu, ifsv, nplu, nplv, mumin, numax, nvmax;
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpknot_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *);
    static integer istart;
    extern /* Subroutine */ int fpopsp_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, integer *);
    static integer nrintu, nrintv;

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpknot,fpopsp */
/*  .. */
/*   set constants */
    /* Parameter adjustments */
    --iopt;
    --ider;
    --nru;
    --u;
    --nrv;
    --v;
    --r__;
    --nrdatu;
    --fpintu;
    --tu;
    --nrdatv;
    --fpintv;
    --tv;
    --c__;
    --dr;
    --step;
    --wrk;

    /* Function Body */
    one = 1.;
    con1 = .1f;
    con9 = .9f;
    con4 = .04f;
/*   initialization */
    ifsu = 0;
    ifsv = 0;
    ifbu = 0;
    ifbv = 0;
    p = -one;
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
    pi = atan2(0., -one);
    per = pi + pi;
    vb = v[1];
    ve = vb + per;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 1: determination of the number of knots and their position.     c */
/* ****************************************************************     c */
/*  given a set of knots we compute the least-squares spline sinf(u,v)  c */
/*  and the corresponding sum of squared residuals fp = f(p=inf).       c */
/*  if iopt(1)=-1  sinf(u,v) is the requested approximation.            c */
/*  if iopt(1)>=0  we check whether we can accept the knots:            c */
/*    if fp <= s we will continue with the current set of knots.        c */
/*    if fp >  s we will increase the number of knots and compute the   c */
/*       corresponding least-squares spline until finally fp <= s.      c */
/*    the initial choice of knots depends on the value of s and iopt.   c */
/*    if s=0 we have spline interpolation; in that case the number of   c */
/*     knots in the u-direction equals nu=numax=mu+6+iopt(2)+iopt(3)    c */
/*     and in the v-direction nv=nvmax=mv+7.                            c */
/*    if s>0 and                                                        c */
/*      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c */
/*       spline without interior knots : nu=8 ; nv=8.                   c */
/*      iopt(1)=1 we start with the set of knots found at the last call c */
/*       of the routine, except for the case that s > fp0; then we      c */
/*       compute the least-squares polynomial directly.                 c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    if (iopt[1] < 0) {
	goto L120;
    }
/*  acc denotes the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
/*  numax and nvmax denote the number of knots needed for interpolation. */
    numax = *mu + 6 + iopt[2] + iopt[3];
    nvmax = *mv + 7;
    nue = min(numax,*nuest);
    nve = min(nvmax,*nvest);
    if (*s > 0.f) {
	goto L100;
    }
/*  if s = 0, s(u,v) is an interpolating spline. */
    *nu = numax;
    *nv = nvmax;
/*  test whether the required storage space exceeds the available one. */
    if (*nu > *nuest || *nv > *nvest) {
	goto L420;
    }
/*  find the position of the knots in the v-direction. */
    i__1 = *mv;
    for (l = 1; l <= i__1; ++l) {
	tv[l + 3] = v[l];
/* L10: */
    }
    tv[*mv + 4] = ve;
    l1 = *mv - 2;
    l2 = *mv + 5;
    for (i__ = 1; i__ <= 3; ++i__) {
	tv[i__] = v[l1] - per;
	tv[l2] = v[i__ + 1] + per;
	++l1;
	++l2;
/* L20: */
    }
/*  if not all the derivative values g(i,j) are given, we will first */
/*  estimate these values by computing a least-squares spline */
    idd[0] = ider[1];
    if (idd[0] == 0) {
	idd[0] = 1;
    }
    if (idd[0] > 0) {
	dr[1] = *r0;
    }
    idd[1] = ider[2];
    idd[2] = ider[3];
    if (idd[2] == 0) {
	idd[2] = 1;
    }
    if (idd[2] > 0) {
	dr[4] = *r1;
    }
    idd[3] = ider[4];
    if (ider[1] < 0 || ider[3] < 0) {
	goto L30;
    }
    if (iopt[2] != 0 && ider[2] == 0) {
	goto L30;
    }
    if (iopt[3] == 0 || ider[4] != 0) {
	goto L70;
    }
/* we set up the knots in the u-direction for computing the least-squares */
/* spline. */
L30:
    i1 = 3;
    i2 = *mu - 2;
    *nu = 4;
    i__1 = *mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i1 > i2) {
	    goto L50;
	}
	++(*nu);
	tu[*nu] = u[i1];
	i1 += 2;
/* L40: */
    }
L50:
    for (i__ = 1; i__ <= 4; ++i__) {
	tu[i__] = 0.f;
	++(*nu);
	tu[*nu] = pi;
/* L60: */
    }
/* we compute the least-squares spline for estimating the derivatives. */
    fpopsp_(&ifsu, &ifsv, &ifbu, &ifbv, &u[1], mu, &v[1], mv, &r__[1], mr, r0,
	     r1, &dr[1], &iopt[1], idd, &tu[1], nu, &tv[1], nv, nuest, nvest, 
	    &p, &step[1], &c__[1], nc, fp, &fpintu[1], &fpintv[1], &nru[1], &
	    nrv[1], &wrk[1], lwrk);
    ifsu = 0;
/* if all the derivatives at the origin are known, we compute the */
/* interpolating spline. */
/* we set up the knots in the u-direction, needed for interpolation. */
L70:
    nn = numax - 8;
    if (nn == 0) {
	goto L95;
    }
    ju = 2 - iopt[2];
    i__1 = nn;
    for (l = 1; l <= i__1; ++l) {
	tu[l + 4] = u[ju];
	++ju;
/* L80: */
    }
    *nu = numax;
    l = *nu;
    for (i__ = 1; i__ <= 4; ++i__) {
	tu[i__] = 0.f;
	tu[l] = pi;
	--l;
/* L90: */
    }
/* we compute the interpolating spline. */
L95:
    fpopsp_(&ifsu, &ifsv, &ifbu, &ifbv, &u[1], mu, &v[1], mv, &r__[1], mr, r0,
	     r1, &dr[1], &iopt[1], idd, &tu[1], nu, &tv[1], nv, nuest, nvest, 
	    &p, &step[1], &c__[1], nc, fp, &fpintu[1], &fpintv[1], &nru[1], &
	    nrv[1], &wrk[1], lwrk);
    goto L430;
/*  if s>0 our initial choice of knots depends on the value of iopt(1). */
L100:
    *ier = 0;
    if (iopt[1] == 0) {
	goto L115;
    }
    step[1] = -step[1];
    step[2] = -step[2];
    if (*fp0 <= *s) {
	goto L115;
    }
/*  if iopt(1)=1 and fp0 > s we start computing the least-squares spline */
/*  according to the set of knots found at the last call of the routine. */
/*  we determine the number of grid coordinates u(i) inside each knot */
/*  interval (tu(l),tu(l+1)). */
    l = 5;
    j = 1;
    nrdatu[1] = 0;
    mu0 = 2 - iopt[2];
    mu1 = *mu - 1 + iopt[3];
    i__1 = mu1;
    for (i__ = mu0; i__ <= i__1; ++i__) {
	++nrdatu[j];
	if (u[i__] < tu[l]) {
	    goto L105;
	}
	--nrdatu[j];
	++l;
	++j;
	nrdatu[j] = 0;
L105:
	;
    }
/*  we determine the number of grid coordinates v(i) inside each knot */
/*  interval (tv(l),tv(l+1)). */
    l = 5;
    j = 1;
    nrdatv[1] = 0;
    i__1 = *mv;
    for (i__ = 2; i__ <= i__1; ++i__) {
	++nrdatv[j];
	if (v[i__] < tv[l]) {
	    goto L110;
	}
	--nrdatv[j];
	++l;
	++j;
	nrdatv[j] = 0;
L110:
	;
    }
    idd[0] = ider[1];
    idd[1] = ider[2];
    idd[2] = ider[3];
    idd[3] = ider[4];
    goto L120;
/*  if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least- */
/*  squares polynomial (which is a spline without interior knots). */
L115:
    *ier = -2;
    idd[0] = ider[1];
    idd[1] = 1;
    idd[2] = ider[3];
    idd[3] = 1;
    *nu = 8;
    *nv = 8;
    nrdatu[1] = *mu - 2 + iopt[2] + iopt[3];
    nrdatv[1] = *mv - 1;
    *lastdi = 0;
    *nplusu = 0;
    *nplusv = 0;
    *fp0 = 0.f;
    *fpold = 0.f;
    *reducu = 0.f;
    *reducv = 0.f;
/*  main loop for the different sets of knots.mpm=mu+mv is a save upper */
/*  bound for the number of trials. */
L120:
    mpm = *mu + *mv;
    i__1 = mpm;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find nrintu (nrintv) which is the number of knot intervals in the */
/*  u-direction (v-direction). */
	nrintu = *nu - 7;
	nrintv = *nv - 7;
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(u,v). */
	i__ = *nu;
	for (j = 1; j <= 4; ++j) {
	    tu[j] = 0.f;
	    tu[i__] = pi;
	    --i__;
/* L125: */
	}
	l1 = 4;
	l2 = l1;
	l3 = *nv - 3;
	l4 = l3;
	tv[l2] = vb;
	tv[l3] = ve;
	for (j = 1; j <= 3; ++j) {
	    ++l1;
	    --l2;
	    ++l3;
	    --l4;
	    tv[l2] = tv[l4] - per;
	    tv[l3] = tv[l1] + per;
/* L130: */
	}
/*  find an estimate of the range of possible values for the optimal */
/*  derivatives at the origin. */
	ktu = nrdatu[1] + 2 - iopt[2];
	if (ktu < mumin) {
	    ktu = mumin;
	}
	if (ktu == *lastu0) {
	    goto L140;
	}
	rmin = *r0;
	rmax = *r0;
	l = *mv * ktu;
	i__2 = l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (r__[i__] < rmin) {
		rmin = r__[i__];
	    }
	    if (r__[i__] > rmax) {
		rmax = r__[i__];
	    }
/* L135: */
	}
	step[1] = rmax - rmin;
	*lastu0 = ktu;
L140:
	ktu = nrdatu[nrintu] + 2 - iopt[3];
	if (ktu < mumin) {
	    ktu = mumin;
	}
	if (ktu == *lastu1) {
	    goto L150;
	}
	rmin = *r1;
	rmax = *r1;
	l = *mv * ktu;
	j = *mr;
	i__2 = l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (r__[j] < rmin) {
		rmin = r__[j];
	    }
	    if (r__[j] > rmax) {
		rmax = r__[j];
	    }
	    --j;
/* L145: */
	}
	step[2] = rmax - rmin;
	*lastu1 = ktu;
/*  find the least-squares spline sinf(u,v). */
L150:
	fpopsp_(&ifsu, &ifsv, &ifbu, &ifbv, &u[1], mu, &v[1], mv, &r__[1], mr,
		 r0, r1, &dr[1], &iopt[1], idd, &tu[1], nu, &tv[1], nv, nuest,
		 nvest, &p, &step[1], &c__[1], nc, fp, &fpintu[1], &fpintv[1],
		 &nru[1], &nrv[1], &wrk[1], lwrk);
	if (step[1] < 0.f) {
	    step[1] = -step[1];
	}
	if (step[2] < 0.f) {
	    step[2] = -step[2];
	}
	if (*ier == -2) {
	    *fp0 = *fp;
	}
/*  test whether the least-squares spline is an acceptable solution. */
	if (iopt[1] < 0) {
	    goto L440;
	}
	fpms = *fp - *s;
	if (abs(fpms) < acc) {
	    goto L440;
	}
/*  if f(p=inf) < s, we accept the choice of knots. */
	if (fpms < 0.f) {
	    goto L300;
	}
/*  if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline */
	if (*nu == numax && *nv == nvmax) {
	    goto L430;
	}
/*  increase the number of knots. */
/*  if nu=nue and nv=nve we cannot further increase the number of knots */
/*  because of the storage capacity limitation. */
	if (*nu == nue && *nv == nve) {
	    goto L420;
	}
	if (ider[1] == 0) {
/* Computing 2nd power */
	    d__1 = *r0 - dr[1];
	    fpintu[1] += d__1 * d__1;
	}
	if (ider[3] == 0) {
/* Computing 2nd power */
	    d__1 = *r1 - dr[4];
	    fpintu[nrintu] += d__1 * d__1;
	}
	*ier = 0;
/*  adjust the parameter reducu or reducv according to the direction */
/*  in which the last added knots were located. */
	if (*lastdi < 0) {
	    goto L160;
	}
	if (*lastdi == 0) {
	    goto L155;
	}
	goto L170;
L155:
	nplv = 3;
	idd[1] = ider[2];
	idd[3] = ider[4];
	*fpold = *fp;
	goto L230;
L160:
	*reducu = *fpold - *fp;
	goto L175;
L170:
	*reducv = *fpold - *fp;
/*  store the sum of squared residuals for the current set of knots. */
L175:
	*fpold = *fp;
/*  find nplu, the number of knots we should add in the u-direction. */
	nplu = 1;
	if (*nu == 8) {
	    goto L180;
	}
	npl1 = *nplusu << 1;
	rn = (doublereal) (*nplusu);
	if (*reducu > acc) {
	    npl1 = (integer) (rn * fpms / *reducu);
	}
/* Computing MIN */
/* Computing MAX */
	i__4 = npl1, i__5 = *nplusu / 2, i__4 = max(i__4,i__5);
	i__2 = *nplusu << 1, i__3 = max(i__4,1);
	nplu = min(i__2,i__3);
/*  find nplv, the number of knots we should add in the v-direction. */
L180:
	nplv = 3;
	if (*nv == 8) {
	    goto L190;
	}
	npl1 = *nplusv << 1;
	rn = (doublereal) (*nplusv);
	if (*reducv > acc) {
	    npl1 = (integer) (rn * fpms / *reducv);
	}
/* Computing MIN */
/* Computing MAX */
	i__4 = npl1, i__5 = *nplusv / 2, i__4 = max(i__4,i__5);
	i__2 = *nplusv << 1, i__3 = max(i__4,1);
	nplv = min(i__2,i__3);
/*  test whether we are going to add knots in the u- or v-direction. */
L190:
	if (nplu < nplv) {
	    goto L210;
	}
	if (nplu == nplv) {
	    goto L200;
	}
	goto L230;
L200:
	if (*lastdi < 0) {
	    goto L230;
	}
L210:
	if (*nu == nue) {
	    goto L230;
	}
/*  addition in the u-direction. */
	*lastdi = -1;
	*nplusu = nplu;
	ifsu = 0;
	istart = 0;
	if (iopt[2] == 0) {
	    istart = 1;
	}
	i__2 = *nplusu;
	for (l = 1; l <= i__2; ++l) {
/*  add a new knot in the u-direction */
	    fpknot_(&u[1], mu, &tu[1], nu, &fpintu[1], &nrdatu[1], &nrintu, 
		    nuest, &istart);
/*  test whether we cannot further increase the number of knots in the */
/*  u-direction. */
	    if (*nu == nue) {
		goto L270;
	    }
/* L220: */
	}
	goto L270;
L230:
	if (*nv == nve) {
	    goto L210;
	}
/*  addition in the v-direction. */
	*lastdi = 1;
	*nplusv = nplv;
	ifsv = 0;
	i__2 = *nplusv;
	for (l = 1; l <= i__2; ++l) {
/*  add a new knot in the v-direction. */
	    fpknot_(&v[1], mv, &tv[1], nv, &fpintv[1], &nrdatv[1], &nrintv, 
		    nvest, &c__1);
/*  test whether we cannot further increase the number of knots in the */
/*  v-direction. */
	    if (*nv == nve) {
		goto L270;
	    }
/* L240: */
	}
/*  restart the computations with the new set of knots. */
L270:
	;
    }
/*  test whether the least-squares polynomial is a solution of our */
/*  approximation problem. */
L300:
    if (*ier == -2) {
	goto L440;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing spline sp(u,v)                c */
/* *****************************************************                c */
/*  we have determined the number of knots and their position. we now   c */
/*  compute the b-spline coefficients of the smoothing spline sp(u,v).  c */
/*  this smoothing spline depends on the parameter p in such a way that c */
/*    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c */
/*  is a continuous, strictly decreasing function of p. moreover the    c */
/*  least-squares polynomial corresponds to p=0 and the least-squares   c */
/*  spline to p=infinity. then iteratively we have to determine the     c */
/*  positive value of p such that f(p)=s. the process which is proposed c */
/*  here makes use of rational interpolation. f(p) is approximated by a c */
/*  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c */
/*  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c */
/*  are used to calculate the new value of p such that r(p)=s.          c */
/*  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  initial value for p. */
    p1 = 0.f;
    f1 = *fp0 - *s;
    p3 = -one;
    f3 = fpms;
    p = one;
    for (i__ = 1; i__ <= 6; ++i__) {
	drr[i__ - 1] = dr[i__];
/* L305: */
    }
    ich1 = 0;
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the smoothing spline sp(u,v) and the corresponding sum f(p). */
	fpopsp_(&ifsu, &ifsv, &ifbu, &ifbv, &u[1], mu, &v[1], mv, &r__[1], mr,
		 r0, r1, drr, &iopt[1], idd, &tu[1], nu, &tv[1], nv, nuest, 
		nvest, &p, &step[1], &c__[1], nc, fp, &fpintu[1], &fpintv[1], 
		&nru[1], &nrv[1], &wrk[1], lwrk);
/*  test whether the approximation sp(u,v) is an acceptable solution. */
	fpms = *fp - *s;
	if (abs(fpms) < acc) {
	    goto L440;
	}
/*  test whether the maximum allowable number of iterations has been */
/*  reached. */
	if (iter == *maxit) {
	    goto L400;
	}
/*  carry out one more step of the iteration process. */
	p2 = p;
	f2 = fpms;
	if (ich3 != 0) {
	    goto L320;
	}
	if (f2 - f3 > acc) {
	    goto L310;
	}
/*  our initial choice of p is too large. */
	p3 = p2;
	f3 = f2;
	p *= con4;
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
	goto L350;
L310:
	if (f2 < 0.f) {
	    ich3 = 1;
	}
L320:
	if (ich1 != 0) {
	    goto L340;
	}
	if (f1 - f2 > acc) {
	    goto L330;
	}
/*  our initial choice of p is too small */
	p1 = p2;
	f1 = f2;
	p /= con4;
	if (p3 < 0.f) {
	    goto L350;
	}
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
	goto L350;
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
L330:
	if (f2 > 0.f) {
	    ich1 = 1;
	}
L340:
	if (f2 >= f1 || f2 <= f3) {
	    goto L410;
	}
/*  find the new value of p. */
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
L350:
	;
    }
/*  error codes and messages. */
L400:
    *ier = 3;
    goto L440;
L410:
    *ier = 2;
    goto L440;
L420:
    *ier = 1;
    goto L440;
L430:
    *ier = -1;
    *fp = 0.f;
L440:
    return 0;
} /* fpspgr_ */

