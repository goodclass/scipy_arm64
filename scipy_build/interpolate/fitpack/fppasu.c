/* fppasu.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fppasu_(integer *iopt, integer *ipar, integer *idim, 
	doublereal *u, integer *mu, doublereal *v, integer *mv, doublereal *
	z__, integer *mz, doublereal *s, integer *nuest, integer *nvest, 
	doublereal *tol, integer *maxit, integer *nc, integer *nu, doublereal 
	*tu, integer *nv, doublereal *tv, doublereal *c__, doublereal *fp, 
	doublereal *fp0, doublereal *fpold, doublereal *reducu, doublereal *
	reducv, doublereal *fpintu, doublereal *fpintv, integer *lastdi, 
	integer *nplusu, integer *nplusv, integer *nru, integer *nrv, integer 
	*nrdatu, integer *nrdatv, doublereal *wrk, integer *lwrk, integer *
	ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, l;
    static doublereal p, f1, f2, f3;
    static integer l1, l2, l3, l4;
    static doublereal p1, p2, p3, ub, vb, ue, ve;
    static integer mm, lq;
    static doublereal rn, acc;
    static integer laa;
    static doublereal one;
    static integer lau, lav, lbu, lbv, lri, nue, nve, mpm, nuk, lsu, lsv, nuu,
	     nvv, ich1, ich3;
    static doublereal con1;
    static integer lau1;
    static doublereal con4;
    static integer lav1;
    static doublereal con9;
    static integer npl1, nk1u, nk1v, ifbu, ifbv, ncof, iter;
    static doublereal fpms;
    static integer ifsu, ifsv;
    static doublereal peru, perv;
    static integer nplu, nplv, mvnu, nminu, nminv, nmaxu, nmaxv;
    extern /* Subroutine */ int fpgrpa_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpknot_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *);
    static integer nrintu, nrintv;

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpgrpa,fpknot */
/*  .. */
/*   set constants */
    /* Parameter adjustments */
    --ipar;
    --nru;
    --u;
    --nrv;
    --v;
    --z__;
    --nrdatu;
    --fpintu;
    --tu;
    --nrdatv;
    --fpintv;
    --tv;
    --c__;
    --wrk;

    /* Function Body */
    one = 1.;
    con1 = .1f;
    con9 = .9f;
    con4 = .04f;
/*  set boundaries of the approximation domain */
    ub = u[1];
    ue = u[*mu];
    vb = v[1];
    ve = v[*mv];
/*  we partition the working space. */
    lsu = 1;
    lsv = lsu + (*mu << 2);
    lri = lsv + (*mv << 2);
    mm = max(*nuest,*mv);
    lq = lri + mm * *idim;
    mvnu = *nuest * *mv * *idim;
    lau = lq + mvnu;
    nuk = *nuest * 5;
    lbu = lau + nuk;
    lav = lbu + nuk;
    nuk = *nvest * 5;
    lbv = lav + nuk;
    laa = lbv + nuk;
    lau1 = lau;
    if (ipar[1] == 0) {
	goto L10;
    }
    peru = ue - ub;
    lau1 = laa;
    laa += *nuest << 2;
L10:
    lav1 = lav;
    if (ipar[2] == 0) {
	goto L20;
    }
    perv = ve - vb;
    lav1 = laa;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 1: determination of the number of knots and their position.     c */
/* ****************************************************************     c */
/*  given a set of knots we compute the least-squares spline sinf(u,v), c */
/*  and the corresponding sum of squared residuals fp=f(p=inf).         c */
/*  if iopt=-1  sinf(u,v) is the requested approximation.               c */
/*  if iopt=0 or iopt=1 we check whether we can accept the knots:       c */
/*    if fp <=s we will continue with the current set of knots.         c */
/*    if fp > s we will increase the number of knots and compute the    c */
/*       corresponding least-squares spline until finally fp<=s.        c */
/*    the initial choice of knots depends on the value of s and iopt.   c */
/*    if s=0 we have spline interpolation; in that case the number of   c */
/*    knots equals nmaxu = mu+4+2*ipar(1) and  nmaxv = mv+4+2*ipar(2)   c */
/*    if s>0 and                                                        c */
/*     *iopt=0 we first compute the least-squares polynomial            c */
/*          nu=nminu=8 and nv=nminv=8                                   c */
/*     *iopt=1 we start with the knots found at the last call of the    c */
/*      routine, except for the case that s > fp0; then we can compute  c */
/*      the least-squares polynomial directly.                          c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  determine the number of knots for polynomial approximation. */
L20:
    nminu = 8;
    nminv = 8;
    if (*iopt < 0) {
	goto L100;
    }
/*  acc denotes the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
/*  find nmaxu and nmaxv which denote the number of knots in u- and v- */
/*  direction in case of spline interpolation. */
    nmaxu = *mu + 4 + (ipar[1] << 1);
    nmaxv = *mv + 4 + (ipar[2] << 1);
/*  find nue and nve which denote the maximum number of knots */
/*  allowed in each direction */
    nue = min(nmaxu,*nuest);
    nve = min(nmaxv,*nvest);
    if (*s > 0.f) {
	goto L60;
    }
/*  if s = 0, s(u,v) is an interpolating spline. */
    *nu = nmaxu;
    *nv = nmaxv;
/*  test whether the required storage space exceeds the available one. */
    if (*nv > *nvest || *nu > *nuest) {
	goto L420;
    }
/*  find the position of the interior knots in case of interpolation. */
/*  the knots in the u-direction. */
    nuu = *nu - 8;
    if (nuu == 0) {
	goto L40;
    }
    i__ = 5;
    j = 3 - ipar[1];
    i__1 = nuu;
    for (l = 1; l <= i__1; ++l) {
	tu[i__] = u[j];
	++i__;
	++j;
/* L30: */
    }
/*  the knots in the v-direction. */
L40:
    nvv = *nv - 8;
    if (nvv == 0) {
	goto L60;
    }
    i__ = 5;
    j = 3 - ipar[2];
    i__1 = nvv;
    for (l = 1; l <= i__1; ++l) {
	tv[i__] = v[j];
	++i__;
	++j;
/* L50: */
    }
    goto L100;
/*  if s > 0 our initial choice of knots depends on the value of iopt. */
L60:
    if (*iopt == 0) {
	goto L90;
    }
    if (*fp0 <= *s) {
	goto L90;
    }
/*  if iopt=1 and fp0 > s we start computing the least- squares spline */
/*  according to the set of knots found at the last call of the routine. */
/*  we determine the number of grid coordinates u(i) inside each knot */
/*  interval (tu(l),tu(l+1)). */
    l = 5;
    j = 1;
    nrdatu[1] = 0;
    mpm = *mu - 1;
    i__1 = mpm;
    for (i__ = 2; i__ <= i__1; ++i__) {
	++nrdatu[j];
	if (u[i__] < tu[l]) {
	    goto L70;
	}
	--nrdatu[j];
	++l;
	++j;
	nrdatu[j] = 0;
L70:
	;
    }
/*  we determine the number of grid coordinates v(i) inside each knot */
/*  interval (tv(l),tv(l+1)). */
    l = 5;
    j = 1;
    nrdatv[1] = 0;
    mpm = *mv - 1;
    i__1 = mpm;
    for (i__ = 2; i__ <= i__1; ++i__) {
	++nrdatv[j];
	if (v[i__] < tv[l]) {
	    goto L80;
	}
	--nrdatv[j];
	++l;
	++j;
	nrdatv[j] = 0;
L80:
	;
    }
    goto L100;
/*  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares */
/*  polynomial (which is a spline without interior knots). */
L90:
    *nu = nminu;
    *nv = nminv;
    nrdatu[1] = *mu - 2;
    nrdatv[1] = *mv - 2;
    *lastdi = 0;
    *nplusu = 0;
    *nplusv = 0;
    *fp0 = 0.f;
    *fpold = 0.f;
    *reducu = 0.f;
    *reducv = 0.f;
L100:
    mpm = *mu + *mv;
    ifsu = 0;
    ifsv = 0;
    ifbu = 0;
    ifbv = 0;
    p = -one;
/*  main loop for the different sets of knots.mpm=mu+mv is a save upper */
/*  bound for the number of trials. */
    i__1 = mpm;
    for (iter = 1; iter <= i__1; ++iter) {
	if (*nu == nminu && *nv == nminv) {
	    *ier = -2;
	}
/*  find nrintu (nrintv) which is the number of knot intervals in the */
/*  u-direction (v-direction). */
	nrintu = *nu - nminu + 1;
	nrintv = *nv - nminv + 1;
/*  find ncof, the number of b-spline coefficients for the current set */
/*  of knots. */
	nk1u = *nu - 4;
	nk1v = *nv - 4;
	ncof = nk1u * nk1v;
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(u,v). */
	if (ipar[1] != 0) {
	    goto L110;
	}
	i__ = *nu;
	for (j = 1; j <= 4; ++j) {
	    tu[j] = ub;
	    tu[i__] = ue;
	    --i__;
/* L105: */
	}
	goto L120;
L110:
	l1 = 4;
	l2 = l1;
	l3 = *nu - 3;
	l4 = l3;
	tu[l2] = ub;
	tu[l3] = ue;
	for (j = 1; j <= 3; ++j) {
	    ++l1;
	    --l2;
	    ++l3;
	    --l4;
	    tu[l2] = tu[l4] - peru;
	    tu[l3] = tu[l1] + peru;
/* L115: */
	}
L120:
	if (ipar[2] != 0) {
	    goto L130;
	}
	i__ = *nv;
	for (j = 1; j <= 4; ++j) {
	    tv[j] = vb;
	    tv[i__] = ve;
	    --i__;
/* L125: */
	}
	goto L140;
L130:
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
	    tv[l2] = tv[l4] - perv;
	    tv[l3] = tv[l1] + perv;
/* L135: */
	}
/*  find the least-squares spline sinf(u,v) and calculate for each knot */
/*  interval tu(j+3)<=u<=tu(j+4) (tv(j+3)<=v<=tv(j+4)) the sum */
/*  of squared residuals fpintu(j),j=1,2,...,nu-7 (fpintv(j),j=1,2,... */
/*  ,nv-7) for the data points having their absciss (ordinate)-value */
/*  belonging to that interval. */
/*  fp gives the total sum of squared residuals. */
L140:
	fpgrpa_(&ifsu, &ifsv, &ifbu, &ifbv, idim, &ipar[1], &u[1], mu, &v[1], 
		mv, &z__[1], mz, &tu[1], nu, &tv[1], nv, &p, &c__[1], nc, fp, 
		&fpintu[1], &fpintv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv], &
		wrk[lri], &wrk[lq], &wrk[lau], &wrk[lau1], &wrk[lav], &wrk[
		lav1], &wrk[lbu], &wrk[lbv], &nru[1], &nrv[1]);
	if (*ier == -2) {
	    *fp0 = *fp;
	}
/*  test whether the least-squares spline is an acceptable solution. */
	if (*iopt < 0) {
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
/*  if nu=nmaxu and nv=nmaxv, sinf(u,v) is an interpolating spline. */
	if (*nu == nmaxu && *nv == nmaxv) {
	    goto L430;
	}
/*  increase the number of knots. */
/*  if nu=nue and nv=nve we cannot further increase the number of knots */
/*  because of the storage capacity limitation. */
	if (*nu == nue && *nv == nve) {
	    goto L420;
	}
	*ier = 0;
/*  adjust the parameter reducu or reducv according to the direction */
/*  in which the last added knots were located. */
	if (*lastdi < 0) {
	    goto L150;
	}
	if (*lastdi == 0) {
	    goto L170;
	}
	goto L160;
L150:
	*reducu = *fpold - *fp;
	goto L170;
L160:
	*reducv = *fpold - *fp;
/*  store the sum of squared residuals for the current set of knots. */
L170:
	*fpold = *fp;
/*  find nplu, the number of knots we should add in the u-direction. */
	nplu = 1;
	if (*nu == nminu) {
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
	nplv = 1;
	if (*nv == nminv) {
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
	i__2 = *nplusu;
	for (l = 1; l <= i__2; ++l) {
/*  add a new knot in the u-direction */
	    fpknot_(&u[1], mu, &tu[1], nu, &fpintu[1], &nrdatu[1], &nrintu, 
		    nuest, &c__1);
/*  test whether we cannot further increase the number of knots in the */
/*  u-direction. */
	    if (*nu == nue) {
		goto L250;
	    }
/* L220: */
	}
	goto L250;
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
		goto L250;
	    }
/* L240: */
	}
/*  restart the computations with the new set of knots. */
L250:
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
/*  this smoothing spline varies with the parameter p in such a way thatc */
/*  f(p)=suml=1,idim(sumi=1,mu(sumj=1,mv((z(i,j,l)-sp(u(i),v(j),l))**2) c */
/*  is a continuous, strictly decreasing function of p. moreover the    c */
/*  least-squares polynomial corresponds to p=0 and the least-squares   c */
/*  spline to p=infinity. iteratively we then have to determine the     c */
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
    ich1 = 0;
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the smoothing spline sp(u,v) and the corresponding sum of */
/*  squared residuals fp. */
	fpgrpa_(&ifsu, &ifsv, &ifbu, &ifbv, idim, &ipar[1], &u[1], mu, &v[1], 
		mv, &z__[1], mz, &tu[1], nu, &tv[1], nv, &p, &c__[1], nc, fp, 
		&fpintu[1], &fpintv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv], &
		wrk[lri], &wrk[lq], &wrk[lau], &wrk[lau1], &wrk[lav], &wrk[
		lav1], &wrk[lbu], &wrk[lbv], &nru[1], &nrv[1]);
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
} /* fppasu_ */

