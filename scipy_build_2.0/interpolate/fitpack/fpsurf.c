/* fpsurf.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpsurf_(integer *iopt, integer *m, doublereal *x, 
	doublereal *y, doublereal *z__, doublereal *w, doublereal *xb, 
	doublereal *xe, doublereal *yb, doublereal *ye, integer *kxx, integer 
	*kyy, doublereal *s, integer *nxest, integer *nyest, doublereal *eta, 
	doublereal *tol, integer *maxit, integer *nmax, integer *km1, integer 
	*km2, integer *ib1, integer *ib3, integer *nc, integer *intest, 
	integer *nrest, integer *nx0, doublereal *tx, integer *ny0, 
	doublereal *ty, doublereal *c__, doublereal *fp, doublereal *fp0, 
	doublereal *fpint, doublereal *coord, doublereal *f, doublereal *ff, 
	doublereal *a, doublereal *q, doublereal *bx, doublereal *by, 
	doublereal *spx, doublereal *spy, doublereal *h__, integer *index, 
	integer *nummer, doublereal *wrk, integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, bx_dim1, bx_offset, by_dim1, 
	    by_offset, spx_dim1, spx_offset, spy_dim1, spy_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l, n;
    static doublereal p, f1, f2, f3;
    static integer i1, i2, i3, j1, l1, l2, n1;
    static doublereal p1, p2, p3, x0, x1, y0, y1;
    static integer la, ii, lf, lh, in;
    static doublereal wi, rn, hx[6], zi, sq;
    static integer kx, ky, lx, ly, nx, ny;
    static doublereal hy[6];
    static integer kx1, kx2, ky1, ky2;
    static doublereal acc;
    static integer ibb;
    static doublereal arg, one, cos__, ten, eps, hxi, sin__;
    static integer nxe, nye;
    static doublereal piv;
    static integer num;
    static doublereal fac1, fac2;
    static integer jxy, nxx, nyy, ich1, ich3;
    static doublereal con1, con4, con9;
    static integer num1, nk1x, nk1y;
    static doublereal half;
    static integer ncof;
    static doublereal dmax__;
    static integer nreg, rank, iter;
    static doublereal fpms, pinv;
    static integer irot, jrot, iband;
    static doublereal sigma, fpmax;
    static integer nminx, nminy;
    static doublereal store;
    static integer nrint, iband1, lwest, iband3, iband4;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static integer ichang;
    extern /* Subroutine */ int fpdisc_(doublereal *, integer *, integer *, 
	    doublereal *, integer *), fporde_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *), 
	    fprank_(doublereal *, doublereal *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), fprota_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpgivs_(doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --nummer;
    --w;
    --z__;
    --y;
    --x;
    --ty;
    --tx;
    spy_dim1 = *m;
    spy_offset = 1 + spy_dim1;
    spy -= spy_offset;
    spx_dim1 = *m;
    spx_offset = 1 + spx_dim1;
    spx -= spx_offset;
    by_dim1 = *nmax;
    by_offset = 1 + by_dim1;
    by -= by_offset;
    bx_dim1 = *nmax;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
    --h__;
    q_dim1 = *nc;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *nc;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ff;
    --f;
    --c__;
    --coord;
    --fpint;
    --index;
    --wrk;

    /* Function Body */
    one = 1.f;
    con1 = .1f;
    con9 = .9f;
    con4 = .04f;
    half = .5f;
    ten = 10.f;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 1: determination of the number of knots and their position.     c */
/* ****************************************************************     c */
/* given a set of knots we compute the least-squares spline sinf(x,y),  c */
/* and the corresponding weighted sum of squared residuals fp=f(p=inf). c */
/* if iopt=-1  sinf(x,y) is the requested approximation.                c */
/* if iopt=0 or iopt=1 we check whether we can accept the knots:        c */
/*   if fp <=s we will continue with the current set of knots.          c */
/*   if fp > s we will increase the number of knots and compute the     c */
/*      corresponding least-squares spline until finally  fp<=s.        c */
/* the initial choice of knots depends on the value of s and iopt.      c */
/*   if iopt=0 we first compute the least-squares polynomial of degree  c */
/*     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c */
/*     fp0=f(0) denotes the corresponding weighted sum of squared       c */
/*     residuals                                                        c */
/*   if iopt=1 we start with the knots found at the last call of the    c */
/*     routine, except for the case that s>=fp0; then we can compute    c */
/*     the least-squares polynomial directly.                           c */
/* eventually the independent variables x and y (and the corresponding  c */
/* parameters) will be switched if this can reduce the bandwidth of the c */
/* system to be solved.                                                 c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  ichang denotes whether(1) or not(-1) the directions have been inter- */
/*  changed. */
    ichang = -1;
    x0 = *xb;
    x1 = *xe;
    y0 = *yb;
    y1 = *ye;
    kx = *kxx;
    ky = *kyy;
    kx1 = kx + 1;
    ky1 = ky + 1;
    nxe = *nxest;
    nye = *nyest;
    eps = sqrt(*eta);
    if (*iopt < 0) {
	goto L20;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
    if (*iopt == 0) {
	goto L10;
    }
    if (*fp0 > *s) {
	goto L20;
    }
/*  initialization for the least-squares polynomial. */
L10:
    nminx = kx1 << 1;
    nminy = ky1 << 1;
    nx = nminx;
    ny = nminy;
    *ier = -2;
    goto L30;
L20:
    nx = *nx0;
    ny = *ny0;
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
L30:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(x,y). */
	l = nx;
	i__2 = kx1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tx[i__] = x0;
	    tx[l] = x1;
	    --l;
/* L40: */
	}
	l = ny;
	i__2 = ky1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ty[i__] = y0;
	    ty[l] = y1;
	    --l;
/* L50: */
	}
/*  find nrint, the total number of knot intervals and nreg, the number */
/*  of panels in which the approximation domain is subdivided by the */
/*  intersection of knots. */
	nxx = nx - (kx1 << 1) + 1;
	nyy = ny - (ky1 << 1) + 1;
	nrint = nxx + nyy;
	nreg = nxx * nyy;
/*  find the bandwidth of the observation matrix a. */
/*  if necessary, interchange the variables x and y, in order to obtain */
/*  a minimal bandwidth. */
	iband1 = kx * (ny - ky1) + ky;
	l = ky * (nx - kx1) + kx;
	if (iband1 <= l) {
	    goto L130;
	}
	iband1 = l;
	ichang = -ichang;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    store = x[i__];
	    x[i__] = y[i__];
	    y[i__] = store;
/* L60: */
	}
	store = x0;
	x0 = y0;
	y0 = store;
	store = x1;
	x1 = y1;
	y1 = store;
	n = min(nx,ny);
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    store = tx[i__];
	    tx[i__] = ty[i__];
	    ty[i__] = store;
/* L70: */
	}
	n1 = n + 1;
	if (nx < ny) {
	    goto L80;
	}
	if (nx == ny) {
	    goto L120;
	}
	goto L100;
L80:
	i__2 = ny;
	for (i__ = n1; i__ <= i__2; ++i__) {
	    tx[i__] = ty[i__];
/* L90: */
	}
	goto L120;
L100:
	i__2 = nx;
	for (i__ = n1; i__ <= i__2; ++i__) {
	    ty[i__] = tx[i__];
/* L110: */
	}
L120:
	l = nx;
	nx = ny;
	ny = l;
	l = nxe;
	nxe = nye;
	nye = l;
	l = nxx;
	nxx = nyy;
	nyy = l;
	l = kx;
	kx = ky;
	ky = l;
	kx1 = kx + 1;
	ky1 = ky + 1;
L130:
	iband = iband1 + 1;
/*  arrange the data points according to the panel they belong to. */
	fporde_(&x[1], &y[1], m, &kx, &ky, &tx[1], &nx, &ty[1], &ny, &nummer[
		1], &index[1], &nreg);
/*  find ncof, the number of b-spline coefficients. */
	nk1x = nx - kx1;
	nk1y = ny - ky1;
	ncof = nk1x * nk1y;
/*  initialize the observation matrix a. */
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[i__] = 0.f;
	    i__3 = iband;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * a_dim1] = 0.f;
/* L140: */
	    }
	}
/*  initialize the sum of squared residuals. */
	*fp = 0.f;
/*  fetch the data points in the new order. main loop for the */
/*  different panels. */
	i__3 = nreg;
	for (num = 1; num <= i__3; ++num) {
/*  fix certain constants for the current panel; jrot records the column */
/*  number of the first non-zero element in a row of the observation */
/*  matrix according to a data point of the panel. */
	    num1 = num - 1;
	    lx = num1 / nyy;
	    l1 = lx + kx1;
	    ly = num1 - lx * nyy;
	    l2 = ly + ky1;
	    jrot = lx * nk1y + ly;
/*  test whether there are still data points in the panel. */
	    in = index[num];
L150:
	    if (in == 0) {
		goto L250;
	    }
/*  fetch a new data point. */
	    wi = w[in];
	    zi = z__[in] * wi;
/*  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in). */
	    fpbspl_(&tx[1], &nx, &kx, &x[in], &l1, hx);
/*  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in). */
	    fpbspl_(&ty[1], &ny, &ky, &y[in], &l2, hy);
/*  store the value of these b-splines in spx and spy respectively. */
	    i__2 = kx1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		spx[in + i__ * spx_dim1] = hx[i__ - 1];
/* L160: */
	    }
	    i__2 = ky1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		spy[in + i__ * spy_dim1] = hy[i__ - 1];
/* L170: */
	    }
/*  initialize the new row of observation matrix. */
	    i__2 = iband;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		h__[i__] = 0.f;
/* L180: */
	    }
/*  calculate the non-zero elements of the new row by making the cross */
/*  products of the non-zero b-splines in x- and y-direction. */
	    i1 = 0;
	    i__2 = kx1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		hxi = hx[i__ - 1];
		j1 = i1;
		i__4 = ky1;
		for (j = 1; j <= i__4; ++j) {
		    ++j1;
		    h__[j1] = hxi * hy[j - 1] * wi;
/* L190: */
		}
		i1 += nk1y;
/* L200: */
	    }
/*  rotate the row into triangle by givens transformations . */
	    irot = jrot;
	    i__2 = iband;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++irot;
		piv = h__[i__];
		if (piv == 0.f) {
		    goto L220;
		}
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &a[irot + a_dim1], &cos__, &sin__);
/*  apply that transformation to the right hand side. */
		fprota_(&cos__, &sin__, &zi, &f[irot]);
		if (i__ == iband) {
		    goto L230;
		}
/*  apply that transformation to the left hand side. */
		i2 = 1;
		i3 = i__ + 1;
		i__4 = iband;
		for (j = i3; j <= i__4; ++j) {
		    ++i2;
		    fprota_(&cos__, &sin__, &h__[j], &a[irot + i2 * a_dim1]);
/* L210: */
		}
L220:
		;
	    }
/*  add the contribution of the row to the sum of squares of residual */
/*  right hand sides. */
L230:
/* Computing 2nd power */
	    d__1 = zi;
	    *fp += d__1 * d__1;
/*  find the number of the next data point in the panel. */
	    in = nummer[in];
	    goto L150;
L250:
	    ;
	}
/*  find dmax, the maximum value for the diagonal elements in the reduced */
/*  triangle. */
	dmax__ = 0.f;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (a[i__ + a_dim1] <= dmax__) {
		goto L260;
	    }
	    dmax__ = a[i__ + a_dim1];
L260:
	    ;
	}
/*  check whether the observation matrix is rank deficient. */
	sigma = eps * dmax__;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (a[i__ + a_dim1] <= sigma) {
		goto L280;
	    }
/* L270: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&a[a_offset], &f[1], &ncof, &iband, &c__[1], nc);
	rank = ncof;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + q_dim1] = a[i__ + a_dim1] / dmax__;
/* L275: */
	}
	goto L300;
/*  in case of rank deficiency, find the minimum norm solution. */
/*  check whether there is sufficient working space */
L280:
	lwest = ncof * iband + ncof + iband;
	if (*lwrk < lwest) {
	    goto L780;
	}
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ff[i__] = f[i__];
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L290: */
	    }
	}
	lf = 1;
	lh = lf + ncof;
	la = lh + iband;
	fprank_(&q[q_offset], &ff[1], &ncof, &iband, nc, &sigma, &c__[1], &sq,
		 &rank, &wrk[la], &wrk[lf], &wrk[lh]);
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L295: */
	}
/*  add to the sum of squared residuals, the contribution of reducing */
/*  the rank. */
	*fp += sq;
L300:
	if (*ier == -2) {
	    *fp0 = *fp;
	}
/*  test whether the least-squares spline is an acceptable solution. */
	if (*iopt < 0) {
	    goto L820;
	}
	fpms = *fp - *s;
	if (abs(fpms) <= acc) {
	    if (*fp <= 0.) {
		goto L815;
	    }
	    goto L820;
	}
/*  test whether we can accept the choice of knots. */
	if (fpms < 0.f) {
	    goto L430;
	}
/*  test whether we cannot further increase the number of knots. */
	if (ncof > *m) {
	    goto L790;
	}
	*ier = 0;
/*  search where to add a new knot. */
/*  find for each interval the sum of squared residuals fpint for the */
/*  data points having the coordinate belonging to that knot interval. */
/*  calculate also coord which is the same sum, weighted by the position */
/*  of the data points considered. */
	i__2 = nrint;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fpint[i__] = 0.f;
	    coord[i__] = 0.f;
/* L320: */
	}
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
	    num1 = num - 1;
	    lx = num1 / nyy;
	    l1 = lx + 1;
	    ly = num1 - lx * nyy;
	    l2 = ly + 1 + nxx;
	    jrot = lx * nk1y + ly;
	    in = index[num];
L330:
	    if (in == 0) {
		goto L360;
	    }
	    store = 0.f;
	    i1 = jrot;
	    i__3 = kx1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		hxi = spx[in + i__ * spx_dim1];
		j1 = i1;
		i__4 = ky1;
		for (j = 1; j <= i__4; ++j) {
		    ++j1;
		    store += hxi * spy[in + j * spy_dim1] * c__[j1];
/* L340: */
		}
		i1 += nk1y;
/* L350: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    store = d__1 * d__1;
	    fpint[l1] += store;
	    coord[l1] += store * x[in];
	    fpint[l2] += store;
	    coord[l2] += store * y[in];
	    in = nummer[in];
	    goto L330;
L360:
	    ;
	}
/*  find the interval for which fpint is maximal on the condition that */
/*  there still can be added a knot. */
L370:
	l = 0;
	fpmax = 0.f;
	l1 = 1;
	l2 = nrint;
	if (nx == nxe) {
	    l1 = nxx + 1;
	}
	if (ny == nye) {
	    l2 = nxx;
	}
	if (l1 > l2) {
	    goto L810;
	}
	i__2 = l2;
	for (i__ = l1; i__ <= i__2; ++i__) {
	    if (fpmax >= fpint[i__]) {
		goto L380;
	    }
	    l = i__;
	    fpmax = fpint[i__];
L380:
	    ;
	}
/*  test whether we cannot further increase the number of knots. */
	if (l == 0) {
	    goto L785;
	}
/*  calculate the position of the new knot. */
	arg = coord[l] / fpint[l];
/*  test in what direction the new knot is going to be added. */
	if (l > nxx) {
	    goto L400;
	}
/*  addition in the x-direction. */
	jxy = l + kx1;
	fpint[l] = 0.f;
	fac1 = tx[jxy] - arg;
	fac2 = arg - tx[jxy - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L370;
	}
	j = nx;
	i__2 = nx;
	for (i__ = jxy; i__ <= i__2; ++i__) {
	    tx[j + 1] = tx[j];
	    --j;
/* L390: */
	}
	tx[jxy] = arg;
	++nx;
	goto L420;
/*  addition in the y-direction. */
L400:
	jxy = l + ky1 - nxx;
	fpint[l] = 0.f;
	fac1 = ty[jxy] - arg;
	fac2 = arg - ty[jxy - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L370;
	}
	j = ny;
	i__2 = ny;
	for (i__ = jxy; i__ <= i__2; ++i__) {
	    ty[j + 1] = ty[j];
	    --j;
/* L410: */
	}
	ty[jxy] = arg;
	++ny;
/*  restart the computations with the new set of knots. */
L420:
	;
    }
/*  test whether the least-squares polynomial is a solution of our */
/*  approximation problem. */
L430:
    if (*ier == -2) {
	goto L830;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing spline sp(x,y)                c */
/* *****************************************************                c */
/* we have determined the number of knots and their position. we now    c */
/* compute the b-spline coefficients of the smoothing spline sp(x,y).   c */
/* the observation matrix a is extended by the rows of a matrix,        c */
/* expressing that sp(x,y) must be a polynomial of degree kx in x and   c */
/* ky in y. the corresponding weights of these additional rows are set  c */
/* to 1./p.  iteratively we than have to determine the value of p       c */
/* such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c */
/* we already know that the least-squares polynomial corresponds to     c */
/* p=0  and that the least-squares spline corresponds to p=infinity.    c */
/* the iteration process which is proposed here makes use of rational   c */
/* interpolation. since f(p) is a convex and strictly decreasing        c */
/* function of p, it can be approximated by a rational function r(p)=   c */
/* (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c */
/* of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c */
/* new value of p such that r(p)=s. convergence is guaranteed by taking c */
/* f1 > 0 and f3 < 0.                                                   c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    kx2 = kx1 + 1;
/*  test whether there are interior knots in the x-direction. */
    if (nk1x == kx1) {
	goto L440;
    }
/*  evaluate the discotinuity jumps of the kx-th order derivative of */
/*  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1. */
    fpdisc_(&tx[1], &nx, &kx2, &bx[bx_offset], nmax);
L440:
    ky2 = ky1 + 1;
/*  test whether there are interior knots in the y-direction. */
    if (nk1y == ky1) {
	goto L450;
    }
/*  evaluate the discontinuity jumps of the ky-th order derivative of */
/*  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1. */
    fpdisc_(&ty[1], &ny, &ky2, &by[by_offset], nmax);
/*  initial value for p. */
L450:
    p1 = 0.f;
    f1 = *fp0 - *s;
    p3 = -one;
    f3 = fpms;
    p = 0.f;
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p += a[i__ + a_dim1];
/* L460: */
    }
    rn = (doublereal) ncof;
    p = rn / p;
/*  find the bandwidth of the extended observation matrix. */
    iband3 = kx1 * nk1y;
    iband4 = iband3 + 1;
    ich1 = 0;
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	pinv = one / p;
/*  store the triangularized observation matrix into q. */
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ff[i__] = f[i__];
	    i__3 = iband;
	    for (j = 1; j <= i__3; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L470: */
	    }
	    ibb = iband + 1;
	    i__3 = iband4;
	    for (j = ibb; j <= i__3; ++j) {
		q[i__ + j * q_dim1] = 0.f;
/* L480: */
	    }
	}
	if (nk1y == ky1) {
	    goto L560;
	}
/*  extend the observation matrix with the rows of a matrix, expressing */
/*  that for x=cst. sp(x,y) must be a polynomial in y of degree ky. */
	i__3 = nk1y;
	for (i__ = ky2; i__ <= i__3; ++i__) {
	    ii = i__ - ky1;
	    i__2 = nk1x;
	    for (j = 1; j <= i__2; ++j) {
/*  initialize the new row. */
		i__4 = iband;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = 0.f;
/* L490: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		i__4 = ky2;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = by[ii + l * by_dim1] * pinv;
/* L500: */
		}
		zi = 0.f;
		jrot = (j - 1) * nk1y + ii;
/*  rotate the new row into triangle by givens transformations without */
/*  square roots. */
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
		    piv = h__[1];
/* Computing MIN */
		    i__5 = iband1, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
		    if (piv == 0.f) {
			if (i2 <= 0) {
			    goto L550;
			}
			goto L520;
		    }
/*  calculate the parameters of the givens transformation. */
		    fpgivs_(&piv, &q[irot + q_dim1], &cos__, &sin__);
/*  apply that givens transformation to the right hand side. */
		    fprota_(&cos__, &sin__, &zi, &ff[irot]);
		    if (i2 == 0) {
			goto L550;
		    }
/*  apply that givens transformation to the left hand side. */
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			l1 = l + 1;
			fprota_(&cos__, &sin__, &h__[l1], &q[irot + l1 * 
				q_dim1]);
/* L510: */
		    }
L520:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			h__[l] = h__[l + 1];
/* L530: */
		    }
		    h__[i2 + 1] = 0.f;
/* L540: */
		}
L550:
		;
	    }
	}
L560:
	if (nk1x == kx1) {
	    goto L640;
	}
/*  extend the observation matrix with the rows of a matrix expressing */
/*  that for y=cst. sp(x,y) must be a polynomial in x of degree kx. */
	i__2 = nk1x;
	for (i__ = kx2; i__ <= i__2; ++i__) {
	    ii = i__ - kx1;
	    i__3 = nk1y;
	    for (j = 1; j <= i__3; ++j) {
/*  initialize the new row */
		i__4 = iband4;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = 0.f;
/* L570: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		j1 = 1;
		i__4 = kx2;
		for (l = 1; l <= i__4; ++l) {
		    h__[j1] = bx[ii + l * bx_dim1] * pinv;
		    j1 += nk1y;
/* L580: */
		}
		zi = 0.f;
		jrot = (i__ - kx2) * nk1y + j;
/*  rotate the new row into triangle by givens transformations . */
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
		    piv = h__[1];
/* Computing MIN */
		    i__5 = iband3, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
		    if (piv == 0.f) {
			if (i2 <= 0) {
			    goto L630;
			}
			goto L600;
		    }
/*  calculate the parameters of the givens transformation. */
		    fpgivs_(&piv, &q[irot + q_dim1], &cos__, &sin__);
/*  apply that givens transformation to the right hand side. */
		    fprota_(&cos__, &sin__, &zi, &ff[irot]);
		    if (i2 == 0) {
			goto L630;
		    }
/*  apply that givens transformation to the left hand side. */
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			l1 = l + 1;
			fprota_(&cos__, &sin__, &h__[l1], &q[irot + l1 * 
				q_dim1]);
/* L590: */
		    }
L600:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			h__[l] = h__[l + 1];
/* L610: */
		    }
		    h__[i2 + 1] = 0.f;
/* L620: */
		}
L630:
		;
	    }
	}
/*  find dmax, the maximum value for the diagonal elements in the */
/*  reduced triangle. */
L640:
	dmax__ = 0.f;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (q[i__ + q_dim1] <= dmax__) {
		goto L650;
	    }
	    dmax__ = q[i__ + q_dim1];
L650:
	    ;
	}
/*  check whether the matrix is rank deficient. */
	sigma = eps * dmax__;
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (q[i__ + q_dim1] <= sigma) {
		goto L670;
	    }
/* L660: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&q[q_offset], &ff[1], &ncof, &iband4, &c__[1], nc);
	rank = ncof;
	goto L675;
/*  in case of rank deficiency, find the minimum norm solution. */
L670:
	lwest = ncof * iband4 + ncof + iband4;
	if (*lwrk < lwest) {
	    goto L780;
	}
	lf = 1;
	lh = lf + ncof;
	la = lh + iband4;
	fprank_(&q[q_offset], &ff[1], &ncof, &iband4, nc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
L675:
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L680: */
	}
/*  compute f(p). */
	*fp = 0.f;
	i__3 = nreg;
	for (num = 1; num <= i__3; ++num) {
	    num1 = num - 1;
	    lx = num1 / nyy;
	    ly = num1 - lx * nyy;
	    jrot = lx * nk1y + ly;
	    in = index[num];
L690:
	    if (in == 0) {
		goto L720;
	    }
	    store = 0.f;
	    i1 = jrot;
	    i__2 = kx1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		hxi = spx[in + i__ * spx_dim1];
		j1 = i1;
		i__4 = ky1;
		for (j = 1; j <= i__4; ++j) {
		    ++j1;
		    store += hxi * spy[in + j * spy_dim1] * c__[j1];
/* L700: */
		}
		i1 += nk1y;
/* L710: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    *fp += d__1 * d__1;
	    in = nummer[in];
	    goto L690;
L720:
	    ;
	}
/*  test whether the approximation sp(x,y) is an acceptable solution. */
	fpms = *fp - *s;
	if (abs(fpms) <= acc) {
	    goto L820;
	}
/*  test whether the maximum allowable number of iterations has been */
/*  reached. */
	if (iter == *maxit) {
	    goto L795;
	}
/*  carry out one more step of the iteration process. */
	p2 = p;
	f2 = fpms;
	if (ich3 != 0) {
	    goto L740;
	}
	if (f2 - f3 > acc) {
	    goto L730;
	}
/*  our initial choice of p is too large. */
	p3 = p2;
	f3 = f2;
	p *= con4;
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
	goto L770;
L730:
	if (f2 < 0.f) {
	    ich3 = 1;
	}
L740:
	if (ich1 != 0) {
	    goto L760;
	}
	if (f1 - f2 > acc) {
	    goto L750;
	}
/*  our initial choice of p is too small */
	p1 = p2;
	f1 = f2;
	p /= con4;
	if (p3 < 0.f) {
	    goto L770;
	}
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
	goto L770;
L750:
	if (f2 > 0.f) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
L760:
	if (f2 >= f1 || f2 <= f3) {
	    goto L800;
	}
/*  find the new value of p. */
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
L770:
	;
    }
/*  error codes and messages. */
L780:
    *ier = lwest;
    goto L830;
L785:
    *ier = 5;
    goto L830;
L790:
    *ier = 4;
    goto L830;
L795:
    *ier = 3;
    goto L830;
L800:
    *ier = 2;
    goto L830;
L810:
    *ier = 1;
    goto L830;
L815:
    *ier = -1;
    *fp = 0.f;
L820:
    if (ncof != rank) {
	*ier = -rank;
    }
/*  test whether x and y are in the original order. */
L830:
    if (ichang < 0) {
	goto L930;
    }
/*  if not, interchange x and y once more. */
    l1 = 1;
    i__1 = nk1x;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2 = i__;
	i__3 = nk1y;
	for (j = 1; j <= i__3; ++j) {
	    f[l2] = c__[l1];
	    ++l1;
	    l2 += nk1x;
/* L840: */
	}
    }
    i__3 = ncof;
    for (i__ = 1; i__ <= i__3; ++i__) {
	c__[i__] = f[i__];
/* L850: */
    }
    i__3 = *m;
    for (i__ = 1; i__ <= i__3; ++i__) {
	store = x[i__];
	x[i__] = y[i__];
	y[i__] = store;
/* L860: */
    }
    n = min(nx,ny);
    i__3 = n;
    for (i__ = 1; i__ <= i__3; ++i__) {
	store = tx[i__];
	tx[i__] = ty[i__];
	ty[i__] = store;
/* L870: */
    }
    n1 = n + 1;
    if (nx < ny) {
	goto L880;
    }
    if (nx == ny) {
	goto L920;
    }
    goto L900;
L880:
    i__3 = ny;
    for (i__ = n1; i__ <= i__3; ++i__) {
	tx[i__] = ty[i__];
/* L890: */
    }
    goto L920;
L900:
    i__3 = nx;
    for (i__ = n1; i__ <= i__3; ++i__) {
	ty[i__] = tx[i__];
/* L910: */
    }
L920:
    l = nx;
    nx = ny;
    ny = l;
L930:
    if (*iopt < 0) {
	goto L940;
    }
    *nx0 = nx;
    *ny0 = ny;
L940:
    return 0;
} /* fpsurf_ */

