/* fpsphe.f -- translated by f2c (version 20190311).
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
static integer c__5 = 5;

/* Subroutine */ int fpsphe_(integer *iopt, integer *m, doublereal *teta, 
	doublereal *phi, doublereal *r__, doublereal *w, doublereal *s, 
	integer *ntest, integer *npest, doublereal *eta, doublereal *tol, 
	integer *maxit, integer *ib1, integer *ib3, integer *nc, integer *ncc,
	 integer *intest, integer *nrest, integer *nt, doublereal *tt, 
	integer *np, doublereal *tp, doublereal *c__, doublereal *fp, 
	doublereal *sup, doublereal *fpint, doublereal *coord, doublereal *f, 
	doublereal *ff, doublereal *row, doublereal *coco, doublereal *cosi, 
	doublereal *a, doublereal *q, doublereal *bt, doublereal *bp, 
	doublereal *spt, doublereal *spp, doublereal *h__, integer *index, 
	integer *nummer, doublereal *wrk, integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, bt_dim1, bt_offset, bp_dim1, 
	    bp_offset, spt_dim1, spt_offset, spp_dim1, spp_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static integer i__, j, l;
    static doublereal p, c1, d1, d2, f1, f2, f3;
    static integer i1, i2, i3, j1, j2, l1, l2;
    static doublereal p1, p2, p3;
    static integer l3, l4;
    static doublereal aa;
    static integer la;
    static doublereal cn, co, fn;
    static integer ii;
    static doublereal pi;
    static integer ij;
    static doublereal ri, si;
    static integer il, in;
    static doublereal wi, rn;
    static integer lf;
    static doublereal sq;
    static integer lh, ll, lp, lt;
    static doublereal ht[4], hp[4], pi2;
    static integer nr1, np4, nt4, nt6;
    static doublereal acc, arg, one, hti, htj, ten, eps;
    static integer jlt, npp;
    static doublereal piv;
    static integer num, nrr, ntt;
    static doublereal fac1, fac2;
    static integer ich1, ich3;
    static doublereal con1, con4, con9;
    static integer num1;
    static doublereal facc, half, facs;
    static integer ncof;
    static doublereal dmax__;
    static integer nreg, rank, iter;
    static doublereal fpms, pinv;
    static integer irot, jrot, iband, ncoff;
    static doublereal sigma, fpmax;
    static integer nrint;
    static doublereal store;
    static integer iband1, lwest, iband3, iband4;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), fpdisc_(doublereal *, 
	    integer *, integer *, doublereal *, integer *), fporde_(
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *), fprank_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    ;
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), fprota_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpgivs_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fprpsp_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *);

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*   fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota,fprpsp */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --nummer;
    spp_dim1 = *m;
    spp_offset = 1 + spp_dim1;
    spp -= spp_offset;
    spt_dim1 = *m;
    spt_offset = 1 + spt_dim1;
    spt -= spt_offset;
    --w;
    --r__;
    --phi;
    --teta;
    bt_dim1 = *ntest;
    bt_offset = 1 + bt_dim1;
    bt -= bt_offset;
    --tt;
    bp_dim1 = *npest;
    bp_offset = 1 + bp_dim1;
    bp -= bp_offset;
    --cosi;
    --coco;
    --row;
    --tp;
    --h__;
    --ff;
    --c__;
    q_dim1 = *ncc;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *ncc;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --f;
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
    pi = atan(one) * 4;
    pi2 = pi + pi;
    eps = sqrt(*eta);
    if (*iopt < 0) {
	goto L70;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
    if (*iopt == 0) {
	goto L10;
    }
    if (*s < *sup) {
	if (*np < 11) {
	    goto L60;
	}
	goto L70;
    }
/*  if iopt=0 we begin by computing the weighted least-squares polynomial */
/*  of the form */
/*     s(teta,phi) = c1*f1(teta) + cn*fn(teta) */
/*  where f1(teta) and fn(teta) are the cubic polynomials satisfying */
/*     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta). */
/*  the corresponding weighted sum of squared residuals gives the upper */
/*  bound sup for the smoothing factor s. */
L10:
    *sup = 0.f;
    d1 = 0.f;
    d2 = 0.f;
    c1 = 0.f;
    cn = 0.f;
    fac1 = pi * (one + half);
/* Computing 3rd power */
    d__1 = pi;
    fac2 = (one + one) / (d__1 * (d__1 * d__1));
    aa = 0.f;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__];
	ri = r__[i__] * wi;
	arg = teta[i__];
	fn = fac2 * arg * arg * (fac1 - arg);
	f1 = (one - fn) * wi;
	fn *= wi;
	if (fn == 0.f) {
	    goto L20;
	}
	fpgivs_(&fn, &d1, &co, &si);
	fprota_(&co, &si, &f1, &aa);
	fprota_(&co, &si, &ri, &cn);
L20:
	if (f1 == 0.f) {
	    goto L30;
	}
	fpgivs_(&f1, &d2, &co, &si);
	fprota_(&co, &si, &ri, &c1);
L30:
	*sup += ri * ri;
/* L40: */
    }
    if (d2 != 0.f) {
	c1 /= d2;
    }
    if (d1 != 0.f) {
	cn = (cn - aa * c1) / d1;
    }
/*  find the b-spline representation of this least-squares polynomial */
    *nt = 8;
    *np = 8;
    for (i__ = 1; i__ <= 4; ++i__) {
	c__[i__] = c1;
	c__[i__ + 4] = c1;
	c__[i__ + 8] = cn;
	c__[i__ + 12] = cn;
	tt[i__] = 0.f;
	tt[i__ + 4] = pi;
	tp[i__] = 0.f;
	tp[i__ + 4] = pi2;
/* L50: */
    }
    *fp = *sup;
/*  test whether the least-squares polynomial is an acceptable solution */
    fpms = *sup - *s;
    if (fpms < acc) {
	goto L960;
    }
/*  test whether we cannot further increase the number of knots. */
L60:
    if (*npest < 11 || *ntest < 9) {
	goto L950;
    }
/*  find the initial set of interior knots of the spherical spline in */
/*  case iopt = 0. */
    *np = 11;
    tp[5] = pi * half;
    tp[6] = pi;
    tp[7] = tp[5] + pi;
    *nt = 9;
    tt[5] = tp[5];
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 1 : computation of least-squares spherical splines.            c */
/*  ********************************************************            c */
/*  if iopt < 0 we compute the least-squares spherical spline according c */
/*  to the given set of knots.                                          c */
/*  if iopt >=0 we compute least-squares spherical splines with increas-c */
/*  ing numbers of knots until the corresponding sum f(p=inf)<=s.       c */
/*  the initial set of knots then depends on the value of iopt:         c */
/*    if iopt=0 we start with one interior knot in the teta-direction   c */
/*              (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). c */
/*    if iopt>0 we start with the set of knots found at the last call   c */
/*              of the routine.                                         c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
L70:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(teta,phi). */
	l1 = 4;
	l2 = l1;
	l3 = *np - 3;
	l4 = l3;
	tp[l2] = 0.f;
	tp[l3] = pi2;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ++l1;
	    --l2;
	    ++l3;
	    --l4;
	    tp[l2] = tp[l4] - pi2;
	    tp[l3] = tp[l1] + pi2;
/* L80: */
	}
	l = *nt;
	for (i__ = 1; i__ <= 4; ++i__) {
	    tt[i__] = 0.f;
	    tt[l] = pi;
	    --l;
/* L90: */
	}
/*  find nrint, the total number of knot intervals and nreg, the number */
/*  of panels in which the approximation domain is subdivided by the */
/*  intersection of knots. */
	ntt = *nt - 7;
	npp = *np - 7;
	nrr = npp / 2;
	nr1 = nrr + 1;
	nrint = ntt + npp;
	nreg = ntt * npp;
/*  arrange the data points according to the panel they belong to. */
	fporde_(&teta[1], &phi[1], m, &c__3, &c__3, &tt[1], nt, &tp[1], np, &
		nummer[1], &index[1], &nreg);
/*  find the b-spline coefficients coco and cosi of the cubic spline */
/*  approximations sc(phi) and ss(phi) for cos(phi) and sin(phi). */
	i__2 = npp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    coco[i__] = 0.f;
	    cosi[i__] = 0.f;
	    i__3 = npp;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * a_dim1] = 0.f;
/* L100: */
	    }
	}
/*  the coefficients coco and cosi are obtained from the conditions */
/*  sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4. */
	i__3 = npp;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    l2 = i__ + 3;
	    arg = tp[l2];
	    fpbspl_(&tp[1], np, &c__3, &arg, &l2, hp);
	    i__2 = npp;
	    for (j = 1; j <= i__2; ++j) {
		row[j] = 0.f;
/* L110: */
	    }
	    ll = i__;
	    for (j = 1; j <= 3; ++j) {
		if (ll > npp) {
		    ll = 1;
		}
		row[ll] += hp[j - 1];
		++ll;
/* L120: */
	    }
	    facc = cos(arg);
	    facs = sin(arg);
	    i__2 = npp;
	    for (j = 1; j <= i__2; ++j) {
		piv = row[j];
		if (piv == 0.f) {
		    goto L140;
		}
		fpgivs_(&piv, &a[j + a_dim1], &co, &si);
		fprota_(&co, &si, &facc, &coco[j]);
		fprota_(&co, &si, &facs, &cosi[j]);
		if (j == npp) {
		    goto L150;
		}
		j1 = j + 1;
		i2 = 1;
		i__4 = npp;
		for (l = j1; l <= i__4; ++l) {
		    ++i2;
		    fprota_(&co, &si, &row[l], &a[j + i2 * a_dim1]);
/* L130: */
		}
L140:
		;
	    }
L150:
	    ;
	}
	fpback_(&a[a_offset], &coco[1], &npp, &npp, &coco[1], ncc);
	fpback_(&a[a_offset], &cosi[1], &npp, &npp, &cosi[1], ncc);
/*  find ncof, the dimension of the spherical spline and ncoff, the */
/*  number of coefficients in the standard b-spline representation. */
	nt4 = *nt - 4;
	np4 = *np - 4;
	ncoff = nt4 * np4;
	ncof = npp * (ntt - 1) + 6;
/*  find the bandwidth of the observation matrix a. */
	iband = npp << 2;
	if (ntt == 4) {
	    iband = (npp + 1) * 3;
	}
	if (ntt < 4) {
	    iband = ncof;
	}
	iband1 = iband - 1;
/*  initialize the observation matrix a. */
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    f[i__] = 0.f;
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
		a[i__ + j * a_dim1] = 0.f;
/* L160: */
	    }
	}
/*  initialize the sum of squared residuals. */
	*fp = 0.f;
/*  fetch the data points in the new order. main loop for the */
/*  different panels. */
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
/*  fix certain constants for the current panel; jrot records the column */
/*  number of the first non-zero element in a row of the observation */
/*  matrix according to a data point of the panel. */
	    num1 = num - 1;
	    lt = num1 / npp;
	    l1 = lt + 4;
	    lp = num1 - lt * npp + 1;
	    l2 = lp + 3;
	    ++lt;
	    jrot = 0;
	    if (lt > 2) {
		jrot = (lt - 3) * npp + 3;
	    }
/*  test whether there are still data points in the current panel. */
	    in = index[num];
L170:
	    if (in == 0) {
		goto L340;
	    }
/*  fetch a new data point. */
	    wi = w[in];
	    ri = r__[in] * wi;
/*  evaluate for the teta-direction, the 4 non-zero b-splines at teta(in) */
	    fpbspl_(&tt[1], nt, &c__3, &teta[in], &l1, ht);
/*  evaluate for the phi-direction, the 4 non-zero b-splines at phi(in) */
	    fpbspl_(&tp[1], np, &c__3, &phi[in], &l2, hp);
/*  store the value of these b-splines in spt and spp resp. */
	    for (i__ = 1; i__ <= 4; ++i__) {
		spp[in + i__ * spp_dim1] = hp[i__ - 1];
		spt[in + i__ * spt_dim1] = ht[i__ - 1];
/* L180: */
	    }
/*  initialize the new row of observation matrix. */
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		h__[i__] = 0.f;
/* L190: */
	    }
/*  calculate the non-zero elements of the new row by making the cross */
/*  products of the non-zero b-splines in teta- and phi-direction and */
/*  by taking into account the conditions of the spherical splines. */
	    i__3 = npp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		row[i__] = 0.f;
/* L200: */
	    }
/*  take into account the condition (3) of the spherical splines. */
	    ll = lp;
	    for (i__ = 1; i__ <= 4; ++i__) {
		if (ll > npp) {
		    ll = 1;
		}
		row[ll] += hp[i__ - 1];
		++ll;
/* L210: */
	    }
/*  take into account the other conditions of the spherical splines. */
	    if (lt > 2 && lt < ntt - 1) {
		goto L230;
	    }
	    facc = 0.f;
	    facs = 0.f;
	    i__3 = npp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		facc += row[i__] * coco[i__];
		facs += row[i__] * cosi[i__];
/* L220: */
	    }
/*  fill in the non-zero elements of the new row. */
L230:
	    j1 = 0;
	    for (j = 1; j <= 4; ++j) {
		jlt = j + lt;
		htj = ht[j - 1];
		if (jlt > 2 && jlt <= nt4) {
		    goto L240;
		}
		++j1;
		h__[j1] += htj;
		goto L280;
L240:
		if (jlt == 3 || jlt == nt4) {
		    goto L260;
		}
		i__3 = npp;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    ++j1;
		    h__[j1] = row[i__] * htj;
/* L250: */
		}
		goto L280;
L260:
		if (jlt == 3) {
		    goto L270;
		}
		h__[j1 + 1] = facc * htj;
		h__[j1 + 2] = facs * htj;
		h__[j1 + 3] = htj;
		j1 += 2;
		goto L280;
L270:
		h__[1] += htj;
		h__[2] = facc * htj;
		h__[3] = facs * htj;
		j1 = 3;
L280:
		;
	    }
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		h__[i__] *= wi;
/* L290: */
	    }
/*  rotate the row into triangle by givens transformations. */
	    irot = jrot;
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		++irot;
		piv = h__[i__];
		if (piv == 0.f) {
		    goto L310;
		}
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the right hand side. */
		fprota_(&co, &si, &ri, &f[irot]);
		if (i__ == iband) {
		    goto L320;
		}
/*  apply that transformation to the left hand side. */
		i2 = 1;
		i3 = i__ + 1;
		i__4 = iband;
		for (j = i3; j <= i__4; ++j) {
		    ++i2;
		    fprota_(&co, &si, &h__[j], &a[irot + i2 * a_dim1]);
/* L300: */
		}
L310:
		;
	    }
/*  add the contribution of the row to the sum of squares of residual */
/*  right hand sides. */
L320:
/* Computing 2nd power */
	    d__1 = ri;
	    *fp += d__1 * d__1;
/*  find the number of the next data point in the panel. */
	    in = nummer[in];
	    goto L170;
L340:
	    ;
	}
/*  find dmax, the maximum value for the diagonal elements in the reduced */
/*  triangle. */
	dmax__ = 0.f;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (a[i__ + a_dim1] <= dmax__) {
		goto L350;
	    }
	    dmax__ = a[i__ + a_dim1];
L350:
	    ;
	}
/*  check whether the observation matrix is rank deficient. */
	sigma = eps * dmax__;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (a[i__ + a_dim1] <= sigma) {
		goto L370;
	    }
/* L360: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&a[a_offset], &f[1], &ncof, &iband, &c__[1], ncc);
	rank = ncof;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + q_dim1] = a[i__ + a_dim1] / dmax__;
/* L365: */
	}
	goto L390;
/*  in case of rank deficiency, find the minimum norm solution. */
L370:
	lwest = ncof * iband + ncof + iband;
	if (*lwrk < lwest) {
	    goto L925;
	}
	lf = 1;
	lh = lf + ncof;
	la = lh + iband;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ff[i__] = f[i__];
	    i__3 = iband;
	    for (j = 1; j <= i__3; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L380: */
	    }
	}
	fprank_(&q[q_offset], &ff[1], &ncof, &iband, ncc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L385: */
	}
/*  add to the sum of squared residuals, the contribution of reducing */
/*  the rank. */
	*fp += sq;
/*  find the coefficients in the standard b-spline representation of */
/*  the spherical spline. */
L390:
	fprpsp_(nt, np, &coco[1], &cosi[1], &c__[1], &ff[1], &ncoff);
/*  test whether the least-squares spline is an acceptable solution. */
	if (*iopt < 0) {
	    if (*fp <= 0.) {
		goto L970;
	    }
	    goto L980;
	}
	fpms = *fp - *s;
	if (abs(fpms) <= acc) {
	    if (*fp <= 0.) {
		goto L970;
	    }
	    goto L980;
	}
/*  if f(p=inf) < s, accept the choice of knots. */
	if (fpms < 0.f) {
	    goto L580;
	}
/*  test whether we cannot further increase the number of knots. */
	if (ncof > *m) {
	    goto L935;
	}
/*  search where to add a new knot. */
/*  find for each interval the sum of squared residuals fpint for the */
/*  data points having the coordinate belonging to that knot interval. */
/*  calculate also coord which is the same sum, weighted by the position */
/*  of the data points considered. */
	i__3 = nrint;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    fpint[i__] = 0.f;
	    coord[i__] = 0.f;
/* L450: */
	}
	i__3 = nreg;
	for (num = 1; num <= i__3; ++num) {
	    num1 = num - 1;
	    lt = num1 / npp;
	    l1 = lt + 1;
	    lp = num1 - lt * npp;
	    l2 = lp + 1 + ntt;
	    jrot = lt * np4 + lp;
	    in = index[num];
L460:
	    if (in == 0) {
		goto L490;
	    }
	    store = 0.f;
	    i1 = jrot;
	    for (i__ = 1; i__ <= 4; ++i__) {
		hti = spt[in + i__ * spt_dim1];
		j1 = i1;
		for (j = 1; j <= 4; ++j) {
		    ++j1;
		    store += hti * spp[in + j * spp_dim1] * c__[j1];
/* L470: */
		}
		i1 += np4;
/* L480: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (r__[in] - store);
	    store = d__1 * d__1;
	    fpint[l1] += store;
	    coord[l1] += store * teta[in];
	    fpint[l2] += store;
	    coord[l2] += store * phi[in];
	    in = nummer[in];
	    goto L460;
L490:
	    ;
	}
/*  find the interval for which fpint is maximal on the condition that */
/*  there still can be added a knot. */
	l1 = 1;
	l2 = nrint;
	if (*ntest < *nt + 1) {
	    l1 = ntt + 1;
	}
	if (*npest < *np + 2) {
	    l2 = ntt;
	}
/*  test whether we cannot further increase the number of knots. */
	if (l1 > l2) {
	    goto L950;
	}
L500:
	fpmax = 0.f;
	l = 0;
	i__3 = l2;
	for (i__ = l1; i__ <= i__3; ++i__) {
	    if (fpmax >= fpint[i__]) {
		goto L510;
	    }
	    l = i__;
	    fpmax = fpint[i__];
L510:
	    ;
	}
	if (l == 0) {
	    goto L930;
	}
/*  calculate the position of the new knot. */
	arg = coord[l] / fpint[l];
/*  test in what direction the new knot is going to be added. */
	if (l > ntt) {
	    goto L530;
	}
/*  addition in the teta-direction */
	l4 = l + 4;
	fpint[l] = 0.f;
	fac1 = tt[l4] - arg;
	fac2 = arg - tt[l4 - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
	j = *nt;
	i__3 = *nt;
	for (i__ = l4; i__ <= i__3; ++i__) {
	    tt[j + 1] = tt[j];
	    --j;
/* L520: */
	}
	tt[l4] = arg;
	++(*nt);
	goto L570;
/*  addition in the phi-direction */
L530:
	l4 = l + 4 - ntt;
	if (arg < pi) {
	    goto L540;
	}
	arg -= pi;
	l4 -= nrr;
L540:
	fpint[l] = 0.f;
	fac1 = tp[l4] - arg;
	fac2 = arg - tp[l4 - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
	ll = nrr + 4;
	j = ll;
	i__3 = ll;
	for (i__ = l4; i__ <= i__3; ++i__) {
	    tp[j + 1] = tp[j];
	    --j;
/* L550: */
	}
	tp[l4] = arg;
	*np += 2;
	++nrr;
	i__3 = ll;
	for (i__ = 5; i__ <= i__3; ++i__) {
	    j = i__ + nrr;
	    tp[j] = tp[i__] + pi;
/* L560: */
	}
/*  restart the computations with the new set of knots. */
L570:
	;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing spherical spline.             c */
/* ********************************************************             c */
/* we have determined the number of knots and their position. we now    c */
/* compute the coefficients of the smoothing spline sp(teta,phi).       c */
/* the observation matrix a is extended by the rows of a matrix, expres-c */
/* sing that sp(teta,phi) must be a constant function in the variable   c */
/* phi and a cubic polynomial in the variable teta. the corresponding   c */
/* weights of these additional rows are set to 1/(p). iteratively       c */
/* we than have to determine the value of p such that f(p) = sum((w(i)* c */
/* (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              c */
/* we already know that the least-squares polynomial corresponds to p=0,c */
/* and that the least-squares spherical spline corresponds to p=infin.  c */
/* the iteration process makes use of rational interpolation. since f(p)c */
/* is a convex and strictly decreasing function of p, it can be approx- c */
/* imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c */
/* three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c */
/* f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c */
/* of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tt(l),l=5,...,nt-4. */
L580:
    fpdisc_(&tt[1], nt, &c__5, &bt[bt_offset], ntest);
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tp(l),l=5,...,np-4. */
    fpdisc_(&tp[1], np, &c__5, &bp[bp_offset], npest);
/*  initial value for p. */
    p1 = 0.f;
    f1 = *sup - *s;
    p3 = -one;
    f3 = fpms;
    p = 0.f;
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p += a[i__ + a_dim1];
/* L585: */
    }
    rn = (doublereal) ncof;
    p = rn / p;
/*  find the bandwidth of the extended observation matrix. */
    iband4 = iband + 3;
    if (ntt <= 4) {
	iband4 = ncof;
    }
    iband3 = iband4 - 1;
    ich1 = 0;
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	pinv = one / p;
/*  store the triangularized observation matrix into q. */
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ff[i__] = f[i__];
	    i__2 = iband4;
	    for (j = 1; j <= i__2; ++j) {
		q[i__ + j * q_dim1] = 0.f;
/* L590: */
	    }
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L600: */
	    }
	}
/*  extend the observation matrix with the rows of a matrix, expressing */
/*  that for teta=cst. sp(teta,phi) must be a constant function. */
	nt6 = *nt - 6;
	i__2 = np4;
	for (i__ = 5; i__ <= i__2; ++i__) {
	    ii = i__ - 4;
	    i__3 = npp;
	    for (l = 1; l <= i__3; ++l) {
		row[l] = 0.f;
/* L610: */
	    }
	    ll = ii;
	    for (l = 1; l <= 5; ++l) {
		if (ll > npp) {
		    ll = 1;
		}
		row[ll] += bp[ii + l * bp_dim1];
		++ll;
/* L620: */
	    }
	    facc = 0.f;
	    facs = 0.f;
	    i__3 = npp;
	    for (l = 1; l <= i__3; ++l) {
		facc += row[l] * coco[l];
		facs += row[l] * cosi[l];
/* L630: */
	    }
	    i__3 = nt6;
	    for (j = 1; j <= i__3; ++j) {
/*  initialize the new row. */
		i__4 = iband;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = 0.f;
/* L640: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		jrot = (j - 2) * npp + 4;
		if (j > 1 && j < nt6) {
		    goto L650;
		}
		h__[1] = facc;
		h__[2] = facs;
		if (j == 1) {
		    jrot = 2;
		}
		goto L670;
L650:
		i__4 = npp;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = row[l];
/* L660: */
		}
L670:
		i__4 = iband;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] *= pinv;
/* L675: */
		}
		ri = 0.f;
/*  rotate the new row into triangle by givens transformations. */
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
		    piv = h__[1];
/* Computing MIN */
		    i__5 = iband1, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
		    if (piv == 0.f) {
			if (i2 <= 0) {
			    goto L720;
			}
			goto L690;
		    }
/*  calculate the parameters of the givens transformation. */
		    fpgivs_(&piv, &q[irot + q_dim1], &co, &si);
/*  apply that givens transformation to the right hand side. */
		    fprota_(&co, &si, &ri, &ff[irot]);
		    if (i2 == 0) {
			goto L720;
		    }
/*  apply that givens transformation to the left hand side. */
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			l1 = l + 1;
			fprota_(&co, &si, &h__[l1], &q[irot + l1 * q_dim1]);
/* L680: */
		    }
L690:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			h__[l] = h__[l + 1];
/* L700: */
		    }
		    h__[i2 + 1] = 0.f;
/* L710: */
		}
L720:
		;
	    }
	}
/*  extend the observation matrix with the rows of a matrix expressing */
/*  that for phi=cst. sp(teta,phi) must be a cubic polynomial. */
	i__3 = nt4;
	for (i__ = 5; i__ <= i__3; ++i__) {
	    ii = i__ - 4;
	    i__2 = npp;
	    for (j = 1; j <= i__2; ++j) {
/*  initialize the new row */
		i__4 = iband4;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] = 0.f;
/* L730: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		j1 = 1;
		for (l = 1; l <= 5; ++l) {
		    il = ii + l;
		    ij = npp;
		    if (il != 3 && il != nt4) {
			goto L750;
		    }
		    j1 = j1 + 3 - j;
		    j2 = j1 - 2;
		    ij = 0;
		    if (il != 3) {
			goto L740;
		    }
		    j1 = 1;
		    j2 = 2;
		    ij = j + 2;
L740:
		    h__[j2] = bt[ii + l * bt_dim1] * coco[j];
		    h__[j2 + 1] = bt[ii + l * bt_dim1] * cosi[j];
L750:
		    h__[j1] += bt[ii + l * bt_dim1];
		    j1 += ij;
/* L760: */
		}
		i__4 = iband4;
		for (l = 1; l <= i__4; ++l) {
		    h__[l] *= pinv;
/* L765: */
		}
		ri = 0.f;
		jrot = 1;
		if (ii > 2) {
		    jrot = j + 3 + (ii - 3) * npp;
		}
/*  rotate the new row into triangle by givens transformations. */
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
		    piv = h__[1];
/* Computing MIN */
		    i__5 = iband3, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
		    if (piv == 0.f) {
			if (i2 <= 0) {
			    goto L810;
			}
			goto L780;
		    }
/*  calculate the parameters of the givens transformation. */
		    fpgivs_(&piv, &q[irot + q_dim1], &co, &si);
/*  apply that givens transformation to the right hand side. */
		    fprota_(&co, &si, &ri, &ff[irot]);
		    if (i2 == 0) {
			goto L810;
		    }
/*  apply that givens transformation to the left hand side. */
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			l1 = l + 1;
			fprota_(&co, &si, &h__[l1], &q[irot + l1 * q_dim1]);
/* L770: */
		    }
L780:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
			h__[l] = h__[l + 1];
/* L790: */
		    }
		    h__[i2 + 1] = 0.f;
/* L800: */
		}
L810:
		;
	    }
	}
/*  find dmax, the maximum value for the diagonal elements in the */
/*  reduced triangle. */
	dmax__ = 0.f;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (q[i__ + q_dim1] <= dmax__) {
		goto L820;
	    }
	    dmax__ = q[i__ + q_dim1];
L820:
	    ;
	}
/*  check whether the matrix is rank deficient. */
	sigma = eps * dmax__;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (q[i__ + q_dim1] <= sigma) {
		goto L840;
	    }
/* L830: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&q[q_offset], &ff[1], &ncof, &iband4, &c__[1], ncc);
	rank = ncof;
	goto L845;
/*  in case of rank deficiency, find the minimum norm solution. */
L840:
	lwest = ncof * iband4 + ncof + iband4;
	if (*lwrk < lwest) {
	    goto L925;
	}
	lf = 1;
	lh = lf + ncof;
	la = lh + iband4;
	fprank_(&q[q_offset], &ff[1], &ncof, &iband4, ncc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
L845:
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L850: */
	}
/*  find the coefficients in the standard b-spline representation of */
/*  the spherical spline. */
	fprpsp_(nt, np, &coco[1], &cosi[1], &c__[1], &ff[1], &ncoff);
/*  compute f(p). */
	*fp = 0.f;
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
	    num1 = num - 1;
	    lt = num1 / npp;
	    lp = num1 - lt * npp;
	    jrot = lt * np4 + lp;
	    in = index[num];
L860:
	    if (in == 0) {
		goto L890;
	    }
	    store = 0.f;
	    i1 = jrot;
	    for (i__ = 1; i__ <= 4; ++i__) {
		hti = spt[in + i__ * spt_dim1];
		j1 = i1;
		for (j = 1; j <= 4; ++j) {
		    ++j1;
		    store += hti * spp[in + j * spp_dim1] * c__[j1];
/* L870: */
		}
		i1 += np4;
/* L880: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (r__[in] - store);
	    *fp += d__1 * d__1;
	    in = nummer[in];
	    goto L860;
L890:
	    ;
	}
/*  test whether the approximation sp(teta,phi) is an acceptable solution */
	fpms = *fp - *s;
	if (abs(fpms) <= acc) {
	    goto L980;
	}
/*  test whether the maximum allowable number of iterations has been */
/*  reached. */
	if (iter == *maxit) {
	    goto L940;
	}
/*  carry out one more step of the iteration process. */
	p2 = p;
	f2 = fpms;
	if (ich3 != 0) {
	    goto L900;
	}
	if (f2 - f3 > acc) {
	    goto L895;
	}
/*  our initial choice of p is too large. */
	p3 = p2;
	f3 = f2;
	p *= con4;
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
	goto L920;
L895:
	if (f2 < 0.f) {
	    ich3 = 1;
	}
L900:
	if (ich1 != 0) {
	    goto L910;
	}
	if (f1 - f2 > acc) {
	    goto L905;
	}
/*  our initial choice of p is too small */
	p1 = p2;
	f1 = f2;
	p /= con4;
	if (p3 < 0.f) {
	    goto L920;
	}
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
	goto L920;
L905:
	if (f2 > 0.f) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
L910:
	if (f2 >= f1 || f2 <= f3) {
	    goto L945;
	}
/*  find the new value of p. */
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
L920:
	;
    }
/*  error codes and messages. */
L925:
    *ier = lwest;
    goto L990;
L930:
    *ier = 5;
    goto L990;
L935:
    *ier = 4;
    goto L990;
L940:
    *ier = 3;
    goto L990;
L945:
    *ier = 2;
    goto L990;
L950:
    *ier = 1;
    goto L990;
L960:
    *ier = -2;
    goto L990;
L970:
    *ier = -1;
    *fp = 0.f;
L980:
    if (ncof != rank) {
	*ier = -rank;
    }
L990:
    return 0;
} /* fpsphe_ */

