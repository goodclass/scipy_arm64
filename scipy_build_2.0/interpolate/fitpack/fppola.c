/* fppola.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fppola_(integer *iopt1, integer *iopt2, integer *iopt3, 
	integer *m, doublereal *u, doublereal *v, doublereal *z__, doublereal 
	*w, D_fp rad, doublereal *s, integer *nuest, integer *nvest, 
	doublereal *eta, doublereal *tol, integer *maxit, integer *ib1, 
	integer *ib3, integer *nc, integer *ncc, integer *intest, integer *
	nrest, integer *nu, doublereal *tu, integer *nv, doublereal *tv, 
	doublereal *c__, doublereal *fp, doublereal *sup, doublereal *fpint, 
	doublereal *coord, doublereal *f, doublereal *ff, doublereal *row, 
	doublereal *cs, doublereal *cosi, doublereal *a, doublereal *q, 
	doublereal *bu, doublereal *bv, doublereal *spu, doublereal *spv, 
	doublereal *h__, integer *index, integer *nummer, doublereal *wrk, 
	integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, bu_dim1, bu_offset, bv_dim1, 
	    bv_offset, spu_dim1, spu_offset, spv_dim1, spv_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal p, r__, c1, c2, c3, c4, f1, f2, f3;
    static integer i1, i2, i3, j1, j2, l1;
    static doublereal p1, p2, p3;
    static integer l2, l3, l4;
    static doublereal u2, u3;
    static integer la;
    static doublereal co;
    static integer ii, lf, il;
    static doublereal pi;
    static integer in;
    static doublereal si;
    static integer lh, ll;
    static doublereal hu[4], wi, rn;
    static integer lu;
    static doublereal sq, zi;
    static integer lv;
    static doublereal hv[4], uu, pi2;
    static integer nr1, nu4, nv4;
    static doublereal acc, fac, arg, one, hui, huj, eps, ten;
    static integer jlu;
    static doublereal piv;
    static integer num, nrr;
    static doublereal fac1, fac2, two;
    static integer nvv, nuu, ich1, ich3;
    static doublereal con1, con4, con9;
    static integer num1;
    static doublereal half;
    static integer ncof;
    static doublereal dmax__;
    static integer ipar, nreg, rank, iter;
    static doublereal fpms, pinv;
    static integer irot, jrot, ipar1, iband, ncoff;
    static doublereal sigma, three, fpmax, ratio;
    static integer numin, nvmin, nrint;
    static doublereal store;
    static integer iband3, iband4, lwest, iband1;
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
	    doublereal *, doublereal *, doublereal *), fprppo_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..user supplied function.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fporde,fpbspl,fpback,fpgivs,fprota,fprank,fpdisc,fprppo */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --nummer;
    spv_dim1 = *m;
    spv_offset = 1 + spv_dim1;
    spv -= spv_offset;
    spu_dim1 = *m;
    spu_offset = 1 + spu_dim1;
    spu -= spu_offset;
    --w;
    --z__;
    --v;
    --u;
    bu_dim1 = *nuest;
    bu_offset = 1 + bu_dim1;
    bu -= bu_offset;
    --tu;
    bv_dim1 = *nvest;
    bv_offset = 1 + bv_dim1;
    bv -= bv_offset;
    cosi -= 6;
    --cs;
    --row;
    --tv;
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
    one = 1.;
    two = 2.;
    three = 3.;
    ten = 10.;
    half = .5f;
    con1 = .1f;
    con9 = .9f;
    con4 = .04f;
    pi = atan(one) * 4;
    pi2 = pi + pi;
    ipar = *iopt2 * (*iopt2 + 3) / 2;
    ipar1 = ipar + 1;
    eps = sqrt(*eta);
    if (*iopt1 < 0) {
	goto L90;
    }
    numin = 9;
    nvmin = *iopt2 * (*iopt2 + 1) + 9;
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
    if (*iopt1 == 0) {
	goto L10;
    }
    if (*s < *sup) {
	if (*nv < nvmin) {
	    goto L70;
	}
	goto L90;
    }
/*  if iopt1 = 0 we begin by computing the weighted least-squares */
/*  polymomial of the form */
/*     s(u,v) = f(1)*(1-u**3)+f(2)*u**3+f(3)*(u**2-u**3)+f(4)*(u-u**3) */
/*  where f(4) = 0 if iopt2> 0 , f(3) = 0 if iopt2 > 1 and */
/*        f(2) = 0 if iopt3> 0. */
/*  the corresponding weighted sum of squared residuals gives the upper */
/*  bound sup for the smoothing factor s. */
L10:
    *sup = 0.f;
    for (i__ = 1; i__ <= 4; ++i__) {
	f[i__] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    a[i__ + j * a_dim1] = 0.f;
/* L20: */
	}
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__];
	zi = z__[i__] * wi;
	uu = u[i__];
	u2 = uu * uu;
	u3 = uu * u2;
	h__[1] = (one - u3) * wi;
	h__[2] = u3 * wi;
	h__[3] = u2 * (one - uu) * wi;
	h__[4] = uu * (one - u2) * wi;
	if (*iopt3 != 0) {
	    h__[2] = 0.f;
	}
	if (*iopt2 > 1) {
	    h__[3] = 0.f;
	}
	if (*iopt2 > 0) {
	    h__[4] = 0.f;
	}
	for (j = 1; j <= 4; ++j) {
	    piv = h__[j];
	    if (piv == 0.f) {
		goto L40;
	    }
	    fpgivs_(&piv, &a[j + a_dim1], &co, &si);
	    fprota_(&co, &si, &zi, &f[j]);
	    if (j == 4) {
		goto L40;
	    }
	    j1 = j + 1;
	    j2 = 1;
	    for (l = j1; l <= 4; ++l) {
		++j2;
		fprota_(&co, &si, &h__[l], &a[j + j2 * a_dim1]);
/* L30: */
	    }
L40:
	    ;
	}
	*sup += zi * zi;
/* L50: */
    }
    if (a[a_dim1 + 4] != 0.f) {
	f[4] /= a[a_dim1 + 4];
    }
    if (a[a_dim1 + 3] != 0.f) {
	f[3] = (f[3] - a[(a_dim1 << 1) + 3] * f[4]) / a[a_dim1 + 3];
    }
    if (a[a_dim1 + 2] != 0.f) {
	f[2] = (f[2] - a[(a_dim1 << 1) + 2] * f[3] - a[a_dim1 * 3 + 2] * f[4])
		 / a[a_dim1 + 2];
    }
    if (a[a_dim1 + 1] != 0.f) {
	f[1] = (f[1] - a[(a_dim1 << 1) + 1] * f[2] - a[a_dim1 * 3 + 1] * f[3] 
		- a[(a_dim1 << 2) + 1] * f[4]) / a[a_dim1 + 1];
    }
/*  find the b-spline representation of this least-squares polynomial */
    c1 = f[1];
    c4 = f[2];
    c2 = f[4] / three + c1;
    c3 = (f[3] + two * f[4]) / three + c1;
    *nu = 8;
    *nv = 8;
    for (i__ = 1; i__ <= 4; ++i__) {
	c__[i__] = c1;
	c__[i__ + 4] = c2;
	c__[i__ + 8] = c3;
	c__[i__ + 12] = c4;
	tu[i__] = 0.f;
	tu[i__ + 4] = one;
	rn = (doublereal) ((i__ << 1) - 9);
	tv[i__] = rn * pi;
	rn = (doublereal) ((i__ << 1) - 1);
	tv[i__ + 4] = rn * pi;
/* L60: */
    }
    *fp = *sup;
/*  test whether the least-squares polynomial is an acceptable solution */
    fpms = *sup - *s;
    if (fpms < acc) {
	goto L960;
    }
/*  test whether we cannot further increase the number of knots. */
L70:
    if (*nuest < numin || *nvest < nvmin) {
	goto L950;
    }
/*  find the initial set of interior knots of the spline in case iopt1=0. */
    *nu = numin;
    *nv = nvmin;
    tu[5] = half;
    nvv = *nv - 8;
    rn = (doublereal) (nvv + 1);
    fac = pi2 / rn;
    i__1 = nvv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rn = (doublereal) i__;
	tv[i__ + 4] = rn * fac - pi;
/* L80: */
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 1 : computation of least-squares bicubic splines.              c */
/*  ******************************************************              c */
/*  if iopt1<0 we compute the least-squares bicubic spline according    c */
/*  to the given set of knots.                                          c */
/*  if iopt1>=0 we compute least-squares bicubic splines with in-       c */
/*  creasing numbers of knots until the corresponding sum f(p=inf)<=s.  c */
/*  the initial set of knots then depends on the value of iopt1         c */
/*    if iopt1=0 we start with one interior knot in the u-direction     c */
/*              (0.5) and 1+iopt2*(iopt2+1) in the v-direction.         c */
/*    if iopt1>0 we start with the set of knots found at the last       c */
/*              call of the routine.                                    c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
L90:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(u,v). */
	l1 = 4;
	l2 = l1;
	l3 = *nv - 3;
	l4 = l3;
	tv[l2] = -pi;
	tv[l3] = pi;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ++l1;
	    --l2;
	    ++l3;
	    --l4;
	    tv[l2] = tv[l4] - pi2;
	    tv[l3] = tv[l1] + pi2;
/* L120: */
	}
	l = *nu;
	for (i__ = 1; i__ <= 4; ++i__) {
	    tu[i__] = 0.f;
	    tu[l] = one;
	    --l;
/* L130: */
	}
/*  find nrint, the total number of knot intervals and nreg, the number */
/*  of panels in which the approximation domain is subdivided by the */
/*  intersection of knots. */
	nuu = *nu - 7;
	nvv = *nv - 7;
	nrr = nvv / 2;
	nr1 = nrr + 1;
	nrint = nuu + nvv;
	nreg = nuu * nvv;
/*  arrange the data points according to the panel they belong to. */
	fporde_(&u[1], &v[1], m, &c__3, &c__3, &tu[1], nu, &tv[1], nv, &
		nummer[1], &index[1], &nreg);
	if (*iopt2 == 0) {
	    goto L195;
	}
/*  find the b-spline coefficients cosi of the cubic spline */
/*  approximations for cr(v)=rad(v)*cos(v) and sr(v) = rad(v)*sin(v) */
/*  if iopt2=1, and additionally also for cr(v)**2,sr(v)**2 and */
/*  2*cr(v)*sr(v) if iopt2=2 */
	i__2 = nvv;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = ipar;
	    for (j = 1; j <= i__3; ++j) {
		cosi[j + i__ * 5] = 0.f;
/* L135: */
	    }
	    i__3 = nvv;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * a_dim1] = 0.f;
/* L140: */
	    }
	}
/*  the coefficients cosi are obtained from interpolation conditions */
/*  at the knots tv(i),i=4,5,...nv-4. */
	i__3 = nvv;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    l2 = i__ + 3;
	    arg = tv[l2];
	    fpbspl_(&tv[1], nv, &c__3, &arg, &l2, hv);
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
		row[j] = 0.f;
/* L145: */
	    }
	    ll = i__;
	    for (j = 1; j <= 3; ++j) {
		if (ll > nvv) {
		    ll = 1;
		}
		row[ll] += hv[j - 1];
		++ll;
/* L150: */
	    }
	    co = cos(arg);
	    si = sin(arg);
	    r__ = (*rad)(&arg);
	    cs[1] = co * r__;
	    cs[2] = si * r__;
	    if (*iopt2 == 1) {
		goto L155;
	    }
	    cs[3] = cs[1] * cs[1];
	    cs[4] = cs[2] * cs[2];
	    cs[5] = cs[1] * cs[2];
L155:
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
		piv = row[j];
		if (piv == 0.f) {
		    goto L170;
		}
		fpgivs_(&piv, &a[j + a_dim1], &co, &si);
		i__4 = ipar;
		for (l = 1; l <= i__4; ++l) {
		    fprota_(&co, &si, &cs[l], &cosi[l + j * 5]);
/* L160: */
		}
		if (j == nvv) {
		    goto L175;
		}
		j1 = j + 1;
		j2 = 1;
		i__4 = nvv;
		for (l = j1; l <= i__4; ++l) {
		    ++j2;
		    fprota_(&co, &si, &row[l], &a[j + j2 * a_dim1]);
/* L165: */
		}
L170:
		;
	    }
L175:
	    ;
	}
	i__3 = ipar;
	for (l = 1; l <= i__3; ++l) {
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
		cs[j] = cosi[l + j * 5];
/* L180: */
	    }
	    fpback_(&a[a_offset], &cs[1], &nvv, &nvv, &cs[1], ncc);
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
		cosi[l + j * 5] = cs[j];
/* L185: */
	    }
/* L190: */
	}
/*  find ncof, the dimension of the spline and ncoff, the number */
/*  of coefficients in the standard b-spline representation. */
L195:
	nu4 = *nu - 4;
	nv4 = *nv - 4;
	ncoff = nu4 * nv4;
	ncof = ipar1 + nvv * (nu4 - 1 - *iopt2 - *iopt3);
/*  find the bandwidth of the observation matrix a. */
	iband = nvv << 2;
	if (nuu - *iopt2 - *iopt3 <= 1) {
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
/* L200: */
	    }
	}
/*  initialize the sum of squared residuals. */
	*fp = 0.f;
	ratio = one + tu[6] / tu[5];
/*  fetch the data points in the new order. main loop for the */
/*  different panels. */
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
/*  fix certain constants for the current panel; jrot records the column */
/*  number of the first non-zero element in a row of the observation */
/*  matrix according to a data point of the panel. */
	    num1 = num - 1;
	    lu = num1 / nvv;
	    l1 = lu + 4;
	    lv = num1 - lu * nvv + 1;
	    l2 = lv + 3;
	    jrot = 0;
	    if (lu > *iopt2) {
		jrot = ipar1 + (lu - *iopt2 - 1) * nvv;
	    }
	    ++lu;
/*  test whether there are still data points in the current panel. */
	    in = index[num];
L210:
	    if (in == 0) {
		goto L380;
	    }
/*  fetch a new data point. */
	    wi = w[in];
	    zi = z__[in] * wi;
/*  evaluate for the u-direction, the 4 non-zero b-splines at u(in) */
	    fpbspl_(&tu[1], nu, &c__3, &u[in], &l1, hu);
/*  evaluate for the v-direction, the 4 non-zero b-splines at v(in) */
	    fpbspl_(&tv[1], nv, &c__3, &v[in], &l2, hv);
/*  store the value of these b-splines in spu and spv resp. */
	    for (i__ = 1; i__ <= 4; ++i__) {
		spu[in + i__ * spu_dim1] = hu[i__ - 1];
		spv[in + i__ * spv_dim1] = hv[i__ - 1];
/* L220: */
	    }
/*  initialize the new row of observation matrix. */
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		h__[i__] = 0.f;
/* L240: */
	    }
/*  calculate the non-zero elements of the new row by making the cross */
/*  products of the non-zero b-splines in u- and v-direction and */
/*  by taking into account the conditions of the splines. */
	    i__3 = nvv;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		row[i__] = 0.f;
/* L250: */
	    }
/*  take into account the periodicity condition of the bicubic splines. */
	    ll = lv;
	    for (i__ = 1; i__ <= 4; ++i__) {
		if (ll > nvv) {
		    ll = 1;
		}
		row[ll] += hv[i__ - 1];
		++ll;
/* L260: */
	    }
/*  take into account the other conditions of the splines. */
	    if (*iopt2 == 0 || lu > *iopt2 + 1) {
		goto L280;
	    }
	    i__3 = ipar;
	    for (l = 1; l <= i__3; ++l) {
		cs[l] = 0.f;
		i__4 = nvv;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    cs[l] += row[i__] * cosi[l + i__ * 5];
/* L270: */
		}
	    }
/*  fill in the non-zero elements of the new row. */
L280:
	    j1 = 0;
	    for (j = 1; j <= 4; ++j) {
		jlu = j + lu;
		huj = hu[j - 1];
		if (jlu > *iopt2 + 2) {
		    goto L320;
		}
		switch (jlu) {
		    case 1:  goto L290;
		    case 2:  goto L290;
		    case 3:  goto L300;
		    case 4:  goto L310;
		}
L290:
		h__[1] = huj;
		j1 = 1;
		goto L330;
L300:
		h__[1] += huj;
		h__[2] = huj * cs[1];
		h__[3] = huj * cs[2];
		j1 = 3;
		goto L330;
L310:
		h__[1] += huj;
		h__[2] += huj * ratio * cs[1];
		h__[3] += huj * ratio * cs[2];
		h__[4] = huj * cs[3];
		h__[5] = huj * cs[4];
		h__[6] = huj * cs[5];
		j1 = 6;
		goto L330;
L320:
		if (jlu > nu4 && *iopt3 != 0) {
		    goto L330;
		}
		i__4 = nvv;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ++j1;
		    h__[j1] = row[i__] * huj;
/* L325: */
		}
L330:
		;
	    }
	    i__4 = iband;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		h__[i__] *= wi;
/* L335: */
	    }
/*  rotate the row into triangle by givens transformations. */
	    irot = jrot;
	    i__4 = iband;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		++irot;
		piv = h__[i__];
		if (piv == 0.f) {
		    goto L350;
		}
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the right hand side. */
		fprota_(&co, &si, &zi, &f[irot]);
		if (i__ == iband) {
		    goto L360;
		}
/*  apply that transformation to the left hand side. */
		i2 = 1;
		i3 = i__ + 1;
		i__3 = iband;
		for (j = i3; j <= i__3; ++j) {
		    ++i2;
		    fprota_(&co, &si, &h__[j], &a[irot + i2 * a_dim1]);
/* L340: */
		}
L350:
		;
	    }
/*  add the contribution of the row to the sum of squares of residual */
/*  right hand sides. */
L360:
/* Computing 2nd power */
	    d__1 = zi;
	    *fp += d__1 * d__1;
/*  find the number of the next data point in the panel. */
	    in = nummer[in];
	    goto L210;
L380:
	    ;
	}
/*  find dmax, the maximum value for the diagonal elements in the reduced */
/*  triangle. */
	dmax__ = 0.f;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (a[i__ + a_dim1] <= dmax__) {
		goto L390;
	    }
	    dmax__ = a[i__ + a_dim1];
L390:
	    ;
	}
/*  check whether the observation matrix is rank deficient. */
	sigma = eps * dmax__;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (a[i__ + a_dim1] <= sigma) {
		goto L410;
	    }
/* L400: */
	}
/*  backward substitution in case of full rank. */
	fpback_(&a[a_offset], &f[1], &ncof, &iband, &c__[1], ncc);
	rank = ncof;
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + q_dim1] = a[i__ + a_dim1] / dmax__;
/* L405: */
	}
	goto L430;
/*  in case of rank deficiency, find the minimum norm solution. */
L410:
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
	    i__4 = iband;
	    for (j = 1; j <= i__4; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L420: */
	    }
	}
	fprank_(&q[q_offset], &ff[1], &ncof, &iband, ncc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
	i__4 = ncof;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    q[i__ + q_dim1] /= dmax__;
/* L425: */
	}
/*  add to the sum of squared residuals, the contribution of reducing */
/*  the rank. */
	*fp += sq;
/*  find the coefficients in the standard b-spline representation of */
/*  the spline. */
L430:
	fprppo_(nu, nv, iopt2, iopt3, &cosi[6], &ratio, &c__[1], &ff[1], &
		ncoff);
/*  test whether the least-squares spline is an acceptable solution. */
	if (*iopt1 < 0) {
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
/*  test whether we cannot further increase the number of knots */
	if (*m < ncof) {
	    goto L935;
	}
/*  search where to add a new knot. */
/*  find for each interval the sum of squared residuals fpint for the */
/*  data points having the coordinate belonging to that knot interval. */
/*  calculate also coord which is the same sum, weighted by the position */
/*  of the data points considered. */
	i__4 = nrint;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    fpint[i__] = 0.f;
	    coord[i__] = 0.f;
/* L450: */
	}
	i__4 = nreg;
	for (num = 1; num <= i__4; ++num) {
	    num1 = num - 1;
	    lu = num1 / nvv;
	    l1 = lu + 1;
	    lv = num1 - lu * nvv;
	    l2 = lv + 1 + nuu;
	    jrot = lu * nv4 + lv;
	    in = index[num];
L460:
	    if (in == 0) {
		goto L490;
	    }
	    store = 0.f;
	    i1 = jrot;
	    for (i__ = 1; i__ <= 4; ++i__) {
		hui = spu[in + i__ * spu_dim1];
		j1 = i1;
		for (j = 1; j <= 4; ++j) {
		    ++j1;
		    store += hui * spv[in + j * spv_dim1] * c__[j1];
/* L470: */
		}
		i1 += nv4;
/* L480: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    store = d__1 * d__1;
	    fpint[l1] += store;
	    coord[l1] += store * u[in];
	    fpint[l2] += store;
	    coord[l2] += store * v[in];
	    in = nummer[in];
	    goto L460;
L490:
	    ;
	}
/* bring together the information concerning knot panels which are */
/* symmetric with respect to the origin. */
	i__4 = nrr;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    l1 = nuu + i__;
	    l2 = l1 + nrr;
	    fpint[l1] += fpint[l2];
	    coord[l1] = coord[l1] + coord[l2] - pi * fpint[l2];
/* L495: */
	}
/*  find the interval for which fpint is maximal on the condition that */
/*  there still can be added a knot. */
	l1 = 1;
	l2 = nuu + nrr;
	if (*nuest < *nu + 1) {
	    l1 = nuu + 1;
	}
	if (*nvest < *nv + 2) {
	    l2 = nuu;
	}
/*  test whether we cannot further increase the number of knots. */
	if (l1 > l2) {
	    goto L950;
	}
L500:
	fpmax = 0.f;
	l = 0;
	i__4 = l2;
	for (i__ = l1; i__ <= i__4; ++i__) {
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
	if (l > nuu) {
	    goto L530;
	}
/*  addition in the u-direction */
	l4 = l + 4;
	fpint[l] = 0.f;
	fac1 = tu[l4] - arg;
	fac2 = arg - tu[l4 - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
	j = *nu;
	i__4 = *nu;
	for (i__ = l4; i__ <= i__4; ++i__) {
	    tu[j + 1] = tu[j];
	    --j;
/* L520: */
	}
	tu[l4] = arg;
	++(*nu);
	goto L570;
/*  addition in the v-direction */
L530:
	l4 = l + 4 - nuu;
	fpint[l] = 0.f;
	fac1 = tv[l4] - arg;
	fac2 = arg - tv[l4 - 1];
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
	ll = nrr + 4;
	j = ll;
	i__4 = ll;
	for (i__ = l4; i__ <= i__4; ++i__) {
	    tv[j + 1] = tv[j];
	    --j;
/* L550: */
	}
	tv[l4] = arg;
	*nv += 2;
	++nrr;
	i__4 = ll;
	for (i__ = 5; i__ <= i__4; ++i__) {
	    j = i__ + nrr;
	    tv[j] = tv[i__] + pi;
/* L560: */
	}
/*  restart the computations with the new set of knots. */
L570:
	;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing bicubic spline.               c */
/* ******************************************************               c */
/* we have determined the number of knots and their position. we now    c */
/* compute the coefficients of the smoothing spline sp(u,v).            c */
/* the observation matrix a is extended by the rows of a matrix, expres-c */
/* sing that sp(u,v) must be a constant function in the variable        c */
/* v and a cubic polynomial in the variable u. the corresponding        c */
/* weights of these additional rows are set to 1/(p). iteratively       c */
/* we than have to determine the value of p such that f(p) = sum((w(i)* c */
/* (z(i)-sp(u(i),v(i))))**2)  be = s.                                   c */
/* we already know that the least-squares polynomial corresponds to p=0,c */
/* and that the least-squares bicubic spline corresponds to p=infin.    c */
/* the iteration process makes use of rational interpolation. since f(p)c */
/* is a convex and strictly decreasing function of p, it can be approx- c */
/* imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c */
/* three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c */
/* f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c */
/* of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tu(l),l=5,...,nu-4. */
L580:
    fpdisc_(&tu[1], nu, &c__5, &bu[bu_offset], nuest);
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tv(l),l=5,...,nv-4. */
    fpdisc_(&tv[1], nv, &c__5, &bv[bv_offset], nvest);
/*  initial value for p. */
    p1 = 0.f;
    f1 = *sup - *s;
    p3 = -one;
    f3 = fpms;
    p = 0.f;
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p += a[i__ + a_dim1];
/* L590: */
    }
    rn = (doublereal) ncof;
    p = rn / p;
/*  find the bandwidth of the extended observation matrix. */
    iband4 = iband + ipar1;
    if (iband4 > ncof) {
	iband4 = ncof;
    }
    iband3 = iband4 - 1;
    ich1 = 0;
    ich3 = 0;
    nuu = nu4 - *iopt3 - 1;
/*  iteration process to find the root of f(p)=s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	pinv = one / p;
/*  store the triangularized observation matrix into q. */
	i__4 = ncof;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    ff[i__] = f[i__];
	    i__2 = iband4;
	    for (j = 1; j <= i__2; ++j) {
		q[i__ + j * q_dim1] = 0.f;
/* L620: */
	    }
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/* L630: */
	    }
	}
/*  extend the observation matrix with the rows of a matrix, expressing */
/*  that for u=constant sp(u,v) must be a constant function. */
	i__2 = nv4;
	for (i__ = 5; i__ <= i__2; ++i__) {
	    ii = i__ - 4;
	    i__4 = nvv;
	    for (l = 1; l <= i__4; ++l) {
		row[l] = 0.f;
/* L635: */
	    }
	    ll = ii;
	    for (l = 1; l <= 5; ++l) {
		if (ll > nvv) {
		    ll = 1;
		}
		row[ll] += bv[ii + l * bv_dim1];
		++ll;
/* L640: */
	    }
	    i__4 = nuu;
	    for (j = 1; j <= i__4; ++j) {
/*  initialize the new row. */
		i__3 = iband;
		for (l = 1; l <= i__3; ++l) {
		    h__[l] = 0.f;
/* L645: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		if (j > *iopt2) {
		    goto L665;
		}
		if (j == 2) {
		    goto L655;
		}
		for (k = 1; k <= 2; ++k) {
		    cs[k] = 0.f;
		    i__3 = nvv;
		    for (l = 1; l <= i__3; ++l) {
			cs[k] += cosi[k + l * 5] * row[l];
/* L650: */
		    }
		}
		h__[1] = cs[1];
		h__[2] = cs[2];
		jrot = 2;
		goto L675;
L655:
		for (k = 3; k <= 5; ++k) {
		    cs[k] = 0.f;
		    i__3 = nvv;
		    for (l = 1; l <= i__3; ++l) {
			cs[k] += cosi[k + l * 5] * row[l];
/* L660: */
		    }
		}
		h__[1] = cs[1] * ratio;
		h__[2] = cs[2] * ratio;
		h__[3] = cs[3];
		h__[4] = cs[4];
		h__[5] = cs[5];
		jrot = 2;
		goto L675;
L665:
		i__3 = nvv;
		for (l = 1; l <= i__3; ++l) {
		    h__[l] = row[l];
/* L670: */
		}
		jrot = ipar1 + 1 + (j - *iopt2 - 1) * nvv;
L675:
		i__3 = iband;
		for (l = 1; l <= i__3; ++l) {
		    h__[l] *= pinv;
/* L677: */
		}
		zi = 0.f;
/*  rotate the new row into triangle by givens transformations. */
		i__3 = ncof;
		for (irot = jrot; irot <= i__3; ++irot) {
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
		    fprota_(&co, &si, &zi, &ff[irot]);
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
/*  that for v=constant. sp(u,v) must be a cubic polynomial. */
	i__4 = nu4;
	for (i__ = 5; i__ <= i__4; ++i__) {
	    ii = i__ - 4;
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
/*  initialize the new row */
		i__3 = iband4;
		for (l = 1; l <= i__3; ++l) {
		    h__[l] = 0.f;
/* L730: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
		j1 = 1;
		for (l = 1; l <= 5; ++l) {
		    il = ii + l - 1;
		    if (il == nu4 && *iopt3 != 0) {
			goto L760;
		    }
		    if (il > *iopt2 + 1) {
			goto L750;
		    }
		    switch (il) {
			case 1:  goto L735;
			case 2:  goto L740;
			case 3:  goto L745;
		    }
L735:
		    h__[1] = bu[ii + l * bu_dim1];
		    j1 = j + 1;
		    goto L760;
L740:
		    h__[1] += bu[ii + l * bu_dim1];
		    h__[2] = bu[ii + l * bu_dim1] * cosi[j * 5 + 1];
		    h__[3] = bu[ii + l * bu_dim1] * cosi[j * 5 + 2];
		    j1 = j + 3;
		    goto L760;
L745:
		    h__[1] += bu[ii + l * bu_dim1];
		    h__[2] = bu[ii + l * bu_dim1] * cosi[j * 5 + 1] * ratio;
		    h__[3] = bu[ii + l * bu_dim1] * cosi[j * 5 + 2] * ratio;
		    h__[4] = bu[ii + l * bu_dim1] * cosi[j * 5 + 3];
		    h__[5] = bu[ii + l * bu_dim1] * cosi[j * 5 + 4];
		    h__[6] = bu[ii + l * bu_dim1] * cosi[j * 5 + 5];
		    j1 = j + 6;
		    goto L760;
L750:
		    h__[j1] = bu[ii + l * bu_dim1];
		    j1 += nvv;
L760:
		    ;
		}
		i__3 = iband4;
		for (l = 1; l <= i__3; ++l) {
		    h__[l] *= pinv;
/* L765: */
		}
		zi = 0.f;
		jrot = 1;
		if (ii > *iopt2 + 1) {
		    jrot = ipar1 + (ii - *iopt2 - 2) * nvv + j;
		}
/*  rotate the new row into triangle by givens transformations. */
		i__3 = ncof;
		for (irot = jrot; irot <= i__3; ++irot) {
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
		    fprota_(&co, &si, &zi, &ff[irot]);
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
/*  the polar spline. */
	fprppo_(nu, nv, iopt2, iopt3, &cosi[6], &ratio, &c__[1], &ff[1], &
		ncoff);
/*  compute f(p). */
	*fp = 0.f;
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
	    num1 = num - 1;
	    lu = num1 / nvv;
	    lv = num1 - lu * nvv;
	    jrot = lu * nv4 + lv;
	    in = index[num];
L860:
	    if (in == 0) {
		goto L890;
	    }
	    store = 0.f;
	    i1 = jrot;
	    for (i__ = 1; i__ <= 4; ++i__) {
		hui = spu[in + i__ * spu_dim1];
		j1 = i1;
		for (j = 1; j <= 4; ++j) {
		    ++j1;
		    store += hui * spv[in + j * spv_dim1] * c__[j1];
/* L870: */
		}
		i1 += nv4;
/* L880: */
	    }
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    *fp += d__1 * d__1;
	    in = nummer[in];
	    goto L860;
L890:
	    ;
	}
/*  test whether the approximation sp(u,v) is an acceptable solution */
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
} /* fppola_ */

