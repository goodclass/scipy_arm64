/* specfun.f -- translated by f2c (version 20190311).
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

static doublereal c_b4 = 2.;
static doublecomplex c_b6 = {1.,0.};
static doublecomplex c_b15 = {2.,0.};
static doublereal c_b50 = .5;
static doublereal c_b51 = 1.5;
static integer c_n1 = -1;
static doublereal c_b149 = 3.;
static doublereal c_b221 = .33333333333333331;
static integer c__200 = 200;
static integer c__15 = 15;
static real c_b317 = 2.f;
static doublereal c_b323 = 1.27;
static doublecomplex c_b340 = {.5,0.};
static integer c__0 = 0;
static doublereal c_b524 = 1.;
static doublecomplex c_b772 = {-.125,0.};
static doublecomplex c_b776 = {.375,0.};
static integer c__1 = 1;
static doublecomplex c_b1162 = {-2.,0.};
static doublereal c_b1169 = .33333;
static doublecomplex c_b1387 = {-1.,0.};

/*       COMPUTATION OF SPECIAL FUNCTIONS */

/*          Shanjie Zhang and Jianming Jin */

/*       Copyrighted but permission granted to use code in programs. */
/*       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc. */

/*       Scipy changes: */
/*       - Compiled into a single source file and changed REAL To DBLE throughout. */
/*       - Changed according to ERRATA. */
/*       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts. */
/*       - Made functions return sf_error codes in ISFER variables instead */
/*         of printing warnings. The codes are */
/*         - SF_ERROR_OK        = 0: no error */
/*         - SF_ERROR_SINGULAR  = 1: singularity encountered */
/*         - SF_ERROR_UNDERFLOW = 2: floating point underflow */
/*         - SF_ERROR_OVERFLOW  = 3: floating point overflow */
/*         - SF_ERROR_SLOW      = 4: too many iterations required */
/*         - SF_ERROR_LOSS      = 5: loss of precision */
/*         - SF_ERROR_NO_RESULT = 6: no result obtained */
/*         - SF_ERROR_DOMAIN    = 7: out of domain */
/*         - SF_ERROR_ARG       = 8: invalid input parameter */
/*         - SF_ERROR_OTHER     = 9: unclassified error */

doublereal dnan_(void)
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = 0.;
    ret_val = 0. / ret_val;
    return ret_val;
} /* dnan_ */

doublereal dinf_(void)
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = 1e300;
    ret_val *= ret_val;
    return ret_val;
} /* dinf_ */

/* Subroutine */ int cpdsa_(integer *n, doublecomplex *z__, doublecomplex *
	cdn)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double sqrt(doublereal);
    void z_exp(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer m;
    static doublereal g0, g1, pd;
    static doublecomplex cr;
    static doublereal gm, pi, vm, xn, vt;
    static doublecomplex ca0, cb0;
    static doublereal ga0, va0, sq2;
    static doublecomplex cdw;
    static doublereal eps;
    extern /* Subroutine */ int gaih_(doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute complex parabolic cylinder function Dn(z) */
/*                for small argument */
/*       Input:   z   --- complex argument of D(z) */
/*                n   --- Order of D(z) (n = 0,-1,-2,...) */
/*       Output:  CDN --- Dn(z) */
/*       Routine called: GAIH for computing Г(x), x=n/2 (n=1,2,...) */
/*       =========================================================== */

    eps = 1e-15;
    pi = 3.141592653589793;
    sq2 = sqrt(2.);
    z__3.r = z__->r * -.25, z__3.i = z__->i * -.25;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z_exp(&z__1, &z__2);
    ca0.r = z__1.r, ca0.i = z__1.i;
    va0 = (1. - *n) * .5;
    if ((real) (*n) == 0.f) {
	cdn->r = ca0.r, cdn->i = ca0.i;
    } else {
	if (z_abs(z__) == 0.f) {
	    if (va0 <= 0.f && va0 == (doublereal) ((integer) va0)) {
		cdn->r = 0., cdn->i = 0.;
	    } else {
		gaih_(&va0, &ga0);
		d__1 = *n * -.5;
		pd = sqrt(pi) / (pow_dd(&c_b4, &d__1) * ga0);
		z__1.r = pd, z__1.i = 0.;
		cdn->r = z__1.r, cdn->i = z__1.i;
	    }
	} else {
	    xn = (doublereal) (-(*n));
	    gaih_(&xn, &g1);
	    d__2 = *n * -.5 - 1.;
	    d__1 = pow_dd(&c_b4, &d__2);
	    z__2.r = d__1 * ca0.r, z__2.i = d__1 * ca0.i;
	    z__1.r = z__2.r / g1, z__1.i = z__2.i / g1;
	    cb0.r = z__1.r, cb0.i = z__1.i;
	    vt = *n * -.5;
	    gaih_(&vt, &g0);
	    z__1.r = g0, z__1.i = 0.;
	    cdn->r = z__1.r, cdn->i = z__1.i;
	    cr.r = 1., cr.i = 0.;
	    for (m = 1; m <= 250; ++m) {
		vm = (m - *n) * .5;
		gaih_(&vm, &gm);
		z__4.r = -cr.r, z__4.i = -cr.i;
		z__3.r = sq2 * z__4.r, z__3.i = sq2 * z__4.i;
		z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * 
			z__->i + z__3.i * z__->r;
		d__1 = (doublereal) m;
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr.r = z__1.r, cr.i = z__1.i;
		z__1.r = gm * cr.r, z__1.i = gm * cr.i;
		cdw.r = z__1.r, cdw.i = z__1.i;
		z__1.r = cdn->r + cdw.r, z__1.i = cdn->i + cdw.i;
		cdn->r = z__1.r, cdn->i = z__1.i;
		if (z_abs(&cdw) < z_abs(cdn) * eps) {
		    goto L20;
		}
/* L10: */
	    }
L20:
	    z__1.r = cb0.r * cdn->r - cb0.i * cdn->i, z__1.i = cb0.r * cdn->i 
		    + cb0.i * cdn->r;
	    cdn->r = z__1.r, cdn->i = z__1.i;
	}
    }
    return 0;
} /* cpdsa_ */

/*       ********************************** */
/* Subroutine */ int cfs_(doublecomplex *z__, doublecomplex *zf, 
	doublecomplex *zd)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_sqrt(
	    doublecomplex *, doublecomplex *), z_sin(doublecomplex *, 
	    doublecomplex *), z_cos(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, m;
    static doublecomplex s;
    static doublereal w0;
    static doublecomplex z0, cf, cg, cr;
    static doublereal pi, wb;
    static doublecomplex zp, cf0, cf1;
    static doublereal wb0;
    static doublecomplex zp2;
    static doublereal eps;


/*       ========================================================= */
/*       Purpose: Compute complex Fresnel Integral S(z) and S'(z) */
/*       Input :  z  --- Argument of S(z) */
/*       Output:  ZF --- S(z) */
/*                ZD --- S'(z) */
/*       ========================================================= */

    eps = 1e-14;
    pi = 3.141592653589793;
    w0 = z_abs(z__);
    d__1 = pi * .5;
    z__2.r = d__1 * z__->r, z__2.i = d__1 * z__->i;
    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * z__->i + 
	    z__2.i * z__->r;
    zp.r = z__1.r, zp.i = z__1.i;
    z__1.r = zp.r * zp.r - zp.i * zp.i, z__1.i = zp.r * zp.i + zp.i * zp.r;
    zp2.r = z__1.r, zp2.i = z__1.i;
    z0.r = 0., z0.i = 0.;
    if (z__->r == z0.r && z__->i == z0.i) {
	s.r = z0.r, s.i = z0.i;
    } else if (w0 <= 2.5f) {
	z__2.r = z__->r * zp.r - z__->i * zp.i, z__2.i = z__->r * zp.i + 
		z__->i * zp.r;
	z__1.r = z__2.r / 3., z__1.i = z__2.i / 3.;
	s.r = z__1.r, s.i = z__1.i;
	cr.r = s.r, cr.i = s.i;
	wb0 = 0.;
	for (k = 1; k <= 80; ++k) {
	    z__6.r = cr.r * -.5, z__6.i = cr.i * -.5;
	    d__1 = k * 4. - 1.;
	    z__5.r = d__1 * z__6.r, z__5.i = d__1 * z__6.i;
	    d__2 = (doublereal) k;
	    z__4.r = z__5.r / d__2, z__4.i = z__5.i / d__2;
	    d__3 = k * 2. + 1.;
	    z__3.r = z__4.r / d__3, z__3.i = z__4.i / d__3;
	    d__4 = k * 4. + 3.;
	    z__2.r = z__3.r / d__4, z__2.i = z__3.i / d__4;
	    z__1.r = z__2.r * zp2.r - z__2.i * zp2.i, z__1.i = z__2.r * zp2.i 
		    + z__2.i * zp2.r;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = s.r + cr.r, z__1.i = s.i + cr.i;
	    s.r = z__1.r, s.i = z__1.i;
	    wb = z_abs(&s);
	    if ((d__1 = wb - wb0, abs(d__1)) < eps && k > 10) {
		goto L30;
	    }
/* L10: */
	    wb0 = wb;
	}
    } else if (w0 > 2.5f && w0 < 4.5f) {
	m = 85;
	s.r = z0.r, s.i = z0.i;
	cf1.r = z0.r, cf1.i = z0.i;
	cf0.r = 1e-100, cf0.i = 0.;
	for (k = m; k >= 0; --k) {
	    d__1 = k * 2. + 3.;
	    z__3.r = d__1 * cf0.r, z__3.i = d__1 * cf0.i;
	    z_div(&z__2, &z__3, &zp);
	    z__1.r = z__2.r - cf1.r, z__1.i = z__2.i - cf1.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	    if (k != k / 2 << 1) {
		z__1.r = s.r + cf.r, z__1.i = s.i + cf.i;
		s.r = z__1.r, s.i = z__1.i;
	    }
	    cf1.r = cf0.r, cf1.i = cf0.i;
/* L15: */
	    cf0.r = cf.r, cf0.i = cf.i;
	}
	z__6.r = pi * zp.r, z__6.i = pi * zp.i;
	z_div(&z__5, &c_b15, &z__6);
	z_sqrt(&z__4, &z__5);
	z_sin(&z__7, &zp);
	z__3.r = z__4.r * z__7.r - z__4.i * z__7.i, z__3.i = z__4.r * z__7.i 
		+ z__4.i * z__7.r;
	z_div(&z__2, &z__3, &cf);
	z__1.r = z__2.r * s.r - z__2.i * s.i, z__1.i = z__2.r * s.i + z__2.i *
		 s.r;
	s.r = z__1.r, s.i = z__1.i;
    } else {
	cr.r = 1., cr.i = 0.;
	cf.r = 1., cf.i = 0.;
	for (k = 1; k <= 20; ++k) {
	    z__4.r = cr.r * -.25, z__4.i = cr.i * -.25;
	    d__1 = k * 4. - 1.;
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    d__2 = k * 4. - 3.;
	    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
	    z_div(&z__1, &z__2, &zp2);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L20: */
	    z__1.r = cf.r + cr.r, z__1.i = cf.i + cr.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	}
	cr.r = 1., cr.i = 0.;
	cg.r = cr.r, cg.i = cr.i;
	for (k = 1; k <= 12; ++k) {
	    z__4.r = cr.r * -.25, z__4.i = cr.i * -.25;
	    d__1 = k * 4. + 1.;
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    d__2 = k * 4. - 1.;
	    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
	    z_div(&z__1, &z__2, &zp2);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L25: */
	    z__1.r = cg.r + cr.r, z__1.i = cg.i + cr.i;
	    cg.r = z__1.r, cg.i = z__1.i;
	}
	z__3.r = pi * z__->r, z__3.i = pi * z__->i;
	z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i 
		+ z__3.i * z__->r;
	z_div(&z__1, &cg, &z__2);
	cg.r = z__1.r, cg.i = z__1.i;
	z_cos(&z__5, &zp);
	z__4.r = cf.r * z__5.r - cf.i * z__5.i, z__4.i = cf.r * z__5.i + cf.i 
		* z__5.r;
	z_sin(&z__7, &zp);
	z__6.r = cg.r * z__7.r - cg.i * z__7.i, z__6.i = cg.r * z__7.i + cg.i 
		* z__7.r;
	z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
	z__8.r = pi * z__->r, z__8.i = pi * z__->i;
	z_div(&z__2, &z__3, &z__8);
	z__1.r = .5 - z__2.r, z__1.i = -z__2.i;
	s.r = z__1.r, s.i = z__1.i;
    }
L30:
    zf->r = s.r, zf->i = s.i;
    d__1 = pi * .5f;
    z__3.r = d__1 * z__->r, z__3.i = d__1 * z__->i;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z_sin(&z__1, &z__2);
    zd->r = z__1.r, zd->i = z__1.i;
    return 0;
} /* cfs_ */

/*       ********************************** */
/* Subroutine */ int lqmn_(integer *mm, integer *m, integer *n, doublereal *x,
	 doublereal *qm, doublereal *qd)
{
    /* System generated locals */
    integer qm_dim1, qm_offset, qd_dim1, qd_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal q0, q1, q10, qf;
    static integer km, ls;
    static doublereal xq, xs, qf0, qf1, qf2;


/*       ========================================================== */
/*       Purpose: Compute the associated Legendre functions of the */
/*                second kind, Qmn(x) and Qmn'(x) */
/*       Input :  x  --- Argument of Qmn(x) */
/*                m  --- Order of Qmn(x)  ( m = 0,1,2,… ) */
/*                n  --- Degree of Qmn(x) ( n = 0,1,2,… ) */
/*                mm --- Physical dimension of QM and QD */
/*       Output:  QM(m,n) --- Qmn(x) */
/*                QD(m,n) --- Qmn'(x) */
/*       ========================================================== */

    /* Parameter adjustments */
    qd_dim1 = *mm - 0 + 1;
    qd_offset = 0 + qd_dim1 * 0;
    qd -= qd_offset;
    qm_dim1 = *mm - 0 + 1;
    qm_offset = 0 + qm_dim1 * 0;
    qm -= qm_offset;

    /* Function Body */
    if (abs(*x) == 1.) {
	i__1 = *m;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = 0; j <= i__2; ++j) {
		qm[i__ + j * qm_dim1] = 1e300;
		qd[i__ + j * qd_dim1] = 1e300;
/* L10: */
	    }
	}
	return 0;
    }
    ls = 1;
    if (abs(*x) > 1.) {
	ls = -1;
    }
    xs = ls * (1. - *x * *x);
    xq = sqrt(xs);
    q0 = log((d__1 = (*x + 1.) / (*x - 1.), abs(d__1))) * .5;
    if (abs(*x) < 1.0001) {
	qm[0] = q0;
	qm[qm_dim1] = *x * q0 - 1.;
	qm[1] = -1. / xq;
	qm[qm_dim1 + 1] = -ls * xq * (q0 + *x / (1. - *x * *x));
	for (i__ = 0; i__ <= 1; ++i__) {
	    i__2 = *n;
	    for (j = 2; j <= i__2; ++j) {
		qm[i__ + j * qm_dim1] = ((j * 2. - 1.) * *x * qm[i__ + (j - 1)
			 * qm_dim1] - (j + i__ - 1.) * qm[i__ + (j - 2) * 
			qm_dim1]) / (j - i__);
/* L15: */
	    }
	}
	i__2 = *n;
	for (j = 0; j <= i__2; ++j) {
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		qm[i__ + j * qm_dim1] = (i__ - 1.) * -2. * *x / xq * qm[i__ - 
			1 + j * qm_dim1] - ls * (j + i__ - 1.) * (j - i__ + 
			2.) * qm[i__ - 2 + j * qm_dim1];
/* L20: */
	    }
	}
    } else {
	if (abs(*x) > 1.1) {
	    km = *m + 40 + *n;
	} else {
	    km = (*m + 40 + *n) * (integer) (-1.f - log(*x - 1.f) * 1.8f);
	}
	qf2 = 0.;
	qf1 = 1.;
	qf0 = 0.;
	for (k = km; k >= 0; --k) {
	    qf0 = (((k << 1) + 3.) * *x * qf1 - (k + 2.) * qf2) / (k + 1.);
	    if (k <= *n) {
		qm[k * qm_dim1] = qf0;
	    }
	    qf2 = qf1;
/* L25: */
	    qf1 = qf0;
	}
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
/* L30: */
	    qm[k * qm_dim1] = q0 * qm[k * qm_dim1] / qf0;
	}
	qf2 = 0.;
	qf1 = 1.;
	for (k = km; k >= 0; --k) {
	    qf0 = (((k << 1) + 3.) * *x * qf1 - (k + 1.) * qf2) / (k + 2.);
	    if (k <= *n) {
		qm[k * qm_dim1 + 1] = qf0;
	    }
	    qf2 = qf1;
/* L35: */
	    qf1 = qf0;
	}
	q10 = -1. / xq;
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
/* L40: */
	    qm[k * qm_dim1 + 1] = q10 * qm[k * qm_dim1 + 1] / qf0;
	}
	i__1 = *n;
	for (j = 0; j <= i__1; ++j) {
	    q0 = qm[j * qm_dim1];
	    q1 = qm[j * qm_dim1 + 1];
	    i__2 = *m - 2;
	    for (i__ = 0; i__ <= i__2; ++i__) {
		qf = (i__ + 1) * -2. * *x / xq * q1 + (j - i__) * (j + i__ + 
			1.) * q0;
		qm[i__ + 2 + j * qm_dim1] = qf;
		q0 = q1;
		q1 = qf;
/* L45: */
	    }
	}
    }
    qd[0] = ls / xs;
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/* L50: */
	qd[j * qd_dim1] = ls * j * (qm[(j - 1) * qm_dim1] - *x * qm[j * 
		qm_dim1]) / xs;
    }
    i__2 = *n;
    for (j = 0; j <= i__2; ++j) {
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    qd[i__ + j * qd_dim1] = ls * i__ * *x / xs * qm[i__ + j * qm_dim1]
		     + (i__ + j) * (j - i__ + 1.) / xq * qm[i__ - 1 + j * 
		    qm_dim1];
/* L55: */
	}
    }
    return 0;
} /* lqmn_ */

/*       ********************************** */
/* Subroutine */ int clpmn_(integer *mm, integer *m, integer *n, doublereal *
	x, doublereal *y, integer *ntype, doublecomplex *cpm, doublecomplex *
	cpd)
{
    /* System generated locals */
    integer cpm_dim1, cpm_offset, cpd_dim1, cpd_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    void z_sqrt(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex z__;
    static integer ls;
    static doublecomplex zq, zs;
    extern doublereal dinf_(void);


/*       ========================================================= */
/*       Purpose: Compute the associated Legendre functions Pmn(z) */
/*                and their derivatives Pmn'(z) for a complex */
/*                argument */
/*       Input :  x     --- Real part of z */
/*                y     --- Imaginary part of z */
/*                m     --- Order of Pmn(z),  m = 0,1,2,...,n */
/*                n     --- Degree of Pmn(z), n = 0,1,2,...,N */
/*                mm    --- Physical dimension of CPM and CPD */
/*                ntype --- type of cut, either 2 or 3 */
/*       Output:  CPM(m,n) --- Pmn(z) */
/*                CPD(m,n) --- Pmn'(z) */
/*       ========================================================= */

    /* Parameter adjustments */
    cpd_dim1 = *mm - 0 + 1;
    cpd_offset = 0 + cpd_dim1 * 0;
    cpd -= cpd_offset;
    cpm_dim1 = *mm - 0 + 1;
    cpm_offset = 0 + cpm_dim1 * 0;
    cpm -= cpm_offset;

    /* Function Body */
    z__1.r = *x, z__1.i = *y;
    z__.r = z__1.r, z__.i = z__1.i;
    i__1 = *n;
    for (i__ = 0; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 0; j <= i__2; ++j) {
	    i__3 = j + i__ * cpm_dim1;
	    cpm[i__3].r = 0., cpm[i__3].i = 0.;
/* L10: */
	    i__3 = j + i__ * cpd_dim1;
	    cpd[i__3].r = 0., cpd[i__3].i = 0.;
	}
    }
    cpm[0].r = 1., cpm[0].i = 0.;
    if (*n == 0) {
	return 0;
    }
    if (abs(*x) == 1. && *y == 0.) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__2 = i__ * cpm_dim1;
	    d__1 = pow_di(x, &i__);
	    cpm[i__2].r = d__1, cpm[i__2].i = 0.;
/* L15: */
	    i__2 = i__ * cpd_dim1;
	    i__1 = i__ + 1;
	    d__1 = i__ * .5 * (i__ + 1) * pow_di(x, &i__1);
	    cpd[i__2].r = d__1, cpd[i__2].i = 0.;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (i__ == 1) {
		    i__3 = i__ + j * cpd_dim1;
		    d__1 = dinf_();
		    cpd[i__3].r = d__1, cpd[i__3].i = 0.;
		} else if (i__ == 2) {
		    i__3 = i__ + j * cpd_dim1;
		    i__4 = j + 1;
		    d__1 = (j + 2) * -.25 * (j + 1) * j * (j - 1) * pow_di(x, 
			    &i__4);
		    cpd[i__3].r = d__1, cpd[i__3].i = 0.;
		}
/* L20: */
	    }
	}
	return 0;
    }
    if (*ntype == 2) {
/*       sqrt(1 - z^2) with branch cut on |x|>1 */
	z__2.r = z__.r * z__.r - z__.i * z__.i, z__2.i = z__.r * z__.i + 
		z__.i * z__.r;
	z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
	zs.r = z__1.r, zs.i = z__1.i;
	z_sqrt(&z__2, &zs);
	z__1.r = -z__2.r, z__1.i = -z__2.i;
	zq.r = z__1.r, zq.i = z__1.i;
	ls = -1;
    } else {
/*       sqrt(z^2 - 1) with branch cut between [-1, 1] */
	z__2.r = z__.r * z__.r - z__.i * z__.i, z__2.i = z__.r * z__.i + 
		z__.i * z__.r;
	z__1.r = z__2.r - 1., z__1.i = z__2.i;
	zs.r = z__1.r, zs.i = z__1.i;
	z_sqrt(&z__1, &zs);
	zq.r = z__1.r, zq.i = z__1.i;
	if (*x < 0.) {
	    z__1.r = -zq.r, z__1.i = -zq.i;
	    zq.r = z__1.r, zq.i = z__1.i;
	}
	ls = 1;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       DLMF 14.7.15 */
/* L25: */
	i__2 = i__ + i__ * cpm_dim1;
	d__1 = i__ * 2. - 1.;
	z__2.r = d__1 * zq.r, z__2.i = d__1 * zq.i;
	i__3 = i__ - 1 + (i__ - 1) * cpm_dim1;
	z__1.r = z__2.r * cpm[i__3].r - z__2.i * cpm[i__3].i, z__1.i = z__2.r 
		* cpm[i__3].i + z__2.i * cpm[i__3].r;
	cpm[i__2].r = z__1.r, cpm[i__2].i = z__1.i;
    }
/* Computing MIN */
    i__3 = *m, i__1 = *n - 1;
    i__2 = min(i__3,i__1);
    for (i__ = 0; i__ <= i__2; ++i__) {
/*       DLMF 14.10.7 */
/* L30: */
	i__3 = i__ + (i__ + 1) * cpm_dim1;
	d__1 = i__ * 2. + 1.;
	z__2.r = d__1 * z__.r, z__2.i = d__1 * z__.i;
	i__1 = i__ + i__ * cpm_dim1;
	z__1.r = z__2.r * cpm[i__1].r - z__2.i * cpm[i__1].i, z__1.i = z__2.r 
		* cpm[i__1].i + z__2.i * cpm[i__1].r;
	cpm[i__3].r = z__1.r, cpm[i__3].i = z__1.i;
    }
    i__3 = *m;
    for (i__ = 0; i__ <= i__3; ++i__) {
	i__1 = *n;
	for (j = i__ + 2; j <= i__1; ++j) {
/*       DLMF 14.10.3 */
	    i__2 = i__ + j * cpm_dim1;
	    d__1 = j * 2. - 1.;
	    z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	    i__4 = i__ + (j - 1) * cpm_dim1;
	    z__3.r = z__4.r * cpm[i__4].r - z__4.i * cpm[i__4].i, z__3.i = 
		    z__4.r * cpm[i__4].i + z__4.i * cpm[i__4].r;
	    d__2 = i__ + j - 1.;
	    i__5 = i__ + (j - 2) * cpm_dim1;
	    z__5.r = d__2 * cpm[i__5].r, z__5.i = d__2 * cpm[i__5].i;
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	    i__6 = j - i__;
	    d__3 = (doublereal) i__6;
	    z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
	    cpm[i__2].r = z__1.r, cpm[i__2].i = z__1.i;
/* L35: */
	}
    }
    cpd[0].r = 0., cpd[0].i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*       DLMF 14.10.5 */
/* L40: */
	i__3 = j * cpd_dim1;
	i__2 = ls * j;
	i__4 = j * cpm_dim1;
	z__4.r = z__.r * cpm[i__4].r - z__.i * cpm[i__4].i, z__4.i = z__.r * 
		cpm[i__4].i + z__.i * cpm[i__4].r;
	i__5 = (j - 1) * cpm_dim1;
	z__3.r = z__4.r - cpm[i__5].r, z__3.i = z__4.i - cpm[i__5].i;
	d__1 = (doublereal) i__2;
	z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	z_div(&z__1, &z__2, &zs);
	cpd[i__3].r = z__1.r, cpd[i__3].i = z__1.i;
    }
    i__3 = *m;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
/*       derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3 */
/*       derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2 */
	    i__4 = i__ + j * cpd_dim1;
	    i__5 = -i__;
	    d__1 = (doublereal) i__5;
	    z__5.r = d__1 * z__.r, z__5.i = d__1 * z__.i;
	    i__1 = i__ + j * cpm_dim1;
	    z__4.r = z__5.r * cpm[i__1].r - z__5.i * cpm[i__1].i, z__4.i = 
		    z__5.r * cpm[i__1].i + z__5.i * cpm[i__1].r;
	    z_div(&z__3, &z__4, &zs);
	    d__2 = (j + i__) * (j - i__ + 1.);
	    z__8.r = d__2, z__8.i = 0.;
	    z_div(&z__7, &z__8, &zq);
	    i__6 = i__ - 1 + j * cpm_dim1;
	    z__6.r = z__7.r * cpm[i__6].r - z__7.i * cpm[i__6].i, z__6.i = 
		    z__7.r * cpm[i__6].i + z__7.i * cpm[i__6].r;
	    z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
	    d__3 = (doublereal) ls;
	    z__1.r = d__3 * z__2.r, z__1.i = d__3 * z__2.i;
	    cpd[i__4].r = z__1.r, cpd[i__4].i = z__1.i;
/* L45: */
	}
    }
    return 0;
} /* clpmn_ */

/*       ********************************** */
/* Subroutine */ int vvsa_(doublereal *va, doublereal *x, doublereal *pv)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sin(doublereal), pow_dd(doublereal *, doublereal *
	    ), sqrt(doublereal);

    /* Local variables */
    static integer m;
    static doublereal r__, a0, g1, r1, v1, gm, ep, pi, gw, vm, sv, ga0, va0, 
	    vb0, sq2, sv0, fac, eps;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute parabolic cylinder function Vv(x) */
/*                for small argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PV --- Vv(x) */
/*       Routine called : GAMMA2 for computing Г(x) */
/*       =================================================== */

    eps = 1e-15;
    pi = 3.141592653589793;
    ep = exp(*x * -.25 * *x);
    va0 = *va * .5 + 1.;
    if (*x == 0.f) {
	if (va0 <= 0.f && va0 == (doublereal) ((integer) va0) || *va == 0.f) {
	    *pv = 0.;
	} else {
	    vb0 = *va * -.5;
	    sv0 = sin(va0 * pi);
	    gamma2_(&va0, &ga0);
	    *pv = pow_dd(&c_b4, &vb0) * sv0 / ga0;
	}
    } else {
	sq2 = sqrt(2.);
	d__1 = *va * -.5;
	a0 = pow_dd(&c_b4, &d__1) * ep / (pi * 2.);
	sv = sin(-(*va + .5) * pi);
	v1 = *va * -.5;
	gamma2_(&v1, &g1);
	*pv = (sv + 1.) * g1;
	r__ = 1.;
	fac = 1.;
	for (m = 1; m <= 250; ++m) {
	    vm = (m - *va) * .5;
	    gamma2_(&vm, &gm);
	    r__ = r__ * sq2 * *x / m;
	    fac = -fac;
	    gw = fac * sv + 1.;
	    r1 = gw * r__ * gm;
	    *pv += r1;
	    if ((d__1 = r1 / *pv, abs(d__1)) < eps && gw != 0.f) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*pv = a0 * *pv;
    }
    return 0;
} /* vvsa_ */

/*       ********************************** */
/*       SciPy: Changed P from a character array to an integer array. */
/* Subroutine */ int jdzo_(integer *nt, integer *n, integer *m, integer *p, 
	doublereal *zo)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal x;
    static integer l0, l1, m1[70], n1[70], l2, p1[70];
    static doublereal x0, x1, x2, bj[101], dj[101], fj[101];
    static integer mm, nm;
    static doublereal xm, zoc[71];
    extern /* Subroutine */ int bjndd_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute the zeros of Bessel functions Jn(x) and */
/*                Jn'(x), and arrange them in the order of their */
/*                magnitudes */
/*       Input :  NT    --- Number of total zeros ( NT ≤ 1200 ) */
/*       Output:  ZO(L) --- Value of the L-th zero of Jn(x) */
/*                          and Jn'(x) */
/*                N(L)  --- n, order of Jn(x) or Jn'(x) associated */
/*                          with the L-th zero */
/*                M(L)  --- m, serial number of the zeros of Jn(x) */
/*                          or Jn'(x) associated with the L-th zero */
/*                          ( L is the serial number of all the */
/*                            zeros of Jn(x) and Jn'(x) ) */
/*                P(L)  --- 0 (TM) or 1 (TE), a code for designating the */
/*                          zeros of Jn(x)  or Jn'(x). */
/*                          In the waveguide applications, the zeros */
/*                          of Jn(x) correspond to TM modes and */
/*                          those of Jn'(x) correspond to TE modes */
/*       Routine called:    BJNDD for computing Jn(x), Jn'(x) and */
/*                          Jn''(x) */
/*       ============================================================= */

    /* Parameter adjustments */
    --p;
    --m;
    --n;

    /* Function Body */
    x = 0.;
    zoc[0] = 0.;
    if (*nt < 600) {
	d__1 = (doublereal) (*nt);
	d__2 = (doublereal) (*nt);
	xm = pow_dd(&d__1, &c_b50) * 2.248485f - 1.f - *nt * .0159382f + 
		pow_dd(&d__2, &c_b51) * 3.208775e-4f;
	nm = (integer) (*nt * .05875f + 14.5f);
	mm = (integer) (*nt * .02f) + 6;
    } else {
	d__1 = (doublereal) (*nt);
	d__2 = (doublereal) (*nt);
	xm = pow_dd(&d__1, &c_b50) * 1.445389f + 5.f + *nt * .01889876f - 
		pow_dd(&d__2, &c_b51) * 2.147763e-4f;
	nm = (integer) (*nt * .0327f + 27.8f);
	mm = (integer) (*nt * .01088f) + 10;
    }
    l0 = 0;
    i__1 = nm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = (doublereal) (i__ - 1);
	x1 = pow_dd(&d__1, &c_b50) * .4795504f + .407658f + (i__ - 1) * 
		.983618f;
	d__1 = (doublereal) (i__ - 1);
	x2 = pow_dd(&d__1, &c_b50) * .8333883f + 1.99535f + (i__ - 1) * 
		.984584f;
	l1 = 0;
	i__2 = mm;
	for (j = 1; j <= i__2; ++j) {
	    if (i__ == 1 && j == 1) {
		goto L15;
	    }
	    x = x1;
L10:
	    bjndd_(&i__, &x, bj, dj, fj);
	    x0 = x;
	    x -= dj[i__ - 1] / fj[i__ - 1];
	    if (x1 > xm) {
		goto L20;
	    }
	    if ((d__1 = x - x0, abs(d__1)) > 1e-10) {
		goto L10;
	    }
L15:
	    ++l1;
	    n1[l1 - 1] = i__ - 1;
	    m1[l1 - 1] = j;
	    if (i__ == 1) {
		m1[l1 - 1] = j - 1;
	    }
	    p1[l1 - 1] = 1;
	    zoc[l1] = x;
	    if (i__ <= 15) {
/* Computing 2nd power */
		i__3 = j + 1;
		x1 = x + 3.057f + (i__ - 1) * .0122f + ((i__ - 1) * .41575f + 
			1.555f) / (i__3 * i__3);
	    } else {
/* Computing 2nd power */
		i__3 = j + 1;
		x1 = x + 2.918f + (i__ - 1) * .01924f + ((i__ - 1) * .13205f 
			+ 6.26f) / (i__3 * i__3);
	    }
L20:
	    x = x2;
L25:
	    bjndd_(&i__, &x, bj, dj, fj);
	    x0 = x;
	    x -= bj[i__ - 1] / dj[i__ - 1];
	    if (x > xm) {
		goto L30;
	    }
	    if ((d__1 = x - x0, abs(d__1)) > 1e-10) {
		goto L25;
	    }
	    ++l1;
	    n1[l1 - 1] = i__ - 1;
	    m1[l1 - 1] = j;
	    p1[l1 - 1] = 0;
	    zoc[l1] = x;
	    if (i__ <= 15) {
/* Computing 2nd power */
		i__3 = j + 1;
		x2 = x + 3.11f + (i__ - 1) * .0138f + ((i__ - 1) * .2804f + 
			.04832f) / (i__3 * i__3);
	    } else {
/* Computing 2nd power */
		i__3 = j + 3;
		x2 = x + 3.001f + (i__ - 1) * .0105f + ((i__ - 1) * .48525f + 
			11.52f) / (i__3 * i__3);
	    }
L30:
	    ;
	}
	l = l0 + l1;
	l2 = l;
L35:
	if (l0 == 0) {
	    i__2 = l;
	    for (k = 1; k <= i__2; ++k) {
		zo[k] = zoc[k];
		n[k] = n1[k - 1];
		m[k] = m1[k - 1];
/* L40: */
		p[k] = p1[k - 1];
	    }
	    l1 = 0;
	} else if (l0 != 0) {
	    if (zo[l0] >= zoc[l1]) {
		zo[l0 + l1] = zo[l0];
		n[l0 + l1] = n[l0];
		m[l0 + l1] = m[l0];
		p[l0 + l1] = p[l0];
		--l0;
	    } else {
		zo[l0 + l1] = zoc[l1];
		n[l0 + l1] = n1[l1 - 1];
		m[l0 + l1] = m1[l1 - 1];
		p[l0 + l1] = p1[l1 - 1];
		--l1;
	    }
	}
	if (l1 != 0) {
	    goto L35;
	}
/* L45: */
	l0 = l2;
    }
    return 0;
} /* jdzo_ */

/*       ********************************** */
/* Subroutine */ int cbk_(integer *m, integer *n, doublereal *c__, doublereal 
	*cv, doublereal *qt, doublereal *ck, doublereal *bk)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t, u[200], v[200], w[200];
    static integer i1, n2;
    static doublereal r1, s1;
    static integer ip, nm;
    static doublereal sw, eps;


/*       ===================================================== */
/*       Purpose: Compute coefficient Bk's for oblate radial */
/*                functions with a small argument */
/*       ===================================================== */

    /* Parameter adjustments */
    --bk;
    --ck;

    /* Function Body */
    eps = 1e-14;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    nm = (integer) ((*n - *m) * .5f + *c__) + 25;
    u[0] = 0.;
    n2 = nm - 2;
    i__1 = n2;
    for (j = 2; j <= i__1; ++j) {
/* L10: */
	u[j - 1] = *c__ * *c__;
    }
    i__1 = n2;
    for (j = 1; j <= i__1; ++j) {
/* L15: */
	v[j - 1] = (j * 2.f - 1.f - ip) * ((j - *m) * 2.f - ip) + *m * (*m - 
		1.f) - *cv;
    }
    i__1 = nm - 1;
    for (j = 1; j <= i__1; ++j) {
/* L20: */
	w[j - 1] = (j * 2.f - ip) * (j * 2.f + 1.f - ip);
    }
    if (ip == 0) {
	sw = 0.;
	i__1 = n2 - 1;
	for (k = 0; k <= i__1; ++k) {
	    s1 = 0.;
	    i1 = k - *m + 1;
	    i__2 = nm;
	    for (i__ = i1; i__ <= i__2; ++i__) {
		if (i__ < 0) {
		    goto L30;
		}
		r1 = 1.;
		i__3 = k;
		for (j = 1; j <= i__3; ++j) {
/* L25: */
		    r1 = r1 * (i__ + *m - j) / j;
		}
		s1 += ck[i__ + 1] * (i__ * 2.f + *m) * r1;
		if ((d__1 = s1 - sw, abs(d__1)) < abs(s1) * eps) {
		    goto L35;
		}
		sw = s1;
L30:
		;
	    }
L35:
	    bk[k + 1] = *qt * s1;
/* L40: */
	}
    } else if (ip == 1) {
	sw = 0.;
	i__1 = n2 - 1;
	for (k = 0; k <= i__1; ++k) {
	    s1 = 0.;
	    i1 = k - *m + 1;
	    i__2 = nm;
	    for (i__ = i1; i__ <= i__2; ++i__) {
		if (i__ < 0) {
		    goto L50;
		}
		r1 = 1.;
		i__3 = k;
		for (j = 1; j <= i__3; ++j) {
/* L45: */
		    r1 = r1 * (i__ + *m - j) / j;
		}
		if (i__ > 0) {
		    s1 += ck[i__] * (i__ * 2.f + *m - 1) * r1;
		}
		s1 -= ck[i__ + 1] * (i__ * 2.f + *m) * r1;
		if ((d__1 = s1 - sw, abs(d__1)) < abs(s1) * eps) {
		    goto L55;
		}
		sw = s1;
L50:
		;
	    }
L55:
	    bk[k + 1] = *qt * s1;
/* L60: */
	}
    }
    w[0] /= v[0];
    bk[1] /= v[0];
    i__1 = n2;
    for (k = 2; k <= i__1; ++k) {
	t = v[k - 1] - w[k - 2] * u[k - 1];
	w[k - 1] /= t;
/* L65: */
	bk[k] = (bk[k] - bk[k - 1] * u[k - 1]) / t;
    }
    for (k = n2 - 1; k >= 1; --k) {
/* L70: */
	bk[k] -= w[k - 1] * bk[k + 1];
    }
    return 0;
} /* cbk_ */

/*       ********************************** */
/* Subroutine */ int rmn2sp_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, doublereal *cv, doublereal *df, integer *kd, 
	doublereal *r2f, doublereal *r2d)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer j, k, j1, j2, l1;
    static doublereal r1, r2, r3, r4, ga, gb, gc, dn[200], pd[252], qd[252];
    static integer ki;
    static doublereal sd;
    static integer ip, nm;
    static doublereal sf, pm[252], qm[252], sw, ck1, ck2, sd0, sd1, sd2;
    static integer nm1, nm2, nm3;
    static doublereal su0, su1, su2, sdm;
    extern /* Subroutine */ int kmn_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *);
    static doublereal eps, spl, sum, spd1, spd2;
    extern /* Subroutine */ int lpmns_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), lqmns_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute prolate spheroidal radial function */
/*                of the second kind with a small argument */
/*       Routines called: */
/*            (1) LPMNS for computing the associated Legendre */
/*                functions of the first kind */
/*            (2) LQMNS for computing the associated Legendre */
/*                functions of the second kind */
/*            (3) KMN for computing expansion coefficients */
/*                and joining factors */
/*       ====================================================== */

    /* Parameter adjustments */
    --df;

    /* Function Body */
    if (abs(df[1]) < 1e-280) {
	*r2f = 1e300;
	*r2d = 1e300;
	return 0;
    }
    eps = 1e-14;
    ip = 1;
    nm1 = (*n - *m) / 2;
    if (*n - *m == nm1 << 1) {
	ip = 0;
    }
    nm = nm1 + 25 + (integer) (*c__);
    nm2 = (nm << 1) + *m;
    kmn_(m, n, c__, cv, kd, &df[1], dn, &ck1, &ck2);
    lpmns_(m, &nm2, x, pm, pd);
    lqmns_(m, &nm2, x, qm, qd);
    su0 = 0.;
    sw = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	j = (k << 1) - 2 + *m + ip;
	su0 += df[k] * qm[j];
	if (k > nm1 && (d__1 = su0 - sw, abs(d__1)) < abs(su0) * eps) {
	    goto L15;
	}
/* L10: */
	sw = su0;
    }
L15:
    sd0 = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	j = (k << 1) - 2 + *m + ip;
	sd0 += df[k] * qd[j];
	if (k > nm1 && (d__1 = sd0 - sw, abs(d__1)) < abs(sd0) * eps) {
	    goto L25;
	}
/* L20: */
	sw = sd0;
    }
L25:
    su1 = 0.;
    sd1 = 0.;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	j = *m - (k << 1) + ip;
	if (j < 0) {
	    j = -j - 1;
	}
	su1 += dn[k - 1] * qm[j];
/* L30: */
	sd1 += dn[k - 1] * qd[j];
    }
    d__1 = (*x - 1.) / (*x + 1.);
    d__2 = *m * .5;
    ga = pow_dd(&d__1, &d__2);
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	j = *m - (k << 1) + ip;
	if (j >= 0) {
	    goto L55;
	}
	if (j < 0) {
	    j = -j - 1;
	}
	r1 = 1.;
	i__2 = j;
	for (j1 = 1; j1 <= i__2; ++j1) {
/* L35: */
	    r1 = (*m + j1) * r1;
	}
	r2 = 1.;
	i__2 = *m - j - 2;
	for (j2 = 1; j2 <= i__2; ++j2) {
/* L40: */
	    r2 = j2 * r2;
	}
	r3 = 1.;
	sf = 1.;
	i__2 = j;
	for (l1 = 1; l1 <= i__2; ++l1) {
	    r3 = r3 * .5 * (-j + l1 - 1.f) * (j + l1) / ((*m + l1) * l1) * (
		    1.f - *x);
/* L45: */
	    sf += r3;
	}
	if (*m - j >= 2) {
	    gb = (*m - j - 1.) * r2;
	}
	if (*m - j <= 1) {
	    gb = 1.;
	}
	spl = r1 * ga * gb * sf;
	i__2 = j + *m;
	su1 += pow_ii(&c_n1, &i__2) * dn[k - 1] * spl;
	spd1 = *m / (*x * *x - 1.) * spl;
	gc = j * .5 * (j + 1.f) / (*m + 1.f);
	sd = 1.;
	r4 = 1.;
	i__2 = j - 1;
	for (l1 = 1; l1 <= i__2; ++l1) {
	    r4 = r4 * .5 * (-j + l1) * (j + l1 + 1.f) / ((*m + l1 + 1.f) * l1)
		     * (1.f - *x);
/* L50: */
	    sd += r4;
	}
	spd2 = r1 * ga * gb * gc * sd;
	i__2 = j + *m;
	sd1 += pow_ii(&c_n1, &i__2) * dn[k - 1] * (spd1 + spd2);
L55:
	;
    }
    su2 = 0.;
    ki = ((*m << 1) + 1 + ip) / 2;
    nm3 = nm + ki;
    i__1 = nm3;
    for (k = ki; k <= i__1; ++k) {
	j = (k << 1) - 1 - *m - ip;
	su2 += dn[k - 1] * pm[j];
	if (j > *m && (d__1 = su2 - sw, abs(d__1)) < abs(su2) * eps) {
	    goto L65;
	}
/* L60: */
	sw = su2;
    }
L65:
    sd2 = 0.;
    i__1 = nm3;
    for (k = ki; k <= i__1; ++k) {
	j = (k << 1) - 1 - *m - ip;
	sd2 += dn[k - 1] * pd[j];
	if (j > *m && (d__1 = sd2 - sw, abs(d__1)) < abs(sd2) * eps) {
	    goto L75;
	}
/* L70: */
	sw = sd2;
    }
L75:
    sum = su0 + su1 + su2;
    sdm = sd0 + sd1 + sd2;
    *r2f = sum / ck2;
    *r2d = sdm / ck2;
    return 0;
} /* rmn2sp_ */

/*       ********************************** */
/* Subroutine */ int bernob_(integer *n, doublereal *bn)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer k, m;
    static doublereal s, r1, r2, tpi;


/*       ====================================== */
/*       Purpose: Compute Bernoulli number Bn */
/*       Input :  n --- Serial number */
/*       Output:  BN(n) --- Bn */
/*       ====================================== */

    tpi = 6.283185307179586;
    bn[0] = 1.;
    bn[1] = -.5;
    bn[2] = .16666666666666666;
/* Computing 2nd power */
    d__1 = 2. / tpi;
    r1 = d__1 * d__1;
    i__1 = *n;
    for (m = 4; m <= i__1; m += 2) {
	r1 = -r1 * (m - 1) * m / (tpi * tpi);
	r2 = 1.;
	for (k = 2; k <= 10000; ++k) {
	    d__1 = 1. / k;
	    s = pow_di(&d__1, &m);
	    r2 += s;
	    if (s < 1e-15) {
		goto L20;
	    }
/* L10: */
	}
L20:
	bn[m] = r1 * r2;
    }
    return 0;
} /* bernob_ */

/*       ********************************** */
/* Subroutine */ int bernoa_(integer *n, doublereal *bn)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, m;
    static doublereal r__, s;


/*       ====================================== */
/*       Purpose: Compute Bernoulli number Bn */
/*       Input :  n --- Serial number */
/*       Output:  BN(n) --- Bn */
/*       ====================================== */

    bn[0] = 1.;
    bn[1] = -.5;
    i__1 = *n;
    for (m = 2; m <= i__1; ++m) {
	s = -(1. / (m + 1.) - .5);
	i__2 = m - 1;
	for (k = 2; k <= i__2; ++k) {
	    r__ = 1.;
	    i__3 = k;
	    for (j = 2; j <= i__3; ++j) {
/* L10: */
		r__ = r__ * (j + m - k) / j;
	    }
/* L20: */
	    s -= r__ * bn[k];
	}
/* L30: */
	bn[m] = s;
    }
    i__1 = *n;
    for (m = 3; m <= i__1; m += 2) {
/* L40: */
	bn[m] = 0.;
    }
    return 0;
} /* bernoa_ */

/*       ********************************** */
/* Subroutine */ int qstar_(integer *m, integer *n, doublereal *c__, 
	doublereal *ck, doublereal *ck1, doublereal *qs, doublereal *qt)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, k, l;
    static doublereal r__, s, ap[200];
    static integer ip;
    static doublereal sk, qs0;


/*       ========================================================= */
/*       Purpose: Compute Q*mn(-ic) for oblate radial functions */
/*                with a small argument */
/*       ========================================================= */

    /* Parameter adjustments */
    --ck;

    /* Function Body */
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
/* Computing 2nd power */
    d__1 = ck[1];
    r__ = 1. / (d__1 * d__1);
    ap[0] = r__;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = 0.;
	i__2 = i__;
	for (l = 1; l <= i__2; ++l) {
	    sk = 0.;
	    i__3 = l;
	    for (k = 0; k <= i__3; ++k) {
/* L10: */
		sk += ck[k + 1] * ck[l - k + 1];
	    }
/* L15: */
	    s += sk * ap[i__ - l];
	}
/* L20: */
	ap[i__] = -r__ * s;
    }
    qs0 = ap[*m];
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
	r__ = 1.;
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L25: */
/* Computing 2nd power */
	    d__1 = k * 2.;
	    r__ = r__ * (k * 2. + ip) * (k * 2. - 1. + ip) / (d__1 * d__1);
	}
/* L30: */
	qs0 += ap[*m - l] * r__;
    }
    *qs = pow_ii(&c_n1, &ip) * *ck1 * (*ck1 * qs0) / *c__;
    *qt = -2. / *ck1 * *qs;
    return 0;
} /* qstar_ */

/*       ********************************** */
/* Subroutine */ int cv0_(integer *kd, integer *m, doublereal *q, doublereal *
	a0)
{
    static doublereal q2;
    extern /* Subroutine */ int cvql_(integer *, integer *, doublereal *, 
	    doublereal *), cvqm_(integer *, doublereal *, doublereal *);


/*       ===================================================== */
/*       Purpose: Compute the initial characteristic value of */
/*                Mathieu functions for m ≤ 12  or q ≤ 300 or */
/*                q ≥ m*m */
/*       Input :  m  --- Order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions */
/*       Output:  A0 --- Characteristic value */
/*       Routines called: */
/*             (1) CVQM for computing initial characteristic */
/*                 value for q ≤ 3*m */
/*             (2) CVQL for computing initial characteristic */
/*                 value for q ≥ m*m */
/*       ==================================================== */

    q2 = *q * *q;
    if (*m == 0) {
	if (*q <= 1.f) {
	    *a0 = (((q2 * .0036392f - .0125868f) * q2 + .0546875f) * q2 - .5f)
		     * q2;
	} else if (*q <= 10.f) {
	    *a0 = ((*q * .003999267 - .09638957) * *q - .88297f) * *q + 
		    .5542818f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 1) {
	if (*q <= 1.f && *kd == 2) {
	    *a0 = (((*q * -6.51e-4f - .015625f) * *q - .125f) * *q + 1.f) * *
		    q + 1.f;
	} else if (*q <= 1.f && *kd == 3) {
	    *a0 = (((*q * -6.51e-4f + .015625f) * *q - .125f) * *q - 1.f) * *
		    q + 1.f;
	} else if (*q <= 10.f && *kd == 2) {
	    *a0 = (((*q * -4.94603e-4 + .0192917) * *q - .3089229f) * *q + 
		    1.33372f) * *q + .811752f;
	} else if (*q <= 10.f && *kd == 3) {
	    *a0 = ((*q * .001971096 - .05482465) * *q - 1.152218f) * *q + 
		    1.10427f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 2) {
	if (*q <= 1.f && *kd == 1) {
	    *a0 = (((q2 * -.0036391f + .0125888f) * q2 - .0551939f) * q2 + 
		    .416667f) * q2 + 4.f;
	} else if (*q <= 1.f && *kd == 4) {
	    *a0 = (q2 * 3.617e-4f - .0833333f) * q2 + 4.f;
	} else if (*q <= 15. && *kd == 1) {
	    *a0 = (((*q * 3.200972e-4 - .008667445) * *q - 1.829032e-4) * *q 
		    + .9919999f) * *q + 3.3290504f;
	} else if (*q <= 10.f && *kd == 4) {
	    *a0 = ((*q * .00238446 - .08725329f) * *q - .004732542) * *q + 
		    4.00909f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 3) {
	if (*q <= 1.f && *kd == 2) {
	    *a0 = ((*q * 6.348e-4f + .015625f) * *q + .0625f) * q2 + 9.f;
	} else if (*q <= 1.f && *kd == 3) {
	    *a0 = ((*q * 6.348e-4f - .015625f) * *q + .0625f) * q2 + 9.f;
	} else if (*q <= 20.f && *kd == 2) {
	    *a0 = (((*q * 3.035731e-4 - .01453021) * *q + .19069602f) * *q - 
		    .1039356f) * *q + 8.9449274f;
	} else if (*q <= 15.f && *kd == 3) {
	    *a0 = ((*q * 9.369364e-5 - .03569325f) * *q + .2689874f) * *q + 
		    8.771735f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 4) {
	if (*q <= 1.f && *kd == 1) {
	    *a0 = ((q2 * -2.1e-6f + 5.012e-4f) * q2 + .0333333f) * q2 + 16.f;
	} else if (*q <= 1.f && *kd == 4) {
	    *a0 = ((q2 * 3.7e-6f - 3.669e-4f) * q2 + .0333333f) * q2 + 16.f;
	} else if (*q <= 25.f && *kd == 1) {
	    *a0 = (((*q * 1.076676e-4 - .0079684875) * *q + .17344854f) * *q 
		    - .5924058f) * *q + 16.620847f;
	} else if (*q <= 20.f && *kd == 4) {
	    *a0 = ((*q * -7.08719e-4 + .0038216144) * *q + .1907493f) * *q + 
		    15.744f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 5) {
	if (*q <= 1.f && *kd == 2) {
	    *a0 = ((*q * 6.8e-6f + 1.42e-5f) * q2 + .0208333f) * q2 + 25.f;
	} else if (*q <= 1.f && *kd == 3) {
	    *a0 = ((*q * -6.8e-6f + 1.42e-5f) * q2 + .0208333f) * q2 + 25.f;
	} else if (*q <= 35.f && *kd == 2) {
	    *a0 = (((*q * 2.238231e-5 - .002983416) * *q + .10706975f) * *q - 
		    .600205f) * *q + 25.93515f;
	} else if (*q <= 25.f && *kd == 3) {
	    *a0 = ((*q * -7.425364e-4 + .0218225) * *q + .0416399) * *q + 
		    24.897f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 6) {
	if (*q <= 1.f) {
	    *a0 = (q2 * 4e-7 + .0142857f) * q2 + 36.f;
	} else if (*q <= 40.f && *kd == 1) {
	    *a0 = (((*q * -1.66846e-5 + 4.80263e-4) * *q + .0253998) * *q - 
		    .181233f) * *q + 36.423f;
	} else if (*q <= 35.f && *kd == 4) {
	    *a0 = ((*q * -4.57146e-4 + .0216609) * *q - .02349616) * *q + 
		    35.99251f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m == 7) {
	if (*q <= 10.f) {
	    cvqm_(m, q, a0);
	} else if (*q <= 50.f && *kd == 2) {
	    *a0 = (((*q * -1.411114e-5 + 9.730514e-4) * *q - .003097887) * *q 
		    + .03533597) * *q + 49.0547f;
	} else if (*q <= 40.f && *kd == 3) {
	    *a0 = ((*q * -3.043872e-4 + .0205511) * *q - .0916292) * *q + 
		    49.19035f;
	} else {
	    cvql_(kd, m, q, a0);
	}
    } else if (*m >= 8) {
	if (*q <= *m * 3.f) {
	    cvqm_(m, q, a0);
	} else if (*q > (doublereal) (*m * *m)) {
	    cvql_(kd, m, q, a0);
	} else {
	    if (*m == 8 && *kd == 1) {
		*a0 = (((*q * 8.634308e-6 - .002100289) * *q + .169072f) * *q 
			- 4.64336f) * *q + 109.4211f;
	    } else if (*m == 8 && *kd == 4) {
		*a0 = ((*q * -6.7842e-5 + .0022057) * *q + .48296f) * *q + 
			56.59f;
	    } else if (*m == 9 && *kd == 2) {
		*a0 = (((*q * 2.906435e-6 - .001019893) * *q + .1101965f) * *
			q - 3.821851f) * *q + 127.6098f;
	    } else if (*m == 9 && *kd == 3) {
		*a0 = ((*q * -9.577289e-5 + .01043839f) * *q + .06588934f) * *
			q + 78.0198f;
	    } else if (*m == 10 && *kd == 1) {
		*a0 = (((*q * 5.44927e-7 - 3.926119e-4) * *q + .0612099f) * *
			q - 2.600805f) * *q + 138.1923f;
	    } else if (*m == 10 && *kd == 4) {
		*a0 = ((*q * -7.660143e-5 + .01132506f) * *q - .09746023f) * *
			q + 99.29494f;
	    } else if (*m == 11 && *kd == 2) {
		*a0 = (((*q * -5.67615e-7 + 7.152722e-6) * *q + .01920291f) * 
			*q - 1.081583f) * *q + 140.88f;
	    } else if (*m == 11 && *kd == 3) {
		*a0 = ((*q * -6.310551e-5 + .0119247f) * *q - .2681195f) * *q 
			+ 123.667f;
	    } else if (*m == 12 && *kd == 1) {
		*a0 = (((*q * -2.38351e-7 - 2.90139e-5) * *q + .02023088f) * *
			q - 1.289f) * *q + 171.2723f;
	    } else if (*m == 12 && *kd == 4) {
		*a0 = (((*q * 3.08902e-7 - 1.577869e-4) * *q + .0247911f) * *
			q - 1.05454f) * *q + 161.471f;
	    }
	}
    }
    return 0;
} /* cv0_ */

/*       ********************************** */
/* Subroutine */ int cvqm_(integer *m, doublereal *q, doublereal *a0)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal hm1, hm3, hm5;


/*       ===================================================== */
/*       Purpose: Compute the characteristic value of Mathieu */
/*                functions for q ≤ m*m */
/*       Input :  m  --- Order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions */
/*       Output:  A0 --- Initial characteristic value */
/*       ===================================================== */

    hm1 = *q * .5f / (*m * *m - 1.f);
/* Computing 3rd power */
    d__1 = hm1;
    hm3 = d__1 * (d__1 * d__1) * .25f / (*m * *m - 4.f);
    hm5 = hm1 * hm3 * *q / ((*m * *m - 1.f) * (*m * *m - 9.f));
/* Computing 4th power */
    i__1 = *m, i__1 *= i__1;
    *a0 = *m * *m + *q * (hm1 + (*m * 5.f * *m + 7.f) * hm3 + (i__1 * i__1 * 
	    9.f + *m * 58.f * *m + 29.f) * hm5);
    return 0;
} /* cvqm_ */

/*       ********************************** */
/* Subroutine */ int cvql_(integer *kd, integer *m, doublereal *q, doublereal 
	*a0)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal w, c1, d1, d2, d3, d4, p1, p2, w2, w3, w4, w6, cv1, cv2;


/*       ======================================================== */
/*       Purpose: Compute the characteristic value of Mathieu */
/*                functions  for q ≥ 3m */
/*       Input :  m  --- Order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions */
/*       Output:  A0 --- Initial characteristic value */
/*       ======================================================== */

    w = 0.;
    if (*kd == 1 || *kd == 2) {
	w = *m * 2. + 1.;
    }
    if (*kd == 3 || *kd == 4) {
	w = *m * 2. - 1.;
    }
    w2 = w * w;
    w3 = w * w2;
    w4 = w2 * w2;
    w6 = w2 * w4;
    d1 = 34.f / w2 + 5.f + 9.f / w4;
    d2 = (410.f / w2 + 33.f + 405.f / w4) / w;
    d3 = (1260.f / w2 + 63.f + 2943.f / w4 + 486.f / w6) / w2;
    d4 = (15617.f / w2 + 527.f + 69001.f / w4 + 41607.f / w6) / w3;
    c1 = 128.f;
    p2 = *q / w4;
    p1 = sqrt(p2);
    cv1 = *q * -2.f + w * 2.f * sqrt(*q) - (w2 + 1.f) / 8.f;
    cv2 = w + 3.f / w + d1 / (p1 * 32.f) + d2 / (c1 * 8.f * p2);
    cv2 = cv2 + d3 / (c1 * 64.f * p1 * p2) + d4 / (c1 * 16.f * c1 * p2 * p2);
    *a0 = cv1 - cv2 / (c1 * p1);
    return 0;
} /* cvql_ */

integer msta1_(doublereal *x, integer *mp)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static doublereal f, a0, f0, f1;
    static integer n0, n1, nn, it;
    extern doublereal envj_(integer *, doublereal *);


/*       =================================================== */
/*       Purpose: Determine the starting point for backward */
/*                recurrence such that the magnitude of */
/*                Jn(x) at that point is about 10^(-MP) */
/*       Input :  x     --- Argument of Jn(x) */
/*                MP    --- Value of magnitude */
/*       Output:  MSTA1 --- Starting point */
/*       =================================================== */

    a0 = abs(*x);
    n0 = (integer) (a0 * 1.1) + 1;
    f0 = envj_(&n0, &a0) - *mp;
    n1 = n0 + 5;
    f1 = envj_(&n1, &a0) - *mp;
    for (it = 1; it <= 20; ++it) {
	nn = (integer) (n1 - (n1 - n0) / (1. - f0 / f1));
	f = envj_(&nn, &a0) - *mp;
	if ((i__1 = nn - n1, abs(i__1)) < 1) {
	    goto L20;
	}
	n0 = n1;
	f0 = f1;
	n1 = nn;
/* L10: */
	f1 = f;
    }
L20:
    ret_val = nn;
    return ret_val;
} /* msta1_ */

integer msta2_(doublereal *x, integer *n, integer *mp)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static doublereal f, a0, f0, f1;
    static integer n0, n1, nn, it;
    static doublereal obj, ejn, hmp;
    extern doublereal envj_(integer *, doublereal *);


/*       =================================================== */
/*       Purpose: Determine the starting point for backward */
/*                recurrence such that all Jn(x) has MP */
/*                significant digits */
/*       Input :  x  --- Argument of Jn(x) */
/*                n  --- Order of Jn(x) */
/*                MP --- Significant digit */
/*       Output:  MSTA2 --- Starting point */
/*       =================================================== */

    a0 = abs(*x);
    hmp = *mp * .5;
    ejn = envj_(n, &a0);
    if (ejn <= hmp) {
	obj = (doublereal) (*mp);
	n0 = (integer) (a0 * 1.1f) + 1;
    } else {
	obj = hmp + ejn;
	n0 = *n;
    }
    f0 = envj_(&n0, &a0) - obj;
    n1 = n0 + 5;
    f1 = envj_(&n1, &a0) - obj;
    for (it = 1; it <= 20; ++it) {
	nn = (integer) (n1 - (n1 - n0) / (1. - f0 / f1));
	f = envj_(&nn, &a0) - obj;
	if ((i__1 = nn - n1, abs(i__1)) < 1) {
	    goto L20;
	}
	n0 = n1;
	f0 = f1;
	n1 = nn;
/* L10: */
	f1 = f;
    }
L20:
    ret_val = nn + 10;
    return ret_val;
} /* msta2_ */

doublereal envj_(integer *n, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_lg10(doublereal *);

    d__1 = *n * 6.28;
    d__2 = *x * 1.36 / *n;
    ret_val = d_lg10(&d__1) * .5 - *n * d_lg10(&d__2);
    return ret_val;
} /* envj_ */

/*       ********************************** */
/* Subroutine */ int ittjyb_(doublereal *x, doublereal *ttj, doublereal *tty)
{
    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal), sqrt(doublereal)
	    ;

    /* Local variables */
    static doublereal t, e0, f0, g0, t1, x1, el, pi, xt;


/*       ========================================================== */
/*       Purpose: Integrate [1-J0(t)]/t with respect to t from 0 */
/*                to x, and Y0(t)/t with respect to t from x to ∞ */
/*       Input :  x   --- Variable in the limits  ( x ≥ 0 ) */
/*       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x */
/*                TTY --- Integration of Y0(t)/t from x to ∞ */
/*       ========================================================== */

    pi = 3.141592653589793;
    el = .5772156649015329;
    if (*x == 0.) {
	*ttj = 0.;
	*tty = -1e300;
    } else if (*x <= 4.) {
	x1 = *x / 4.;
	t = x1 * x1;
	*ttj = ((((((t * 3.5817e-5 - 6.39765e-4) * t + .007092535) * t - 
		.055544803) * t + .296292677) * t - .999999326) * t + 
		1.999999936) * t;
	*tty = (((((((t * -3.546e-6 + 7.6217e-5) * t - .001059499) * t + 
		.010787555) * t - .07810271) * t + .377255736) * t - 
		1.114084491) * t + 1.909859297) * t;
	e0 = el + log(*x / 2.);
	*tty = pi / 6. + e0 / pi * (*ttj * 2. - e0) - *tty;
    } else if (*x <= 8.) {
	xt = *x + pi * .25;
	t1 = 4. / *x;
	t = t1 * t1;
	f0 = (((((t * .0145369 - .0666297) * t + .1341551) * t - .1647797) * 
		t + .1608874) * t - .2021547) * t + .7977506;
	g0 = ((((((t * .0160672 - .0759339) * t + .1576116) * t - .1960154) * 
		t + .1797457) * t - .1702778) * t + .3235819) * t1;
	*ttj = (f0 * cos(xt) + g0 * sin(xt)) / (sqrt(*x) * *x);
	*ttj = *ttj + el + log(*x / 2.);
	*tty = (f0 * sin(xt) - g0 * cos(xt)) / (sqrt(*x) * *x);
    } else {
	t = 8. / *x;
	xt = *x + pi * .25;
	f0 = (((((t * .0018118 - .0091909) * t + .017033) * t - 9.394e-4) * t 
		- .051445) * t - 1.1e-6) * t + .7978846;
	g0 = (((((t * -.0023731 + .0059842) * t + .0024437) * t - .0233178) * 
		t + 5.95e-5) * t + .1620695) * t;
	*ttj = (f0 * cos(xt) + g0 * sin(xt)) / (sqrt(*x) * *x) + el + log(*x /
		 2.);
	*tty = (f0 * sin(xt) - g0 * cos(xt)) / (sqrt(*x) * *x);
    }
    return 0;
} /* ittjyb_ */

/*       ********************************** */
/* Subroutine */ int ittjya_(doublereal *x, doublereal *ttj, doublereal *tty)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer k, l;
    static doublereal r__, t, a0, b1, e0, g0, g1, r0, r1, r2, el, pi, xk, rs, 
	    px, qx, vt, bj0, bj1, by0, by1;


/*       ========================================================= */
/*       Purpose: Integrate [1-J0(t)]/t with respect to t from 0 */
/*                to x, and Y0(t)/t with respect to t from x to ∞ */
/*       Input :  x   --- Variable in the limits  ( x ≥ 0 ) */
/*       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x */
/*                TTY --- Integration of Y0(t)/t from x to ∞ */
/*       ========================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    if (*x == 0.) {
	*ttj = 0.;
	*tty = -1e300;
    } else if (*x <= 20.) {
	*ttj = 1.;
	r__ = 1.;
	for (k = 2; k <= 100; ++k) {
	    r__ = r__ * -.25 * (k - 1.) / (k * k * k) * *x * *x;
	    *ttj += r__;
	    if (abs(r__) < abs(*ttj) * 1e-12) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*ttj = *ttj * .125 * *x * *x;
	e0 = (pi * pi / 6. - el * el) * .5 - (log(*x / 2.) * .5 + el) * log(*
		x / 2.);
	b1 = el + log(*x / 2.) - 1.5;
	rs = 1.;
	r__ = -1.;
	for (k = 2; k <= 100; ++k) {
	    r__ = r__ * -.25 * (k - 1.) / (k * k * k) * *x * *x;
	    rs += 1. / k;
	    r2 = r__ * (rs + 1. / (k * 2.) - (el + log(*x / 2.)));
	    b1 += r2;
	    if (abs(r2) < abs(b1) * 1e-12) {
		goto L25;
	    }
/* L20: */
	}
L25:
	*tty = 2. / pi * (e0 + *x * .125 * *x * b1);
    } else {
	a0 = sqrt(2. / (pi * *x));
	bj0 = 0.;
	by0 = 0.;
	bj1 = 0.;
	for (l = 0; l <= 1; ++l) {
	    vt = l * 4. * l;
	    px = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 14; ++k) {
/* Computing 2nd power */
		d__1 = k * 4. - 3.;
/* Computing 2nd power */
		d__2 = k * 4. - 1.;
		r__ = r__ * -.0078125 * (vt - d__1 * d__1) / (*x * k) * (vt - 
			d__2 * d__2) / ((k * 2. - 1.) * *x);
		px += r__;
		if (abs(r__) < abs(px) * 1e-12) {
		    goto L35;
		}
/* L30: */
	    }
L35:
	    qx = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 14; ++k) {
/* Computing 2nd power */
		d__1 = k * 4. - 1.;
/* Computing 2nd power */
		d__2 = k * 4. + 1.;
		r__ = r__ * -.0078125 * (vt - d__1 * d__1) / (*x * k) * (vt - 
			d__2 * d__2) / (k * 2. + 1.) / *x;
		qx += r__;
		if (abs(r__) < abs(qx) * 1e-12) {
		    goto L45;
		}
/* L40: */
	    }
L45:
	    qx = (vt - 1.) * .125 / *x * qx;
	    xk = *x - (l * .5 + .25) * pi;
	    bj1 = a0 * (px * cos(xk) - qx * sin(xk));
	    by1 = a0 * (px * sin(xk) + qx * cos(xk));
	    if (l == 0) {
		bj0 = bj1;
		by0 = by1;
	    }
/* L50: */
	}
	t = 2. / *x;
	g0 = 1.;
	r0 = 1.;
	for (k = 1; k <= 10; ++k) {
	    r0 = -k * k * t * t * r0;
/* L55: */
	    g0 += r0;
	}
	g1 = 1.;
	r1 = 1.;
	for (k = 1; k <= 10; ++k) {
	    r1 = -k * (k + 1.) * t * t * r1;
/* L60: */
	    g1 += r1;
	}
	*ttj = g1 * 2. * bj0 / (*x * *x) - g0 * bj1 / *x + el + log(*x / 2.);
	*tty = g1 * 2. * by0 / (*x * *x) - g0 * by1 / *x;
    }
    return 0;
} /* ittjya_ */

/*       ********************************** */
/* Subroutine */ int cjylv_(doublereal *v, doublecomplex *z__, doublecomplex *
	cbjv, doublecomplex *cdjv, doublecomplex *cbyv, doublecomplex *cdyv)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    void z_sqrt(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *), z_log(doublecomplex *, 
	    doublecomplex *), pow_zi(doublecomplex *, doublecomplex *, 
	    integer *);
    double pow_di(doublereal *, integer *);
    void z_exp(doublecomplex *, doublecomplex *);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal a[91];
    static integer i__, k, l, l0;
    static doublereal v0;
    static doublecomplex cf[12];
    static integer lf;
    static doublecomplex ct;
    static integer km;
    static doublereal pi, vr;
    static doublecomplex ct2, cfj;
    extern /* Subroutine */ int cjk_(integer *, doublereal *);
    static doublecomplex csj, cfy, cws, csy, ceta;


/*       =================================================== */
/*       Purpose: Compute Bessel functions Jv(z) and Yv(z) */
/*                and their derivatives with a complex */
/*                argument and a large order */
/*       Input:   v --- Order of Jv(z) and Yv(z) */
/*                z --- Complex argument */
/*       Output:  CBJV --- Jv(z) */
/*                CDJV --- Jv'(z) */
/*                CBYV --- Yv(z) */
/*                CDYV --- Yv'(z) */
/*       Routine called: */
/*                CJK to compute the expansion coefficients */
/*       =================================================== */

    km = 12;
    cjk_(&km, a);
    pi = 3.141592653589793;
    for (l = 1; l >= 0; --l) {
	v0 = *v - l;
	z__4.r = z__->r / v0, z__4.i = z__->i / v0;
	z__5.r = z__->r / v0, z__5.i = z__->i / v0;
	z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, z__3.i = z__4.r * z__5.i 
		+ z__4.i * z__5.r;
	z__2.r = 1. - z__3.r, z__2.i = -z__3.i;
	z_sqrt(&z__1, &z__2);
	cws.r = z__1.r, cws.i = z__1.i;
	z__4.r = z__->r / v0, z__4.i = z__->i / v0;
	z__5.r = cws.r + 1., z__5.i = cws.i;
	z_div(&z__3, &z__4, &z__5);
	z_log(&z__2, &z__3);
	z__1.r = cws.r + z__2.r, z__1.i = cws.i + z__2.i;
	ceta.r = z__1.r, ceta.i = z__1.i;
	z_div(&z__1, &c_b6, &cws);
	ct.r = z__1.r, ct.i = z__1.i;
	z__1.r = ct.r * ct.r - ct.i * ct.i, z__1.i = ct.r * ct.i + ct.i * 
		ct.r;
	ct2.r = z__1.r, ct2.i = z__1.i;
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
	    l0 = k * (k + 1) / 2 + 1;
	    lf = l0 + k;
	    i__2 = k - 1;
	    i__3 = lf - 1;
	    cf[i__2].r = a[i__3], cf[i__2].i = 0.;
	    i__2 = l0;
	    for (i__ = lf - 1; i__ >= i__2; --i__) {
/* L10: */
		i__3 = k - 1;
		i__4 = k - 1;
		z__2.r = cf[i__4].r * ct2.r - cf[i__4].i * ct2.i, z__2.i = cf[
			i__4].r * ct2.i + cf[i__4].i * ct2.r;
		i__5 = i__ - 1;
		z__1.r = z__2.r + a[i__5], z__1.i = z__2.i;
		cf[i__3].r = z__1.r, cf[i__3].i = z__1.i;
	    }
/* L15: */
	    i__3 = k - 1;
	    i__4 = k - 1;
	    pow_zi(&z__2, &ct, &k);
	    z__1.r = cf[i__4].r * z__2.r - cf[i__4].i * z__2.i, z__1.i = cf[
		    i__4].r * z__2.i + cf[i__4].i * z__2.r;
	    cf[i__3].r = z__1.r, cf[i__3].i = z__1.i;
	}
	vr = 1. / v0;
	csj.r = 1., csj.i = 0.;
	i__3 = km;
	for (k = 1; k <= i__3; ++k) {
/* L20: */
	    i__4 = k - 1;
	    d__1 = pow_di(&vr, &k);
	    z__2.r = d__1 * cf[i__4].r, z__2.i = d__1 * cf[i__4].i;
	    z__1.r = csj.r + z__2.r, z__1.i = csj.i + z__2.i;
	    csj.r = z__1.r, csj.i = z__1.i;
	}
	d__1 = pi * 2. * v0;
	z__4.r = ct.r / d__1, z__4.i = ct.i / d__1;
	z_sqrt(&z__3, &z__4);
	z__6.r = v0 * ceta.r, z__6.i = v0 * ceta.i;
	z_exp(&z__5, &z__6);
	z__2.r = z__3.r * z__5.r - z__3.i * z__5.i, z__2.i = z__3.r * z__5.i 
		+ z__3.i * z__5.r;
	z__1.r = z__2.r * csj.r - z__2.i * csj.i, z__1.i = z__2.r * csj.i + 
		z__2.i * csj.r;
	cbjv->r = z__1.r, cbjv->i = z__1.i;
	if (l == 1) {
	    cfj.r = cbjv->r, cfj.i = cbjv->i;
	}
	csy.r = 1., csy.i = 0.;
	i__4 = km;
	for (k = 1; k <= i__4; ++k) {
/* L25: */
	    i__3 = pow_ii(&c_n1, &k);
	    i__1 = k - 1;
	    d__1 = (doublereal) i__3;
	    z__3.r = d__1 * cf[i__1].r, z__3.i = d__1 * cf[i__1].i;
	    d__2 = pow_di(&vr, &k);
	    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
	    z__1.r = csy.r + z__2.r, z__1.i = csy.i + z__2.i;
	    csy.r = z__1.r, csy.i = z__1.i;
	}
	z__6.r = ct.r * 2., z__6.i = ct.i * 2.;
	d__1 = pi * v0;
	z__5.r = z__6.r / d__1, z__5.i = z__6.i / d__1;
	z_sqrt(&z__4, &z__5);
	z__3.r = -z__4.r, z__3.i = -z__4.i;
	d__2 = -v0;
	z__8.r = d__2 * ceta.r, z__8.i = d__2 * ceta.i;
	z_exp(&z__7, &z__8);
	z__2.r = z__3.r * z__7.r - z__3.i * z__7.i, z__2.i = z__3.r * z__7.i 
		+ z__3.i * z__7.r;
	z__1.r = z__2.r * csy.r - z__2.i * csy.i, z__1.i = z__2.r * csy.i + 
		z__2.i * csy.r;
	cbyv->r = z__1.r, cbyv->i = z__1.i;
	if (l == 1) {
	    cfy.r = cbyv->r, cfy.i = cbyv->i;
	}
/* L30: */
    }
    d__1 = -(*v);
    z__4.r = d__1, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbjv->r - z__3.i * cbjv->i, z__2.i = z__3.r * cbjv->i + 
	    z__3.i * cbjv->r;
    z__1.r = z__2.r + cfj.r, z__1.i = z__2.i + cfj.i;
    cdjv->r = z__1.r, cdjv->i = z__1.i;
    d__1 = -(*v);
    z__4.r = d__1, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbyv->r - z__3.i * cbyv->i, z__2.i = z__3.r * cbyv->i + 
	    z__3.i * cbyv->r;
    z__1.r = z__2.r + cfy.r, z__1.i = z__2.i + cfy.i;
    cdyv->r = z__1.r, cdyv->i = z__1.i;
    return 0;
} /* cjylv_ */

/*       ********************************** */
/* Subroutine */ int rmn2l_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, doublereal *df, integer *kd, doublereal *r2f, 
	doublereal *r2d, integer *id)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), d_lg10(doublereal *);

    /* Local variables */
    static integer j, k, l;
    static doublereal r__, a0, b0, r0;
    static integer lg, ip, nm;
    static doublereal cx, dy[252];
    static integer np;
    static doublereal sw, sy[252];
    static integer id1, id2, nm1, nm2;
    static doublereal reg, eps, suc, sud, eps1, eps2;
    extern /* Subroutine */ int sphy_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*       ======================================================== */
/*       Purpose: Compute prolate and oblate spheroidal radial */
/*                functions of the second kind for given m, n, */
/*                c and a large cx */
/*       Routine called: */
/*                SPHY for computing the spherical Bessel */
/*                functions of the second kind */
/*       ======================================================== */

    /* Parameter adjustments */
    --df;

    /* Function Body */
    eps = 1e-14;
    ip = 1;
    nm1 = (*n - *m) / 2;
    if (*n - *m == nm1 << 1) {
	ip = 0;
    }
    nm = nm1 + 25 + (integer) (*c__);
    reg = 1.;
    if (*m + nm > 80) {
	reg = 1e-200;
    }
    nm2 = (nm << 1) + *m;
    cx = *c__ * *x;
    sphy_(&nm2, &cx, &nm2, sy, dy);
    r0 = reg;
    i__1 = (*m << 1) + ip;
    for (j = 1; j <= i__1; ++j) {
/* L10: */
	r0 *= j;
    }
    r__ = r0;
    suc = r__ * df[1];
    sw = 0.;
    i__1 = nm;
    for (k = 2; k <= i__1; ++k) {
	r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k + ip 
		- 1.5);
	suc += r__ * df[k];
	if (k > nm1 && (d__1 = suc - sw, abs(d__1)) < abs(suc) * eps) {
	    goto L20;
	}
/* L15: */
	sw = suc;
    }
L20:
    d__1 = 1. - *kd / (*x * *x);
    d__2 = *m * .5;
    a0 = pow_dd(&d__1, &d__2) / suc;
    *r2f = 0.;
    eps1 = 0.;
    np = 0;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	l = (k << 1) + *m - *n - 2 + ip;
	lg = 1;
	if (l != l / 4 << 2) {
	    lg = -1;
	}
	if (k == 1) {
	    r__ = r0;
	} else {
	    r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k 
		    + ip - 1.5);
	}
	np = *m + (k << 1) - 2 + ip;
	*r2f += lg * r__ * (df[k] * sy[np]);
	eps1 = (d__1 = *r2f - sw, abs(d__1));
	if (k > nm1 && eps1 < abs(*r2f) * eps) {
	    goto L55;
	}
/* L50: */
	sw = *r2f;
    }
L55:
    d__1 = eps1 / abs(*r2f) + eps;
    id1 = (integer) d_lg10(&d__1);
    *r2f *= a0;
    if (np >= nm2) {
	*id = 10;
	return 0;
    }
    b0 = *kd * *m / pow_dd(x, &c_b149) / (1.f - *kd / (*x * *x)) * *r2f;
    sud = 0.;
    eps2 = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	l = (k << 1) + *m - *n - 2 + ip;
	lg = 1;
	if (l != l / 4 << 2) {
	    lg = -1;
	}
	if (k == 1) {
	    r__ = r0;
	} else {
	    r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k 
		    + ip - 1.5);
	}
	np = *m + (k << 1) - 2 + ip;
	sud += lg * r__ * (df[k] * dy[np]);
	eps2 = (d__1 = sud - sw, abs(d__1));
	if (k > nm1 && eps2 < abs(sud) * eps) {
	    goto L65;
	}
/* L60: */
	sw = sud;
    }
L65:
    *r2d = b0 + a0 * *c__ * sud;
    d__1 = eps2 / abs(sud) + eps;
    id2 = (integer) d_lg10(&d__1);
    *id = max(id1,id2);
    return 0;
} /* rmn2l_ */

/*       ********************************** */
/* Subroutine */ int psi_spec__(doublereal *x, doublereal *ps)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, n;
    static doublereal s, a1, a2, a3, a4, a5, a6, a7, a8, x2, el, xa, pi;


/*       ====================================== */
/*       Purpose: Compute Psi function */
/*       Input :  x  --- Argument of psi(x) */
/*       Output:  PS --- psi(x) */
/*       ====================================== */

    xa = abs(*x);
    pi = 3.141592653589793;
    el = .5772156649015329;
    s = 0.;
    if (*x == (doublereal) ((integer) (*x)) && *x <= 0.f) {
	*ps = 1e300;
	return 0;
    } else if (xa == (doublereal) ((integer) xa)) {
	n = (integer) xa;
	i__1 = n - 1;
	for (k = 1; k <= i__1; ++k) {
/* L10: */
	    s += 1. / k;
	}
	*ps = -el + s;
    } else if (xa + .5f == (doublereal) ((integer) (xa + .5f))) {
	n = (integer) (xa - .5f);
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
/* L20: */
	    s += 1.f / (k * 2. - 1.);
	}
	*ps = -el + s * 2. - 1.386294361119891;
    } else {
	if (xa < 10.f) {
	    n = 10 - (integer) xa;
	    i__1 = n - 1;
	    for (k = 0; k <= i__1; ++k) {
/* L30: */
		s += 1. / (xa + k);
	    }
	    xa += n;
	}
	x2 = 1. / (xa * xa);
	a1 = -.08333333333333;
	a2 = .0083333333333333333;
	a3 = -.0039682539682539683;
	a4 = .0041666666666666667;
	a5 = -.0075757575757575758;
	a6 = .021092796092796093;
	a7 = -.083333333333333333;
	a8 = .4432598039215686;
	*ps = log(xa) - .5 / xa + x2 * (((((((a8 * x2 + a7) * x2 + a6) * x2 + 
		a5) * x2 + a4) * x2 + a3) * x2 + a2) * x2 + a1);
	*ps -= s;
    }
    if (*x < 0.f) {
	*ps = *ps - pi * cos(pi * *x) / sin(pi * *x) - 1. / *x;
    }
    return 0;
} /* psi_spec__ */

/*       ********************************** */
/* Subroutine */ int cva2_(integer *kd, integer *m, doublereal *q, doublereal 
	*a)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal a1, a2, q1, q2;
    static integer nn;
    static doublereal qq;
    extern /* Subroutine */ int cv0_(integer *, integer *, doublereal *, 
	    doublereal *);
    static integer ndiv;
    extern /* Subroutine */ int cvql_(integer *, integer *, doublereal *, 
	    doublereal *), cvqm_(integer *, doublereal *, doublereal *);
    static integer iflag;
    static doublereal delta;
    extern /* Subroutine */ int refine_(integer *, integer *, doublereal *, 
	    doublereal *);


/*       ====================================================== */
/*       Purpose: Calculate a specific characteristic value of */
/*                Mathieu functions */
/*       Input :  m  --- Order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions */
/*                KD --- Case code */
/*                       KD=1 for cem(x,q)  ( m = 0,2,4,...) */
/*                       KD=2 for cem(x,q)  ( m = 1,3,5,...) */
/*                       KD=3 for sem(x,q)  ( m = 1,3,5,...) */
/*                       KD=4 for sem(x,q)  ( m = 2,4,6,...) */
/*       Output:  A  --- Characteristic value */
/*       Routines called: */
/*             (1) REFINE for finding accurate characteristic */
/*                 value using an iteration method */
/*             (2) CV0 for finding initial characteristic */
/*                 values using polynomial approximation */
/*             (3) CVQM for computing initial characteristic */
/*                 values for q ≤ 3*m */
/*             (3) CVQL for computing initial characteristic */
/*                 values for q ≥ m*m */
/*       ====================================================== */

    if (*m <= 12 || *q <= *m * 3.f || *q > (doublereal) (*m * *m)) {
	cv0_(kd, m, q, a);
	if (*q != 0. && *m != 2) {
	    refine_(kd, m, q, a);
	}
	if (*q > .002 && *m == 2) {
	    refine_(kd, m, q, a);
	}
    } else {
	ndiv = 10;
	delta = (*m - 3.f) * *m / ndiv;
	if (*q - *m * 3.f <= *m * *m - *q) {
L5:
	    nn = (integer) ((*q - *m * 3.f) / delta) + 1;
	    delta = (*q - *m * 3.f) / nn;
	    q1 = *m * 2.f;
	    cvqm_(m, &q1, &a1);
	    q2 = *m * 3.f;
	    cvqm_(m, &q2, &a2);
	    qq = *m * 3.f;
	    i__1 = nn;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		qq += delta;
		*a = (a1 * q2 - a2 * q1 + (a2 - a1) * qq) / (q2 - q1);
		iflag = 1;
		if (i__ == nn) {
		    iflag = -1;
		}
		refine_(kd, m, &qq, a);
		q1 = q2;
		q2 = qq;
		a1 = a2;
		a2 = *a;
/* L10: */
	    }
	    if (iflag == -10) {
		ndiv <<= 1;
		delta = (*m - 3.f) * *m / ndiv;
		goto L5;
	    }
	} else {
L15:
	    nn = (integer) ((*m * *m - *q) / delta) + 1;
	    delta = (*m * *m - *q) / nn;
	    q1 = *m * (*m - 1.f);
	    cvql_(kd, m, &q1, &a1);
	    q2 = (doublereal) (*m * *m);
	    cvql_(kd, m, &q2, &a2);
	    qq = (doublereal) (*m * *m);
	    i__1 = nn;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		qq -= delta;
		*a = (a1 * q2 - a2 * q1 + (a2 - a1) * qq) / (q2 - q1);
		iflag = 1;
		if (i__ == nn) {
		    iflag = -1;
		}
		refine_(kd, m, &qq, a);
		q1 = q2;
		q2 = qq;
		a1 = a2;
		a2 = *a;
/* L20: */
	    }
	    if (iflag == -10) {
		ndiv <<= 1;
		delta = (*m - 3.f) * *m / ndiv;
		goto L15;
	    }
	}
    }
    return 0;
} /* cva2_ */

/*       ********************************** */
/* Subroutine */ int lpmns_(integer *m, integer *n, doublereal *x, doublereal 
	*pm, doublereal *pd)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal x0, pm0, pm1, pm2, pmk;


/*       ======================================================== */
/*       Purpose: Compute associated Legendre functions Pmn(x) */
/*                and Pmn'(x) for a given order */
/*       Input :  x --- Argument of Pmn(x) */
/*                m --- Order of Pmn(x),  m = 0,1,2,...,n */
/*                n --- Degree of Pmn(x), n = 0,1,2,...,N */
/*       Output:  PM(n) --- Pmn(x) */
/*                PD(n) --- Pmn'(x) */
/*       ======================================================== */

    i__1 = *n;
    for (k = 0; k <= i__1; ++k) {
	pm[k] = 0.;
/* L10: */
	pd[k] = 0.;
    }
    if (abs(*x) == 1.) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    if (*m == 0) {
		pm[k] = 1.;
		pd[k] = k * .5 * (k + 1.f);
		if (*x < 0.f) {
		    pm[k] = pow_ii(&c_n1, &k) * pm[k];
		    i__2 = k + 1;
		    pd[k] = pow_ii(&c_n1, &i__2) * pd[k];
		}
	    } else if (*m == 1) {
		pd[k] = 1e300;
	    } else if (*m == 2) {
		pd[k] = (k + 2.f) * -.25 * (k + 1.f) * k * (k - 1.f);
		if (*x < 0.f) {
		    i__2 = k + 1;
		    pd[k] = pow_ii(&c_n1, &i__2) * pd[k];
		}
	    }
/* L15: */
	}
	return 0;
    }
    x0 = (d__1 = 1. - *x * *x, abs(d__1));
    pm0 = 1.;
    pmk = pm0;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	pmk = (k * 2. - 1.) * sqrt(x0) * pm0;
/* L20: */
	pm0 = pmk;
    }
    pm1 = (*m * 2. + 1.) * *x * pm0;
    pm[*m] = pmk;
    pm[*m + 1] = pm1;
    i__1 = *n;
    for (k = *m + 2; k <= i__1; ++k) {
	pm2 = ((k * 2. - 1.) * *x * pm1 - (k + *m - 1.) * pmk) / (k - *m);
	pm[k] = pm2;
	pmk = pm1;
/* L25: */
	pm1 = pm2;
    }
    pd[0] = ((1. - *m) * pm[1] - *x * pm[0]) / (*x * *x - 1.f);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* L30: */
	pd[k] = (k * *x * pm[k] - (k + *m) * pm[k - 1]) / (*x * *x - 1.);
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	pm[k] = pow_ii(&c_n1, m) * pm[k];
/* L35: */
	pd[k] = pow_ii(&c_n1, m) * pd[k];
    }
    return 0;
} /* lpmns_ */

/*       ********************************** */
/* Subroutine */ int cerf_(doublecomplex *z__, doublecomplex *cer, 
	doublecomplex *cder)
{
    /* System generated locals */
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal), exp(doublereal), cos(
	    doublereal), sin(doublereal), cosh(doublereal), sinh(doublereal);
    void z_exp(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, n;
    static doublereal r__, w, x, y, c0, w1, w2, x2, cs, er, pi, ss, ei1, ei2, 
	    er0, er1, er2, eri, eps, err;


/*       ========================================================== */
/*       Purpose: Compute complex Error function erf(z) & erf'(z) */
/*       Input:   z   --- Complex argument of erf(z) */
/*                x   --- Real part of z */
/*                y   --- Imaginary part of z */
/*       Output:  CER --- erf(z) */
/*                CDER --- erf'(z) */
/*       ========================================================== */
    eps = 1e-12;
    pi = 3.141592653589793;
    x = z__->r;
    y = d_imag(z__);
    x2 = x * x;
    if (x <= 3.5) {
	er = 1.;
	r__ = 1.;
	w = 0.;
	for (k = 1; k <= 100; ++k) {
	    r__ = r__ * x2 / (k + .5);
	    er += r__;
	    if ((d__1 = er - w, abs(d__1)) <= eps * abs(er)) {
		goto L15;
	    }
/* L10: */
	    w = er;
	}
L15:
	c0 = 2. / sqrt(pi) * x * exp(-x2);
	er0 = c0 * er;
    } else {
	er = 1.;
	r__ = 1.;
	for (k = 1; k <= 12; ++k) {
	    r__ = -r__ * (k - .5) / x2;
/* L20: */
	    er += r__;
	}
	c0 = exp(-x2) / (x * sqrt(pi));
	er0 = 1. - c0 * er;
    }
    if (y == 0.) {
	err = er0;
	eri = 0.;
    } else {
	cs = cos(x * 2. * y);
	ss = sin(x * 2. * y);
	er1 = exp(-x2) * (1. - cs) / (pi * 2. * x);
	ei1 = exp(-x2) * ss / (pi * 2. * x);
	er2 = 0.;
	w1 = 0.;
	for (n = 1; n <= 100; ++n) {
	    er2 += exp(n * -.25 * n) / (n * n + x2 * 4.) * (x * 2. - x * 2. * 
		    cosh(n * y) * cs + n * sinh(n * y) * ss);
	    if ((d__1 = (er2 - w1) / er2, abs(d__1)) < eps) {
		goto L30;
	    }
/* L25: */
	    w1 = er2;
	}
L30:
	c0 = exp(-x2) * 2. / pi;
	err = er0 + er1 + c0 * er2;
	ei2 = 0.;
	w2 = 0.;
	for (n = 1; n <= 100; ++n) {
	    ei2 += exp(n * -.25 * n) / (n * n + x2 * 4.) * (x * 2. * cosh(n * 
		    y) * ss + n * sinh(n * y) * cs);
	    if ((d__1 = (ei2 - w2) / ei2, abs(d__1)) < eps) {
		goto L40;
	    }
/* L35: */
	    w2 = ei2;
	}
L40:
	eri = ei1 + c0 * ei2;
    }
    z__1.r = err, z__1.i = eri;
    cer->r = z__1.r, cer->i = z__1.i;
    d__1 = 2. / sqrt(pi);
    z__4.r = -z__->r, z__4.i = -z__->i;
    z__3.r = z__4.r * z__->r - z__4.i * z__->i, z__3.i = z__4.r * z__->i + 
	    z__4.i * z__->r;
    z_exp(&z__2, &z__3);
    z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
    cder->r = z__1.r, cder->i = z__1.i;
    return 0;
} /* cerf_ */

/*       ********************************** */
/* Subroutine */ int rswfp_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, doublereal *cv, integer *kf, doublereal *r1f, 
	doublereal *r1d, doublereal *r2f, doublereal *r2d)
{
    static doublereal df[200];
    static integer id, kd;
    extern /* Subroutine */ int rmn1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *)
	    , sdmn_(integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *), rmn2l_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *), rmn2sp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);


/*       ============================================================== */
/*       Purpose: Compute prolate spheriodal radial functions of the */
/*                first and second kinds, and their derivatives */
/*       Input :  m  --- Mode parameter, m = 0,1,2,... */
/*                n  --- Mode parameter, n = m,m+1,m+2,... */
/*                c  --- Spheroidal parameter */
/*                x  --- Argument of radial function ( x > 1.0 ) */
/*                cv --- Characteristic value */
/*                KF --- Function code */
/*                       KF=1 for the first kind */
/*                       KF=2 for the second kind */
/*                       KF=3 for both the first and second kinds */
/*       Output:  R1F --- Radial function of the first kind */
/*                R1D --- Derivative of the radial function of */
/*                        the first kind */
/*                R2F --- Radial function of the second kind */
/*                R2D --- Derivative of the radial function of */
/*                        the second kind */
/*       Routines called: */
/*            (1) SDMN for computing expansion coefficients dk */
/*            (2) RMN1 for computing prolate and oblate radial */
/*                functions of the first kind */
/*            (3) RMN2L for computing prolate and oblate radial */
/*                functions of the second kind for a large argument */
/*            (4) RMN2SP for computing the prolate radial function */
/*                of the second kind for a small argument */
/*       ============================================================== */

    kd = 1;
    sdmn_(m, n, c__, cv, &kd, df);
    if (*kf != 2) {
	rmn1_(m, n, c__, x, df, &kd, r1f, r1d);
    }
    if (*kf > 1) {
	rmn2l_(m, n, c__, x, df, &kd, r2f, r2d, &id);
	if (id > -8) {
	    rmn2sp_(m, n, c__, x, cv, df, &kd, r2f, r2d);
	}
    }
    return 0;
} /* rswfp_ */

/*       ********************************** */
/* Subroutine */ int jyndd_(integer *n, doublereal *x, doublereal *bjn, 
	doublereal *djn, doublereal *fjn, doublereal *byn, doublereal *dyn, 
	doublereal *fyn)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal bj[2], by[2];
    static integer nm;
    extern /* Subroutine */ int jynbh_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute Bessel functions Jn(x) and Yn(x), and */
/*                their first and second derivatives */
/*       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 ) */
/*                n   ---  Order of Jn(x) and Yn(x) */
/*       Output:  BJN ---  Jn(x) */
/*                DJN ---  Jn'(x) */
/*                FJN ---  Jn"(x) */
/*                BYN ---  Yn(x) */
/*                DYN ---  Yn'(x) */
/*                FYN ---  Yn"(x) */
/*       Routines called: */
/*                JYNBH to compute Jn and Yn */
/*       =========================================================== */

    i__1 = *n + 1;
    jynbh_(&i__1, n, x, &nm, bj, by);
/*       Compute derivatives by differentiation formulas */
    *bjn = bj[0];
    *byn = by[0];
    *djn = -bj[1] + *n * bj[0] / *x;
    *dyn = -by[1] + *n * by[0] / *x;
    *fjn = (*n * *n / (*x * *x) - 1.) * *bjn - *djn / *x;
    *fyn = (*n * *n / (*x * *x) - 1.) * *byn - *dyn / *x;
    return 0;
} /* jyndd_ */

/*       ********************************** */
/* Subroutine */ int gam0_(doublereal *x, doublereal *ga)
{
    /* Initialized data */

    static doublereal g[25] = { 1.,.5772156649015329,-.6558780715202538,
	    -.0420026350340952,.1665386113822915,-.0421977345555443,
	    -.009621971527877,.007218943246663,-.0011651675918591,
	    -2.152416741149e-4,1.280502823882e-4,-2.01348547807e-5,
	    -1.2504934821e-6,1.133027232e-6,-2.056338417e-7,6.116095e-9,
	    5.0020075e-9,-1.1812746e-9,1.043427e-10,7.7823e-12,-3.6968e-12,
	    5.1e-13,-2.06e-14,-5.4e-15,1.4e-15 };

    static integer k;
    static doublereal gr;


/*       ================================================ */
/*       Purpose: Compute gamma function Г(x) */
/*       Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 ) */
/*       Output:  GA --- Г(x) */
/*       ================================================ */

    gr = 25.;
    for (k = 24; k >= 1; --k) {
/* L20: */
	gr = gr * *x + g[k - 1];
    }
    *ga = 1. / (gr * *x);
    return 0;
} /* gam0_ */

/*       ********************************** */
/* Subroutine */ int cisib_(doublereal *x, doublereal *ci, doublereal *si)
{
    /* Builtin functions */
    double log(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal x2, fx, gx;


/*       ============================================= */
/*       Purpose: Compute cosine and sine integrals */
/*                Si(x) and Ci(x) ( x ≥ 0 ) */
/*       Input :  x  --- Argument of Ci(x) and Si(x) */
/*       Output:  CI --- Ci(x) */
/*                SI --- Si(x) */
/*       ============================================= */

    x2 = *x * *x;
    if (*x == 0.f) {
	*ci = -1e300;
	*si = 0.;
    } else if (*x <= 1.) {
	*ci = ((((x2 * -3e-8 + 3.1e-6) * x2 - 2.3148e-4) * x2 + .01041667) * 
		x2 - .25f) * x2 + .577215665 + log(*x);
	*si = ((((x2 * 3.1e-7 - 2.834e-5) * x2 + .00166667) * x2 - .05555556) 
		* x2 + 1.f) * *x;
    } else {
	fx = ((((x2 + 38.027264) * x2 + 265.187033) * x2 + 335.67732) * x2 + 
		38.102495) / ((((x2 + 40.021433) * x2 + 322.624911) * x2 + 
		570.23628) * x2 + 157.105423);
	gx = ((((x2 + 42.242855) * x2 + 302.757865) * x2 + 352.018498) * x2 + 
		21.821899) / ((((x2 + 48.196927) * x2 + 482.485984) * x2 + 
		1114.978885) * x2 + 449.690326) / *x;
	*ci = fx * sin(*x) / *x - gx * cos(*x) / *x;
	*si = 1.570796327 - fx * cos(*x) / *x - gx * sin(*x) / *x;
    }
    return 0;
} /* cisib_ */

/*       ********************************** */
/* Subroutine */ int eulera_(integer *n, doublereal *en)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, m;
    static doublereal r__, s;


/*       ====================================== */
/*       Purpose: Compute Euler number En */
/*       Input :  n --- Serial number */
/*       Output:  EN(n) --- En */
/*       ====================================== */

    en[0] = 1.;
    i__1 = *n / 2;
    for (m = 1; m <= i__1; ++m) {
	s = 1.;
	i__2 = m - 1;
	for (k = 1; k <= i__2; ++k) {
	    r__ = 1.;
	    i__3 = k << 1;
	    for (j = 1; j <= i__3; ++j) {
/* L10: */
		r__ = r__ * (m * 2. - k * 2. + j) / j;
	    }
/* L20: */
	    s += r__ * en[k * 2];
	}
/* L30: */
	en[m * 2] = -s;
    }
    return 0;
} /* eulera_ */

/*       ********************************** */
/* Subroutine */ int refine_(integer *kd, integer *m, doublereal *q, 
	doublereal *a)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal f, x, f0, f1, x0, x1, ca;
    static integer mj, it;
    extern /* Subroutine */ int cvf_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    static doublereal eps, delta;


/*       ===================================================== */
/*       Purpose: calculate the accurate characteristic value */
/*                by the secant method */
/*       Input :  m --- Order of Mathieu functions */
/*                q --- Parameter of Mathieu functions */
/*                A --- Initial characteristic value */
/*       Output:  A --- Refineed characteristic value */
/*       Routine called:  CVF for computing the value of F for */
/*                        characteristic equation */
/*       ======================================================== */

    eps = 1e-14;
    mj = *m + 10;
    ca = *a;
    delta = 0.;
    x0 = *a;
    cvf_(kd, m, q, &x0, &mj, &f0);
    x1 = *a * 1.002f;
    cvf_(kd, m, q, &x1, &mj, &f1);
    for (it = 1; it <= 100; ++it) {
	++mj;
	x = x1 - (x1 - x0) / (1. - f0 / f1);
	cvf_(kd, m, q, &x, &mj, &f);
	if ((d__1 = 1.f - x1 / x, abs(d__1)) < eps || f == 0.f) {
	    goto L15;
	}
	x0 = x1;
	f0 = f1;
	x1 = x;
/* L10: */
	f1 = f;
    }
L15:
    *a = x;
    return 0;
} /* refine_ */

/*       ********************************** */
/* Subroutine */ int cisia_(doublereal *x, doublereal *ci, doublereal *si)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, m;
    static doublereal p2, x2, bj[101], el, xa, xf, xg, xr, xs, xa0, xa1, xg1, 
	    xg2, eps, xcs, xss;


/*       ============================================= */
/*       Purpose: Compute cosine and sine integrals */
/*                Si(x) and Ci(x)  ( x ≥ 0 ) */
/*       Input :  x  --- Argument of Ci(x) and Si(x) */
/*       Output:  CI --- Ci(x) */
/*                SI --- Si(x) */
/*       ============================================= */

    p2 = 1.570796326794897;
    el = .5772156649015329;
    eps = 1e-15;
    x2 = *x * *x;
    if (*x == 0.) {
	*ci = -1e300;
	*si = 0.;
    } else if (*x <= 16.) {
	xr = x2 * -.25;
	*ci = el + log(*x) + xr;
	for (k = 2; k <= 40; ++k) {
	    xr = xr * -.5 * (k - 1) / (k * k * ((k << 1) - 1)) * x2;
	    *ci += xr;
	    if (abs(xr) < abs(*ci) * eps) {
		goto L15;
	    }
/* L10: */
	}
L15:
	xr = *x;
	*si = *x;
	for (k = 1; k <= 40; ++k) {
	    xr = xr * -.5 * ((k << 1) - 1) / k / ((k << 2) * k + (k << 2) + 1)
		     * x2;
	    *si += xr;
	    if (abs(xr) < abs(*si) * eps) {
		return 0;
	    }
/* L20: */
	}
    } else if (*x <= 32.) {
	m = (integer) (*x * .82f + 47.2f);
	xa1 = 0.;
	xa0 = 1e-100;
	for (k = m; k >= 1; --k) {
	    xa = k * 4. * xa0 / *x - xa1;
	    bj[k - 1] = xa;
	    xa1 = xa0;
/* L25: */
	    xa0 = xa;
	}
	xs = bj[0];
	i__1 = m;
	for (k = 3; k <= i__1; k += 2) {
/* L30: */
	    xs += bj[k - 1] * 2.;
	}
	bj[0] /= xs;
	i__1 = m;
	for (k = 2; k <= i__1; ++k) {
/* L35: */
	    bj[k - 1] /= xs;
	}
	xr = 1.;
	xg1 = bj[0];
	i__1 = m;
	for (k = 2; k <= i__1; ++k) {
/* Computing 2nd power */
	    r__1 = k * 2.f - 3.f;
/* Computing 2nd power */
	    r__2 = k * 2.f - 1.f;
	    xr = xr * .25 * (r__1 * r__1) / ((k - 1.f) * (r__2 * r__2)) * *x;
/* L40: */
	    xg1 += bj[k - 1] * xr;
	}
	xr = 1.;
	xg2 = bj[0];
	i__1 = m;
	for (k = 2; k <= i__1; ++k) {
/* Computing 2nd power */
	    r__1 = k * 2.f - 5.f;
/* Computing 2nd power */
	    r__2 = k * 2.f - 3.f;
	    xr = xr * .25 * (r__1 * r__1) / ((k - 1.f) * (r__2 * r__2)) * *x;
/* L45: */
	    xg2 += bj[k - 1] * xr;
	}
	xcs = cos(*x / 2.);
	xss = sin(*x / 2.);
	*ci = el + log(*x) - *x * xss * xg1 + xcs * 2 * xg2 - xcs * 2 * xcs;
	*si = *x * xcs * xg1 + xss * 2 * xg2 - sin(*x);
    } else {
	xr = 1.;
	xf = 1.;
	for (k = 1; k <= 9; ++k) {
	    xr = xr * -2. * k * ((k << 1) - 1) / x2;
/* L50: */
	    xf += xr;
	}
	xr = 1. / *x;
	xg = xr;
	for (k = 1; k <= 8; ++k) {
	    xr = xr * -2. * ((k << 1) + 1) * k / x2;
/* L55: */
	    xg += xr;
	}
	*ci = xf * sin(*x) / *x - xg * cos(*x) / *x;
	*si = p2 - xf * cos(*x) / *x - xg * sin(*x) / *x;
    }
    return 0;
} /* cisia_ */

/*       ********************************** */
/* Subroutine */ int itsl0_(doublereal *x, doublereal *tl0)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal a[18];
    static integer k;
    static doublereal r__, s, a0, a1, s0, af, el, rd, pi, ti;


/*       =========================================================== */
/*       Purpose: Evaluate the integral of modified Struve function */
/*                L0(t) with respect to t from 0 to x */
/*       Input :  x   --- Upper limit  ( x ≥ 0 ) */
/*       Output:  TL0 --- Integration of L0(t) from 0 to x */
/*       =========================================================== */

    pi = 3.141592653589793;
    r__ = 1.;
    if (*x <= 20.f) {
	s = .5;
	for (k = 1; k <= 100; ++k) {
	    rd = 1.;
	    if (k == 1) {
		rd = .5;
	    }
/* Computing 2nd power */
	    d__1 = *x / (k * 2. + 1.);
	    r__ = r__ * rd * k / (k + 1.) * (d__1 * d__1);
	    s += r__;
	    if ((d__1 = r__ / s, abs(d__1)) < 1e-12) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*tl0 = 2. / pi * *x * *x * s;
    } else {
	s = 1.;
	for (k = 1; k <= 10; ++k) {
/* Computing 2nd power */
	    d__1 = (k * 2. + 1.) / *x;
	    r__ = r__ * k / (k + 1.) * (d__1 * d__1);
	    s += r__;
	    if ((d__1 = r__ / s, abs(d__1)) < 1e-12) {
		goto L25;
	    }
/* L20: */
	}
L25:
	el = .57721566490153;
	s0 = -s / (pi * *x * *x) + 2. / pi * (log(*x * 2.) + el);
	a0 = 1.;
	a1 = .625;
	a[0] = a1;
	for (k = 1; k <= 10; ++k) {
/* Computing 2nd power */
	    d__1 = k + .5;
	    af = ((k + .5) * 1.5 * (k + .83333333333333337) * a1 - d__1 * 
		    d__1 * .5 * (k - .5) * a0) / (k + 1.);
	    a[k] = af;
	    a0 = a1;
/* L30: */
	    a1 = af;
	}
	ti = 1.;
	r__ = 1.;
	for (k = 1; k <= 11; ++k) {
	    r__ /= *x;
/* L35: */
	    ti += a[k - 1] * r__;
	}
	*tl0 = ti / sqrt(pi * 2 * *x) * exp(*x) + s0;
    }
    return 0;
} /* itsl0_ */

/*       ********************************** */
/* Subroutine */ int clqn_(integer *n, doublereal *x, doublereal *y, 
	doublecomplex *cqn, doublecomplex *cqd)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_log(
	    doublecomplex *, doublecomplex *);
    double log(doublereal);

    /* Local variables */
    static integer k;
    static doublecomplex z__;
    static integer km, ls;
    static doublecomplex cq0, cq1, cqf0, cqf1, cqf2;


/*       ================================================== */
/*       Purpose: Compute the Legendre functions Qn(z) and */
/*                their derivatives Qn'(z) for a complex */
/*                argument */
/*       Input :  x --- Real part of z */
/*                y --- Imaginary part of z */
/*                n --- Degree of Qn(z), n = 0,1,2,... */
/*       Output:  CQN(n) --- Qn(z) */
/*                CQD(n) --- Qn'(z) */
/*       ================================================== */

    z__1.r = *x, z__1.i = *y;
    z__.r = z__1.r, z__.i = z__1.i;
    if (z__.r == 1. && z__.i == 0.) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cqn[i__2].r = 1e300, cqn[i__2].i = 0.;
/* L10: */
	    i__2 = k;
	    cqd[i__2].r = 1e300, cqd[i__2].i = 0.;
	}
	return 0;
    }
    ls = 1;
    if (z_abs(&z__) > 1.) {
	ls = -1;
    }
    z__5.r = z__.r + 1., z__5.i = z__.i;
    d__1 = (doublereal) ls;
    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
    z__6.r = 1. - z__.r, z__6.i = -z__.i;
    z_div(&z__3, &z__4, &z__6);
    z_log(&z__2, &z__3);
    z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
    cq0.r = z__1.r, cq0.i = z__1.i;
    z__2.r = z__.r * cq0.r - z__.i * cq0.i, z__2.i = z__.r * cq0.i + z__.i * 
	    cq0.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    cq1.r = z__1.r, cq1.i = z__1.i;
    cqn[0].r = cq0.r, cqn[0].i = cq0.i;
    cqn[1].r = cq1.r, cqn[1].i = cq1.i;
    if (z_abs(&z__) < 1.0001) {
	cqf0.r = cq0.r, cqf0.i = cq0.i;
	cqf1.r = cq1.r, cqf1.i = cq1.i;
	i__2 = *n;
	for (k = 2; k <= i__2; ++k) {
	    d__1 = k * 2. - 1.;
	    z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	    z__3.r = z__4.r * cqf1.r - z__4.i * cqf1.i, z__3.i = z__4.r * 
		    cqf1.i + z__4.i * cqf1.r;
	    d__2 = k - 1.;
	    z__5.r = d__2 * cqf0.r, z__5.i = d__2 * cqf0.i;
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	    d__3 = (doublereal) k;
	    z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
	    cqf2.r = z__1.r, cqf2.i = z__1.i;
	    i__1 = k;
	    cqn[i__1].r = cqf2.r, cqn[i__1].i = cqf2.i;
	    cqf0.r = cqf1.r, cqf0.i = cqf1.i;
/* L15: */
	    cqf1.r = cqf2.r, cqf1.i = cqf2.i;
	}
    } else {
	if (z_abs(&z__) > 1.1) {
	    km = *n + 40;
	} else {
	    z__1.r = z__.r - 1.f, z__1.i = z__.i;
	    km = (*n + 40) * (integer) (-1.f - log(z_abs(&z__1)) * 1.8f);
	}
	cqf2.r = 0., cqf2.i = 0.;
	cqf1.r = 1., cqf1.i = 0.;
	for (k = km; k >= 0; --k) {
	    d__1 = (k << 1) + 3.;
	    z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	    z__3.r = z__4.r * cqf1.r - z__4.i * cqf1.i, z__3.i = z__4.r * 
		    cqf1.i + z__4.i * cqf1.r;
	    d__2 = k + 2.;
	    z__5.r = d__2 * cqf2.r, z__5.i = d__2 * cqf2.i;
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	    d__3 = k + 1.;
	    z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
	    cqf0.r = z__1.r, cqf0.i = z__1.i;
	    if (k <= *n) {
		i__2 = k;
		cqn[i__2].r = cqf0.r, cqn[i__2].i = cqf0.i;
	    }
	    cqf2.r = cqf1.r, cqf2.i = cqf1.i;
/* L20: */
	    cqf1.r = cqf0.r, cqf1.i = cqf0.i;
	}
	i__2 = *n;
	for (k = 0; k <= i__2; ++k) {
/* L25: */
	    i__1 = k;
	    i__3 = k;
	    z__2.r = cqn[i__3].r * cq0.r - cqn[i__3].i * cq0.i, z__2.i = cqn[
		    i__3].r * cq0.i + cqn[i__3].i * cq0.r;
	    z_div(&z__1, &z__2, &cqf0);
	    cqn[i__1].r = z__1.r, cqn[i__1].i = z__1.i;
	}
    }
    z__3.r = z__.r * cqn[0].r - z__.i * cqn[0].i, z__3.i = z__.r * cqn[0].i + 
	    z__.i * cqn[0].r;
    z__2.r = cqn[1].r - z__3.r, z__2.i = cqn[1].i - z__3.i;
    z__5.r = z__.r * z__.r - z__.i * z__.i, z__5.i = z__.r * z__.i + z__.i * 
	    z__.r;
    z__4.r = z__5.r - 1., z__4.i = z__5.i;
    z_div(&z__1, &z__2, &z__4);
    cqd[0].r = z__1.r, cqd[0].i = z__1.i;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* L30: */
	i__3 = k;
	d__1 = (doublereal) k;
	z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	i__2 = k;
	z__3.r = z__4.r * cqn[i__2].r - z__4.i * cqn[i__2].i, z__3.i = z__4.r 
		* cqn[i__2].i + z__4.i * cqn[i__2].r;
	i__4 = k - 1;
	d__2 = (doublereal) k;
	z__5.r = d__2 * cqn[i__4].r, z__5.i = d__2 * cqn[i__4].i;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	z__7.r = z__.r * z__.r - z__.i * z__.i, z__7.i = z__.r * z__.i + 
		z__.i * z__.r;
	z__6.r = z__7.r - 1., z__6.i = z__7.i;
	z_div(&z__1, &z__2, &z__6);
	cqd[i__3].r = z__1.r, cqd[i__3].i = z__1.i;
    }
    return 0;
} /* clqn_ */

/*       ********************************** */
/* Subroutine */ int airyzo_(integer *nt, integer *kf, doublereal *xa, 
	doublereal *xb, doublereal *xc, doublereal *xd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal u, x, u1, ad, bd, ai, bi, pi, rt, rt0, err;
    extern /* Subroutine */ int airyb_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);


/*       ======================================================== */
/*       Purpose: Compute the first NT zeros of Airy functions */
/*                Ai(x) and Ai'(x), a and a', and the associated */
/*                values of Ai(a') and Ai'(a); and the first NT */
/*                zeros of Airy functions Bi(x) and Bi'(x), b and */
/*                b', and the associated values of Bi(b') and */
/*                Bi'(b) */
/*       Input :  NT    --- Total number of zeros */
/*                KF    --- Function code */
/*                          KF=1 for Ai(x) and Ai'(x) */
/*                          KF=2 for Bi(x) and Bi'(x) */
/*       Output:  XA(m) --- a, the m-th zero of Ai(x) or */
/*                          b, the m-th zero of Bi(x) */
/*                XB(m) --- a', the m-th zero of Ai'(x) or */
/*                          b', the m-th zero of Bi'(x) */
/*                XC(m) --- Ai(a') or Bi(b') */
/*                XD(m) --- Ai'(a) or Bi'(b) */
/*                          ( m --- Serial number of zeros ) */
/*       Routine called: AIRYB for computing Airy functions and */
/*                       their derivatives */
/*       ======================================================= */

    /* Parameter adjustments */
    --xd;
    --xc;
    --xb;
    --xa;

    /* Function Body */
    pi = 3.141592653589793;
    rt = 0.;
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rt0 = 0.;
	if (*kf == 1) {
	    u = pi * 3. * (i__ * 4. - 1) / 8.;
	    u1 = 1 / (u * u);
	} else if (*kf == 2) {
	    if (i__ == 1) {
		rt0 = -1.17371;
	    } else {
		u = pi * 3. * (i__ * 4. - 3.) / 8.;
		u1 = 1 / (u * u);
	    }
	}
	if (rt0 == 0.) {
/*             DLMF 9.9.18 */
	    d__1 = u * u;
	    rt0 = -pow_dd(&d__1, &c_b221) * (u1 * (u1 * (u1 * (u1 * 
		    -15.509155201673648 + .92984423225308643) - 
		    .1388888888888889) + .10416666666666667) + 1.);
	}
L10:
	x = rt0;
	airyb_(&x, &ai, &bi, &ad, &bd);
	if (*kf == 1) {
	    rt = rt0 - ai / ad;
	}
	if (*kf == 2) {
	    rt = rt0 - bi / bd;
	}
	err = (d__1 = (rt - rt0) / rt, abs(d__1));
	if (err > 1e-12) {
	    rt0 = rt;
	    goto L10;
	} else {
	    xa[i__] = rt;
	    if (err > 1e-14) {
		airyb_(&rt, &ai, &bi, &ad, &bd);
	    }
	    if (*kf == 1) {
		xd[i__] = ad;
	    }
	    if (*kf == 2) {
		xd[i__] = bd;
	    }
	}
/* L15: */
    }
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rt0 = 0.;
	if (*kf == 1) {
	    if (i__ == 1) {
		rt0 = -1.01879;
	    } else {
		u = pi * 3. * (i__ * 4. - 3.) / 8.;
		u1 = 1 / (u * u);
	    }
	} else if (*kf == 2) {
	    if (i__ == 1) {
		rt0 = -2.29444;
	    } else {
		u = pi * 3. * (i__ * 4. - 1.) / 8.;
		u1 = 1 / (u * u);
	    }
	}
	if (rt0 == 0.) {
/*             DLMF 9.9.19 */
	    d__1 = u * u;
	    rt0 = -pow_dd(&d__1, &c_b221) * (u1 * (u1 * (u1 * (u1 * 
		    15.016855549125514 - .87395351080246919) + 
		    .12152777777777778) - .14583333333333334) + 1.);
	}
L20:
	x = rt0;
	airyb_(&x, &ai, &bi, &ad, &bd);
	if (*kf == 1) {
	    rt = rt0 - ad / (ai * x);
	}
	if (*kf == 2) {
	    rt = rt0 - bd / (bi * x);
	}
	err = (d__1 = (rt - rt0) / rt, abs(d__1));
	if (err > 1e-12) {
	    rt0 = rt;
	    goto L20;
	} else {
	    xb[i__] = rt;
	    if (err > 1e-14) {
		airyb_(&rt, &ai, &bi, &ad, &bd);
	    }
	    if (*kf == 1) {
		xc[i__] = ai;
	    }
	    if (*kf == 2) {
		xc[i__] = bi;
	    }
	}
/* L25: */
    }
    return 0;
} /* airyzo_ */

/*       ********************************** */
/* Subroutine */ int error_(doublereal *x, doublereal *err)
{
    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, c0, x2, er, pi, eps;


/*       ========================================= */
/*       Purpose: Compute error function erf(x) */
/*       Input:   x   --- Argument of erf(x) */
/*       Output:  ERR --- erf(x) */
/*       ========================================= */

    eps = 1e-15;
    pi = 3.141592653589793;
    x2 = *x * *x;
    if (abs(*x) < 3.5) {
	er = 1.;
	r__ = 1.;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * x2 / (k + .5);
	    er += r__;
	    if (abs(r__) <= abs(er) * eps) {
		goto L15;
	    }
/* L10: */
	}
L15:
	c0 = 2. / sqrt(pi) * *x * exp(-x2);
	*err = c0 * er;
    } else {
	er = 1.;
	r__ = 1.;
	for (k = 1; k <= 12; ++k) {
	    r__ = -r__ * (k - .5) / x2;
/* L20: */
	    er += r__;
	}
	c0 = exp(-x2) / (abs(*x) * sqrt(pi));
	*err = 1. - c0 * er;
	if (*x < 0.f) {
	    *err = -(*err);
	}
    }
    return 0;
} /* error_ */

/*       ********************************** */
/* Subroutine */ int cerror_(doublecomplex *z__, doublecomplex *cer)
{
    /* System generated locals */
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_exp(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal a0;
    static doublecomplex c0, z1, cl, cr, cs;
    static doublereal pi;


/*       ==================================================== */
/*       Purpose: Compute error function erf(z) for a complex */
/*                argument (z=x+iy) */
/*       Input :  z   --- Complex argument */
/*       Output:  CER --- erf(z) */
/*       ==================================================== */

    a0 = z_abs(z__);
    z__3.r = -z__->r, z__3.i = -z__->i;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z_exp(&z__1, &z__2);
    c0.r = z__1.r, c0.i = z__1.i;
    pi = 3.141592653589793;
    z1.r = z__->r, z1.i = z__->i;
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }

/*       Cutoff radius R = 4.36; determined by balancing rounding error */
/*       and asymptotic expansion error, see below. */

/*       The resulting maximum global accuracy expected is around 1e-8 */

    if (a0 <= 4.36) {

/*          Rounding error in the Taylor expansion is roughly */

/*          ~ R*R * EPSILON * R**(2 R**2) / (2 R**2 Gamma(R**2 + 1/2)) */

	cs.r = z1.r, cs.i = z1.i;
	cr.r = z1.r, cr.i = z1.i;
	for (k = 1; k <= 120; ++k) {
	    z__3.r = cr.r * z1.r - cr.i * z1.i, z__3.i = cr.r * z1.i + cr.i * 
		    z1.r;
	    z__2.r = z__3.r * z1.r - z__3.i * z1.i, z__2.i = z__3.r * z1.i + 
		    z__3.i * z1.r;
	    d__1 = k + .5;
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cs.r + cr.r, z__1.i = cs.i + cr.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	    z_div(&z__1, &cr, &cs);
	    if (z_abs(&z__1) < 1e-15) {
		goto L15;
	    }
/* L10: */
	}
L15:
	z__3.r = c0.r * 2., z__3.i = c0.i * 2.;
	z__2.r = z__3.r * cs.r - z__3.i * cs.i, z__2.i = z__3.r * cs.i + 
		z__3.i * cs.r;
	d__1 = sqrt(pi);
	z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	cer->r = z__1.r, cer->i = z__1.i;
    } else {
	z_div(&z__1, &c_b6, &z1);
	cl.r = z__1.r, cl.i = z__1.i;
	cr.r = cl.r, cr.i = cl.i;

/*          Asymptotic series; maximum K must be at most ~ R^2. */

/*          The maximum accuracy obtainable from this expansion is roughly */

/*          ~ Gamma(2R**2 + 2) / ( */
/*                   (2 R**2)**(R**2 + 1/2) Gamma(R**2 + 3/2) 2**(R**2 + 1/2)) */

	for (k = 1; k <= 20; ++k) {
	    z__3.r = -cr.r, z__3.i = -cr.i;
	    d__1 = k - .5;
	    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	    z__4.r = z1.r * z1.r - z1.i * z1.i, z__4.i = z1.r * z1.i + z1.i * 
		    z1.r;
	    z_div(&z__1, &z__2, &z__4);
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cl.r + cr.r, z__1.i = cl.i + cr.i;
	    cl.r = z__1.r, cl.i = z__1.i;
	    z_div(&z__1, &cr, &cl);
	    if (z_abs(&z__1) < 1e-15) {
		goto L25;
	    }
/* L20: */
	}
L25:
	z__3.r = c0.r * cl.r - c0.i * cl.i, z__3.i = c0.r * cl.i + c0.i * 
		cl.r;
	d__1 = sqrt(pi);
	z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
	z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
	cer->r = z__1.r, cer->i = z__1.i;
    }
    if (z__->r < 0.f) {
	z__1.r = -cer->r, z__1.i = -cer->i;
	cer->r = z__1.r, cer->i = z__1.i;
    }
    return 0;
} /* cerror_ */

/*       ********************************** */
/* Subroutine */ int eulerb_(integer *n, doublereal *en)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer k, m;
    static doublereal s, r1, r2, hpi;
    static integer isgn;


/*       ====================================== */
/*       Purpose: Compute Euler number En */
/*       Input :  n --- Serial number */
/*       Output:  EN(n) --- En */
/*       ====================================== */

    hpi = .63661977236758138;
    en[0] = 1.;
    en[2] = -1.;
/* Computing 3rd power */
    d__1 = hpi;
    r1 = d__1 * (d__1 * d__1) * -4.;
    i__1 = *n;
    for (m = 4; m <= i__1; m += 2) {
	r1 = -r1 * (m - 1) * m * hpi * hpi;
	r2 = 1.;
	isgn = 1;
	for (k = 3; k <= 1000; k += 2) {
	    isgn = -isgn;
	    d__1 = 1. / k;
	    i__2 = m + 1;
	    s = pow_di(&d__1, &i__2);
	    r2 += isgn * s;
	    if (s < 1e-15) {
		goto L20;
	    }
/* L10: */
	}
L20:
	en[m] = r1 * r2;
    }
    return 0;
} /* eulerb_ */

/*       ********************************** */
/* Subroutine */ int cva1_(integer *kd, integer *m, doublereal *q, doublereal 
	*cv)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal d__[500], e[500], f[500], g[200], h__[200];
    static integer i__, j, k;
    static doublereal s, t;
    static integer k1;
    static doublereal t1, x1;
    static integer ic;
    static doublereal xa;
    static integer nm;
    static doublereal xb;
    static integer nm1, icm;
    static doublereal eps;


/*       ============================================================ */
/*       Purpose: Compute a sequence of characteristic values of */
/*                Mathieu functions */
/*       Input :  M  --- Maximum order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions */
/*                KD --- Case code */
/*                       KD=1 for cem(x,q)  ( m = 0,2,4,… ) */
/*                       KD=2 for cem(x,q)  ( m = 1,3,5,… ) */
/*                       KD=3 for sem(x,q)  ( m = 1,3,5,… ) */
/*                       KD=4 for sem(x,q)  ( m = 2,4,6,… ) */
/*       Output:  CV(I) --- Characteristic values; I = 1,2,3,... */
/*                For KD=1, CV(1), CV(2), CV(3),..., correspond to */
/*                the characteristic values of cem for m = 0,2,4,... */
/*                For KD=2, CV(1), CV(2), CV(3),..., correspond to */
/*                the characteristic values of cem for m = 1,3,5,... */
/*                For KD=3, CV(1), CV(2), CV(3),..., correspond to */
/*                the characteristic values of sem for m = 1,3,5,... */
/*                For KD=4, CV(1), CV(2), CV(3),..., correspond to */
/*                the characteristic values of sem for m = 0,2,4,... */
/*       ============================================================ */

    /* Parameter adjustments */
    --cv;

    /* Function Body */
    eps = 1e-14;
    icm = *m / 2 + 1;
    if (*kd == 4) {
	icm = *m / 2;
    }
    if (*q == 0.) {
	if (*kd == 1) {
	    i__1 = icm;
	    for (ic = 1; ic <= i__1; ++ic) {
/* L10: */
/* Computing 2nd power */
		d__1 = ic - 1.;
		cv[ic] = d__1 * d__1 * 4.;
	    }
	} else if (*kd != 4) {
	    i__1 = icm;
	    for (ic = 1; ic <= i__1; ++ic) {
/* L15: */
/* Computing 2nd power */
		d__1 = ic * 2. - 1.;
		cv[ic] = d__1 * d__1;
	    }
	} else {
	    i__1 = icm;
	    for (ic = 1; ic <= i__1; ++ic) {
/* L20: */
		cv[ic] = ic * 4. * ic;
	    }
	}
    } else {
	nm = (integer) (*m * 1.5f + 10 + *q * .5f);
	e[0] = 0.;
	f[0] = 0.;
	if (*kd == 1) {
	    d__[0] = 0.;
	    i__1 = nm;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = i__ - 1.;
		d__[i__ - 1] = d__1 * d__1 * 4.;
		e[i__ - 1] = *q;
/* L25: */
		f[i__ - 1] = *q * *q;
	    }
	    e[1] = sqrt(2.) * *q;
	    f[1] = *q * 2. * *q;
	} else if (*kd != 4) {
	    d__[0] = pow_ii(&c_n1, kd) * *q + 1.;
	    i__1 = nm;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = i__ * 2. - 1.;
		d__[i__ - 1] = d__1 * d__1;
		e[i__ - 1] = *q;
/* L30: */
		f[i__ - 1] = *q * *q;
	    }
	} else {
	    d__[0] = 4.;
	    i__1 = nm;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		d__[i__ - 1] = i__ * 4. * i__;
		e[i__ - 1] = *q;
/* L35: */
		f[i__ - 1] = *q * *q;
	    }
	}
	xa = d__[nm - 1] + (d__1 = e[nm - 1], abs(d__1));
	xb = d__[nm - 1] - (d__1 = e[nm - 1], abs(d__1));
	nm1 = nm - 1;
	i__1 = nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = (d__1 = e[i__ - 1], abs(d__1)) + (d__2 = e[i__], abs(d__2));
	    t1 = d__[i__ - 1] + t;
	    if (xa < t1) {
		xa = t1;
	    }
	    t1 = d__[i__ - 1] - t;
	    if (t1 < xb) {
		xb = t1;
	    }
/* L40: */
	}
	i__1 = icm;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    g[i__ - 1] = xa;
/* L45: */
	    h__[i__ - 1] = xb;
	}
	i__1 = icm;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = icm;
	    for (k1 = k; k1 <= i__2; ++k1) {
		if (g[k1 - 1] < g[k - 1]) {
		    g[k - 1] = g[k1 - 1];
		    goto L55;
		}
/* L50: */
	    }
L55:
	    if (k != 1 && h__[k - 1] < h__[k - 2]) {
		h__[k - 1] = h__[k - 2];
	    }
L60:
	    x1 = (g[k - 1] + h__[k - 1]) / 2.;
	    cv[k] = x1;
	    if ((d__1 = (g[k - 1] - h__[k - 1]) / x1, abs(d__1)) < eps) {
		goto L70;
	    }
	    j = 0;
	    s = 1.;
	    i__2 = nm;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (s == 0.) {
		    s += 1e-30;
		}
		t = f[i__ - 1] / s;
		s = d__[i__ - 1] - t - x1;
		if (s < 0.f) {
		    ++j;
		}
/* L65: */
	    }
	    if (j < k) {
		h__[k - 1] = x1;
	    } else {
		g[k - 1] = x1;
		if (j >= icm) {
		    g[icm - 1] = x1;
		} else {
		    if (h__[j] < x1) {
			h__[j] = x1;
		    }
		    if (x1 < g[j - 1]) {
			g[j - 1] = x1;
		    }
		}
	    }
	    goto L60;
L70:
	    cv[k] = x1;
/* L75: */
	}
    }
    return 0;
} /* cva1_ */

/*       ********************************** */
/* Subroutine */ int ittikb_(doublereal *x, doublereal *tti, doublereal *ttk)
{
    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), log(doublereal);

    /* Local variables */
    static doublereal t, e0, t1, x1, el, pi;


/*       ========================================================= */
/*       Purpose: Integrate [I0(t)-1]/t with respect to t from 0 */
/*                to x, and K0(t)/t with respect to t from x to ∞ */
/*       Input :  x   --- Variable in the limits  ( x ≥ 0 ) */
/*       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x */
/*                TTK --- Integration of K0(t)/t from x to ∞ */
/*       ========================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    if (*x == 0.) {
	*tti = 0.;
    } else if (*x <= 5.) {
	x1 = *x / 5.;
	t = x1 * x1;
	*tti = (((((((t * 1.263e-4 + 9.6442e-4) * t + .00968217) * t + 
		.06615507) * t + .33116853) * t + 1.13027241) * t + 
		2.44140746) * t + 3.12499991) * t;
    } else {
	t = 5. / *x;
	*tti = (((((((((t * 2.1945464 - 3.5195009) * t - 11.9094395) * t + 
		40.394734) * t - 48.0524115) * t + 28.1221478) * t - 
		8.6556013) * t + 1.4780044) * t - .0493843) * t + .1332055) * 
		t + .3989314;
	*tti = *tti * exp(*x) / (sqrt(*x) * *x);
    }
    if (*x == 0.) {
	*ttk = 1e300;
    } else if (*x <= 2.) {
	t1 = *x / 2.;
	t = t1 * t1;
	*ttk = (((((t * 7.7e-7 + 1.544e-5) * t + 4.8077e-4) * t + .00925821) *
		 t + .10937537) * t + .74999993) * t;
	e0 = el + log(*x / 2.);
	*ttk = pi * pi / 24. + e0 * (e0 * .5 + *tti) - *ttk;
    } else if (*x <= 4.) {
	t = 2. / *x;
	*ttk = (((t * .06084 - .280367) * t + .590944) * t - .850013) * t + 
		1.234684;
	*ttk = *ttk * exp(-(*x)) / (sqrt(*x) * *x);
    } else {
	t = 4. / *x;
	*ttk = (((((t * .02724 - .1110396) * t + .2060126) * t - .2621446) * 
		t + .3219184) * t - .5091339) * t + 1.2533141;
	*ttk = *ttk * exp(-(*x)) / (sqrt(*x) * *x);
    }
    return 0;
} /* ittikb_ */

/*       ********************************** */
/* Subroutine */ int lqnb_(integer *n, doublereal *x, doublereal *qn, 
	doublereal *qd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer j, k, l;
    static doublereal q0, q1, x2, qf;
    static integer nl;
    static doublereal qr, qc1, qc2, qf0, qf1, qf2, eps;


/*       ==================================================== */
/*       Purpose: Compute Legendre functions Qn(x) & Qn'(x) */
/*       Input :  x  --- Argument of Qn(x) */
/*                n  --- Degree of Qn(x)  ( n = 0,1,2,…) */
/*       Output:  QN(n) --- Qn(x) */
/*                QD(n) --- Qn'(x) */
/*       ==================================================== */

    eps = 1e-14;
    if (abs(*x) == 1.) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    qn[k] = 1e300;
/* L10: */
	    qd[k] = 1e300;
	}
	return 0;
    }
    if (*x <= 1.021) {
	x2 = (d__1 = (*x + 1.) / (1. - *x), abs(d__1));
	q0 = log(x2) * .5;
	q1 = *x * q0 - 1.;
	qn[0] = q0;
	qn[1] = q1;
	qd[0] = 1. / (1. - *x * *x);
	qd[1] = qn[0] + *x * qd[0];
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    qf = ((k * 2. - 1.) * *x * q1 - (k - 1.) * q0) / k;
	    qn[k] = qf;
	    qd[k] = (qn[k - 1] - *x * qf) * k / (1. - *x * *x);
	    q0 = q1;
/* L15: */
	    q1 = qf;
	}
    } else {
	qc1 = 0.;
	qc2 = 1. / *x;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    qc2 = qc2 * j / ((j * 2.f + 1.) * *x);
	    if (j == *n - 1) {
		qc1 = qc2;
	    }
/* L20: */
	}
	for (l = 0; l <= 1; ++l) {
	    nl = *n + l;
	    qf = 1.;
	    qr = 1.;
	    for (k = 1; k <= 500; ++k) {
		qr = qr * (nl * .5 + k - 1.) * ((nl - 1) * .5 + k) / ((nl + k 
			- .5) * k * *x * *x);
		qf += qr;
		if ((d__1 = qr / qf, abs(d__1)) < eps) {
		    goto L30;
		}
/* L25: */
	    }
L30:
	    if (l == 0) {
		qn[*n - 1] = qf * qc1;
	    } else {
		qn[*n] = qf * qc2;
	    }
/* L35: */
	}
	qf2 = qn[*n];
	qf1 = qn[*n - 1];
	for (k = *n; k >= 2; --k) {
	    qf0 = (((k << 1) - 1.) * *x * qf1 - k * qf2) / (k - 1.);
	    qn[k - 2] = qf0;
	    qf2 = qf1;
/* L40: */
	    qf1 = qf0;
	}
	qd[0] = 1. / (1. - *x * *x);
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
/* L45: */
	    qd[k] = k * (qn[k - 1] - *x * qn[k]) / (1. - *x * *x);
	}
    }
    return 0;
} /* lqnb_ */

/*       ********************************** */
/* Subroutine */ int cjk_(integer *km, doublereal *a)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal f, g;
    static integer j, k;
    static doublereal f0, g0;
    static integer l1, l2, l3, l4;


/*       ======================================================== */
/*       Purpose: Compute the expansion coefficients for the */
/*                asymptotic expansion of Bessel functions */
/*                with large orders */
/*       Input :  Km   --- Maximum k */
/*       Output:  A(L) --- Cj(k) where j and k are related to L */
/*                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km */
/*       ======================================================== */

    /* Parameter adjustments */
    --a;

    /* Function Body */
    a[1] = 1.;
    f0 = 1.;
    g0 = 1.;
    i__1 = *km - 1;
    for (k = 0; k <= i__1; ++k) {
	l1 = (k + 1) * (k + 2) / 2 + 1;
	l2 = (k + 1) * (k + 2) / 2 + k + 2;
	f = (k * .5 + .125 / (k + 1)) * f0;
	g = -(k * 1.5 + .625 / ((k + 1.) * 3.f)) * g0;
	a[l1] = f;
	a[l2] = g;
	f0 = f;
/* L10: */
	g0 = g;
    }
    i__1 = *km - 1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    l3 = k * (k + 1) / 2 + j + 1;
	    l4 = (k + 1) * (k + 2) / 2 + j + 1;
	    a[l4] = (j + k * .5 + .125 / (j * 2.f + k + 1.f)) * a[l3] - (j + 
		    k * .5 - 1.f + .625 / (j * 2.f + k + 1.f)) * a[l3 - 1];
/* L15: */
	}
    }
    return 0;
} /* cjk_ */

/*       ********************************** */
/* Subroutine */ int ittika_(doublereal *x, doublereal *tti, doublereal *ttk)
{
    /* Initialized data */

    static doublereal c__[8] = { 1.625,4.1328125,14.5380859375,65.53353881835,
	    360.66157150269,2344.8727161884,17588.273098916,149506.39538279 };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), log(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, b1, e0, r2, el, rc, pi, rs;


/*       ========================================================= */
/*       Purpose: Integrate [I0(t)-1]/t with respect to t from 0 */
/*                to x, and K0(t)/t with respect to t from x to ∞ */
/*       Input :  x   --- Variable in the limits  ( x ≥ 0 ) */
/*       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x */
/*                TTK --- Integration of K0(t)/t from x to ∞ */
/*       ========================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    if (*x == 0.) {
	*tti = 0.;
	*ttk = 1e300;
	return 0;
    }
    if (*x < 40.) {
	*tti = 1.;
	r__ = 1.;
	for (k = 2; k <= 50; ++k) {
	    r__ = r__ * .25 * (k - 1.) / (k * k * k) * *x * *x;
	    *tti += r__;
	    if ((d__1 = r__ / *tti, abs(d__1)) < 1e-12) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*tti = *tti * .125 * *x * *x;
    } else {
	*tti = 1.;
	r__ = 1.;
	for (k = 1; k <= 8; ++k) {
	    r__ /= *x;
/* L20: */
	    *tti += c__[k - 1] * r__;
	}
	rc = *x * sqrt(pi * 2. * *x);
	*tti = *tti * exp(*x) / rc;
    }
    if (*x <= 12.) {
	e0 = (log(*x / 2.) * .5 + el) * log(*x / 2.) + pi * pi / 24. + el * 
		.5 * el;
	b1 = 1.5 - (el + log(*x / 2.));
	rs = 1.;
	r__ = 1.;
	for (k = 2; k <= 50; ++k) {
	    r__ = r__ * .25 * (k - 1.) / (k * k * k) * *x * *x;
	    rs += 1. / k;
	    r2 = r__ * (rs + 1. / (k * 2.) - (el + log(*x / 2.)));
	    b1 += r2;
	    if ((d__1 = r2 / b1, abs(d__1)) < 1e-12) {
		goto L30;
	    }
/* L25: */
	}
L30:
	*ttk = e0 - *x * .125 * *x * b1;
    } else {
	*ttk = 1.;
	r__ = 1.;
	for (k = 1; k <= 8; ++k) {
	    r__ = -r__ / *x;
/* L35: */
	    *ttk += c__[k - 1] * r__;
	}
	rc = *x * sqrt(2. / pi * *x);
	*ttk = *ttk * exp(-(*x)) / rc;
    }
    return 0;
} /* ittika_ */

/*       ********************************** */
/* Subroutine */ int lamv_(doublereal *v, doublereal *x, doublereal *vm, 
	doublereal *vl, doublereal *dl)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), cos(
	    doublereal), sin(doublereal);

    /* Local variables */
    static doublereal f;
    static integer i__, j, k, m, n;
    static doublereal r__, a0, f0, f1, f2;
    static integer k0;
    static doublereal r0, v0, x2, ga, bk, ck, rc, cs, pi, sk, uk, vk, rp, rq, 
	    xk, px, qx, vv, rp2, fac;
    extern /* Subroutine */ int gam0_(doublereal *, doublereal *);
    static doublereal bjv0, bjv1;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    static doublereal elsecs;


/*       ========================================================= */
/*       Purpose: Compute lambda function with arbitrary order v, */
/*                and their derivative */
/*       Input :  x --- Argument of lambda function */
/*                v --- Order of lambda function */
/*       Output:  VL(n) --- Lambda function of order n+v0 */
/*                DL(n) --- Derivative of lambda function */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*            (2) GAM0 for computing gamma function (|x| ≤ 1) */
/*       ========================================================= */

    pi = 3.141592653589793;
    rp2 = .63661977236758;
    *x = abs(*x);
    x2 = *x * *x;
    n = (integer) (*v);
    v0 = *v - n;
    *vm = *v;
    if (*x <= 12.) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    vk = v0 + k;
	    bk = 1.;
	    r__ = 1.;
	    for (i__ = 1; i__ <= 50; ++i__) {
		r__ = r__ * -.25 * x2 / (i__ * (i__ + vk));
		bk += r__;
		if (abs(r__) < abs(bk) * 1e-15) {
		    goto L15;
		}
/* L10: */
	    }
L15:
	    vl[k] = bk;
	    uk = 1.;
	    r__ = 1.;
	    for (i__ = 1; i__ <= 50; ++i__) {
		r__ = r__ * -.25 * x2 / (i__ * (i__ + vk + 1.));
		uk += r__;
		if (abs(r__) < abs(uk) * 1e-15) {
		    goto L25;
		}
/* L20: */
	    }
L25:
	    dl[k] = *x * -.5 / (vk + 1.) * uk;
	}
	return 0;
    }
    k0 = 11;
    if (*x >= 35.) {
	k0 = 10;
    }
    if (*x >= 50.) {
	k0 = 8;
    }
    bjv0 = 0.;
    bjv1 = 0.;
    for (j = 0; j <= 1; ++j) {
	vv = (j + v0) * 4. * (j + v0);
	px = 1.;
	rp = 1.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = (doublereal) (k * 4.f - 3.f);
	    d__2 = (doublereal) (k * 4.f - 1.f);
	    rp = rp * -.0078125 * (vv - pow_dd(&d__1, &c_b4)) * (vv - pow_dd(&
		    d__2, &c_b4)) / (k * (k * 2.f - 1.f) * x2);
/* L30: */
	    px += rp;
	}
	qx = 1.;
	rq = 1.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = (doublereal) (k * 4.f - 1.f);
	    d__2 = (doublereal) (k * 4.f + 1.f);
	    rq = rq * -.0078125 * (vv - pow_dd(&d__1, &c_b4)) * (vv - pow_dd(&
		    d__2, &c_b4)) / (k * (k * 2.f + 1.f) * x2);
/* L35: */
	    qx += rq;
	}
	qx = (vv - 1.) * .125 * qx / *x;
	xk = *x - ((j + v0) * .5 + .25) * pi;
	a0 = sqrt(rp2 / *x);
	ck = cos(xk);
	sk = sin(xk);
	if (j == 0) {
	    bjv0 = a0 * (px * ck - qx * sk);
	}
	if (j == 1) {
	    bjv1 = a0 * (px * ck - qx * sk);
	}
/* L40: */
    }
    if (v0 == 0.) {
	ga = 1.;
    } else {
	gam0_(&v0, &ga);
	ga = v0 * ga;
    }
    d__1 = 2. / *x;
    fac = pow_dd(&d__1, &v0) * ga;
    vl[0] = bjv0;
    dl[0] = -bjv1 + v0 / *x * bjv0;
    vl[1] = bjv1;
    dl[1] = bjv0 - (v0 + 1.) / *x * bjv1;
    r0 = (v0 + 1.) * 2. / *x;
    if (n <= 1) {
	vl[0] = fac * vl[0];
	dl[0] = fac * dl[0] - v0 / *x * vl[0];
	vl[1] = fac * r0 * vl[1];
	dl[1] = fac * r0 * dl[1] - (v0 + 1.) / *x * vl[1];
	return 0;
    }
    if (n >= 2 && n <= (integer) (*x * .9f)) {
	f0 = bjv0;
	f1 = bjv1;
	i__1 = n;
	for (k = 2; k <= i__1; ++k) {
	    f = (k + v0 - 1.) * 2. / *x * f1 - f0;
	    f0 = f1;
	    f1 = f;
/* L45: */
	    vl[k] = f;
	}
    } else if (n >= 2) {
	m = msta1_(x, &c__200);
	if (m < n) {
	    n = m;
	} else {
	    m = msta2_(x, &n, &c__15);
	}
	f = 0.;
	f2 = 0.;
	f1 = 1e-100;
	for (k = m; k >= 0; --k) {
	    f = (v0 + k + 1.) * 2. / *x * f1 - f2;
	    if (k <= n) {
		vl[k] = f;
	    }
	    f2 = f1;
/* L50: */
	    f1 = f;
	}
	cs = 0.;
	if (abs(bjv0) > abs(bjv1)) {
	    cs = bjv0 / f;
	}
	elsecs = bjv1 / f2;
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
/* L55: */
	    vl[k] = cs * vl[k];
	}
    }
    vl[0] = fac * vl[0];
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	rc = fac * r0;
	vl[j] = rc * vl[j];
	dl[j - 1] = *x * -.5 / (j + v0) * vl[j];
/* L65: */
	r0 = (j + v0 + 1) * 2. / *x * r0;
    }
    dl[n] = (v0 + n) * 2. * (vl[n - 1] - vl[n]) / *x;
    *vm = n + v0;
    return 0;
} /* lamv_ */

/*       ********************************** */
/* Subroutine */ int chguit_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *hu, integer *id)
{
    /* Initialized data */

    static doublereal t[30] = { .0259597723012478,.0778093339495366,
	    .129449135396945,.180739964873425,.231543551376029,
	    .281722937423262,.331142848268448,.379670056576798,
	    .427173741583078,.473525841761707,.51860140005857,
	    .562278900753945,.60444059704851,.644972828489477,
	    .683766327381356,.72071651335573,.755723775306586,
	    .788693739932264,.819537526162146,.84817198478593,
	    .874519922646898,.898510310810046,.920078476177628,
	    .939166276116423,.955722255839996,.969701788765053,
	    .981067201752598,.989787895222222,.995840525118838,
	    .999210123227436 };
    static doublereal w[30] = { .0519078776312206,.0517679431749102,
	    .051488451500981,.0510701560698557,.0505141845325094,
	    .0498220356905502,.0489955754557568,.0480370318199712,
	    .0469489888489122,.0457343797161145,.0443964787957872,
	    .0429388928359356,.0413655512355848,.0396806954523808,
	    .0378888675692434,.0359948980510845,.0340038927249464,
	    .0319212190192963,.029752491500789,.0275035567499248,
	    .0251804776215213,.0227895169439978,.0203371207294572,
	    .0178299010142074,.0152746185967848,.0126781664768159,
	    .010047557182288,.00738993116334531,.00471272992695363,
	    .00202681196887362 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal c__, d__, g;
    static integer j, k, m;
    static doublereal s, a1, b1, f1, f2, t1, t2, t3, t4, ga, hu0, hu1, hu2;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute hypergeometric function U(a,b,x) by */
/*                using Gaussian-Legendre integration (n=60) */
/*       Input  : a  --- Parameter ( a > 0 ) */
/*                b  --- Parameter */
/*                x  --- Argument ( x > 0 ) */
/*       Output:  HU --- U(a,b,z) */
/*                ID --- Estimated number of significant digits */
/*       Routine called: GAMMA2 for computing Г(x) */
/*       ====================================================== */

    *id = 9;
/*       DLMF 13.4.4, integration up to C=12/X */
    a1 = *a - 1.;
    b1 = *b - *a - 1.;
    c__ = 12. / *x;
    hu0 = 0.;
    for (m = 10; m <= 100; m += 5) {
	hu1 = 0.;
	g = c__ * .5 / m;
	d__ = g;
	i__1 = m;
	for (j = 1; j <= i__1; ++j) {
	    s = 0.;
	    for (k = 1; k <= 30; ++k) {
		t1 = d__ + g * t[k - 1];
		t2 = d__ - g * t[k - 1];
		d__1 = t1 + 1.;
		f1 = exp(-(*x) * t1) * pow_dd(&t1, &a1) * pow_dd(&d__1, &b1);
		d__1 = t2 + 1.;
		f2 = exp(-(*x) * t2) * pow_dd(&t2, &a1) * pow_dd(&d__1, &b1);
		s += w[k - 1] * (f1 + f2);
/* L10: */
	    }
	    hu1 += s * g;
	    d__ += g * 2.;
/* L15: */
	}
	if ((d__1 = 1. - hu0 / hu1, abs(d__1)) < 1e-9) {
	    goto L25;
	}
	hu0 = hu1;
/* L20: */
    }
L25:
    gamma2_(a, &ga);
    hu1 /= ga;
/*       DLMF 13.4.4 with substitution t=C/(1-u) */
/*       integration u from 0 to 1, i.e. t from C=12/X to infinity */
    for (m = 2; m <= 10; m += 2) {
	hu2 = 0.;
	g = .5 / m;
	d__ = g;
	i__1 = m;
	for (j = 1; j <= i__1; ++j) {
	    s = 0.;
	    for (k = 1; k <= 30; ++k) {
		t1 = d__ + g * t[k - 1];
		t2 = d__ - g * t[k - 1];
		t3 = c__ / (1. - t1);
		t4 = c__ / (1. - t2);
		d__1 = t3 + 1.;
		f1 = t3 * t3 / c__ * exp(-(*x) * t3) * pow_dd(&t3, &a1) * 
			pow_dd(&d__1, &b1);
		d__1 = t4 + 1.;
		f2 = t4 * t4 / c__ * exp(-(*x) * t4) * pow_dd(&t4, &a1) * 
			pow_dd(&d__1, &b1);
		s += w[k - 1] * (f1 + f2);
/* L30: */
	    }
	    hu2 += s * g;
	    d__ += g * 2.;
/* L35: */
	}
	if ((d__1 = 1. - hu0 / hu2, abs(d__1)) < 1e-9) {
	    goto L45;
	}
	hu0 = hu2;
/* L40: */
    }
L45:
    gamma2_(a, &ga);
    hu2 /= ga;
    *hu = hu1 + hu2;
    return 0;
} /* chguit_ */

/*       ********************************** */
/* Subroutine */ int kmn_(integer *m, integer *n, doublereal *c__, doublereal 
	*cv, integer *kd, doublereal *df, doublereal *dn, doublereal *ck1, 
	doublereal *ck2)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double pow_ri(real *, integer *), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal r__, t, u[200], v[200], w[200], g0, r1, r2, r3, r4, r5, 
	    cs;
    static integer ip, nm, nn;
    static doublereal rk[200], tp[200], sw, gk0, gk1, gk2, gk3, sa0, sb0;
    static integer nm1;
    static doublereal su0, dnp;


/*       =================================================== */
/*       Purpose: Compute the expansion coefficients of the */
/*                prolate and oblate spheroidal functions */
/*                and joining factors */
/*       =================================================== */

    /* Parameter adjustments */
    --dn;
    --df;

    /* Function Body */
    nm = (integer) ((*n - *m) * .5f + *c__) + 25;
    nn = nm + *m;
    cs = *c__ * *c__ * *kd;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    k = 0;
    i__1 = nn + 3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ip == 0) {
	    k = (i__ - 1) * -2;
	}
	if (ip == 1) {
	    k = -((i__ << 1) - 3);
	}
	gk0 = *m * 2. + k;
	gk1 = (*m + k) * (*m + k + 1.);
	gk2 = (*m + k) * 2. - 1.;
	gk3 = (*m + k) * 2. + 3.;
	u[i__ - 1] = gk0 * (gk0 - 1.) * cs / (gk2 * (gk2 + 2.));
	v[i__ - 1] = gk1 - *cv + ((gk1 - *m * *m) * 2. - 1.) * cs / (gk2 * 
		gk3);
/* L10: */
	w[i__ - 1] = (k + 1.) * (k + 2.) * cs / ((gk2 + 2.) * gk3);
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	t = v[*m];
	i__2 = *m - k - 1;
	for (l = 0; l <= i__2; ++l) {
/* L15: */
	    t = v[*m - l - 1] - w[*m - l] * u[*m - l - 1] / t;
	}
/* L20: */
	rk[k - 1] = -u[k - 1] / t;
    }
    r__ = 1.;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	r__ *= rk[k - 1];
/* L25: */
	dn[k] = df[1] * r__;
    }
    tp[nn - 1] = v[nn];
    i__1 = *m + 1;
    for (k = nn - 1; k >= i__1; --k) {
	tp[k - 1] = v[k] - w[k + 1] * u[k] / tp[k];
	if (k > *m + 1) {
	    rk[k - 1] = -u[k - 1] / tp[k - 1];
	}
/* L30: */
    }
    if (*m == 0) {
	dnp = df[1];
    }
    if (*m != 0) {
	dnp = dn[*m];
    }
    dn[*m + 1] = pow_ii(&c_n1, &ip) * dnp * cs / ((*m * 2.f - 1.f) * (*m * 
	    2.f + 1.f - ip * 4.f) * tp[*m]);
    i__1 = nn;
    for (k = *m + 2; k <= i__1; ++k) {
/* L35: */
	dn[k] = rk[k - 1] * dn[k - 1];
    }
    r1 = 1.;
    i__1 = (*n + *m + ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L40: */
	r1 *= j + (*n + *m + ip) * .5;
    }
    nm1 = (*n - *m) / 2;
    r__ = 1.;
    i__1 = (*m << 1) + ip;
    for (j = 1; j <= i__1; ++j) {
/* L45: */
	r__ *= j;
    }
    su0 = r__ * df[1];
    sw = 0.;
    i__1 = nm;
    for (k = 2; k <= i__1; ++k) {
	r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k + ip 
		- 1.5);
	su0 += r__ * df[k];
	if (k > nm1 && (d__1 = (su0 - sw) / su0, abs(d__1)) < 1e-14) {
	    goto L55;
	}
/* L50: */
	sw = su0;
    }
L55:
    if (*kd == 1) {
	goto L70;
    }
    r2 = 1.;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
/* L60: */
	r2 = *c__ * 2. * r2 * j;
    }
    r3 = 1.;
    i__1 = (*n - *m - ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L65: */
	r3 *= j;
    }
    sa0 = ((*m + ip) * 2.f + 1.f) * r1 / (pow_ri(&c_b317, n) * pow_di(c__, &
	    ip) * r2 * r3 * df[1]);
    *ck1 = sa0 * su0;
    if (*kd == -1) {
	return 0;
    }
L70:
    r4 = 1.;
    i__1 = (*n - *m - ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L75: */
	r4 = r4 * 4. * j;
    }
    r5 = 1.;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
/* L80: */
	r5 = r5 * (j + *m) / *c__;
    }
    g0 = dn[*m];
    if (*m == 0) {
	g0 = df[1];
    }
    i__1 = ip + 1;
    sb0 = (ip + 1.f) * pow_di(c__, &i__1) / (ip * 2.f * (*m - 2.f) + 1.f) / (*
	    m * 2.f - 1.f);
    *ck2 = pow_ii(&c_n1, &ip) * sb0 * r4 * r5 * g0 / r1 * su0;
    return 0;
} /* kmn_ */

/*       ********************************** */
/* Subroutine */ int lagzo_(integer *n, doublereal *x, doublereal *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal p, q, z__, f0, f1, z0, fd, gd, pd, hn, pf;
    static integer it, nr;
    static doublereal wp;


/*       ========================================================= */
/*       Purpose : Compute the zeros of Laguerre polynomial Ln(x) */
/*                 in the interval [0,∞], and the corresponding */
/*                 weighting coefficients for Gauss-Laguerre */
/*                 integration */
/*       Input :   n    --- Order of the Laguerre polynomial */
/*                 X(n) --- Zeros of the Laguerre polynomial */
/*                 W(n) --- Corresponding weighting coefficients */
/*       ========================================================= */

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    hn = 1. / *n;
    pf = 0.;
    pd = 0.;
    i__1 = *n;
    for (nr = 1; nr <= i__1; ++nr) {
	z__ = hn;
	if (nr > 1) {
	    d__1 = (doublereal) nr;
	    z__ = x[nr - 1] + hn * pow_dd(&d__1, &c_b323);
	}
	it = 0;
L10:
	++it;
	z0 = z__;
	p = 1.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L15: */
	    p *= z__ - x[i__];
	}
	f0 = 1.;
	f1 = 1. - z__;
	i__2 = *n;
	for (k = 2; k <= i__2; ++k) {
	    pf = ((k * 2. - 1. - z__) * f1 - (k - 1.) * f0) / k;
	    pd = k / z__ * (pf - f1);
	    f0 = f1;
/* L20: */
	    f1 = pf;
	}
	fd = pf / p;
	q = 0.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wp = 1.;
	    i__3 = nr - 1;
	    for (j = 1; j <= i__3; ++j) {
		if (j == i__) {
		    goto L25;
		}
		wp *= z__ - x[j];
L25:
		;
	    }
	    q += wp;
/* L30: */
	}
	gd = (pd - q * fd) / p;
	z__ -= fd / gd;
	if (it <= 40 && (d__1 = (z__ - z0) / z__, abs(d__1)) > 1e-15) {
	    goto L10;
	}
	x[nr] = z__;
	w[nr] = 1. / (z__ * pd * pd);
/* L35: */
    }
    return 0;
} /* lagzo_ */

/*       ********************************** */
/* Subroutine */ int vvla_(doublereal *va, doublereal *x, doublereal *pv)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, a0, x1, gl, qe, pi, pdl, dsl, eps;
    extern /* Subroutine */ int dvla_(doublereal *, doublereal *, doublereal *
	    ), gamma2_(doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute parabolic cylinder function Vv(x) */
/*                for large argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PV --- Vv(x) */
/*       Routines called: */
/*             (1) DVLA for computing Dv(x) for large |x| */
/*             (2) GAMMA2 for computing Г(x) */
/*       =================================================== */

    pi = 3.141592653589793;
    eps = 1e-12;
    qe = exp(*x * .25f * *x);
    d__1 = abs(*x);
    d__2 = -(*va) - 1.;
    a0 = pow_dd(&d__1, &d__2) * sqrt(2. / pi) * qe;
    r__ = 1.;
    *pv = 1.;
    for (k = 1; k <= 18; ++k) {
	r__ = r__ * .5 * (k * 2.f + *va - 1.f) * (k * 2.f + *va) / (k * *x * *
		x);
	*pv += r__;
	if ((d__1 = r__ / *pv, abs(d__1)) < eps) {
	    goto L15;
	}
/* L10: */
    }
L15:
    *pv = a0 * *pv;
    if (*x < 0.) {
	x1 = -(*x);
	dvla_(va, &x1, &pdl);
	d__1 = -(*va);
	gamma2_(&d__1, &gl);
	dsl = sin(pi * *va) * sin(pi * *va);
	*pv = dsl * gl / pi * pdl - cos(pi * *va) * *pv;
    }
    return 0;
} /* vvla_ */

/*       ********************************** */
/* Subroutine */ int cjyva_(doublereal *v, doublecomplex *z__, doublereal *vm,
	 doublecomplex *cbj, doublecomplex *cdj, doublecomplex *cby, 
	doublecomplex *cdy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *);
    double pow_dd(doublereal *, doublereal *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_sqrt(
	    doublecomplex *, doublecomplex *), z_cos(doublecomplex *, 
	    doublecomplex *), z_sin(doublecomplex *, doublecomplex *);
    double cos(doublereal), sin(doublereal);
    void z_log(doublecomplex *, doublecomplex *), z_exp(doublecomplex *, 
	    doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer j, k, l, m, n;
    static doublereal a0;
    static integer k0;
    static doublereal v0, w0, w1;
    static doublecomplex z1, z2, ca, cb;
    static doublereal ga, gb;
    static doublecomplex cf, ci;
    static integer lb;
    static doublecomplex cr, cs;
    static doublereal wa, pi, vg, vl;
    static doublecomplex zk;
    static doublereal vv;
    static doublecomplex ca0, cf0, cf1, cf2, cg0, cg1;
    static integer lb0;
    static doublecomplex ch2, ch1, ch0, cr0, cs0, cs1, cr1;
    static doublereal ya0, ya1, rp2, pv0, pv1;
    static doublecomplex cec, cck, cp11, cp12, cp22, cp21, csk, crp, crq, cyk;
    static doublereal yak;
    static doublecomplex cpz, cqz, cju0, cjv0, cjv1, cju1, cyl1, cyl2, cyv0, 
	    cyv1, cjvl, cylk, cfac0, cfac1;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute Bessel functions Jv(z), Yv(z) and their */
/*                derivatives for a complex argument */
/*       Input :  z --- Complex argument */
/*                v --- Order of Jv(z) and Yv(z) */
/*                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 ) */
/*       Output:  CBJ(n) --- Jn+v0(z) */
/*                CDJ(n) --- Jn+v0'(z) */
/*                CBY(n) --- Yn+v0(z) */
/*                CDY(n) --- Yn+v0'(z) */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) GAMMA2 for computing the gamma function */
/*            (2) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       =========================================================== */

    pi = 3.141592653589793;
    rp2 = .63661977236758;
    ci.r = 0., ci.i = 1.;
    a0 = z_abs(z__);
    z1.r = z__->r, z1.i = z__->i;
    z__1.r = z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = z__1.r, z2.i = z__1.i;
    n = (integer) (*v);
    v0 = *v - n;
    pv0 = pi * v0;
    pv1 = pi * (v0 + 1.);
    if (a0 < 1e-100) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbj[i__2].r = 0., cbj[i__2].i = 0.;
	    i__2 = k;
	    cdj[i__2].r = 0., cdj[i__2].i = 0.;
	    i__2 = k;
	    cby[i__2].r = -1e300, cby[i__2].i = -0.;
/* L10: */
	    i__2 = k;
	    cdy[i__2].r = 1e300, cdy[i__2].i = 0.;
	}
	if (v0 == 0.f) {
	    cbj[0].r = 1., cbj[0].i = 0.;
	    cdj[1].r = .5, cdj[1].i = 0.;
	} else {
	    cdj[0].r = 1e300, cdj[0].i = 0.;
	}
	*vm = *v;
	return 0;
    }
    lb0 = 0;
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (a0 <= 12.f) {
	for (l = 0; l <= 1; ++l) {
	    vl = v0 + l;
	    cjvl.r = 1., cjvl.i = 0.;
	    cr.r = 1., cr.i = 0.;
	    for (k = 1; k <= 40; ++k) {
		z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
		z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			z2.i + z__3.i * z2.r;
		d__1 = k * (k + vl);
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr.r = z__1.r, cr.i = z__1.i;
		z__1.r = cjvl.r + cr.r, z__1.i = cjvl.i + cr.i;
		cjvl.r = z__1.r, cjvl.i = z__1.i;
		if (z_abs(&cr) < z_abs(&cjvl) * 1e-15) {
		    goto L20;
		}
/* L15: */
	    }
L20:
	    vg = vl + 1.;
	    gamma2_(&vg, &ga);
	    z__3.r = z1.r * .5, z__3.i = z1.i * .5;
	    z__4.r = vl, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z__1.r = z__2.r / ga, z__1.i = z__2.i / ga;
	    ca.r = z__1.r, ca.i = z__1.i;
	    if (l == 0) {
		z__1.r = cjvl.r * ca.r - cjvl.i * ca.i, z__1.i = cjvl.r * 
			ca.i + cjvl.i * ca.r;
		cjv0.r = z__1.r, cjv0.i = z__1.i;
	    }
	    if (l == 1) {
		z__1.r = cjvl.r * ca.r - cjvl.i * ca.i, z__1.i = cjvl.r * 
			ca.i + cjvl.i * ca.r;
		cjv1.r = z__1.r, cjv1.i = z__1.i;
	    }
/* L25: */
	}
    } else {
	k0 = 11;
	if (a0 >= 35.f) {
	    k0 = 10;
	}
	if (a0 >= 50.f) {
	    k0 = 8;
	}
	for (j = 0; j <= 1; ++j) {
	    vv = (j + v0) * 4. * (j + v0);
	    cpz.r = 1., cpz.i = 0.;
	    crp.r = 1., crp.i = 0.;
	    i__2 = k0;
	    for (k = 1; k <= i__2; ++k) {
		z__4.r = crp.r * -.0078125, z__4.i = crp.i * -.0078125;
		d__2 = (doublereal) (k * 4.f - 3.f);
		d__1 = vv - pow_dd(&d__2, &c_b4);
		z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
		d__4 = (doublereal) (k * 4.f - 1.f);
		d__3 = vv - pow_dd(&d__4, &c_b4);
		z__2.r = d__3 * z__3.r, z__2.i = d__3 * z__3.i;
		r__1 = k * (k * 2.f - 1.f);
		z__5.r = r__1 * z2.r, z__5.i = r__1 * z2.i;
		z_div(&z__1, &z__2, &z__5);
		crp.r = z__1.r, crp.i = z__1.i;
/* L30: */
		z__1.r = cpz.r + crp.r, z__1.i = cpz.i + crp.i;
		cpz.r = z__1.r, cpz.i = z__1.i;
	    }
	    cqz.r = 1., cqz.i = 0.;
	    crq.r = 1., crq.i = 0.;
	    i__2 = k0;
	    for (k = 1; k <= i__2; ++k) {
		z__4.r = crq.r * -.0078125, z__4.i = crq.i * -.0078125;
		d__2 = (doublereal) (k * 4.f - 1.f);
		d__1 = vv - pow_dd(&d__2, &c_b4);
		z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
		d__4 = (doublereal) (k * 4.f + 1.f);
		d__3 = vv - pow_dd(&d__4, &c_b4);
		z__2.r = d__3 * z__3.r, z__2.i = d__3 * z__3.i;
		r__1 = k * (k * 2.f + 1.f);
		z__5.r = r__1 * z2.r, z__5.i = r__1 * z2.i;
		z_div(&z__1, &z__2, &z__5);
		crq.r = z__1.r, crq.i = z__1.i;
/* L35: */
		z__1.r = cqz.r + crq.r, z__1.i = cqz.i + crq.i;
		cqz.r = z__1.r, cqz.i = z__1.i;
	    }
	    d__1 = (vv - 1.f) * .125;
	    z__2.r = d__1 * cqz.r, z__2.i = d__1 * cqz.i;
	    z_div(&z__1, &z__2, &z1);
	    cqz.r = z__1.r, cqz.i = z__1.i;
	    d__1 = ((j + v0) * .5 + .25) * pi;
	    z__1.r = z1.r - d__1, z__1.i = z1.i;
	    zk.r = z__1.r, zk.i = z__1.i;
	    z__3.r = rp2, z__3.i = 0.;
	    z_div(&z__2, &z__3, &z1);
	    z_sqrt(&z__1, &z__2);
	    ca0.r = z__1.r, ca0.i = z__1.i;
	    z_cos(&z__1, &zk);
	    cck.r = z__1.r, cck.i = z__1.i;
	    z_sin(&z__1, &zk);
	    csk.r = z__1.r, csk.i = z__1.i;
	    if (j == 0) {
		z__3.r = cpz.r * cck.r - cpz.i * cck.i, z__3.i = cpz.r * 
			cck.i + cpz.i * cck.r;
		z__4.r = cqz.r * csk.r - cqz.i * csk.i, z__4.i = cqz.r * 
			csk.i + cqz.i * csk.r;
		z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
		z__1.r = ca0.r * z__2.r - ca0.i * z__2.i, z__1.i = ca0.r * 
			z__2.i + ca0.i * z__2.r;
		cjv0.r = z__1.r, cjv0.i = z__1.i;
		z__3.r = cpz.r * csk.r - cpz.i * csk.i, z__3.i = cpz.r * 
			csk.i + cpz.i * csk.r;
		z__4.r = cqz.r * cck.r - cqz.i * cck.i, z__4.i = cqz.r * 
			cck.i + cqz.i * cck.r;
		z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
		z__1.r = ca0.r * z__2.r - ca0.i * z__2.i, z__1.i = ca0.r * 
			z__2.i + ca0.i * z__2.r;
		cyv0.r = z__1.r, cyv0.i = z__1.i;
	    } else if (j == 1) {
		z__3.r = cpz.r * cck.r - cpz.i * cck.i, z__3.i = cpz.r * 
			cck.i + cpz.i * cck.r;
		z__4.r = cqz.r * csk.r - cqz.i * csk.i, z__4.i = cqz.r * 
			csk.i + cqz.i * csk.r;
		z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
		z__1.r = ca0.r * z__2.r - ca0.i * z__2.i, z__1.i = ca0.r * 
			z__2.i + ca0.i * z__2.r;
		cjv1.r = z__1.r, cjv1.i = z__1.i;
		z__3.r = cpz.r * csk.r - cpz.i * csk.i, z__3.i = cpz.r * 
			csk.i + cpz.i * csk.r;
		z__4.r = cqz.r * cck.r - cqz.i * cck.i, z__4.i = cqz.r * 
			cck.i + cqz.i * cck.r;
		z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
		z__1.r = ca0.r * z__2.r - ca0.i * z__2.i, z__1.i = ca0.r * 
			z__2.i + ca0.i * z__2.r;
		cyv1.r = z__1.r, cyv1.i = z__1.i;
	    }
/* L40: */
	}
    }
    if (a0 <= 12.f) {
	if (v0 != 0.f) {
	    for (l = 0; l <= 1; ++l) {
		vl = v0 + l;
		cjvl.r = 1., cjvl.i = 0.;
		cr.r = 1., cr.i = 0.;
		for (k = 1; k <= 40; ++k) {
		    z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
		    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			    z2.i + z__3.i * z2.r;
		    d__1 = k * (k - vl);
		    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		    cr.r = z__1.r, cr.i = z__1.i;
		    z__1.r = cjvl.r + cr.r, z__1.i = cjvl.i + cr.i;
		    cjvl.r = z__1.r, cjvl.i = z__1.i;
		    if (z_abs(&cr) < z_abs(&cjvl) * 1e-15) {
			goto L50;
		    }
/* L45: */
		}
L50:
		vg = 1. - vl;
		gamma2_(&vg, &gb);
		z_div(&z__3, &c_b15, &z1);
		z__4.r = vl, z__4.i = 0.;
		pow_zz(&z__2, &z__3, &z__4);
		z__1.r = z__2.r / gb, z__1.i = z__2.i / gb;
		cb.r = z__1.r, cb.i = z__1.i;
		if (l == 0) {
		    z__1.r = cjvl.r * cb.r - cjvl.i * cb.i, z__1.i = cjvl.r * 
			    cb.i + cjvl.i * cb.r;
		    cju0.r = z__1.r, cju0.i = z__1.i;
		}
		if (l == 1) {
		    z__1.r = cjvl.r * cb.r - cjvl.i * cb.i, z__1.i = cjvl.r * 
			    cb.i + cjvl.i * cb.r;
		    cju1.r = z__1.r, cju1.i = z__1.i;
		}
/* L55: */
	    }
	    d__1 = cos(pv0);
	    z__3.r = d__1 * cjv0.r, z__3.i = d__1 * cjv0.i;
	    z__2.r = z__3.r - cju0.r, z__2.i = z__3.i - cju0.i;
	    d__2 = sin(pv0);
	    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	    d__1 = cos(pv1);
	    z__3.r = d__1 * cjv1.r, z__3.i = d__1 * cjv1.i;
	    z__2.r = z__3.r - cju1.r, z__2.i = z__3.i - cju1.i;
	    d__2 = sin(pv1);
	    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	    cyv1.r = z__1.r, cyv1.i = z__1.i;
	} else {
	    z__3.r = z1.r / 2., z__3.i = z1.i / 2.;
	    z_log(&z__2, &z__3);
	    z__1.r = z__2.r + .5772156649015329, z__1.i = z__2.i;
	    cec.r = z__1.r, cec.i = z__1.i;
	    cs0.r = 0., cs0.i = 0.;
	    w0 = 0.;
	    cr0.r = 1., cr0.i = 0.;
	    for (k = 1; k <= 30; ++k) {
		w0 += 1. / k;
		z__3.r = cr0.r * -.25, z__3.i = cr0.i * -.25;
		i__2 = k * k;
		d__1 = (doublereal) i__2;
		z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
		z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * 
			z2.i + z__2.i * z2.r;
		cr0.r = z__1.r, cr0.i = z__1.i;
/* L60: */
		z__2.r = w0 * cr0.r, z__2.i = w0 * cr0.i;
		z__1.r = cs0.r + z__2.r, z__1.i = cs0.i + z__2.i;
		cs0.r = z__1.r, cs0.i = z__1.i;
	    }
	    z__3.r = cec.r * cjv0.r - cec.i * cjv0.i, z__3.i = cec.r * cjv0.i 
		    + cec.i * cjv0.r;
	    z__2.r = z__3.r - cs0.r, z__2.i = z__3.i - cs0.i;
	    z__1.r = rp2 * z__2.r, z__1.i = rp2 * z__2.i;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	    cs1.r = 1., cs1.i = 0.;
	    w1 = 0.;
	    cr1.r = 1., cr1.i = 0.;
	    for (k = 1; k <= 30; ++k) {
		w1 += 1. / k;
		z__3.r = cr1.r * -.25, z__3.i = cr1.i * -.25;
		i__2 = k * (k + 1);
		d__1 = (doublereal) i__2;
		z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
		z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * 
			z2.i + z__2.i * z2.r;
		cr1.r = z__1.r, cr1.i = z__1.i;
/* L65: */
		d__1 = w1 * 2. + 1. / (k + 1.);
		z__2.r = d__1 * cr1.r, z__2.i = d__1 * cr1.i;
		z__1.r = cs1.r + z__2.r, z__1.i = cs1.i + z__2.i;
		cs1.r = z__1.r, cs1.i = z__1.i;
	    }
	    z__4.r = cec.r * cjv1.r - cec.i * cjv1.i, z__4.i = cec.r * cjv1.i 
		    + cec.i * cjv1.r;
	    z_div(&z__5, &c_b6, &z1);
	    z__3.r = z__4.r - z__5.r, z__3.i = z__4.i - z__5.i;
	    z__7.r = z1.r * .25, z__7.i = z1.i * .25;
	    z__6.r = z__7.r * cs1.r - z__7.i * cs1.i, z__6.i = z__7.r * cs1.i 
		    + z__7.i * cs1.r;
	    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
	    z__1.r = rp2 * z__2.r, z__1.i = rp2 * z__2.i;
	    cyv1.r = z__1.r, cyv1.i = z__1.i;
	}
    }
    if (z__->r < 0.) {
	z__2.r = pv0 * ci.r, z__2.i = pv0 * ci.i;
	z_exp(&z__1, &z__2);
	cfac0.r = z__1.r, cfac0.i = z__1.i;
	z__2.r = pv1 * ci.r, z__2.i = pv1 * ci.i;
	z_exp(&z__1, &z__2);
	cfac1.r = z__1.r, cfac1.i = z__1.i;
	if (d_imag(z__) < 0.) {
	    z__2.r = cfac0.r * cyv0.r - cfac0.i * cyv0.i, z__2.i = cfac0.r * 
		    cyv0.i + cfac0.i * cyv0.r;
	    z__5.r = ci.r * 2., z__5.i = ci.i * 2.;
	    d__1 = cos(pv0);
	    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	    z__3.r = z__4.r * cjv0.r - z__4.i * cjv0.i, z__3.i = z__4.r * 
		    cjv0.i + z__4.i * cjv0.r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	    z__2.r = cfac1.r * cyv1.r - cfac1.i * cyv1.i, z__2.i = cfac1.r * 
		    cyv1.i + cfac1.i * cyv1.r;
	    z__5.r = ci.r * 2., z__5.i = ci.i * 2.;
	    d__1 = cos(pv1);
	    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	    z__3.r = z__4.r * cjv1.r - z__4.i * cjv1.i, z__3.i = z__4.r * 
		    cjv1.i + z__4.i * cjv1.r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cyv1.r = z__1.r, cyv1.i = z__1.i;
	    z_div(&z__1, &cjv0, &cfac0);
	    cjv0.r = z__1.r, cjv0.i = z__1.i;
	    z_div(&z__1, &cjv1, &cfac1);
	    cjv1.r = z__1.r, cjv1.i = z__1.i;
	} else if (d_imag(z__) > 0.) {
	    z_div(&z__2, &cyv0, &cfac0);
	    z__5.r = ci.r * 2., z__5.i = ci.i * 2.;
	    d__1 = cos(pv0);
	    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	    z__3.r = z__4.r * cjv0.r - z__4.i * cjv0.i, z__3.i = z__4.r * 
		    cjv0.i + z__4.i * cjv0.r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	    z_div(&z__2, &cyv1, &cfac1);
	    z__5.r = ci.r * 2., z__5.i = ci.i * 2.;
	    d__1 = cos(pv1);
	    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	    z__3.r = z__4.r * cjv1.r - z__4.i * cjv1.i, z__3.i = z__4.r * 
		    cjv1.i + z__4.i * cjv1.r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    cyv1.r = z__1.r, cyv1.i = z__1.i;
	    z__1.r = cfac0.r * cjv0.r - cfac0.i * cjv0.i, z__1.i = cfac0.r * 
		    cjv0.i + cfac0.i * cjv0.r;
	    cjv0.r = z__1.r, cjv0.i = z__1.i;
	    z__1.r = cfac1.r * cjv1.r - cfac1.i * cjv1.i, z__1.i = cfac1.r * 
		    cjv1.i + cfac1.i * cjv1.r;
	    cjv1.r = z__1.r, cjv1.i = z__1.i;
	}
    }
    cbj[0].r = cjv0.r, cbj[0].i = cjv0.i;
    cbj[1].r = cjv1.r, cbj[1].i = cjv1.i;
    if (n >= 2 && n <= (integer) (a0 * .25f)) {
	cf0.r = cjv0.r, cf0.i = cjv0.i;
	cf1.r = cjv1.r, cf1.i = cjv1.i;
	i__2 = n;
	for (k = 2; k <= i__2; ++k) {
	    d__1 = (k + v0 - 1.) * 2.;
	    z__4.r = d__1, z__4.i = 0.;
	    z_div(&z__3, &z__4, z__);
	    z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i 
		    + z__3.i * cf1.r;
	    z__1.r = z__2.r - cf0.r, z__1.i = z__2.i - cf0.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	    i__1 = k;
	    cbj[i__1].r = cf.r, cbj[i__1].i = cf.i;
	    cf0.r = cf1.r, cf0.i = cf1.i;
/* L70: */
	    cf1.r = cf.r, cf1.i = cf.i;
	}
    } else if (n >= 2) {
	m = msta1_(&a0, &c__200);
	if (m < n) {
	    n = m;
	} else {
	    m = msta2_(&a0, &n, &c__15);
	}
	cf2.r = 0., cf2.i = 0.;
	cf1.r = 1e-100, cf1.i = 0.;
	for (k = m; k >= 0; --k) {
	    d__1 = (v0 + k + 1.) * 2.;
	    z__4.r = d__1, z__4.i = 0.;
	    z_div(&z__3, &z__4, z__);
	    z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i 
		    + z__3.i * cf1.r;
	    z__1.r = z__2.r - cf2.r, z__1.i = z__2.i - cf2.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	    if (k <= n) {
		i__2 = k;
		cbj[i__2].r = cf.r, cbj[i__2].i = cf.i;
	    }
	    cf2.r = cf1.r, cf2.i = cf1.i;
/* L75: */
	    cf1.r = cf.r, cf1.i = cf.i;
	}
	if (z_abs(&cjv0) > z_abs(&cjv1)) {
	    z_div(&z__1, &cjv0, &cf);
	    cs.r = z__1.r, cs.i = z__1.i;
	}
	if (z_abs(&cjv0) <= z_abs(&cjv1)) {
	    z_div(&z__1, &cjv1, &cf2);
	    cs.r = z__1.r, cs.i = z__1.i;
	}
	i__2 = n;
	for (k = 0; k <= i__2; ++k) {
/* L80: */
	    i__1 = k;
	    i__3 = k;
	    z__1.r = cs.r * cbj[i__3].r - cs.i * cbj[i__3].i, z__1.i = cs.r * 
		    cbj[i__3].i + cs.i * cbj[i__3].r;
	    cbj[i__1].r = z__1.r, cbj[i__1].i = z__1.i;
	}
    }
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbj[0].r - z__3.i * cbj[0].i, z__2.i = z__3.r * cbj[0]
	    .i + z__3.i * cbj[0].r;
    z__1.r = z__2.r - cbj[1].r, z__1.i = z__2.i - cbj[1].i;
    cdj[0].r = z__1.r, cdj[0].i = z__1.i;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L85: */
	i__3 = k;
	d__1 = -(k + v0);
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = k;
	z__2.r = z__3.r * cbj[i__2].r - z__3.i * cbj[i__2].i, z__2.i = z__3.r 
		* cbj[i__2].i + z__3.i * cbj[i__2].r;
	i__4 = k - 1;
	z__1.r = z__2.r + cbj[i__4].r, z__1.i = z__2.i + cbj[i__4].i;
	cdj[i__3].r = z__1.r, cdj[i__3].i = z__1.i;
    }
    cby[0].r = cyv0.r, cby[0].i = cyv0.i;
    cby[1].r = cyv1.r, cby[1].i = cyv1.i;
    ya0 = z_abs(&cyv0);
    lb = 0;
    cg0.r = cyv0.r, cg0.i = cyv0.i;
    cg1.r = cyv1.r, cg1.i = cyv1.i;
    i__3 = n;
    for (k = 2; k <= i__3; ++k) {
	d__1 = (v0 + k - 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	z__2.r = z__3.r * cg1.r - z__3.i * cg1.i, z__2.i = z__3.r * cg1.i + 
		z__3.i * cg1.r;
	z__1.r = z__2.r - cg0.r, z__1.i = z__2.i - cg0.i;
	cyk.r = z__1.r, cyk.i = z__1.i;
	if (z_abs(&cyk) > 1e290) {
	    goto L90;
	}
	yak = z_abs(&cyk);
	ya1 = z_abs(&cg0);
	if (yak < ya0 && yak < ya1) {
	    lb = k;
	}
	i__2 = k;
	cby[i__2].r = cyk.r, cby[i__2].i = cyk.i;
	cg0.r = cg1.r, cg0.i = cg1.i;
	cg1.r = cyk.r, cg1.i = cyk.i;
L90:
	;
    }
    if (lb <= 4 || d_imag(z__) == 0.) {
	goto L125;
    }
L95:
    if (lb == lb0) {
	goto L125;
    }
    ch2.r = 1., ch2.i = 0.;
    ch1.r = 0., ch1.i = 0.;
    lb0 = lb;
    for (k = lb; k >= 1; --k) {
	d__1 = (k + v0) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	z__2.r = z__3.r * ch1.r - z__3.i * ch1.i, z__2.i = z__3.r * ch1.i + 
		z__3.i * ch1.r;
	z__1.r = z__2.r - ch2.r, z__1.i = z__2.i - ch2.i;
	ch0.r = z__1.r, ch0.i = z__1.i;
	ch2.r = ch1.r, ch2.i = ch1.i;
/* L100: */
	ch1.r = ch0.r, ch1.i = ch0.i;
    }
    cp12.r = ch0.r, cp12.i = ch0.i;
    cp22.r = ch2.r, cp22.i = ch2.i;
    ch2.r = 0., ch2.i = 0.;
    ch1.r = 1., ch1.i = 0.;
    for (k = lb; k >= 1; --k) {
	d__1 = (k + v0) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	z__2.r = z__3.r * ch1.r - z__3.i * ch1.i, z__2.i = z__3.r * ch1.i + 
		z__3.i * ch1.r;
	z__1.r = z__2.r - ch2.r, z__1.i = z__2.i - ch2.i;
	ch0.r = z__1.r, ch0.i = z__1.i;
	ch2.r = ch1.r, ch2.i = ch1.i;
/* L105: */
	ch1.r = ch0.r, ch1.i = ch0.i;
    }
    cp11.r = ch0.r, cp11.i = ch0.i;
    cp21.r = ch2.r, cp21.i = ch2.i;
    if (lb == n) {
	i__3 = lb + 1;
	d__1 = (lb + v0) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = lb;
	z__2.r = z__3.r * cbj[i__2].r - z__3.i * cbj[i__2].i, z__2.i = z__3.r 
		* cbj[i__2].i + z__3.i * cbj[i__2].r;
	i__4 = lb - 1;
	z__1.r = z__2.r - cbj[i__4].r, z__1.i = z__2.i - cbj[i__4].i;
	cbj[i__3].r = z__1.r, cbj[i__3].i = z__1.i;
    }
    if (z_abs(cbj) > z_abs(&cbj[1])) {
	i__3 = lb + 1;
	i__2 = lb + 1;
	z__3.r = cbj[i__2].r * cyv0.r - cbj[i__2].i * cyv0.i, z__3.i = cbj[
		i__2].r * cyv0.i + cbj[i__2].i * cyv0.r;
	z__5.r = cp11.r * 2., z__5.i = cp11.i * 2.;
	z__6.r = pi * z__->r, z__6.i = pi * z__->i;
	z_div(&z__4, &z__5, &z__6);
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z_div(&z__1, &z__2, cbj);
	cby[i__3].r = z__1.r, cby[i__3].i = z__1.i;
	i__3 = lb;
	i__2 = lb;
	z__3.r = cbj[i__2].r * cyv0.r - cbj[i__2].i * cyv0.i, z__3.i = cbj[
		i__2].r * cyv0.i + cbj[i__2].i * cyv0.r;
	z__5.r = cp12.r * 2., z__5.i = cp12.i * 2.;
	z__6.r = pi * z__->r, z__6.i = pi * z__->i;
	z_div(&z__4, &z__5, &z__6);
	z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	z_div(&z__1, &z__2, cbj);
	cby[i__3].r = z__1.r, cby[i__3].i = z__1.i;
    } else {
	i__3 = lb + 1;
	i__2 = lb + 1;
	z__3.r = cbj[i__2].r * cyv1.r - cbj[i__2].i * cyv1.i, z__3.i = cbj[
		i__2].r * cyv1.i + cbj[i__2].i * cyv1.r;
	z__5.r = cp21.r * 2., z__5.i = cp21.i * 2.;
	z__6.r = pi * z__->r, z__6.i = pi * z__->i;
	z_div(&z__4, &z__5, &z__6);
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z_div(&z__1, &z__2, &cbj[1]);
	cby[i__3].r = z__1.r, cby[i__3].i = z__1.i;
	i__3 = lb;
	i__2 = lb;
	z__3.r = cbj[i__2].r * cyv1.r - cbj[i__2].i * cyv1.i, z__3.i = cbj[
		i__2].r * cyv1.i + cbj[i__2].i * cyv1.r;
	z__5.r = cp22.r * 2., z__5.i = cp22.i * 2.;
	z__6.r = pi * z__->r, z__6.i = pi * z__->i;
	z_div(&z__4, &z__5, &z__6);
	z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	z_div(&z__1, &z__2, &cbj[1]);
	cby[i__3].r = z__1.r, cby[i__3].i = z__1.i;
    }
    i__3 = lb + 1;
    cyl2.r = cby[i__3].r, cyl2.i = cby[i__3].i;
    i__3 = lb;
    cyl1.r = cby[i__3].r, cyl1.i = cby[i__3].i;
    for (k = lb - 1; k >= 0; --k) {
	d__1 = (k + v0 + 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	z__2.r = z__3.r * cyl1.r - z__3.i * cyl1.i, z__2.i = z__3.r * cyl1.i 
		+ z__3.i * cyl1.r;
	z__1.r = z__2.r - cyl2.r, z__1.i = z__2.i - cyl2.i;
	cylk.r = z__1.r, cylk.i = z__1.i;
	i__3 = k;
	cby[i__3].r = cylk.r, cby[i__3].i = cylk.i;
	cyl2.r = cyl1.r, cyl2.i = cyl1.i;
/* L110: */
	cyl1.r = cylk.r, cyl1.i = cylk.i;
    }
    i__3 = lb;
    cyl1.r = cby[i__3].r, cyl1.i = cby[i__3].i;
    i__3 = lb + 1;
    cyl2.r = cby[i__3].r, cyl2.i = cby[i__3].i;
    i__3 = n - 1;
    for (k = lb + 1; k <= i__3; ++k) {
	d__1 = (k + v0) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	z__2.r = z__3.r * cyl2.r - z__3.i * cyl2.i, z__2.i = z__3.r * cyl2.i 
		+ z__3.i * cyl2.r;
	z__1.r = z__2.r - cyl1.r, z__1.i = z__2.i - cyl1.i;
	cylk.r = z__1.r, cylk.i = z__1.i;
	i__2 = k + 1;
	cby[i__2].r = cylk.r, cby[i__2].i = cylk.i;
	cyl1.r = cyl2.r, cyl1.i = cyl2.i;
/* L115: */
	cyl2.r = cylk.r, cyl2.i = cylk.i;
    }
    i__3 = n;
    for (k = 2; k <= i__3; ++k) {
	wa = z_abs(&cby[k]);
	if (wa < z_abs(&cby[k - 1])) {
	    lb = k;
	}
/* L120: */
    }
    goto L95;
L125:
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cby[0].r - z__3.i * cby[0].i, z__2.i = z__3.r * cby[0]
	    .i + z__3.i * cby[0].r;
    z__1.r = z__2.r - cby[1].r, z__1.i = z__2.i - cby[1].i;
    cdy[0].r = z__1.r, cdy[0].i = z__1.i;
    i__3 = n;
    for (k = 1; k <= i__3; ++k) {
/* L130: */
	i__2 = k;
	i__4 = k - 1;
	d__1 = k + v0;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__1 = k;
	z__2.r = z__3.r * cby[i__1].r - z__3.i * cby[i__1].i, z__2.i = z__3.r 
		* cby[i__1].i + z__3.i * cby[i__1].r;
	z__1.r = cby[i__4].r - z__2.r, z__1.i = cby[i__4].i - z__2.i;
	cdy[i__2].r = z__1.r, cdy[i__2].i = z__1.i;
    }
    *vm = n + v0;
    return 0;
} /* cjyva_ */

/*       ********************************** */
/* Subroutine */ int cjyvb_(doublereal *v, doublecomplex *z__, doublereal *vm,
	 doublecomplex *cbj, doublecomplex *cdj, doublecomplex *cby, 
	doublecomplex *cdy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *);
    double pow_dd(doublereal *, doublereal *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_sqrt(
	    doublecomplex *, doublecomplex *), z_cos(doublecomplex *, 
	    doublecomplex *), z_sin(doublecomplex *, doublecomplex *);
    double cos(doublereal), sin(doublereal);
    void z_log(doublecomplex *, doublecomplex *), z_exp(doublecomplex *, 
	    doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer k, m, n;
    static doublereal a0;
    static integer k0;
    static doublereal v0, w0;
    static doublecomplex z1, z2, ca, cb;
    static doublereal ga, gb;
    static doublecomplex cf, ci, cr, cs;
    static doublereal pi, vg;
    static doublecomplex zk;
    static doublereal vv;
    static doublecomplex ca0, cf1, cf2, cr0, cs0;
    static doublereal rp2, pv0;
    static doublecomplex cec, cck, csk, crp, crq, cpz, cqz, cyy, cju0, cjv0, 
	    cyv0, cjvn, cfac0;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute Bessel functions Jv(z), Yv(z) and their */
/*                derivatives for a complex argument */
/*       Input :  z --- Complex argument */
/*                v --- Order of Jv(z) and Yv(z) */
/*                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 ) */
/*       Output:  CBJ(n) --- Jn+v0(z) */
/*                CDJ(n) --- Jn+v0'(z) */
/*                CBY(n) --- Yn+v0(z) */
/*                CDY(n) --- Yn+v0'(z) */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) GAMMA2 for computing the gamma function */
/*            (2) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       =========================================================== */

    pi = 3.141592653589793;
    rp2 = .63661977236758;
    ci.r = 0., ci.i = 1.;
    a0 = z_abs(z__);
    z1.r = z__->r, z1.i = z__->i;
    z__1.r = z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = z__1.r, z2.i = z__1.i;
    n = (integer) (*v);
    v0 = *v - n;
    pv0 = pi * v0;
    if (a0 < 1e-100) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbj[i__2].r = 0., cbj[i__2].i = 0.;
	    i__2 = k;
	    cdj[i__2].r = 0., cdj[i__2].i = 0.;
	    i__2 = k;
	    cby[i__2].r = -1e300, cby[i__2].i = -0.;
/* L10: */
	    i__2 = k;
	    cdy[i__2].r = 1e300, cdy[i__2].i = 0.;
	}
	if (v0 == 0.f) {
	    cbj[0].r = 1., cbj[0].i = 0.;
	    cdj[1].r = .5, cdj[1].i = 0.;
	} else {
	    cdj[0].r = 1e300, cdj[0].i = 0.;
	}
	*vm = *v;
	return 0;
    }
    if (z__->r < 0.) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (a0 <= 12.f) {
	cjv0.r = 1., cjv0.i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 40; ++k) {
	    z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    d__1 = k * (k + v0);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cjv0.r + cr.r, z__1.i = cjv0.i + cr.i;
	    cjv0.r = z__1.r, cjv0.i = z__1.i;
	    if (z_abs(&cr) < z_abs(&cjv0) * 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	vg = v0 + 1.;
	gamma2_(&vg, &ga);
	z__3.r = z1.r * .5, z__3.i = z1.i * .5;
	z__4.r = v0, z__4.i = 0.;
	pow_zz(&z__2, &z__3, &z__4);
	z__1.r = z__2.r / ga, z__1.i = z__2.i / ga;
	ca.r = z__1.r, ca.i = z__1.i;
	z__1.r = cjv0.r * ca.r - cjv0.i * ca.i, z__1.i = cjv0.r * ca.i + 
		cjv0.i * ca.r;
	cjv0.r = z__1.r, cjv0.i = z__1.i;
    } else {
	k0 = 11;
	if (a0 >= 35.f) {
	    k0 = 10;
	}
	if (a0 >= 50.f) {
	    k0 = 8;
	}
	vv = v0 * 4. * v0;
	cpz.r = 1., cpz.i = 0.;
	crp.r = 1., crp.i = 0.;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
	    z__4.r = crp.r * -.0078125, z__4.i = crp.i * -.0078125;
	    d__2 = (doublereal) (k * 4.f - 3.f);
	    d__1 = vv - pow_dd(&d__2, &c_b4);
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    d__4 = (doublereal) (k * 4.f - 1.f);
	    d__3 = vv - pow_dd(&d__4, &c_b4);
	    z__2.r = d__3 * z__3.r, z__2.i = d__3 * z__3.i;
	    r__1 = k * (k * 2.f - 1.f);
	    z__5.r = r__1 * z2.r, z__5.i = r__1 * z2.i;
	    z_div(&z__1, &z__2, &z__5);
	    crp.r = z__1.r, crp.i = z__1.i;
/* L25: */
	    z__1.r = cpz.r + crp.r, z__1.i = cpz.i + crp.i;
	    cpz.r = z__1.r, cpz.i = z__1.i;
	}
	cqz.r = 1., cqz.i = 0.;
	crq.r = 1., crq.i = 0.;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
	    z__4.r = crq.r * -.0078125, z__4.i = crq.i * -.0078125;
	    d__2 = (doublereal) (k * 4.f - 1.f);
	    d__1 = vv - pow_dd(&d__2, &c_b4);
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    d__4 = (doublereal) (k * 4.f + 1.f);
	    d__3 = vv - pow_dd(&d__4, &c_b4);
	    z__2.r = d__3 * z__3.r, z__2.i = d__3 * z__3.i;
	    r__1 = k * (k * 2.f + 1.f);
	    z__5.r = r__1 * z2.r, z__5.i = r__1 * z2.i;
	    z_div(&z__1, &z__2, &z__5);
	    crq.r = z__1.r, crq.i = z__1.i;
/* L30: */
	    z__1.r = cqz.r + crq.r, z__1.i = cqz.i + crq.i;
	    cqz.r = z__1.r, cqz.i = z__1.i;
	}
	d__1 = (vv - 1.f) * .125;
	z__2.r = d__1 * cqz.r, z__2.i = d__1 * cqz.i;
	z_div(&z__1, &z__2, &z1);
	cqz.r = z__1.r, cqz.i = z__1.i;
	d__1 = (v0 * .5 + .25) * pi;
	z__1.r = z1.r - d__1, z__1.i = z1.i;
	zk.r = z__1.r, zk.i = z__1.i;
	z__3.r = rp2, z__3.i = 0.;
	z_div(&z__2, &z__3, &z1);
	z_sqrt(&z__1, &z__2);
	ca0.r = z__1.r, ca0.i = z__1.i;
	z_cos(&z__1, &zk);
	cck.r = z__1.r, cck.i = z__1.i;
	z_sin(&z__1, &zk);
	csk.r = z__1.r, csk.i = z__1.i;
	z__3.r = cpz.r * cck.r - cpz.i * cck.i, z__3.i = cpz.r * cck.i + 
		cpz.i * cck.r;
	z__4.r = cqz.r * csk.r - cqz.i * csk.i, z__4.i = cqz.r * csk.i + 
		cqz.i * csk.r;
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z__1.r = ca0.r * z__2.r - ca0.i * z__2.i, z__1.i = ca0.r * z__2.i + 
		ca0.i * z__2.r;
	cjv0.r = z__1.r, cjv0.i = z__1.i;
	z__3.r = cpz.r * csk.r - cpz.i * csk.i, z__3.i = cpz.r * csk.i + 
		cpz.i * csk.r;
	z__4.r = cqz.r * cck.r - cqz.i * cck.i, z__4.i = cqz.r * cck.i + 
		cqz.i * cck.r;
	z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	z__1.r = ca0.r * z__2.r - ca0.i * z__2.i, z__1.i = ca0.r * z__2.i + 
		ca0.i * z__2.r;
	cyv0.r = z__1.r, cyv0.i = z__1.i;
    }
    if (a0 <= 12.f) {
	if (v0 != 0.f) {
	    cjvn.r = 1., cjvn.i = 0.;
	    cr.r = 1., cr.i = 0.;
	    for (k = 1; k <= 40; ++k) {
		z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
		z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			z2.i + z__3.i * z2.r;
		d__1 = k * (k - v0);
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr.r = z__1.r, cr.i = z__1.i;
		z__1.r = cjvn.r + cr.r, z__1.i = cjvn.i + cr.i;
		cjvn.r = z__1.r, cjvn.i = z__1.i;
		if (z_abs(&cr) < z_abs(&cjvn) * 1e-15) {
		    goto L40;
		}
/* L35: */
	    }
L40:
	    vg = 1. - v0;
	    gamma2_(&vg, &gb);
	    z_div(&z__3, &c_b15, &z1);
	    z__4.r = v0, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z__1.r = z__2.r / gb, z__1.i = z__2.i / gb;
	    cb.r = z__1.r, cb.i = z__1.i;
	    z__1.r = cjvn.r * cb.r - cjvn.i * cb.i, z__1.i = cjvn.r * cb.i + 
		    cjvn.i * cb.r;
	    cju0.r = z__1.r, cju0.i = z__1.i;
	    d__1 = cos(pv0);
	    z__3.r = d__1 * cjv0.r, z__3.i = d__1 * cjv0.i;
	    z__2.r = z__3.r - cju0.r, z__2.i = z__3.i - cju0.i;
	    d__2 = sin(pv0);
	    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	} else {
	    z__3.r = z1.r / 2., z__3.i = z1.i / 2.;
	    z_log(&z__2, &z__3);
	    z__1.r = z__2.r + .5772156649015329, z__1.i = z__2.i;
	    cec.r = z__1.r, cec.i = z__1.i;
	    cs0.r = 0., cs0.i = 0.;
	    w0 = 0.;
	    cr0.r = 1., cr0.i = 0.;
	    for (k = 1; k <= 30; ++k) {
		w0 += 1. / k;
		z__3.r = cr0.r * -.25, z__3.i = cr0.i * -.25;
		i__2 = k * k;
		d__1 = (doublereal) i__2;
		z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
		z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * 
			z2.i + z__2.i * z2.r;
		cr0.r = z__1.r, cr0.i = z__1.i;
/* L45: */
		z__2.r = w0 * cr0.r, z__2.i = w0 * cr0.i;
		z__1.r = cs0.r + z__2.r, z__1.i = cs0.i + z__2.i;
		cs0.r = z__1.r, cs0.i = z__1.i;
	    }
	    z__3.r = cec.r * cjv0.r - cec.i * cjv0.i, z__3.i = cec.r * cjv0.i 
		    + cec.i * cjv0.r;
	    z__2.r = z__3.r - cs0.r, z__2.i = z__3.i - cs0.i;
	    z__1.r = rp2 * z__2.r, z__1.i = rp2 * z__2.i;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	}
    }
    if (n == 0) {
	n = 1;
    }
    m = msta1_(&a0, &c__200);
    if (m < n) {
	n = m;
    } else {
	m = msta2_(&a0, &n, &c__15);
    }
    cf2.r = 0., cf2.i = 0.;
    cf1.r = 1e-100, cf1.i = 0.;
    for (k = m; k >= 0; --k) {
	d__1 = (v0 + k + 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, &z1);
	z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i + 
		z__3.i * cf1.r;
	z__1.r = z__2.r - cf2.r, z__1.i = z__2.i - cf2.i;
	cf.r = z__1.r, cf.i = z__1.i;
	if (k <= n) {
	    i__2 = k;
	    cbj[i__2].r = cf.r, cbj[i__2].i = cf.i;
	}
	cf2.r = cf1.r, cf2.i = cf1.i;
/* L50: */
	cf1.r = cf.r, cf1.i = cf.i;
    }
    z_div(&z__1, &cjv0, &cf);
    cs.r = z__1.r, cs.i = z__1.i;
    i__2 = n;
    for (k = 0; k <= i__2; ++k) {
/* L55: */
	i__1 = k;
	i__3 = k;
	z__1.r = cs.r * cbj[i__3].r - cs.i * cbj[i__3].i, z__1.i = cs.r * cbj[
		i__3].i + cs.i * cbj[i__3].r;
	cbj[i__1].r = z__1.r, cbj[i__1].i = z__1.i;
    }
    if (z__->r < 0.) {
	z__2.r = pv0 * ci.r, z__2.i = pv0 * ci.i;
	z_exp(&z__1, &z__2);
	cfac0.r = z__1.r, cfac0.i = z__1.i;
	if (d_imag(z__) < 0.) {
	    z__2.r = cfac0.r * cyv0.r - cfac0.i * cyv0.i, z__2.i = cfac0.r * 
		    cyv0.i + cfac0.i * cyv0.r;
	    z__5.r = ci.r * 2., z__5.i = ci.i * 2.;
	    d__1 = cos(pv0);
	    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	    z__3.r = z__4.r * cjv0.r - z__4.i * cjv0.i, z__3.i = z__4.r * 
		    cjv0.i + z__4.i * cjv0.r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	} else if (d_imag(z__) > 0.) {
	    z_div(&z__2, &cyv0, &cfac0);
	    z__5.r = ci.r * 2., z__5.i = ci.i * 2.;
	    d__1 = cos(pv0);
	    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	    z__3.r = z__4.r * cjv0.r - z__4.i * cjv0.i, z__3.i = z__4.r * 
		    cjv0.i + z__4.i * cjv0.r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    cyv0.r = z__1.r, cyv0.i = z__1.i;
	}
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    if (d_imag(z__) < 0.) {
		i__3 = k;
		d__1 = -pi * (k + v0);
		z__3.r = d__1 * ci.r, z__3.i = d__1 * ci.i;
		z_exp(&z__2, &z__3);
		i__2 = k;
		z__1.r = z__2.r * cbj[i__2].r - z__2.i * cbj[i__2].i, z__1.i =
			 z__2.r * cbj[i__2].i + z__2.i * cbj[i__2].r;
		cbj[i__3].r = z__1.r, cbj[i__3].i = z__1.i;
	    } else if (d_imag(z__) > 0.) {
		i__3 = k;
		d__1 = pi * (k + v0);
		z__3.r = d__1 * ci.r, z__3.i = d__1 * ci.i;
		z_exp(&z__2, &z__3);
		i__2 = k;
		z__1.r = z__2.r * cbj[i__2].r - z__2.i * cbj[i__2].i, z__1.i =
			 z__2.r * cbj[i__2].i + z__2.i * cbj[i__2].r;
		cbj[i__3].r = z__1.r, cbj[i__3].i = z__1.i;
	    }
/* L60: */
	}
	z1.r = z1.r, z1.i = z1.i;
    }
    cby[0].r = cyv0.r, cby[0].i = cyv0.i;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__3 = k;
	i__2 = k - 1;
	z__3.r = cbj[i__3].r * cby[i__2].r - cbj[i__3].i * cby[i__2].i, 
		z__3.i = cbj[i__3].r * cby[i__2].i + cbj[i__3].i * cby[i__2]
		.r;
	z__5.r = pi * z__->r, z__5.i = pi * z__->i;
	z_div(&z__4, &c_b15, &z__5);
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z_div(&z__1, &z__2, &cbj[k - 1]);
	cyy.r = z__1.r, cyy.i = z__1.i;
	i__3 = k;
	cby[i__3].r = cyy.r, cby[i__3].i = cyy.i;
/* L65: */
    }
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbj[0].r - z__3.i * cbj[0].i, z__2.i = z__3.r * cbj[0]
	    .i + z__3.i * cbj[0].r;
    z__1.r = z__2.r - cbj[1].r, z__1.i = z__2.i - cbj[1].i;
    cdj[0].r = z__1.r, cdj[0].i = z__1.i;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L70: */
	i__3 = k;
	d__1 = -(k + v0);
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = k;
	z__2.r = z__3.r * cbj[i__2].r - z__3.i * cbj[i__2].i, z__2.i = z__3.r 
		* cbj[i__2].i + z__3.i * cbj[i__2].r;
	i__4 = k - 1;
	z__1.r = z__2.r + cbj[i__4].r, z__1.i = z__2.i + cbj[i__4].i;
	cdj[i__3].r = z__1.r, cdj[i__3].i = z__1.i;
    }
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cby[0].r - z__3.i * cby[0].i, z__2.i = z__3.r * cby[0]
	    .i + z__3.i * cby[0].r;
    z__1.r = z__2.r - cby[1].r, z__1.i = z__2.i - cby[1].i;
    cdy[0].r = z__1.r, cdy[0].i = z__1.i;
    i__3 = n;
    for (k = 1; k <= i__3; ++k) {
/* L75: */
	i__2 = k;
	i__4 = k - 1;
	d__1 = k + v0;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__1 = k;
	z__2.r = z__3.r * cby[i__1].r - z__3.i * cby[i__1].i, z__2.i = z__3.r 
		* cby[i__1].i + z__3.i * cby[i__1].r;
	z__1.r = cby[i__4].r - z__2.r, z__1.i = cby[i__4].i - z__2.i;
	cdy[i__2].r = z__1.r, cdy[i__2].i = z__1.i;
    }
    *vm = n + v0;
    return 0;
} /* cjyvb_ */

/*       ********************************** */
/* Subroutine */ int jy01a_(doublereal *x, doublereal *bj0, doublereal *dj0, 
	doublereal *bj1, doublereal *dj1, doublereal *by0, doublereal *dy0, 
	doublereal *by1, doublereal *dy1)
{
    /* Initialized data */

    static doublereal a[12] = { -.0703125,.112152099609375,-.5725014209747314,
	    6.074042001273483,-110.0171402692467,3038.090510922384,
	    -118838.4262567832,6252951.493434797,-425939216.5047669,
	    36468400807.06556,-3833534661393.944,485401468685290.1 };
    static doublereal b[12] = { .0732421875,-.2271080017089844,
	    1.727727502584457,-24.38052969955606,551.3358961220206,
	    -18257.75547429318,832859.3040162893,-50069589.53198893,
	    3836255180.230433,-364901081884.9833,42189715702840.96,
	    -5827244631566907. };
    static doublereal a1[12] = { .1171875,-.144195556640625,.6765925884246826,
	    -6.883914268109947,121.5978918765359,-3302.272294480852,
	    127641.2726461746,-6656367.718817688,450278600.3050393,
	    -38338575207.4279,4011838599133.198,-506056850331472.7 };
    static doublereal b1[12] = { -.1025390625,.2775764465332031,
	    -1.993531733751297,27.24882731126854,-603.8440767050702,
	    19718.37591223663,-890297.8767070678,53104110.10968522,
	    -4043620325.107754,382701134659.8605,-44064814178522.78,
	    6065091351222699. };

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), sqrt(doublereal),
	     cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__;
    static integer k0;
    static doublereal p0, q0, r0, r1, p1, t1, t2, w0, w1, q1, x2, ec, cu, pi, 
	    cs0, cs1, rp2;


/*       ======================================================= */
/*       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x), */
/*                Y1(x), and their derivatives */
/*       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ≥ 0 ) */
/*       Output:  BJ0 --- J0(x) */
/*                DJ0 --- J0'(x) */
/*                BJ1 --- J1(x) */
/*                DJ1 --- J1'(x) */
/*                BY0 --- Y0(x) */
/*                DY0 --- Y0'(x) */
/*                BY1 --- Y1(x) */
/*                DY1 --- Y1'(x) */
/*       ======================================================= */

    pi = 3.141592653589793;
    rp2 = .63661977236758;
    x2 = *x * *x;
    if (*x == 0.) {
	*bj0 = 1.;
	*bj1 = 0.;
	*dj0 = 0.;
	*dj1 = .5;
	*by0 = -1e300;
	*by1 = -1e300;
	*dy0 = 1e300;
	*dy1 = 1e300;
	return 0;
    }
    if (*x <= 12.) {
	*bj0 = 1.;
	r__ = 1.;
	for (k = 1; k <= 30; ++k) {
	    r__ = r__ * -.25 * x2 / (k * k);
	    *bj0 += r__;
	    if (abs(r__) < abs(*bj0) * 1e-15) {
		goto L10;
	    }
/* L5: */
	}
L10:
	*bj1 = 1.;
	r__ = 1.;
	for (k = 1; k <= 30; ++k) {
	    r__ = r__ * -.25 * x2 / (k * (k + 1.));
	    *bj1 += r__;
	    if (abs(r__) < abs(*bj1) * 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	*bj1 = *x * .5 * *bj1;
	ec = log(*x / 2.) + .5772156649015329;
	cs0 = 0.;
	w0 = 0.;
	r0 = 1.;
	for (k = 1; k <= 30; ++k) {
	    w0 += 1. / k;
	    r0 = r0 * -.25 / (k * k) * x2;
	    r__ = r0 * w0;
	    cs0 += r__;
	    if (abs(r__) < abs(cs0) * 1e-15) {
		goto L30;
	    }
/* L25: */
	}
L30:
	*by0 = rp2 * (ec * *bj0 - cs0);
	cs1 = 1.;
	w1 = 0.;
	r1 = 1.;
	for (k = 1; k <= 30; ++k) {
	    w1 += 1. / k;
	    r1 = r1 * -.25 / (k * (k + 1)) * x2;
	    r__ = r1 * (w1 * 2. + 1. / (k + 1.));
	    cs1 += r__;
	    if (abs(r__) < abs(cs1) * 1e-15) {
		goto L40;
	    }
/* L35: */
	}
L40:
	*by1 = rp2 * (ec * *bj1 - 1. / *x - *x * .25 * cs1);
    } else {
	k0 = 12;
	if (*x >= 35.f) {
	    k0 = 10;
	}
	if (*x >= 50.f) {
	    k0 = 8;
	}
	t1 = *x - pi * .25;
	p0 = 1.;
	q0 = -.125 / *x;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = k * -2;
	    p0 += a[k - 1] * pow_di(x, &i__2);
/* L45: */
	    i__2 = k * -2 - 1;
	    q0 += b[k - 1] * pow_di(x, &i__2);
	}
	cu = sqrt(rp2 / *x);
	*bj0 = cu * (p0 * cos(t1) - q0 * sin(t1));
	*by0 = cu * (p0 * sin(t1) + q0 * cos(t1));
	t2 = *x - pi * .75;
	p1 = 1.;
	q1 = .375 / *x;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = k * -2;
	    p1 += a1[k - 1] * pow_di(x, &i__1);
/* L50: */
	    i__1 = k * -2 - 1;
	    q1 += b1[k - 1] * pow_di(x, &i__1);
	}
	cu = sqrt(rp2 / *x);
	*bj1 = cu * (p1 * cos(t2) - q1 * sin(t2));
	*by1 = cu * (p1 * sin(t2) + q1 * cos(t2));
    }
    *dj0 = -(*bj1);
    *dj1 = *bj0 - *bj1 / *x;
    *dy0 = -(*by1);
    *dy1 = *by0 - *by1 / *x;
    return 0;
} /* jy01a_ */

/*       ********************************** */
/* Subroutine */ int incog_(doublereal *a, doublereal *x, doublereal *gin, 
	doublereal *gim, doublereal *gip, integer *isfer)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, s, t0, ga, xam;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute the incomplete gamma function */
/*                r(a,x), Г(a,x) and P(a,x) */
/*       Input :  a   --- Parameter ( a ≤ 170 ) */
/*                x   --- Argument */
/*       Output:  GIN --- r(a,x) */
/*                GIM --- Г(a,x) */
/*                GIP --- P(a,x) */
/*                ISFER --- Error flag */
/*       Routine called: GAMMA2 for computing Г(x) */
/*       =================================================== */

    *isfer = 0;
    xam = -(*x) + *a * log(*x);
    if (xam > 700.f || *a > 170.f) {
	*isfer = 6;
	return 0;
    }
    if (*x == 0.f) {
	*gin = 0.f;
	gamma2_(a, &ga);
	*gim = ga;
	*gip = 0.f;
    } else if (*x <= *a + 1.f) {
	s = 1. / *a;
	r__ = s;
	for (k = 1; k <= 60; ++k) {
	    r__ = r__ * *x / (*a + k);
	    s += r__;
	    if ((d__1 = r__ / s, abs(d__1)) < 1e-15) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*gin = exp(xam) * s;
	gamma2_(a, &ga);
	*gip = *gin / ga;
	*gim = ga - *gin;
    } else if (*x > *a + 1.f) {
	t0 = 0.;
	for (k = 60; k >= 1; --k) {
	    t0 = (k - *a) / (k / (*x + t0) + 1.);
/* L20: */
	}
	*gim = exp(xam) / (*x + t0);
	gamma2_(a, &ga);
	*gin = ga - *gim;
	*gip = 1. - *gim / ga;
    }
    return 0;
} /* incog_ */

/*       ********************************** */
/* Subroutine */ int itikb_(doublereal *x, doublereal *ti, doublereal *tk)
{
    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), log(doublereal);

    /* Local variables */
    static doublereal t, t1, pi;


/*       ======================================================= */
/*       Purpose: Integrate Bessel functions I0(t) and K0(t) */
/*                with respect to t from 0 to x */
/*       Input :  x  --- Upper limit of the integral ( x ≥ 0 ) */
/*       Output:  TI --- Integration of I0(t) from 0 to x */
/*                TK --- Integration of K0(t) from 0 to x */
/*       ======================================================= */

    pi = 3.141592653589793;
    if (*x == 0.) {
	*ti = 0.;
    } else if (*x < 5.) {
	t1 = *x / 5.;
	t = t1 * t1;
	*ti = ((((((((t * 5.9434e-4 + .004500642) * t + .044686921) * t + 
		.300704878) * t + 1.471860153) * t + 4.844024624) * t + 
		9.765629849) * t + 10.416666367) * t + 5.) * t1;
    } else if (*x >= 5.f && *x <= 8.) {
	t = 5. / *x;
	*ti = (((t * -.015166 - .0202292) * t + .1294122) * t - .0302912) * t 
		+ .4161224;
	*ti = *ti * exp(*x) / sqrt(*x);
    } else {
	t = 8. / *x;
	*ti = (((((t * -.0073995 + .017744) * t - .0114858) * t + .0055956) * 
		t + .0059191) * t + .0311734) * t + .3989423;
	*ti = *ti * exp(*x) / sqrt(*x);
    }
    if (*x == 0.) {
	*tk = 0.;
    } else if (*x <= 2.) {
	t1 = *x / 2.;
	t = t1 * t1;
	*tk = ((((((t * 1.16e-6 + 2.069e-5) * t + 6.2664e-4) * t + .01110118) 
		* t + .11227902) * t + .50407836) * t + .84556868) * t1;
	*tk -= log(*x / 2.) * *ti;
    } else if (*x > 2.f && *x <= 4.) {
	t = 2. / *x;
	*tk = (((t * .0160395 - .0781715) * t + .185984) * t - .3584641) * t 
		+ 1.2494934;
	*tk = pi / 2. - *tk * exp(-(*x)) / sqrt(*x);
    } else if (*x > 4.f && *x <= 7.) {
	t = 4. / *x;
	*tk = (((((t * .0037128 - .0158449) * t + .0320504) * t - .0481455) * 
		t + .0787284) * t - .1958273) * t + 1.2533141;
	*tk = pi / 2. - *tk * exp(-(*x)) / sqrt(*x);
    } else {
	t = 7. / *x;
	*tk = (((((t * 3.3934e-4 - .00163271) * t + .00417454) * t - 
		.00933944) * t + .02576646) * t - .11190289) * t + 1.25331414;
	*tk = pi / 2. - *tk * exp(-(*x)) / sqrt(*x);
    }
    return 0;
} /* itikb_ */

/*       ********************************** */
/* Subroutine */ int itika_(doublereal *x, doublereal *ti, doublereal *tk)
{
    /* Initialized data */

    static doublereal a[10] = { .625,1.0078125,2.5927734375,9.1868591308594,
	    41.567974090576,229.19635891914,1491.504060477,11192.354495579,
	    95159.39374212,904124.25769041 };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), log(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, b1, b2, e0, x2, el, pi, rs, tw, rc1, rc2;


/*       ======================================================= */
/*       Purpose: Integrate modified Bessel functions I0(t) and */
/*                K0(t) with respect to t from 0 to x */
/*       Input :  x  --- Upper limit of the integral  ( x ≥ 0 ) */
/*       Output:  TI --- Integration of I0(t) from 0 to x */
/*                TK --- Integration of K0(t) from 0 to x */
/*       ======================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    if (*x == 0.) {
	*ti = 0.;
	*tk = 0.;
	return 0;
    } else if (*x < 20.) {
	x2 = *x * *x;
	*ti = 1.;
	r__ = 1.;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * .25 * ((k << 1) - 1.) / ((k << 1) + 1.) / (k * k) * 
		    x2;
	    *ti += r__;
	    if ((d__1 = r__ / *ti, abs(d__1)) < 1e-12) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*ti *= *x;
    } else {
	x2 = 0.;
	*ti = 1.;
	r__ = 1.;
	for (k = 1; k <= 10; ++k) {
	    r__ /= *x;
/* L20: */
	    *ti += a[k - 1] * r__;
	}
	rc1 = 1. / sqrt(pi * 2. * *x);
	*ti = rc1 * exp(*x) * *ti;
    }
    if (*x < 12.) {
	e0 = el + log(*x / 2.);
	b1 = 1. - e0;
	b2 = 0.;
	rs = 0.;
	r__ = 1.;
	tw = 0.;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * .25 * ((k << 1) - 1.) / ((k << 1) + 1.) / (k * k) * 
		    x2;
	    b1 += r__ * (1. / ((k << 1) + 1) - e0);
	    rs += 1. / k;
	    b2 += r__ * rs;
	    *tk = b1 + b2;
	    if ((d__1 = (*tk - tw) / *tk, abs(d__1)) < 1e-12) {
		goto L30;
	    }
/* L25: */
	    tw = *tk;
	}
L30:
	*tk *= *x;
    } else {
	*tk = 1.;
	r__ = 1.;
	for (k = 1; k <= 10; ++k) {
	    r__ = -r__ / *x;
/* L35: */
	    *tk += a[k - 1] * r__;
	}
	rc2 = sqrt(pi / (*x * 2.));
	*tk = pi / 2. - rc2 * *tk * exp(-(*x));
    }
    return 0;
} /* itika_ */

/*       ********************************** */
/* Subroutine */ int jyv_(doublereal *v, doublereal *x, doublereal *vm, 
	doublereal *bj, doublereal *dj, doublereal *by, doublereal *dy)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), cos(
	    doublereal), sin(doublereal), log(doublereal);

    /* Local variables */
    static doublereal a, b, f;
    static integer j, k, l, m, n;
    static doublereal r__, a0, f0, f1, f2;
    static integer k0;
    static doublereal r0, r1, v0, w0, w1, x2, ga, gb, ec, ck, el, cs, pi, vg, 
	    sk, vl, rp, rq, xk, px, qx, vv, cs0, cs1, rp2, pv0, pv1, bju0, 
	    bjv0, bjv1, bju1, byv0, byv1, bjvl, byvk;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ======================================================= */
/*       Purpose: Compute Bessel functions Jv(x) and Yv(x) */
/*                and their derivatives */
/*       Input :  x --- Argument of Jv(x) and Yv(x) */
/*                v --- Order of Jv(x) and Yv(x) */
/*                      ( v = n+v0, 0 ≤ v0 < 1, n = 0,1,2,... ) */
/*       Output:  BJ(n) --- Jn+v0(x) */
/*                DJ(n) --- Jn+v0'(x) */
/*                BY(n) --- Yn+v0(x) */
/*                DY(n) --- Yn+v0'(x) */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) GAMMA2 for computing gamma function */
/*            (2) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       ======================================================= */

    el = .5772156649015329;
    pi = 3.141592653589793;
    rp2 = .63661977236758;
    x2 = *x * *x;
    n = (integer) (*v);
    v0 = *v - n;
    if (*x < 1e-100) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    bj[k] = 0.;
	    dj[k] = 0.;
	    by[k] = -1e300;
/* L10: */
	    dy[k] = 1e300;
	}
	if (v0 == 0.f) {
	    bj[0] = 1.;
	    dj[1] = .5;
	} else {
	    dj[0] = 1e300;
	}
	*vm = *v;
	return 0;
    }
    bjv0 = 0.;
    bjv1 = 0.;
    byv0 = 0.;
    byv1 = 0.;
    if (*x <= 12.f) {
	for (l = 0; l <= 1; ++l) {
	    vl = v0 + l;
	    bjvl = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 40; ++k) {
		r__ = r__ * -.25 * x2 / (k * (k + vl));
		bjvl += r__;
		if (abs(r__) < abs(bjvl) * 1e-15) {
		    goto L20;
		}
/* L15: */
	    }
L20:
	    vg = vl + 1.;
	    gamma2_(&vg, &ga);
	    d__1 = *x * .5;
	    a = pow_dd(&d__1, &vl) / ga;
	    if (l == 0) {
		bjv0 = bjvl * a;
	    }
	    if (l == 1) {
		bjv1 = bjvl * a;
	    }
/* L25: */
	}
    } else {
	k0 = 11;
	if (*x >= 35.f) {
	    k0 = 10;
	}
	if (*x >= 50.f) {
	    k0 = 8;
	}
	for (j = 0; j <= 1; ++j) {
	    vv = (j + v0) * 4. * (j + v0);
	    px = 1.;
	    rp = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = (doublereal) (k * 4.f - 3.f);
		d__2 = (doublereal) (k * 4.f - 1.f);
		rp = rp * -.0078125 * (vv - pow_dd(&d__1, &c_b4)) * (vv - 
			pow_dd(&d__2, &c_b4)) / (k * (k * 2.f - 1.f) * x2);
/* L30: */
		px += rp;
	    }
	    qx = 1.;
	    rq = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = (doublereal) (k * 4.f - 1.f);
		d__2 = (doublereal) (k * 4.f + 1.f);
		rq = rq * -.0078125 * (vv - pow_dd(&d__1, &c_b4)) * (vv - 
			pow_dd(&d__2, &c_b4)) / (k * (k * 2.f + 1.f) * x2);
/* L35: */
		qx += rq;
	    }
	    qx = (vv - 1.f) * .125 * qx / *x;
	    xk = *x - ((j + v0) * .5 + .25) * pi;
	    a0 = sqrt(rp2 / *x);
	    ck = cos(xk);
	    sk = sin(xk);
	    if (j == 0) {
		bjv0 = a0 * (px * ck - qx * sk);
		byv0 = a0 * (px * sk + qx * ck);
	    } else if (j == 1) {
		bjv1 = a0 * (px * ck - qx * sk);
		byv1 = a0 * (px * sk + qx * ck);
	    }
/* L40: */
	}
    }
    bj[0] = bjv0;
    bj[1] = bjv1;
    dj[0] = v0 / *x * bj[0] - bj[1];
    dj[1] = -(v0 + 1.) / *x * bj[1] + bj[0];
    if (n >= 2 && n <= (integer) (*x * .9f)) {
	f0 = bjv0;
	f1 = bjv1;
	i__1 = n;
	for (k = 2; k <= i__1; ++k) {
	    f = (k + v0 - 1.) * 2. / *x * f1 - f0;
	    bj[k] = f;
	    f0 = f1;
/* L45: */
	    f1 = f;
	}
    } else if (n >= 2) {
	m = msta1_(x, &c__200);
	if (m < n) {
	    n = m;
	} else {
	    m = msta2_(x, &n, &c__15);
	}
	f = 0.;
	f2 = 0.;
	f1 = 1e-100;
	for (k = m; k >= 0; --k) {
	    f = (v0 + k + 1.) * 2. / *x * f1 - f2;
	    if (k <= n) {
		bj[k] = f;
	    }
	    f2 = f1;
/* L50: */
	    f1 = f;
	}
	if (abs(bjv0) > abs(bjv1)) {
	    cs = bjv0 / f;
	} else {
	    cs = bjv1 / f2;
	}
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
/* L55: */
	    bj[k] = cs * bj[k];
	}
    }
    i__1 = n;
    for (k = 2; k <= i__1; ++k) {
/* L60: */
	dj[k] = -(k + v0) / *x * bj[k] + bj[k - 1];
    }
    if (*x <= 12.) {
	if (v0 != 0.f) {
	    bju0 = 0.;
	    bju1 = 0.;
	    for (l = 0; l <= 1; ++l) {
		vl = v0 + l;
		bjvl = 1.;
		r__ = 1.;
		for (k = 1; k <= 40; ++k) {
		    r__ = r__ * -.25 * x2 / (k * (k - vl));
		    bjvl += r__;
		    if (abs(r__) < abs(bjvl) * 1e-15) {
			goto L70;
		    }
/* L65: */
		}
L70:
		vg = 1. - vl;
		gamma2_(&vg, &gb);
		d__1 = 2. / *x;
		b = pow_dd(&d__1, &vl) / gb;
		if (l == 0) {
		    bju0 = bjvl * b;
		}
		if (l == 1) {
		    bju1 = bjvl * b;
		}
/* L75: */
	    }
	    pv0 = pi * v0;
	    pv1 = pi * (v0 + 1.);
	    byv0 = (bjv0 * cos(pv0) - bju0) / sin(pv0);
	    byv1 = (bjv1 * cos(pv1) - bju1) / sin(pv1);
	} else {
	    ec = log(*x / 2.) + el;
	    cs0 = 0.;
	    w0 = 0.;
	    r0 = 1.;
	    for (k = 1; k <= 30; ++k) {
		w0 += 1. / k;
		r0 = r0 * -.25 / (k * k) * x2;
/* L80: */
		cs0 += r0 * w0;
	    }
	    byv0 = rp2 * (ec * bjv0 - cs0);
	    cs1 = 1.;
	    w1 = 0.;
	    r1 = 1.;
	    for (k = 1; k <= 30; ++k) {
		w1 += 1. / k;
		r1 = r1 * -.25 / (k * (k + 1)) * x2;
/* L85: */
		cs1 += r1 * (w1 * 2. + 1. / (k + 1.));
	    }
	    byv1 = rp2 * (ec * bjv1 - 1. / *x - *x * .25 * cs1);
	}
    }
    by[0] = byv0;
    by[1] = byv1;
    i__1 = n;
    for (k = 2; k <= i__1; ++k) {
	byvk = (v0 + k - 1.) * 2. / *x * byv1 - byv0;
	by[k] = byvk;
	byv0 = byv1;
/* L90: */
	byv1 = byvk;
    }
    dy[0] = v0 / *x * by[0] - by[1];
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L95: */
	dy[k] = -(k + v0) / *x * by[k] + by[k - 1];
    }
    *vm = n + v0;
    return 0;
} /* jyv_ */

/*       ********************************** */
/* Subroutine */ int jynb_(integer *n, doublereal *x, integer *nm, doublereal 
	*bj, doublereal *dj, doublereal *by, doublereal *dy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int jynbh_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);


/*       ===================================================== */
/*       Purpose: Compute Bessel functions Jn(x), Yn(x) and */
/*                their derivatives */
/*       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 ) */
/*                n --- Order of Jn(x) and Yn(x) */
/*       Output:  BJ(n) --- Jn(x) */
/*                DJ(n) --- Jn'(x) */
/*                BY(n) --- Yn(x) */
/*                DY(n) --- Yn'(x) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                JYNBH to calculate the Jn and Yn */
/*       ===================================================== */

    jynbh_(n, &c__0, x, nm, bj, by);
/*       Compute derivatives by differentiation formulas */
    if (*x < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    dj[k] = 0.;
/* L10: */
	    dy[k] = 1e300;
	}
	dj[1] = .5;
    } else {
	dj[0] = -bj[1];
	i__1 = *nm;
	for (k = 1; k <= i__1; ++k) {
/* L40: */
	    dj[k] = bj[k - 1] - k / *x * bj[k];
	}
	dy[0] = -by[1];
	i__1 = *nm;
	for (k = 1; k <= i__1; ++k) {
/* L50: */
	    dy[k] = by[k - 1] - k * by[k] / *x;
	}
    }
    return 0;
} /* jynb_ */

/*       ********************************** */
/* Subroutine */ int jynbh_(integer *n, integer *nmin, doublereal *x, integer 
	*nm, doublereal *bj, doublereal *by)
{
    /* Initialized data */

    static doublereal a[4] = { -.0703125,.112152099609375,-.5725014209747314,
	    6.074042001273483 };
    static doublereal b[4] = { .0732421875,-.2271080017089844,
	    1.727727502584457,-24.38052969955606 };
    static doublereal a1[4] = { .1171875,-.144195556640625,.6765925884246826,
	    -6.883914268109947 };
    static doublereal b1[4] = { -.1025390625,.2775764465332031,
	    -1.993531733751297,27.24882731126854 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double log(doublereal), pow_di(doublereal *, integer *), sqrt(doublereal),
	     cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, m;
    static doublereal f1, f2, p0, q0, p1, s0, q1, t1, t2, ec, bs, cu, pi;
    static integer ky;
    static doublereal su, sv, bj0, bj1, by0, by1, r2p, bjk, byk;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ===================================================== */
/*       Purpose: Compute Bessel functions Jn(x), Yn(x) */
/*       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 ) */
/*                n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 ) */
/*                nmin -- Lowest order computed  ( nmin ≥ 0 ) */
/*       Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0 */
/*                BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0 */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 to calculate the starting */
/*                point for backward recurrence */
/*       ===================================================== */

    pi = 3.141592653589793;
    r2p = .63661977236758;
    *nm = *n;
    if (*x < 1e-100) {
	i__1 = *n;
	for (k = *nmin; k <= i__1; ++k) {
	    bj[k - *nmin] = 0.;
/* L10: */
	    by[k - *nmin] = -1e300;
	}
	if (*nmin == 0) {
	    bj[0] = 1.;
	}
	return 0;
    }
    if (*x <= 300.f || *n > (integer) (*x * .9f)) {
/*          Backward recurrence for Jn */
	if (*n == 0) {
	    *nm = 1;
	}
	m = msta1_(x, &c__200);
	if (m < *nm) {
	    *nm = m;
	} else {
	    m = msta2_(x, nm, &c__15);
	}
	bs = 0.;
	su = 0.;
	sv = 0.;
	f2 = 0.;
	f1 = 1e-100;
	f = 0.;
	for (k = m; k >= 0; --k) {
	    f = (k + 1.) * 2. / *x * f1 - f2;
	    if (k <= *nm && k >= *nmin) {
		bj[k - *nmin] = f;
	    }
	    if (k == k / 2 << 1 && k != 0) {
		bs += f * 2.;
		i__1 = k / 2;
		su += pow_ii(&c_n1, &i__1) * f / k;
	    } else if (k > 1) {
		i__1 = k / 2;
		sv += pow_ii(&c_n1, &i__1) * k / (k * k - 1.) * f;
	    }
	    f2 = f1;
/* L15: */
	    f1 = f;
	}
	s0 = bs + f;
	i__1 = *nm;
	for (k = *nmin; k <= i__1; ++k) {
/* L20: */
	    bj[k - *nmin] /= s0;
	}
/*          Estimates for Yn at start of recurrence */
	bj0 = f1 / s0;
	bj1 = f2 / s0;
	ec = log(*x / 2.) + .5772156649015329;
	by0 = r2p * (ec * bj0 - su * 4. / s0);
	by1 = r2p * ((ec - 1.) * bj1 - bj0 / *x - sv * 4. / s0);
	if (0 >= *nmin) {
	    by[-(*nmin)] = by0;
	}
	if (1 >= *nmin) {
	    by[1 - *nmin] = by1;
	}
	ky = 2;
    } else {
/*          Hankel expansion */
	t1 = *x - pi * .25;
	p0 = 1.;
	q0 = -.125 / *x;
	for (k = 1; k <= 4; ++k) {
	    i__1 = k * -2;
	    p0 += a[k - 1] * pow_di(x, &i__1);
/* L25: */
	    i__1 = k * -2 - 1;
	    q0 += b[k - 1] * pow_di(x, &i__1);
	}
	cu = sqrt(r2p / *x);
	bj0 = cu * (p0 * cos(t1) - q0 * sin(t1));
	by0 = cu * (p0 * sin(t1) + q0 * cos(t1));
	if (0 >= *nmin) {
	    bj[-(*nmin)] = bj0;
	}
	if (0 >= *nmin) {
	    by[-(*nmin)] = by0;
	}
	t2 = *x - pi * .75;
	p1 = 1.;
	q1 = .375 / *x;
	for (k = 1; k <= 4; ++k) {
	    i__1 = k * -2;
	    p1 += a1[k - 1] * pow_di(x, &i__1);
/* L30: */
	    i__1 = k * -2 - 1;
	    q1 += b1[k - 1] * pow_di(x, &i__1);
	}
	bj1 = cu * (p1 * cos(t2) - q1 * sin(t2));
	by1 = cu * (p1 * sin(t2) + q1 * cos(t2));
	if (1 >= *nmin) {
	    bj[1 - *nmin] = bj1;
	}
	if (1 >= *nmin) {
	    by[1 - *nmin] = by1;
	}
	i__1 = *nm;
	for (k = 2; k <= i__1; ++k) {
	    bjk = (k - 1.) * 2. / *x * bj1 - bj0;
	    if (k >= *nmin) {
		bj[k - *nmin] = bjk;
	    }
	    bj0 = bj1;
/* L35: */
	    bj1 = bjk;
	}
	ky = 2;
    }
/*       Forward recurrence for Yn */
    i__1 = *nm;
    for (k = ky; k <= i__1; ++k) {
	byk = (k - 1.) * 2. * by1 / *x - by0;
	if (k >= *nmin) {
	    by[k - *nmin] = byk;
	}
	by0 = by1;
/* L45: */
	by1 = byk;
    }
    return 0;
} /* jynbh_ */

/*       ********************************** */
/* Subroutine */ int legzo_(integer *n, doublereal *x, doublereal *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal p, q, z__, f0, f1;
    static integer n0;
    static doublereal z0, fd, gd, pd, pf;
    static integer nr;
    static doublereal wp;


/*       ========================================================= */
/*       Purpose : Compute the zeros of Legendre polynomial Pn(x) */
/*                 in the interval [-1,1], and the corresponding */
/*                 weighting coefficients for Gauss-Legendre */
/*                 integration */
/*       Input :   n    --- Order of the Legendre polynomial */
/*       Output:   X(n) --- Zeros of the Legendre polynomial */
/*                 W(n) --- Corresponding weighting coefficients */
/*       ========================================================= */

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    n0 = (*n + 1) / 2;
    pf = 0.;
    pd = 0.;
    i__1 = n0;
    for (nr = 1; nr <= i__1; ++nr) {
	z__ = cos((nr - .25) * 3.1415926 / *n);
L10:
	z0 = z__;
	p = 1.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L15: */
	    p *= z__ - x[i__];
	}
	f0 = 1.;
	if (nr == n0 && *n != *n / 2 << 1) {
	    z__ = 0.;
	}
	f1 = z__;
	i__2 = *n;
	for (k = 2; k <= i__2; ++k) {
	    pf = (2. - 1. / k) * z__ * f1 - (1. - 1. / k) * f0;
	    pd = k * (f1 - z__ * pf) / (1. - z__ * z__);
	    f0 = f1;
/* L20: */
	    f1 = pf;
	}
	if (z__ == 0.f) {
	    goto L40;
	}
	fd = pf / p;
	q = 0.;
	i__2 = nr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wp = 1.;
	    i__3 = nr;
	    for (j = 1; j <= i__3; ++j) {
		if (j != i__) {
		    wp *= z__ - x[j];
		}
/* L30: */
	    }
/* L35: */
	    q += wp;
	}
	gd = (pd - q * fd) / p;
	z__ -= fd / gd;
	if ((d__1 = z__ - z0, abs(d__1)) > abs(z__) * 1e-15) {
	    goto L10;
	}
L40:
	x[nr] = z__;
	x[*n + 1 - nr] = -z__;
	w[nr] = 2. / ((1. - z__ * z__) * pd * pd);
/* L45: */
	w[*n + 1 - nr] = w[nr];
    }
    return 0;
} /* legzo_ */

/*       ********************************** */
/* Subroutine */ int aswfa_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, integer *kd, doublereal *cv, doublereal *s1f, 
	doublereal *s1d)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;

    /* Local variables */
    static integer k;
    static doublereal r__, a0, d0, d1, x0, x1, df[200], ck[200];
    static integer ip, nm, nm2;
    static doublereal su1, su2, eps;
    extern /* Subroutine */ int sckb_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), sdmn_(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute the prolate and oblate spheroidal angular */
/*                functions of the first kind and their derivatives */
/*       Input :  m  --- Mode parameter,  m = 0,1,2,... */
/*                n  --- Mode parameter,  n = m,m+1,... */
/*                c  --- Spheroidal parameter */
/*                x  --- Argument of angular function, |x| < 1.0 */
/*                KD --- Function code */
/*                       KD=1 for prolate;  KD=-1 for oblate */
/*                cv --- Characteristic value */
/*       Output:  S1F --- Angular function of the first kind */
/*                S1D --- Derivative of the angular function of */
/*                        the first kind */
/*       Routine called: */
/*                SCKB for computing expansion coefficients ck */
/*       =========================================================== */

    eps = 1e-14;
    x0 = *x;
    *x = abs(*x);
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    nm = (integer) ((*n - *m) / 2 + *c__) + 40;
    nm2 = nm / 2 - 2;
    sdmn_(m, n, c__, cv, kd, df);
    sckb_(m, n, c__, df, ck);
    x1 = 1. - *x * *x;
    if (*m == 0 && x1 == 0.) {
	a0 = 1.;
    } else {
	d__1 = *m * .5;
	a0 = pow_dd(&x1, &d__1);
    }
    su1 = ck[0];
    i__1 = nm2;
    for (k = 1; k <= i__1; ++k) {
	r__ = ck[k] * pow_di(&x1, &k);
	su1 += r__;
	if (k >= 10 && (d__1 = r__ / su1, abs(d__1)) < eps) {
	    goto L15;
	}
/* L10: */
    }
L15:
    *s1f = a0 * pow_di(x, &ip) * su1;
    if (*x == 1.) {
	if (*m == 0) {
	    *s1d = ip * ck[0] - ck[1] * 2.;
	}
	if (*m == 1) {
	    *s1d = -1e100;
	}
	if (*m == 2) {
	    *s1d = ck[0] * -2.;
	}
	if (*m >= 3) {
	    *s1d = 0.;
	}
    } else {
	d__1 = ip + 1.;
	d0 = ip - *m / x1 * pow_dd(x, &d__1);
	d__1 = ip + 1.;
	d1 = a0 * -2. * pow_dd(x, &d__1);
	su2 = ck[1];
	i__1 = nm2;
	for (k = 2; k <= i__1; ++k) {
	    d__1 = k - 1.;
	    r__ = k * ck[k] * pow_dd(&x1, &d__1);
	    su2 += r__;
	    if (k >= 10 && (d__1 = r__ / su2, abs(d__1)) < eps) {
		goto L25;
	    }
/* L20: */
	}
L25:
	*s1d = d0 * a0 * su1 + d1 * su2;
    }
    if (x0 < 0. && ip == 0) {
	*s1d = -(*s1d);
    }
    if (x0 < 0. && ip == 1) {
	*s1f = -(*s1f);
    }
    *x = x0;
    return 0;
} /* aswfa_ */

/*       ********************************** */
/* Subroutine */ int jyna_(integer *n, doublereal *x, integer *nm, doublereal 
	*bj, doublereal *dj, doublereal *by, doublereal *dy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal f;
    static integer k, m;
    static doublereal f0, f1, f2, cs, bj0, bj1, dj0, dj1, by0, by1, dy0, dy1, 
	    bjk;
    extern /* Subroutine */ int jy01b_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ========================================================== */
/*       Purpose: Compute Bessel functions Jn(x) & Yn(x) and */
/*                their derivatives */
/*       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ≥ 0 ) */
/*                n --- Order of Jn(x) & Yn(x) */
/*       Output:  BJ(n) --- Jn(x) */
/*                DJ(n) --- Jn'(x) */
/*                BY(n) --- Yn(x) */
/*                DY(n) --- Yn'(x) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*            (1) JY01B to calculate J0(x), J1(x), Y0(x) & Y1(x) */
/*            (2) MSTA1 and MSTA2 to calculate the starting */
/*                point for backward recurrence */
/*       ========================================================= */

    *nm = *n;
    if (*x < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    bj[k] = 0.;
	    dj[k] = 0.;
	    by[k] = -1e300;
/* L10: */
	    dy[k] = 1e300;
	}
	bj[0] = 1.;
	dj[1] = .5;
	return 0;
    }
    jy01b_(x, &bj0, &dj0, &bj1, &dj1, &by0, &dy0, &by1, &dy1);
    bj[0] = bj0;
    bj[1] = bj1;
    by[0] = by0;
    by[1] = by1;
    dj[0] = dj0;
    dj[1] = dj1;
    dy[0] = dy0;
    dy[1] = dy1;
    if (*n <= 1) {
	return 0;
    }
    if (*n < (integer) (*x * .9f)) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    bjk = (k - 1.) * 2. / *x * bj1 - bj0;
	    bj[k] = bjk;
	    bj0 = bj1;
/* L20: */
	    bj1 = bjk;
	}
    } else {
	m = msta1_(x, &c__200);
	if (m < *n) {
	    *nm = m;
	} else {
	    m = msta2_(x, n, &c__15);
	}
	f2 = 0.;
	f1 = 1e-100;
	f = 0.;
	for (k = m; k >= 0; --k) {
	    f = (k + 1.) * 2. / *x * f1 - f2;
	    if (k <= *nm) {
		bj[k] = f;
	    }
	    f2 = f1;
/* L30: */
	    f1 = f;
	}
	if (abs(bj0) > abs(bj1)) {
	    cs = bj0 / f;
	} else {
	    cs = bj1 / f2;
	}
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
/* L40: */
	    bj[k] = cs * bj[k];
	}
    }
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
/* L50: */
	dj[k] = bj[k - 1] - k / *x * bj[k];
    }
    f0 = by[0];
    f1 = by[1];
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
	f = (k - 1.) * 2. / *x * f1 - f0;
	by[k] = f;
	f0 = f1;
/* L60: */
	f1 = f;
    }
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
/* L70: */
	dy[k] = by[k - 1] - k * by[k] / *x;
    }
    return 0;
} /* jyna_ */

/*       ********************************** */
/* Subroutine */ int pbdv_(doublereal *v, doublereal *x, doublereal *dv, 
	doublereal *dp, doublereal *pdf, doublereal *pdd)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), exp(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, l, m;
    static doublereal f0, f1, s0, v0, v1, v2;
    static integer ja, na;
    static doublereal ep, pd, xa;
    static integer nk;
    static doublereal vh;
    static integer nv;
    static doublereal pd0, pd1;
    extern /* Subroutine */ int dvla_(doublereal *, doublereal *, doublereal *
	    ), dvsa_(doublereal *, doublereal *, doublereal *);


/*       ==================================================== */
/*       Purpose: Compute parabolic cylinder functions Dv(x) */
/*                and their derivatives */
/*       Input:   x --- Argument of Dv(x) */
/*                v --- Order of Dv(x) */
/*       Output:  DV(na) --- Dn+v0(x) */
/*                DP(na) --- Dn+v0'(x) */
/*                ( na = |n|, v0 = v-n, |v0| < 1, */
/*                  n = 0,±1,±2,… ) */
/*                PDF --- Dv(x) */
/*                PDD --- Dv'(x) */
/*       Routines called: */
/*             (1) DVSA for computing Dv(x) for small |x| */
/*             (2) DVLA for computing Dv(x) for large |x| */
/*       ==================================================== */

    xa = abs(*x);
    vh = *v;
    *v += d_sign(&c_b524, v);
    nv = (integer) (*v);
    v0 = *v - nv;
    na = abs(nv);
    ep = exp(*x * -.25 * *x);
    ja = 0;
    if (na >= 1) {
	ja = 1;
    }
    if (*v >= 0.f) {
	if (v0 == 0.f) {
	    pd0 = ep;
	    pd1 = *x * ep;
	} else {
	    i__1 = ja;
	    for (l = 0; l <= i__1; ++l) {
		v1 = v0 + l;
		if (xa <= 5.8f) {
		    dvsa_(&v1, x, &pd1);
		}
		if (xa > 5.8f) {
		    dvla_(&v1, x, &pd1);
		}
		if (l == 0) {
		    pd0 = pd1;
		}
/* L10: */
	    }
	}
	dv[0] = pd0;
	dv[1] = pd1;
	i__1 = na;
	for (k = 2; k <= i__1; ++k) {
	    *pdf = *x * pd1 - (k + v0 - 1.) * pd0;
	    dv[k] = *pdf;
	    pd0 = pd1;
/* L15: */
	    pd1 = *pdf;
	}
    } else {
	if (*x <= 0.f) {
	    if (xa <= 5.8) {
		dvsa_(&v0, x, &pd0);
		v1 = v0 - 1.;
		dvsa_(&v1, x, &pd1);
	    } else {
		dvla_(&v0, x, &pd0);
		v1 = v0 - 1.;
		dvla_(&v1, x, &pd1);
	    }
	    dv[0] = pd0;
	    dv[1] = pd1;
	    i__1 = na;
	    for (k = 2; k <= i__1; ++k) {
		pd = (-(*x) * pd1 + pd0) / (k - 1. - v0);
		dv[k] = pd;
		pd0 = pd1;
/* L20: */
		pd1 = pd;
	    }
	} else if (*x <= 2.f) {
	    v2 = nv + v0;
	    if (nv == 0) {
		v2 += -1.;
	    }
	    nk = (integer) (-v2);
	    dvsa_(&v2, x, &f1);
	    v1 = v2 + 1.;
	    dvsa_(&v1, x, &f0);
	    dv[nk] = f1;
	    dv[nk - 1] = f0;
	    for (k = nk - 2; k >= 0; --k) {
		f = *x * f0 + (k - v0 + 1.) * f1;
		dv[k] = f;
		f1 = f0;
/* L25: */
		f0 = f;
	    }
	} else {
	    if (xa <= 5.8f) {
		dvsa_(&v0, x, &pd0);
	    }
	    if (xa > 5.8f) {
		dvla_(&v0, x, &pd0);
	    }
	    dv[0] = pd0;
	    m = na + 100;
	    f1 = 0.;
	    f0 = 1e-30;
	    f = 0.;
	    for (k = m; k >= 0; --k) {
		f = *x * f0 + (k - v0 + 1.) * f1;
		if (k <= na) {
		    dv[k] = f;
		}
		f1 = f0;
/* L30: */
		f0 = f;
	    }
	    s0 = pd0 / f;
	    i__1 = na;
	    for (k = 0; k <= i__1; ++k) {
/* L35: */
		dv[k] = s0 * dv[k];
	    }
	}
    }
    i__1 = na - 1;
    for (k = 0; k <= i__1; ++k) {
	v1 = abs(v0) + k;
	if (*v >= 0.) {
	    dp[k] = *x * .5 * dv[k] - dv[k + 1];
	} else {
	    dp[k] = *x * -.5 * dv[k] - v1 * dv[k + 1];
	}
/* L40: */
    }
    *pdf = dv[na - 1];
    *pdd = dp[na - 1];
    *v = vh;
    return 0;
} /* pbdv_ */

/*       ********************************** */
/* Subroutine */ int itsh0_(doublereal *x, doublereal *th0)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal), sqrt(doublereal)
	    ;

    /* Local variables */
    static doublereal a[25];
    static integer k;
    static doublereal r__, s, a0, a1, s0, af, bf, bg, el, rd, pi, xp, ty;


/*       =================================================== */
/*       Purpose: Evaluate the integral of Struve function */
/*                H0(t) with respect to t from 0 and x */
/*       Input :  x   --- Upper limit  ( x ≥ 0 ) */
/*       Output:  TH0 --- Integration of H0(t) from 0 and x */
/*       =================================================== */

    pi = 3.141592653589793;
    r__ = 1.;
    if (*x <= 30.f) {
	s = .5;
	for (k = 1; k <= 100; ++k) {
	    rd = 1.;
	    if (k == 1) {
		rd = .5;
	    }
/* Computing 2nd power */
	    d__1 = *x / (k * 2. + 1.);
	    r__ = -r__ * rd * k / (k + 1.) * (d__1 * d__1);
	    s += r__;
	    if (abs(r__) < abs(s) * 1e-12) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*th0 = 2. / pi * *x * *x * s;
    } else {
	s = 1.;
	for (k = 1; k <= 12; ++k) {
/* Computing 2nd power */
	    d__1 = (k * 2. + 1.) / *x;
	    r__ = -r__ * k / (k + 1.) * (d__1 * d__1);
	    s += r__;
	    if (abs(r__) < abs(s) * 1e-12) {
		goto L25;
	    }
/* L20: */
	}
L25:
	el = .57721566490153;
	s0 = s / (pi * *x * *x) + 2. / pi * (log(*x * 2.) + el);
	a0 = 1.;
	a1 = .625;
	a[0] = a1;
	for (k = 1; k <= 20; ++k) {
	    af = ((k + .5) * 1.5 * (k + .83333333333333337) * a1 - (k + .5) * 
		    .5 * (k + .5) * (k - .5) * a0) / (k + 1.);
	    a[k] = af;
	    a0 = a1;
/* L30: */
	    a1 = af;
	}
	bf = 1.;
	r__ = 1.;
	for (k = 1; k <= 10; ++k) {
	    r__ = -r__ / (*x * *x);
/* L35: */
	    bf += a[(k << 1) - 1] * r__;
	}
	bg = a[0] / *x;
	r__ = 1. / *x;
	for (k = 1; k <= 10; ++k) {
	    r__ = -r__ / (*x * *x);
/* L40: */
	    bg += a[k * 2] * r__;
	}
	xp = *x + pi * .25;
	ty = sqrt(2. / (pi * *x)) * (bg * cos(xp) - bf * sin(xp));
	*th0 = ty + s0;
    }
    return 0;
} /* itsh0_ */

/*       ********************************** */
/* Subroutine */ int cerzo_(integer *nt, doublecomplex *zo)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal w;
    static doublecomplex z__;
    static doublereal w0, pi;
    static integer it;
    static doublecomplex zd;
    static integer nr;
    static doublecomplex zf;
    static doublereal pu, pv, px, py;
    static doublecomplex zp, zq, zw, zfd, zgd;
    extern /* Subroutine */ int cerf_(doublecomplex *, doublecomplex *, 
	    doublecomplex *);


/*       =============================================================== */
/*       Purpose : Evaluate the complex zeros of error function erf(z) */
/*                 using the modified Newton's iteration method */
/*       Input :   NT --- Total number of zeros */
/*       Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT */
/*       Routine called: CERF for computing erf(z) and erf'(z) */
/*       =============================================================== */

    /* Parameter adjustments */
    --zo;

    /* Function Body */
    pi = 3.141592653589793;
    w = 0.;
    i__1 = *nt;
    for (nr = 1; nr <= i__1; ++nr) {
	pu = sqrt(pi * (nr * 4. - .5));
	pv = pi * sqrt(nr * 2. - .25);
	px = pu * .5f - log(pv) * .5f / pu;
	py = pu * .5f + log(pv) * .5f / pu;
	z__1.r = px, z__1.i = py;
	z__.r = z__1.r, z__.i = z__1.i;
	it = 0;
L15:
	++it;
	cerf_(&z__, &zf, &zd);
	zp.r = 1., zp.i = 0.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    i__3 = i__;
	    z__2.r = z__.r - zo[i__3].r, z__2.i = z__.i - zo[i__3].i;
	    z__1.r = zp.r * z__2.r - zp.i * z__2.i, z__1.i = zp.r * z__2.i + 
		    zp.i * z__2.r;
	    zp.r = z__1.r, zp.i = z__1.i;
	}
	z_div(&z__1, &zf, &zp);
	zfd.r = z__1.r, zfd.i = z__1.i;
	zq.r = 0., zq.i = 0.;
	i__3 = nr - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    zw.r = 1., zw.i = 0.;
	    i__2 = nr - 1;
	    for (j = 1; j <= i__2; ++j) {
		if (j == i__) {
		    goto L25;
		}
		i__4 = j;
		z__2.r = z__.r - zo[i__4].r, z__2.i = z__.i - zo[i__4].i;
		z__1.r = zw.r * z__2.r - zw.i * z__2.i, z__1.i = zw.r * 
			z__2.i + zw.i * z__2.r;
		zw.r = z__1.r, zw.i = z__1.i;
L25:
		;
	    }
/* L30: */
	    z__1.r = zq.r + zw.r, z__1.i = zq.i + zw.i;
	    zq.r = z__1.r, zq.i = z__1.i;
	}
	z__3.r = zq.r * zfd.r - zq.i * zfd.i, z__3.i = zq.r * zfd.i + zq.i * 
		zfd.r;
	z__2.r = zd.r - z__3.r, z__2.i = zd.i - z__3.i;
	z_div(&z__1, &z__2, &zp);
	zgd.r = z__1.r, zgd.i = z__1.i;
	z_div(&z__2, &zfd, &zgd);
	z__1.r = z__.r - z__2.r, z__1.i = z__.i - z__2.i;
	z__.r = z__1.r, z__.i = z__1.i;
	w0 = w;
	w = z_abs(&z__);
	if (it <= 50 && (d__1 = (w - w0) / w, abs(d__1)) > 1e-11) {
	    goto L15;
	}
/* L35: */
	i__3 = nr;
	zo[i__3].r = z__.r, zo[i__3].i = z__.i;
    }
    return 0;
} /* cerzo_ */

/*       ********************************** */
/* Subroutine */ int gamma2_(doublereal *x, doublereal *ga)
{
    /* Initialized data */

    static doublereal g[26] = { 1.,.5772156649015329,-.6558780715202538,
	    -.0420026350340952,.1665386113822915,-.0421977345555443,
	    -.009621971527877,.007218943246663,-.0011651675918591,
	    -2.152416741149e-4,1.280502823882e-4,-2.01348547807e-5,
	    -1.2504934821e-6,1.133027232e-6,-2.056338417e-7,6.116095e-9,
	    5.0020075e-9,-1.1812746e-9,1.043427e-10,7.7823e-12,-3.6968e-12,
	    5.1e-13,-2.06e-14,-5.4e-15,1.4e-15,1e-16 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static integer k, m;
    static doublereal r__, z__;
    static integer m1;
    static doublereal pi, gr;


/*       ================================================== */
/*       Purpose: Compute gamma function Г(x) */
/*       Input :  x  --- Argument of Г(x) */
/*                       ( x is not equal to 0,-1,-2,…) */
/*       Output:  GA --- Г(x) */
/*       ================================================== */

    pi = 3.141592653589793;
    if (*x == (doublereal) ((integer) (*x))) {
	if (*x > 0.) {
	    *ga = 1.;
	    m1 = (integer) (*x - 1);
	    i__1 = m1;
	    for (k = 2; k <= i__1; ++k) {
/* L10: */
		*ga *= k;
	    }
	} else {
	    *ga = 1e300;
	}
    } else {
	r__ = 1.;
	if (abs(*x) > 1.) {
	    z__ = abs(*x);
	    m = (integer) z__;
	    i__1 = m;
	    for (k = 1; k <= i__1; ++k) {
/* L15: */
		r__ *= z__ - k;
	    }
	    z__ -= m;
	} else {
	    z__ = *x;
	}
	gr = g[25];
	for (k = 25; k >= 1; --k) {
/* L20: */
	    gr = gr * z__ + g[k - 1];
	}
	*ga = 1. / (gr * z__);
	if (abs(*x) > 1.) {
	    *ga *= r__;
	    if (*x < 0.) {
		*ga = -pi / (*x * *ga * sin(pi * *x));
	    }
	}
    }
    return 0;
} /* gamma2_ */

/*       ********************************** */
/* Subroutine */ int chgu_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *hu, integer *md, integer *isfer)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal a00, aa, b00;
    static integer id;
    static logical bn;
    static integer id1;
    static logical bl1, bl2, bl3, il1, il2, il3;
    static doublereal hu1;
    extern /* Subroutine */ int chgul_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *), chgus_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), chgubi_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), chguit_(
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;


/*       ======================================================= */
/*       Purpose: Compute the confluent hypergeometric function */
/*                U(a,b,x) */
/*       Input  : a  --- Parameter */
/*                b  --- Parameter */
/*                x  --- Argument  ( x > 0 ) */
/*       Output:  HU --- U(a,b,x) */
/*                MD --- Method code */
/*                ISFER --- Error flag */
/*       Routines called: */
/*            (1) CHGUS for small x ( MD=1 ) */
/*            (2) CHGUL for large x ( MD=2 ) */
/*            (3) CHGUBI for integer b ( MD=3 ) */
/*            (4) CHGUIT for numerical integration ( MD=4 ) */
/*       ======================================================= */

    aa = *a - *b + 1.;
    *isfer = 0;
    il1 = *a == (doublereal) ((integer) (*a)) && *a <= 0.f;
    il2 = aa == (doublereal) ((integer) aa) && aa <= 0.f;
    il3 = (d__1 = *a * (*a - *b + 1.f), abs(d__1)) / *x <= 2.f;
    bl1 = *x <= 5.f || *x <= 10.f && *a <= 2.f;
    bl2 = *x > 5.f && *x <= 12.5f && (*a >= 1.f && *b >= *a + 4.f);
    bl3 = *x > 12.5f && *a >= 5.f && *b >= *a + 5.f;
    bn = *b == (doublereal) ((integer) (*b)) && *b != 0.f;
    id1 = -100;
    hu1 = 0.;
    if (*b != (doublereal) ((integer) (*b))) {
	chgus_(a, b, x, hu, &id1);
	*md = 1;
	if (id1 >= 9) {
	    return 0;
	}
	hu1 = *hu;
    }
    if (il1 || il2 || il3) {
	chgul_(a, b, x, hu, &id);
	*md = 2;
	if (id >= 9) {
	    return 0;
	}
	if (id1 > id) {
	    *md = 1;
	    id = id1;
	    *hu = hu1;
	}
    }
    if (*a >= 1.f) {
	if (bn && (bl1 || bl2 || bl3)) {
	    chgubi_(a, b, x, hu, &id);
	    *md = 3;
	} else {
	    chguit_(a, b, x, hu, &id);
	    *md = 4;
	}
    } else {
	if (*b <= *a) {
	    a00 = *a;
	    b00 = *b;
	    *a = *a - *b + 1.;
	    *b = 2. - *b;
	    chguit_(a, b, x, hu, &id);
	    d__1 = 1. - b00;
	    *hu = pow_dd(x, &d__1) * *hu;
	    *a = a00;
	    *b = b00;
	    *md = 4;
	} else if (bn && ! il1) {
	    chgubi_(a, b, x, hu, &id);
	    *md = 3;
	}
    }
    if (id < 6) {
	*isfer = 6;
    }
    return 0;
} /* chgu_ */

/*       ********************************** */
/* Subroutine */ int lamn_(integer *n, doublereal *x, integer *nm, doublereal 
	*bl, doublereal *dl)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal f;
    static integer i__, k, m;
    static doublereal r__, f0, f1, r0, x2, bg, bk, bs, uk;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ========================================================= */
/*       Purpose: Compute lambda functions and their derivatives */
/*       Input:   x --- Argument of lambda function */
/*                n --- Order of lambda function */
/*       Output:  BL(n) --- Lambda function of order n */
/*                DL(n) --- Derivative of lambda function */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 for computing the start */
/*                point for backward recurrence */
/*       ========================================================= */

    *nm = *n;
    if (abs(*x) < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    bl[k] = 0.;
/* L10: */
	    dl[k] = 0.;
	}
	bl[0] = 1.;
	dl[1] = .5;
	return 0;
    }
    if (*x <= 12.) {
	x2 = *x * *x;
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    bk = 1.;
	    r__ = 1.;
	    for (i__ = 1; i__ <= 50; ++i__) {
		r__ = r__ * -.25 * x2 / (i__ * (i__ + k));
		bk += r__;
		if (abs(r__) < abs(bk) * 1e-15) {
		    goto L20;
		}
/* L15: */
	    }
L20:
	    bl[k] = bk;
/* L25: */
	    if (k >= 1) {
		dl[k - 1] = *x * -.5 / k * bk;
	    }
	}
	uk = 1.;
	r__ = 1.;
	for (i__ = 1; i__ <= 50; ++i__) {
	    r__ = r__ * -.25 * x2 / (i__ * (i__ + *n + 1.));
	    uk += r__;
	    if (abs(r__) < abs(uk) * 1e-15) {
		goto L35;
	    }
/* L30: */
	}
L35:
	dl[*n] = *x * -.5 / (*n + 1.) * uk;
	return 0;
    }
    if (*n == 0) {
	*nm = 1;
    }
    m = msta1_(x, &c__200);
    if (m < *nm) {
	*nm = m;
    } else {
	m = msta2_(x, nm, &c__15);
    }
    bs = 0.;
    f = 0.;
    f0 = 0.;
    f1 = 1e-100;
    for (k = m; k >= 0; --k) {
	f = (k + 1.) * 2. * f1 / *x - f0;
	if (k <= *nm) {
	    bl[k] = f;
	}
	if (k == k / 2 << 1) {
	    bs += f * 2.;
	}
	f0 = f1;
/* L40: */
	f1 = f;
    }
    bg = bs - f;
    i__1 = *nm;
    for (k = 0; k <= i__1; ++k) {
/* L45: */
	bl[k] /= bg;
    }
    r0 = 1.;
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
	r0 = r0 * 2. * k / *x;
/* L50: */
	bl[k] = r0 * bl[k];
    }
    dl[0] = *x * -.5 * bl[1];
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L55: */
	dl[k] = k * 2. / *x * (bl[k - 1] - bl[k]);
    }
    return 0;
} /* lamn_ */

/*       ********************************** */
/* Subroutine */ int comelp_(doublereal *hk, doublereal *ck, doublereal *ce)
{
    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal ae, be, ak, bk, pk;


/*       ================================================== */
/*       Purpose: Compute complete elliptic integrals K(k) */
/*                and E(k) */
/*       Input  : K  --- Modulus k ( 0 ≤ k ≤ 1 ) */
/*       Output : CK --- K(k) */
/*                CE --- E(k) */
/*       ================================================== */

    pk = 1. - *hk * *hk;
    if (*hk == 1.f) {
	*ck = 1e300;
	*ce = 1.;
    } else {
	ak = (((pk * .01451196212 + .03742563713) * pk + .03590092383) * pk + 
		.09666344259) * pk + 1.38629436112;
	bk = (((pk * .00441787012 + .03328355346) * pk + .06880248576) * pk + 
		.12498593597) * pk + .5;
	*ck = ak - bk * log(pk);
	ae = (((pk * .01736506451 + .04757383546) * pk + .0626060122) * pk + 
		.44325141463) * pk + 1.;
	be = (((pk * .00526449639 + .04069697526) * pk + .09200180037) * pk + 
		.2499836831) * pk;
	*ce = ae - be * log(pk);
    }
    return 0;
} /* comelp_ */

/*       ********************************** */
/* Subroutine */ int incob_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *bix)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer k;
    static doublereal s0, t1, t2, dk[51], fk[51], ta, bt, tb;
    extern /* Subroutine */ int beta_(doublereal *, doublereal *, doublereal *
	    );


/*       ======================================================== */
/*       Purpose: Compute the incomplete beta function Ix(a,b) */
/*       Input :  a --- Parameter */
/*                b --- Parameter */
/*                x --- Argument ( 0 ≤ x ≤ 1 ) */
/*       Output:  BIX --- Ix(a,b) */
/*       Routine called: BETA for computing beta function B(p,q) */
/*       ======================================================== */

    s0 = (*a + 1.) / (*a + *b + 2.);
    beta_(a, b, &bt);
    if (*x <= s0) {
	for (k = 1; k <= 20; ++k) {
/* L10: */
	    dk[(k << 1) - 1] = k * (*b - k) * *x / (*a + k * 2. - 1.) / (*a + 
		    k * 2.);
	}
	for (k = 0; k <= 20; ++k) {
/* L15: */
	    dk[k * 2] = -(*a + k) * (*a + *b + k) * *x / (*a + k * 2.) / (*a 
		    + k * 2.f + 1.f);
	}
	t1 = 0.;
	for (k = 20; k >= 1; --k) {
/* L20: */
	    t1 = dk[k - 1] / (t1 + 1.);
	}
	ta = 1. / (t1 + 1.);
	d__1 = 1. - *x;
	*bix = pow_dd(x, a) * pow_dd(&d__1, b) / (*a * bt) * ta;
    } else {
	for (k = 1; k <= 20; ++k) {
/* L25: */
	    fk[(k << 1) - 1] = k * (*a - k) * (1. - *x) / (*b + k * 2.f - 1.f)
		     / (*b + k * 2.f);
	}
	for (k = 0; k <= 20; ++k) {
/* L30: */
	    fk[k * 2] = -(*b + k) * (*a + *b + k) * (1. - *x) / (*b + k * 2.) 
		    / (*b + k * 2. + 1.);
	}
	t2 = 0.;
	for (k = 20; k >= 1; --k) {
/* L35: */
	    t2 = fk[k - 1] / (t2 + 1.);
	}
	tb = 1. / (t2 + 1.);
	d__1 = 1. - *x;
	*bix = 1. - pow_dd(x, a) * pow_dd(&d__1, b) / (*b * bt) * tb;
    }
    return 0;
} /* incob_ */

/*       ********************************** */
/* Subroutine */ int cvf_(integer *kd, integer *m, doublereal *q, doublereal *
	a, integer *mj, doublereal *f)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal b;
    static integer j, l, j0, l0;
    static doublereal t0, t1, t2;
    static integer ic, jf;


/*       ====================================================== */
/*       Purpose: Compute the value of F for characteristic */
/*                equation of Mathieu functions */
/*       Input :  m --- Order of Mathieu functions */
/*                q --- Parameter of Mathieu functions */
/*                A --- Characteristic value */
/*       Output:  F --- Value of F for characteristic equation */
/*       ====================================================== */

    b = *a;
    ic = *m / 2;
    l = 0;
    l0 = 0;
    j0 = 2;
    jf = ic;
    if (*kd == 1) {
	l0 = 2;
    }
    if (*kd == 1) {
	j0 = 3;
    }
    if (*kd == 2 || *kd == 3) {
	l = 1;
    }
    if (*kd == 4) {
	jf = ic - 1;
    }
    t1 = 0.;
    i__1 = ic + 1;
    for (j = *mj; j >= i__1; --j) {
/* L10: */
/* Computing 2nd power */
	d__1 = j * 2. + l;
	t1 = -(*q) * *q / (d__1 * d__1 - b + t1);
    }
    if (*m <= 2) {
	t2 = 0.;
	if (*kd == 1 && *m == 0) {
	    t1 += t1;
	}
	if (*kd == 1 && *m == 2) {
	    t1 = *q * -2. * *q / (4. - b + t1) - 4.;
	}
	if (*kd == 2 && *m == 1) {
	    t1 += *q;
	}
	if (*kd == 3 && *m == 1) {
	    t1 -= *q;
	}
    } else {
	t0 = 0.;
	if (*kd == 1) {
	    t0 = 4. - b + *q * 2. * *q / b;
	}
	if (*kd == 2) {
	    t0 = 1. - b + *q;
	}
	if (*kd == 3) {
	    t0 = 1. - b - *q;
	}
	if (*kd == 4) {
	    t0 = 4. - b;
	}
	t2 = -(*q) * *q / t0;
	i__1 = jf;
	for (j = j0; j <= i__1; ++j) {
/* L15: */
/* Computing 2nd power */
	    d__1 = j * 2. - l - l0;
	    t2 = -(*q) * *q / (d__1 * d__1 - b + t2);
	}
    }
/* Computing 2nd power */
    d__1 = ic * 2. + l;
    *f = d__1 * d__1 + t1 + t2 - b;
    return 0;
} /* cvf_ */

/*       ********************************** */
/* Subroutine */ int clpn_(integer *n, doublereal *x, doublereal *y, 
	doublecomplex *cpn, doublecomplex *cpd)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex z__, cp0, cp1, cpf;


/*       ================================================== */
/*       Purpose: Compute Legendre polynomials Pn(z) and */
/*                their derivatives Pn'(z) for a complex */
/*                argument */
/*       Input :  x --- Real part of z */
/*                y --- Imaginary part of z */
/*                n --- Degree of Pn(z), n = 0,1,2,... */
/*       Output:  CPN(n) --- Pn(z) */
/*                CPD(n) --- Pn'(z) */
/*       ================================================== */

    z__1.r = *x, z__1.i = *y;
    z__.r = z__1.r, z__.i = z__1.i;
    cpn[0].r = 1., cpn[0].i = 0.;
    cpn[1].r = z__.r, cpn[1].i = z__.i;
    cpd[0].r = 0., cpd[0].i = 0.;
    cpd[1].r = 1., cpd[1].i = 0.;
    cp0.r = 1., cp0.i = 0.;
    cp1.r = z__.r, cp1.i = z__.i;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	d__1 = (k * 2. - 1.) / k;
	z__3.r = d__1 * z__.r, z__3.i = d__1 * z__.i;
	z__2.r = z__3.r * cp1.r - z__3.i * cp1.i, z__2.i = z__3.r * cp1.i + 
		z__3.i * cp1.r;
	d__2 = (k - 1.) / k;
	z__4.r = d__2 * cp0.r, z__4.i = d__2 * cp0.i;
	z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
	cpf.r = z__1.r, cpf.i = z__1.i;
	i__2 = k;
	cpn[i__2].r = cpf.r, cpn[i__2].i = cpf.i;
	if (abs(*x) == 1. && *y == 0.) {
	    i__2 = k;
	    i__3 = k + 1;
	    d__1 = pow_di(x, &i__3) * .5 * k * (k + 1.);
	    cpd[i__2].r = d__1, cpd[i__2].i = 0.;
	} else {
	    i__2 = k;
	    z__4.r = z__.r * cpf.r - z__.i * cpf.i, z__4.i = z__.r * cpf.i + 
		    z__.i * cpf.r;
	    z__3.r = cp1.r - z__4.r, z__3.i = cp1.i - z__4.i;
	    d__1 = (doublereal) k;
	    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	    z__6.r = z__.r * z__.r - z__.i * z__.i, z__6.i = z__.r * z__.i + 
		    z__.i * z__.r;
	    z__5.r = 1. - z__6.r, z__5.i = -z__6.i;
	    z_div(&z__1, &z__2, &z__5);
	    cpd[i__2].r = z__1.r, cpd[i__2].i = z__1.i;
	}
	cp0.r = cp1.r, cp0.i = cp1.i;
/* L10: */
	cp1.r = cpf.r, cp1.i = cpf.i;
    }
    return 0;
} /* clpn_ */

/*       ********************************** */
/* Subroutine */ int lqmns_(integer *m, integer *n, doublereal *x, doublereal 
	*qm, doublereal *qd)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer k, l;
    static doublereal q0, q00, q10, q01, q11;
    static integer km, ls;
    static doublereal xq, qf0, qf1, qf2, qg0, qg1, qh0, qh1, qm0, qm1, qh2, 
	    q0l, q1l, qmk;


/*       ======================================================== */
/*       Purpose: Compute associated Legendre functions Qmn(x) */
/*                and Qmn'(x) for a given order */
/*       Input :  x --- Argument of Qmn(x) */
/*                m --- Order of Qmn(x),  m = 0,1,2,... */
/*                n --- Degree of Qmn(x), n = 0,1,2,... */
/*       Output:  QM(n) --- Qmn(x) */
/*                QD(n) --- Qmn'(x) */
/*       ======================================================== */

    i__1 = *n;
    for (k = 0; k <= i__1; ++k) {
	qm[k] = 0.;
/* L10: */
	qd[k] = 0.;
    }
    if (abs(*x) == 1.) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    qm[k] = 1e300;
/* L15: */
	    qd[k] = 1e300;
	}
	return 0;
    }
    ls = 1;
    if (abs(*x) > 1.) {
	ls = -1;
    }
    xq = sqrt(ls * (1. - *x * *x));
    q0 = log((d__1 = (*x + 1.f) / (*x - 1.f), abs(d__1))) * .5;
    q00 = q0;
    q10 = -1. / xq;
    q01 = *x * q0 - 1.;
    q11 = -ls * xq * (q0 + *x / (1. - *x * *x));
    qf0 = q00;
    qf1 = q10;
    qm0 = 0.;
    qm1 = 0.;
    i__1 = *m;
    for (k = 2; k <= i__1; ++k) {
	qm0 = (k - 1.f) * -2. / xq * *x * qf1 - ls * (k - 1.f) * (2.f - k) * 
		qf0;
	qf0 = qf1;
/* L20: */
	qf1 = qm0;
    }
    if (*m == 0) {
	qm0 = q00;
    }
    if (*m == 1) {
	qm0 = q10;
    }
    qm[0] = qm0;
    if (abs(*x) < 1.0001) {
	if (*m == 0 && *n > 0) {
	    qf0 = q00;
	    qf1 = q01;
	    i__1 = *n;
	    for (k = 2; k <= i__1; ++k) {
		qf2 = ((k * 2.f - 1.) * *x * qf1 - (k - 1.f) * qf0) / k;
		qm[k] = qf2;
		qf0 = qf1;
/* L25: */
		qf1 = qf2;
	    }
	}
	qg0 = q01;
	qg1 = q11;
	i__1 = *m;
	for (k = 2; k <= i__1; ++k) {
	    qm1 = (k - 1.f) * -2. / xq * *x * qg1 - ls * k * (3.f - k) * qg0;
	    qg0 = qg1;
/* L30: */
	    qg1 = qm1;
	}
	if (*m == 0) {
	    qm1 = q01;
	}
	if (*m == 1) {
	    qm1 = q11;
	}
	qm[1] = qm1;
	if (*m == 1 && *n > 1) {
	    qh0 = q10;
	    qh1 = q11;
	    i__1 = *n;
	    for (k = 2; k <= i__1; ++k) {
		qh2 = ((k * 2.f - 1.) * *x * qh1 - k * qh0) / (k - 1.f);
		qm[k] = qh2;
		qh0 = qh1;
/* L35: */
		qh1 = qh2;
	    }
	} else if (*m >= 2) {
	    qg0 = q00;
	    qg1 = q01;
	    qh0 = q10;
	    qh1 = q11;
	    qmk = 0.;
	    i__1 = *n;
	    for (l = 2; l <= i__1; ++l) {
		q0l = ((l * 2. - 1.) * *x * qg1 - (l - 1.) * qg0) / l;
		q1l = ((l * 2.f - 1.) * *x * qh1 - l * qh0) / (l - 1.);
		qf0 = q0l;
		qf1 = q1l;
		i__2 = *m;
		for (k = 2; k <= i__2; ++k) {
		    qmk = (k - 1.f) * -2. / xq * *x * qf1 - ls * (k + l - 1.f)
			     * (l + 2.f - k) * qf0;
		    qf0 = qf1;
/* L40: */
		    qf1 = qmk;
		}
		qm[l] = qmk;
		qg0 = qg1;
		qg1 = q0l;
		qh0 = qh1;
/* L45: */
		qh1 = q1l;
	    }
	}
    } else {
	if (abs(*x) > 1.1f) {
	    km = *m + 40 + *n;
	} else {
	    km = (*m + 40 + *n) * (integer) (-1.f - log(*x - 1.f) * 1.8f);
	}
	qf2 = 0.;
	qf1 = 1.;
	for (k = km; k >= 0; --k) {
	    qf0 = ((k * 2.f + 3.) * *x * qf1 - (k + 2.f - *m) * qf2) / (k + *
		    m + 1.f);
	    if (k <= *n) {
		qm[k] = qf0;
	    }
	    qf2 = qf1;
/* L50: */
	    qf1 = qf0;
	}
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
/* L55: */
	    qm[k] = qm[k] * qm0 / qf0;
	}
    }
    if (abs(*x) < 1.) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
/* L60: */
	    qm[k] = pow_ii(&c_n1, m) * qm[k];
	}
    }
    qd[0] = ((1. - *m) * qm[1] - *x * qm[0]) / (*x * *x - 1.f);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* L65: */
	qd[k] = (k * *x * qm[k] - (k + *m) * qm[k - 1]) / (*x * *x - 1.f);
    }
    return 0;
} /* lqmns_ */

/*       ********************************** */
/* Subroutine */ int ciklv_(doublereal *v, doublecomplex *z__, doublecomplex *
	cbiv, doublecomplex *cdiv, doublecomplex *cbkv, doublecomplex *cdkv)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    void z_sqrt(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *), z_log(doublecomplex *, 
	    doublecomplex *), pow_zi(doublecomplex *, doublecomplex *, 
	    integer *);
    double pow_di(doublereal *, integer *);
    void z_exp(doublecomplex *, doublecomplex *);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal a[91];
    static integer i__, k, l, l0;
    static doublereal v0;
    static doublecomplex cf[12];
    static integer lf;
    static doublecomplex ct;
    static integer km;
    static doublereal pi, vr;
    static doublecomplex ct2, cfi, cfk;
    extern /* Subroutine */ int cjk_(integer *, doublereal *);
    static doublecomplex csi, csk, cws, ceta;


/*       ===================================================== */
/*       Purpose: Compute modified Bessel functions Iv(z) and */
/*                Kv(z) and their derivatives with a complex */
/*                argument and a large order */
/*       Input:   v --- Order of Iv(z) and Kv(z) */
/*                z --- Complex argument */
/*       Output:  CBIV --- Iv(z) */
/*                CDIV --- Iv'(z) */
/*                CBKV --- Kv(z) */
/*                CDKV --- Kv'(z) */
/*       Routine called: */
/*                CJK to compute the expansion coefficients */
/*       ==================================================== */

    pi = 3.141592653589793;
    km = 12;
    cjk_(&km, a);
    for (l = 1; l >= 0; --l) {
	v0 = *v - l;
	z__4.r = z__->r / v0, z__4.i = z__->i / v0;
	z__5.r = z__->r / v0, z__5.i = z__->i / v0;
	z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, z__3.i = z__4.r * z__5.i 
		+ z__4.i * z__5.r;
	z__2.r = z__3.r + 1., z__2.i = z__3.i;
	z_sqrt(&z__1, &z__2);
	cws.r = z__1.r, cws.i = z__1.i;
	z__4.r = z__->r / v0, z__4.i = z__->i / v0;
	z__5.r = cws.r + 1., z__5.i = cws.i;
	z_div(&z__3, &z__4, &z__5);
	z_log(&z__2, &z__3);
	z__1.r = cws.r + z__2.r, z__1.i = cws.i + z__2.i;
	ceta.r = z__1.r, ceta.i = z__1.i;
	z_div(&z__1, &c_b6, &cws);
	ct.r = z__1.r, ct.i = z__1.i;
	z__1.r = ct.r * ct.r - ct.i * ct.i, z__1.i = ct.r * ct.i + ct.i * 
		ct.r;
	ct2.r = z__1.r, ct2.i = z__1.i;
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
	    l0 = k * (k + 1) / 2 + 1;
	    lf = l0 + k;
	    i__2 = k - 1;
	    i__3 = lf - 1;
	    cf[i__2].r = a[i__3], cf[i__2].i = 0.;
	    i__2 = l0;
	    for (i__ = lf - 1; i__ >= i__2; --i__) {
/* L10: */
		i__3 = k - 1;
		i__4 = k - 1;
		z__2.r = cf[i__4].r * ct2.r - cf[i__4].i * ct2.i, z__2.i = cf[
			i__4].r * ct2.i + cf[i__4].i * ct2.r;
		i__5 = i__ - 1;
		z__1.r = z__2.r + a[i__5], z__1.i = z__2.i;
		cf[i__3].r = z__1.r, cf[i__3].i = z__1.i;
	    }
/* L15: */
	    i__3 = k - 1;
	    i__4 = k - 1;
	    pow_zi(&z__2, &ct, &k);
	    z__1.r = cf[i__4].r * z__2.r - cf[i__4].i * z__2.i, z__1.i = cf[
		    i__4].r * z__2.i + cf[i__4].i * z__2.r;
	    cf[i__3].r = z__1.r, cf[i__3].i = z__1.i;
	}
	vr = 1. / v0;
	csi.r = 1., csi.i = 0.;
	i__3 = km;
	for (k = 1; k <= i__3; ++k) {
/* L20: */
	    i__4 = k - 1;
	    d__1 = pow_di(&vr, &k);
	    z__2.r = d__1 * cf[i__4].r, z__2.i = d__1 * cf[i__4].i;
	    z__1.r = csi.r + z__2.r, z__1.i = csi.i + z__2.i;
	    csi.r = z__1.r, csi.i = z__1.i;
	}
	d__1 = pi * 2. * v0;
	z__4.r = ct.r / d__1, z__4.i = ct.i / d__1;
	z_sqrt(&z__3, &z__4);
	z__6.r = v0 * ceta.r, z__6.i = v0 * ceta.i;
	z_exp(&z__5, &z__6);
	z__2.r = z__3.r * z__5.r - z__3.i * z__5.i, z__2.i = z__3.r * z__5.i 
		+ z__3.i * z__5.r;
	z__1.r = z__2.r * csi.r - z__2.i * csi.i, z__1.i = z__2.r * csi.i + 
		z__2.i * csi.r;
	cbiv->r = z__1.r, cbiv->i = z__1.i;
	if (l == 1) {
	    cfi.r = cbiv->r, cfi.i = cbiv->i;
	}
	csk.r = 1., csk.i = 0.;
	i__4 = km;
	for (k = 1; k <= i__4; ++k) {
/* L25: */
	    i__3 = pow_ii(&c_n1, &k);
	    i__1 = k - 1;
	    d__1 = (doublereal) i__3;
	    z__3.r = d__1 * cf[i__1].r, z__3.i = d__1 * cf[i__1].i;
	    d__2 = pow_di(&vr, &k);
	    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
	    z__1.r = csk.r + z__2.r, z__1.i = csk.i + z__2.i;
	    csk.r = z__1.r, csk.i = z__1.i;
	}
	z__5.r = pi * ct.r, z__5.i = pi * ct.i;
	d__1 = v0 * 2.;
	z__4.r = z__5.r / d__1, z__4.i = z__5.i / d__1;
	z_sqrt(&z__3, &z__4);
	d__2 = -v0;
	z__7.r = d__2 * ceta.r, z__7.i = d__2 * ceta.i;
	z_exp(&z__6, &z__7);
	z__2.r = z__3.r * z__6.r - z__3.i * z__6.i, z__2.i = z__3.r * z__6.i 
		+ z__3.i * z__6.r;
	z__1.r = z__2.r * csk.r - z__2.i * csk.i, z__1.i = z__2.r * csk.i + 
		z__2.i * csk.r;
	cbkv->r = z__1.r, cbkv->i = z__1.i;
	if (l == 1) {
	    cfk.r = cbkv->r, cfk.i = cbkv->i;
	}
/* L30: */
    }
    z__4.r = *v, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbiv->r - z__3.i * cbiv->i, z__2.i = z__3.r * cbiv->i + 
	    z__3.i * cbiv->r;
    z__1.r = cfi.r - z__2.r, z__1.i = cfi.i - z__2.i;
    cdiv->r = z__1.r, cdiv->i = z__1.i;
    z__2.r = -cfk.r, z__2.i = -cfk.i;
    z__5.r = *v, z__5.i = 0.;
    z_div(&z__4, &z__5, z__);
    z__3.r = z__4.r * cbkv->r - z__4.i * cbkv->i, z__3.i = z__4.r * cbkv->i + 
	    z__4.i * cbkv->r;
    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
    cdkv->r = z__1.r, cdkv->i = z__1.i;
    return 0;
} /* ciklv_ */

/*       ********************************** */
/* Subroutine */ int elit_(doublereal *hk, doublereal *phi, doublereal *fe, 
	doublereal *ee)
{
    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal), cos(doublereal), log(doublereal)
	    , tan(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal a, b, c__, d__, g;
    static integer n;
    static doublereal r__, a0, b0, d0, ce, ck, pi, fac;


/*       ================================================== */
/*       Purpose: Compute complete and incomplete elliptic */
/*                integrals F(k,phi) and E(k,phi) */
/*       Input  : HK  --- Modulus k ( 0 ≤ k ≤ 1 ) */
/*                Phi --- Argument ( in degrees ) */
/*       Output : FE  --- F(k,phi) */
/*                EE  --- E(k,phi) */
/*       ================================================== */

    g = 0.;
    pi = 3.14159265358979;
    a0 = 1.;
    b0 = sqrt(1. - *hk * *hk);
    d0 = pi / 180. * *phi;
    r__ = *hk * *hk;
    if (*hk == 1. && *phi == 90.) {
	*fe = 1e300;
	*ee = 1.;
    } else if (*hk == 1.) {
	*fe = log((sin(d0) + 1.) / cos(d0));
	*ee = sin(d0);
    } else {
	fac = 1.;
	d__ = 0.;
	for (n = 1; n <= 40; ++n) {
	    a = (a0 + b0) / 2.;
	    b = sqrt(a0 * b0);
	    c__ = (a0 - b0) / 2.;
	    fac *= 2.;
	    r__ += fac * c__ * c__;
	    if (*phi != 90.) {
		d__ = d0 + atan(b0 / a0 * tan(d0));
		g += c__ * sin(d__);
		d0 = d__ + pi * (integer) (d__ / pi + .5);
	    }
	    a0 = a;
	    b0 = b;
	    if (c__ < 1e-7) {
		goto L15;
	    }
/* L10: */
	}
L15:
	ck = pi / (a * 2.);
	ce = pi * (2. - r__) / (a * 4.);
	if (*phi == 90.) {
	    *fe = ck;
	    *ee = ce;
	} else {
	    *fe = d__ / (fac * a);
	    *ee = *fe * ce / ck + g;
	}
    }
    return 0;
} /* elit_ */

/*       ********************************** */
/* Subroutine */ int elit3_(doublereal *phi, doublereal *hk, doublereal *c__, 
	doublereal *el3)
{
    /* Initialized data */

    static doublereal t[10] = { .9931285991850949,.9639719272779138,
	    .9122344282513259,.8391169718222188,.7463319064601508,
	    .636053680726515,.5108670019508271,.3737060887154195,
	    .2277858511416451,.07652652113349734 };
    static doublereal w[10] = { .01761400713915212,.04060142980038694,
	    .06267204833410907,.08327674157670475,.1019301198172404,
	    .1181945319615184,.1316886384491766,.142096109318382,
	    .1491729864726037,.1527533871307258 };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal c0, c1, c2, f1, f2, t1, t2;
    static logical lb1, lb2;


/*       ========================================================= */
/*       Purpose: Compute the elliptic integral of the third kind */
/*                using Gauss-Legendre quadrature */
/*       Input :  Phi --- Argument ( in degrees ) */
/*                 k  --- Modulus   ( 0 ≤ k ≤ 1.0 ) */
/*                 c  --- Parameter ( 0 ≤ c ≤ 1.0 ) */
/*       Output:  EL3 --- Value of the elliptic integral of the */
/*                        third kind */
/*       ========================================================= */

    lb1 = *hk == 1. && (d__1 = *phi - 90.f, abs(d__1)) <= 1e-8;
    lb2 = *c__ == 1. && (d__1 = *phi - 90.f, abs(d__1)) <= 1e-8;
    if (lb1 || lb2) {
	*el3 = 1e300;
	return 0;
    }
    c1 = *phi * .0087266462599716;
    c2 = c1;
    *el3 = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	c0 = c2 * t[i__ - 1];
	t1 = c1 + c0;
	t2 = c1 - c0;
	f1 = 1. / ((1. - *c__ * sin(t1) * sin(t1)) * sqrt(1. - *hk * *hk * 
		sin(t1) * sin(t1)));
	f2 = 1. / ((1. - *c__ * sin(t2) * sin(t2)) * sqrt(1. - *hk * *hk * 
		sin(t2) * sin(t2)));
/* L10: */
	*el3 += w[i__ - 1] * (f1 + f2);
    }
    *el3 = c1 * *el3;
    return 0;
} /* elit3_ */

/*       ********************************** */
/* Subroutine */ int eix_(doublereal *x, doublereal *ei)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, ga;
    extern /* Subroutine */ int e1xb_(doublereal *, doublereal *);


/*       ============================================ */
/*       Purpose: Compute exponential integral Ei(x) */
/*       Input :  x  --- Argument of Ei(x) */
/*       Output:  EI --- Ei(x) */
/*       ============================================ */

    if (*x == 0.f) {
	*ei = -1e300;
    } else if (*x < 0.) {
	d__1 = -(*x);
	e1xb_(&d__1, ei);
	*ei = -(*ei);
    } else if (abs(*x) <= 40.f) {
/*          Power series around x=0 */
	*ei = 1.;
	r__ = 1.;
	for (k = 1; k <= 100; ++k) {
/* Computing 2nd power */
	    d__1 = k + 1.;
	    r__ = r__ * k * *x / (d__1 * d__1);
	    *ei += r__;
	    if ((d__1 = r__ / *ei, abs(d__1)) <= 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	ga = .5772156649015328;
	*ei = ga + log(*x) + *x * *ei;
    } else {
/*          Asymptotic expansion (the series is not convergent) */
	*ei = 1.;
	r__ = 1.;
	for (k = 1; k <= 20; ++k) {
	    r__ = r__ * k / *x;
/* L25: */
	    *ei += r__;
	}
	*ei = exp(*x) / *x * *ei;
    }
    return 0;
} /* eix_ */

/*       ********************************** */
/* Subroutine */ int eixz_(doublecomplex *z__, doublecomplex *cei)
{
    /* System generated locals */
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal pi;
    extern /* Subroutine */ int e1z_(doublecomplex *, doublecomplex *);


/*       ============================================ */
/*       Purpose: Compute exponential integral Ei(x) */
/*       Input :  x  --- Complex argument of Ei(x) */
/*       Output:  EI --- Ei(x) */
/*       ============================================ */

    pi = 3.141592653589793;
    z__1.r = -z__->r, z__1.i = -z__->i;
    e1z_(&z__1, cei);
    z__1.r = -cei->r, z__1.i = -cei->i;
    cei->r = z__1.r, cei->i = z__1.i;
    if (d_imag(z__) > 0.) {
	z__2.r = pi * 0., z__2.i = pi * 1.;
	z__1.r = cei->r + z__2.r, z__1.i = cei->i + z__2.i;
	cei->r = z__1.r, cei->i = z__1.i;
    } else if (d_imag(z__) < 0.) {
	z__2.r = pi * 0., z__2.i = pi * 1.;
	z__1.r = cei->r - z__2.r, z__1.i = cei->i - z__2.i;
	cei->r = z__1.r, cei->i = z__1.i;
    } else if (d_imag(z__) == 0.) {
	if (z__->r > 0.) {
	    z__2.r = pi * 0., z__2.i = pi * 1.;
	    z__1.r = cei->r - z__2.r, z__1.i = cei->i - z__2.i;
	    cei->r = z__1.r, cei->i = z__1.i;
	}
    }
    return 0;
} /* eixz_ */

/*       ********************************** */
/* Subroutine */ int e1xb_(doublereal *x, doublereal *e1)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer k, m;
    static doublereal r__, t, t0, ga;


/*       ============================================ */
/*       Purpose: Compute exponential integral E1(x) */
/*       Input :  x  --- Argument of E1(x) */
/*       Output:  E1 --- E1(x)  ( x > 0 ) */
/*       ============================================ */

    if (*x == 0.f) {
	*e1 = 1e300;
    } else if (*x <= 1.f) {
	*e1 = 1.;
	r__ = 1.;
	for (k = 1; k <= 25; ++k) {
/* Computing 2nd power */
	    d__1 = k + 1.;
	    r__ = -r__ * k * *x / (d__1 * d__1);
	    *e1 += r__;
	    if (abs(r__) <= abs(*e1) * 1e-15) {
		goto L15;
	    }
/* L10: */
	}
L15:
	ga = .5772156649015328;
	*e1 = -ga - log(*x) + *x * *e1;
    } else {
	m = (integer) (80.f / *x) + 20;
	t0 = 0.;
	for (k = m; k >= 1; --k) {
	    t0 = k / (k / (*x + t0) + 1.);
/* L20: */
	}
	t = 1. / (*x + t0);
	*e1 = exp(-(*x)) * t;
    }
    return 0;
} /* e1xb_ */

/*       ********************************** */
/* Subroutine */ int chgm_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *hg)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double exp(doublereal);
    void z_exp(doublecomplex *, doublecomplex *);
    double cos(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, m, n;
    static doublereal r__, y, a0, a1, r1, r2, x0, y0, y1;
    static integer la;
    static doublereal pi;
    static integer nl;
    static doublereal rg, xg, hg1, hg2;
    static doublecomplex cta, ctb;
    static doublereal tai, tbi, tar, tbr, sum1, sum2;
    static doublecomplex ctba;
    static doublereal tbai, tbar;
    extern /* Subroutine */ int cgama_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute confluent hypergeometric function */
/*                M(a,b,x) */
/*       Input  : a  --- Parameter */
/*                b  --- Parameter ( b <> 0,-1,-2,... ) */
/*                x  --- Argument */
/*       Output:  HG --- M(a,b,x) */
/*       Routine called: CGAMA for computing complex ln[Г(x)] */
/*       =================================================== */

    pi = 3.141592653589793;
    a0 = *a;
    a1 = *a;
    x0 = *x;
    *hg = 0.;
    if (*b == 0. || *b == (doublereal) (-(i__1 = (integer) (*b), abs(i__1)))) 
	    {
	*hg = 1e300;
    } else if (*a == 0. || *x == 0.) {
	*hg = 1.;
    } else if (*a == -1.) {
	*hg = 1. - *x / *b;
    } else if (*a == *b) {
	*hg = exp(*x);
    } else if (*a - *b == 1.) {
	*hg = (*x / *b + 1.) * exp(*x);
    } else if (*a == 1. && *b == 2.) {
	*hg = (exp(*x) - 1.) / *x;
    } else if (*a == (doublereal) ((integer) (*a)) && *a < 0.) {
	m = (integer) (-(*a));
	r__ = 1.;
	*hg = 1.;
	i__1 = m;
	for (k = 1; k <= i__1; ++k) {
	    r__ = r__ * (*a + k - 1.) / k / (*b + k - 1.) * *x;
/* L10: */
	    *hg += r__;
	}
    }
    if (*hg != 0.) {
	return 0;
    }
/*       DLMF 13.2.39 */
    if (*x < 0.) {
	*a = *b - *a;
	a0 = *a;
	*x = abs(*x);
    }
    nl = 0;
    la = 0;
    if (*a >= 2.) {
/*       preparing terms for DLMF 13.3.1 */
	nl = 1;
	la = (integer) (*a);
	*a = *a - la - 1.;
    }
    y0 = 0.;
    y1 = 0.;
    i__1 = nl;
    for (n = 0; n <= i__1; ++n) {
	if (a0 >= 2.) {
	    *a += 1.;
	}
	if (*x <= abs(*b) + 30. || *a < 0.) {
	    *hg = 1.;
	    rg = 1.;
	    for (j = 1; j <= 500; ++j) {
		rg = rg * (*a + j - 1.) / (j * (*b + j - 1.)) * *x;
		*hg += rg;
		if (*hg != 0. && (d__1 = rg / *hg, abs(d__1)) < 1e-15) {
/*       DLMF 13.2.39 (cf. above) */
		    if (x0 < 0.) {
			*hg *= exp(x0);
		    }
		    goto L25;
		}
/* L15: */
	    }
	} else {
/*       DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum */
	    y = 0.;
	    cgama_(a, &y, &c__0, &tar, &tai);
	    z__1.r = tar, z__1.i = tai;
	    cta.r = z__1.r, cta.i = z__1.i;
	    y = 0.;
	    cgama_(b, &y, &c__0, &tbr, &tbi);
	    z__1.r = tbr, z__1.i = tbi;
	    ctb.r = z__1.r, ctb.i = z__1.i;
	    xg = *b - *a;
	    y = 0.;
	    cgama_(&xg, &y, &c__0, &tbar, &tbai);
	    z__1.r = tbar, z__1.i = tbai;
	    ctba.r = z__1.r, ctba.i = z__1.i;
	    sum1 = 1.;
	    sum2 = 1.;
	    r1 = 1.;
	    r2 = 1.;
	    for (i__ = 1; i__ <= 8; ++i__) {
		r1 = -r1 * (*a + i__ - 1.) * (*a - *b + i__) / (*x * i__);
		r2 = -r2 * (*b - *a + i__ - 1.) * (*a - i__) / (*x * i__);
		sum1 += r1;
/* L20: */
		sum2 += r2;
	    }
	    if (x0 >= 0.) {
		z__2.r = ctb.r - ctba.r, z__2.i = ctb.i - ctba.i;
		z_exp(&z__1, &z__2);
		d__1 = -(*a);
		hg1 = z__1.r * pow_dd(x, &d__1) * cos(pi * *a) * sum1;
		z__3.r = ctb.r - cta.r, z__3.i = ctb.i - cta.i;
		z__2.r = z__3.r + *x, z__2.i = z__3.i;
		z_exp(&z__1, &z__2);
		d__1 = *a - *b;
		hg2 = z__1.r * pow_dd(x, &d__1) * sum2;
	    } else {
/*       DLMF 13.2.39 (cf. above) */
		z__3.r = ctb.r - ctba.r, z__3.i = ctb.i - ctba.i;
		z__2.r = z__3.r + x0, z__2.i = z__3.i;
		z_exp(&z__1, &z__2);
		d__1 = -(*a);
		hg1 = z__1.r * pow_dd(x, &d__1) * cos(pi * *a) * sum1;
		z__2.r = ctb.r - cta.r, z__2.i = ctb.i - cta.i;
		z_exp(&z__1, &z__2);
		d__1 = *a - *b;
		hg2 = z__1.r * pow_dd(x, &d__1) * sum2;
	    }
	    *hg = hg1 + hg2;
	}
L25:
	if (n == 0) {
	    y0 = *hg;
	}
	if (n == 1) {
	    y1 = *hg;
	}
/* L30: */
    }
    if (a0 >= 2.) {
/*       DLMF 13.3.1 */
	i__1 = la - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *hg = ((*a * 2. - *b + *x) * y1 + (*b - *a) * y0) / *a;
	    y0 = y1;
	    y1 = *hg;
/* L35: */
	    *a += 1.;
	}
    }
    *a = a1;
    *x = x0;
    return 0;
} /* chgm_ */

/*       ********************************** */
/* Subroutine */ int hygfx_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *x, doublereal *hf, integer *isfer)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), pow_di(
	    doublereal *, integer *), log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    extern /* Subroutine */ int psi_spec__(doublereal *, doublereal *);
    static integer j, k, m;
    static doublereal r__, a0, c0, c1, f0, g0, g1, g2, g3, f1;
    static logical l0, l1, l2, l3, l4, l5;
    static doublereal r0, r1, x1, aa, bb, ga, gb, gc, el, pa, pb, gm, pi;
    static integer nm;
    static doublereal hw, rm, sm, rp, sp, sp0, gca, gcb, gam, gbm, eps, gcab, 
	    gabc;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ==================================================== */
/*       Purpose: Compute hypergeometric function F(a,b,c,x) */
/*       Input :  a --- Parameter */
/*                b --- Parameter */
/*                c --- Parameter, c <> 0,-1,-2,... */
/*                x --- Argument   ( x < 1 ) */
/*       Output:  HF --- F(a,b,c,x) */
/*                ISFER --- Error flag */
/*       Routines called: */
/*            (1) GAMMA2 for computing gamma function */
/*            (2) PSI_SPEC for computing psi function */
/*       ==================================================== */

    pi = 3.141592653589793;
    el = .5772156649015329;
    *isfer = 0;
    l0 = *c__ == (doublereal) ((integer) (*c__)) && *c__ < 0.f;
    l1 = 1. - *x < 1e-15 && *c__ - *a - *b <= 0.f;
    l2 = *a == (doublereal) ((integer) (*a)) && *a < 0.f;
    l3 = *b == (doublereal) ((integer) (*b)) && *b < 0.f;
    l4 = *c__ - *a == (doublereal) ((integer) (*c__ - *a)) && *c__ - *a <= 
	    0.f;
    l5 = *c__ - *b == (doublereal) ((integer) (*c__ - *b)) && *c__ - *b <= 
	    0.f;
    if (l0 || l1) {
	*isfer = 3;
	return 0;
    }
    eps = 1e-15;
    if (*x > .95f) {
	eps = 1e-8;
    }
    if (*x == 0.f || *a == 0.f || *b == 0.f) {
	*hf = 1.;
	return 0;
    } else if (1. - *x == eps && *c__ - *a - *b > 0.f) {
	gamma2_(c__, &gc);
	d__1 = *c__ - *a - *b;
	gamma2_(&d__1, &gcab);
	d__1 = *c__ - *a;
	gamma2_(&d__1, &gca);
	d__1 = *c__ - *b;
	gamma2_(&d__1, &gcb);
	*hf = gc * gcab / (gca * gcb);
	return 0;
    } else if (*x + 1. <= eps && (d__1 = *c__ - *a + *b - 1.f, abs(d__1)) <= 
	    eps) {
	d__1 = -(*a);
	g0 = sqrt(pi) * pow_dd(&c_b4, &d__1);
	gamma2_(c__, &g1);
	d__1 = *a / 2.f + 1. - *b;
	gamma2_(&d__1, &g2);
	d__1 = *a * .5f + .5;
	gamma2_(&d__1, &g3);
	*hf = g0 * g1 / (g2 * g3);
	return 0;
    } else if (l2 || l3) {
	if (l2) {
	    nm = (integer) abs(*a);
	}
	if (l3) {
	    nm = (integer) abs(*b);
	}
	*hf = 1.;
	r__ = 1.;
	i__1 = nm;
	for (k = 1; k <= i__1; ++k) {
	    r__ = r__ * (*a + k - 1.) * (*b + k - 1.) / (k * (*c__ + k - 1.)) 
		    * *x;
/* L10: */
	    *hf += r__;
	}
	return 0;
    } else if (l4 || l5) {
	if (l4) {
	    nm = (d__1 = *c__ - *a, (integer) abs(d__1));
	}
	if (l5) {
	    nm = (d__1 = *c__ - *b, (integer) abs(d__1));
	}
	*hf = 1.;
	r__ = 1.;
	i__1 = nm;
	for (k = 1; k <= i__1; ++k) {
	    r__ = r__ * (*c__ - *a + k - 1.) * (*c__ - *b + k - 1.) / (k * (*
		    c__ + k - 1.)) * *x;
/* L15: */
	    *hf += r__;
	}
	d__1 = 1. - *x;
	d__2 = *c__ - *a - *b;
	*hf = pow_dd(&d__1, &d__2) * *hf;
	return 0;
    }
    aa = *a;
    bb = *b;
    x1 = *x;
    if (*x < 0.) {
	*x /= *x - 1.;
	if (*c__ > *a && *b < *a && *b > 0.f) {
	    *a = bb;
	    *b = aa;
	}
	*b = *c__ - *b;
    }
    hw = 0.;
    if (*x >= .75) {
	gm = 0.;
	if ((d__1 = *c__ - *a - *b - (integer) (*c__ - *a - *b), abs(d__1)) < 
		1e-15) {
	    m = (integer) (*c__ - *a - *b);
	    gamma2_(a, &ga);
	    gamma2_(b, &gb);
	    gamma2_(c__, &gc);
	    d__1 = *a + m;
	    gamma2_(&d__1, &gam);
	    d__1 = *b + m;
	    gamma2_(&d__1, &gbm);
	    psi_spec__(a, &pa);
	    psi_spec__(b, &pb);
	    if (m != 0) {
		gm = 1.;
	    }
	    i__1 = abs(m) - 1;
	    for (j = 1; j <= i__1; ++j) {
/* L30: */
		gm *= j;
	    }
	    rm = 1.;
	    i__1 = abs(m);
	    for (j = 1; j <= i__1; ++j) {
/* L35: */
		rm *= j;
	    }
	    f0 = 1.;
	    r0 = 1.;
	    r1 = 1.;
	    sp0 = 0.;
	    sp = 0.;
	    if (m >= 0) {
		c0 = gm * gc / (gam * gbm);
		d__1 = *x - 1.;
		c1 = -gc * pow_di(&d__1, &m) / (ga * gb * rm);
		i__1 = m - 1;
		for (k = 1; k <= i__1; ++k) {
		    r0 = r0 * (*a + k - 1.) * (*b + k - 1.f) / (k * (k - m)) *
			     (1.f - *x);
/* L40: */
		    f0 += r0;
		}
		i__1 = m;
		for (k = 1; k <= i__1; ++k) {
/* L45: */
		    sp0 = sp0 + 1. / (*a + k - 1.f) + 1.f / (*b + k - 1.f) - 
			    1.f / k;
		}
		f1 = pa + pb + sp0 + el * 2. + log(1. - *x);
		for (k = 1; k <= 250; ++k) {
		    sp = sp + (1. - *a) / (k * (*a + k - 1.f)) + (1.f - *b) / 
			    (k * (*b + k - 1.f));
		    sm = 0.;
		    i__1 = m;
		    for (j = 1; j <= i__1; ++j) {
/* L50: */
			sm = sm + (1. - *a) / ((j + k) * (*a + j + k - 1.f)) 
				+ 1.f / (*b + j + k - 1.f);
		    }
		    rp = pa + pb + el * 2. + sp + sm + log(1. - *x);
		    r1 = r1 * (*a + m + k - 1.) * (*b + m + k - 1.f) / (k * (
			    m + k)) * (1.f - *x);
		    f1 += r1 * rp;
		    if ((d__1 = f1 - hw, abs(d__1)) < abs(f1) * eps) {
			goto L60;
		    }
/* L55: */
		    hw = f1;
		}
L60:
		*hf = f0 * c0 + f1 * c1;
	    } else if (m < 0) {
		m = -m;
		d__1 = 1. - *x;
		c0 = gm * gc / (ga * gb * pow_di(&d__1, &m));
		c1 = -pow_ii(&c_n1, &m) * gc / (gam * gbm * rm);
		i__1 = m - 1;
		for (k = 1; k <= i__1; ++k) {
		    r0 = r0 * (*a - m + k - 1.) * (*b - m + k - 1.f) / (k * (
			    k - m)) * (1.f - *x);
/* L65: */
		    f0 += r0;
		}
		i__1 = m;
		for (k = 1; k <= i__1; ++k) {
/* L70: */
		    sp0 += 1. / k;
		}
		f1 = pa + pb - sp0 + el * 2. + log(1. - *x);
		for (k = 1; k <= 250; ++k) {
		    sp = sp + (1. - *a) / (k * (*a + k - 1.f)) + (1.f - *b) / 
			    (k * (*b + k - 1.f));
		    sm = 0.;
		    i__1 = m;
		    for (j = 1; j <= i__1; ++j) {
/* L75: */
			sm += 1. / (j + k);
		    }
		    rp = pa + pb + el * 2. + sp - sm + log(1. - *x);
		    r1 = r1 * (*a + k - 1.) * (*b + k - 1.f) / (k * (m + k)) *
			     (1.f - *x);
		    f1 += r1 * rp;
		    if ((d__1 = f1 - hw, abs(d__1)) < abs(f1) * eps) {
			goto L85;
		    }
/* L80: */
		    hw = f1;
		}
L85:
		*hf = f0 * c0 + f1 * c1;
	    }
	} else {
	    gamma2_(a, &ga);
	    gamma2_(b, &gb);
	    gamma2_(c__, &gc);
	    d__1 = *c__ - *a;
	    gamma2_(&d__1, &gca);
	    d__1 = *c__ - *b;
	    gamma2_(&d__1, &gcb);
	    d__1 = *c__ - *a - *b;
	    gamma2_(&d__1, &gcab);
	    d__1 = *a + *b - *c__;
	    gamma2_(&d__1, &gabc);
	    c0 = gc * gcab / (gca * gcb);
	    d__1 = 1. - *x;
	    d__2 = *c__ - *a - *b;
	    c1 = gc * gabc / (ga * gb) * pow_dd(&d__1, &d__2);
	    *hf = 0.;
	    r0 = c0;
	    r1 = c1;
	    for (k = 1; k <= 250; ++k) {
		r0 = r0 * (*a + k - 1.) * (*b + k - 1.f) / (k * (*a + *b - *
			c__ + k)) * (1.f - *x);
		r1 = r1 * (*c__ - *a + k - 1.) * (*c__ - *b + k - 1.f) / (k * 
			(*c__ - *a - *b + k)) * (1.f - *x);
		*hf = *hf + r0 + r1;
		if ((d__1 = *hf - hw, abs(d__1)) < abs(*hf) * eps) {
		    goto L95;
		}
/* L90: */
		hw = *hf;
	    }
L95:
	    *hf = *hf + c0 + c1;
	}
    } else {
	a0 = 1.;
	if (*c__ > *a && *c__ < *a * 2. && *c__ > *b && *c__ < *b * 2.) {
	    d__1 = 1. - *x;
	    d__2 = *c__ - *a - *b;
	    a0 = pow_dd(&d__1, &d__2);
	    *a = *c__ - *a;
	    *b = *c__ - *b;
	}
	*hf = 1.;
	r__ = 1.;
	for (k = 1; k <= 250; ++k) {
	    r__ = r__ * (*a + k - 1.) * (*b + k - 1.) / (k * (*c__ + k - 1.)) 
		    * *x;
	    *hf += r__;
	    if ((d__1 = *hf - hw, abs(d__1)) <= abs(*hf) * eps) {
		goto L105;
	    }
/* L100: */
	    hw = *hf;
	}
L105:
	*hf = a0 * *hf;
    }
    if (x1 < 0.) {
	*x = x1;
	d__1 = 1. - *x;
	c0 = 1. / pow_dd(&d__1, &aa);
	*hf = c0 * *hf;
    }
    *a = aa;
    *b = bb;
    if (k > 120) {
	*isfer = 5;
    }
    return 0;
} /* hygfx_ */

/*       ********************************** */
/* Subroutine */ int cchg_(doublereal *a, doublereal *b, doublecomplex *z__, 
	doublecomplex *chg)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    void z_exp(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), d_imag(doublecomplex *), atan(doublereal), 
	    cos(doublereal);
    void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, m, n;
    static doublereal x, y, a0, a1, x0;
    static doublecomplex z0;
    static doublereal ba;
    static doublecomplex ci;
    static integer la;
    static doublecomplex cr;
    static doublereal pi;
    static integer nl, ns;
    static doublecomplex cg1, cg2, cg3;
    static doublereal g1i, g2i, g3i;
    static doublecomplex cr1, cs1, cs2, cr2;
    static doublereal g1r, g2r, g3r;
    static doublecomplex cy0, cy1, crg;
    static doublereal phi;
    static doublecomplex chw, chg1, chg2, cfac;
    extern /* Subroutine */ int cgama_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute confluent hypergeometric function */
/*                M(a,b,z) with real parameters a, b and a */
/*                complex argument z */
/*       Input :  a --- Parameter */
/*                b --- Parameter */
/*                z --- Complex argument */
/*       Output:  CHG --- M(a,b,z) */
/*       Routine called: CGAMA for computing complex ln[Г(x)] */
/*       =================================================== */

    pi = 3.141592653589793;
    ci.r = 0., ci.i = 1.;
    a0 = *a;
    a1 = *a;
    z0.r = z__->r, z0.i = z__->i;
    if (*b == 0.f || *b == (doublereal) (-((integer) abs(*b)))) {
	chg->r = 1e300, chg->i = 0.;
    } else if (*a == 0. || z__->r == 0. && z__->i == 0.) {
	chg->r = 1., chg->i = 0.;
    } else if (*a == -1.) {
	z__2.r = z__->r / *b, z__2.i = z__->i / *b;
	z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
	chg->r = z__1.r, chg->i = z__1.i;
    } else if (*a == *b) {
	z_exp(&z__1, z__);
	chg->r = z__1.r, chg->i = z__1.i;
    } else if (*a - *b == 1.) {
	z__3.r = z__->r / *b, z__3.i = z__->i / *b;
	z__2.r = z__3.r + 1., z__2.i = z__3.i;
	z_exp(&z__4, z__);
	z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * z__4.i 
		+ z__2.i * z__4.r;
	chg->r = z__1.r, chg->i = z__1.i;
    } else if (*a == 1. && *b == 2.) {
	z_exp(&z__3, z__);
	z__2.r = z__3.r - 1., z__2.i = z__3.i;
	z_div(&z__1, &z__2, z__);
	chg->r = z__1.r, chg->i = z__1.i;
    } else if (*a == (doublereal) ((integer) (*a)) && *a < 0.) {
	m = (integer) (-(*a));
	cr.r = 1., cr.i = 0.;
	chg->r = 1., chg->i = 0.;
	i__1 = m;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = *a + k - 1.;
	    z__4.r = d__1 * cr.r, z__4.i = d__1 * cr.i;
	    d__2 = (doublereal) k;
	    z__3.r = z__4.r / d__2, z__3.i = z__4.i / d__2;
	    d__3 = *b + k - 1.;
	    z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
	    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * 
		    z__->i + z__2.i * z__->r;
	    cr.r = z__1.r, cr.i = z__1.i;
/* L10: */
	    z__1.r = chg->r + cr.r, z__1.i = chg->i + cr.i;
	    chg->r = z__1.r, chg->i = z__1.i;
	}
    } else {
	x0 = z__->r;
	if (x0 < 0.) {
	    *a = *b - *a;
	    a0 = *a;
	    z__1.r = -z__->r, z__1.i = -z__->i;
	    z__->r = z__1.r, z__->i = z__1.i;
	}
	nl = 0;
	la = 0;
	if (*a >= 2.) {
	    nl = 1;
	    la = (integer) (*a);
	    *a = *a - la - 1.;
	}
	ns = 0;
	i__1 = nl;
	for (n = 0; n <= i__1; ++n) {
	    if (a0 >= 2.) {
		*a += 1.;
	    }
	    if (z_abs(z__) < abs(*b) + 20. || *a < 0.) {
		chg->r = 1., chg->i = 0.;
		crg.r = 1., crg.i = 0.;
		for (j = 1; j <= 500; ++j) {
		    d__1 = *a + j - 1.;
		    z__3.r = d__1 * crg.r, z__3.i = d__1 * crg.i;
		    d__2 = j * (*b + j - 1.);
		    z__2.r = z__3.r / d__2, z__2.i = z__3.i / d__2;
		    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = 
			    z__2.r * z__->i + z__2.i * z__->r;
		    crg.r = z__1.r, crg.i = z__1.i;
		    z__1.r = chg->r + crg.r, z__1.i = chg->i + crg.i;
		    chg->r = z__1.r, chg->i = z__1.i;
		    z__2.r = chg->r - chw.r, z__2.i = chg->i - chw.i;
		    z_div(&z__1, &z__2, chg);
		    if (z_abs(&z__1) < 1e-15) {
			goto L25;
		    }
		    chw.r = chg->r, chw.i = chg->i;
/* L15: */
		}
	    } else {
		y = 0.;
		cgama_(a, &y, &c__0, &g1r, &g1i);
		z__1.r = g1r, z__1.i = g1i;
		cg1.r = z__1.r, cg1.i = z__1.i;
		y = 0.;
		cgama_(b, &y, &c__0, &g2r, &g2i);
		z__1.r = g2r, z__1.i = g2i;
		cg2.r = z__1.r, cg2.i = z__1.i;
		ba = *b - *a;
		y = 0.;
		cgama_(&ba, &y, &c__0, &g3r, &g3i);
		z__1.r = g3r, z__1.i = g3i;
		cg3.r = z__1.r, cg3.i = z__1.i;
		cs1.r = 1., cs1.i = 0.;
		cs2.r = 1., cs2.i = 0.;
		cr1.r = 1., cr1.i = 0.;
		cr2.r = 1., cr2.i = 0.;
		for (i__ = 1; i__ <= 8; ++i__) {
		    z__4.r = -cr1.r, z__4.i = -cr1.i;
		    d__1 = *a + i__ - 1.;
		    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
		    d__2 = *a - *b + i__;
		    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
		    d__3 = (doublereal) i__;
		    z__5.r = d__3 * z__->r, z__5.i = d__3 * z__->i;
		    z_div(&z__1, &z__2, &z__5);
		    cr1.r = z__1.r, cr1.i = z__1.i;
		    d__1 = *b - *a + i__ - 1.;
		    z__3.r = d__1 * cr2.r, z__3.i = d__1 * cr2.i;
		    d__2 = i__ - *a;
		    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
		    d__3 = (doublereal) i__;
		    z__4.r = d__3 * z__->r, z__4.i = d__3 * z__->i;
		    z_div(&z__1, &z__2, &z__4);
		    cr2.r = z__1.r, cr2.i = z__1.i;
		    z__1.r = cs1.r + cr1.r, z__1.i = cs1.i + cr1.i;
		    cs1.r = z__1.r, cs1.i = z__1.i;
/* L20: */
		    z__1.r = cs2.r + cr2.r, z__1.i = cs2.i + cr2.i;
		    cs2.r = z__1.r, cs2.i = z__1.i;
		}
		x = z__->r;
		y = d_imag(z__);
		if (x == 0.f && y >= 0.f) {
		    phi = pi * .5;
		} else if (x == 0.f && y <= 0.f) {
		    phi = pi * -.5;
		} else {
		    phi = atan(y / x);
		}
		if (phi > pi * -.5f && phi < pi * 1.5f) {
		    ns = 1;
		}
		if (phi > pi * -1.5f && phi <= pi * -.5f) {
		    ns = -1;
		}
		d__1 = (doublereal) ns;
		z__4.r = d__1 * ci.r, z__4.i = d__1 * ci.i;
		z__3.r = pi * z__4.r, z__3.i = pi * z__4.i;
		z__2.r = *a * z__3.r, z__2.i = *a * z__3.i;
		z_exp(&z__1, &z__2);
		cfac.r = z__1.r, cfac.i = z__1.i;
		if (y == 0.) {
		    d__1 = cos(pi * *a);
		    cfac.r = d__1, cfac.i = 0.;
		}
		z__5.r = cg2.r - cg3.r, z__5.i = cg2.i - cg3.i;
		z_exp(&z__4, &z__5);
		d__1 = -(*a);
		z__7.r = d__1, z__7.i = 0.;
		pow_zz(&z__6, z__, &z__7);
		z__3.r = z__4.r * z__6.r - z__4.i * z__6.i, z__3.i = z__4.r * 
			z__6.i + z__4.i * z__6.r;
		z__2.r = z__3.r * cfac.r - z__3.i * cfac.i, z__2.i = z__3.r * 
			cfac.i + z__3.i * cfac.r;
		z__1.r = z__2.r * cs1.r - z__2.i * cs1.i, z__1.i = z__2.r * 
			cs1.i + z__2.i * cs1.r;
		chg1.r = z__1.r, chg1.i = z__1.i;
		z__5.r = cg2.r - cg1.r, z__5.i = cg2.i - cg1.i;
		z__4.r = z__5.r + z__->r, z__4.i = z__5.i + z__->i;
		z_exp(&z__3, &z__4);
		d__1 = *a - *b;
		z__7.r = d__1, z__7.i = 0.;
		pow_zz(&z__6, z__, &z__7);
		z__2.r = z__3.r * z__6.r - z__3.i * z__6.i, z__2.i = z__3.r * 
			z__6.i + z__3.i * z__6.r;
		z__1.r = z__2.r * cs2.r - z__2.i * cs2.i, z__1.i = z__2.r * 
			cs2.i + z__2.i * cs2.r;
		chg2.r = z__1.r, chg2.i = z__1.i;
		z__1.r = chg1.r + chg2.r, z__1.i = chg1.i + chg2.i;
		chg->r = z__1.r, chg->i = z__1.i;
	    }
L25:
	    if (n == 0) {
		cy0.r = chg->r, cy0.i = chg->i;
	    }
	    if (n == 1) {
		cy1.r = chg->r, cy1.i = chg->i;
	    }
/* L30: */
	}
	if (a0 >= 2.) {
	    i__1 = la - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = *a * 2. - *b;
		z__4.r = d__1 + z__->r, z__4.i = z__->i;
		z__3.r = z__4.r * cy1.r - z__4.i * cy1.i, z__3.i = z__4.r * 
			cy1.i + z__4.i * cy1.r;
		d__2 = *b - *a;
		z__5.r = d__2 * cy0.r, z__5.i = d__2 * cy0.i;
		z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
		z__1.r = z__2.r / *a, z__1.i = z__2.i / *a;
		chg->r = z__1.r, chg->i = z__1.i;
		cy0.r = cy1.r, cy0.i = cy1.i;
		cy1.r = chg->r, cy1.i = chg->i;
/* L35: */
		*a += 1.;
	    }
	}
	if (x0 < 0.) {
	    z__3.r = -z__->r, z__3.i = -z__->i;
	    z_exp(&z__2, &z__3);
	    z__1.r = chg->r * z__2.r - chg->i * z__2.i, z__1.i = chg->r * 
		    z__2.i + chg->i * z__2.r;
	    chg->r = z__1.r, chg->i = z__1.i;
	}
    }
    *a = a1;
    z__->r = z0.r, z__->i = z0.i;
    return 0;
} /* cchg_ */

/*       ********************************** */
/* Subroutine */ int hygfz_(doublereal *a, doublereal *b, doublereal *c__, 
	doublecomplex *z__, doublecomplex *zhf, integer *isfer)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *), sqrt(doublereal), 
	    pow_dd(doublereal *, doublereal *);
    void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *), z_div(
	    doublecomplex *, doublecomplex *, doublecomplex *);
    double d_sign(doublereal *, doublereal *);
    void pow_zi(doublecomplex *, doublecomplex *, integer *), z_log(
	    doublecomplex *, doublecomplex *);
    integer pow_ii(integer *, integer *);
    double tan(doublereal);

    /* Local variables */
    extern /* Subroutine */ int psi_spec__(doublereal *, doublereal *);
    static integer j, k, m;
    static doublereal x, y, a0, g0, g1, g2, g3;
    static logical l0, l1, l2, l3, l4, l5, l6;
    static doublereal t0, w0;
    static doublecomplex z1;
    static doublereal aa, bb, ca, cb, ga, gb, gc, el, pa, pb, gm, pi;
    static doublecomplex z00;
    static integer nm;
    static doublereal rm, sm, sp, sq;
    static doublecomplex zp;
    static doublereal ws;
    static doublecomplex zr, zw, zc0, zc1;
    static doublereal rk1;
    static doublecomplex zf0, zf1;
    static doublereal sj1, sp0, rk2, sj2;
    static doublecomplex zp0, zr0, zr1;
    static doublereal gab, gca, gcb, gba;
    static integer mab, nca, ncb;
    static doublereal pca, gam, gbm, pac, eps, gcab, gabc;
    static integer mcab;
    static doublereal gcbk;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute the hypergeometric function for a */
/*                complex argument, F(a,b,c,z) */
/*       Input :  a --- Parameter */
/*                b --- Parameter */
/*                c --- Parameter,  c <> 0,-1,-2,... */
/*                z --- Complex argument */
/*       Output:  ZHF --- F(a,b,c,z) */
/*                ISFER --- Error flag */
/*       Routines called: */
/*            (1) GAMMA2 for computing gamma function */
/*            (2) PSI_SPEC for computing psi function */
/*       ====================================================== */

    x = z__->r;
    y = d_imag(z__);
    eps = 1e-15;
    *isfer = 0;
    l0 = *c__ == (doublereal) ((integer) (*c__)) && *c__ < 0.;
    l1 = (d__1 = 1. - x, abs(d__1)) < eps && y == 0. && *c__ - *a - *b <= 0.;
    z__1.r = z__->r + 1., z__1.i = z__->i;
    l2 = z_abs(&z__1) < eps && (d__1 = *c__ - *a + *b - 1., abs(d__1)) < eps;
    l3 = *a == (doublereal) ((integer) (*a)) && *a < 0.;
    l4 = *b == (doublereal) ((integer) (*b)) && *b < 0.;
    l5 = *c__ - *a == (doublereal) ((integer) (*c__ - *a)) && *c__ - *a <= 0.;
    l6 = *c__ - *b == (doublereal) ((integer) (*c__ - *b)) && *c__ - *b <= 0.;
    aa = *a;
    bb = *b;
    a0 = z_abs(z__);
    if (a0 > .95) {
	eps = 1e-8;
    }
    pi = 3.141592653589793;
    el = .5772156649015329;
    if (l0 || l1) {
	*isfer = 3;
	return 0;
    }
    nm = 0;
    if (a0 == 0. || *a == 0. || *b == 0.) {
	zhf->r = 1., zhf->i = 0.;
    } else if (z__->r == 1. && z__->i == 0. && *c__ - *a - *b > 0.) {
	gamma2_(c__, &gc);
	d__1 = *c__ - *a - *b;
	gamma2_(&d__1, &gcab);
	d__1 = *c__ - *a;
	gamma2_(&d__1, &gca);
	d__1 = *c__ - *b;
	gamma2_(&d__1, &gcb);
	d__1 = gc * gcab / (gca * gcb);
	zhf->r = d__1, zhf->i = 0.;
    } else if (l2) {
	d__1 = -(*a);
	g0 = sqrt(pi) * pow_dd(&c_b4, &d__1);
	gamma2_(c__, &g1);
	d__1 = *a / 2. + 1. - *b;
	gamma2_(&d__1, &g2);
	d__1 = *a * .5 + .5;
	gamma2_(&d__1, &g3);
	d__1 = g0 * g1 / (g2 * g3);
	zhf->r = d__1, zhf->i = 0.;
    } else if (l3 || l4) {
	if (l3) {
	    nm = (integer) abs(*a);
	}
	if (l4) {
	    nm = (integer) abs(*b);
	}
	zhf->r = 1., zhf->i = 0.;
	zr.r = 1., zr.i = 0.;
	i__1 = nm;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = *a + k - 1.;
	    z__4.r = d__1 * zr.r, z__4.i = d__1 * zr.i;
	    d__2 = *b + k - 1.;
	    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
	    d__3 = k * (*c__ + k - 1.);
	    z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
	    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * 
		    z__->i + z__2.i * z__->r;
	    zr.r = z__1.r, zr.i = z__1.i;
/* L10: */
	    z__1.r = zhf->r + zr.r, z__1.i = zhf->i + zr.i;
	    zhf->r = z__1.r, zhf->i = z__1.i;
	}
    } else if (l5 || l6) {
	if (l5) {
	    nm = (d__1 = *c__ - *a, (integer) abs(d__1));
	}
	if (l6) {
	    nm = (d__1 = *c__ - *b, (integer) abs(d__1));
	}
	zhf->r = 1., zhf->i = 0.;
	zr.r = 1., zr.i = 0.;
	i__1 = nm;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = *c__ - *a + k - 1.;
	    z__4.r = d__1 * zr.r, z__4.i = d__1 * zr.i;
	    d__2 = *c__ - *b + k - 1.;
	    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
	    d__3 = k * (*c__ + k - 1.);
	    z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
	    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * 
		    z__->i + z__2.i * z__->r;
	    zr.r = z__1.r, zr.i = z__1.i;
/* L15: */
	    z__1.r = zhf->r + zr.r, z__1.i = zhf->i + zr.i;
	    zhf->r = z__1.r, zhf->i = z__1.i;
	}
	z__3.r = 1. - z__->r, z__3.i = -z__->i;
	d__1 = *c__ - *a - *b;
	z__4.r = d__1, z__4.i = 0.;
	pow_zz(&z__2, &z__3, &z__4);
	z__1.r = z__2.r * zhf->r - z__2.i * zhf->i, z__1.i = z__2.r * zhf->i 
		+ z__2.i * zhf->r;
	zhf->r = z__1.r, zhf->i = z__1.i;
    } else if (a0 <= 1.) {
	if (x < 0.) {
	    z__2.r = z__->r - 1., z__2.i = z__->i;
	    z_div(&z__1, z__, &z__2);
	    z1.r = z__1.r, z1.i = z__1.i;
	    if (*c__ > *a && *b < *a && *b > 0.f) {
		*a = bb;
		*b = aa;
	    }
	    z__3.r = 1. - z__->r, z__3.i = -z__->i;
	    z__4.r = *a, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z_div(&z__1, &c_b6, &z__2);
	    zc0.r = z__1.r, zc0.i = z__1.i;
	    zhf->r = 1., zhf->i = 0.;
	    zr0.r = 1., zr0.i = 0.;
	    for (k = 1; k <= 500; ++k) {
		d__1 = *a + k - 1.;
		z__4.r = d__1 * zr0.r, z__4.i = d__1 * zr0.i;
		d__2 = *c__ - *b + k - 1.;
		z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
		d__3 = k * (*c__ + k - 1.);
		z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
		z__1.r = z__2.r * z1.r - z__2.i * z1.i, z__1.i = z__2.r * 
			z1.i + z__2.i * z1.r;
		zr0.r = z__1.r, zr0.i = z__1.i;
		z__1.r = zhf->r + zr0.r, z__1.i = zhf->i + zr0.i;
		zhf->r = z__1.r, zhf->i = z__1.i;
		z__1.r = zhf->r - zw.r, z__1.i = zhf->i - zw.i;
		if (z_abs(&z__1) < z_abs(zhf) * eps) {
		    goto L25;
		}
/* L20: */
		zw.r = zhf->r, zw.i = zhf->i;
	    }
L25:
	    z__1.r = zc0.r * zhf->r - zc0.i * zhf->i, z__1.i = zc0.r * zhf->i 
		    + zc0.i * zhf->r;
	    zhf->r = z__1.r, zhf->i = z__1.i;
	} else if (a0 >= .9) {
	    gm = 0.;
	    d__1 = *c__ - *a - *b;
	    mcab = (integer) (*c__ - *a - *b + eps * d_sign(&c_b524, &d__1));
	    if ((d__1 = *c__ - *a - *b - mcab, abs(d__1)) < eps) {
		m = (integer) (*c__ - *a - *b);
		gamma2_(a, &ga);
		gamma2_(b, &gb);
		gamma2_(c__, &gc);
		d__1 = *a + m;
		gamma2_(&d__1, &gam);
		d__1 = *b + m;
		gamma2_(&d__1, &gbm);
		psi_spec__(a, &pa);
		psi_spec__(b, &pb);
		if (m != 0) {
		    gm = 1.;
		}
		i__1 = abs(m) - 1;
		for (j = 1; j <= i__1; ++j) {
/* L30: */
		    gm *= j;
		}
		rm = 1.;
		i__1 = abs(m);
		for (j = 1; j <= i__1; ++j) {
/* L35: */
		    rm *= j;
		}
		zf0.r = 1., zf0.i = 0.;
		zr0.r = 1., zr0.i = 0.;
		zr1.r = 1., zr1.i = 0.;
		sp0 = 0.;
		sp = 0.;
		if (m >= 0) {
		    d__1 = gm * gc / (gam * gbm);
		    zc0.r = d__1, zc0.i = 0.;
		    d__1 = -gc;
		    z__4.r = z__->r - 1., z__4.i = z__->i;
		    pow_zi(&z__3, &z__4, &m);
		    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
		    d__2 = ga * gb * rm;
		    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
		    zc1.r = z__1.r, zc1.i = z__1.i;
		    i__1 = m - 1;
		    for (k = 1; k <= i__1; ++k) {
			d__1 = *a + k - 1.;
			z__4.r = d__1 * zr0.r, z__4.i = d__1 * zr0.i;
			d__2 = *b + k - 1.;
			z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
			i__2 = k * (k - m);
			d__3 = (doublereal) i__2;
			z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
			z__5.r = 1. - z__->r, z__5.i = -z__->i;
			z__1.r = z__2.r * z__5.r - z__2.i * z__5.i, z__1.i = 
				z__2.r * z__5.i + z__2.i * z__5.r;
			zr0.r = z__1.r, zr0.i = z__1.i;
/* L40: */
			z__1.r = zf0.r + zr0.r, z__1.i = zf0.i + zr0.i;
			zf0.r = z__1.r, zf0.i = z__1.i;
		    }
		    i__1 = m;
		    for (k = 1; k <= i__1; ++k) {
/* L45: */
			sp0 = sp0 + 1. / (*a + k - 1.) + 1.f / (*b + k - 1.) 
				- 1. / k;
		    }
		    d__1 = pa + pb + sp0 + el * 2.;
		    z__3.r = 1. - z__->r, z__3.i = -z__->i;
		    z_log(&z__2, &z__3);
		    z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
		    zf1.r = z__1.r, zf1.i = z__1.i;
		    for (k = 1; k <= 500; ++k) {
			sp = sp + (1. - *a) / (k * (*a + k - 1.)) + (1. - *b) 
				/ (k * (*b + k - 1.));
			sm = 0.;
			i__1 = m;
			for (j = 1; j <= i__1; ++j) {
			    sm = sm + (1. - *a) / ((j + k) * (*a + j + k - 1.)
				    ) + 1. / (*b + j + k - 1.);
/* L50: */
			}
			d__1 = pa + pb + el * 2. + sp + sm;
			z__3.r = 1. - z__->r, z__3.i = -z__->i;
			z_log(&z__2, &z__3);
			z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
			zp.r = z__1.r, zp.i = z__1.i;
			d__1 = *a + m + k - 1.;
			z__4.r = d__1 * zr1.r, z__4.i = d__1 * zr1.i;
			d__2 = *b + m + k - 1.;
			z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
			i__1 = k * (m + k);
			d__3 = (doublereal) i__1;
			z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
			z__5.r = 1. - z__->r, z__5.i = -z__->i;
			z__1.r = z__2.r * z__5.r - z__2.i * z__5.i, z__1.i = 
				z__2.r * z__5.i + z__2.i * z__5.r;
			zr1.r = z__1.r, zr1.i = z__1.i;
			z__2.r = zr1.r * zp.r - zr1.i * zp.i, z__2.i = zr1.r *
				 zp.i + zr1.i * zp.r;
			z__1.r = zf1.r + z__2.r, z__1.i = zf1.i + z__2.i;
			zf1.r = z__1.r, zf1.i = z__1.i;
			z__1.r = zf1.r - zw.r, z__1.i = zf1.i - zw.i;
			if (z_abs(&z__1) < z_abs(&zf1) * eps) {
			    goto L60;
			}
/* L55: */
			zw.r = zf1.r, zw.i = zf1.i;
		    }
L60:
		    z__2.r = zf0.r * zc0.r - zf0.i * zc0.i, z__2.i = zf0.r * 
			    zc0.i + zf0.i * zc0.r;
		    z__3.r = zf1.r * zc1.r - zf1.i * zc1.i, z__3.i = zf1.r * 
			    zc1.i + zf1.i * zc1.r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    zhf->r = z__1.r, zhf->i = z__1.i;
		} else if (m < 0) {
		    m = -m;
		    d__1 = gm * gc;
		    z__2.r = d__1, z__2.i = 0.;
		    d__2 = ga * gb;
		    z__5.r = 1. - z__->r, z__5.i = -z__->i;
		    pow_zi(&z__4, &z__5, &m);
		    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
		    z_div(&z__1, &z__2, &z__3);
		    zc0.r = z__1.r, zc0.i = z__1.i;
		    d__1 = -pow_ii(&c_n1, &m) * gc / (gam * gbm * rm);
		    zc1.r = d__1, zc1.i = 0.;
		    i__1 = m - 1;
		    for (k = 1; k <= i__1; ++k) {
			d__1 = *a - m + k - 1.;
			z__4.r = d__1 * zr0.r, z__4.i = d__1 * zr0.i;
			d__2 = *b - m + k - 1.;
			z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
			i__2 = k * (k - m);
			d__3 = (doublereal) i__2;
			z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
			z__5.r = 1. - z__->r, z__5.i = -z__->i;
			z__1.r = z__2.r * z__5.r - z__2.i * z__5.i, z__1.i = 
				z__2.r * z__5.i + z__2.i * z__5.r;
			zr0.r = z__1.r, zr0.i = z__1.i;
/* L65: */
			z__1.r = zf0.r + zr0.r, z__1.i = zf0.i + zr0.i;
			zf0.r = z__1.r, zf0.i = z__1.i;
		    }
		    i__1 = m;
		    for (k = 1; k <= i__1; ++k) {
/* L70: */
			sp0 += 1. / k;
		    }
		    d__1 = pa + pb - sp0 + el * 2.;
		    z__3.r = 1. - z__->r, z__3.i = -z__->i;
		    z_log(&z__2, &z__3);
		    z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
		    zf1.r = z__1.r, zf1.i = z__1.i;
		    for (k = 1; k <= 500; ++k) {
			sp = sp + (1. - *a) / (k * (*a + k - 1.)) + (1. - *b) 
				/ (k * (*b + k - 1.));
			sm = 0.;
			i__1 = m;
			for (j = 1; j <= i__1; ++j) {
/* L75: */
			    sm += 1. / (j + k);
			}
			d__1 = pa + pb + el * 2. + sp - sm;
			z__3.r = 1. - z__->r, z__3.i = -z__->i;
			z_log(&z__2, &z__3);
			z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
			zp.r = z__1.r, zp.i = z__1.i;
			d__1 = *a + k - 1.;
			z__4.r = d__1 * zr1.r, z__4.i = d__1 * zr1.i;
			d__2 = *b + k - 1.;
			z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
			i__1 = k * (m + k);
			d__3 = (doublereal) i__1;
			z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
			z__5.r = 1. - z__->r, z__5.i = -z__->i;
			z__1.r = z__2.r * z__5.r - z__2.i * z__5.i, z__1.i = 
				z__2.r * z__5.i + z__2.i * z__5.r;
			zr1.r = z__1.r, zr1.i = z__1.i;
			z__2.r = zr1.r * zp.r - zr1.i * zp.i, z__2.i = zr1.r *
				 zp.i + zr1.i * zp.r;
			z__1.r = zf1.r + z__2.r, z__1.i = zf1.i + z__2.i;
			zf1.r = z__1.r, zf1.i = z__1.i;
			z__1.r = zf1.r - zw.r, z__1.i = zf1.i - zw.i;
			if (z_abs(&z__1) < z_abs(&zf1) * eps) {
			    goto L85;
			}
/* L80: */
			zw.r = zf1.r, zw.i = zf1.i;
		    }
L85:
		    z__2.r = zf0.r * zc0.r - zf0.i * zc0.i, z__2.i = zf0.r * 
			    zc0.i + zf0.i * zc0.r;
		    z__3.r = zf1.r * zc1.r - zf1.i * zc1.i, z__3.i = zf1.r * 
			    zc1.i + zf1.i * zc1.r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    zhf->r = z__1.r, zhf->i = z__1.i;
		}
	    } else {
		gamma2_(a, &ga);
		gamma2_(b, &gb);
		gamma2_(c__, &gc);
		d__1 = *c__ - *a;
		gamma2_(&d__1, &gca);
		d__1 = *c__ - *b;
		gamma2_(&d__1, &gcb);
		d__1 = *c__ - *a - *b;
		gamma2_(&d__1, &gcab);
		d__1 = *a + *b - *c__;
		gamma2_(&d__1, &gabc);
		d__1 = gc * gcab / (gca * gcb);
		zc0.r = d__1, zc0.i = 0.;
		d__1 = gc * gabc / (ga * gb);
		z__3.r = 1. - z__->r, z__3.i = -z__->i;
		d__2 = *c__ - *a - *b;
		z__4.r = d__2, z__4.i = 0.;
		pow_zz(&z__2, &z__3, &z__4);
		z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
		zc1.r = z__1.r, zc1.i = z__1.i;
		zhf->r = 0., zhf->i = 0.;
		zr0.r = zc0.r, zr0.i = zc0.i;
		zr1.r = zc1.r, zr1.i = zc1.i;
		for (k = 1; k <= 500; ++k) {
		    d__1 = *a + k - 1.;
		    z__4.r = d__1 * zr0.r, z__4.i = d__1 * zr0.i;
		    d__2 = *b + k - 1.;
		    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
		    d__3 = k * (*a + *b - *c__ + k);
		    z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
		    z__5.r = 1. - z__->r, z__5.i = -z__->i;
		    z__1.r = z__2.r * z__5.r - z__2.i * z__5.i, z__1.i = 
			    z__2.r * z__5.i + z__2.i * z__5.r;
		    zr0.r = z__1.r, zr0.i = z__1.i;
		    d__1 = *c__ - *a + k - 1.;
		    z__4.r = d__1 * zr1.r, z__4.i = d__1 * zr1.i;
		    d__2 = *c__ - *b + k - 1.;
		    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
		    d__3 = k * (*c__ - *a - *b + k);
		    z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
		    z__5.r = 1. - z__->r, z__5.i = -z__->i;
		    z__1.r = z__2.r * z__5.r - z__2.i * z__5.i, z__1.i = 
			    z__2.r * z__5.i + z__2.i * z__5.r;
		    zr1.r = z__1.r, zr1.i = z__1.i;
		    z__2.r = zhf->r + zr0.r, z__2.i = zhf->i + zr0.i;
		    z__1.r = z__2.r + zr1.r, z__1.i = z__2.i + zr1.i;
		    zhf->r = z__1.r, zhf->i = z__1.i;
		    z__1.r = zhf->r - zw.r, z__1.i = zhf->i - zw.i;
		    if (z_abs(&z__1) < z_abs(zhf) * eps) {
			goto L95;
		    }
/* L90: */
		    zw.r = zhf->r, zw.i = zhf->i;
		}
L95:
		z__2.r = zhf->r + zc0.r, z__2.i = zhf->i + zc0.i;
		z__1.r = z__2.r + zc1.r, z__1.i = z__2.i + zc1.i;
		zhf->r = z__1.r, zhf->i = z__1.i;
	    }
	} else {
	    z00.r = 1., z00.i = 0.;
	    if (*c__ - *a < *a && *c__ - *b < *b) {
		z__2.r = 1. - z__->r, z__2.i = -z__->i;
		d__1 = *c__ - *a - *b;
		z__3.r = d__1, z__3.i = 0.;
		pow_zz(&z__1, &z__2, &z__3);
		z00.r = z__1.r, z00.i = z__1.i;
		*a = *c__ - *a;
		*b = *c__ - *b;
	    }
	    zhf->r = 1., zhf->i = 0.;
	    zr.r = 1., zr.i = 0.;
	    for (k = 1; k <= 1500; ++k) {
		d__1 = *a + k - 1.;
		z__4.r = d__1 * zr.r, z__4.i = d__1 * zr.i;
		d__2 = *b + k - 1.;
		z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
		d__3 = k * (*c__ + k - 1.);
		z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
		z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * 
			z__->i + z__2.i * z__->r;
		zr.r = z__1.r, zr.i = z__1.i;
		z__1.r = zhf->r + zr.r, z__1.i = zhf->i + zr.i;
		zhf->r = z__1.r, zhf->i = z__1.i;
		z__1.r = zhf->r - zw.r, z__1.i = zhf->i - zw.i;
		if (z_abs(&z__1) <= z_abs(zhf) * eps) {
		    goto L105;
		}
/* L100: */
		zw.r = zhf->r, zw.i = zhf->i;
	    }
L105:
	    z__1.r = z00.r * zhf->r - z00.i * zhf->i, z__1.i = z00.r * zhf->i 
		    + z00.i * zhf->r;
	    zhf->r = z__1.r, zhf->i = z__1.i;
	}
    } else if (a0 > 1.) {
	d__1 = *a - *b;
	mab = (integer) (*a - *b + eps * d_sign(&c_b524, &d__1));
	if ((d__1 = *a - *b - mab, abs(d__1)) < eps && a0 <= 1.1) {
	    *b += eps;
	}
	if ((d__1 = *a - *b - mab, abs(d__1)) > eps) {
	    gamma2_(a, &ga);
	    gamma2_(b, &gb);
	    gamma2_(c__, &gc);
	    d__1 = *a - *b;
	    gamma2_(&d__1, &gab);
	    d__1 = *b - *a;
	    gamma2_(&d__1, &gba);
	    d__1 = *c__ - *a;
	    gamma2_(&d__1, &gca);
	    d__1 = *c__ - *b;
	    gamma2_(&d__1, &gcb);
	    d__1 = gc * gba;
	    z__2.r = d__1, z__2.i = 0.;
	    d__2 = gca * gb;
	    z__5.r = -z__->r, z__5.i = -z__->i;
	    z__6.r = *a, z__6.i = 0.;
	    pow_zz(&z__4, &z__5, &z__6);
	    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
	    z_div(&z__1, &z__2, &z__3);
	    zc0.r = z__1.r, zc0.i = z__1.i;
	    d__1 = gc * gab;
	    z__2.r = d__1, z__2.i = 0.;
	    d__2 = gcb * ga;
	    z__5.r = -z__->r, z__5.i = -z__->i;
	    z__6.r = *b, z__6.i = 0.;
	    pow_zz(&z__4, &z__5, &z__6);
	    z__3.r = d__2 * z__4.r, z__3.i = d__2 * z__4.i;
	    z_div(&z__1, &z__2, &z__3);
	    zc1.r = z__1.r, zc1.i = z__1.i;
	    zr0.r = zc0.r, zr0.i = zc0.i;
	    zr1.r = zc1.r, zr1.i = zc1.i;
	    zhf->r = 0., zhf->i = 0.;
	    for (k = 1; k <= 500; ++k) {
		d__1 = *a + k - 1.;
		z__3.r = d__1 * zr0.r, z__3.i = d__1 * zr0.i;
		d__2 = *a - *c__ + k;
		z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
		d__3 = (*a - *b + k) * k;
		z__4.r = d__3 * z__->r, z__4.i = d__3 * z__->i;
		z_div(&z__1, &z__2, &z__4);
		zr0.r = z__1.r, zr0.i = z__1.i;
		d__1 = *b + k - 1.;
		z__3.r = d__1 * zr1.r, z__3.i = d__1 * zr1.i;
		d__2 = *b - *c__ + k;
		z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
		d__3 = (*b - *a + k) * k;
		z__4.r = d__3 * z__->r, z__4.i = d__3 * z__->i;
		z_div(&z__1, &z__2, &z__4);
		zr1.r = z__1.r, zr1.i = z__1.i;
		z__2.r = zhf->r + zr0.r, z__2.i = zhf->i + zr0.i;
		z__1.r = z__2.r + zr1.r, z__1.i = z__2.i + zr1.i;
		zhf->r = z__1.r, zhf->i = z__1.i;
		z__2.r = zhf->r - zw.r, z__2.i = zhf->i - zw.i;
		z_div(&z__1, &z__2, zhf);
		if (z_abs(&z__1) <= eps) {
		    goto L115;
		}
/* L110: */
		zw.r = zhf->r, zw.i = zhf->i;
	    }
L115:
	    z__2.r = zhf->r + zc0.r, z__2.i = zhf->i + zc0.i;
	    z__1.r = z__2.r + zc1.r, z__1.i = z__2.i + zc1.i;
	    zhf->r = z__1.r, zhf->i = z__1.i;
	} else {
	    if (*a - *b < 0.) {
		*a = bb;
		*b = aa;
	    }
	    ca = *c__ - *a;
	    cb = *c__ - *b;
	    nca = (integer) (ca + eps * d_sign(&c_b524, &ca));
	    ncb = (integer) (cb + eps * d_sign(&c_b524, &cb));
	    if ((d__1 = ca - nca, abs(d__1)) < eps || (d__2 = cb - ncb, abs(
		    d__2)) < eps) {
		*c__ += eps;
	    }
	    gamma2_(a, &ga);
	    gamma2_(c__, &gc);
	    d__1 = *c__ - *b;
	    gamma2_(&d__1, &gcb);
	    psi_spec__(a, &pa);
	    d__1 = *c__ - *a;
	    psi_spec__(&d__1, &pca);
	    d__1 = *a - *c__;
	    psi_spec__(&d__1, &pac);
	    mab = (integer) (*a - *b + eps);
	    z__2.r = gc, z__2.i = 0.;
	    z__5.r = -z__->r, z__5.i = -z__->i;
	    z__6.r = *b, z__6.i = 0.;
	    pow_zz(&z__4, &z__5, &z__6);
	    z__3.r = ga * z__4.r, z__3.i = ga * z__4.i;
	    z_div(&z__1, &z__2, &z__3);
	    zc0.r = z__1.r, zc0.i = z__1.i;
	    d__1 = *a - *b;
	    gamma2_(&d__1, &gm);
	    d__1 = gm / gcb;
	    z__1.r = d__1 * zc0.r, z__1.i = d__1 * zc0.i;
	    zf0.r = z__1.r, zf0.i = z__1.i;
	    zr.r = zc0.r, zr.i = zc0.i;
	    i__1 = mab - 1;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = *b + k - 1.;
		z__2.r = d__1 * zr.r, z__2.i = d__1 * zr.i;
		d__2 = (doublereal) k;
		z__3.r = d__2 * z__->r, z__3.i = d__2 * z__->i;
		z_div(&z__1, &z__2, &z__3);
		zr.r = z__1.r, zr.i = z__1.i;
		t0 = *a - *b - k;
		gamma2_(&t0, &g0);
		d__1 = *c__ - *b - k;
		gamma2_(&d__1, &gcbk);
/* L120: */
		z__3.r = g0 * zr.r, z__3.i = g0 * zr.i;
		z__2.r = z__3.r / gcbk, z__2.i = z__3.i / gcbk;
		z__1.r = zf0.r + z__2.r, z__1.i = zf0.i + z__2.i;
		zf0.r = z__1.r, zf0.i = z__1.i;
	    }
	    if (mab == 0) {
		zf0.r = 0., zf0.i = 0.;
	    }
	    z__2.r = gc, z__2.i = 0.;
	    d__1 = ga * gcb;
	    z__5.r = -z__->r, z__5.i = -z__->i;
	    z__6.r = *a, z__6.i = 0.;
	    pow_zz(&z__4, &z__5, &z__6);
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    z_div(&z__1, &z__2, &z__3);
	    zc1.r = z__1.r, zc1.i = z__1.i;
	    sp = el * -2. - pa - pca;
	    i__1 = mab;
	    for (j = 1; j <= i__1; ++j) {
/* L125: */
		sp += 1. / j;
	    }
	    z__3.r = -z__->r, z__3.i = -z__->i;
	    z_log(&z__2, &z__3);
	    z__1.r = sp + z__2.r, z__1.i = z__2.i;
	    zp0.r = z__1.r, zp0.i = z__1.i;
	    sq = 1.;
	    i__1 = mab;
	    for (j = 1; j <= i__1; ++j) {
/* L130: */
		sq = sq * (*b + j - 1.) * (*b - *c__ + j) / j;
	    }
	    z__2.r = sq * zp0.r, z__2.i = sq * zp0.i;
	    z__1.r = z__2.r * zc1.r - z__2.i * zc1.i, z__1.i = z__2.r * zc1.i 
		    + z__2.i * zc1.r;
	    zf1.r = z__1.r, zf1.i = z__1.i;
	    zr.r = zc1.r, zr.i = zc1.i;
	    rk1 = 1.;
	    sj1 = 0.;
	    w0 = 0.;
	    for (k = 1; k <= 10000; ++k) {
		z_div(&z__1, &zr, z__);
		zr.r = z__1.r, zr.i = z__1.i;
		rk1 = rk1 * (*b + k - 1.) * (*b - *c__ + k) / (k * k);
		rk2 = rk1;
		i__1 = k + mab;
		for (j = k + 1; j <= i__1; ++j) {
/* L135: */
		    rk2 = rk2 * (*b + j - 1.) * (*b - *c__ + j) / j;
		}
		sj1 = sj1 + (*a - 1.) / (k * (*a + k - 1.)) + (*a - *c__ - 1.)
			 / (k * (*a - *c__ + k - 1.));
		sj2 = sj1;
		i__1 = k + mab;
		for (j = k + 1; j <= i__1; ++j) {
/* L140: */
		    sj2 += 1. / j;
		}
		d__1 = el * -2. - pa - pac + sj2 - 1. / (k + *a - *c__) - pi /
			 tan(pi * (k + *a - *c__));
		z__3.r = -z__->r, z__3.i = -z__->i;
		z_log(&z__2, &z__3);
		z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
		zp.r = z__1.r, zp.i = z__1.i;
		z__3.r = rk2 * zr.r, z__3.i = rk2 * zr.i;
		z__2.r = z__3.r * zp.r - z__3.i * zp.i, z__2.i = z__3.r * 
			zp.i + z__3.i * zp.r;
		z__1.r = zf1.r + z__2.r, z__1.i = zf1.i + z__2.i;
		zf1.r = z__1.r, zf1.i = z__1.i;
		ws = z_abs(&zf1);
		if ((d__1 = (ws - w0) / ws, abs(d__1)) < eps) {
		    goto L150;
		}
/* L145: */
		w0 = ws;
	    }
L150:
	    z__1.r = zf0.r + zf1.r, z__1.i = zf0.i + zf1.i;
	    zhf->r = z__1.r, zhf->i = z__1.i;
	}
    }
    *a = aa;
    *b = bb;
    if (k > 150) {
	*isfer = 5;
    }
    return 0;
} /* hygfz_ */

/*       ********************************** */
/* Subroutine */ int itairy_(doublereal *x, doublereal *apt, doublereal *bpt, 
	doublereal *ant, doublereal *bnt)
{
    /* Initialized data */

    static doublereal a[16] = { .569444444444444,.891300154320988,
	    2.26624344493027,7.98950124766861,36.0688546785343,
	    198.670292131169,1292.23456582211,9694.838696696,82418.4704952483,
	    783031.092490225,8222104.93622814,94555739.9360556,
	    1181955956.4073,15956465304.0121,231369166433.05,3586225227969.69 
	    };

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal), exp(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer k, l;
    static doublereal r__, c1, c2, q0, q1, q2, pi, xe, fx, gx, sr3, su1, su2, 
	    xr1, su3, xr2, xp6, su4, su5, su6, eps;


/*       ====================================================== */
/*       Purpose: Compute the integrals of Airy fnctions with */
/*                respect to t from 0 and x ( x ≥ 0 ) */
/*       Input  : x   --- Upper limit of the integral */
/*       Output : APT --- Integration of Ai(t) from 0 and x */
/*                BPT --- Integration of Bi(t) from 0 and x */
/*                ANT --- Integration of Ai(-t) from 0 and x */
/*                BNT --- Integration of Bi(-t) from 0 and x */
/*       ====================================================== */

    eps = 1e-15;
    pi = 3.141592653589793;
    c1 = .355028053887817;
    c2 = .258819403792807;
    sr3 = 1.732050807568877;
    if (*x == 0.) {
	*apt = 0.;
	*bpt = 0.;
	*ant = 0.;
	*bnt = 0.;
    } else {
	if (abs(*x) <= 9.25) {
	    for (l = 0; l <= 1; ++l) {
		*x = pow_ii(&c_n1, &l) * *x;
		fx = *x;
		r__ = *x;
		for (k = 1; k <= 40; ++k) {
		    r__ = r__ * (k * 3.f - 2.) / (k * 3.f + 1.) * *x / (k * 
			    3.f) * *x / (k * 3.f - 1.) * *x;
		    fx += r__;
		    if (abs(r__) < abs(fx) * eps) {
			goto L15;
		    }
/* L10: */
		}
L15:
		gx = *x * .5 * *x;
		r__ = gx;
		for (k = 1; k <= 40; ++k) {
		    r__ = r__ * (k * 3.f - 1.) / (k * 3.f + 2.) * *x / (k * 
			    3.f) * *x / (k * 3.f + 1.) * *x;
		    gx += r__;
		    if (abs(r__) < abs(gx) * eps) {
			goto L25;
		    }
/* L20: */
		}
L25:
		*ant = c1 * fx - c2 * gx;
		*bnt = sr3 * (c1 * fx + c2 * gx);
		if (l == 0) {
		    *apt = *ant;
		    *bpt = *bnt;
		} else {
		    *ant = -(*ant);
		    *bnt = -(*bnt);
		    *x = -(*x);
		}
/* L30: */
	    }
	} else {
	    q2 = 1.414213562373095;
	    q0 = .3333333333333333;
	    q1 = .6666666666666667;
	    xe = *x * sqrt(*x) / 1.5;
	    xp6 = 1. / sqrt(pi * 6. * xe);
	    su1 = 1.;
	    r__ = 1.;
	    xr1 = 1. / xe;
	    for (k = 1; k <= 16; ++k) {
		r__ = -r__ * xr1;
/* L35: */
		su1 += a[k - 1] * r__;
	    }
	    su2 = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 16; ++k) {
		r__ *= xr1;
/* L40: */
		su2 += a[k - 1] * r__;
	    }
	    *apt = q0 - exp(-xe) * xp6 * su1;
	    *bpt = exp(xe) * 2. * xp6 * su2;
	    su3 = 1.;
	    r__ = 1.;
	    xr2 = 1. / (xe * xe);
	    for (k = 1; k <= 8; ++k) {
		r__ = -r__ * xr2;
/* L45: */
		su3 += a[(k << 1) - 1] * r__;
	    }
	    su4 = a[0] * xr1;
	    r__ = xr1;
	    for (k = 1; k <= 7; ++k) {
		r__ = -r__ * xr2;
/* L50: */
		su4 += a[k * 2] * r__;
	    }
	    su5 = su3 + su4;
	    su6 = su3 - su4;
	    *ant = q1 - q2 * xp6 * (su5 * cos(xe) - su6 * sin(xe));
	    *bnt = q2 * xp6 * (su5 * sin(xe) + su6 * cos(xe));
	}
    }
    return 0;
} /* itairy_ */

/*       ********************************** */
/* Subroutine */ int ikna_(integer *n, doublereal *x, integer *nm, doublereal 
	*bi, doublereal *di, doublereal *bk, doublereal *dk)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal f, g, h__;
    static integer k, m;
    static doublereal f0, f1, h0, h1, g0, g1, s0, bi0, bi1, di0, di1, bk0, 
	    dk0, bk1, dk1;
    extern /* Subroutine */ int ik01a_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ======================================================== */
/*       Purpose: Compute modified Bessel functions In(x) and */
/*                Kn(x), and their derivatives */
/*       Input:   x --- Argument of In(x) and Kn(x) ( x ≥ 0 ) */
/*                n --- Order of In(x) and Kn(x) */
/*       Output:  BI(n) --- In(x) */
/*                DI(n) --- In'(x) */
/*                BK(n) --- Kn(x) */
/*                DK(n) --- Kn'(x) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*            (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x) */
/*            (2) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       ======================================================== */

    *nm = *n;
    if (*x <= 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    bi[k] = 0.;
	    di[k] = 0.;
	    bk[k] = 1e300;
/* L10: */
	    dk[k] = -1e300;
	}
	bi[0] = 1.;
	di[1] = .5;
	return 0;
    }
    ik01a_(x, &bi0, &di0, &bi1, &di1, &bk0, &dk0, &bk1, &dk1);
    bi[0] = bi0;
    bi[1] = bi1;
    bk[0] = bk0;
    bk[1] = bk1;
    di[0] = di0;
    di[1] = di1;
    dk[0] = dk0;
    dk[1] = dk1;
    if (*n <= 1) {
	return 0;
    }
    if (*x > 40.f && *n < (integer) (*x * .25f)) {
	h0 = bi0;
	h1 = bi1;
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    h__ = (k - 1.) * -2. / *x * h1 + h0;
	    bi[k] = h__;
	    h0 = h1;
/* L15: */
	    h1 = h__;
	}
    } else {
	m = msta1_(x, &c__200);
	if (m < *n) {
	    *nm = m;
	} else {
	    m = msta2_(x, n, &c__15);
	}
	f0 = 0.;
	f1 = 1e-100;
	f = 0.;
	for (k = m; k >= 0; --k) {
	    f = (k + 1.) * 2. * f1 / *x + f0;
	    if (k <= *nm) {
		bi[k] = f;
	    }
	    f0 = f1;
/* L20: */
	    f1 = f;
	}
	s0 = bi0 / f;
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
/* L25: */
	    bi[k] = s0 * bi[k];
	}
    }
    g0 = bk0;
    g1 = bk1;
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
	g = (k - 1.) * 2. / *x * g1 + g0;
	bk[k] = g;
	g0 = g1;
/* L30: */
	g1 = g;
    }
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
	di[k] = bi[k - 1] - k / *x * bi[k];
/* L40: */
	dk[k] = -bk[k - 1] - k / *x * bk[k];
    }
    return 0;
} /* ikna_ */

/*       ********************************** */
/* Subroutine */ int cjynb_(integer *n, doublecomplex *z__, integer *nm, 
	doublecomplex *cbj, doublecomplex *cdj, doublecomplex *cby, 
	doublecomplex *cdy)
{
    /* Initialized data */

    static doublereal a[4] = { -.0703125,.112152099609375,-.5725014209747314,
	    6.074042001273483 };
    static doublereal b[4] = { .0732421875,-.2271080017089844,
	    1.727727502584457,-24.38052969955606 };
    static doublereal a1[4] = { .1171875,-.144195556640625,.6765925884246826,
	    -6.883914268109947 };
    static doublereal b1[4] = { -.1025390625,.2775764465332031,
	    -1.993531733751297,27.24882731126854 };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    integer pow_ii(integer *, integer *);
    void z_cos(doublecomplex *, doublecomplex *), z_log(doublecomplex *, 
	    doublecomplex *), pow_zi(doublecomplex *, doublecomplex *, 
	    integer *), z_sqrt(doublecomplex *, doublecomplex *), z_sin(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, m;
    static doublereal a0, y0;
    static doublecomplex ce, cf;
    static doublereal el;
    static doublecomplex cu;
    static doublereal pi;
    static doublecomplex cf1, cf2, cp0, cq0, cp1, cs0, cq1, ct1, ct2;
    static doublereal r2p;
    static doublecomplex cbs, csu, csv, cyy, cbj0, cbj1, cby0, cby1, cbjk;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ======================================================= */
/*       Purpose: Compute Bessel functions Jn(z), Yn(z) and */
/*                their derivatives for a complex argument */
/*       Input :  z --- Complex argument of Jn(z) and Yn(z) */
/*                n --- Order of Jn(z) and Yn(z) */
/*       Output:  CBJ(n) --- Jn(z) */
/*                CDJ(n) --- Jn'(z) */
/*                CBY(n) --- Yn(z) */
/*                CDY(n) --- Yn'(z) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 to calculate the starting */
/*                point for backward recurrence */
/*       ======================================================= */

    el = .5772156649015329;
    pi = 3.141592653589793;
    r2p = .63661977236758;
    y0 = (d__1 = d_imag(z__), abs(d__1));
    a0 = z_abs(z__);
    *nm = *n;
    if (a0 < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbj[i__2].r = 0., cbj[i__2].i = 0.;
	    i__2 = k;
	    cdj[i__2].r = 0., cdj[i__2].i = 0.;
	    i__2 = k;
	    cby[i__2].r = -1e300, cby[i__2].i = -0.;
/* L10: */
	    i__2 = k;
	    cdy[i__2].r = 1e300, cdy[i__2].i = 0.;
	}
	cbj[0].r = 1., cbj[0].i = 0.;
	cdj[1].r = .5, cdj[1].i = 0.;
	return 0;
    }
    if (a0 <= 300. || *n > 80) {
	if (*n == 0) {
	    *nm = 1;
	}
	m = msta1_(&a0, &c__200);
	if (m < *nm) {
	    *nm = m;
	} else {
	    m = msta2_(&a0, nm, &c__15);
	}
	cbs.r = 0., cbs.i = 0.;
	csu.r = 0., csu.i = 0.;
	csv.r = 0., csv.i = 0.;
	cf2.r = 0., cf2.i = 0.;
	cf1.r = 1e-100, cf1.i = 0.;
	for (k = m; k >= 0; --k) {
	    d__1 = (k + 1.) * 2.;
	    z__4.r = d__1, z__4.i = 0.;
	    z_div(&z__3, &z__4, z__);
	    z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i 
		    + z__3.i * cf1.r;
	    z__1.r = z__2.r - cf2.r, z__1.i = z__2.i - cf2.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	    if (k <= *nm) {
		i__2 = k;
		cbj[i__2].r = cf.r, cbj[i__2].i = cf.i;
	    }
	    if (k == k / 2 << 1 && k != 0) {
		if (y0 <= 1.) {
		    z__2.r = cf.r * 2., z__2.i = cf.i * 2.;
		    z__1.r = cbs.r + z__2.r, z__1.i = cbs.i + z__2.i;
		    cbs.r = z__1.r, cbs.i = z__1.i;
		} else {
		    i__2 = k / 2;
		    d__1 = pow_ii(&c_n1, &i__2) * 2.;
		    z__2.r = d__1 * cf.r, z__2.i = d__1 * cf.i;
		    z__1.r = cbs.r + z__2.r, z__1.i = cbs.i + z__2.i;
		    cbs.r = z__1.r, cbs.i = z__1.i;
		}
		i__1 = k / 2;
		i__2 = pow_ii(&c_n1, &i__1);
		d__1 = (doublereal) i__2;
		z__3.r = d__1 * cf.r, z__3.i = d__1 * cf.i;
		d__2 = (doublereal) k;
		z__2.r = z__3.r / d__2, z__2.i = z__3.i / d__2;
		z__1.r = csu.r + z__2.r, z__1.i = csu.i + z__2.i;
		csu.r = z__1.r, csu.i = z__1.i;
	    } else if (k > 1) {
		i__2 = k / 2;
		d__1 = pow_ii(&c_n1, &i__2) * k / (k * k - 1.);
		z__2.r = d__1 * cf.r, z__2.i = d__1 * cf.i;
		z__1.r = csv.r + z__2.r, z__1.i = csv.i + z__2.i;
		csv.r = z__1.r, csv.i = z__1.i;
	    }
	    cf2.r = cf1.r, cf2.i = cf1.i;
/* L15: */
	    cf1.r = cf.r, cf1.i = cf.i;
	}
	if (y0 <= 1.) {
	    z__1.r = cbs.r + cf.r, z__1.i = cbs.i + cf.i;
	    cs0.r = z__1.r, cs0.i = z__1.i;
	} else {
	    z__2.r = cbs.r + cf.r, z__2.i = cbs.i + cf.i;
	    z_cos(&z__3, z__);
	    z_div(&z__1, &z__2, &z__3);
	    cs0.r = z__1.r, cs0.i = z__1.i;
	}
	i__2 = *nm;
	for (k = 0; k <= i__2; ++k) {
/* L20: */
	    i__1 = k;
	    z_div(&z__1, &cbj[k], &cs0);
	    cbj[i__1].r = z__1.r, cbj[i__1].i = z__1.i;
	}
	z__3.r = z__->r / 2., z__3.i = z__->i / 2.;
	z_log(&z__2, &z__3);
	z__1.r = z__2.r + el, z__1.i = z__2.i;
	ce.r = z__1.r, ce.i = z__1.i;
	z__3.r = ce.r * cbj[0].r - ce.i * cbj[0].i, z__3.i = ce.r * cbj[0].i 
		+ ce.i * cbj[0].r;
	z__5.r = csu.r * 4., z__5.i = csu.i * 4.;
	z_div(&z__4, &z__5, &cs0);
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z__1.r = r2p * z__2.r, z__1.i = r2p * z__2.i;
	cby[0].r = z__1.r, cby[0].i = z__1.i;
	z__5.r = -cbj[0].r, z__5.i = -cbj[0].i;
	z_div(&z__4, &z__5, z__);
	z__7.r = ce.r - 1., z__7.i = ce.i;
	z__6.r = z__7.r * cbj[1].r - z__7.i * cbj[1].i, z__6.i = z__7.r * cbj[
		1].i + z__7.i * cbj[1].r;
	z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
	z__9.r = csv.r * 4., z__9.i = csv.i * 4.;
	z_div(&z__8, &z__9, &cs0);
	z__2.r = z__3.r - z__8.r, z__2.i = z__3.i - z__8.i;
	z__1.r = r2p * z__2.r, z__1.i = r2p * z__2.i;
	cby[1].r = z__1.r, cby[1].i = z__1.i;
    } else {
	d__1 = pi * .25;
	z__1.r = z__->r - d__1, z__1.i = z__->i;
	ct1.r = z__1.r, ct1.i = z__1.i;
	cp0.r = 1., cp0.i = 0.;
	for (k = 1; k <= 4; ++k) {
/* L25: */
	    i__1 = k - 1;
	    i__2 = k * -2;
	    pow_zi(&z__3, z__, &i__2);
	    z__2.r = a[i__1] * z__3.r, z__2.i = a[i__1] * z__3.i;
	    z__1.r = cp0.r + z__2.r, z__1.i = cp0.i + z__2.i;
	    cp0.r = z__1.r, cp0.i = z__1.i;
	}
	z_div(&z__1, &c_b772, z__);
	cq0.r = z__1.r, cq0.i = z__1.i;
	for (k = 1; k <= 4; ++k) {
/* L30: */
	    i__1 = k - 1;
	    i__2 = k * -2 - 1;
	    pow_zi(&z__3, z__, &i__2);
	    z__2.r = b[i__1] * z__3.r, z__2.i = b[i__1] * z__3.i;
	    z__1.r = cq0.r + z__2.r, z__1.i = cq0.i + z__2.i;
	    cq0.r = z__1.r, cq0.i = z__1.i;
	}
	z__3.r = r2p, z__3.i = 0.;
	z_div(&z__2, &z__3, z__);
	z_sqrt(&z__1, &z__2);
	cu.r = z__1.r, cu.i = z__1.i;
	z_cos(&z__4, &ct1);
	z__3.r = cp0.r * z__4.r - cp0.i * z__4.i, z__3.i = cp0.r * z__4.i + 
		cp0.i * z__4.r;
	z_sin(&z__6, &ct1);
	z__5.r = cq0.r * z__6.r - cq0.i * z__6.i, z__5.i = cq0.r * z__6.i + 
		cq0.i * z__6.r;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cbj0.r = z__1.r, cbj0.i = z__1.i;
	z_sin(&z__4, &ct1);
	z__3.r = cp0.r * z__4.r - cp0.i * z__4.i, z__3.i = cp0.r * z__4.i + 
		cp0.i * z__4.r;
	z_cos(&z__6, &ct1);
	z__5.r = cq0.r * z__6.r - cq0.i * z__6.i, z__5.i = cq0.r * z__6.i + 
		cq0.i * z__6.r;
	z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cby0.r = z__1.r, cby0.i = z__1.i;
	cbj[0].r = cbj0.r, cbj[0].i = cbj0.i;
	cby[0].r = cby0.r, cby[0].i = cby0.i;
	d__1 = pi * .75;
	z__1.r = z__->r - d__1, z__1.i = z__->i;
	ct2.r = z__1.r, ct2.i = z__1.i;
	cp1.r = 1., cp1.i = 0.;
	for (k = 1; k <= 4; ++k) {
/* L35: */
	    i__1 = k - 1;
	    i__2 = k * -2;
	    pow_zi(&z__3, z__, &i__2);
	    z__2.r = a1[i__1] * z__3.r, z__2.i = a1[i__1] * z__3.i;
	    z__1.r = cp1.r + z__2.r, z__1.i = cp1.i + z__2.i;
	    cp1.r = z__1.r, cp1.i = z__1.i;
	}
	z_div(&z__1, &c_b776, z__);
	cq1.r = z__1.r, cq1.i = z__1.i;
	for (k = 1; k <= 4; ++k) {
/* L40: */
	    i__1 = k - 1;
	    i__2 = k * -2 - 1;
	    pow_zi(&z__3, z__, &i__2);
	    z__2.r = b1[i__1] * z__3.r, z__2.i = b1[i__1] * z__3.i;
	    z__1.r = cq1.r + z__2.r, z__1.i = cq1.i + z__2.i;
	    cq1.r = z__1.r, cq1.i = z__1.i;
	}
	z_cos(&z__4, &ct2);
	z__3.r = cp1.r * z__4.r - cp1.i * z__4.i, z__3.i = cp1.r * z__4.i + 
		cp1.i * z__4.r;
	z_sin(&z__6, &ct2);
	z__5.r = cq1.r * z__6.r - cq1.i * z__6.i, z__5.i = cq1.r * z__6.i + 
		cq1.i * z__6.r;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cbj1.r = z__1.r, cbj1.i = z__1.i;
	z_sin(&z__4, &ct2);
	z__3.r = cp1.r * z__4.r - cp1.i * z__4.i, z__3.i = cp1.r * z__4.i + 
		cp1.i * z__4.r;
	z_cos(&z__6, &ct2);
	z__5.r = cq1.r * z__6.r - cq1.i * z__6.i, z__5.i = cq1.r * z__6.i + 
		cq1.i * z__6.r;
	z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cby1.r = z__1.r, cby1.i = z__1.i;
	cbj[1].r = cbj1.r, cbj[1].i = cbj1.i;
	cby[1].r = cby1.r, cby[1].i = cby1.i;
	i__1 = *nm;
	for (k = 2; k <= i__1; ++k) {
	    d__1 = (k - 1.) * 2.;
	    z__4.r = d__1, z__4.i = 0.;
	    z_div(&z__3, &z__4, z__);
	    z__2.r = z__3.r * cbj1.r - z__3.i * cbj1.i, z__2.i = z__3.r * 
		    cbj1.i + z__3.i * cbj1.r;
	    z__1.r = z__2.r - cbj0.r, z__1.i = z__2.i - cbj0.i;
	    cbjk.r = z__1.r, cbjk.i = z__1.i;
	    i__2 = k;
	    cbj[i__2].r = cbjk.r, cbj[i__2].i = cbjk.i;
	    cbj0.r = cbj1.r, cbj0.i = cbj1.i;
/* L45: */
	    cbj1.r = cbjk.r, cbj1.i = cbjk.i;
	}
    }
    z__1.r = -cbj[1].r, z__1.i = -cbj[1].i;
    cdj[0].r = z__1.r, cdj[0].i = z__1.i;
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L50: */
	i__2 = k;
	i__3 = k - 1;
	z__4.r = (doublereal) k, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__4 = k;
	z__2.r = z__3.r * cbj[i__4].r - z__3.i * cbj[i__4].i, z__2.i = z__3.r 
		* cbj[i__4].i + z__3.i * cbj[i__4].r;
	z__1.r = cbj[i__3].r - z__2.r, z__1.i = cbj[i__3].i - z__2.i;
	cdj[i__2].r = z__1.r, cdj[i__2].i = z__1.i;
    }
    if (z_abs(cbj) > 1.) {
	z__3.r = cbj[1].r * cby[0].r - cbj[1].i * cby[0].i, z__3.i = cbj[1].r 
		* cby[0].i + cbj[1].i * cby[0].r;
	z__5.r = pi * z__->r, z__5.i = pi * z__->i;
	z_div(&z__4, &c_b15, &z__5);
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z_div(&z__1, &z__2, cbj);
	cby[1].r = z__1.r, cby[1].i = z__1.i;
    }
    i__2 = *nm;
    for (k = 2; k <= i__2; ++k) {
	if (z_abs(&cbj[k - 1]) >= z_abs(&cbj[k - 2])) {
	    i__3 = k;
	    i__4 = k - 1;
	    z__3.r = cbj[i__3].r * cby[i__4].r - cbj[i__3].i * cby[i__4].i, 
		    z__3.i = cbj[i__3].r * cby[i__4].i + cbj[i__3].i * cby[
		    i__4].r;
	    z__5.r = pi * z__->r, z__5.i = pi * z__->i;
	    z_div(&z__4, &c_b15, &z__5);
	    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	    z_div(&z__1, &z__2, &cbj[k - 1]);
	    cyy.r = z__1.r, cyy.i = z__1.i;
	} else {
	    i__3 = k;
	    i__4 = k - 2;
	    z__3.r = cbj[i__3].r * cby[i__4].r - cbj[i__3].i * cby[i__4].i, 
		    z__3.i = cbj[i__3].r * cby[i__4].i + cbj[i__3].i * cby[
		    i__4].r;
	    d__1 = (k - 1.) * 4.;
	    z__5.r = d__1, z__5.i = 0.;
	    z__7.r = pi * z__->r, z__7.i = pi * z__->i;
	    z__6.r = z__7.r * z__->r - z__7.i * z__->i, z__6.i = z__7.r * 
		    z__->i + z__7.i * z__->r;
	    z_div(&z__4, &z__5, &z__6);
	    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	    z_div(&z__1, &z__2, &cbj[k - 2]);
	    cyy.r = z__1.r, cyy.i = z__1.i;
	}
	i__3 = k;
	cby[i__3].r = cyy.r, cby[i__3].i = cyy.i;
/* L55: */
    }
    z__1.r = -cby[1].r, z__1.i = -cby[1].i;
    cdy[0].r = z__1.r, cdy[0].i = z__1.i;
    i__2 = *nm;
    for (k = 1; k <= i__2; ++k) {
/* L60: */
	i__3 = k;
	i__4 = k - 1;
	z__4.r = (doublereal) k, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__1 = k;
	z__2.r = z__3.r * cby[i__1].r - z__3.i * cby[i__1].i, z__2.i = z__3.r 
		* cby[i__1].i + z__3.i * cby[i__1].r;
	z__1.r = cby[i__4].r - z__2.r, z__1.i = cby[i__4].i - z__2.i;
	cdy[i__3].r = z__1.r, cdy[i__3].i = z__1.i;
    }
    return 0;
} /* cjynb_ */

/*       ********************************** */
/* Subroutine */ int iknb_(integer *n, doublereal *x, integer *nm, doublereal 
	*bi, doublereal *di, doublereal *bk, doublereal *dk)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double exp(doublereal), log(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal f, g;
    static integer k, l, m;
    static doublereal r__, a0, f0, f1, g0, g1;
    static integer k0;
    static doublereal s0, el, bs, pi, vt, sk0, bkl;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ============================================================ */
/*       Purpose: Compute modified Bessel functions In(x) and Kn(x), */
/*                and their derivatives */
/*       Input:   x --- Argument of In(x) and Kn(x) ( 0 ≤ x ≤ 700 ) */
/*                n --- Order of In(x) and Kn(x) */
/*       Output:  BI(n) --- In(x) */
/*                DI(n) --- In'(x) */
/*                BK(n) --- Kn(x) */
/*                DK(n) --- Kn'(x) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 for computing the starting point */
/*                for backward recurrence */
/*       =========================================================== */

    pi = 3.141592653589793;
    el = .5772156649015329;
    *nm = *n;
    if (*x <= 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    bi[k] = 0.;
	    di[k] = 0.;
	    bk[k] = 1e300;
/* L10: */
	    dk[k] = -1e300;
	}
	bi[0] = 1.;
	di[1] = .5;
	return 0;
    }
    if (*n == 0) {
	*nm = 1;
    }
    m = msta1_(x, &c__200);
    if (m < *nm) {
	*nm = m;
    } else {
	m = msta2_(x, nm, &c__15);
    }
    bs = 0.;
    sk0 = 0.;
    f = 0.;
    f0 = 0.;
    f1 = 1e-100;
    for (k = m; k >= 0; --k) {
	f = (k + 1.) * 2. / *x * f1 + f0;
	if (k <= *nm) {
	    bi[k] = f;
	}
	if (k != 0 && k == k / 2 << 1) {
	    sk0 += f * 4. / k;
	}
	bs += f * 2.;
	f0 = f1;
/* L15: */
	f1 = f;
    }
    s0 = exp(*x) / (bs - f);
    i__1 = *nm;
    for (k = 0; k <= i__1; ++k) {
/* L20: */
	bi[k] = s0 * bi[k];
    }
    if (*x <= 8.) {
	bk[0] = -(log(*x * .5) + el) * bi[0] + s0 * sk0;
	bk[1] = (1. / *x - bi[1] * bk[0]) / bi[0];
    } else {
	a0 = sqrt(pi / (*x * 2.)) * exp(-(*x));
	k0 = 16;
	if (*x >= 25.f) {
	    k0 = 10;
	}
	if (*x >= 80.f) {
	    k0 = 8;
	}
	if (*x >= 200.f) {
	    k0 = 6;
	}
	for (l = 0; l <= 1; ++l) {
	    bkl = 1.;
	    vt = l * 4.;
	    r__ = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
/* Computing 2nd power */
		r__1 = k * 2.f - 1.f;
		r__ = r__ * .125 * (vt - r__1 * r__1) / (k * *x);
/* L25: */
		bkl += r__;
	    }
	    bk[l] = a0 * bkl;
/* L30: */
	}
    }
    g0 = bk[0];
    g1 = bk[1];
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
	g = (k - 1.) * 2. / *x * g1 + g0;
	bk[k] = g;
	g0 = g1;
/* L35: */
	g1 = g;
    }
    di[0] = bi[1];
    dk[0] = -bk[1];
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
	di[k] = bi[k - 1] - k / *x * bi[k];
/* L40: */
	dk[k] = -bk[k - 1] - k / *x * bk[k];
    }
    return 0;
} /* iknb_ */

/*       ********************************** */
/* Subroutine */ int lpmn_(integer *mm, integer *m, integer *n, doublereal *x,
	 doublereal *pm, doublereal *pd)
{
    /* System generated locals */
    integer pm_dim1, pm_offset, pd_dim1, pd_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, ls;
    static doublereal xq, xs;
    extern doublereal dinf_(void);


/*       ===================================================== */
/*       Purpose: Compute the associated Legendre functions */
/*                Pmn(x) and their derivatives Pmn'(x) for */
/*                real argument */
/*       Input :  x  --- Argument of Pmn(x) */
/*                m  --- Order of Pmn(x),  m = 0,1,2,...,n */
/*                n  --- Degree of Pmn(x), n = 0,1,2,...,N */
/*                mm --- Physical dimension of PM and PD */
/*       Output:  PM(m,n) --- Pmn(x) */
/*                PD(m,n) --- Pmn'(x) */
/*       ===================================================== */

    /* Parameter adjustments */
    pd_dim1 = *mm - 0 + 1;
    pd_offset = 0 + pd_dim1 * 0;
    pd -= pd_offset;
    pm_dim1 = *mm - 0 + 1;
    pm_offset = 0 + pm_dim1 * 0;
    pm -= pm_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 0; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 0; j <= i__2; ++j) {
	    pm[j + i__ * pm_dim1] = 0.;
/* L10: */
	    pd[j + i__ * pd_dim1] = 0.;
	}
    }
    pm[0] = 1.;
    if (*n == 0) {
	return 0;
    }
    if (abs(*x) == 1.) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    pm[i__ * pm_dim1] = pow_di(x, &i__);
/* L15: */
	    i__1 = i__ + 1;
	    pd[i__ * pd_dim1] = i__ * .5 * (i__ + 1.) * pow_di(x, &i__1);
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (i__ == 1) {
		    pd[i__ + j * pd_dim1] = dinf_();
		} else if (i__ == 2) {
		    i__3 = j + 1;
		    pd[i__ + j * pd_dim1] = (j + 2) * -.25 * (j + 1) * j * (j 
			    - 1) * pow_di(x, &i__3);
		}
/* L20: */
	    }
	}
	return 0;
    }
    ls = 1;
    if (abs(*x) > 1.) {
	ls = -1;
    }
    xq = sqrt(ls * (1. - *x * *x));
/*       Ensure connection to the complex-valued function for |x| > 1 */
    if (*x < -1.) {
	xq = -xq;
    }
    xs = ls * (1. - *x * *x);
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L30: */
	pm[i__ + i__ * pm_dim1] = -ls * (i__ * 2. - 1.) * xq * pm[i__ - 1 + (
		i__ - 1) * pm_dim1];
    }
/* Computing MIN */
    i__1 = *m, i__3 = *n - 1;
    i__2 = min(i__1,i__3);
    for (i__ = 0; i__ <= i__2; ++i__) {
/* L35: */
	pm[i__ + (i__ + 1) * pm_dim1] = (i__ * 2. + 1.) * *x * pm[i__ + i__ * 
		pm_dim1];
    }
    i__2 = *m;
    for (i__ = 0; i__ <= i__2; ++i__) {
	i__1 = *n;
	for (j = i__ + 2; j <= i__1; ++j) {
	    pm[i__ + j * pm_dim1] = ((j * 2. - 1.) * *x * pm[i__ + (j - 1) * 
		    pm_dim1] - (i__ + j - 1.) * pm[i__ + (j - 2) * pm_dim1]) /
		     (j - i__);
/* L40: */
	}
    }
    pd[0] = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L45: */
	pd[j * pd_dim1] = ls * j * (pm[(j - 1) * pm_dim1] - *x * pm[j * 
		pm_dim1]) / xs;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    pd[i__ + j * pd_dim1] = ls * i__ * *x * pm[i__ + j * pm_dim1] / 
		    xs + (j + i__) * (j - i__ + 1.) / xq * pm[i__ - 1 + j * 
		    pm_dim1];
/* L50: */
	}
    }
    return 0;
} /* lpmn_ */

/*       ********************************** */
/* Subroutine */ int mtu0_(integer *kf, integer *m, doublereal *q, doublereal 
	*x, doublereal *csf, doublereal *csd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal a;
    static integer k, ic;
    static doublereal fg[251];
    static integer kd;
    static doublereal rd;
    static integer km;
    static doublereal qm, xr, eps;
    extern /* Subroutine */ int cva2_(integer *, integer *, doublereal *, 
	    doublereal *);
    extern doublereal dnan_(void);
    extern /* Subroutine */ int fcoef_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);


/*       =============================================================== */
/*       Purpose: Compute Mathieu functions cem(x,q) and sem(x,q) */
/*                and their derivatives ( q ≥ 0 ) */
/*       Input :  KF  --- Function code */
/*                        KF=1 for computing cem(x,q) and cem'(x,q) */
/*                        KF=2 for computing sem(x,q) and sem'(x,q) */
/*                m   --- Order of Mathieu functions */
/*                q   --- Parameter of Mathieu functions */
/*                x   --- Argument of Mathieu functions (in degrees) */
/*       Output:  CSF --- cem(x,q) or sem(x,q) */
/*                CSD --- cem'x,q) or sem'x,q) */
/*       Routines called: */
/*            (1) CVA2 for computing the characteristic values */
/*            (2) FCOEF for computing the expansion coefficients */
/*       =============================================================== */

    eps = 1e-14;
    if (*kf == 1 && *m == *m / 2 << 1) {
	kd = 1;
    }
    if (*kf == 1 && *m != *m / 2 << 1) {
	kd = 2;
    }
    if (*kf == 2 && *m != *m / 2 << 1) {
	kd = 3;
    }
    if (*kf == 2 && *m == *m / 2 << 1) {
	kd = 4;
    }
    cva2_(&kd, m, q, &a);
    if (*q <= 1.) {
	qm = sqrt(*q) * 56.1f + 7.5f - *q * 134.7f + sqrt(*q) * 90.7f * *q;
    } else {
	qm = sqrt(*q) * 3.1f + 17.f - *q * .126f + sqrt(*q) * .0037f * *q;
    }
    km = (integer) (qm + *m * .5f);
    if (km > 251) {
	*csf = dnan_();
	*csd = dnan_();
	return 0;
    }
    fcoef_(&kd, m, q, &a, fg);
    ic = *m / 2 + 1;
    rd = .0174532925199433;
    xr = *x * rd;
    *csf = 0.;
    i__1 = km;
    for (k = 1; k <= i__1; ++k) {
	if (kd == 1) {
	    *csf += fg[k - 1] * cos(((k << 1) - 2) * xr);
	} else if (kd == 2) {
	    *csf += fg[k - 1] * cos(((k << 1) - 1) * xr);
	} else if (kd == 3) {
	    *csf += fg[k - 1] * sin(((k << 1) - 1) * xr);
	} else if (kd == 4) {
	    *csf += fg[k - 1] * sin((k << 1) * xr);
	}
	if (k >= ic && (d__1 = fg[k - 1], abs(d__1)) < abs(*csf) * eps) {
	    goto L15;
	}
/* L10: */
    }
L15:
    *csd = 0.;
    i__1 = km;
    for (k = 1; k <= i__1; ++k) {
	if (kd == 1) {
	    *csd -= ((k << 1) - 2) * fg[k - 1] * sin(((k << 1) - 2) * xr);
	} else if (kd == 2) {
	    *csd -= ((k << 1) - 1) * fg[k - 1] * sin(((k << 1) - 1) * xr);
	} else if (kd == 3) {
	    *csd += ((k << 1) - 1) * fg[k - 1] * cos(((k << 1) - 1) * xr);
	} else if (kd == 4) {
	    *csd += k * 2. * fg[k - 1] * cos((k << 1) * xr);
	}
	if (k >= ic && (d__1 = fg[k - 1], abs(d__1)) < abs(*csd) * eps) {
	    goto L25;
	}
/* L20: */
    }
L25:
    return 0;
} /* mtu0_ */

/*       ********************************** */
/* Subroutine */ int cy01_(integer *kf, doublecomplex *z__, doublecomplex *zf,
	 doublecomplex *zd)
{
    /* Initialized data */

    static doublereal a[12] = { -.0703125,.112152099609375,-.5725014209747314,
	    6.074042001273483,-110.0171402692467,3038.090510922384,
	    -118838.4262567832,6252951.493434797,-425939216.5047669,
	    36468400807.06556,-3833534661393.944,485401468685290.1 };
    static doublereal b[12] = { .0732421875,-.2271080017089844,
	    1.727727502584457,-24.38052969955606,551.3358961220206,
	    -18257.75547429318,832859.3040162893,-50069589.53198893,
	    3836255180.230433,-364901081884.9833,42189715702840.96,
	    -5827244631566907. };
    static doublereal a1[12] = { .1171875,-.144195556640625,.6765925884246826,
	    -6.883914268109947,121.5978918765359,-3302.272294480852,
	    127641.2726461746,-6656367.718817688,450278600.3050393,
	    -38338575207.4279,4011838599133.198,-506056850331472.7 };
    static doublereal b1[12] = { -.1025390625,.2775764465332031,
	    -1.993531733751297,27.24882731126854,-603.8440767050702,
	    19718.37591223663,-890297.8767070678,53104110.10968522,
	    -4043620325.107754,382701134659.8605,-44064814178522.78,
	    6065091351222699. };

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9, z__10;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_log(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *), pow_zi(doublecomplex *, 
	    doublecomplex *, integer *), z_sqrt(doublecomplex *, 
	    doublecomplex *), z_cos(doublecomplex *, doublecomplex *), z_sin(
	    doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer k;
    static doublereal a0;
    static integer k0;
    static doublereal w0, w1;
    static doublecomplex z1, z2, ci;
    static doublereal el;
    static doublecomplex cp, cr, cs, cu;
    static doublereal pi;
    static doublecomplex cp0, cq0, cp1, cq1, ct1, ct2;
    static doublereal rp2;
    static doublecomplex cbj0, cbj1, cby0, cby1, cdy0, cdy1;


/*       =========================================================== */
/*       Purpose: Compute complex Bessel functions Y0(z), Y1(z) */
/*                and their derivatives */
/*       Input :  z  --- Complex argument of Yn(z) ( n=0,1 ) */
/*                KF --- Function choice code */
/*                    KF=0 for ZF=Y0(z) and ZD=Y0'(z) */
/*                    KF=1 for ZF=Y1(z) and ZD=Y1'(z) */
/*                    KF=2 for ZF=Y1'(z) and ZD=Y1''(z) */
/*       Output:  ZF --- Y0(z) or Y1(z) or Y1'(z) */
/*                ZD --- Y0'(z) or Y1'(z) or Y1''(z) */
/*       =========================================================== */

    pi = 3.141592653589793;
    el = .5772156649015329;
    rp2 = 2. / pi;
    ci.r = 0., ci.i = 1.;
    a0 = z_abs(z__);
    z__1.r = z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = z__1.r, z2.i = z__1.i;
    z1.r = z__->r, z1.i = z__->i;
    if (a0 == 0.) {
	cbj0.r = 1., cbj0.i = 0.;
	cbj1.r = 0., cbj1.i = 0.;
	cby0.r = -1e300, cby0.i = -0.;
	cby1.r = -1e300, cby1.i = -0.;
	cdy0.r = 1e300, cdy0.i = 0.;
	cdy1.r = 1e300, cdy1.i = 0.;
	goto L70;
    }
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (a0 <= 12.f) {
	cbj0.r = 1., cbj0.i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 40; ++k) {
	    z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    i__1 = k * k;
	    d__1 = (doublereal) i__1;
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cbj0.r + cr.r, z__1.i = cbj0.i + cr.i;
	    cbj0.r = z__1.r, cbj0.i = z__1.i;
	    if (z_abs(&cr) < z_abs(&cbj0) * 1e-15) {
		goto L15;
	    }
/* L10: */
	}
L15:
	cbj1.r = 1., cbj1.i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 40; ++k) {
	    z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    d__1 = k * (k + 1.);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cbj1.r + cr.r, z__1.i = cbj1.i + cr.i;
	    cbj1.r = z__1.r, cbj1.i = z__1.i;
	    if (z_abs(&cr) < z_abs(&cbj1) * 1e-15) {
		goto L25;
	    }
/* L20: */
	}
L25:
	z__2.r = z1.r * .5, z__2.i = z1.i * .5;
	z__1.r = z__2.r * cbj1.r - z__2.i * cbj1.i, z__1.i = z__2.r * cbj1.i 
		+ z__2.i * cbj1.r;
	cbj1.r = z__1.r, cbj1.i = z__1.i;
	w0 = 0.;
	cr.r = 1., cr.i = 0.;
	cs.r = 0., cs.i = 0.;
	for (k = 1; k <= 40; ++k) {
	    w0 += 1. / k;
	    z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
	    i__1 = k * k;
	    d__1 = (doublereal) i__1;
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
	    z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * z2.i + 
		    z__2.i * z2.r;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = w0 * cr.r, z__1.i = w0 * cr.i;
	    cp.r = z__1.r, cp.i = z__1.i;
	    z__1.r = cs.r + cp.r, z__1.i = cs.i + cp.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	    if (z_abs(&cp) < z_abs(&cs) * 1e-15) {
		goto L35;
	    }
/* L30: */
	}
L35:
	z__6.r = z1.r / 2., z__6.i = z1.i / 2.;
	z_log(&z__5, &z__6);
	z__4.r = z__5.r + el, z__4.i = z__5.i;
	z__3.r = rp2 * z__4.r, z__3.i = rp2 * z__4.i;
	z__2.r = z__3.r * cbj0.r - z__3.i * cbj0.i, z__2.i = z__3.r * cbj0.i 
		+ z__3.i * cbj0.r;
	z__7.r = rp2 * cs.r, z__7.i = rp2 * cs.i;
	z__1.r = z__2.r - z__7.r, z__1.i = z__2.i - z__7.i;
	cby0.r = z__1.r, cby0.i = z__1.i;
	w1 = 0.;
	cr.r = 1., cr.i = 0.;
	cs.r = 1., cs.i = 0.;
	for (k = 1; k <= 40; ++k) {
	    w1 += 1. / k;
	    z__3.r = cr.r * -.25, z__3.i = cr.i * -.25;
	    i__1 = k * (k + 1);
	    d__1 = (doublereal) i__1;
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
	    z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * z2.i + 
		    z__2.i * z2.r;
	    cr.r = z__1.r, cr.i = z__1.i;
	    d__1 = w1 * 2. + 1. / (k + 1.);
	    z__1.r = d__1 * cr.r, z__1.i = d__1 * cr.i;
	    cp.r = z__1.r, cp.i = z__1.i;
	    z__1.r = cs.r + cp.r, z__1.i = cs.i + cp.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	    if (z_abs(&cp) < z_abs(&cs) * 1e-15) {
		goto L45;
	    }
/* L40: */
	}
L45:
	z__7.r = z1.r / 2., z__7.i = z1.i / 2.;
	z_log(&z__6, &z__7);
	z__5.r = z__6.r + el, z__5.i = z__6.i;
	z__4.r = z__5.r * cbj1.r - z__5.i * cbj1.i, z__4.i = z__5.r * cbj1.i 
		+ z__5.i * cbj1.r;
	z_div(&z__8, &c_b6, &z1);
	z__3.r = z__4.r - z__8.r, z__3.i = z__4.i - z__8.i;
	z__10.r = z1.r * .25, z__10.i = z1.i * .25;
	z__9.r = z__10.r * cs.r - z__10.i * cs.i, z__9.i = z__10.r * cs.i + 
		z__10.i * cs.r;
	z__2.r = z__3.r - z__9.r, z__2.i = z__3.i - z__9.i;
	z__1.r = rp2 * z__2.r, z__1.i = rp2 * z__2.i;
	cby1.r = z__1.r, cby1.i = z__1.i;
    } else {
	k0 = 12;
	if (a0 >= 35.f) {
	    k0 = 10;
	}
	if (a0 >= 50.f) {
	    k0 = 8;
	}
	d__1 = pi * .25;
	z__1.r = z1.r - d__1, z__1.i = z1.i;
	ct1.r = z__1.r, ct1.i = z__1.i;
	cp0.r = 1., cp0.i = 0.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
/* L50: */
	    i__2 = k - 1;
	    i__3 = k * -2;
	    pow_zi(&z__3, &z1, &i__3);
	    z__2.r = a[i__2] * z__3.r, z__2.i = a[i__2] * z__3.i;
	    z__1.r = cp0.r + z__2.r, z__1.i = cp0.i + z__2.i;
	    cp0.r = z__1.r, cp0.i = z__1.i;
	}
	z_div(&z__1, &c_b772, &z1);
	cq0.r = z__1.r, cq0.i = z__1.i;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
/* L55: */
	    i__3 = k - 1;
	    i__1 = k * -2 - 1;
	    pow_zi(&z__3, &z1, &i__1);
	    z__2.r = b[i__3] * z__3.r, z__2.i = b[i__3] * z__3.i;
	    z__1.r = cq0.r + z__2.r, z__1.i = cq0.i + z__2.i;
	    cq0.r = z__1.r, cq0.i = z__1.i;
	}
	z__3.r = rp2, z__3.i = 0.;
	z_div(&z__2, &z__3, &z1);
	z_sqrt(&z__1, &z__2);
	cu.r = z__1.r, cu.i = z__1.i;
	z_cos(&z__4, &ct1);
	z__3.r = cp0.r * z__4.r - cp0.i * z__4.i, z__3.i = cp0.r * z__4.i + 
		cp0.i * z__4.r;
	z_sin(&z__6, &ct1);
	z__5.r = cq0.r * z__6.r - cq0.i * z__6.i, z__5.i = cq0.r * z__6.i + 
		cq0.i * z__6.r;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cbj0.r = z__1.r, cbj0.i = z__1.i;
	z_sin(&z__4, &ct1);
	z__3.r = cp0.r * z__4.r - cp0.i * z__4.i, z__3.i = cp0.r * z__4.i + 
		cp0.i * z__4.r;
	z_cos(&z__6, &ct1);
	z__5.r = cq0.r * z__6.r - cq0.i * z__6.i, z__5.i = cq0.r * z__6.i + 
		cq0.i * z__6.r;
	z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cby0.r = z__1.r, cby0.i = z__1.i;
	d__1 = pi * .75;
	z__1.r = z1.r - d__1, z__1.i = z1.i;
	ct2.r = z__1.r, ct2.i = z__1.i;
	cp1.r = 1., cp1.i = 0.;
	i__3 = k0;
	for (k = 1; k <= i__3; ++k) {
/* L60: */
	    i__1 = k - 1;
	    i__2 = k * -2;
	    pow_zi(&z__3, &z1, &i__2);
	    z__2.r = a1[i__1] * z__3.r, z__2.i = a1[i__1] * z__3.i;
	    z__1.r = cp1.r + z__2.r, z__1.i = cp1.i + z__2.i;
	    cp1.r = z__1.r, cp1.i = z__1.i;
	}
	z_div(&z__1, &c_b776, &z1);
	cq1.r = z__1.r, cq1.i = z__1.i;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
/* L65: */
	    i__2 = k - 1;
	    i__3 = k * -2 - 1;
	    pow_zi(&z__3, &z1, &i__3);
	    z__2.r = b1[i__2] * z__3.r, z__2.i = b1[i__2] * z__3.i;
	    z__1.r = cq1.r + z__2.r, z__1.i = cq1.i + z__2.i;
	    cq1.r = z__1.r, cq1.i = z__1.i;
	}
	z_cos(&z__4, &ct2);
	z__3.r = cp1.r * z__4.r - cp1.i * z__4.i, z__3.i = cp1.r * z__4.i + 
		cp1.i * z__4.r;
	z_sin(&z__6, &ct2);
	z__5.r = cq1.r * z__6.r - cq1.i * z__6.i, z__5.i = cq1.r * z__6.i + 
		cq1.i * z__6.r;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cbj1.r = z__1.r, cbj1.i = z__1.i;
	z_sin(&z__4, &ct2);
	z__3.r = cp1.r * z__4.r - cp1.i * z__4.i, z__3.i = cp1.r * z__4.i + 
		cp1.i * z__4.r;
	z_cos(&z__6, &ct2);
	z__5.r = cq1.r * z__6.r - cq1.i * z__6.i, z__5.i = cq1.r * z__6.i + 
		cq1.i * z__6.r;
	z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
	z__1.r = cu.r * z__2.r - cu.i * z__2.i, z__1.i = cu.r * z__2.i + cu.i 
		* z__2.r;
	cby1.r = z__1.r, cby1.i = z__1.i;
    }
    if (z__->r < 0.f) {
	if (d_imag(z__) < 0.f) {
	    z__3.r = ci.r * 2., z__3.i = ci.i * 2.;
	    z__2.r = z__3.r * cbj0.r - z__3.i * cbj0.i, z__2.i = z__3.r * 
		    cbj0.i + z__3.i * cbj0.r;
	    z__1.r = cby0.r - z__2.r, z__1.i = cby0.i - z__2.i;
	    cby0.r = z__1.r, cby0.i = z__1.i;
	}
	if (d_imag(z__) > 0.f) {
	    z__3.r = ci.r * 2., z__3.i = ci.i * 2.;
	    z__2.r = z__3.r * cbj0.r - z__3.i * cbj0.i, z__2.i = z__3.r * 
		    cbj0.i + z__3.i * cbj0.r;
	    z__1.r = cby0.r + z__2.r, z__1.i = cby0.i + z__2.i;
	    cby0.r = z__1.r, cby0.i = z__1.i;
	}
	if (d_imag(z__) < 0.f) {
	    z__4.r = ci.r * 2., z__4.i = ci.i * 2.;
	    z__3.r = z__4.r * cbj1.r - z__4.i * cbj1.i, z__3.i = z__4.r * 
		    cbj1.i + z__4.i * cbj1.r;
	    z__2.r = cby1.r - z__3.r, z__2.i = cby1.i - z__3.i;
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    cby1.r = z__1.r, cby1.i = z__1.i;
	}
	if (d_imag(z__) > 0.f) {
	    z__4.r = ci.r * 2., z__4.i = ci.i * 2.;
	    z__3.r = z__4.r * cbj1.r - z__4.i * cbj1.i, z__3.i = z__4.r * 
		    cbj1.i + z__4.i * cbj1.r;
	    z__2.r = cby1.r + z__3.r, z__2.i = cby1.i + z__3.i;
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    cby1.r = z__1.r, cby1.i = z__1.i;
	}
	z__1.r = -cbj1.r, z__1.i = -cbj1.i;
	cbj1.r = z__1.r, cbj1.i = z__1.i;
    }
    z__1.r = -cby1.r, z__1.i = -cby1.i;
    cdy0.r = z__1.r, cdy0.i = z__1.i;
    z_div(&z__3, &c_b6, z__);
    z__2.r = z__3.r * cby1.r - z__3.i * cby1.i, z__2.i = z__3.r * cby1.i + 
	    z__3.i * cby1.r;
    z__1.r = cby0.r - z__2.r, z__1.i = cby0.i - z__2.i;
    cdy1.r = z__1.r, cdy1.i = z__1.i;
L70:
    if (*kf == 0) {
	zf->r = cby0.r, zf->i = cby0.i;
	zd->r = cdy0.r, zd->i = cdy0.i;
    } else if (*kf == 1) {
	zf->r = cby1.r, zf->i = cby1.i;
	zd->r = cdy1.r, zd->i = cdy1.i;
    } else if (*kf == 2) {
	zf->r = cdy1.r, zf->i = cdy1.i;
	z__3.r = -cdy1.r, z__3.i = -cdy1.i;
	z_div(&z__2, &z__3, z__);
	z__7.r = z__->r * z__->r - z__->i * z__->i, z__7.i = z__->r * z__->i 
		+ z__->i * z__->r;
	z_div(&z__6, &c_b6, &z__7);
	z__5.r = 1. - z__6.r, z__5.i = -z__6.i;
	z__4.r = z__5.r * cby1.r - z__5.i * cby1.i, z__4.i = z__5.r * cby1.i 
		+ z__5.i * cby1.r;
	z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
	zd->r = z__1.r, zd->i = z__1.i;
    }
    return 0;
} /* cy01_ */

/*       ********************************** */
/* Subroutine */ int ffk_(integer *ks, doublereal *x, doublereal *fr, 
	doublereal *fi, doublereal *fm, doublereal *fa, doublereal *gr, 
	doublereal *gi, doublereal *gm, doublereal *ga)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer pow_ii(integer *, integer *);
    double sin(doublereal), cos(doublereal), atan(doublereal);

    /* Local variables */
    static integer k, m;
    static doublereal c1, s1, x2, x4, cs, pi, xa, xc, xf, xg, ss, xp, xq, xr, 
	    xs, xw, fi0, xf0, xf1, pp2, p2p, xq2, eps, srd, xsu;


/*       ======================================================= */
/*       Purpose: Compute modified Fresnel integrals F±(x) */
/*                and K±(x) */
/*       Input :  x   --- Argument of F±(x) and K±(x) */
/*                KS  --- Sign code */
/*                        KS=0 for calculating F+(x) and K+(x) */
/*                        KS=1 for calculating F_(x) and K_(x) */
/*       Output:  FR  --- Re[F±(x)] */
/*                FI  --- Im[F±(x)] */
/*                FM  --- |F±(x)| */
/*                FA  --- Arg[F±(x)]  (Degs.) */
/*                GR  --- Re[K±(x)] */
/*                GI  --- Im[K±(x)] */
/*                GM  --- |K±(x)| */
/*                GA  --- Arg[K±(x)]  (Degs.) */
/*       ====================================================== */

    srd = 57.29577951308233;
    eps = 1e-15;
    pi = 3.141592653589793;
    pp2 = 1.2533141373155;
    p2p = .7978845608028654;
    xa = abs(*x);
    x2 = *x * *x;
    x4 = x2 * x2;
    if (*x == 0.) {
	*fr = sqrt(pi * .5) * .5;
	*fi = pow_ii(&c_n1, ks) * *fr;
	*fm = sqrt(pi * .25);
	*fa = pow_ii(&c_n1, ks) * 45.;
	*gr = .5;
	*gi = 0.;
	*gm = .5;
	*ga = 0.;
    } else {
	if (xa <= 2.5) {
	    xr = p2p * xa;
	    c1 = xr;
	    for (k = 1; k <= 50; ++k) {
		xr = xr * -.5 * (k * 4. - 3.) / k / (k * 2. - 1.) / (k * 4. + 
			1.) * x4;
		c1 += xr;
		if ((d__1 = xr / c1, abs(d__1)) < eps) {
		    goto L15;
		}
/* L10: */
	    }
L15:
	    s1 = p2p * xa * xa * xa / 3.;
	    xr = s1;
	    for (k = 1; k <= 50; ++k) {
		xr = xr * -.5 * (k * 4. - 1.) / k / (k * 2. + 1.) / (k * 4. + 
			3.) * x4;
		s1 += xr;
		if ((d__1 = xr / s1, abs(d__1)) < eps) {
		    goto L40;
		}
/* L20: */
	    }
	} else if (xa < 5.5) {
	    m = (integer) (x2 * 1.75f + 42);
	    xsu = 0.;
	    xc = 0.;
	    xs = 0.;
	    xf1 = 0.;
	    xf0 = 1e-100;
	    for (k = m; k >= 0; --k) {
		xf = (k * 2. + 3.) * xf0 / x2 - xf1;
		if (k == k / 2 << 1) {
		    xc += xf;
		} else {
		    xs += xf;
		}
		xsu += (k * 2. + 1.) * xf * xf;
		xf1 = xf0;
/* L25: */
		xf0 = xf;
	    }
	    xq = sqrt(xsu);
	    xw = p2p * xa / xq;
	    c1 = xc * xw;
	    s1 = xs * xw;
	} else {
	    xr = 1.;
	    xf = 1.;
	    for (k = 1; k <= 12; ++k) {
		xr = xr * -.25 * (k * 4. - 1.) * (k * 4. - 3.) / x4;
/* L30: */
		xf += xr;
	    }
	    xr = 1. / (xa * 2. * xa);
	    xg = xr;
	    for (k = 1; k <= 12; ++k) {
		xr = xr * -.25 * (k * 4. + 1.) * (k * 4. - 1.) / x4;
/* L35: */
		xg += xr;
	    }
	    c1 = (xf * sin(x2) - xg * cos(x2)) / sqrt(pi * 2.) / xa + .5;
	    s1 = .5 - (xf * cos(x2) + xg * sin(x2)) / sqrt(pi * 2.) / xa;
	}
L40:
	*fr = pp2 * (.5 - c1);
	fi0 = pp2 * (.5 - s1);
	*fi = pow_ii(&c_n1, ks) * fi0;
	*fm = sqrt(*fr * *fr + *fi * *fi);
	if (*fr >= 0.f) {
	    *fa = srd * atan(*fi / *fr);
	} else if (*fi > 0.f) {
	    *fa = srd * (atan(*fi / *fr) + pi);
	} else if (*fi < 0.f) {
	    *fa = srd * (atan(*fi / *fr) - pi);
	}
	xp = *x * *x + pi / 4.;
	cs = cos(xp);
	ss = sin(xp);
	xq2 = 1. / sqrt(pi);
	*gr = xq2 * (*fr * cs + fi0 * ss);
	*gi = pow_ii(&c_n1, ks) * xq2 * (fi0 * cs - *fr * ss);
	*gm = sqrt(*gr * *gr + *gi * *gi);
	if (*gr >= 0.f) {
	    *ga = srd * atan(*gi / *gr);
	} else if (*gi > 0.f) {
	    *ga = srd * (atan(*gi / *gr) + pi);
	} else if (*gi < 0.f) {
	    *ga = srd * (atan(*gi / *gr) - pi);
	}
	if (*x < 0.) {
	    *fr = pp2 - *fr;
	    *fi = pow_ii(&c_n1, ks) * pp2 - *fi;
	    *fm = sqrt(*fr * *fr + *fi * *fi);
	    *fa = srd * atan(*fi / *fr);
	    *gr = cos(*x * *x) - *gr;
	    *gi = -pow_ii(&c_n1, ks) * sin(*x * *x) - *gi;
	    *gm = sqrt(*gr * *gr + *gi * *gi);
	    *ga = srd * atan(*gi / *gr);
	}
    }
    return 0;
} /* ffk_ */

/*       ********************************** */
/* Subroutine */ int airya_(doublereal *x, doublereal *ai, doublereal *bi, 
	doublereal *ad, doublereal *bd)
{
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal z__, c1, c2, xa, xq, vi1, vj1, vj2, vi2, vk1, vk2, sr3, 
	    vy1, vy2, pir;
    extern /* Subroutine */ int ajyik_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute Airy functions and their derivatives */
/*       Input:   x  --- Argument of Airy function */
/*       Output:  AI --- Ai(x) */
/*                BI --- Bi(x) */
/*                AD --- Ai'(x) */
/*                BD --- Bi'(x) */
/*       Routine called: */
/*                AJYIK for computing Jv(x), Yv(x), Iv(x) and */
/*                Kv(x) with v=1/3 and 2/3 */
/*       ====================================================== */

    xa = abs(*x);
    pir = .318309886183891;
    c1 = .355028053887817;
    c2 = .258819403792807;
    sr3 = 1.732050807568877;
    z__ = pow_dd(&xa, &c_b51) / 1.5;
    xq = sqrt(xa);
    ajyik_(&z__, &vj1, &vj2, &vy1, &vy2, &vi1, &vi2, &vk1, &vk2);
    if (*x == 0.) {
	*ai = c1;
	*bi = sr3 * c1;
	*ad = -c2;
	*bd = sr3 * c2;
    } else if (*x > 0.) {
	*ai = pir * xq / sr3 * vk1;
	*bi = xq * (pir * vk1 + 2. / sr3 * vi1);
	*ad = -xa / sr3 * pir * vk2;
	*bd = xa * (pir * vk2 + 2. / sr3 * vi2);
    } else {
	*ai = xq * .5 * (vj1 - vy1 / sr3);
	*bi = xq * -.5 * (vj1 / sr3 + vy1);
	*ad = xa * .5 * (vj2 + vy2 / sr3);
	*bd = xa * .5 * (vj2 / sr3 - vy2);
    }
    return 0;
} /* airya_ */

/*       ********************************** */
/* Subroutine */ int airyb_(doublereal *x, doublereal *ai, doublereal *bi, 
	doublereal *ad, doublereal *bd)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer k;
    static doublereal r__, c1, c2, df, dg, ck[51], dk[51];
    static integer km;
    static doublereal pi, xa, xe, fx, gx, xf, rp, xm, xq;
    static integer km2;
    static doublereal sr3, xp1, xr1, xr2, sad, sbd, sda, sdb, sai, sbi, ssa, 
	    eps, ssb, xar, xcs, xss;
    static integer kmax;


/*       ======================================================= */
/*       Purpose: Compute Airy functions and their derivatives */
/*       Input:   x  --- Argument of Airy function */
/*       Output:  AI --- Ai(x) */
/*                BI --- Bi(x) */
/*                AD --- Ai'(x) */
/*                BD --- Bi'(x) */
/*       ======================================================= */

    eps = 1e-15;
    pi = 3.141592653589793;
    c1 = .355028053887817;
    c2 = .258819403792807;
    sr3 = 1.732050807568877;
    xa = abs(*x);
    xq = sqrt(xa);
    xm = 8.;
    if (*x > 0.) {
	xm = 5.;
    }
    if (*x == 0.) {
	*ai = c1;
	*bi = sr3 * c1;
	*ad = -c2;
	*bd = sr3 * c2;
	return 0;
    }
    if (xa <= xm) {
	fx = 1.;
	r__ = 1.;
	for (k = 1; k <= 40; ++k) {
	    r__ = r__ * *x / (k * 3.) * *x / (k * 3. - 1.) * *x;
	    fx += r__;
	    if (abs(r__) < abs(fx) * eps) {
		goto L15;
	    }
/* L10: */
	}
L15:
	gx = *x;
	r__ = *x;
	for (k = 1; k <= 40; ++k) {
	    r__ = r__ * *x / (k * 3.) * *x / (k * 3. + 1.) * *x;
	    gx += r__;
	    if (abs(r__) < abs(gx) * eps) {
		goto L25;
	    }
/* L20: */
	}
L25:
	*ai = c1 * fx - c2 * gx;
	*bi = sr3 * (c1 * fx + c2 * gx);
	df = *x * .5 * *x;
	r__ = df;
	for (k = 1; k <= 40; ++k) {
	    r__ = r__ * *x / (k * 3.) * *x / (k * 3. + 2.) * *x;
	    df += r__;
	    if (abs(r__) < abs(df) * eps) {
		goto L35;
	    }
/* L30: */
	}
L35:
	dg = 1.;
	r__ = 1.;
	for (k = 1; k <= 40; ++k) {
	    r__ = r__ * *x / (k * 3.) * *x / (k * 3. - 2.) * *x;
	    dg += r__;
	    if (abs(r__) < abs(dg) * eps) {
		goto L45;
	    }
/* L40: */
	}
L45:
	*ad = c1 * df - c2 * dg;
	*bd = sr3 * (c1 * df + c2 * dg);
    } else {
	km = (integer) (24.5f - xa);
	if (xa < 6.f) {
	    km = 14;
	}
	if (xa > 15.f) {
	    km = 10;
	}
	if (*x > 0.) {
	    kmax = km;
	} else {
/*             Choose cutoffs so that the remainder term in asymptotic */
/*             expansion is epsilon size. The X<0 branch needs to be fast */
/*             in order to make AIRYZO efficient */
	    if (xa > 70.f) {
		km = 3;
	    }
	    if (xa > 500.f) {
		km = 2;
	    }
	    if (xa > 1e3f) {
		km = 1;
	    }
	    km2 = km;
	    if (xa > 150.f) {
		km2 = 1;
	    }
	    if (xa > 3e3f) {
		km2 = 0;
	    }
	    kmax = (km << 1) + 1;
	}
	xe = xa * xq / 1.5;
	xr1 = 1. / xe;
	xar = 1. / xq;
	xf = sqrt(xar);
	rp = .5641895835477563;
	r__ = 1.;
	i__1 = kmax;
	for (k = 1; k <= i__1; ++k) {
	    r__ = r__ * (k * 6. - 1.) / 216. * (k * 6. - 3.) / k * (k * 6. - 
		    5.) / (k * 2. - 1.);
	    ck[k - 1] = r__;
/* L50: */
	    dk[k - 1] = -(k * 6. + 1.) / (k * 6. - 1.) * ck[k - 1];
	}
	if (*x > 0.) {
	    sai = 1.;
	    sad = 1.;
	    r__ = 1.;
	    i__1 = km;
	    for (k = 1; k <= i__1; ++k) {
		r__ = -r__ * xr1;
		sai += ck[k - 1] * r__;
/* L55: */
		sad += dk[k - 1] * r__;
	    }
	    sbi = 1.;
	    sbd = 1.;
	    r__ = 1.;
	    i__1 = km;
	    for (k = 1; k <= i__1; ++k) {
		r__ *= xr1;
		sbi += ck[k - 1] * r__;
/* L60: */
		sbd += dk[k - 1] * r__;
	    }
	    xp1 = exp(-xe);
	    *ai = rp * .5 * xf * xp1 * sai;
	    *bi = rp * xf / xp1 * sbi;
	    *ad = rp * -.5 / xf * xp1 * sad;
	    *bd = rp / xf / xp1 * sbd;
	} else {
	    xcs = cos(xe + pi / 4.);
	    xss = sin(xe + pi / 4.);
	    ssa = 1.;
	    sda = 1.;
	    r__ = 1.;
	    xr2 = 1. / (xe * xe);
	    i__1 = km;
	    for (k = 1; k <= i__1; ++k) {
		r__ = -r__ * xr2;
		ssa += ck[(k << 1) - 1] * r__;
/* L65: */
		sda += dk[(k << 1) - 1] * r__;
	    }
	    ssb = ck[0] * xr1;
	    sdb = dk[0] * xr1;
	    r__ = xr1;
	    i__1 = km2;
	    for (k = 1; k <= i__1; ++k) {
		r__ = -r__ * xr2;
		ssb += ck[k * 2] * r__;
/* L70: */
		sdb += dk[k * 2] * r__;
	    }
	    *ai = rp * xf * (xss * ssa - xcs * ssb);
	    *bi = rp * xf * (xcs * ssa + xss * ssb);
	    *ad = -rp / xf * (xcs * sda + xss * sdb);
	    *bd = rp / xf * (xss * sda - xcs * sdb);
	}
    }
    return 0;
} /* airyb_ */

/*       ********************************** */
/* Subroutine */ int scka_(integer *m, integer *n, doublereal *c__, 
	doublereal *cv, integer *kd, doublereal *ck)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal f;
    static integer j, k;
    static doublereal f0, f1, f2;
    static integer k1;
    static doublereal r1, r2, s0;
    static integer kb;
    static doublereal fl, cs;
    static integer ip, nm;
    static doublereal fs, su1, su2;


/*       ====================================================== */
/*       Purpose: Compute the expansion coefficients of the */
/*                prolate and oblate spheroidal functions, c2k */
/*       Input :  m  --- Mode parameter */
/*                n  --- Mode parameter */
/*                c  --- Spheroidal parameter */
/*                cv --- Characteristic value */
/*                KD --- Function code */
/*                       KD=1 for prolate; KD=-1 for oblate */
/*       Output:  CK(k) --- Expansion coefficients ck; */
/*                          CK(1), CK(2),... correspond to */
/*                          c0, c2,... */
/*       ====================================================== */

    /* Parameter adjustments */
    --ck;

    /* Function Body */
    if (*c__ <= 1e-10) {
	*c__ = 1e-10;
    }
    nm = (integer) ((*n - *m) / 2 + *c__) + 25;
    cs = *c__ * *c__ * *kd;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    fs = 1.;
    f1 = 0.;
    f0 = 1e-100;
    kb = 0;
    ck[nm + 1] = 0.;
    fl = 0.;
    for (k = nm; k >= 1; --k) {
	f = (((k * 2. + *m + ip) * (k * 2. + *m + 1. + ip) - *cv + cs) * f0 - 
		(k + 1.) * 4. * (k + *m + 1.) * f1) / cs;
	if (abs(f) > (d__1 = ck[k + 1], abs(d__1))) {
	    ck[k] = f;
	    f1 = f0;
	    f0 = f;
	    if (abs(f) > 1e100) {
		i__1 = k;
		for (k1 = nm; k1 >= i__1; --k1) {
/* L5: */
		    ck[k1] *= 1e-100;
		}
		f1 *= 1e-100;
		f0 *= 1e-100;
	    }
	} else {
	    kb = k;
	    fl = ck[k + 1];
	    f1 = 1.;
	    f2 = ((*m + ip) * (*m + ip + 1.f) - *cv + cs) * .25 / (*m + 1.f) *
		     f1;
	    ck[1] = f1;
	    if (kb == 1) {
		fs = f2;
	    } else if (kb == 2) {
		ck[2] = f2;
		fs = (((*m + ip + 2.f) * (*m + ip + 3.f) - *cv + cs) * f2 - 
			cs * f1) * .125 / (*m + 2.f);
	    } else {
		ck[2] = f2;
		i__1 = kb + 1;
		for (j = 3; j <= i__1; ++j) {
		    f = (((j * 2.f + *m + ip - 4.f) * (j * 2.f + *m + ip - 
			    3.f) - *cv + cs) * f2 - cs * f1) * .25 / ((j - 
			    1.f) * (j + *m - 1.f));
		    if (j <= kb) {
			ck[j] = f;
		    }
		    f1 = f2;
/* L10: */
		    f2 = f;
		}
		fs = f;
	    }
	    goto L20;
	}
/* L15: */
    }
L20:
    su1 = 0.;
    i__1 = kb;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	su1 += ck[k];
    }
    su2 = 0.;
    i__1 = nm;
    for (k = kb + 1; k <= i__1; ++k) {
/* L30: */
	su2 += ck[k];
    }
    r1 = 1.;
    i__1 = (*n + *m + ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L35: */
	r1 *= j + (*n + *m + ip) * .5;
    }
    r2 = 1.;
    i__1 = (*n - *m - ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L40: */
	r2 = -r2 * j;
    }
    if (kb == 0) {
	s0 = r1 / (pow_di(&c_b4, n) * r2 * su2);
    } else {
	s0 = r1 / (pow_di(&c_b4, n) * r2 * (fl / fs * su1 + su2));
    }
    i__1 = kb;
    for (k = 1; k <= i__1; ++k) {
/* L45: */
	ck[k] = fl / fs * s0 * ck[k];
    }
    i__1 = nm;
    for (k = kb + 1; k <= i__1; ++k) {
/* L50: */
	ck[k] = s0 * ck[k];
    }
    return 0;
} /* scka_ */

/*       ********************************** */
/* Subroutine */ int sckb_(integer *m, integer *n, doublereal *c__, 
	doublereal *df, doublereal *ck)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, k;
    static doublereal r__, d1, d2, d3;
    static integer i1, i2;
    static doublereal r1;
    static integer ip, nm;
    static doublereal sw, fac, reg, sum;


/*       ====================================================== */
/*       Purpose: Compute the expansion coefficients of the */
/*                prolate and oblate spheroidal functions */
/*       Input :  m  --- Mode parameter */
/*                n  --- Mode parameter */
/*                c  --- Spheroidal parameter */
/*                DF(k) --- Expansion coefficients dk */
/*       Output:  CK(k) --- Expansion coefficients ck; */
/*                          CK(1), CK(2), ... correspond to */
/*                          c0, c2, ... */
/*       ====================================================== */

    /* Parameter adjustments */
    --ck;
    --df;

    /* Function Body */
    if (*c__ <= 1e-10) {
	*c__ = 1e-10;
    }
    nm = (integer) ((*n - *m) * .5f + *c__) + 25;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    reg = 1.;
    if (*m + nm > 80) {
	reg = 1e-200;
    }
    fac = -pow_di(&c_b50, m);
    sw = 0.;
    i__1 = nm - 1;
    for (k = 0; k <= i__1; ++k) {
	fac = -fac;
	i1 = (k << 1) + ip + 1;
	r__ = reg;
	i__2 = i1 + (*m << 1) - 1;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L10: */
	    r__ *= i__;
	}
	i2 = k + *m + ip;
	i__2 = i2 + k - 1;
	for (i__ = i2; i__ <= i__2; ++i__) {
/* L15: */
	    r__ *= i__ + .5;
	}
	sum = r__ * df[k + 1];
	i__2 = nm;
	for (i__ = k + 1; i__ <= i__2; ++i__) {
	    d1 = i__ * 2. + ip;
	    d2 = *m * 2. + d1;
	    d3 = i__ + *m + ip - .5;
	    r__ = r__ * d2 * (d2 - 1.) * i__ * (d3 + k) / (d1 * (d1 - 1.) * (
		    i__ - k) * d3);
	    sum += r__ * df[i__ + 1];
	    if ((d__1 = sw - sum, abs(d__1)) < abs(sum) * 1e-14) {
		goto L25;
	    }
/* L20: */
	    sw = sum;
	}
L25:
	r1 = reg;
	i__2 = *m + k;
	for (i__ = 2; i__ <= i__2; ++i__) {
/* L30: */
	    r1 *= i__;
	}
/* L35: */
	ck[k + 1] = fac * sum / r1;
    }
    return 0;
} /* sckb_ */

/*       ********************************** */
/* Subroutine */ int cpdla_(integer *n, doublecomplex *z__, doublecomplex *
	cdn)
{
    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    void z_exp(doublecomplex *, doublecomplex *), pow_zi(doublecomplex *, 
	    doublecomplex *, integer *), z_div(doublecomplex *, doublecomplex 
	    *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex cr, cb0;


/*       =========================================================== */
/*       Purpose: Compute complex parabolic cylinder function Dn(z) */
/*                for large argument */
/*       Input:   z   --- Complex argument of Dn(z) */
/*                n   --- Order of Dn(z) (n = 0,±1,±2,…) */
/*       Output:  CDN --- Dn(z) */
/*       =========================================================== */

    pow_zi(&z__2, z__, n);
    z__5.r = z__->r * -.25, z__5.i = z__->i * -.25;
    z__4.r = z__5.r * z__->r - z__5.i * z__->i, z__4.i = z__5.r * z__->i + 
	    z__5.i * z__->r;
    z_exp(&z__3, &z__4);
    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * z__3.i + 
	    z__2.i * z__3.r;
    cb0.r = z__1.r, cb0.i = z__1.i;
    cr.r = 1., cr.i = 0.;
    cdn->r = 1., cdn->i = 0.;
    for (k = 1; k <= 16; ++k) {
	z__4.r = cr.r * -.5, z__4.i = cr.i * -.5;
	r__1 = k * 2.f - *n - 1.f;
	z__3.r = r__1 * z__4.r, z__3.i = r__1 * z__4.i;
	r__2 = k * 2.f - *n - 2.f;
	z__2.r = r__2 * z__3.r, z__2.i = r__2 * z__3.i;
	d__1 = (doublereal) k;
	z__6.r = d__1 * z__->r, z__6.i = d__1 * z__->i;
	z__5.r = z__6.r * z__->r - z__6.i * z__->i, z__5.i = z__6.r * z__->i 
		+ z__6.i * z__->r;
	z_div(&z__1, &z__2, &z__5);
	cr.r = z__1.r, cr.i = z__1.i;
	z__1.r = cdn->r + cr.r, z__1.i = cdn->i + cr.i;
	cdn->r = z__1.r, cdn->i = z__1.i;
	if (z_abs(&cr) < z_abs(cdn) * 1e-12) {
	    goto L15;
	}
/* L10: */
    }
L15:
    z__1.r = cb0.r * cdn->r - cb0.i * cdn->i, z__1.i = cb0.r * cdn->i + cb0.i 
	    * cdn->r;
    cdn->r = z__1.r, cdn->i = z__1.i;
    return 0;
} /* cpdla_ */

/*       ********************************** */
/* Subroutine */ int fcszo_(integer *kf, integer *nt, doublecomplex *zo)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), log(
	    doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal w;
    static doublecomplex z__;
    static doublereal w0, pi;
    static integer it;
    static doublecomplex zd;
    static integer nr;
    static doublecomplex zf;
    static doublereal px, py;
    static doublecomplex zp, zq, zw;
    extern /* Subroutine */ int cfc_(doublecomplex *, doublecomplex *, 
	    doublecomplex *), cfs_(doublecomplex *, doublecomplex *, 
	    doublecomplex *);
    static doublecomplex zfd, zgd;
    static doublereal psq;


/*       =============================================================== */
/*       Purpose: Compute the complex zeros of Fresnel integral C(z) */
/*                or S(z) using modified Newton's iteration method */
/*       Input :  KF  --- Function code */
/*                        KF=1 for C(z) or KF=2 for S(z) */
/*                NT  --- Total number of zeros */
/*       Output:  ZO(L) --- L-th zero of C(z) or S(z) */
/*       Routines called: */
/*            (1) CFC for computing Fresnel integral C(z) */
/*            (2) CFS for computing Fresnel integral S(z) */
/*       ============================================================== */

    /* Parameter adjustments */
    --zo;

    /* Function Body */
    pi = 3.141592653589793;
    psq = 0.;
    w = 0.;
    i__1 = *nt;
    for (nr = 1; nr <= i__1; ++nr) {
	if (*kf == 1) {
	    psq = sqrt(nr * 4. - 1.);
	}
	if (*kf == 2) {
	    d__1 = (doublereal) nr;
	    psq = pow_dd(&d__1, &c_b50) * 2.;
	}
	px = psq - log(pi * psq) / (pi * pi * pow_dd(&psq, &c_b149));
	py = log(pi * psq) / (pi * psq);
	z__1.r = px, z__1.i = py;
	z__.r = z__1.r, z__.i = z__1.i;
	if (*kf == 2) {
	    if (nr == 2) {
		z__.r = 2.8334f, z__.i = .2443f;
	    }
	    if (nr == 3) {
		z__.r = 3.4674f, z__.i = .2185f;
	    }
	    if (nr == 4) {
		z__.r = 4.0025f, z__.i = .2008f;
	    }
	}
	it = 0;
L15:
	++it;
	if (*kf == 1) {
	    cfc_(&z__, &zf, &zd);
	}
	if (*kf == 2) {
	    cfs_(&z__, &zf, &zd);
	}
	zp.r = 1., zp.i = 0.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    i__3 = i__;
	    z__2.r = z__.r - zo[i__3].r, z__2.i = z__.i - zo[i__3].i;
	    z__1.r = zp.r * z__2.r - zp.i * z__2.i, z__1.i = zp.r * z__2.i + 
		    zp.i * z__2.r;
	    zp.r = z__1.r, zp.i = z__1.i;
	}
	z_div(&z__1, &zf, &zp);
	zfd.r = z__1.r, zfd.i = z__1.i;
	zq.r = 0., zq.i = 0.;
	i__3 = nr - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    zw.r = 1., zw.i = 0.;
	    i__2 = nr - 1;
	    for (j = 1; j <= i__2; ++j) {
		if (j == i__) {
		    goto L25;
		}
		i__4 = j;
		z__2.r = z__.r - zo[i__4].r, z__2.i = z__.i - zo[i__4].i;
		z__1.r = zw.r * z__2.r - zw.i * z__2.i, z__1.i = zw.r * 
			z__2.i + zw.i * z__2.r;
		zw.r = z__1.r, zw.i = z__1.i;
L25:
		;
	    }
/* L30: */
	    z__1.r = zq.r + zw.r, z__1.i = zq.i + zw.i;
	    zq.r = z__1.r, zq.i = z__1.i;
	}
	z__3.r = zq.r * zfd.r - zq.i * zfd.i, z__3.i = zq.r * zfd.i + zq.i * 
		zfd.r;
	z__2.r = zd.r - z__3.r, z__2.i = zd.i - z__3.i;
	z_div(&z__1, &z__2, &zp);
	zgd.r = z__1.r, zgd.i = z__1.i;
	z_div(&z__2, &zfd, &zgd);
	z__1.r = z__.r - z__2.r, z__1.i = z__.i - z__2.i;
	z__.r = z__1.r, z__.i = z__1.i;
	w0 = w;
	w = z_abs(&z__);
	if (it <= 50 && (d__1 = (w - w0) / w, abs(d__1)) > 1e-12) {
	    goto L15;
	}
/* L35: */
	i__3 = nr;
	zo[i__3].r = z__.r, zo[i__3].i = z__.i;
    }
    return 0;
} /* fcszo_ */

/*       ********************************** */
/* Subroutine */ int e1xa_(doublereal *x, doublereal *e1)
{
    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal es1, es2;


/*       ============================================ */
/*       Purpose: Compute exponential integral E1(x) */
/*       Input :  x  --- Argument of E1(x) */
/*       Output:  E1 --- E1(x) ( x > 0 ) */
/*       ============================================ */

    if (*x == 0.f) {
	*e1 = 1e300;
    } else if (*x <= 1.f) {
	*e1 = -log(*x) + ((((*x * .00107857 - .00976004) * *x + .05519968) * *
		x - .24991055) * *x + .99999193) * *x - .57721566;
    } else {
	es1 = (((*x + 8.5733287401) * *x + 18.059016973) * *x + 8.6347608925) 
		* *x + .2677737343;
	es2 = (((*x + 9.5733223454) * *x + 25.6329561486) * *x + 
		21.0996530827) * *x + 3.9584969228;
	*e1 = exp(-(*x)) / *x * es1 / es2;
    }
    return 0;
} /* e1xa_ */

/*       ********************************** */
/* Subroutine */ int lpmv0_(doublereal *v, integer *m, doublereal *x, 
	doublereal *pmv)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer pow_ii(integer *, integer *);
    double sin(doublereal), tan(doublereal), log(doublereal);

    /* Local variables */
    extern /* Subroutine */ int psi_spec__(doublereal *, doublereal *);
    static integer j, k;
    static doublereal r__, s, c0, r0, s0, r2, r1, v0, s1, s2, el, pa, pi, rg, 
	    qr;
    static integer nv;
    static doublereal xq, vs, pv0, eps, pss, psv;


/*       ======================================================= */
/*       Purpose: Compute the associated Legendre function */
/*                Pmv(x) with an integer order and an arbitrary */
/*                nonnegative degree v */
/*       Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 ) */
/*                m   --- Order of Pmv(x) */
/*                v   --- Degree of Pmv(x) */
/*       Output:  PMV --- Pmv(x) */
/*       Routine called:  PSI_SPEC for computing Psi function */
/*       ======================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    eps = 1e-14;
    nv = (integer) (*v);
    v0 = *v - nv;
    if (*x == -1. && *v != (doublereal) nv) {
	if (*m == 0) {
	    *pmv = -1e300;
	}
	if (*m != 0) {
	    *pmv = 1e300;
	}
	return 0;
    }
    c0 = 1.;
    if (*m != 0) {
	rg = *v * (*v + *m);
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
/* L10: */
	    rg *= *v * *v - j * j;
	}
	xq = sqrt(1. - *x * *x);
	r0 = 1.;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
/* L15: */
	    r0 = r0 * .5 * xq / j;
	}
	c0 = r0 * rg;
    }
    if (v0 == 0.) {
/*          DLMF 14.3.4, 14.7.17, 15.2.4 */
	*pmv = 1.;
	r__ = 1.;
	i__1 = nv - *m;
	for (k = 1; k <= i__1; ++k) {
	    r__ = r__ * .5 * (-nv + *m + k - 1.) * (nv + *m + k) / (k * (k + *
		    m)) * (*x + 1.);
/* L20: */
	    *pmv += r__;
	}
	*pmv = pow_ii(&c_n1, &nv) * c0 * *pmv;
    } else {
	if (*x >= -.35) {
/*             DLMF 14.3.4, 15.2.1 */
	    *pmv = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 100; ++k) {
		r__ = r__ * .5 * (-(*v) + *m + k - 1.) * (*v + *m + k) / (k * 
			(*m + k)) * (1. - *x);
		*pmv += r__;
		if (k > 12 && (d__1 = r__ / *pmv, abs(d__1)) < eps) {
		    goto L30;
		}
/* L25: */
	    }
L30:
	    *pmv = pow_ii(&c_n1, m) * c0 * *pmv;
	} else {
/*             DLMF 14.3.5, 15.8.10 */
	    vs = sin(*v * pi) / pi;
	    pv0 = 0.;
	    if (*m != 0) {
		qr = sqrt((1. - *x) / (*x + 1.));
		r2 = 1.;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
/* L35: */
		    r2 = r2 * qr * j;
		}
		s0 = 1.;
		r1 = 1.;
		i__1 = *m - 1;
		for (k = 1; k <= i__1; ++k) {
		    r1 = r1 * .5 * (-(*v) + k - 1) * (*v + k) / (k * (k - *m))
			     * (*x + 1.);
/* L40: */
		    s0 += r1;
		}
		pv0 = -vs * r2 / *m * s0;
	    }
	    psi_spec__(v, &psv);
	    pa = (psv + el) * 2. + pi / tan(pi * *v) + 1. / *v;
	    s1 = 0.;
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
/* L45: */
		s1 += (j * j + *v * *v) / (j * (j * j - *v * *v));
	    }
	    *pmv = pa + s1 - 1. / (*m - *v) + log((*x + 1.) * .5);
	    r__ = 1.;
	    for (k = 1; k <= 100; ++k) {
		r__ = r__ * .5 * (-(*v) + *m + k - 1.) * (*v + *m + k) / (k * 
			(k + *m)) * (*x + 1.);
		s = 0.;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
/* L50: */
/* Computing 2nd power */
		    i__2 = k + j;
/* Computing 2nd power */
		    i__3 = k + j;
		    s += (i__2 * i__2 + *v * *v) / ((k + j) * (i__3 * i__3 - *
			    v * *v));
		}
		s2 = 0.;
		i__2 = k;
		for (j = 1; j <= i__2; ++j) {
/* L55: */
		    s2 += 1. / (j * (j * j - *v * *v));
		}
		pss = pa + s + *v * 2. * *v * s2 - 1. / (*m + k - *v) + log((*
			x + 1.) * .5);
		r2 = pss * r__;
		*pmv += r2;
		if ((d__1 = r2 / *pmv, abs(d__1)) < eps) {
		    goto L65;
		}
/* L60: */
	    }
L65:
	    *pmv = pv0 + *pmv * vs * c0;
	}
    }
    return 0;
} /* lpmv0_ */

/*       ********************************** */
/* Subroutine */ int lpmv_(doublereal *v, integer *m, doublereal *x, 
	doublereal *pmv)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer j;
    static doublereal g1, g2, p0, p1, v0;
    static integer nv, mx;
    static doublereal vx;
    extern doublereal dinf_(void), dnan_(void);
    extern /* Subroutine */ int lpmv0_(doublereal *, integer *, doublereal *, 
	    doublereal *);
    static integer neg_m__;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ======================================================= */
/*       Purpose: Compute the associated Legendre function */
/*                Pmv(x) with an integer order and an arbitrary */
/*                degree v, using recursion for large degrees */
/*       Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 ) */
/*                m   --- Order of Pmv(x) */
/*                v   --- Degree of Pmv(x) */
/*       Output:  PMV --- Pmv(x) */
/*       Routine called:  LPMV0 */
/*       ======================================================= */

    if (*x == -1. && *v != (doublereal) ((integer) (*v))) {
	if (*m == 0) {
	    *pmv = -dinf_();
	}
	if (*m != 0) {
	    *pmv = dinf_();
	}
	return 0;
    }
    vx = *v;
    mx = *m;
/*       DLMF 14.9.5 */
    if (*v < 0.) {
	vx = -vx - 1;
    }
    neg_m__ = 0;
    if (*m < 0) {
	if (vx + *m + 1 > 0. || vx != (doublereal) ((integer) vx)) {
	    neg_m__ = 1;
	    mx = -(*m);
	} else {
/*             We don't handle cases where DLMF 14.9.3 doesn't help */
	    *pmv = dnan_();
	    return 0;
	}
    }
    nv = (integer) vx;
    v0 = vx - nv;
    if (nv > 2 && nv > mx) {
/*          Up-recursion on degree, AMS 8.5.3 / DLMF 14.10.3 */
	d__1 = v0 + mx;
	lpmv0_(&d__1, &mx, x, &p0);
	d__1 = v0 + mx + 1;
	lpmv0_(&d__1, &mx, x, &p1);
	*pmv = p1;
	i__1 = nv;
	for (j = mx + 2; j <= i__1; ++j) {
	    *pmv = (((v0 + j) * 2 - 1) * *x * p1 - (v0 + j - 1 + mx) * p0) / (
		    v0 + j - mx);
	    p0 = p1;
	    p1 = *pmv;
/* L10: */
	}
    } else {
	lpmv0_(&vx, &mx, x, pmv);
    }
    if (neg_m__ != 0 && abs(*pmv) < 1e300) {
/*          DLMF 14.9.3 */
	d__1 = vx - mx + 1;
	gamma2_(&d__1, &g1);
	d__1 = vx + mx + 1;
	gamma2_(&d__1, &g2);
	*pmv = *pmv * g1 / g2 * pow_ii(&c_n1, &mx);
    }
    return 0;
} /* lpmv_ */

/*       ********************************** */
/* Subroutine */ int cgama_(doublereal *x, doublereal *y, integer *kf, 
	doublereal *gr, doublereal *gi)
{
    /* Initialized data */

    static doublereal a[10] = { .08333333333333333,-.002777777777777778,
	    7.936507936507937e-4,-5.952380952380952e-4,8.417508417508418e-4,
	    -.001917526917526918,.00641025641025641,-.02955065359477124,
	    .1796443723688307,-1.3924322169059 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal), log(doublereal), pow_di(
	    doublereal *, integer *), cos(doublereal), sin(doublereal), cosh(
	    doublereal), sinh(doublereal), exp(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal t, g0, x0, x1, y1, z1, z2;
    static integer na;
    static doublereal pi, th, si, sr, gi1, gr1, th1, th2;


/*       ========================================================= */
/*       Purpose: Compute the gamma function Г(z) or ln[Г(z)] */
/*                for a complex argument */
/*       Input :  x  --- Real part of z */
/*                y  --- Imaginary part of z */
/*                KF --- Function code */
/*                       KF=0 for ln[Г(z)] */
/*                       KF=1 for Г(z) */
/*       Output:  GR --- Real part of ln[Г(z)] or Г(z) */
/*                GI --- Imaginary part of ln[Г(z)] or Г(z) */
/*       ======================================================== */

    pi = 3.141592653589793;
    if (*y == 0. && *x == (doublereal) ((integer) (*x)) && *x <= 0.) {
	*gr = 1e300;
	*gi = 0.;
	return 0;
    } else if (*x < 0.) {
	x1 = *x;
	y1 = *y;
	*x = -(*x);
	*y = -(*y);
    } else {
	y1 = 0.;
	x1 = *x;
    }
    x0 = *x;
    na = 0;
    if (*x <= 7.f) {
	na = (integer) (7 - *x);
	x0 = *x + na;
    }
    z1 = sqrt(x0 * x0 + *y * *y);
    th = atan(*y / x0);
    *gr = (x0 - .5) * log(z1) - th * *y - x0 + log(pi * 2.) * .5;
    *gi = th * (x0 - .5) + *y * log(z1) - *y;
    for (k = 1; k <= 10; ++k) {
	i__1 = 1 - (k << 1);
	t = pow_di(&z1, &i__1);
	*gr += a[k - 1] * t * cos((k * 2. - 1.) * th);
/* L10: */
	*gi -= a[k - 1] * t * sin((k * 2. - 1.) * th);
    }
    if (*x <= 7.f) {
	gr1 = 0.;
	gi1 = 0.;
	i__1 = na - 1;
	for (j = 0; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = *x + j;
	    gr1 += log(d__1 * d__1 + *y * *y) * .5;
/* L15: */
	    gi1 += atan(*y / (*x + j));
	}
	*gr -= gr1;
	*gi -= gi1;
    }
    if (x1 < 0.) {
	z1 = sqrt(*x * *x + *y * *y);
	th1 = atan(*y / *x);
	sr = -sin(pi * *x) * cosh(pi * *y);
	si = -cos(pi * *x) * sinh(pi * *y);
	z2 = sqrt(sr * sr + si * si);
	th2 = atan(si / sr);
	if (sr < 0.) {
	    th2 = pi + th2;
	}
	*gr = log(pi / (z1 * z2)) - *gr;
	*gi = -th1 - th2 - *gi;
	*x = x1;
	*y = y1;
    }
    if (*kf == 1) {
	g0 = exp(*gr);
	*gr = g0 * cos(*gi);
	*gi = g0 * sin(*gi);
    }
    return 0;
} /* cgama_ */

/*       ********************************** */
/* Subroutine */ int aswfb_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, integer *kd, doublereal *cv, doublereal *s1f, 
	doublereal *s1d)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer k;
    static doublereal df[200], pd[252];
    static integer mk, ip, nm;
    static doublereal pm[252], sw;
    static integer nm2;
    static doublereal su1, eps;
    extern /* Subroutine */ int sdmn_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *), lpmns_(integer *, integer 
	    *, doublereal *, doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute the prolate and oblate spheroidal angular */
/*                functions of the first kind and their derivatives */
/*       Input :  m  --- Mode parameter,  m = 0,1,2,... */
/*                n  --- Mode parameter,  n = m,m+1,... */
/*                c  --- Spheroidal parameter */
/*                x  --- Argument of angular function, |x| < 1.0 */
/*                KD --- Function code */
/*                       KD=1 for prolate;  KD=-1 for oblate */
/*                cv --- Characteristic value */
/*       Output:  S1F --- Angular function of the first kind */
/*                S1D --- Derivative of the angular function of */
/*                        the first kind */
/*       Routines called: */
/*            (1) SDMN for computing expansion coefficients dk */
/*            (2) LPMNS for computing associated Legendre function */
/*                of the first kind Pmn(x) */
/*       =========================================================== */

    eps = 1e-14;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    nm = (integer) ((*n - *m) / 2 + *c__) + 25;
    nm2 = (nm << 1) + *m;
    sdmn_(m, n, c__, cv, kd, df);
    lpmns_(m, &nm2, x, pm, pd);
    sw = 0.;
    su1 = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	mk = *m + (k - 1 << 1) + ip;
	su1 += df[k - 1] * pm[mk];
	if ((d__1 = sw - su1, abs(d__1)) < abs(su1) * eps) {
	    goto L15;
	}
/* L10: */
	sw = su1;
    }
L15:
    *s1f = pow_ii(&c_n1, m) * su1;
    su1 = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	mk = *m + (k - 1 << 1) + ip;
	su1 += df[k - 1] * pd[mk];
	if ((d__1 = sw - su1, abs(d__1)) < abs(su1) * eps) {
	    goto L25;
	}
/* L20: */
	sw = su1;
    }
L25:
    *s1d = pow_ii(&c_n1, m) * su1;
    return 0;
} /* aswfb_ */

/*       ********************************** */
/* Subroutine */ int chgus_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *hu, integer *id)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), pow_dd(doublereal *, doublereal *), d_lg10(
	    doublereal *);

    /* Local variables */
    static integer j;
    static doublereal d1, d2, h0, r1, r2, ga, gb, pi, gb2, hu0, xg1, xg2, gab,
	     hua, hmin, hmax;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute confluent hypergeometric function */
/*                U(a,b,x) for small argument x */
/*       Input  : a  --- Parameter */
/*                b  --- Parameter ( b <> 0,-1,-2,...) */
/*                x  --- Argument */
/*       Output:  HU --- U(a,b,x) */
/*                ID --- Estimated number of significant digits */
/*       Routine called: GAMMA2 for computing gamma function */
/*       ====================================================== */

/*       DLMF 13.2.42 with prefactors rewritten according to */
/*       DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2 */

    *id = -100;
    pi = 3.141592653589793;
    gamma2_(a, &ga);
    gamma2_(b, &gb);
    xg1 = *a + 1. - *b;
    gamma2_(&xg1, &gab);
    xg2 = 2. - *b;
    gamma2_(&xg2, &gb2);
    hu0 = pi / sin(pi * *b);
    r1 = hu0 / (gab * gb);
    d__1 = 1. - *b;
    r2 = hu0 * pow_dd(x, &d__1) / (ga * gb2);
    *hu = r1 - r2;
    hmax = 0.;
    hmin = 1e300;
    h0 = 0.;
    for (j = 1; j <= 150; ++j) {
	r1 = r1 * (*a + j - 1.) / (j * (*b + j - 1.)) * *x;
	r2 = r2 * (*a - *b + j) / (j * (1. - *b + j)) * *x;
	*hu = *hu + r1 - r2;
	hua = abs(*hu);
	if (hua > hmax) {
	    hmax = hua;
	}
	if (hua < hmin) {
	    hmin = hua;
	}
	if ((d__1 = *hu - h0, abs(d__1)) < abs(*hu) * 1e-15) {
	    goto L15;
	}
/* L10: */
	h0 = *hu;
    }
L15:
    d1 = d_lg10(&hmax);
    d2 = 0.;
    if (hmin != 0.f) {
	d2 = d_lg10(&hmin);
    }
    *id = (integer) (15 - (d__1 = d1 - d2, abs(d__1)));
    return 0;
} /* chgus_ */

/*       ********************************** */
/* Subroutine */ int itth0_(doublereal *x, doublereal *tth)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, s, t, f0, g0, pi, xt, tty;


/*       =========================================================== */
/*       Purpose: Evaluate the integral H0(t)/t with respect to t */
/*                from x to infinity */
/*       Input :  x   --- Lower limit  ( x ≥ 0 ) */
/*       Output:  TTH --- Integration of H0(t)/t from x to infinity */
/*       =========================================================== */

    pi = 3.141592653589793;
    s = 1.;
    r__ = 1.;
    if (*x < 24.5) {
	for (k = 1; k <= 60; ++k) {
/* Computing 3rd power */
	    d__1 = k * 2.f + 1.;
	    r__ = -r__ * *x * *x * (k * 2.f - 1.) / (d__1 * (d__1 * d__1));
	    s += r__;
	    if (abs(r__) < abs(s) * 1e-12) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*tth = pi / 2. - 2. / pi * *x * s;
    } else {
	for (k = 1; k <= 10; ++k) {
/* Computing 3rd power */
	    d__1 = k * 2.f - 1.;
	    r__ = -r__ * (d__1 * (d__1 * d__1)) / ((k * 2.f + 1.) * *x * *x);
	    s += r__;
	    if (abs(r__) < abs(s) * 1e-12) {
		goto L25;
	    }
/* L20: */
	}
L25:
	*tth = 2. / (pi * *x) * s;
	t = 8. / *x;
	xt = *x + pi * .25;
	f0 = (((((t * .0018118 - .0091909) * t + .017033) * t - 9.394e-4) * t 
		- .051445) * t - 1.1e-6) * t + .7978846;
	g0 = (((((t * -.0023731 + .0059842) * t + .0024437) * t - .0233178) * 
		t + 5.95e-5) * t + .1620695) * t;
	tty = (f0 * sin(xt) - g0 * cos(xt)) / (sqrt(*x) * *x);
	*tth += tty;
    }
    return 0;
} /* itth0_ */

/*       ********************************** */
/* Subroutine */ int lgama_(integer *kf, doublereal *x, doublereal *gl)
{
    /* Initialized data */

    static doublereal a[10] = { .08333333333333333,-.002777777777777778,
	    7.936507936507937e-4,-5.952380952380952e-4,8.417508417508418e-4,
	    -.001917526917526918,.00641025641025641,-.02955065359477124,
	    .1796443723688307,-1.3924322169059 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer k, n;
    static doublereal x0, x2, xp, gl0;


/*       ================================================== */
/*       Purpose: Compute gamma function Г(x) or ln[Г(x)] */
/*       Input:   x  --- Argument of Г(x) ( x > 0 ) */
/*                KF --- Function code */
/*                       KF=1 for Г(x); KF=0 for ln[Г(x)] */
/*       Output:  GL --- Г(x) or ln[Г(x)] */
/*       ================================================== */

    x0 = *x;
    n = 0;
    if (*x == 1.f || *x == 2.f) {
	*gl = 0.;
	goto L20;
    } else if (*x <= 7.f) {
	n = (integer) (7 - *x);
	x0 = *x + n;
    }
    x2 = 1. / (x0 * x0);
    xp = 6.283185307179586477;
    gl0 = a[9];
    for (k = 9; k >= 1; --k) {
/* L10: */
	gl0 = gl0 * x2 + a[k - 1];
    }
    *gl = gl0 / x0 + log(xp) * .5 + (x0 - .5) * log(x0) - x0;
    if (*x <= 7.f) {
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
	    *gl -= log(x0 - 1.);
/* L15: */
	    x0 += -1.;
	}
    }
L20:
    if (*kf == 1) {
	*gl = exp(*gl);
    }
    return 0;
} /* lgama_ */

/*       ********************************** */
/* Subroutine */ int lqna_(integer *n, doublereal *x, doublereal *qn, 
	doublereal *qd)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer k;
    static doublereal q0, q1, qf;


/*       ===================================================== */
/*       Purpose: Compute Legendre functions Qn(x) and Qn'(x) */
/*       Input :  x  --- Argument of Qn(x) ( -1 ≤ x ≤ 1 ) */
/*                n  --- Degree of Qn(x) ( n = 0,1,2,… ) */
/*       Output:  QN(n) --- Qn(x) */
/*                QD(n) --- Qn'(x) */
/*                ( 1.0D+300 stands for infinity ) */
/*       ===================================================== */

    if (abs(*x) == 1.) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    qn[k] = 1e300;
	    qd[k] = -1e300;
/* L10: */
	}
    } else if (abs(*x) < 1.) {
	q0 = log((*x + 1.) / (1. - *x)) * .5;
	q1 = *x * q0 - 1.;
	qn[0] = q0;
	qn[1] = q1;
	qd[0] = 1. / (1. - *x * *x);
	qd[1] = qn[0] + *x * qd[0];
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    qf = (((k << 1) - 1) * *x * q1 - (k - 1) * q0) / k;
	    qn[k] = qf;
	    qd[k] = (qn[k - 1] - *x * qf) * k / (1. - *x * *x);
	    q0 = q1;
/* L15: */
	    q1 = qf;
	}
    }
    return 0;
} /* lqna_ */

/*       ********************************** */
/* Subroutine */ int dvla_(doublereal *va, doublereal *x, doublereal *pd)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), pow_dd(doublereal *, doublereal *), cos(
	    doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__, a0, x1, gl, ep, pi, vl, eps;
    extern /* Subroutine */ int vvla_(doublereal *, doublereal *, doublereal *
	    ), gamma2_(doublereal *, doublereal *);


/*       ==================================================== */
/*       Purpose: Compute parabolic cylinder functions Dv(x) */
/*                for large argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PD --- Dv(x) */
/*       Routines called: */
/*             (1) VVLA for computing Vv(x) for large |x| */
/*             (2) GAMMA2 for computing Г(x) */
/*       ==================================================== */

    pi = 3.141592653589793;
    eps = 1e-12;
    ep = exp(*x * -.25f * *x);
    d__1 = abs(*x);
    a0 = pow_dd(&d__1, va) * ep;
    r__ = 1.;
    *pd = 1.;
    for (k = 1; k <= 16; ++k) {
	r__ = r__ * -.5 * (k * 2.f - *va - 1.f) * (k * 2.f - *va - 2.f) / (k *
		 *x * *x);
	*pd += r__;
	if ((d__1 = r__ / *pd, abs(d__1)) < eps) {
	    goto L15;
	}
/* L10: */
    }
L15:
    *pd = a0 * *pd;
    if (*x < 0.) {
	x1 = -(*x);
	vvla_(va, &x1, &vl);
	d__1 = -(*va);
	gamma2_(&d__1, &gl);
	*pd = pi * vl / gl + cos(pi * *va) * *pd;
    }
    return 0;
} /* dvla_ */

/*       ********************************** */
/* Subroutine */ int ik01a_(doublereal *x, doublereal *bi0, doublereal *di0, 
	doublereal *bi1, doublereal *di1, doublereal *bk0, doublereal *dk0, 
	doublereal *bk1, doublereal *dk1)
{
    /* Initialized data */

    static doublereal a[12] = { .125,.0703125,.0732421875,.11215209960938,
	    .22710800170898,.57250142097473,1.7277275025845,6.0740420012735,
	    24.380529699556,110.01714026925,551.33589612202,3038.0905109224 };
    static doublereal b[12] = { -.375,-.1171875,-.1025390625,-.14419555664063,
	    -.2775764465332,-.67659258842468,-1.9935317337513,
	    -6.8839142681099,-27.248827311269,-121.59789187654,
	    -603.84407670507,-3302.2722944809 };
    static doublereal a1[8] = { .125,.2109375,1.0986328125,11.775970458984,
	    214.61706161499,5951.1522710323,233476.45606175,12312234.987631 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), pow_di(doublereal *, integer *),
	     log(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r__;
    static integer k0;
    static doublereal w0, x2, ca, cb, el, ct, pi, xr, ww, xr2;


/*       ========================================================= */
/*       Purpose: Compute modified Bessel functions I0(x), I1(1), */
/*                K0(x) and K1(x), and their derivatives */
/*       Input :  x   --- Argument ( x ≥ 0 ) */
/*       Output:  BI0 --- I0(x) */
/*                DI0 --- I0'(x) */
/*                BI1 --- I1(x) */
/*                DI1 --- I1'(x) */
/*                BK0 --- K0(x) */
/*                DK0 --- K0'(x) */
/*                BK1 --- K1(x) */
/*                DK1 --- K1'(x) */
/*       ========================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    x2 = *x * *x;
    if (*x == 0.) {
	*bi0 = 1.;
	*bi1 = 0.;
	*bk0 = 1e300;
	*bk1 = 1e300;
	*di0 = 0.;
	*di1 = .5;
	*dk0 = -1e300;
	*dk1 = -1e300;
	return 0;
    } else if (*x <= 18.) {
	*bi0 = 1.;
	r__ = 1.;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * .25 * x2 / (k * k);
	    *bi0 += r__;
	    if ((d__1 = r__ / *bi0, abs(d__1)) < 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	*bi1 = 1.;
	r__ = 1.;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * .25 * x2 / (k * (k + 1));
	    *bi1 += r__;
	    if ((d__1 = r__ / *bi1, abs(d__1)) < 1e-15) {
		goto L30;
	    }
/* L25: */
	}
L30:
	*bi1 = *x * .5 * *bi1;
    } else {
	k0 = 12;
	if (*x >= 35.f) {
	    k0 = 9;
	}
	if (*x >= 50.f) {
	    k0 = 7;
	}
	ca = exp(*x) / sqrt(pi * 2. * *x);
	*bi0 = 1.;
	xr = 1. / *x;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
/* L35: */
	    *bi0 += a[k - 1] * pow_di(&xr, &k);
	}
	*bi0 = ca * *bi0;
	*bi1 = 1.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
/* L40: */
	    *bi1 += b[k - 1] * pow_di(&xr, &k);
	}
	*bi1 = ca * *bi1;
    }
    ww = 0.;
    if (*x <= 9.) {
	ct = -(log(*x / 2.) + el);
	*bk0 = 0.;
	w0 = 0.;
	r__ = 1.;
	for (k = 1; k <= 50; ++k) {
	    w0 += 1. / k;
	    r__ = r__ * .25 / (k * k) * x2;
	    *bk0 += r__ * (w0 + ct);
	    if ((d__1 = (*bk0 - ww) / *bk0, abs(d__1)) < 1e-15) {
		goto L70;
	    }
/* L65: */
	    ww = *bk0;
	}
L70:
	*bk0 += ct;
    } else {
	cb = .5 / *x;
	xr2 = 1. / x2;
	*bk0 = 1.;
	for (k = 1; k <= 8; ++k) {
/* L75: */
	    *bk0 += a1[k - 1] * pow_di(&xr2, &k);
	}
	*bk0 = cb * *bk0 / *bi0;
    }
    *bk1 = (1. / *x - *bi1 * *bk0) / *bi0;
    *di0 = *bi1;
    *di1 = *bi0 - *bi1 / *x;
    *dk0 = -(*bk1);
    *dk1 = -(*bk0) - *bk1 / *x;
    return 0;
} /* ik01a_ */

/*       ********************************** */
/* Subroutine */ int cpbdn_(integer *n, doublecomplex *z__, doublecomplex *
	cpb, doublecomplex *cpd)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_exp(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, m;
    static doublereal x, a0;
    static doublecomplex c0;
    static integer n0, n1;
    static doublecomplex z1, cf;
    static doublereal pi;
    static doublecomplex ca0, cf0, cf1, cs0;
    static integer nm1;
    static doublecomplex cfa, cfb;
    extern /* Subroutine */ int cpdla_(integer *, doublecomplex *, 
	    doublecomplex *), cpdsa_(integer *, doublecomplex *, 
	    doublecomplex *);


/*       ================================================== */
/*       Purpose: Compute the parabolic cylinder functions */
/*                 Dn(z) and Dn'(z) for a complex argument */
/*       Input:   z --- Complex argument of Dn(z) */
/*                n --- Order of Dn(z)  ( n=0,±1,±2,… ) */
/*       Output:  CPB(|n|) --- Dn(z) */
/*                CPD(|n|) --- Dn'(z) */
/*       Routines called: */
/*            (1) CPDSA for computing Dn(z) for a small |z| */
/*            (2) CPDLA for computing Dn(z) for a large |z| */
/*       ================================================== */

    pi = 3.141592653589793;
    x = z__->r;
    a0 = z_abs(z__);
    c0.r = 0., c0.i = 0.;
    z__3.r = z__->r * -.25, z__3.i = z__->i * -.25;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z_exp(&z__1, &z__2);
    ca0.r = z__1.r, ca0.i = z__1.i;
    n0 = 0;
    if (*n >= 0) {
	cf0.r = ca0.r, cf0.i = ca0.i;
	z__1.r = z__->r * ca0.r - z__->i * ca0.i, z__1.i = z__->r * ca0.i + 
		z__->i * ca0.r;
	cf1.r = z__1.r, cf1.i = z__1.i;
	cpb[0].r = cf0.r, cpb[0].i = cf0.i;
	cpb[1].r = cf1.r, cpb[1].i = cf1.i;
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    z__2.r = z__->r * cf1.r - z__->i * cf1.i, z__2.i = z__->r * cf1.i 
		    + z__->i * cf1.r;
	    d__1 = k - 1.;
	    z__3.r = d__1 * cf0.r, z__3.i = d__1 * cf0.i;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	    i__2 = k;
	    cpb[i__2].r = cf.r, cpb[i__2].i = cf.i;
	    cf0.r = cf1.r, cf0.i = cf1.i;
/* L10: */
	    cf1.r = cf.r, cf1.i = cf.i;
	}
    } else {
	n0 = -(*n);
	if (x <= 0.f || z_abs(z__) == 0.f) {
	    cf0.r = ca0.r, cf0.i = ca0.i;
	    cpb[0].r = cf0.r, cpb[0].i = cf0.i;
	    z__1.r = -z__->r, z__1.i = -z__->i;
	    z1.r = z__1.r, z1.i = z__1.i;
	    if (a0 <= 7.f) {
		cpdsa_(&c_n1, &z1, &cf1);
	    } else {
		cpdla_(&c_n1, &z1, &cf1);
	    }
	    d__1 = sqrt(pi * 2.);
	    z__3.r = d__1, z__3.i = 0.;
	    z_div(&z__2, &z__3, &ca0);
	    z__1.r = z__2.r - cf1.r, z__1.i = z__2.i - cf1.i;
	    cf1.r = z__1.r, cf1.i = z__1.i;
	    cpb[1].r = cf1.r, cpb[1].i = cf1.i;
	    i__1 = n0;
	    for (k = 2; k <= i__1; ++k) {
		z__4.r = -z__->r, z__4.i = -z__->i;
		z__3.r = z__4.r * cf1.r - z__4.i * cf1.i, z__3.i = z__4.r * 
			cf1.i + z__4.i * cf1.r;
		z__2.r = z__3.r + cf0.r, z__2.i = z__3.i + cf0.i;
		d__1 = k - 1.;
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cf.r = z__1.r, cf.i = z__1.i;
		i__2 = k;
		cpb[i__2].r = cf.r, cpb[i__2].i = cf.i;
		cf0.r = cf1.r, cf0.i = cf1.i;
/* L15: */
		cf1.r = cf.r, cf1.i = cf.i;
	    }
	} else {
	    if (a0 <= 3.f) {
		i__1 = -n0;
		cpdsa_(&i__1, z__, &cfa);
		i__1 = n0;
		cpb[i__1].r = cfa.r, cpb[i__1].i = cfa.i;
		n1 = n0 + 1;
		i__1 = -n1;
		cpdsa_(&i__1, z__, &cfb);
		i__1 = n1;
		cpb[i__1].r = cfb.r, cpb[i__1].i = cfb.i;
		nm1 = n0 - 1;
		for (k = nm1; k >= 0; --k) {
		    z__2.r = z__->r * cfa.r - z__->i * cfa.i, z__2.i = z__->r 
			    * cfa.i + z__->i * cfa.r;
		    d__1 = k + 1.;
		    z__3.r = d__1 * cfb.r, z__3.i = d__1 * cfb.i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    cf.r = z__1.r, cf.i = z__1.i;
		    i__1 = k;
		    cpb[i__1].r = cf.r, cpb[i__1].i = cf.i;
		    cfb.r = cfa.r, cfb.i = cfa.i;
/* L20: */
		    cfa.r = cf.r, cfa.i = cf.i;
		}
	    } else {
		m = abs(*n) + 100;
		cfa.r = c0.r, cfa.i = c0.i;
		cfb.r = 1e-30, cfb.i = 0.;
		for (k = m; k >= 0; --k) {
		    z__2.r = z__->r * cfb.r - z__->i * cfb.i, z__2.i = z__->r 
			    * cfb.i + z__->i * cfb.r;
		    d__1 = k + 1.;
		    z__3.r = d__1 * cfa.r, z__3.i = d__1 * cfa.i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    cf.r = z__1.r, cf.i = z__1.i;
		    if (k <= n0) {
			i__1 = k;
			cpb[i__1].r = cf.r, cpb[i__1].i = cf.i;
		    }
		    cfa.r = cfb.r, cfa.i = cfb.i;
/* L25: */
		    cfb.r = cf.r, cfb.i = cf.i;
		}
		z_div(&z__1, &ca0, &cf);
		cs0.r = z__1.r, cs0.i = z__1.i;
		i__1 = n0;
		for (k = 0; k <= i__1; ++k) {
/* L30: */
		    i__2 = k;
		    i__3 = k;
		    z__1.r = cs0.r * cpb[i__3].r - cs0.i * cpb[i__3].i, 
			    z__1.i = cs0.r * cpb[i__3].i + cs0.i * cpb[i__3]
			    .r;
		    cpb[i__2].r = z__1.r, cpb[i__2].i = z__1.i;
		}
	    }
	}
    }
    z__2.r = z__->r * -.5, z__2.i = z__->i * -.5;
    z__1.r = z__2.r * cpb[0].r - z__2.i * cpb[0].i, z__1.i = z__2.r * cpb[0]
	    .i + z__2.i * cpb[0].r;
    cpd[0].r = z__1.r, cpd[0].i = z__1.i;
    if (*n >= 0) {
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
/* L35: */
	    i__3 = k;
	    z__3.r = z__->r * -.5, z__3.i = z__->i * -.5;
	    i__1 = k;
	    z__2.r = z__3.r * cpb[i__1].r - z__3.i * cpb[i__1].i, z__2.i = 
		    z__3.r * cpb[i__1].i + z__3.i * cpb[i__1].r;
	    i__4 = k - 1;
	    d__1 = (doublereal) k;
	    z__4.r = d__1 * cpb[i__4].r, z__4.i = d__1 * cpb[i__4].i;
	    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
	    cpd[i__3].r = z__1.r, cpd[i__3].i = z__1.i;
	}
    } else {
	i__3 = n0;
	for (k = 1; k <= i__3; ++k) {
/* L40: */
	    i__1 = k;
	    z__3.r = z__->r * .5, z__3.i = z__->i * .5;
	    i__4 = k;
	    z__2.r = z__3.r * cpb[i__4].r - z__3.i * cpb[i__4].i, z__2.i = 
		    z__3.r * cpb[i__4].i + z__3.i * cpb[i__4].r;
	    i__2 = k - 1;
	    z__1.r = z__2.r - cpb[i__2].r, z__1.i = z__2.i - cpb[i__2].i;
	    cpd[i__1].r = z__1.r, cpd[i__1].i = z__1.i;
	}
    }
    return 0;
} /* cpbdn_ */

/*       ********************************** */
/* Subroutine */ int ik01b_(doublereal *x, doublereal *bi0, doublereal *di0, 
	doublereal *bi1, doublereal *di1, doublereal *bk0, doublereal *dk0, 
	doublereal *bk1, doublereal *dk1)
{
    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), log(doublereal);

    /* Local variables */
    static doublereal t, t2;


/*       ========================================================= */
/*       Purpose: Compute modified Bessel functions I0(x), I1(1), */
/*                K0(x) and K1(x), and their derivatives */
/*       Input :  x   --- Argument ( x ≥ 0 ) */
/*       Output:  BI0 --- I0(x) */
/*                DI0 --- I0'(x) */
/*                BI1 --- I1(x) */
/*                DI1 --- I1'(x) */
/*                BK0 --- K0(x) */
/*                DK0 --- K0'(x) */
/*                BK1 --- K1(x) */
/*                DK1 --- K1'(x) */
/*       ========================================================= */

    if (*x == 0.) {
	*bi0 = 1.;
	*bi1 = 0.;
	*bk0 = 1e300;
	*bk1 = 1e300;
	*di0 = 0.;
	*di1 = .5;
	*dk0 = -1e300;
	*dk1 = -1e300;
	return 0;
    } else if (*x <= 3.75) {
	t = *x / 3.75;
	t2 = t * t;
	*bi0 = (((((t2 * .0045813 + .0360768) * t2 + .2659732) * t2 + 
		1.2067492) * t2 + 3.0899424) * t2 + 3.5156229) * t2 + 1.;
	*bi1 = *x * ((((((t2 * 3.2411e-4 + .00301532) * t2 + .02658733) * t2 
		+ .15084934) * t2 + .51498869) * t2 + .87890594) * t2 + .5);
    } else {
	t = 3.75 / *x;
	*bi0 = ((((((((t * .00392377 - .01647633) * t + .02635537) * t - 
		.02057706) * t + .00916281) * t - .00157565) * t + .00225319) 
		* t + .01328592) * t + .39894228) * exp(*x) / sqrt(*x);
	*bi1 = ((((((((t * -.00420059 + .01787654) * t - .02895312) * t + 
		.02282967) * t - .01031555) * t + .00163801) * t - .00362018) 
		* t - .03988024) * t + .39894228) * exp(*x) / sqrt(*x);
    }
    if (*x <= 2.) {
	t = *x / 2.;
	t2 = t * t;
	*bk0 = (((((t2 * 7.4e-6 + 1.075e-4) * t2 + .00262698) * t2 + .0348859)
		 * t2 + .23069756) * t2 + .4227842) * t2 - .57721566 - *bi0 * 
		log(t);
	*bk1 = ((((((t2 * -4.686e-5 - .00110404) * t2 - .01919402) * t2 - 
		.18156897) * t2 - .67278579) * t2 + .15443144) * t2 + 1.) / *
		x + *bi1 * log(t);
    } else {
	t = 2. / *x;
	t2 = t * t;
	*bk0 = ((((((t * 5.3208e-4 - .0025154) * t + .00587872) * t - 
		.01062446) * t + .02189568) * t - .07832358) * t + 1.25331414)
		 * exp(-(*x)) / sqrt(*x);
	*bk1 = ((((((t * -6.8245e-4 + .00325614) * t - .00780353) * t + 
		.01504268) * t - .0365562) * t + .23498619) * t + 1.25331414) 
		* exp(-(*x)) / sqrt(*x);
    }
    *di0 = *bi1;
    *di1 = *bi0 - *bi1 / *x;
    *dk0 = -(*bk1);
    *dk1 = -(*bk0) - *bk1 / *x;
    return 0;
} /* ik01b_ */

/*       ********************************** */
/* Subroutine */ int beta_(doublereal *p, doublereal *q, doublereal *bt)
{
    static doublereal gp, gq, gpq, ppq;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ========================================== */
/*       Purpose: Compute the beta function B(p,q) */
/*       Input :  p  --- Parameter  ( p > 0 ) */
/*                q  --- Parameter  ( q > 0 ) */
/*       Output:  BT --- B(p,q) */
/*       Routine called: GAMMA2 for computing Г(x) */
/*       ========================================== */

    gamma2_(p, &gp);
    gamma2_(q, &gq);
    ppq = *p + *q;
    gamma2_(&ppq, &gpq);
    *bt = gp * gq / gpq;
    return 0;
} /* beta_ */

/*       ********************************** */
/* Subroutine */ int lpn_(integer *n, doublereal *x, doublereal *pn, 
	doublereal *pd)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer k;
    static doublereal p0, p1, pf;


/*       =============================================== */
/*       Purpose: Compute Legendre polynomials Pn(x) */
/*                and their derivatives Pn'(x) */
/*       Input :  x --- Argument of Pn(x) */
/*                n --- Degree of Pn(x) ( n = 0,1,...) */
/*       Output:  PN(n) --- Pn(x) */
/*                PD(n) --- Pn'(x) */
/*       =============================================== */

    pn[0] = 1.;
    pn[1] = *x;
    pd[0] = 0.;
    pd[1] = 1.;
    p0 = 1.;
    p1 = *x;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	pf = (k * 2. - 1.) / k * *x * p1 - (k - 1.) / k * p0;
	pn[k] = pf;
	if (abs(*x) == 1.) {
	    i__2 = k + 1;
	    pd[k] = pow_di(x, &i__2) * .5 * k * (k + 1.);
	} else {
	    pd[k] = k * (p1 - *x * pf) / (1. - *x * *x);
	}
	p0 = p1;
/* L10: */
	p1 = pf;
    }
    return 0;
} /* lpn_ */

/*       ********************************** */
/* Subroutine */ int fcoef_(integer *kd, integer *m, doublereal *q, 
	doublereal *a, doublereal *fc)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal f;
    static integer i__, j, k;
    static doublereal s, u, v, f1, f2, f3, s0;
    static integer kb, jm, km;
    static doublereal qm, sp, ss;
    extern doublereal dnan_(void);
    static doublereal fnan;


/*       ===================================================== */
/*       Purpose: Compute expansion coefficients for Mathieu */
/*                functions and modified Mathieu functions */
/*       Input :  m  --- Order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions */
/*                KD --- Case code */
/*                       KD=1 for cem(x,q)  ( m = 0,2,4,...) */
/*                       KD=2 for cem(x,q)  ( m = 1,3,5,...) */
/*                       KD=3 for sem(x,q)  ( m = 1,3,5,...) */
/*                       KD=4 for sem(x,q)  ( m = 2,4,6,...) */
/*                A  --- Characteristic value of Mathieu */
/*                       functions for given m and q */
/*       Output:  FC(k) --- Expansion coefficients of Mathieu */
/*                       functions ( k= 1,2,...,KM ) */
/*                       FC(1),FC(2),FC(3),... correspond to */
/*                       A0,A2,A4,... for KD=1 case, A1,A3, */
/*                       A5,... for KD=2 case, B1,B3,B5,... */
/*                       for KD=3 case and B2,B4,B6,... for */
/*                       KD=4 case */
/*       ===================================================== */

    /* Parameter adjustments */
    --fc;

    /* Function Body */
    for (i__ = 1; i__ <= 251; ++i__) {
/* L5: */
	fc[i__] = 0.;
    }
    if (abs(*q) <= 1e-7) {
/*          Expansion up to order Q^1 (Abramowitz & Stegun 20.2.27-28) */
	if (*kd == 1) {
	    jm = *m / 2 + 1;
	} else if (*kd == 2 || *kd == 3) {
	    jm = (*m - 1) / 2 + 1;
	} else if (*kd == 4) {
	    jm = *m / 2;
	}
/*          Check for overflow */
	if (jm + 1 > 251) {
	    goto L6;
	}
/*          Proceed using the simplest expansion */
	if (*kd == 1 || *kd == 2) {
	    if (*m == 0) {
		fc[1] = 1 / sqrt(2.);
		fc[2] = -(*q) / 2. / sqrt(2.);
	    } else if (*m == 1) {
		fc[1] = 1.;
		fc[2] = -(*q) / 8.;
	    } else if (*m == 2) {
		fc[1] = *q / 4.;
		fc[2] = 1.;
		fc[3] = -(*q) / 12.;
	    } else {
		fc[jm] = 1.;
		fc[jm + 1] = -(*q) / ((*m + 1) * 4.);
		fc[jm - 1] = *q / ((*m - 1) * 4.);
	    }
	} else if (*kd == 3 || *kd == 4) {
	    if (*m == 1) {
		fc[1] = 1.;
		fc[2] = -(*q) / 8.;
	    } else if (*m == 2) {
		fc[1] = 1.;
		fc[2] = -(*q) / 12.;
	    } else {
		fc[jm] = 1.;
		fc[jm + 1] = -(*q) / ((*m + 1) * 4.);
		fc[jm - 1] = *q / ((*m - 1) * 4.);
	    }
	}
	return 0;
    } else if (*q <= 1.) {
	qm = sqrt(*q) * 56.1f + 7.5f - *q * 134.7f + sqrt(*q) * 90.7f * *q;
    } else {
	qm = sqrt(*q) * 3.1f + 17.f - *q * .126f + sqrt(*q) * .0037f * *q;
    }
    km = (integer) (qm + *m * .5f);
    if (km > 251) {
/*          Overflow, generate NaNs */
L6:
	fnan = dnan_();
	for (i__ = 1; i__ <= 251; ++i__) {
/* L7: */
	    fc[i__] = fnan;
	}
	return 0;
    }
    kb = 0;
    s = 0.;
    f = 1e-100;
    u = 0.;
    fc[km] = 0.;
    f2 = 0.;
    if (*kd == 1) {
	for (k = km; k >= 3; --k) {
	    v = u;
	    u = f;
	    f = (*a - k * 4. * k) * u / *q - v;
	    if (abs(f) < (d__1 = fc[k + 1], abs(d__1))) {
		kb = k;
		fc[1] = 1e-100;
		sp = 0.;
		f3 = fc[k + 1];
		fc[2] = *a / *q * fc[1];
		fc[3] = (*a - 4.) * fc[2] / *q - fc[1] * 2.;
		u = fc[2];
		f1 = fc[3];
		i__1 = kb;
		for (i__ = 3; i__ <= i__1; ++i__) {
		    v = u;
		    u = f1;
/* Computing 2nd power */
		    d__1 = i__ - 1.;
		    f1 = (*a - d__1 * d__1 * 4.) * u / *q - v;
		    fc[i__ + 1] = f1;
		    if (i__ == kb) {
			f2 = f1;
		    }
		    if (i__ != kb) {
			sp += f1 * f1;
		    }
/* L15: */
		}
/* Computing 2nd power */
		d__1 = fc[1];
/* Computing 2nd power */
		d__2 = fc[2];
/* Computing 2nd power */
		d__3 = fc[3];
		sp = sp + d__1 * d__1 * 2. + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
		d__1 = f3 / f2;
		ss = s + sp * (d__1 * d__1);
		s0 = sqrt(1. / ss);
		i__1 = km;
		for (j = 1; j <= i__1; ++j) {
		    if (j <= kb + 1) {
			fc[j] = s0 * fc[j] * f3 / f2;
		    } else {
			fc[j] = s0 * fc[j];
		    }
/* L20: */
		}
		goto L85;
	    } else {
		fc[k] = f;
		s += f * f;
	    }
/* L25: */
	}
	fc[2] = *q * fc[3] / (*a - 4. - *q * 2. * *q / *a);
	fc[1] = *q / *a * fc[2];
/* Computing 2nd power */
	d__1 = fc[1];
/* Computing 2nd power */
	d__2 = fc[2];
	s = s + d__1 * d__1 * 2. + d__2 * d__2;
	s0 = sqrt(1. / s);
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
/* L30: */
	    fc[k] = s0 * fc[k];
	}
    } else if (*kd == 2 || *kd == 3) {
	for (k = km; k >= 3; --k) {
	    v = u;
	    u = f;
/* Computing 2nd power */
	    d__1 = k * 2. - 1;
	    f = (*a - d__1 * d__1) * u / *q - v;
	    if (abs(f) >= (d__1 = fc[k], abs(d__1))) {
		fc[k - 1] = f;
		s += f * f;
	    } else {
		kb = k;
		f3 = fc[k];
		goto L45;
	    }
/* L35: */
	}
	fc[1] = *q / (*a - 1. - pow_ii(&c_n1, kd) * *q) * fc[2];
	s += fc[1] * fc[1];
	s0 = sqrt(1. / s);
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
/* L40: */
	    fc[k] = s0 * fc[k];
	}
	goto L85;
L45:
	fc[1] = 1e-100;
	fc[2] = (*a - 1. - pow_ii(&c_n1, kd) * *q) / *q * fc[1];
	sp = 0.;
	u = fc[1];
	f1 = fc[2];
	i__1 = kb - 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    v = u;
	    u = f1;
/* Computing 2nd power */
	    d__1 = i__ * 2. - 1.;
	    f1 = (*a - d__1 * d__1) * u / *q - v;
	    if (i__ != kb - 1) {
		fc[i__ + 1] = f1;
		sp += f1 * f1;
	    } else {
		f2 = f1;
	    }
/* L50: */
	}
/* Computing 2nd power */
	d__1 = fc[1];
/* Computing 2nd power */
	d__2 = fc[2];
	sp = sp + d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
	d__1 = f3 / f2;
	ss = s + sp * (d__1 * d__1);
	s0 = 1. / sqrt(ss);
	i__1 = km;
	for (j = 1; j <= i__1; ++j) {
	    if (j < kb) {
		fc[j] = s0 * fc[j] * f3 / f2;
	    }
	    if (j >= kb) {
		fc[j] = s0 * fc[j];
	    }
/* L55: */
	}
    } else if (*kd == 4) {
	for (k = km; k >= 3; --k) {
	    v = u;
	    u = f;
	    f = (*a - k * 4. * k) * u / *q - v;
	    if (abs(f) >= (d__1 = fc[k], abs(d__1))) {
		fc[k - 1] = f;
		s += f * f;
	    } else {
		kb = k;
		f3 = fc[k];
		goto L70;
	    }
/* L60: */
	}
	fc[1] = *q / (*a - 4.) * fc[2];
	s += fc[1] * fc[1];
	s0 = sqrt(1. / s);
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
/* L65: */
	    fc[k] = s0 * fc[k];
	}
	goto L85;
L70:
	fc[1] = 1e-100;
	fc[2] = (*a - 4.) / *q * fc[1];
	sp = 0.;
	u = fc[1];
	f1 = fc[2];
	i__1 = kb - 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    v = u;
	    u = f1;
	    f1 = (*a - i__ * 4. * i__) * u / *q - v;
	    if (i__ != kb - 1) {
		fc[i__ + 1] = f1;
		sp += f1 * f1;
	    } else {
		f2 = f1;
	    }
/* L75: */
	}
/* Computing 2nd power */
	d__1 = fc[1];
/* Computing 2nd power */
	d__2 = fc[2];
	sp = sp + d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
	d__1 = f3 / f2;
	ss = s + sp * (d__1 * d__1);
	s0 = 1. / sqrt(ss);
	i__1 = km;
	for (j = 1; j <= i__1; ++j) {
	    if (j < kb) {
		fc[j] = s0 * fc[j] * f3 / f2;
	    }
	    if (j >= kb) {
		fc[j] = s0 * fc[j];
	    }
/* L80: */
	}
    }
L85:
    if (fc[1] < 0.) {
	i__1 = km;
	for (j = 1; j <= i__1; ++j) {
/* L90: */
	    fc[j] = -fc[j];
	}
    }
    return 0;
} /* fcoef_ */

/*       ********************************** */
/* Subroutine */ int sphi_(integer *n, doublereal *x, integer *nm, doublereal 
	*si, doublereal *di)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sinh(doublereal), cosh(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, m;
    static doublereal f0, f1, cs, si0;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ======================================================== */
/*       Purpose: Compute modified spherical Bessel functions */
/*                of the first kind, in(x) and in'(x) */
/*       Input :  x --- Argument of in(x) */
/*                n --- Order of in(x) ( n = 0,1,2,... ) */
/*       Output:  SI(n) --- in(x) */
/*                DI(n) --- in'(x) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       ======================================================== */

    *nm = *n;
    if (abs(*x) < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    si[k] = 0.;
/* L10: */
	    di[k] = 0.;
	}
	si[0] = 1.;
	di[1] = .333333333333333;
	return 0;
    }
    si[0] = sinh(*x) / *x;
    si[1] = -(sinh(*x) / *x - cosh(*x)) / *x;
    si0 = si[0];
    if (*n >= 2) {
	m = msta1_(x, &c__200);
	if (m < *n) {
	    *nm = m;
	} else {
	    m = msta2_(x, n, &c__15);
	}
	f = 0.;
	f0 = 0.;
	f1 = -99.;
	for (k = m; k >= 0; --k) {
	    f = (k * 2. + 3.) * f1 / *x + f0;
	    if (k <= *nm) {
		si[k] = f;
	    }
	    f0 = f1;
/* L15: */
	    f1 = f;
	}
	cs = si0 / f;
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
/* L20: */
	    si[k] = cs * si[k];
	}
    }
    di[0] = si[1];
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	di[k] = si[k - 1] - (k + 1.) / *x * si[k];
    }
    return 0;
} /* sphi_ */

/*       ********************************** */
/* Subroutine */ int pbwa_(doublereal *a, doublereal *x, doublereal *w1f, 
	doublereal *w1d, doublereal *w2f, doublereal *w2d)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__[80], h__[100];
    static integer k, m;
    static doublereal r__, d1, d2, f1, g1, g2, f2, h0, h1;
    static integer l1, l2;
    static doublereal p0, r1, x1, y1, x2, dl, hl, y1f, y1d, y2f, y2d, ugi, 
	    vgi, eps, ugr, vgr;
    extern /* Subroutine */ int cgama_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute parabolic cylinder functions W(a,±x) */
/*                and their derivatives */
/*       Input  : a --- Parameter  ( 0 ≤ |a| ≤ 5 ) */
/*                x --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 ) */
/*       Output : W1F --- W(a,x) */
/*                W1D --- W'(a,x) */
/*                W2F --- W(a,-x) */
/*                W2D --- W'(a,-x) */
/*       Routine called: */
/*               CGAMA for computing complex gamma function */
/*       ====================================================== */

    eps = 1e-15;
    p0 = .59460355750136;
    if (*a == 0.) {
	g1 = 3.625609908222;
	g2 = 1.225416702465;
    } else {
	x1 = .25;
	y1 = *a * .5;
	cgama_(&x1, &y1, &c__1, &ugr, &ugi);
	g1 = sqrt(ugr * ugr + ugi * ugi);
	x2 = .75;
	cgama_(&x2, &y1, &c__1, &vgr, &vgi);
	g2 = sqrt(vgr * vgr + vgi * vgi);
    }
    f1 = sqrt(g1 / g2);
    f2 = sqrt(g2 * 2. / g1);
    h0 = 1.;
    h1 = *a;
    h__[0] = *a;
    for (l1 = 4; l1 <= 200; l1 += 2) {
	m = l1 / 2;
	hl = *a * h1 - (l1 - 2.) * .25 * (l1 - 3.) * h0;
	h__[m - 1] = hl;
	h0 = h1;
/* L10: */
	h1 = hl;
    }
    y1f = 1.;
    r__ = 1.;
    for (k = 1; k <= 100; ++k) {
	r__ = r__ * .5 * *x * *x / (k * (k * 2. - 1.));
	r1 = h__[k - 1] * r__;
	y1f += r1;
	if (abs(r1) <= eps * abs(y1f) && k > 30) {
	    goto L20;
	}
/* L15: */
    }
L20:
    y1d = *a;
    r__ = 1.;
    for (k = 1; k <= 99; ++k) {
	r__ = r__ * .5 * *x * *x / (k * (k * 2. + 1.));
	r1 = h__[k] * r__;
	y1d += r1;
	if (abs(r1) <= eps * abs(y1d) && k > 30) {
	    goto L30;
	}
/* L25: */
    }
L30:
    y1d = *x * y1d;
    d1 = 1.;
    d2 = *a;
    d__[0] = 1.;
    d__[1] = *a;
    for (l2 = 5; l2 <= 160; l2 += 2) {
	m = (l2 + 1) / 2;
	dl = *a * d2 - (l2 - 2.) * .25 * (l2 - 3.) * d1;
	d__[m - 1] = dl;
	d1 = d2;
/* L40: */
	d2 = dl;
    }
    y2f = 1.;
    r__ = 1.;
    for (k = 1; k <= 79; ++k) {
	r__ = r__ * .5 * *x * *x / (k * (k * 2. + 1.));
	r1 = d__[k] * r__;
	y2f += r1;
	if (abs(r1) <= eps * abs(y2f) && k > 30) {
	    goto L50;
	}
/* L45: */
    }
L50:
    y2f = *x * y2f;
    y2d = 1.;
    r__ = 1.;
    for (k = 1; k <= 79; ++k) {
	r__ = r__ * .5 * *x * *x / (k * (k * 2. - 1.));
	r1 = d__[k] * r__;
	y2d += r1;
	if (abs(r1) <= eps * abs(y2f) && k > 30) {
	    goto L60;
	}
/* L55: */
    }
L60:
    *w1f = p0 * (f1 * y1f - f2 * y2f);
    *w2f = p0 * (f1 * y1f + f2 * y2f);
    *w1d = p0 * (f1 * y1d - f2 * y2d);
    *w2d = p0 * (f1 * y1d + f2 * y2d);
    return 0;
} /* pbwa_ */

/*       ********************************** */
/* Subroutine */ int rmn1_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, doublereal *df, integer *kd, doublereal *r1f, 
	doublereal *r1d)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_ri(real *, integer *), pow_di(doublereal *, integer *), pow_dd(
	    doublereal *, doublereal *);

    /* Local variables */
    static integer j, k, l;
    static doublereal r__, a0, b0, r0, r1, r2, r3, ck[200], dj[252];
    static integer lg, ip, nm;
    static doublereal cx, sj[252];
    static integer np;
    static doublereal sw, sa0;
    static integer nm1, nm2;
    static doublereal sw1, reg, eps, suc, sud, sum;
    extern /* Subroutine */ int sckb_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), sphj_(integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);


/*       ======================================================= */
/*       Purpose: Compute prolate and oblate spheroidal radial */
/*                functions of the first kind for given m, n, */
/*                c and x */
/*       Routines called: */
/*            (1) SCKB for computing expansion coefficients c2k */
/*            (2) SPHJ for computing the spherical Bessel */
/*                functions of the first kind */
/*       ======================================================= */

    /* Parameter adjustments */
    --df;

    /* Function Body */
    eps = 1e-14;
    ip = 1;
    nm1 = (*n - *m) / 2;
    if (*n - *m == nm1 << 1) {
	ip = 0;
    }
    nm = nm1 + 25 + (integer) (*c__);
    reg = 1.;
    if (*m + nm > 80) {
	reg = 1e-200;
    }
    r0 = reg;
    i__1 = (*m << 1) + ip;
    for (j = 1; j <= i__1; ++j) {
/* L10: */
	r0 *= j;
    }
    r__ = r0;
    suc = r__ * df[1];
    sw = 0.;
    i__1 = nm;
    for (k = 2; k <= i__1; ++k) {
	r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k + ip 
		- 1.5);
	suc += r__ * df[k];
	if (k > nm1 && (d__1 = suc - sw, abs(d__1)) < abs(suc) * eps) {
	    goto L20;
	}
/* L15: */
	sw = suc;
    }
L20:
    if (*x == 0.f) {
	sckb_(m, n, c__, &df[1], ck);
	sum = 0.;
	sw1 = 0.;
	i__1 = nm;
	for (j = 1; j <= i__1; ++j) {
	    sum += ck[j - 1];
	    if ((d__1 = sum - sw1, abs(d__1)) < abs(sum) * eps) {
		goto L30;
	    }
/* L25: */
	    sw1 = sum;
	}
L30:
	r1 = 1.;
	i__1 = (*n + *m + ip) / 2;
	for (j = 1; j <= i__1; ++j) {
/* L35: */
	    r1 *= j + (*n + *m + ip) * .5;
	}
	r2 = 1.;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
/* L40: */
	    r2 = *c__ * 2. * r2 * j;
	}
	r3 = 1.;
	i__1 = (*n - *m - ip) / 2;
	for (j = 1; j <= i__1; ++j) {
/* L45: */
	    r3 *= j;
	}
	sa0 = ((*m + ip) * 2.f + 1.f) * r1 / (pow_ri(&c_b317, n) * pow_di(c__,
		 &ip) * r2 * r3);
	if (ip == 0) {
	    *r1f = sum / (sa0 * suc) * df[1] * reg;
	    *r1d = 0.;
	} else if (ip == 1) {
	    *r1f = 0.;
	    *r1d = sum / (sa0 * suc) * df[1] * reg;
	}
	return 0;
    }
    cx = *c__ * *x;
    nm2 = (nm << 1) + *m;
    sphj_(&nm2, &cx, &nm2, sj, dj);
    d__1 = 1. - *kd / (*x * *x);
    d__2 = *m * .5;
    a0 = pow_dd(&d__1, &d__2) / suc;
    *r1f = 0.;
    sw = 0.;
    lg = 0;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	l = (k << 1) + *m - *n - 2 + ip;
	if (l == l / 4 << 2) {
	    lg = 1;
	}
	if (l != l / 4 << 2) {
	    lg = -1;
	}
	if (k == 1) {
	    r__ = r0;
	} else {
	    r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k 
		    + ip - 1.5);
	}
	np = *m + (k << 1) - 2 + ip;
	*r1f += lg * r__ * df[k] * sj[np];
	if (k > nm1 && (d__1 = *r1f - sw, abs(d__1)) < abs(*r1f) * eps) {
	    goto L55;
	}
/* L50: */
	sw = *r1f;
    }
L55:
    *r1f *= a0;
    b0 = *kd * *m / pow_dd(x, &c_b149) / (1.f - *kd / (*x * *x)) * *r1f;
    sud = 0.;
    sw = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	l = (k << 1) + *m - *n - 2 + ip;
	if (l == l / 4 << 2) {
	    lg = 1;
	}
	if (l != l / 4 << 2) {
	    lg = -1;
	}
	if (k == 1) {
	    r__ = r0;
	} else {
	    r__ = r__ * (*m + k - 1.f) * (*m + k + ip - 1.5) / (k - 1.) / (k 
		    + ip - 1.5);
	}
	np = *m + (k << 1) - 2 + ip;
	sud += lg * r__ * df[k] * dj[np];
	if (k > nm1 && (d__1 = sud - sw, abs(d__1)) < abs(sud) * eps) {
	    goto L65;
	}
/* L60: */
	sw = sud;
    }
L65:
    *r1d = b0 + a0 * *c__ * sud;
    return 0;
} /* rmn1_ */

/*       ********************************** */
/* Subroutine */ int dvsa_(doublereal *va, doublereal *x, doublereal *pd)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static integer m;
    static doublereal r__, a0, g0, g1, r1, ep, gm, pi, vm, vt, ga0, va0, sq2, 
	    eps;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute parabolic cylinder function Dv(x) */
/*                for small argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PD --- Dv(x) */
/*       Routine called: GAMMA2 for computing Г(x) */
/*       =================================================== */

    eps = 1e-15;
    pi = 3.141592653589793;
    sq2 = sqrt(2.);
    ep = exp(*x * -.25 * *x);
    va0 = (1. - *va) * .5;
    if (*va == 0.f) {
	*pd = ep;
    } else {
	if (*x == 0.f) {
	    if (va0 <= 0.f && va0 == (doublereal) ((integer) va0)) {
		*pd = 0.;
	    } else {
		gamma2_(&va0, &ga0);
		d__1 = *va * -.5;
		*pd = sqrt(pi) / (pow_dd(&c_b4, &d__1) * ga0);
	    }
	} else {
	    d__1 = -(*va);
	    gamma2_(&d__1, &g1);
	    d__1 = *va * -.5 - 1.;
	    a0 = pow_dd(&c_b4, &d__1) * ep / g1;
	    vt = *va * -.5;
	    gamma2_(&vt, &g0);
	    *pd = g0;
	    r__ = 1.;
	    for (m = 1; m <= 250; ++m) {
		vm = (m - *va) * .5;
		gamma2_(&vm, &gm);
		r__ = -r__ * sq2 * *x / m;
		r1 = gm * r__;
		*pd += r1;
		if (abs(r1) < abs(*pd) * eps) {
		    goto L15;
		}
/* L10: */
	    }
L15:
	    *pd = a0 * *pd;
	}
    }
    return 0;
} /* dvsa_ */

/*       ********************************** */
/* Subroutine */ int e1z_(doublecomplex *z__, doublecomplex *ce1)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);
    void z_log(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *), z_exp(doublecomplex *, 
	    doublecomplex *);

    /* Local variables */
    static integer k;
    static doublereal x, a0, el;
    static doublecomplex cr;
    static doublereal pi;
    static doublecomplex zc, zd;
    static doublereal xt;
    static doublecomplex zdc;


/*       ==================================================== */
/*       Purpose: Compute complex exponential integral E1(z) */
/*       Input :  z   --- Argument of E1(z) */
/*       Output:  CE1 --- E1(z) */
/*       ==================================================== */

    pi = 3.141592653589793;
    el = .5772156649015328;
    x = z__->r;
    a0 = z_abs(z__);
/*       Continued fraction converges slowly near negative real axis, */
/*       so use power series in a wedge around it until radius 40.0 */
    xt = (d__1 = d_imag(z__), abs(d__1)) * -2;
    if (a0 == 0.) {
	ce1->r = 1e300, ce1->i = 0.;
    } else if (a0 <= 5.f || x < xt && a0 < 40.f) {
/*          Power series */
	ce1->r = 1., ce1->i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 500; ++k) {
	    z__4.r = -cr.r, z__4.i = -cr.i;
	    d__1 = (doublereal) k;
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * 
		    z__->i + z__3.i * z__->r;
/* Computing 2nd power */
	    d__3 = k + 1.;
	    d__2 = d__3 * d__3;
	    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = ce1->r + cr.r, z__1.i = ce1->i + cr.i;
	    ce1->r = z__1.r, ce1->i = z__1.i;
	    if (z_abs(&cr) <= z_abs(ce1) * 1e-15) {
		goto L15;
	    }
/* L10: */
	}
L15:
	if (x <= 0.f && d_imag(z__) == 0.f) {
/*             Careful on the branch cut -- avoid signed zeros */
	    d__1 = -el;
	    z__5.r = -z__->r, z__5.i = -z__->i;
	    z_log(&z__4, &z__5);
	    z__3.r = d__1 - z__4.r, z__3.i = -z__4.i;
	    z__6.r = z__->r * ce1->r - z__->i * ce1->i, z__6.i = z__->r * 
		    ce1->i + z__->i * ce1->r;
	    z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
	    z__7.r = pi * 0., z__7.i = pi * 1.;
	    z__1.r = z__2.r - z__7.r, z__1.i = z__2.i - z__7.i;
	    ce1->r = z__1.r, ce1->i = z__1.i;
	} else {
	    d__1 = -el;
	    z_log(&z__3, z__);
	    z__2.r = d__1 - z__3.r, z__2.i = -z__3.i;
	    z__4.r = z__->r * ce1->r - z__->i * ce1->i, z__4.i = z__->r * 
		    ce1->i + z__->i * ce1->r;
	    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
	    ce1->r = z__1.r, ce1->i = z__1.i;
	}
    } else {
/*          Continued fraction https://dlmf.nist.gov/6.9 */

/*                           1     1     1     2     2     3     3 */
/*          E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ... */
/*                         Z +   1 +   Z +   1 +   Z +   1 +   Z + */
	zc.r = 0., zc.i = 0.;
	z_div(&z__1, &c_b6, z__);
	zd.r = z__1.r, zd.i = z__1.i;
	zdc.r = zd.r, zdc.i = zd.i;
	z__1.r = zc.r + zdc.r, z__1.i = zc.i + zdc.i;
	zc.r = z__1.r, zc.i = z__1.i;
	for (k = 1; k <= 500; ++k) {
	    d__1 = (doublereal) k;
	    z__3.r = d__1 * zd.r, z__3.i = d__1 * zd.i;
	    z__2.r = z__3.r + 1, z__2.i = z__3.i;
	    z_div(&z__1, &c_b6, &z__2);
	    zd.r = z__1.r, zd.i = z__1.i;
	    z__2.r = zd.r - 1, z__2.i = zd.i;
	    z__1.r = z__2.r * zdc.r - z__2.i * zdc.i, z__1.i = z__2.r * zdc.i 
		    + z__2.i * zdc.r;
	    zdc.r = z__1.r, zdc.i = z__1.i;
	    z__1.r = zc.r + zdc.r, z__1.i = zc.i + zdc.i;
	    zc.r = z__1.r, zc.i = z__1.i;
	    d__1 = (doublereal) k;
	    z__3.r = d__1 * zd.r, z__3.i = d__1 * zd.i;
	    z__2.r = z__3.r + z__->r, z__2.i = z__3.i + z__->i;
	    z_div(&z__1, &c_b6, &z__2);
	    zd.r = z__1.r, zd.i = z__1.i;
	    z__3.r = z__->r * zd.r - z__->i * zd.i, z__3.i = z__->r * zd.i + 
		    z__->i * zd.r;
	    z__2.r = z__3.r - 1, z__2.i = z__3.i;
	    z__1.r = z__2.r * zdc.r - z__2.i * zdc.i, z__1.i = z__2.r * zdc.i 
		    + z__2.i * zdc.r;
	    zdc.r = z__1.r, zdc.i = z__1.i;
	    z__1.r = zc.r + zdc.r, z__1.i = zc.i + zdc.i;
	    zc.r = z__1.r, zc.i = z__1.i;
	    if (z_abs(&zdc) <= z_abs(&zc) * 1e-15 && k > 20) {
		goto L25;
	    }
/* L20: */
	}
L25:
	z__3.r = -z__->r, z__3.i = -z__->i;
	z_exp(&z__2, &z__3);
	z__1.r = z__2.r * zc.r - z__2.i * zc.i, z__1.i = z__2.r * zc.i + 
		z__2.i * zc.r;
	ce1->r = z__1.r, ce1->i = z__1.i;
	if (x <= 0.f && d_imag(z__) == 0.f) {
	    z__2.r = pi * 0., z__2.i = pi * 1.;
	    z__1.r = ce1->r - z__2.r, z__1.i = ce1->i - z__2.i;
	    ce1->r = z__1.r, ce1->i = z__1.i;
	}
    }
    return 0;
} /* e1z_ */

/*       ********************************** */
/* Subroutine */ int itjyb_(doublereal *x, doublereal *tj, doublereal *ty)
{
    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal), sqrt(doublereal)
	    ;

    /* Local variables */
    static doublereal t, f0, g0, x1, pi, xt;


/*       ======================================================= */
/*       Purpose: Integrate Bessel functions J0(t) and Y0(t) */
/*                with respect to t from 0 to x ( x ≥ 0 ) */
/*       Input :  x  --- Upper limit of the integral */
/*       Output:  TJ --- Integration of J0(t) from 0 to x */
/*                TY --- Integration of Y0(t) from 0 to x */
/*       ======================================================= */

    pi = 3.141592653589793;
    if (*x == 0.) {
	*tj = 0.;
	*ty = 0.;
    } else if (*x <= 4.) {
	x1 = *x / 4.;
	t = x1 * x1;
	*tj = (((((((t * -1.33718e-4 + .002362211) * t - .025791036) * t + 
		.197492634) * t - 1.015860606) * t + 3.199997842) * t - 
		5.333333161) * t + 4.) * x1;
	*ty = ((((((((t * 1.3351e-5 - 2.35002e-4) * t + .003034322) * t - 
		.029600855) * t + .203380298) * t - .904755062) * t + 
		2.287317974) * t - 2.567250468) * t + 1.076611469) * x1;
	*ty = 2. / pi * log(*x / 2.) * *tj - *ty;
    } else if (*x <= 8.) {
	xt = *x - pi * .25;
	t = 16. / (*x * *x);
	f0 = ((((((t * .001496119 - .00739083) * t + .016236617) * t - 
		.022007499) * t + .023644978) * t - .031280848) * t + 
		.124611058) * 4. / *x;
	g0 = (((((t * .001076103 - .005434851) * t + .01242264) * t - 
		.018255209f) * t + .023664841) * t - .049635633) * t + 
		.79784879;
	*tj = 1. - (f0 * cos(xt) - g0 * sin(xt)) / sqrt(*x);
	*ty = -(f0 * sin(xt) + g0 * cos(xt)) / sqrt(*x);
    } else {
	t = 64. / (*x * *x);
	xt = *x - pi * .25;
	f0 = (((((((t * -2.68482e-5 + 1.270039e-4) * t - 2.755037e-4) * t + 
		3.992825e-4) * t - 5.366169e-4) * t + .0010089872) * t - 
		.0040403539) * t + .0623347304) * 8. / *x;
	g0 = ((((((t * -2.26238e-5 + 1.107299e-4) * t - 2.543955e-4) * t + 
		4.100676e-4) * t - 6.740148e-4) * t + .0017870944) * t - 
		.01256424405) * t + .79788456;
	*tj = 1. - (f0 * cos(xt) - g0 * sin(xt)) / sqrt(*x);
	*ty = -(f0 * sin(xt) + g0 * cos(xt)) / sqrt(*x);
    }
    return 0;
} /* itjyb_ */

/*       ********************************** */
/* Subroutine */ int chgul_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *hu, integer *id)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), d_lg10(doublereal *);

    /* Local variables */
    static integer k;
    static doublereal r__, r0, aa, ra;
    static integer nm;
    static logical il1, il2;


/*       ======================================================= */
/*       Purpose: Compute the confluent hypergeometric function */
/*                U(a,b,x) for large argument x */
/*       Input  : a  --- Parameter */
/*                b  --- Parameter */
/*                x  --- Argument */
/*       Output:  HU --- U(a,b,x) */
/*                ID --- Estimated number of significant digits */
/*       ======================================================= */

    *id = -100;
    aa = *a - *b + 1.;
    il1 = *a == (doublereal) ((integer) (*a)) && *a <= 0.f;
    il2 = aa == (doublereal) ((integer) aa) && aa <= 0.f;
    nm = 0;
    if (il1) {
	nm = (integer) abs(*a);
    }
    if (il2) {
	nm = (integer) abs(aa);
    }
/*       IL1: DLMF 13.2.7 with k=-s-a */
/*       IL2: DLMF 13.2.8 */
    if (il1 || il2) {
	*hu = 1.;
	r__ = 1.;
	i__1 = nm;
	for (k = 1; k <= i__1; ++k) {
	    r__ = -r__ * (*a + k - 1.) * (*a - *b + k) / (k * *x);
	    *hu += r__;
/* L10: */
	}
	d__1 = -(*a);
	*hu = pow_dd(x, &d__1) * *hu;
	*id = 10;
    } else {
/*       DLMF 13.7.3 */
	*hu = 1.;
	r__ = 1.;
	for (k = 1; k <= 25; ++k) {
	    r__ = -r__ * (*a + k - 1.) * (*a - *b + k) / (k * *x);
	    ra = abs(r__);
	    if (k > 5 && ra >= r0 || ra < 1e-15) {
		goto L20;
	    }
	    r0 = ra;
/* L15: */
	    *hu += r__;
	}
L20:
	*id = (d__1 = d_lg10(&ra), (integer) abs(d__1));
	d__1 = -(*a);
	*hu = pow_dd(x, &d__1) * *hu;
    }
    return 0;
} /* chgul_ */

/*       ********************************** */
/* Subroutine */ int gmn_(integer *m, integer *n, doublereal *c__, doublereal 
	*x, doublereal *bk, doublereal *gf, doublereal *gd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;

    /* Local variables */
    static integer k, ip, nm;
    static doublereal gw, xm, gd0, gd1, gf0, eps;


/*       =========================================================== */
/*       Purpose: Compute gmn(-ic,ix) and its derivative for oblate */
/*                radial functions with a small argument */
/*       =========================================================== */

    /* Parameter adjustments */
    --bk;

    /* Function Body */
    eps = 1e-14;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    nm = (integer) ((*n - *m) * .5f + *c__) + 25;
    d__1 = *x * *x + 1.;
    d__2 = *m * -.5;
    xm = pow_dd(&d__1, &d__2);
    gf0 = 0.;
    gw = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	d__1 = (doublereal) (k * 2.f - 2.f);
	gf0 += bk[k] * pow_dd(x, &d__1);
	if ((d__1 = (gf0 - gw) / gf0, abs(d__1)) < eps && k >= 10) {
	    goto L15;
	}
/* L10: */
	gw = gf0;
    }
L15:
    i__1 = 1 - ip;
    *gf = xm * gf0 * pow_di(x, &i__1);
    gd1 = -(*m) * *x / (*x * *x + 1.) * *gf;
    gd0 = 0.;
    i__1 = nm;
    for (k = 1; k <= i__1; ++k) {
	if (ip == 0) {
	    d__1 = (doublereal) (k * 2.f - 2.f);
	    gd0 += (k * 2. - 1.f) * bk[k] * pow_dd(x, &d__1);
	} else {
	    d__1 = (doublereal) (k * 2.f - 1.f);
	    gd0 += k * 2. * bk[k + 1] * pow_dd(x, &d__1);
	}
	if ((d__1 = (gd0 - gw) / gd0, abs(d__1)) < eps && k >= 10) {
	    goto L25;
	}
/* L20: */
	gw = gd0;
    }
L25:
    *gd = gd1 + xm * gd0;
    return 0;
} /* gmn_ */

/*       ********************************** */
/* Subroutine */ int itjya_(doublereal *x, doublereal *tj, doublereal *ty)
{
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static doublereal a[18];
    static integer k;
    static doublereal r__, a0, a1, r2, x2, af, bf, bg, el, rc, pi, rs, xp, 
	    ty1, ty2, eps;


/*       ========================================================== */
/*       Purpose: Integrate Bessel functions J0(t) & Y0(t) with */
/*                respect to t from 0 to x */
/*       Input :  x  --- Upper limit of the integral ( x >= 0 ) */
/*       Output:  TJ --- Integration of J0(t) from 0 to x */
/*                TY --- Integration of Y0(t) from 0 to x */
/*       ======================================================= */

    pi = 3.141592653589793;
    el = .5772156649015329;
    eps = 1e-12;
    if (*x == 0.) {
	*tj = 0.;
	*ty = 0.;
    } else if (*x <= 20.) {
	x2 = *x * *x;
	*tj = *x;
	r__ = *x;
	for (k = 1; k <= 60; ++k) {
	    r__ = r__ * -.25 * ((k << 1) - 1.) / ((k << 1) + 1.) / (k * k) * 
		    x2;
	    *tj += r__;
	    if (abs(r__) < abs(*tj) * eps) {
		goto L15;
	    }
/* L10: */
	}
L15:
	ty1 = (el + log(*x / 2.)) * *tj;
	rs = 0.;
	ty2 = 1.;
	r__ = 1.;
	for (k = 1; k <= 60; ++k) {
	    r__ = r__ * -.25 * ((k << 1) - 1.) / ((k << 1) + 1.) / (k * k) * 
		    x2;
	    rs += 1. / k;
	    r2 = r__ * (rs + 1. / (k * 2. + 1.));
	    ty2 += r2;
	    if (abs(r2) < abs(ty2) * eps) {
		goto L25;
	    }
/* L20: */
	}
L25:
	*ty = (ty1 - *x * ty2) * 2. / pi;
    } else {
	a0 = 1.;
	a1 = .625;
	a[0] = a1;
	for (k = 1; k <= 16; ++k) {
	    af = ((k + .5) * 1.5 * (k + .83333333333333337) * a1 - (k + .5) * 
		    .5 * (k + .5) * (k - .5) * a0) / (k + 1.);
	    a[k] = af;
	    a0 = a1;
/* L30: */
	    a1 = af;
	}
	bf = 1.;
	r__ = 1.;
	for (k = 1; k <= 8; ++k) {
	    r__ = -r__ / (*x * *x);
/* L35: */
	    bf += a[(k << 1) - 1] * r__;
	}
	bg = a[0] / *x;
	r__ = 1. / *x;
	for (k = 1; k <= 8; ++k) {
	    r__ = -r__ / (*x * *x);
/* L40: */
	    bg += a[k * 2] * r__;
	}
	xp = *x + pi * .25;
	rc = sqrt(2. / (pi * *x));
	*tj = 1. - rc * (bf * cos(xp) + bg * sin(xp));
	*ty = rc * (bg * cos(xp) - bf * sin(xp));
    }
    return 0;
} /* itjya_ */

/*       ********************************** */
/* Subroutine */ int rcty_(integer *n, doublereal *x, integer *nm, doublereal 
	*ry, doublereal *dy)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal rf0, rf1, rf2;


/*       ======================================================== */
/*       Purpose: Compute Riccati-Bessel functions of the second */
/*                kind and their derivatives */
/*       Input:   x --- Argument of Riccati-Bessel function */
/*                n --- Order of yn(x) */
/*       Output:  RY(n) --- x·yn(x) */
/*                DY(n) --- [x·yn(x)]' */
/*                NM --- Highest order computed */
/*       ======================================================== */

    *nm = *n;
    if (*x < 1e-60) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    ry[k] = -1e300;
/* L10: */
	    dy[k] = 1e300;
	}
	ry[0] = -1.;
	dy[0] = 0.;
	return 0;
    }
    ry[0] = -cos(*x);
    ry[1] = ry[0] / *x - sin(*x);
    rf0 = ry[0];
    rf1 = ry[1];
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	rf2 = (k * 2. - 1.) * rf1 / *x - rf0;
	if (abs(rf2) > 1e300) {
	    goto L20;
	}
	ry[k] = rf2;
	rf0 = rf1;
/* L15: */
	rf1 = rf2;
    }
L20:
    *nm = k - 1;
    dy[0] = sin(*x);
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	dy[k] = -k * ry[k] / *x + ry[k - 1];
    }
    return 0;
} /* rcty_ */

/*       ********************************** */
/* Subroutine */ int lpni_(integer *n, doublereal *x, doublereal *pn, 
	doublereal *pd, doublereal *pl)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer j, k;
    static doublereal r__;
    static integer n1;
    static doublereal p0, p1, pf;


/*       ===================================================== */
/*       Purpose: Compute Legendre polynomials Pn(x), Pn'(x) */
/*                and the integral of Pn(t) from 0 to x */
/*       Input :  x --- Argument of Pn(x) */
/*                n --- Degree of Pn(x) ( n = 0,1,... ) */
/*       Output:  PN(n) --- Pn(x) */
/*                PD(n) --- Pn'(x) */
/*                PL(n) --- Integral of Pn(t) from 0 to x */
/*       ===================================================== */

    pn[0] = 1.;
    pn[1] = *x;
    pd[0] = 0.;
    pd[1] = 1.;
    pl[0] = *x;
    pl[1] = *x * .5 * *x;
    p0 = 1.;
    p1 = *x;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	pf = (k * 2. - 1.) / k * *x * p1 - (k - 1.) / k * p0;
	pn[k] = pf;
	if (abs(*x) == 1.) {
	    i__2 = k + 1;
	    pd[k] = pow_di(x, &i__2) * .5 * k * (k + 1.);
	} else {
	    pd[k] = k * (p1 - *x * pf) / (1. - *x * *x);
	}
	pl[k] = (*x * pn[k] - pn[k - 1]) / (k + 1.);
	p0 = p1;
	p1 = pf;
	if (k == k / 2 << 1) {
	    goto L15;
	}
	r__ = 1. / (k + 1.);
	n1 = (k - 1) / 2;
	i__2 = n1;
	for (j = 1; j <= i__2; ++j) {
/* L10: */
	    r__ = (.5 / j - 1.) * r__;
	}
	pl[k] += r__;
L15:
	;
    }
    return 0;
} /* lpni_ */

/*       ********************************** */
/* Subroutine */ int klvna_(doublereal *x, doublereal *ber, doublereal *bei, 
	doublereal *ger, doublereal *gei, doublereal *der, doublereal *dei, 
	doublereal *her, doublereal *hei)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), cos(doublereal), sin(doublereal), sqrt(doublereal)
	    , exp(doublereal);

    /* Local variables */
    static integer k, m;
    static doublereal r__, r0, r1, x2, x4, el, rc, cs;
    static integer km;
    static doublereal pi, gs, xd, rs, ss, xt, cn0, cp0, xc1, xc2, pn0, qn0, 
	    pp0, qp0, xe1, xe2, sp0, sn0, pp1, pn1, qp1, qn1, fac, eps;


/*       ====================================================== */
/*       Purpose: Compute Kelvin functions ber x, bei x, ker x */
/*                and kei x, and their derivatives  ( x > 0 ) */
/*       Input :  x   --- Argument of Kelvin functions */
/*       Output:  BER --- ber x */
/*                BEI --- bei x */
/*                GER --- ker x */
/*                GEI --- kei x */
/*                DER --- ber'x */
/*                DEI --- bei'x */
/*                HER --- ker'x */
/*                HEI --- kei'x */
/*       ================================================ */

    pi = 3.141592653589793;
    el = .5772156649015329;
    eps = 1e-15;
    if (*x == 0.) {
	*ber = 1.;
	*bei = 0.;
	*ger = 1e300;
	*gei = pi * -.25;
	*der = 0.;
	*dei = 0.;
	*her = -1e300;
	*hei = 0.;
	return 0;
    }
    x2 = *x * .25 * *x;
    x4 = x2 * x2;
    if (abs(*x) < 10.) {
	*ber = 1.;
	r__ = 1.;
	for (m = 1; m <= 60; ++m) {
/* Computing 2nd power */
	    d__1 = m * 2. - 1.;
	    r__ = r__ * -.25 / (m * m) / (d__1 * d__1) * x4;
	    *ber += r__;
	    if (abs(r__) < abs(*ber) * eps) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*bei = x2;
	r__ = x2;
	for (m = 1; m <= 60; ++m) {
/* Computing 2nd power */
	    d__1 = m * 2. + 1.;
	    r__ = r__ * -.25 / (m * m) / (d__1 * d__1) * x4;
	    *bei += r__;
	    if (abs(r__) < abs(*bei) * eps) {
		goto L25;
	    }
/* L20: */
	}
L25:
	*ger = -(log(*x / 2.) + el) * *ber + pi * .25 * *bei;
	r__ = 1.;
	gs = 0.;
	for (m = 1; m <= 60; ++m) {
/* Computing 2nd power */
	    d__1 = m * 2. - 1.;
	    r__ = r__ * -.25 / (m * m) / (d__1 * d__1) * x4;
	    gs = gs + 1. / (m * 2. - 1.) + 1. / (m * 2.);
	    *ger += r__ * gs;
	    if ((d__1 = r__ * gs, abs(d__1)) < abs(*ger) * eps) {
		goto L35;
	    }
/* L30: */
	}
L35:
	*gei = x2 - (log(*x / 2.) + el) * *bei - pi * .25 * *ber;
	r__ = x2;
	gs = 1.;
	for (m = 1; m <= 60; ++m) {
/* Computing 2nd power */
	    d__1 = m * 2. + 1.;
	    r__ = r__ * -.25 / (m * m) / (d__1 * d__1) * x4;
	    gs = gs + 1. / (m * 2.) + 1. / (m * 2. + 1.);
	    *gei += r__ * gs;
	    if ((d__1 = r__ * gs, abs(d__1)) < abs(*gei) * eps) {
		goto L45;
	    }
/* L40: */
	}
L45:
	*der = *x * -.25 * x2;
	r__ = *der;
	for (m = 1; m <= 60; ++m) {
/* Computing 2nd power */
	    d__1 = m * 2. + 1.;
	    r__ = r__ * -.25 / m / (m + 1.) / (d__1 * d__1) * x4;
	    *der += r__;
	    if (abs(r__) < abs(*der) * eps) {
		goto L55;
	    }
/* L50: */
	}
L55:
	*dei = *x * .5;
	r__ = *dei;
	for (m = 1; m <= 60; ++m) {
	    r__ = r__ * -.25 / (m * m) / (m * 2. - 1.) / (m * 2. + 1.) * x4;
	    *dei += r__;
	    if (abs(r__) < abs(*dei) * eps) {
		goto L65;
	    }
/* L60: */
	}
L65:
	r__ = *x * -.25 * x2;
	gs = 1.5;
	*her = r__ * 1.5 - *ber / *x - (log(*x / 2.) + el) * *der + pi * .25f 
		* *dei;
	for (m = 1; m <= 60; ++m) {
/* Computing 2nd power */
	    d__1 = m * 2. + 1.;
	    r__ = r__ * -.25 / m / (m + 1.) / (d__1 * d__1) * x4;
	    gs = gs + 1. / ((m << 1) + 1.) + 1. / ((m << 1) + 2.);
	    *her += r__ * gs;
	    if ((d__1 = r__ * gs, abs(d__1)) < abs(*her) * eps) {
		goto L75;
	    }
/* L70: */
	}
L75:
	r__ = *x * .5;
	gs = 1.;
	*hei = *x * .5 - *bei / *x - (log(*x / 2.) + el) * *dei - pi * .25f * 
		*der;
	for (m = 1; m <= 60; ++m) {
	    r__ = r__ * -.25 / (m * m) / ((m << 1) - 1.) / ((m << 1) + 1.) * 
		    x4;
	    gs = gs + 1. / (m * 2.) + 1. / ((m << 1) + 1.);
	    *hei += r__ * gs;
	    if ((d__1 = r__ * gs, abs(d__1)) < abs(*hei) * eps) {
		return 0;
	    }
/* L80: */
	}
    } else {
	pp0 = 1.;
	pn0 = 1.;
	qp0 = 0.;
	qn0 = 0.;
	r0 = 1.;
	km = 18;
	if (abs(*x) >= 40.f) {
	    km = 10;
	}
	fac = 1.;
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
	    fac = -fac;
	    xt = k * .25 * pi - (integer) (k * .125) * 2. * pi;
	    cs = cos(xt);
	    ss = sin(xt);
/* Computing 2nd power */
	    d__1 = k * 2. - 1.;
	    r0 = r0 * .125 * (d__1 * d__1) / k / *x;
	    rc = r0 * cs;
	    rs = r0 * ss;
	    pp0 += rc;
	    pn0 += fac * rc;
	    qp0 += rs;
/* L85: */
	    qn0 += fac * rs;
	}
	xd = *x / sqrt(2.);
	xe1 = exp(xd);
	xe2 = exp(-xd);
	xc1 = 1. / sqrt(pi * 2. * *x);
	xc2 = sqrt(pi * .5 / *x);
	cp0 = cos(xd + pi * .125);
	cn0 = cos(xd - pi * .125);
	sp0 = sin(xd + pi * .125);
	sn0 = sin(xd - pi * .125);
	*ger = xc2 * xe2 * (pn0 * cp0 - qn0 * sp0);
	*gei = xc2 * xe2 * (-pn0 * sp0 - qn0 * cp0);
	*ber = xc1 * xe1 * (pp0 * cn0 + qp0 * sn0) - *gei / pi;
	*bei = xc1 * xe1 * (pp0 * sn0 - qp0 * cn0) + *ger / pi;
	pp1 = 1.;
	pn1 = 1.;
	qp1 = 0.;
	qn1 = 0.;
	r1 = 1.;
	fac = 1.;
	i__1 = km;
	for (k = 1; k <= i__1; ++k) {
	    fac = -fac;
	    xt = k * .25 * pi - (integer) (k * .125) * 2. * pi;
	    cs = cos(xt);
	    ss = sin(xt);
/* Computing 2nd power */
	    d__1 = k * 2. - 1.;
	    r1 = r1 * .125 * (4. - d__1 * d__1) / k / *x;
	    rc = r1 * cs;
	    rs = r1 * ss;
	    pp1 += fac * rc;
	    pn1 += rc;
	    qp1 += fac * rs;
	    qn1 += rs;
/* L90: */
	}
	*her = xc2 * xe2 * (-pn1 * cn0 + qn1 * sn0);
	*hei = xc2 * xe2 * (pn1 * sn0 + qn1 * cn0);
	*der = xc1 * xe1 * (pp1 * cp0 + qp1 * sp0) - *hei / pi;
	*dei = xc1 * xe1 * (pp1 * sp0 - qp1 * cp0) + *her / pi;
    }
    return 0;
} /* klvna_ */

/*       ********************************** */
/* Subroutine */ int chgubi_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *hu, integer *id)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double pow_di(doublereal *, integer *), d_lg10(doublereal *), log(
	    doublereal);

    /* Local variables */
    extern /* Subroutine */ int psi_spec__(doublereal *, doublereal *);
    static integer j, k, m, n;
    static doublereal r__, a0, a1, a2, h0, s0, s1, s2, ga, el, sa, sb, ua, ub,
	     hw, rn, ps, da1, da2, db1, ga1, db2;
    static integer id1, id2;
    static doublereal hm1, hm2, hm3, hu1, hu2, rn1, hmin, hmax;
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ====================================================== */
/*       Purpose: Compute confluent hypergeometric function */
/*                U(a,b,x) with integer b ( b = ±1,±2,... ) */
/*       Input  : a  --- Parameter */
/*                b  --- Parameter */
/*                x  --- Argument */
/*       Output:  HU --- U(a,b,x) */
/*                ID --- Estimated number of significant digits */
/*       Routines called: */
/*            (1) GAMMA2 for computing gamma function Г(x) */
/*            (2) PSI_SPEC for computing psi function */
/*       ====================================================== */

    *id = -100;
    el = .5772156649015329;
    n = (d__1 = *b - 1, (integer) abs(d__1));
    rn1 = 1.;
    rn = 1.;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	rn *= j;
	if (j == n - 1) {
	    rn1 = rn;
	}
/* L10: */
    }
    psi_spec__(a, &ps);
    gamma2_(a, &ga);
    if (*b > 0.f) {
	a0 = *a;
	a1 = *a - n;
	a2 = a1;
	gamma2_(&a1, &ga1);
	i__1 = n - 1;
	ua = pow_ii(&c_n1, &i__1) / (rn * ga1);
	i__1 = -n;
	ub = rn1 / ga * pow_di(x, &i__1);
    } else {
	a0 = *a + n;
	a1 = a0;
	a2 = *a;
	gamma2_(&a1, &ga1);
	i__1 = n - 1;
	ua = pow_ii(&c_n1, &i__1) / (rn * ga) * pow_di(x, &n);
	ub = rn1 / ga1;
    }
    hm1 = 1.;
    r__ = 1.;
    hmax = 0.;
    hmin = 1e300;
    h0 = 0.;
    for (k = 1; k <= 150; ++k) {
	r__ = r__ * (a0 + k - 1.) * *x / ((n + k) * k);
	hm1 += r__;
	hu1 = abs(hm1);
	if (hu1 > hmax) {
	    hmax = hu1;
	}
	if (hu1 < hmin) {
	    hmin = hu1;
	}
	if ((d__1 = hm1 - h0, abs(d__1)) < abs(hm1) * 1e-15) {
	    goto L20;
	}
/* L15: */
	h0 = hm1;
    }
L20:
    da1 = d_lg10(&hmax);
    da2 = 0.;
    if (hmin != 0.f) {
	da2 = d_lg10(&hmin);
    }
    *id = (integer) (15 - (d__1 = da1 - da2, abs(d__1)));
    hm1 *= log(*x);
    s0 = 0.;
    i__1 = n;
    for (m = 1; m <= i__1; ++m) {
	if (*b >= 0.f) {
	    s0 -= 1. / m;
	}
/* L25: */
	if (*b < 0.f) {
	    s0 += (1. - *a) / (m * (*a + m - 1.));
	}
    }
    hm2 = ps + el * 2. + s0;
    r__ = 1.;
    hmax = 0.;
    hmin = 1e300;
    for (k = 1; k <= 150; ++k) {
	s1 = 0.;
	s2 = 0.;
	if (*b > 0.f) {
	    i__1 = k;
	    for (m = 1; m <= i__1; ++m) {
/* L30: */
		s1 -= (m + *a * 2. - 2.) / (m * (m + *a - 1.));
	    }
	    i__1 = n;
	    for (m = 1; m <= i__1; ++m) {
/* L35: */
		s2 += 1. / (k + m);
	    }
	} else {
	    i__1 = k + n;
	    for (m = 1; m <= i__1; ++m) {
/* L40: */
		s1 += (1. - *a) / (m * (m + *a - 1.));
	    }
	    i__1 = k;
	    for (m = 1; m <= i__1; ++m) {
/* L45: */
		s2 += 1. / m;
	    }
	}
	hw = el * 2. + ps + s1 - s2;
	r__ = r__ * (a0 + k - 1.) * *x / ((n + k) * k);
	hm2 += r__ * hw;
	hu2 = abs(hm2);
	if (hu2 > hmax) {
	    hmax = hu2;
	}
	if (hu2 < hmin) {
	    hmin = hu2;
	}
	if ((d__1 = (hm2 - h0) / hm2, abs(d__1)) < 1e-15) {
	    goto L55;
	}
/* L50: */
	h0 = hm2;
    }
L55:
    db1 = d_lg10(&hmax);
    db2 = 0.;
    if (hmin != 0.f) {
	db2 = d_lg10(&hmin);
    }
    id1 = (integer) (15 - (d__1 = db1 - db2, abs(d__1)));
    if (id1 < *id) {
	*id = id1;
    }
    hm3 = 1.;
    if (n == 0) {
	hm3 = 0.;
    }
    r__ = 1.;
    i__1 = n - 1;
    for (k = 1; k <= i__1; ++k) {
	r__ = r__ * (a2 + k - 1.) / ((k - n) * k) * *x;
/* L60: */
	hm3 += r__;
    }
    sa = ua * (hm1 + hm2);
    sb = ub * hm3;
    *hu = sa + sb;
    id2 = 0;
    if (sa != 0.f) {
	d__1 = abs(sa);
	id1 = (integer) d_lg10(&d__1);
    }
    if (*hu != 0.f) {
	d__1 = abs(*hu);
	id2 = (integer) d_lg10(&d__1);
    }
    if (sa * sb < 0.f) {
	*id -= (i__1 = id1 - id2, abs(i__1));
    }
    return 0;
} /* chgubi_ */

/*       ********************************** */
/* Subroutine */ int cyzo_(integer *nt, integer *kf, integer *kc, 
	doublecomplex *zo, doublecomplex *zv)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static doublereal h__;
    static integer i__, j;
    static doublereal w, x, y;
    static doublecomplex z__;
    static doublereal w0;
    static integer it;
    static doublecomplex zd;
    static integer nr;
    static doublecomplex zf, zp, zq, zw;
    extern /* Subroutine */ int cy01_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    static doublecomplex zfd, zgd, zero;


/*       =========================================================== */
/*       Purpose : Compute the complex zeros of Y0(z), Y1(z) and */
/*                 Y1'(z), and their associated values at the zeros */
/*                 using the modified Newton's iteration method */
/*       Input:    NT --- Total number of zeros/roots */
/*                 KF --- Function choice code */
/*                        KF=0 for  Y0(z) & Y1(z0) */
/*                        KF=1 for  Y1(z) & Y0(z1) */
/*                        KF=2 for  Y1'(z) & Y1(z1') */
/*                 KC --- Choice code */
/*                        KC=0 for complex roots */
/*                        KC=1 for real roots */
/*       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z) */
/*                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z) */
/*                           at the L-th zero */
/*       Routine called: CY01 for computing Y0(z) and Y1(z), and */
/*                       their derivatives */
/*       =========================================================== */
    /* Parameter adjustments */
    --zv;
    --zo;

    /* Function Body */
    x = 0.;
    y = 0.;
    h__ = 0.;
    if (*kc == 0) {
	x = -2.4;
	y = .54;
	h__ = 3.14;
    } else if (*kc == 1) {
	x = .89f;
	y = 0.f;
	h__ = -3.14f;
    }
    if (*kf == 1) {
	x = -.503f;
    }
    if (*kf == 2) {
	x = .577f;
    }
    z__1.r = x, z__1.i = y;
    zero.r = z__1.r, zero.i = z__1.i;
    z__.r = zero.r, z__.i = zero.i;
    w = 0.;
    i__1 = *nt;
    for (nr = 1; nr <= i__1; ++nr) {
	if (nr != 1) {
	    i__2 = nr - 1;
	    z__1.r = zo[i__2].r - h__, z__1.i = zo[i__2].i;
	    z__.r = z__1.r, z__.i = z__1.i;
	}
	it = 0;
L15:
	++it;
	cy01_(kf, &z__, &zf, &zd);
	zp.r = 1., zp.i = 0.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    i__3 = i__;
	    z__2.r = z__.r - zo[i__3].r, z__2.i = z__.i - zo[i__3].i;
	    z__1.r = zp.r * z__2.r - zp.i * z__2.i, z__1.i = zp.r * z__2.i + 
		    zp.i * z__2.r;
	    zp.r = z__1.r, zp.i = z__1.i;
	}
	z_div(&z__1, &zf, &zp);
	zfd.r = z__1.r, zfd.i = z__1.i;
	zq.r = 0., zq.i = 0.;
	i__3 = nr - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    zw.r = 1., zw.i = 0.;
	    i__2 = nr - 1;
	    for (j = 1; j <= i__2; ++j) {
		if (j == i__) {
		    goto L25;
		}
		i__4 = j;
		z__2.r = z__.r - zo[i__4].r, z__2.i = z__.i - zo[i__4].i;
		z__1.r = zw.r * z__2.r - zw.i * z__2.i, z__1.i = zw.r * 
			z__2.i + zw.i * z__2.r;
		zw.r = z__1.r, zw.i = z__1.i;
L25:
		;
	    }
	    z__1.r = zq.r + zw.r, z__1.i = zq.i + zw.i;
	    zq.r = z__1.r, zq.i = z__1.i;
/* L30: */
	}
	z__3.r = zq.r * zfd.r - zq.i * zfd.i, z__3.i = zq.r * zfd.i + zq.i * 
		zfd.r;
	z__2.r = zd.r - z__3.r, z__2.i = zd.i - z__3.i;
	z_div(&z__1, &z__2, &zp);
	zgd.r = z__1.r, zgd.i = z__1.i;
	z_div(&z__2, &zfd, &zgd);
	z__1.r = z__.r - z__2.r, z__1.i = z__.i - z__2.i;
	z__.r = z__1.r, z__.i = z__1.i;
	w0 = w;
	w = z_abs(&z__);
	if (it <= 50 && (d__1 = (w - w0) / w, abs(d__1)) > 1e-12) {
	    goto L15;
	}
	i__3 = nr;
	zo[i__3].r = z__.r, zo[i__3].i = z__.i;
/* L35: */
    }
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = i__;
	z__.r = zo[i__3].r, z__.i = zo[i__3].i;
	if (*kf == 0 || *kf == 2) {
	    cy01_(&c__1, &z__, &zf, &zd);
	    i__3 = i__;
	    zv[i__3].r = zf.r, zv[i__3].i = zf.i;
	} else if (*kf == 1) {
	    cy01_(&c__0, &z__, &zf, &zd);
	    i__3 = i__;
	    zv[i__3].r = zf.r, zv[i__3].i = zf.i;
	}
/* L40: */
    }
    return 0;
} /* cyzo_ */

/*       ********************************** */
/* Subroutine */ int klvnb_(doublereal *x, doublereal *ber, doublereal *bei, 
	doublereal *ger, doublereal *gei, doublereal *der, doublereal *dei, 
	doublereal *her, doublereal *hei)
{
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal), exp(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer l;
    static doublereal t, u, v, t2, pi, yd, yc1, yc2, ye1, ye2, csn, csp, fxi, 
	    pni, ppi, tni, tpi, fxr, pnr, ppr, tnr, ssn, tpr, ssp;


/*       ====================================================== */
/*       Purpose: Compute Kelvin functions ber x, bei x, ker x */
/*                and kei x, and their derivatives  ( x > 0 ) */
/*       Input :  x   --- Argument of Kelvin functions */
/*       Output:  BER --- ber x */
/*                BEI --- bei x */
/*                GER --- ker x */
/*                GEI --- kei x */
/*                DER --- ber'x */
/*                DEI --- bei'x */
/*                HER --- ker'x */
/*                HEI --- kei'x */
/*       ================================================ */

    pi = 3.141592653589793;
    if (*x == 0.) {
	*ber = 1.;
	*bei = 0.;
	*ger = 1e300;
	*gei = pi * -.25;
	*der = 0.;
	*dei = 0.;
	*her = -1e300;
	*hei = 0.;
    } else if (*x < 8.) {
	t = *x / 8.;
	t2 = t * t;
	u = t2 * t2;
	*ber = ((((((u * -9.01e-6 + .00122552) * u - .08349609) * u + 
		2.64191397) * u - 32.36345652) * u + 113.77777774) * u - 64.) 
		* u + 1.;
	*bei = t * t * ((((((u * 1.1346e-4 - .01103667) * u + .52185615) * u 
		- 10.56765779) * u + 72.81777742) * u - 113.77777774) * u + 
		16.);
	*ger = ((((((u * -2.458e-5 + .00309699) * u - .19636347) * u + 
		5.65539121) * u - 60.60977451) * u + 171.36272133) * u - 
		59.05819744) * u - .57721566;
	*ger = *ger - log(*x * .5) * *ber + pi * .25 * *bei;
	*gei = t2 * ((((((u * 2.9532e-4 - .02695875) * u + 1.17509064) * u - 
		21.30060904) * u + 124.2356965) * u - 142.91827687) * u + 
		6.76454936);
	*gei = *gei - log(*x * .5) * *bei - pi * .25 * *ber;
	*der = *x * t2 * ((((((u * -3.94e-6 + 4.5957e-4) * u - .02609253) * u 
		+ .66047849) * u - 6.0681481) * u + 14.22222222) * u - 4.);
	*dei = *x * ((((((u * 4.609e-5 - .00379386) * u + .14677204) * u - 
		2.31167514) * u + 11.37777772) * u - 10.66666666) * u + .5);
	*her = *x * t2 * ((((((u * -1.075e-5 + .00116137) * u - .06136358) * 
		u + 1.4138478) * u - 11.36433272) * u + 21.42034017) * u - 
		3.69113734);
	*her = *her - log(*x * .5) * *der - *ber / *x + pi * .25 * *dei;
	*hei = *x * ((((((u * 1.1997e-4 - .00926707) * u + .33049424) * u - 
		4.65950823) * u + 19.41182758) * u - 13.39858846) * u + 
		.21139217);
	*hei = *hei - log(*x * .5) * *dei - *bei / *x - pi * .25 * *der;
    } else {
	t = 8. / *x;
	tnr = 0.;
	tni = 0.;
	for (l = 1; l <= 2; ++l) {
	    v = pow_ii(&c_n1, &l) * t;
	    tpr = ((((v * 6e-7 - 3.4e-6) * v - 2.52e-5) * v - 9.06e-5) * v * 
		    v + .0110486) * v;
	    tpi = ((((v * 1.9e-6 + 5.1e-6) * v * v - 9.01e-5) * v - 9.765e-4) 
		    * v - .0110485) * v - .3926991;
	    if (l == 1) {
		tnr = tpr;
		tni = tpi;
	    }
/* L10: */
	}
	yd = *x / sqrt(2.);
	ye1 = exp(yd + tpr);
	ye2 = exp(-yd + tnr);
	yc1 = 1. / sqrt(pi * 2. * *x);
	yc2 = sqrt(pi / (*x * 2.));
	csp = cos(yd + tpi);
	ssp = sin(yd + tpi);
	csn = cos(-yd + tni);
	ssn = sin(-yd + tni);
	*ger = yc2 * ye2 * csn;
	*gei = yc2 * ye2 * ssn;
	fxr = yc1 * ye1 * csp;
	fxi = yc1 * ye1 * ssp;
	*ber = fxr - *gei / pi;
	*bei = fxi + *ger / pi;
	pnr = 0.;
	pni = 0.;
	for (l = 1; l <= 2; ++l) {
	    v = pow_ii(&c_n1, &l) * t;
	    ppr = (((((v * 1.6e-6 + 1.17e-5) * v + 3.46e-5) * v + 5e-7) * v - 
		    .0013813) * v - .0625001) * v + .7071068;
	    ppi = (((((v * -3.2e-6 - 2.4e-6) * v + 3.38e-5) * v + 2.452e-4) * 
		    v + .0013811) * v - 1e-7) * v + .7071068;
	    if (l == 1) {
		pnr = ppr;
		pni = ppi;
	    }
/* L15: */
	}
	*her = *gei * pni - *ger * pnr;
	*hei = -(*gei * pnr + *ger * pni);
	*der = fxr * ppr - fxi * ppi - *hei / pi;
	*dei = fxi * ppr + fxr * ppi + *her / pi;
    }
    return 0;
} /* klvnb_ */

/*       ********************************** */
/* Subroutine */ int rmn2so_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, doublereal *cv, doublereal *df, integer *kd, 
	doublereal *r2f, doublereal *r2d)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer j;
    static doublereal h0, gd, bk[200], ck[200], gf, dn[200], pi;
    static integer nm, ip;
    static doublereal qs, qt, sw, ck1, ck2, r1d, r1f;
    extern /* Subroutine */ int cbk_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), gmn_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), kmn_(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal eps, sum;
    extern /* Subroutine */ int rmn1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *)
	    , sckb_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *), qstar_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/*       ============================================================= */
/*       Purpose: Compute oblate radial functions of the second kind */
/*                with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix) */
/*       Routines called: */
/*            (1) SCKB for computing the expansion coefficients c2k */
/*            (2) KMN for computing the joining factors */
/*            (3) QSTAR for computing the factor defined in (15.7.3) */
/*            (4) CBK for computing the the expansion coefficient */
/*                defined in (15.7.6) */
/*            (5) GMN for computing the function defined in (15.7.4) */
/*            (6) RMN1 for computing the radial function of the first */
/*                kind */
/*       ============================================================= */

    /* Parameter adjustments */
    --df;

    /* Function Body */
    if (abs(df[1]) <= 1e-280) {
	*r2f = 1e300;
	*r2d = 1e300;
	return 0;
    }
    eps = 1e-14;
    pi = 3.141592653589793;
    nm = (integer) ((*n - *m) / 2 + *c__) + 25;
    ip = 1;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    sckb_(m, n, c__, &df[1], ck);
    kmn_(m, n, c__, cv, kd, &df[1], dn, &ck1, &ck2);
    qstar_(m, n, c__, ck, &ck1, &qs, &qt);
    cbk_(m, n, c__, cv, &qt, ck, bk);
    if (*x == 0.) {
	sum = 0.;
	sw = 0.;
	i__1 = nm;
	for (j = 1; j <= i__1; ++j) {
	    sum += ck[j - 1];
	    if ((d__1 = sum - sw, abs(d__1)) < abs(sum) * eps) {
		goto L15;
	    }
/* L10: */
	    sw = sum;
	}
L15:
	if (ip == 0) {
	    r1f = sum / ck1;
	    *r2f = pi * -.5 * qs * r1f;
	    *r2d = qs * r1f + bk[0];
	} else if (ip == 1) {
	    r1d = sum / ck1;
	    *r2f = bk[0];
	    *r2d = pi * -.5 * qs * r1d;
	}
	return 0;
    } else {
	gmn_(m, n, c__, x, bk, &gf, &gd);
	rmn1_(m, n, c__, x, &df[1], kd, &r1f, &r1d);
	h0 = atan(*x) - pi * .5;
	*r2f = qs * r1f * h0 + gf;
	*r2d = qs * (r1d * h0 + r1f / (*x * *x + 1.)) + gd;
    }
    return 0;
} /* rmn2so_ */

/*       ********************************** */
/* Subroutine */ int bjndd_(integer *n, doublereal *x, doublereal *bj, 
	doublereal *dj, doublereal *fj)
{
    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    double r_lg10(real *), d_lg10(doublereal *);

    /* Local variables */
    static doublereal f;
    static integer k, m;
    static doublereal f0, f1, bs;
    static integer mt, nt;


/*       ===================================================== */
/*       Purpose: Compute Bessel functions Jn(x) and their */
/*                first and second derivatives ( n= 0,1,… ) */
/*       Input:   x ---  Argument of Jn(x)  ( x ≥ 0 ) */
/*                n ---  Order of Jn(x) */
/*       Output:  BJ(n+1) ---  Jn(x) */
/*                DJ(n+1) ---  Jn'(x) */
/*                FJ(n+1) ---  Jn"(x) */
/*       ===================================================== */

    /* Parameter adjustments */
    --fj;
    --dj;
    --bj;

    /* Function Body */
    for (nt = 1; nt <= 900; ++nt) {
	r__1 = nt * 6.28f;
	d__1 = abs(*x) * 1.36f / nt;
	mt = (integer) (r_lg10(&r__1) * .5f - nt * d_lg10(&d__1));
	if (mt > 20) {
	    goto L15;
	}
/* L10: */
    }
L15:
    m = nt;
    bs = 0.;
    f = 0.;
    f0 = 0.;
    f1 = 1e-35;
    for (k = m; k >= 0; --k) {
	f = (k + 1.) * 2. * f1 / *x - f0;
	if (k <= *n) {
	    bj[k + 1] = f;
	}
	if (k == k / 2 << 1) {
	    bs += f * 2.;
	}
	f0 = f1;
/* L20: */
	f1 = f;
    }
    i__1 = *n;
    for (k = 0; k <= i__1; ++k) {
/* L25: */
	bj[k + 1] /= bs - f;
    }
    dj[1] = -bj[2];
    fj[1] = bj[1] * -1. - dj[1] / *x;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	dj[k + 1] = bj[k] - k * bj[k + 1] / *x;
/* L30: */
	fj[k + 1] = (k * k / (*x * *x) - 1.) * bj[k + 1] - dj[k + 1] / *x;
    }
    return 0;
} /* bjndd_ */

/*       ********************************** */
/* Subroutine */ int sphj_(integer *n, doublereal *x, integer *nm, doublereal 
	*sj, doublereal *dj)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, m;
    static doublereal f0, f1, sa, sb, cs;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);

/*       MODIFIED to ALLOW N=0 CASE (ALSO IN SPHY) */

/*       ======================================================= */
/*       Purpose: Compute spherical Bessel functions jn(x) and */
/*                their derivatives */
/*       Input :  x --- Argument of jn(x) */
/*                n --- Order of jn(x)  ( n = 0,1,… ) */
/*       Output:  SJ(n) --- jn(x) */
/*                DJ(n) --- jn'(x) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       ======================================================= */

    *nm = *n;
    if (abs(*x) < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    sj[k] = 0.;
/* L10: */
	    dj[k] = 0.;
	}
	sj[0] = 1.;
	if (*n > 0) {
	    dj[1] = .3333333333333333;
	}
	return 0;
    }
    sj[0] = sin(*x) / *x;
    dj[0] = (cos(*x) - sin(*x) / *x) / *x;
    if (*n < 1) {
	return 0;
    }
    sj[1] = (sj[0] - cos(*x)) / *x;
    if (*n >= 2) {
	sa = sj[0];
	sb = sj[1];
	m = msta1_(x, &c__200);
	if (m < *n) {
	    *nm = m;
	} else {
	    m = msta2_(x, n, &c__15);
	}
	f = 0.;
	f0 = 0.;
	f1 = -99.;
	for (k = m; k >= 0; --k) {
	    f = (k * 2. + 3.) * f1 / *x - f0;
	    if (k <= *nm) {
		sj[k] = f;
	    }
	    f0 = f1;
/* L15: */
	    f1 = f;
	}
	cs = 0.;
	if (abs(sa) > abs(sb)) {
	    cs = sa / f;
	}
	if (abs(sa) <= abs(sb)) {
	    cs = sb / f0;
	}
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
/* L20: */
	    sj[k] = cs * sj[k];
	}
    }
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	dj[k] = sj[k - 1] - (k + 1.) * sj[k] / *x;
    }
    return 0;
} /* sphj_ */

/*       ********************************** */
/* Subroutine */ int othpl_(integer *kf, integer *n, doublereal *x, 
	doublereal *pl, doublereal *dpl)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal a, b, c__;
    static integer k;
    static doublereal y0, y1, yn, dy0, dy1, dyn;


/*       ========================================================== */
/*       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x), */
/*                or Ln(x) or Hn(x), and their derivatives */
/*       Input :  KF --- Function code */
/*                       KF=1 for Chebyshev polynomial Tn(x) */
/*                       KF=2 for Chebyshev polynomial Un(x) */
/*                       KF=3 for Laguerre polynomial Ln(x) */
/*                       KF=4 for Hermite polynomial Hn(x) */
/*                n ---  Order of orthogonal polynomials */
/*                x ---  Argument of orthogonal polynomials */
/*       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x) */
/*                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x) */
/*       ========================================================= */

    a = 2.;
    b = 0.;
    c__ = 1.;
    y0 = 1.;
    y1 = *x * 2.;
    dy0 = 0.;
    dy1 = 2.;
    pl[0] = 1.;
    pl[1] = *x * 2.;
    dpl[0] = 0.;
    dpl[1] = 2.;
    if (*kf == 1) {
	y1 = *x;
	dy1 = 1.;
	pl[1] = *x;
	dpl[1] = 1.;
    } else if (*kf == 3) {
	y1 = 1. - *x;
	dy1 = -1.;
	pl[1] = 1. - *x;
	dpl[1] = -1.;
    }
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	if (*kf == 3) {
	    a = -1. / k;
	    b = a + 2.;
	    c__ = a + 1.;
	} else if (*kf == 4) {
	    c__ = (k - 1.) * 2.;
	}
	yn = (a * *x + b) * y1 - c__ * y0;
	dyn = a * y1 + (a * *x + b) * dy1 - c__ * dy0;
	pl[k] = yn;
	dpl[k] = dyn;
	y0 = y1;
	y1 = yn;
	dy0 = dy1;
/* L10: */
	dy1 = dyn;
    }
    return 0;
} /* othpl_ */

/*       ********************************** */
/* Subroutine */ int klvnzo_(integer *nt, integer *kd, doublereal *zo)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer m;
    static doublereal rt, rt0[8], bei, ddi, dei, gdi, gei, hei, ber, ddr, der,
	     gdr, ger, her;
    extern /* Subroutine */ int klvna_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);


/*       ==================================================== */
/*       Purpose: Compute the zeros of Kelvin functions */
/*       Input :  NT  --- Total number of zeros */
/*                KD  --- Function code */
/*                KD=1 to 8 for ber x, bei x, ker x, kei x, */
/*                          ber'x, bei'x, ker'x and kei'x, */
/*                          respectively. */
/*       Output:  ZO(M) --- the M-th zero of Kelvin function */
/*                          for code KD */
/*       Routine called: */
/*                KLVNA for computing Kelvin functions and */
/*                their derivatives */
/*       ==================================================== */

    /* Parameter adjustments */
    --zo;

    /* Function Body */
    rt0[0] = 2.84891f;
    rt0[1] = 5.02622f;
    rt0[2] = 1.71854f;
    rt0[3] = 3.91467f;
    rt0[4] = 6.03871f;
    rt0[5] = 3.77268f;
    rt0[6] = 2.66584f;
    rt0[7] = 4.93181f;
    rt = rt0[*kd - 1];
    i__1 = *nt;
    for (m = 1; m <= i__1; ++m) {
L10:
	klvna_(&rt, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
	if (*kd == 1) {
	    rt -= ber / der;
	} else if (*kd == 2) {
	    rt -= bei / dei;
	} else if (*kd == 3) {
	    rt -= ger / her;
	} else if (*kd == 4) {
	    rt -= gei / hei;
	} else if (*kd == 5) {
	    ddr = -bei - der / rt;
	    rt -= der / ddr;
	} else if (*kd == 6) {
	    ddi = ber - dei / rt;
	    rt -= dei / ddi;
	} else if (*kd == 7) {
	    gdr = -gei - her / rt;
	    rt -= her / gdr;
	} else {
	    gdi = ger - hei / rt;
	    rt -= hei / gdi;
	}
	if ((d__1 = rt - rt0[*kd - 1], abs(d__1)) > 5e-10) {
	    rt0[*kd - 1] = rt;
	    goto L10;
	}
	zo[m] = rt;
/* L15: */
	rt += 4.44;
    }
    return 0;
} /* klvnzo_ */

/*       ********************************** */
/* Subroutine */ int rswfo_(integer *m, integer *n, doublereal *c__, 
	doublereal *x, doublereal *cv, integer *kf, doublereal *r1f, 
	doublereal *r1d, doublereal *r2f, doublereal *r2d)
{
    static doublereal df[200];
    static integer id, kd;
    extern /* Subroutine */ int rmn1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *)
	    , sdmn_(integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *), rmn2l_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *), rmn2so_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);


/*       ========================================================== */
/*       Purpose: Compute oblate radial functions of the first */
/*                and second kinds, and their derivatives */
/*       Input :  m  --- Mode parameter,  m = 0,1,2,... */
/*                n  --- Mode parameter,  n = m,m+1,m+2,... */
/*                c  --- Spheroidal parameter */
/*                x  --- Argument (x ≥ 0) */
/*                cv --- Characteristic value */
/*                KF --- Function code */
/*                       KF=1 for the first kind */
/*                       KF=2 for the second kind */
/*                       KF=3 for both the first and second kinds */
/*       Output:  R1F --- Radial function of the first kind */
/*                R1D --- Derivative of the radial function of */
/*                        the first kind */
/*                R2F --- Radial function of the second kind */
/*                R2D --- Derivative of the radial function of */
/*                        the second kind */
/*       Routines called: */
/*            (1) SDMN for computing expansion coefficients dk */
/*            (2) RMN1 for computing prolate or oblate radial */
/*                function of the first kind */
/*            (3) RMN2L for computing prolate or oblate radial */
/*                function of the second kind for a large argument */
/*            (4) RMN2SO for computing oblate radial functions of */
/*                the second kind for a small argument */
/*       ========================================================== */

    kd = -1;
    sdmn_(m, n, c__, cv, &kd, df);
    if (*kf != 2) {
	rmn1_(m, n, c__, x, df, &kd, r1f, r1d);
    }
    if (*kf > 1) {
	id = 10;
	if (*x > 1e-8) {
	    rmn2l_(m, n, c__, x, df, &kd, r2f, r2d, &id);
	}
	if (id > -1) {
	    rmn2so_(m, n, c__, x, cv, df, &kd, r2f, r2d);
	}
    }
    return 0;
} /* rswfo_ */

/*       ********************************** */
/* Subroutine */ int ch12n_(integer *n, doublecomplex *z__, integer *nm, 
	doublecomplex *chf1, doublecomplex *chd1, doublecomplex *chf2, 
	doublecomplex *chd2)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex ci;
    static doublereal pi;
    static doublecomplex zi, cf1, cbi[251], cbj[251], cdi[251], cdj[251], cbk[
	    251], cdk[251], cby[251], cdy[251], cfac;
    extern /* Subroutine */ int ciknb_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    ), cjynb_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *);


/*       ==================================================== */
/*       Purpose: Compute Hankel functions of the first and */
/*                second kinds and their derivatives for a */
/*                complex argument */
/*       Input :  z --- Complex argument */
/*                n --- Order of Hn(1)(z) and Hn(2)(z) */
/*       Output:  CHF1(n) --- Hn(1)(z) */
/*                CHD1(n) --- Hn(1)'(z) */
/*                CHF2(n) --- Hn(2)(z) */
/*                CHD2(n) --- Hn(2)'(z) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*             (1) CJYNB for computing Jn(z) and Yn(z) */
/*             (2) CIKNB for computing In(z) and Kn(z) */
/*       ==================================================== */

    ci.r = 0., ci.i = 1.;
    pi = 3.141592653589793;
    if (d_imag(z__) < 0.) {
	cjynb_(n, z__, nm, cbj, cdj, cby, cdy);
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    i__3 = k;
	    i__4 = k;
	    z__2.r = ci.r * cby[i__4].r - ci.i * cby[i__4].i, z__2.i = ci.r * 
		    cby[i__4].i + ci.i * cby[i__4].r;
	    z__1.r = cbj[i__3].r + z__2.r, z__1.i = cbj[i__3].i + z__2.i;
	    chf1[i__2].r = z__1.r, chf1[i__2].i = z__1.i;
/* L10: */
	    i__2 = k;
	    i__3 = k;
	    i__4 = k;
	    z__2.r = ci.r * cdy[i__4].r - ci.i * cdy[i__4].i, z__2.i = ci.r * 
		    cdy[i__4].i + ci.i * cdy[i__4].r;
	    z__1.r = cdj[i__3].r + z__2.r, z__1.i = cdj[i__3].i + z__2.i;
	    chd1[i__2].r = z__1.r, chd1[i__2].i = z__1.i;
	}
	z__1.r = ci.r * z__->r - ci.i * z__->i, z__1.i = ci.r * z__->i + ci.i 
		* z__->r;
	zi.r = z__1.r, zi.i = z__1.i;
	ciknb_(n, &zi, nm, cbi, cdi, cbk, cdk);
	z__2.r = pi * ci.r, z__2.i = pi * ci.i;
	z_div(&z__1, &c_b1162, &z__2);
	cfac.r = z__1.r, cfac.i = z__1.i;
	i__2 = *nm;
	for (k = 0; k <= i__2; ++k) {
	    i__3 = k;
	    i__4 = k;
	    z__1.r = cfac.r * cbk[i__4].r - cfac.i * cbk[i__4].i, z__1.i = 
		    cfac.r * cbk[i__4].i + cfac.i * cbk[i__4].r;
	    chf2[i__3].r = z__1.r, chf2[i__3].i = z__1.i;
	    i__3 = k;
	    z__2.r = cfac.r * ci.r - cfac.i * ci.i, z__2.i = cfac.r * ci.i + 
		    cfac.i * ci.r;
	    i__4 = k;
	    z__1.r = z__2.r * cdk[i__4].r - z__2.i * cdk[i__4].i, z__1.i = 
		    z__2.r * cdk[i__4].i + z__2.i * cdk[i__4].r;
	    chd2[i__3].r = z__1.r, chd2[i__3].i = z__1.i;
/* L15: */
	    z__1.r = cfac.r * ci.r - cfac.i * ci.i, z__1.i = cfac.r * ci.i + 
		    cfac.i * ci.r;
	    cfac.r = z__1.r, cfac.i = z__1.i;
	}
    } else if (d_imag(z__) > 0.) {
	z__2.r = -ci.r, z__2.i = -ci.i;
	z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * z__->i 
		+ z__2.i * z__->r;
	zi.r = z__1.r, zi.i = z__1.i;
	ciknb_(n, &zi, nm, cbi, cdi, cbk, cdk);
	z__1.r = -ci.r, z__1.i = -ci.i;
	cf1.r = z__1.r, cf1.i = z__1.i;
	z__2.r = pi * ci.r, z__2.i = pi * ci.i;
	z_div(&z__1, &c_b15, &z__2);
	cfac.r = z__1.r, cfac.i = z__1.i;
	i__2 = *nm;
	for (k = 0; k <= i__2; ++k) {
	    i__3 = k;
	    i__4 = k;
	    z__1.r = cfac.r * cbk[i__4].r - cfac.i * cbk[i__4].i, z__1.i = 
		    cfac.r * cbk[i__4].i + cfac.i * cbk[i__4].r;
	    chf1[i__3].r = z__1.r, chf1[i__3].i = z__1.i;
	    i__3 = k;
	    z__3.r = -cfac.r, z__3.i = -cfac.i;
	    z__2.r = z__3.r * ci.r - z__3.i * ci.i, z__2.i = z__3.r * ci.i + 
		    z__3.i * ci.r;
	    i__4 = k;
	    z__1.r = z__2.r * cdk[i__4].r - z__2.i * cdk[i__4].i, z__1.i = 
		    z__2.r * cdk[i__4].i + z__2.i * cdk[i__4].r;
	    chd1[i__3].r = z__1.r, chd1[i__3].i = z__1.i;
/* L20: */
	    z__1.r = cfac.r * cf1.r - cfac.i * cf1.i, z__1.i = cfac.r * cf1.i 
		    + cfac.i * cf1.r;
	    cfac.r = z__1.r, cfac.i = z__1.i;
	}
	cjynb_(n, z__, nm, cbj, cdj, cby, cdy);
	i__2 = *nm;
	for (k = 0; k <= i__2; ++k) {
	    i__3 = k;
	    i__4 = k;
	    i__1 = k;
	    z__2.r = ci.r * cby[i__1].r - ci.i * cby[i__1].i, z__2.i = ci.r * 
		    cby[i__1].i + ci.i * cby[i__1].r;
	    z__1.r = cbj[i__4].r - z__2.r, z__1.i = cbj[i__4].i - z__2.i;
	    chf2[i__3].r = z__1.r, chf2[i__3].i = z__1.i;
/* L25: */
	    i__3 = k;
	    i__4 = k;
	    i__1 = k;
	    z__2.r = ci.r * cdy[i__1].r - ci.i * cdy[i__1].i, z__2.i = ci.r * 
		    cdy[i__1].i + ci.i * cdy[i__1].r;
	    z__1.r = cdj[i__4].r - z__2.r, z__1.i = cdj[i__4].i - z__2.i;
	    chd2[i__3].r = z__1.r, chd2[i__3].i = z__1.i;
	}
    } else {
	cjynb_(n, z__, nm, cbj, cdj, cby, cdy);
	i__3 = *nm;
	for (k = 0; k <= i__3; ++k) {
	    i__4 = k;
	    i__1 = k;
	    i__2 = k;
	    z__2.r = ci.r * cby[i__2].r - ci.i * cby[i__2].i, z__2.i = ci.r * 
		    cby[i__2].i + ci.i * cby[i__2].r;
	    z__1.r = cbj[i__1].r + z__2.r, z__1.i = cbj[i__1].i + z__2.i;
	    chf1[i__4].r = z__1.r, chf1[i__4].i = z__1.i;
	    i__4 = k;
	    i__1 = k;
	    i__2 = k;
	    z__2.r = ci.r * cdy[i__2].r - ci.i * cdy[i__2].i, z__2.i = ci.r * 
		    cdy[i__2].i + ci.i * cdy[i__2].r;
	    z__1.r = cdj[i__1].r + z__2.r, z__1.i = cdj[i__1].i + z__2.i;
	    chd1[i__4].r = z__1.r, chd1[i__4].i = z__1.i;
	    i__4 = k;
	    i__1 = k;
	    i__2 = k;
	    z__2.r = ci.r * cby[i__2].r - ci.i * cby[i__2].i, z__2.i = ci.r * 
		    cby[i__2].i + ci.i * cby[i__2].r;
	    z__1.r = cbj[i__1].r - z__2.r, z__1.i = cbj[i__1].i - z__2.i;
	    chf2[i__4].r = z__1.r, chf2[i__4].i = z__1.i;
/* L30: */
	    i__4 = k;
	    i__1 = k;
	    i__2 = k;
	    z__2.r = ci.r * cdy[i__2].r - ci.i * cdy[i__2].i, z__2.i = ci.r * 
		    cdy[i__2].i + ci.i * cdy[i__2].r;
	    z__1.r = cdj[i__1].r - z__2.r, z__1.i = cdj[i__1].i - z__2.i;
	    chd2[i__4].r = z__1.r, chd2[i__4].i = z__1.i;
	}
    }
    return 0;
} /* ch12n_ */

/*       ********************************** */
/* Subroutine */ int jyzo_(integer *n, integer *nt, doublereal *rj0, 
	doublereal *rj1, doublereal *ry0, doublereal *ry1)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer l;
    static doublereal x, x0, pi, bjn, djn, fjn, byn, dyn, fyn;
    extern /* Subroutine */ int jyndd_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal xguess;


/*       ====================================================== */
/*       Purpose: Compute the zeros of Bessel functions Jn(x), */
/*                Yn(x), and their derivatives */
/*       Input :  n  --- Order of Bessel functions  (n >= 0) */
/*                NT --- Number of zeros (roots) */
/*       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT */
/*                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT */
/*                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT */
/*                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT */
/*       Routine called: JYNDD for computing Jn(x), Yn(x), and */
/*                       their first and second derivatives */
/*       ====================================================== */

    /* Parameter adjustments */
    --ry1;
    --ry0;
    --rj1;
    --rj0;

    /* Function Body */
    pi = 3.141592653589793;
/*       -- Newton method for j_{N,L} */
/*       1) initial guess for j_{N,1} */
    if (*n <= 20) {
	x = *n * 1.15859f + 2.82141f;
    } else {
/*          Abr & Stg (9.5.14) */
	d__1 = (doublereal) (*n);
	d__2 = (doublereal) (*n);
	x = *n + pow_dd(&d__1, &c_b1169) * 1.85576f + 1.03315f / pow_dd(&d__2,
		 &c_b1169);
    }
    l = 0;
/*       2) iterate */
    xguess = x;
L10:
    x0 = x;
    jyndd_(n, &x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= bjn / djn;
    if (x - x0 < -1.) {
	x = x0 - 1;
    }
    if (x - x0 > 1.) {
	x = x0 + 1;
    }
    if ((d__1 = x - x0, abs(d__1)) > 1e-11) {
	goto L10;
    }
/*       3) initial guess for j_{N,L+1} */
    if (l >= 1) {
	if (x <= rj0[l] + .5f) {
	    x = xguess + pi;
	    xguess = x;
	    goto L10;
	}
    }
    ++l;
    rj0[l] = x;
/*       XXX: should have a better initial guess for large N ~> 100 here */
/* Computing MAX */
/* Computing 2nd power */
    i__1 = *n;
    d__1 = (*n * .0679f + .0972 - i__1 * i__1 * 3.54e-4f) / l;
    x = x + pi + max(d__1,0.);
    if (l < *nt) {
	goto L10;
    }
/*       -- Newton method for j_{N,L}' */
    if (*n <= 20) {
	x = *n * 1.07703f + .961587f;
    } else {
	d__1 = (doublereal) (*n);
	d__2 = (doublereal) (*n);
	x = *n + pow_dd(&d__1, &c_b1169) * .80861f + .07249f / pow_dd(&d__2, &
		c_b1169);
    }
    if (*n == 0) {
	x = 3.8317f;
    }
    l = 0;
    xguess = x;
L15:
    x0 = x;
    jyndd_(n, &x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= djn / fjn;
    if (x - x0 < -1.) {
	x = x0 - 1;
    }
    if (x - x0 > 1.) {
	x = x0 + 1;
    }
    if ((d__1 = x - x0, abs(d__1)) > 1e-11) {
	goto L15;
    }
    if (l >= 1) {
	if (x <= rj1[l] + .5f) {
	    x = xguess + pi;
	    xguess = x;
	    goto L15;
	}
    }
    ++l;
    rj1[l] = x;
/*       XXX: should have a better initial guess for large N ~> 100 here */
/* Computing MAX */
/* Computing 2nd power */
    i__1 = *n;
    d__1 = (*n * .0915f + .4955 - i__1 * i__1 * 4.35e-4f) / l;
    x = x + pi + max(d__1,0.);
    if (l < *nt) {
	goto L15;
    }
/*       -- Newton method for y_{N,L} */
    if (*n <= 20) {
	x = *n * 1.08933f + 1.19477f;
    } else {
	d__1 = (doublereal) (*n);
	d__2 = (doublereal) (*n);
	x = *n + pow_dd(&d__1, &c_b1169) * .93158f + .26035f / pow_dd(&d__2, &
		c_b1169);
    }
    l = 0;
    xguess = x;
L20:
    x0 = x;
    jyndd_(n, &x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= byn / dyn;
    if (x - x0 < -1.) {
	x = x0 - 1;
    }
    if (x - x0 > 1.) {
	x = x0 + 1;
    }
    if ((d__1 = x - x0, abs(d__1)) > 1e-11) {
	goto L20;
    }
    if (l >= 1) {
	if (x <= ry0[l] + .5f) {
	    x = xguess + pi;
	    xguess = x;
	    goto L20;
	}
    }
    ++l;
    ry0[l] = x;
/*       XXX: should have a better initial guess for large N ~> 100 here */
/* Computing MAX */
/* Computing 2nd power */
    i__1 = *n;
    d__1 = (*n * .0852f + .312 - i__1 * i__1 * 4.03e-4f) / l;
    x = x + pi + max(d__1,0.);
    if (l < *nt) {
	goto L20;
    }
/*       -- Newton method for y_{N,L}' */
    if (*n <= 20) {
	x = *n * 1.16099f + 2.67257f;
    } else {
	d__1 = (doublereal) (*n);
	d__2 = (doublereal) (*n);
	x = *n + pow_dd(&d__1, &c_b1169) * 1.8211f + .94001f / pow_dd(&d__2, &
		c_b1169);
    }
    l = 0;
    xguess = x;
L25:
    x0 = x;
    jyndd_(n, &x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= dyn / fyn;
    if ((d__1 = x - x0, abs(d__1)) > 1e-11) {
	goto L25;
    }
    if (l >= 1) {
	if (x <= ry1[l] + .5f) {
	    x = xguess + pi;
	    xguess = x;
	    goto L25;
	}
    }
    ++l;
    ry1[l] = x;
/*       XXX: should have a better initial guess for large N ~> 100 here */
/* Computing MAX */
/* Computing 2nd power */
    i__1 = *n;
    d__1 = (*n * .0643f + .197 - i__1 * i__1 * 2.86e-4f) / l;
    x = x + pi + max(d__1,0.);
    if (l < *nt) {
	goto L25;
    }
    return 0;
} /* jyzo_ */

/*       ********************************** */
/* Subroutine */ int ikv_(doublereal *v, doublereal *x, doublereal *vm, 
	doublereal *bi, doublereal *di, doublereal *bk, doublereal *dk)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), exp(doublereal), sqrt(
	    doublereal), log(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, m, n;
    static doublereal r__, a1, a2, f1, f2;
    static integer k0;
    static doublereal r1, r2, v0, w0, x2, ca, cb, cs, ct, wa, pi, vt, ww, bi0,
	     bk0, bk1, bk2, v0n, v0p, gan, gap, piv, sum;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ======================================================= */
/*       Purpose: Compute modified Bessel functions Iv(x) and */
/*                Kv(x), and their derivatives */
/*       Input :  x --- Argument ( x ≥ 0 ) */
/*                v --- Order of Iv(x) and Kv(x) */
/*                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 ) */
/*       Output:  BI(n) --- In+v0(x) */
/*                DI(n) --- In+v0'(x) */
/*                BK(n) --- Kn+v0(x) */
/*                DK(n) --- Kn+v0'(x) */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) GAMMA2 for computing the gamma function */
/*            (2) MSTA1 and MSTA2 to compute the starting */
/*                point for backward recurrence */
/*       ======================================================= */

    pi = 3.141592653589793;
    x2 = *x * *x;
    n = (integer) (*v);
    v0 = *v - n;
    if (n == 0) {
	n = 1;
    }
    if (*x < 1e-100) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    bi[k] = 0.;
	    di[k] = 0.;
	    bk[k] = -1e300;
/* L10: */
	    dk[k] = 1e300;
	}
	if (*v == 0.f) {
	    bi[0] = 1.;
	    di[1] = .5;
	}
	*vm = *v;
	return 0;
    }
    piv = pi * v0;
    vt = v0 * 4. * v0;
    if (v0 == 0.) {
	a1 = 1.;
    } else {
	v0p = v0 + 1.;
	gamma2_(&v0p, &gap);
	d__1 = *x * .5;
	a1 = pow_dd(&d__1, &v0) / gap;
    }
    k0 = 14;
    if (*x >= 35.f) {
	k0 = 10;
    }
    if (*x >= 50.f) {
	k0 = 8;
    }
    if (*x <= 18.f) {
	bi0 = 1.;
	r__ = 1.;
	for (k = 1; k <= 30; ++k) {
	    r__ = r__ * .25 * x2 / (k * (k + v0));
	    bi0 += r__;
	    if ((d__1 = r__ / bi0, abs(d__1)) < 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	bi0 *= a1;
    } else {
	ca = exp(*x) / sqrt(pi * 2. * *x);
	sum = 1.;
	r__ = 1.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = k * 2. - 1.;
	    r__ = r__ * -.125 * (vt - pow_dd(&d__1, &c_b4)) / (k * *x);
/* L25: */
	    sum += r__;
	}
	bi0 = ca * sum;
    }
    m = msta1_(x, &c__200);
    if (m < n) {
	n = m;
    } else {
	m = msta2_(x, &n, &c__15);
    }
    f = 0.;
    f2 = 0.;
    f1 = 1e-100;
    ww = 0.;
    for (k = m; k >= 0; --k) {
	f = (v0 + k + 1.) * 2. / *x * f1 + f2;
	if (k <= n) {
	    bi[k] = f;
	}
	f2 = f1;
/* L30: */
	f1 = f;
    }
    cs = bi0 / f;
    i__1 = n;
    for (k = 0; k <= i__1; ++k) {
/* L35: */
	bi[k] = cs * bi[k];
    }
    di[0] = v0 / *x * bi[0] + bi[1];
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L40: */
	di[k] = -(k + v0) / *x * bi[k] + bi[k - 1];
    }
    if (*x <= 9.) {
	if (v0 == 0.) {
	    ct = -log(*x * .5) - .5772156649015329;
	    cs = 0.;
	    w0 = 0.;
	    r__ = 1.;
	    for (k = 1; k <= 50; ++k) {
		w0 += 1. / k;
		r__ = r__ * .25 / (k * k) * x2;
		cs += r__ * (w0 + ct);
		wa = abs(cs);
		if ((d__1 = (wa - ww) / wa, abs(d__1)) < 1e-15) {
		    goto L50;
		}
/* L45: */
		ww = wa;
	    }
L50:
	    bk0 = ct + cs;
	} else {
	    v0n = 1. - v0;
	    gamma2_(&v0n, &gan);
	    d__1 = *x * .5;
	    a2 = 1. / (gan * pow_dd(&d__1, &v0));
	    d__1 = *x * .5;
	    a1 = pow_dd(&d__1, &v0) / gap;
	    sum = a2 - a1;
	    r1 = 1.;
	    r2 = 1.;
	    for (k = 1; k <= 120; ++k) {
		r1 = r1 * .25 * x2 / (k * (k - v0));
		r2 = r2 * .25 * x2 / (k * (k + v0));
		sum = sum + a2 * r1 - a1 * r2;
		wa = abs(sum);
		if ((d__1 = (wa - ww) / wa, abs(d__1)) < 1e-15) {
		    goto L60;
		}
/* L55: */
		ww = wa;
	    }
L60:
	    bk0 = pi * .5 * sum / sin(piv);
	}
    } else {
	cb = exp(-(*x)) * sqrt(pi * .5 / *x);
	sum = 1.;
	r__ = 1.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    d__1 = (doublereal) (k * 2.f - 1.f);
	    r__ = r__ * .125 * (vt - pow_dd(&d__1, &c_b4)) / (k * *x);
/* L65: */
	    sum += r__;
	}
	bk0 = cb * sum;
    }
    bk1 = (1. / *x - bi[1] * bk0) / bi[0];
    bk[0] = bk0;
    bk[1] = bk1;
    i__1 = n;
    for (k = 2; k <= i__1; ++k) {
	bk2 = (v0 + k - 1.) * 2. / *x * bk1 + bk0;
	bk[k] = bk2;
	bk0 = bk1;
/* L70: */
	bk1 = bk2;
    }
    dk[0] = v0 / *x * bk[0] - bk[1];
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L80: */
	dk[k] = -(k + v0) / *x * bk[k] - bk[k - 1];
    }
    *vm = n + v0;
    return 0;
} /* ikv_ */

/*       ********************************** */
/* Subroutine */ int sdmn_(integer *m, integer *n, doublereal *c__, 
	doublereal *cv, integer *kd, doublereal *df)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal a[200], d__[200], f, g[200];
    static integer i__, j, k;
    static doublereal f0, f1, f2;
    static integer k1;
    static doublereal r1, s0, r3, r4;
    static integer kb;
    static doublereal fl, cs;
    static integer ip, nm;
    static doublereal fs, sw, dk0, dk1, dk2, d2k, su1, su2;


/*       ===================================================== */
/*       Purpose: Compute the expansion coefficients of the */
/*                prolate and oblate spheroidal functions, dk */
/*       Input :  m  --- Mode parameter */
/*                n  --- Mode parameter */
/*                c  --- Spheroidal parameter */
/*                cv --- Characteristic value */
/*                KD --- Function code */
/*                       KD=1 for prolate; KD=-1 for oblate */
/*       Output:  DF(k) --- Expansion coefficients dk; */
/*                          DF(1), DF(2), ... correspond to */
/*                          d0, d2, ... for even n-m and d1, */
/*                          d3, ... for odd n-m */
/*       ===================================================== */

    /* Parameter adjustments */
    --df;

    /* Function Body */
    nm = (integer) ((*n - *m) * .5f + *c__) + 25;
    if (*c__ < 1e-10) {
	i__1 = nm;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	    df[i__] = 0.;
	}
	df[(*n - *m) / 2 + 1] = 1.;
	return 0;
    }
    cs = *c__ * *c__ * *kd;
    ip = 1;
    k = 0;
    if (*n - *m == (*n - *m) / 2 << 1) {
	ip = 0;
    }
    i__1 = nm + 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ip == 0) {
	    k = i__ - 1 << 1;
	}
	if (ip == 1) {
	    k = (i__ << 1) - 1;
	}
	dk0 = (doublereal) (*m + k);
	dk1 = (doublereal) (*m + k + 1);
	dk2 = (doublereal) (*m + k << 1);
	d2k = (doublereal) ((*m << 1) + k);
	a[i__ - 1] = (d2k + 2.f) * (d2k + 1.f) / ((dk2 + 3.f) * (dk2 + 5.f)) *
		 cs;
	d__[i__ - 1] = dk0 * dk1 + (dk0 * 2.f * dk1 - *m * 2.f * *m - 1.f) / (
		(dk2 - 1.f) * (dk2 + 3.f)) * cs;
	g[i__ - 1] = k * (k - 1.f) / ((dk2 - 3.f) * (dk2 - 1.f)) * cs;
/* L10: */
    }
    fs = 1.;
    f1 = 0.;
    f0 = 1e-100;
    kb = 0;
    df[nm + 1] = 0.;
    fl = 0.;
    for (k = nm; k >= 1; --k) {
	f = -((d__[k] - *cv) * f0 + a[k] * f1) / g[k];
	if (abs(f) > (d__1 = df[k + 1], abs(d__1))) {
	    df[k] = f;
	    f1 = f0;
	    f0 = f;
	    if (abs(f) > 1e100) {
		i__1 = nm;
		for (k1 = k; k1 <= i__1; ++k1) {
/* L12: */
		    df[k1] *= 1e-100;
		}
		f1 *= 1e-100;
		f0 *= 1e-100;
	    }
	} else {
	    kb = k;
	    fl = df[k + 1];
	    f1 = 1e-100;
	    f2 = -(d__[0] - *cv) / a[0] * f1;
	    df[1] = f1;
	    if (kb == 1) {
		fs = f2;
	    } else if (kb == 2) {
		df[2] = f2;
		fs = -((d__[1] - *cv) * f2 + g[1] * f1) / a[1];
	    } else {
		df[2] = f2;
		i__1 = kb + 1;
		for (j = 3; j <= i__1; ++j) {
		    f = -((d__[j - 2] - *cv) * f2 + g[j - 2] * f1) / a[j - 2];
		    if (j <= kb) {
			df[j] = f;
		    }
		    if (abs(f) > 1e100) {
			i__2 = j;
			for (k1 = 1; k1 <= i__2; ++k1) {
/* L15: */
			    df[k1] *= 1e-100;
			}
			f *= 1e-100;
			f2 *= 1e-100;
		    }
		    f1 = f2;
/* L20: */
		    f2 = f;
		}
		fs = f;
	    }
	    goto L35;
	}
/* L30: */
    }
L35:
    su1 = 0.;
    r1 = 1.;
    i__1 = *m + ip << 1;
    for (j = *m + ip + 1; j <= i__1; ++j) {
/* L40: */
	r1 *= j;
    }
    su1 = df[1] * r1;
    i__1 = kb;
    for (k = 2; k <= i__1; ++k) {
	r1 = -r1 * (k + *m + ip - 1.5) / (k - 1.);
/* L45: */
	su1 += r1 * df[k];
    }
    su2 = 0.;
    sw = 0.;
    i__1 = nm;
    for (k = kb + 1; k <= i__1; ++k) {
	if (k != 1) {
	    r1 = -r1 * (k + *m + ip - 1.5) / (k - 1.);
	}
	su2 += r1 * df[k];
	if ((d__1 = sw - su2, abs(d__1)) < abs(su2) * 1e-14) {
	    goto L55;
	}
/* L50: */
	sw = su2;
    }
L55:
    r3 = 1.;
    i__1 = (*m + *n + ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L60: */
	r3 *= j + (*n + *m + ip) * .5;
    }
    r4 = 1.;
    i__1 = (*n - *m - ip) / 2;
    for (j = 1; j <= i__1; ++j) {
/* L65: */
	r4 = r4 * -4. * j;
    }
    s0 = r3 / (fl * (su1 / fs) + su2) / r4;
    i__1 = kb;
    for (k = 1; k <= i__1; ++k) {
/* L70: */
	df[k] = fl / fs * s0 * df[k];
    }
    i__1 = nm;
    for (k = kb + 1; k <= i__1; ++k) {
/* L75: */
	df[k] = s0 * df[k];
    }
    return 0;
} /* sdmn_ */

/*       ********************************** */
/* Subroutine */ int ajyik_(doublereal *x, doublereal *vj1, doublereal *vj2, 
	doublereal *vy1, doublereal *vy2, doublereal *vi1, doublereal *vi2, 
	doublereal *vk1, doublereal *vk2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), cos(
	    doublereal), sin(doublereal), exp(doublereal);

    /* Local variables */
    static integer k, l;
    static doublereal r__, a0, b0, c0;
    static integer k0;
    static doublereal x2, ck, gn, pi, sk, vl, rp, rq, xk, px, qx, vv, gn1, 
	    gn2, gp1, gp2, rp2, uj1, uj2, pv1, uu0, pv2, vv0, vil, vjl, vsl, 
	    sum;


/*       ======================================================= */
/*       Purpose: Compute Bessel functions Jv(x) and Yv(x), */
/*                and modified Bessel functions Iv(x) and */
/*                Kv(x), and their derivatives with v=1/3,2/3 */
/*       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and */
/*                      Kv(x) ( x ≥ 0 ) */
/*       Output:  VJ1 --- J1/3(x) */
/*                VJ2 --- J2/3(x) */
/*                VY1 --- Y1/3(x) */
/*                VY2 --- Y2/3(x) */
/*                VI1 --- I1/3(x) */
/*                VI2 --- I2/3(x) */
/*                VK1 --- K1/3(x) */
/*                VK2 --- K2/3(x) */
/*       ======================================================= */

    if (*x == 0.) {
	*vj1 = 0.;
	*vj2 = 0.;
	*vy1 = -1e300;
	*vy2 = 1e300;
	*vi1 = 0.;
	*vi2 = 0.;
	*vk1 = -1e300;
	*vk2 = -1e300;
	return 0;
    }
    pi = 3.141592653589793;
    rp2 = .63661977236758;
    gp1 = .892979511569249;
    gp2 = .902745292950934;
    gn1 = 1.3541179394264;
    gn2 = 2.678938534707747;
    vv0 = .444444444444444;
    uu0 = 1.1547005383793;
    x2 = *x * *x;
    k0 = 12;
    if (*x >= 35.f) {
	k0 = 10;
    }
    if (*x >= 50.f) {
	k0 = 8;
    }
    if (*x <= 12.f) {
	for (l = 1; l <= 2; ++l) {
	    vl = l / 3.;
	    vjl = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 40; ++k) {
		r__ = r__ * -.25 * x2 / (k * (k + vl));
		vjl += r__;
		if (abs(r__) < 1e-15) {
		    goto L20;
		}
/* L15: */
	    }
L20:
	    d__1 = *x * .5;
	    a0 = pow_dd(&d__1, &vl);
	    if (l == 1) {
		*vj1 = a0 / gp1 * vjl;
	    }
	    if (l == 2) {
		*vj2 = a0 / gp2 * vjl;
	    }
/* L25: */
	}
    } else {
	for (l = 1; l <= 2; ++l) {
	    vv = vv0 * l * l;
	    px = 1.;
	    rp = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = (doublereal) (k * 4.f - 3.f);
		d__2 = (doublereal) (k * 4.f - 1.f);
		rp = rp * -.0078125 * (vv - pow_dd(&d__1, &c_b4)) * (vv - 
			pow_dd(&d__2, &c_b4)) / (k * (k * 2.f - 1.f) * x2);
/* L30: */
		px += rp;
	    }
	    qx = 1.;
	    rq = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = (doublereal) (k * 4.f - 1.f);
		d__2 = (doublereal) (k * 4.f + 1.f);
		rq = rq * -.0078125 * (vv - pow_dd(&d__1, &c_b4)) * (vv - 
			pow_dd(&d__2, &c_b4)) / (k * (k * 2.f + 1.f) * x2);
/* L35: */
		qx += rq;
	    }
	    qx = (vv - 1.f) * .125 * qx / *x;
	    xk = *x - (l * .5 / 3. + .25) * pi;
	    a0 = sqrt(rp2 / *x);
	    ck = cos(xk);
	    sk = sin(xk);
	    if (l == 1) {
		*vj1 = a0 * (px * ck - qx * sk);
		*vy1 = a0 * (px * sk + qx * ck);
	    } else if (l == 2) {
		*vj2 = a0 * (px * ck - qx * sk);
		*vy2 = a0 * (px * sk + qx * ck);
	    }
/* L40: */
	}
    }
    if (*x <= 12.) {
	uj1 = 0.;
	uj2 = 0.;
	for (l = 1; l <= 2; ++l) {
	    vl = l / 3.;
	    vjl = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 40; ++k) {
		r__ = r__ * -.25 * x2 / (k * (k - vl));
		vjl += r__;
		if (abs(r__) < 1e-15) {
		    goto L50;
		}
/* L45: */
	    }
L50:
	    d__1 = 2. / *x;
	    b0 = pow_dd(&d__1, &vl);
	    if (l == 1) {
		uj1 = b0 * vjl / gn1;
	    }
	    if (l == 2) {
		uj2 = b0 * vjl / gn2;
	    }
/* L55: */
	}
	pv1 = pi / 3.;
	pv2 = pi / 1.5;
	*vy1 = uu0 * (*vj1 * cos(pv1) - uj1);
	*vy2 = uu0 * (*vj2 * cos(pv2) - uj2);
    }
    if (*x <= 18.f) {
	for (l = 1; l <= 2; ++l) {
	    vl = l / 3.;
	    vil = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 40; ++k) {
		r__ = r__ * .25 * x2 / (k * (k + vl));
		vil += r__;
		if (abs(r__) < 1e-15) {
		    goto L65;
		}
/* L60: */
	    }
L65:
	    d__1 = *x * .5;
	    a0 = pow_dd(&d__1, &vl);
	    if (l == 1) {
		*vi1 = a0 / gp1 * vil;
	    }
	    if (l == 2) {
		*vi2 = a0 / gp2 * vil;
	    }
/* L70: */
	}
    } else {
	c0 = exp(*x) / sqrt(pi * 2. * *x);
	for (l = 1; l <= 2; ++l) {
	    vv = vv0 * l * l;
	    vsl = 1.;
	    r__ = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = k * 2. - 1.;
		r__ = r__ * -.125 * (vv - pow_dd(&d__1, &c_b4)) / (k * *x);
/* L75: */
		vsl += r__;
	    }
	    if (l == 1) {
		*vi1 = c0 * vsl;
	    }
	    if (l == 2) {
		*vi2 = c0 * vsl;
	    }
/* L80: */
	}
    }
    if (*x <= 9.) {
	gn = 0.;
	for (l = 1; l <= 2; ++l) {
	    vl = l / 3.;
	    if (l == 1) {
		gn = gn1;
	    }
	    if (l == 2) {
		gn = gn2;
	    }
	    d__1 = 2. / *x;
	    a0 = pow_dd(&d__1, &vl) / gn;
	    sum = 1.;
	    r__ = 1.;
	    for (k = 1; k <= 60; ++k) {
		r__ = r__ * .25 * x2 / (k * (k - vl));
		sum += r__;
		if (abs(r__) < 1e-15) {
		    goto L90;
		}
/* L85: */
	    }
L90:
	    if (l == 1) {
		*vk1 = uu0 * .5 * pi * (sum * a0 - *vi1);
	    }
	    if (l == 2) {
		*vk2 = uu0 * .5 * pi * (sum * a0 - *vi2);
	    }
/* L95: */
	}
    } else {
	c0 = exp(-(*x)) * sqrt(pi * .5 / *x);
	for (l = 1; l <= 2; ++l) {
	    vv = vv0 * l * l;
	    sum = 1.;
	    r__ = 1.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = (doublereal) (k * 2.f - 1.f);
		r__ = r__ * .125 * (vv - pow_dd(&d__1, &c_b4)) / (k * *x);
/* L100: */
		sum += r__;
	    }
	    if (l == 1) {
		*vk1 = c0 * sum;
	    }
	    if (l == 2) {
		*vk2 = c0 * sum;
	    }
/* L105: */
	}
    }
    return 0;
} /* ajyik_ */

/*       ********************************** */
/* Subroutine */ int cikvb_(doublereal *v, doublecomplex *z__, doublereal *vm,
	 doublecomplex *cbi, doublecomplex *cdi, doublecomplex *cbk, 
	doublecomplex *cdk)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *), z_div(
	    doublecomplex *, doublecomplex *, doublecomplex *), z_exp(
	    doublecomplex *, doublecomplex *), z_sqrt(doublecomplex *, 
	    doublecomplex *);
    double pow_dd(doublereal *, doublereal *);
    void z_log(doublecomplex *, doublecomplex *);
    double sin(doublereal), d_imag(doublecomplex *);

    /* Local variables */
    static integer k, m, n;
    static doublereal a0;
    static integer k0;
    static doublereal v0, w0;
    static doublecomplex z1, z2, ca, cb, cf, ci, cp, cr, cs, ct;
    static doublereal pi, vt;
    static doublecomplex ca1, ca2, cf1, cf2, ci0, cr1, cr2;
    static doublereal v0n, v0p, gan, gap;
    static doublecomplex ckk, cvk, csu;
    static doublereal piv;
    static doublecomplex cbi0, cbk0;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       =========================================================== */
/*       Purpose: Compute the modified Bessel functions Iv(z), Kv(z) */
/*                and their derivatives for an arbitrary order and */
/*                complex argument */
/*       Input :  z --- Complex argument z */
/*                v --- Real order of Iv(z) and Kv(z) */
/*                      ( v =n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 ) */
/*       Output:  CBI(n) --- In+v0(z) */
/*                CDI(n) --- In+v0'(z) */
/*                CBK(n) --- Kn+v0(z) */
/*                CDK(n) --- Kn+v0'(z) */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) GAMMA2 for computing the gamma function */
/*            (2) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       =========================================================== */

    z1.r = z__->r, z1.i = z__->i;
    z__1.r = z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = z__1.r, z2.i = z__1.i;
    a0 = z_abs(z__);
    pi = 3.141592653589793;
    ci.r = 0., ci.i = 1.;
    n = (integer) (*v);
    v0 = *v - n;
    piv = pi * v0;
    vt = v0 * 4. * v0;
    if (n == 0) {
	n = 1;
    }
    if (a0 < 1e-100) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbi[i__2].r = 0., cbi[i__2].i = 0.;
	    i__2 = k;
	    cdi[i__2].r = 0., cdi[i__2].i = 0.;
	    i__2 = k;
	    cbk[i__2].r = -1e300, cbk[i__2].i = 0.;
/* L10: */
	    i__2 = k;
	    cdk[i__2].r = 1e300, cdk[i__2].i = 0.;
	}
	if (v0 == 0.f) {
	    cbi[0].r = 1., cbi[0].i = 0.;
	    cdi[1].r = .5, cdi[1].i = 0.;
	}
	*vm = *v;
	return 0;
    }
    k0 = 14;
    if (a0 >= 35.f) {
	k0 = 10;
    }
    if (a0 >= 50.f) {
	k0 = 8;
    }
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (a0 < 18.f) {
	if (v0 == 0.f) {
	    ca1.r = 1., ca1.i = 0.;
	} else {
	    v0p = v0 + 1.;
	    gamma2_(&v0p, &gap);
	    z__3.r = z1.r * .5, z__3.i = z1.i * .5;
	    z__4.r = v0, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z__1.r = z__2.r / gap, z__1.i = z__2.i / gap;
	    ca1.r = z__1.r, ca1.i = z__1.i;
	}
	ci0.r = 1., ci0.i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 50; ++k) {
	    z__3.r = cr.r * .25, z__3.i = cr.i * .25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    d__1 = k * (k + v0);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = ci0.r + cr.r, z__1.i = ci0.i + cr.i;
	    ci0.r = z__1.r, ci0.i = z__1.i;
	    z_div(&z__1, &cr, &ci0);
	    if (z_abs(&z__1) < 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	z__1.r = ci0.r * ca1.r - ci0.i * ca1.i, z__1.i = ci0.r * ca1.i + 
		ci0.i * ca1.r;
	cbi0.r = z__1.r, cbi0.i = z__1.i;
    } else {
	z_exp(&z__2, &z1);
	d__1 = pi * 2.;
	z__4.r = d__1 * z1.r, z__4.i = d__1 * z1.i;
	z_sqrt(&z__3, &z__4);
	z_div(&z__1, &z__2, &z__3);
	ca.r = z__1.r, ca.i = z__1.i;
	cs.r = 1., cs.i = 0.;
	cr.r = 1., cr.i = 0.;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
	    z__3.r = cr.r * -.125, z__3.i = cr.i * -.125;
	    d__2 = k * 2. - 1.;
	    d__1 = vt - pow_dd(&d__2, &c_b4);
	    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	    d__3 = (doublereal) k;
	    z__4.r = d__3 * z1.r, z__4.i = d__3 * z1.i;
	    z_div(&z__1, &z__2, &z__4);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L25: */
	    z__1.r = cs.r + cr.r, z__1.i = cs.i + cr.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	}
	z__1.r = ca.r * cs.r - ca.i * cs.i, z__1.i = ca.r * cs.i + ca.i * 
		cs.r;
	cbi0.r = z__1.r, cbi0.i = z__1.i;
    }
    m = msta1_(&a0, &c__200);
    if (m < n) {
	n = m;
    } else {
	m = msta2_(&a0, &n, &c__15);
    }
    cf2.r = 0., cf2.i = 0.;
    cf1.r = 1e-100, cf1.i = 0.;
    for (k = m; k >= 0; --k) {
	d__1 = (v0 + k + 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, &z1);
	z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i + 
		z__3.i * cf1.r;
	z__1.r = z__2.r + cf2.r, z__1.i = z__2.i + cf2.i;
	cf.r = z__1.r, cf.i = z__1.i;
	if (k <= n) {
	    i__2 = k;
	    cbi[i__2].r = cf.r, cbi[i__2].i = cf.i;
	}
	cf2.r = cf1.r, cf2.i = cf1.i;
/* L30: */
	cf1.r = cf.r, cf1.i = cf.i;
    }
    z_div(&z__1, &cbi0, &cf);
    cs.r = z__1.r, cs.i = z__1.i;
    i__2 = n;
    for (k = 0; k <= i__2; ++k) {
/* L35: */
	i__1 = k;
	i__3 = k;
	z__1.r = cs.r * cbi[i__3].r - cs.i * cbi[i__3].i, z__1.i = cs.r * cbi[
		i__3].i + cs.i * cbi[i__3].r;
	cbi[i__1].r = z__1.r, cbi[i__1].i = z__1.i;
    }
    if (a0 <= 9.f) {
	if (v0 == 0.f) {
	    z__4.r = z1.r * .5, z__4.i = z1.i * .5;
	    z_log(&z__3, &z__4);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z__1.r = z__2.r - .5772156649015329, z__1.i = z__2.i;
	    ct.r = z__1.r, ct.i = z__1.i;
	    cs.r = 0., cs.i = 0.;
	    w0 = 0.;
	    cr.r = 1., cr.i = 0.;
	    for (k = 1; k <= 50; ++k) {
		w0 += 1. / k;
		z__3.r = cr.r * .25, z__3.i = cr.i * .25;
		i__1 = k * k;
		d__1 = (doublereal) i__1;
		z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
		z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * 
			z2.i + z__2.i * z2.r;
		cr.r = z__1.r, cr.i = z__1.i;
		z__2.r = w0 + ct.r, z__2.i = ct.i;
		z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * 
			z__2.i + cr.i * z__2.r;
		cp.r = z__1.r, cp.i = z__1.i;
		z__1.r = cs.r + cp.r, z__1.i = cs.i + cp.i;
		cs.r = z__1.r, cs.i = z__1.i;
		z_div(&z__1, &cp, &cs);
		if (k >= 10 && z_abs(&z__1) < 1e-15) {
		    goto L45;
		}
/* L40: */
	    }
L45:
	    z__1.r = ct.r + cs.r, z__1.i = ct.i + cs.i;
	    cbk0.r = z__1.r, cbk0.i = z__1.i;
	} else {
	    v0n = 1. - v0;
	    gamma2_(&v0n, &gan);
	    z__4.r = z1.r * .5, z__4.i = z1.i * .5;
	    z__5.r = v0, z__5.i = 0.;
	    pow_zz(&z__3, &z__4, &z__5);
	    z__2.r = gan * z__3.r, z__2.i = gan * z__3.i;
	    z_div(&z__1, &c_b6, &z__2);
	    ca2.r = z__1.r, ca2.i = z__1.i;
	    z__3.r = z1.r * .5, z__3.i = z1.i * .5;
	    z__4.r = v0, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z__1.r = z__2.r / gap, z__1.i = z__2.i / gap;
	    ca1.r = z__1.r, ca1.i = z__1.i;
	    z__1.r = ca2.r - ca1.r, z__1.i = ca2.i - ca1.i;
	    csu.r = z__1.r, csu.i = z__1.i;
	    cr1.r = 1., cr1.i = 0.;
	    cr2.r = 1., cr2.i = 0.;
	    for (k = 1; k <= 50; ++k) {
		z__3.r = cr1.r * .25, z__3.i = cr1.i * .25;
		z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			z2.i + z__3.i * z2.r;
		d__1 = k * (k - v0);
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr1.r = z__1.r, cr1.i = z__1.i;
		z__3.r = cr2.r * .25, z__3.i = cr2.i * .25;
		z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			z2.i + z__3.i * z2.r;
		d__1 = k * (k + v0);
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr2.r = z__1.r, cr2.i = z__1.i;
		z__2.r = ca2.r * cr1.r - ca2.i * cr1.i, z__2.i = ca2.r * 
			cr1.i + ca2.i * cr1.r;
		z__3.r = ca1.r * cr2.r - ca1.i * cr2.i, z__3.i = ca1.r * 
			cr2.i + ca1.i * cr2.r;
		z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		cp.r = z__1.r, cp.i = z__1.i;
		z__1.r = csu.r + cp.r, z__1.i = csu.i + cp.i;
		csu.r = z__1.r, csu.i = z__1.i;
		z_div(&z__1, &cp, &csu);
		if (k >= 10 && z_abs(&z__1) < 1e-15) {
		    goto L55;
		}
/* L50: */
	    }
L55:
	    d__1 = pi * .5;
	    z__2.r = d__1 * csu.r, z__2.i = d__1 * csu.i;
	    d__2 = sin(piv);
	    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	    cbk0.r = z__1.r, cbk0.i = z__1.i;
	}
    } else {
	z__3.r = -z1.r, z__3.i = -z1.i;
	z_exp(&z__2, &z__3);
	d__1 = pi * .5;
	z__6.r = d__1, z__6.i = 0.;
	z_div(&z__5, &z__6, &z1);
	z_sqrt(&z__4, &z__5);
	z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * z__4.i 
		+ z__2.i * z__4.r;
	cb.r = z__1.r, cb.i = z__1.i;
	cs.r = 1., cs.i = 0.;
	cr.r = 1., cr.i = 0.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    z__3.r = cr.r * .125, z__3.i = cr.i * .125;
	    d__2 = k * 2. - 1.;
	    d__1 = vt - pow_dd(&d__2, &c_b4);
	    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	    d__3 = (doublereal) k;
	    z__4.r = d__3 * z1.r, z__4.i = d__3 * z1.i;
	    z_div(&z__1, &z__2, &z__4);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L60: */
	    z__1.r = cs.r + cr.r, z__1.i = cs.i + cr.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	}
	z__1.r = cb.r * cs.r - cb.i * cs.i, z__1.i = cb.r * cs.i + cb.i * 
		cs.r;
	cbk0.r = z__1.r, cbk0.i = z__1.i;
    }
    cbk[0].r = cbk0.r, cbk[0].i = cbk0.i;
    if (z__->r < 0.f) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    d__1 = (k + v0) * pi;
	    z__2.r = d__1 * ci.r, z__2.i = d__1 * ci.i;
	    z_exp(&z__1, &z__2);
	    cvk.r = z__1.r, cvk.i = z__1.i;
	    if (d_imag(z__) < 0.) {
		i__3 = k;
		i__2 = k;
		z__2.r = cvk.r * cbk[i__2].r - cvk.i * cbk[i__2].i, z__2.i = 
			cvk.r * cbk[i__2].i + cvk.i * cbk[i__2].r;
		z__4.r = pi * ci.r, z__4.i = pi * ci.i;
		i__4 = k;
		z__3.r = z__4.r * cbi[i__4].r - z__4.i * cbi[i__4].i, z__3.i =
			 z__4.r * cbi[i__4].i + z__4.i * cbi[i__4].r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		cbk[i__3].r = z__1.r, cbk[i__3].i = z__1.i;
		i__3 = k;
		z_div(&z__1, &cbi[k], &cvk);
		cbi[i__3].r = z__1.r, cbi[i__3].i = z__1.i;
	    } else if (d_imag(z__) > 0.f) {
		i__3 = k;
		z_div(&z__2, &cbk[k], &cvk);
		z__4.r = pi * ci.r, z__4.i = pi * ci.i;
		i__2 = k;
		z__3.r = z__4.r * cbi[i__2].r - z__4.i * cbi[i__2].i, z__3.i =
			 z__4.r * cbi[i__2].i + z__4.i * cbi[i__2].r;
		z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		cbk[i__3].r = z__1.r, cbk[i__3].i = z__1.i;
		i__3 = k;
		i__2 = k;
		z__1.r = cvk.r * cbi[i__2].r - cvk.i * cbi[i__2].i, z__1.i = 
			cvk.r * cbi[i__2].i + cvk.i * cbi[i__2].r;
		cbi[i__3].r = z__1.r, cbi[i__3].i = z__1.i;
	    }
/* L65: */
	}
    }
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	z_div(&z__3, &c_b6, z__);
	i__3 = k;
	i__2 = k - 1;
	z__4.r = cbi[i__3].r * cbk[i__2].r - cbi[i__3].i * cbk[i__2].i, 
		z__4.i = cbi[i__3].r * cbk[i__2].i + cbi[i__3].i * cbk[i__2]
		.r;
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z_div(&z__1, &z__2, &cbi[k - 1]);
	ckk.r = z__1.r, ckk.i = z__1.i;
	i__3 = k;
	cbk[i__3].r = ckk.r, cbk[i__3].i = ckk.i;
/* L70: */
    }
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbi[0].r - z__3.i * cbi[0].i, z__2.i = z__3.r * cbi[0]
	    .i + z__3.i * cbi[0].r;
    z__1.r = z__2.r + cbi[1].r, z__1.i = z__2.i + cbi[1].i;
    cdi[0].r = z__1.r, cdi[0].i = z__1.i;
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbk[0].r - z__3.i * cbk[0].i, z__2.i = z__3.r * cbk[0]
	    .i + z__3.i * cbk[0].r;
    z__1.r = z__2.r - cbk[1].r, z__1.i = z__2.i - cbk[1].i;
    cdk[0].r = z__1.r, cdk[0].i = z__1.i;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__3 = k;
	d__1 = -(k + v0);
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = k;
	z__2.r = z__3.r * cbi[i__2].r - z__3.i * cbi[i__2].i, z__2.i = z__3.r 
		* cbi[i__2].i + z__3.i * cbi[i__2].r;
	i__4 = k - 1;
	z__1.r = z__2.r + cbi[i__4].r, z__1.i = z__2.i + cbi[i__4].i;
	cdi[i__3].r = z__1.r, cdi[i__3].i = z__1.i;
/* L80: */
	i__3 = k;
	d__1 = -(k + v0);
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = k;
	z__2.r = z__3.r * cbk[i__2].r - z__3.i * cbk[i__2].i, z__2.i = z__3.r 
		* cbk[i__2].i + z__3.i * cbk[i__2].r;
	i__4 = k - 1;
	z__1.r = z__2.r - cbk[i__4].r, z__1.i = z__2.i - cbk[i__4].i;
	cdk[i__3].r = z__1.r, cdk[i__3].i = z__1.i;
    }
    *vm = n + v0;
    return 0;
} /* cikvb_ */

/*       ********************************** */
/* Subroutine */ int cikva_(doublereal *v, doublecomplex *z__, doublereal *vm,
	 doublecomplex *cbi, doublecomplex *cdi, doublecomplex *cbk, 
	doublecomplex *cdk)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *), z_exp(
	    doublecomplex *, doublecomplex *), z_sqrt(doublecomplex *, 
	    doublecomplex *), z_div(doublecomplex *, doublecomplex *, 
	    doublecomplex *);
    double pow_dd(doublereal *, doublereal *);
    void z_log(doublecomplex *, doublecomplex *);
    double sin(doublereal), d_imag(doublecomplex *);

    /* Local variables */
    static integer k, m, n;
    static doublereal a0;
    static integer k0;
    static doublereal v0, w0;
    static doublecomplex z1, z2, ca, cb, cf, ci, cp, cr, cs, ct;
    static doublereal pi, vt, ws;
    static doublecomplex ca1, ca2, cf1, cf2, ci0, cg0, cg1, cr1, cr2;
    static doublereal v0n, v0p, ws0;
    static doublecomplex cgk;
    static doublereal gan, gap;
    static doublecomplex cvk, csu;
    static doublereal piv;
    static doublecomplex cbi0, cbk0, cbk1;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int gamma2_(doublereal *, doublereal *);


/*       ============================================================ */
/*       Purpose: Compute the modified Bessel functions Iv(z), Kv(z) */
/*                and their derivatives for an arbitrary order and */
/*                complex argument */
/*       Input :  z --- Complex argument */
/*                v --- Real order of Iv(z) and Kv(z) */
/*                      ( v = n+v0, n = 0,1,2,…, 0 ≤ v0 < 1 ) */
/*       Output:  CBI(n) --- In+v0(z) */
/*                CDI(n) --- In+v0'(z) */
/*                CBK(n) --- Kn+v0(z) */
/*                CDK(n) --- Kn+v0'(z) */
/*                VM --- Highest order computed */
/*       Routines called: */
/*            (1) GAMMA2 for computing the gamma function */
/*            (2) MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       ============================================================ */

    pi = 3.141592653589793;
    ci.r = 0., ci.i = 1.;
    a0 = z_abs(z__);
    z1.r = z__->r, z1.i = z__->i;
    z__1.r = z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = z__1.r, z2.i = z__1.i;
    n = (integer) (*v);
    v0 = *v - n;
    piv = pi * v0;
    vt = v0 * 4. * v0;
    if (n == 0) {
	n = 1;
    }
    if (a0 < 1e-100) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbi[i__2].r = 0., cbi[i__2].i = 0.;
	    i__2 = k;
	    cdi[i__2].r = 0., cdi[i__2].i = 0.;
	    i__2 = k;
	    cbk[i__2].r = -1e300, cbk[i__2].i = 0.;
/* L10: */
	    i__2 = k;
	    cdk[i__2].r = 1e300, cdk[i__2].i = 0.;
	}
	if (v0 == 0.f) {
	    cbi[0].r = 1., cbi[0].i = 0.;
	    cdi[1].r = .5, cdi[1].i = 0.;
	}
	*vm = *v;
	return 0;
    }
    k0 = 14;
    if (a0 >= 35.f) {
	k0 = 10;
    }
    if (a0 >= 50.f) {
	k0 = 8;
    }
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (a0 < 18.f) {
	if (v0 == 0.f) {
	    ca1.r = 1., ca1.i = 0.;
	} else {
	    v0p = v0 + 1.;
	    gamma2_(&v0p, &gap);
	    z__3.r = z1.r * .5, z__3.i = z1.i * .5;
	    z__4.r = v0, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z__1.r = z__2.r / gap, z__1.i = z__2.i / gap;
	    ca1.r = z__1.r, ca1.i = z__1.i;
	}
	ci0.r = 1., ci0.i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 50; ++k) {
	    z__3.r = cr.r * .25, z__3.i = cr.i * .25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    d__1 = k * (k + v0);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = ci0.r + cr.r, z__1.i = ci0.i + cr.i;
	    ci0.r = z__1.r, ci0.i = z__1.i;
	    if (z_abs(&cr) < z_abs(&ci0) * 1e-15) {
		goto L20;
	    }
/* L15: */
	}
L20:
	z__1.r = ci0.r * ca1.r - ci0.i * ca1.i, z__1.i = ci0.r * ca1.i + 
		ci0.i * ca1.r;
	cbi0.r = z__1.r, cbi0.i = z__1.i;
    } else {
	z_exp(&z__2, &z1);
	d__1 = pi * 2.;
	z__4.r = d__1 * z1.r, z__4.i = d__1 * z1.i;
	z_sqrt(&z__3, &z__4);
	z_div(&z__1, &z__2, &z__3);
	ca.r = z__1.r, ca.i = z__1.i;
	cs.r = 1., cs.i = 0.;
	cr.r = 1., cr.i = 0.;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
	    z__3.r = cr.r * -.125, z__3.i = cr.i * -.125;
	    d__2 = k * 2. - 1.;
	    d__1 = vt - pow_dd(&d__2, &c_b4);
	    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	    d__3 = (doublereal) k;
	    z__4.r = d__3 * z1.r, z__4.i = d__3 * z1.i;
	    z_div(&z__1, &z__2, &z__4);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L25: */
	    z__1.r = cs.r + cr.r, z__1.i = cs.i + cr.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	}
	z__1.r = ca.r * cs.r - ca.i * cs.i, z__1.i = ca.r * cs.i + ca.i * 
		cs.r;
	cbi0.r = z__1.r, cbi0.i = z__1.i;
    }
    m = msta1_(&a0, &c__200);
    if (m < n) {
	n = m;
    } else {
	m = msta2_(&a0, &n, &c__15);
    }
    cf2.r = 0., cf2.i = 0.;
    cf1.r = 1e-100, cf1.i = 0.;
    for (k = m; k >= 0; --k) {
	d__1 = (v0 + k + 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, &z1);
	z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i + 
		z__3.i * cf1.r;
	z__1.r = z__2.r + cf2.r, z__1.i = z__2.i + cf2.i;
	cf.r = z__1.r, cf.i = z__1.i;
	if (k <= n) {
	    i__2 = k;
	    cbi[i__2].r = cf.r, cbi[i__2].i = cf.i;
	}
	cf2.r = cf1.r, cf2.i = cf1.i;
/* L30: */
	cf1.r = cf.r, cf1.i = cf.i;
    }
    z_div(&z__1, &cbi0, &cf);
    cs.r = z__1.r, cs.i = z__1.i;
    i__2 = n;
    for (k = 0; k <= i__2; ++k) {
/* L35: */
	i__1 = k;
	i__3 = k;
	z__1.r = cs.r * cbi[i__3].r - cs.i * cbi[i__3].i, z__1.i = cs.r * cbi[
		i__3].i + cs.i * cbi[i__3].r;
	cbi[i__1].r = z__1.r, cbi[i__1].i = z__1.i;
    }
    if (a0 <= 9.f) {
	if (v0 == 0.f) {
	    z__4.r = z1.r * .5, z__4.i = z1.i * .5;
	    z_log(&z__3, &z__4);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z__1.r = z__2.r - .5772156649015329, z__1.i = z__2.i;
	    ct.r = z__1.r, ct.i = z__1.i;
	    cs.r = 0., cs.i = 0.;
	    w0 = 0.;
	    cr.r = 1., cr.i = 0.;
	    for (k = 1; k <= 50; ++k) {
		w0 += 1. / k;
		z__3.r = cr.r * .25, z__3.i = cr.i * .25;
		i__1 = k * k;
		d__1 = (doublereal) i__1;
		z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
		z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * 
			z2.i + z__2.i * z2.r;
		cr.r = z__1.r, cr.i = z__1.i;
		z__2.r = w0 + ct.r, z__2.i = ct.i;
		z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * 
			z__2.i + cr.i * z__2.r;
		cp.r = z__1.r, cp.i = z__1.i;
		z__1.r = cs.r + cp.r, z__1.i = cs.i + cp.i;
		cs.r = z__1.r, cs.i = z__1.i;
		z_div(&z__1, &cp, &cs);
		if (k >= 10 && z_abs(&z__1) < 1e-15) {
		    goto L45;
		}
/* L40: */
	    }
L45:
	    z__1.r = ct.r + cs.r, z__1.i = ct.i + cs.i;
	    cbk0.r = z__1.r, cbk0.i = z__1.i;
	} else {
	    v0n = 1. - v0;
	    gamma2_(&v0n, &gan);
	    z__4.r = z1.r * .5, z__4.i = z1.i * .5;
	    z__5.r = v0, z__5.i = 0.;
	    pow_zz(&z__3, &z__4, &z__5);
	    z__2.r = gan * z__3.r, z__2.i = gan * z__3.i;
	    z_div(&z__1, &c_b6, &z__2);
	    ca2.r = z__1.r, ca2.i = z__1.i;
	    z__3.r = z1.r * .5, z__3.i = z1.i * .5;
	    z__4.r = v0, z__4.i = 0.;
	    pow_zz(&z__2, &z__3, &z__4);
	    z__1.r = z__2.r / gap, z__1.i = z__2.i / gap;
	    ca1.r = z__1.r, ca1.i = z__1.i;
	    z__1.r = ca2.r - ca1.r, z__1.i = ca2.i - ca1.i;
	    csu.r = z__1.r, csu.i = z__1.i;
	    cr1.r = 1., cr1.i = 0.;
	    cr2.r = 1., cr2.i = 0.;
	    ws0 = 0.;
	    for (k = 1; k <= 50; ++k) {
		z__3.r = cr1.r * .25, z__3.i = cr1.i * .25;
		z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			z2.i + z__3.i * z2.r;
		d__1 = k * (k - v0);
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr1.r = z__1.r, cr1.i = z__1.i;
		z__3.r = cr2.r * .25, z__3.i = cr2.i * .25;
		z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * 
			z2.i + z__3.i * z2.r;
		d__1 = k * (k + v0);
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
		cr2.r = z__1.r, cr2.i = z__1.i;
		z__3.r = ca2.r * cr1.r - ca2.i * cr1.i, z__3.i = ca2.r * 
			cr1.i + ca2.i * cr1.r;
		z__2.r = csu.r + z__3.r, z__2.i = csu.i + z__3.i;
		z__4.r = ca1.r * cr2.r - ca1.i * cr2.i, z__4.i = ca1.r * 
			cr2.i + ca1.i * cr2.r;
		z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
		csu.r = z__1.r, csu.i = z__1.i;
		ws = z_abs(&csu);
		if (k >= 10 && (d__1 = ws - ws0, abs(d__1)) / ws < 1e-15) {
		    goto L55;
		}
		ws0 = ws;
/* L50: */
	    }
L55:
	    d__1 = pi * .5;
	    z__2.r = d__1 * csu.r, z__2.i = d__1 * csu.i;
	    d__2 = sin(piv);
	    z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	    cbk0.r = z__1.r, cbk0.i = z__1.i;
	}
    } else {
	z__3.r = -z1.r, z__3.i = -z1.i;
	z_exp(&z__2, &z__3);
	d__1 = pi * .5;
	z__6.r = d__1, z__6.i = 0.;
	z_div(&z__5, &z__6, &z1);
	z_sqrt(&z__4, &z__5);
	z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * z__4.i 
		+ z__2.i * z__4.r;
	cb.r = z__1.r, cb.i = z__1.i;
	cs.r = 1., cs.i = 0.;
	cr.r = 1., cr.i = 0.;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
	    z__3.r = cr.r * .125, z__3.i = cr.i * .125;
	    d__2 = k * 2. - 1.;
	    d__1 = vt - pow_dd(&d__2, &c_b4);
	    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	    d__3 = (doublereal) k;
	    z__4.r = d__3 * z1.r, z__4.i = d__3 * z1.i;
	    z_div(&z__1, &z__2, &z__4);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L60: */
	    z__1.r = cs.r + cr.r, z__1.i = cs.i + cr.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	}
	z__1.r = cb.r * cs.r - cb.i * cs.i, z__1.i = cb.r * cs.i + cb.i * 
		cs.r;
	cbk0.r = z__1.r, cbk0.i = z__1.i;
    }
    z_div(&z__3, &c_b6, &z1);
    z__4.r = cbi[1].r * cbk0.r - cbi[1].i * cbk0.i, z__4.i = cbi[1].r * 
	    cbk0.i + cbi[1].i * cbk0.r;
    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
    z_div(&z__1, &z__2, cbi);
    cbk1.r = z__1.r, cbk1.i = z__1.i;
    cbk[0].r = cbk0.r, cbk[0].i = cbk0.i;
    cbk[1].r = cbk1.r, cbk[1].i = cbk1.i;
    cg0.r = cbk0.r, cg0.i = cbk0.i;
    cg1.r = cbk1.r, cg1.i = cbk1.i;
    i__1 = n;
    for (k = 2; k <= i__1; ++k) {
	d__1 = (v0 + k - 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, &z1);
	z__2.r = z__3.r * cg1.r - z__3.i * cg1.i, z__2.i = z__3.r * cg1.i + 
		z__3.i * cg1.r;
	z__1.r = z__2.r + cg0.r, z__1.i = z__2.i + cg0.i;
	cgk.r = z__1.r, cgk.i = z__1.i;
	i__3 = k;
	cbk[i__3].r = cgk.r, cbk[i__3].i = cgk.i;
	cg0.r = cg1.r, cg0.i = cg1.i;
/* L65: */
	cg1.r = cgk.r, cg1.i = cgk.i;
    }
    if (z__->r < 0.f) {
	i__1 = n;
	for (k = 0; k <= i__1; ++k) {
	    d__1 = (k + v0) * pi;
	    z__2.r = d__1 * ci.r, z__2.i = d__1 * ci.i;
	    z_exp(&z__1, &z__2);
	    cvk.r = z__1.r, cvk.i = z__1.i;
	    if (d_imag(z__) < 0.) {
		i__3 = k;
		i__2 = k;
		z__2.r = cvk.r * cbk[i__2].r - cvk.i * cbk[i__2].i, z__2.i = 
			cvk.r * cbk[i__2].i + cvk.i * cbk[i__2].r;
		z__4.r = pi * ci.r, z__4.i = pi * ci.i;
		i__4 = k;
		z__3.r = z__4.r * cbi[i__4].r - z__4.i * cbi[i__4].i, z__3.i =
			 z__4.r * cbi[i__4].i + z__4.i * cbi[i__4].r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		cbk[i__3].r = z__1.r, cbk[i__3].i = z__1.i;
		i__3 = k;
		z_div(&z__1, &cbi[k], &cvk);
		cbi[i__3].r = z__1.r, cbi[i__3].i = z__1.i;
	    } else if (d_imag(z__) > 0.f) {
		i__3 = k;
		z_div(&z__2, &cbk[k], &cvk);
		z__4.r = pi * ci.r, z__4.i = pi * ci.i;
		i__2 = k;
		z__3.r = z__4.r * cbi[i__2].r - z__4.i * cbi[i__2].i, z__3.i =
			 z__4.r * cbi[i__2].i + z__4.i * cbi[i__2].r;
		z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		cbk[i__3].r = z__1.r, cbk[i__3].i = z__1.i;
		i__3 = k;
		i__2 = k;
		z__1.r = cvk.r * cbi[i__2].r - cvk.i * cbi[i__2].i, z__1.i = 
			cvk.r * cbi[i__2].i + cvk.i * cbi[i__2].r;
		cbi[i__3].r = z__1.r, cbi[i__3].i = z__1.i;
	    }
/* L70: */
	}
    }
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbi[0].r - z__3.i * cbi[0].i, z__2.i = z__3.r * cbi[0]
	    .i + z__3.i * cbi[0].r;
    z__1.r = z__2.r + cbi[1].r, z__1.i = z__2.i + cbi[1].i;
    cdi[0].r = z__1.r, cdi[0].i = z__1.i;
    z__4.r = v0, z__4.i = 0.;
    z_div(&z__3, &z__4, z__);
    z__2.r = z__3.r * cbk[0].r - z__3.i * cbk[0].i, z__2.i = z__3.r * cbk[0]
	    .i + z__3.i * cbk[0].r;
    z__1.r = z__2.r - cbk[1].r, z__1.i = z__2.i - cbk[1].i;
    cdk[0].r = z__1.r, cdk[0].i = z__1.i;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__3 = k;
	d__1 = -(k + v0);
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = k;
	z__2.r = z__3.r * cbi[i__2].r - z__3.i * cbi[i__2].i, z__2.i = z__3.r 
		* cbi[i__2].i + z__3.i * cbi[i__2].r;
	i__4 = k - 1;
	z__1.r = z__2.r + cbi[i__4].r, z__1.i = z__2.i + cbi[i__4].i;
	cdi[i__3].r = z__1.r, cdi[i__3].i = z__1.i;
/* L75: */
	i__3 = k;
	d__1 = -(k + v0);
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__2 = k;
	z__2.r = z__3.r * cbk[i__2].r - z__3.i * cbk[i__2].i, z__2.i = z__3.r 
		* cbk[i__2].i + z__3.i * cbk[i__2].r;
	i__4 = k - 1;
	z__1.r = z__2.r - cbk[i__4].r, z__1.i = z__2.i - cbk[i__4].i;
	cdk[i__3].r = z__1.r, cdk[i__3].i = z__1.i;
    }
    *vm = n + v0;
    return 0;
} /* cikva_ */

/*       ********************************** */
/* Subroutine */ int cfc_(doublecomplex *z__, doublecomplex *zf, 
	doublecomplex *zd)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_sqrt(
	    doublecomplex *, doublecomplex *), z_sin(doublecomplex *, 
	    doublecomplex *), z_cos(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex c__;
    static integer k, m;
    static doublereal w0;
    static doublecomplex z0, cf, cg, cr;
    static doublereal wa, pi;
    static doublecomplex zp, cf0, cf1;
    static doublereal wa0;
    static doublecomplex zp2;
    static doublereal eps;


/*       ========================================================= */
/*       Purpose: Compute complex Fresnel integral C(z) and C'(z) */
/*       Input :  z --- Argument of C(z) */
/*       Output:  ZF --- C(z) */
/*                ZD --- C'(z) */
/*       ========================================================= */

    eps = 1e-14;
    pi = 3.141592653589793;
    w0 = z_abs(z__);
    d__1 = pi * .5;
    z__2.r = d__1 * z__->r, z__2.i = d__1 * z__->i;
    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * z__->i + 
	    z__2.i * z__->r;
    zp.r = z__1.r, zp.i = z__1.i;
    z__1.r = zp.r * zp.r - zp.i * zp.i, z__1.i = zp.r * zp.i + zp.i * zp.r;
    zp2.r = z__1.r, zp2.i = z__1.i;
    z0.r = 0., z0.i = 0.;
    if (z__->r == z0.r && z__->i == z0.i) {
	c__.r = z0.r, c__.i = z0.i;
    } else if (w0 <= 2.5f) {
	cr.r = z__->r, cr.i = z__->i;
	c__.r = cr.r, c__.i = cr.i;
	wa0 = 0.;
	for (k = 1; k <= 80; ++k) {
	    z__6.r = cr.r * -.5, z__6.i = cr.i * -.5;
	    d__1 = k * 4. - 3.;
	    z__5.r = d__1 * z__6.r, z__5.i = d__1 * z__6.i;
	    d__2 = (doublereal) k;
	    z__4.r = z__5.r / d__2, z__4.i = z__5.i / d__2;
	    d__3 = k * 2. - 1.;
	    z__3.r = z__4.r / d__3, z__3.i = z__4.i / d__3;
	    d__4 = k * 4. + 1.;
	    z__2.r = z__3.r / d__4, z__2.i = z__3.i / d__4;
	    z__1.r = z__2.r * zp2.r - z__2.i * zp2.i, z__1.i = z__2.r * zp2.i 
		    + z__2.i * zp2.r;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = c__.r + cr.r, z__1.i = c__.i + cr.i;
	    c__.r = z__1.r, c__.i = z__1.i;
	    wa = z_abs(&c__);
	    if ((d__1 = (wa - wa0) / wa, abs(d__1)) < eps && k > 10) {
		goto L30;
	    }
/* L10: */
	    wa0 = wa;
	}
    } else if (w0 > 2.5f && w0 < 4.5f) {
	m = 85;
	c__.r = z0.r, c__.i = z0.i;
	cf1.r = z0.r, cf1.i = z0.i;
	cf0.r = 1e-100, cf0.i = 0.;
	for (k = m; k >= 0; --k) {
	    d__1 = k * 2. + 3.;
	    z__3.r = d__1 * cf0.r, z__3.i = d__1 * cf0.i;
	    z_div(&z__2, &z__3, &zp);
	    z__1.r = z__2.r - cf1.r, z__1.i = z__2.i - cf1.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	    if (k == k / 2 << 1) {
		z__1.r = c__.r + cf.r, z__1.i = c__.i + cf.i;
		c__.r = z__1.r, c__.i = z__1.i;
	    }
	    cf1.r = cf0.r, cf1.i = cf0.i;
/* L15: */
	    cf0.r = cf.r, cf0.i = cf.i;
	}
	z__6.r = pi * zp.r, z__6.i = pi * zp.i;
	z_div(&z__5, &c_b15, &z__6);
	z_sqrt(&z__4, &z__5);
	z_sin(&z__7, &zp);
	z__3.r = z__4.r * z__7.r - z__4.i * z__7.i, z__3.i = z__4.r * z__7.i 
		+ z__4.i * z__7.r;
	z_div(&z__2, &z__3, &cf);
	z__1.r = z__2.r * c__.r - z__2.i * c__.i, z__1.i = z__2.r * c__.i + 
		z__2.i * c__.r;
	c__.r = z__1.r, c__.i = z__1.i;
    } else {
	cr.r = 1., cr.i = 0.;
	cf.r = 1., cf.i = 0.;
	for (k = 1; k <= 20; ++k) {
	    z__4.r = cr.r * -.25, z__4.i = cr.i * -.25;
	    d__1 = k * 4. - 1.;
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    d__2 = k * 4. - 3.;
	    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
	    z_div(&z__1, &z__2, &zp2);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L20: */
	    z__1.r = cf.r + cr.r, z__1.i = cf.i + cr.i;
	    cf.r = z__1.r, cf.i = z__1.i;
	}
	z__3.r = pi * z__->r, z__3.i = pi * z__->i;
	z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i 
		+ z__3.i * z__->r;
	z_div(&z__1, &c_b6, &z__2);
	cr.r = z__1.r, cr.i = z__1.i;
	cg.r = cr.r, cg.i = cr.i;
	for (k = 1; k <= 12; ++k) {
	    z__4.r = cr.r * -.25, z__4.i = cr.i * -.25;
	    d__1 = k * 4. + 1.;
	    z__3.r = d__1 * z__4.r, z__3.i = d__1 * z__4.i;
	    d__2 = k * 4. - 1.;
	    z__2.r = d__2 * z__3.r, z__2.i = d__2 * z__3.i;
	    z_div(&z__1, &z__2, &zp2);
	    cr.r = z__1.r, cr.i = z__1.i;
/* L25: */
	    z__1.r = cg.r + cr.r, z__1.i = cg.i + cr.i;
	    cg.r = z__1.r, cg.i = z__1.i;
	}
	z_sin(&z__5, &zp);
	z__4.r = cf.r * z__5.r - cf.i * z__5.i, z__4.i = cf.r * z__5.i + cf.i 
		* z__5.r;
	z_cos(&z__7, &zp);
	z__6.r = cg.r * z__7.r - cg.i * z__7.i, z__6.i = cg.r * z__7.i + cg.i 
		* z__7.r;
	z__3.r = z__4.r - z__6.r, z__3.i = z__4.i - z__6.i;
	z__8.r = pi * z__->r, z__8.i = pi * z__->i;
	z_div(&z__2, &z__3, &z__8);
	z__1.r = z__2.r + .5, z__1.i = z__2.i;
	c__.r = z__1.r, c__.i = z__1.i;
    }
L30:
    zf->r = c__.r, zf->i = c__.i;
    d__1 = pi * .5f;
    z__3.r = d__1 * z__->r, z__3.i = d__1 * z__->i;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z_cos(&z__1, &z__2);
    zd->r = z__1.r, zd->i = z__1.i;
    return 0;
} /* cfc_ */

/*       ********************************** */
/* Subroutine */ int fcs_(doublereal *x, doublereal *c__, doublereal *s)
{
    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal f, g;
    static integer k, m;
    static doublereal q, r__, t, f0, f1, t0, t2, pi, xa, px, su, eps;


/*       ================================================= */
/*       Purpose: Compute Fresnel integrals C(x) and S(x) */
/*       Input :  x --- Argument of C(x) and S(x) */
/*       Output:  C --- C(x) */
/*                S --- S(x) */
/*       ================================================= */

    eps = 1e-15;
    pi = 3.141592653589793;
    xa = abs(*x);
    px = pi * xa;
    t = px * .5 * xa;
    t2 = t * t;
    if (xa == 0.f) {
	*c__ = 0.;
	*s = 0.;
    } else if (xa < 2.5) {
	r__ = xa;
	*c__ = r__;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * -.5 * (k * 4. - 3.) / k / (k * 2. - 1.) / (k * 4. + 
		    1.) * t2;
	    *c__ += r__;
	    if (abs(r__) < abs(*c__) * eps) {
		goto L15;
	    }
/* L10: */
	}
L15:
	*s = xa * t / 3.;
	r__ = *s;
	for (k = 1; k <= 50; ++k) {
	    r__ = r__ * -.5 * (k * 4. - 1.) / k / (k * 2. + 1.) / (k * 4. + 
		    3.) * t2;
	    *s += r__;
	    if (abs(r__) < abs(*s) * eps) {
		goto L40;
	    }
/* L20: */
	}
    } else if (xa < 4.5) {
	m = (integer) (t * 1.75f + 42.f);
	su = 0.;
	*c__ = 0.;
	*s = 0.;
	f1 = 0.;
	f0 = 1e-100;
	for (k = m; k >= 0; --k) {
	    f = (k * 2. + 3.) * f0 / t - f1;
	    if (k == k / 2 << 1) {
		*c__ += f;
	    } else {
		*s += f;
	    }
	    su += (k * 2. + 1.) * f * f;
	    f1 = f0;
/* L25: */
	    f0 = f;
	}
	q = sqrt(su);
	*c__ = *c__ * xa / q;
	*s = *s * xa / q;
    } else {
	r__ = 1.;
	f = 1.;
	for (k = 1; k <= 20; ++k) {
	    r__ = r__ * -.25 * (k * 4. - 1.) * (k * 4. - 3.) / t2;
/* L30: */
	    f += r__;
	}
	r__ = 1. / (px * xa);
	g = r__;
	for (k = 1; k <= 12; ++k) {
	    r__ = r__ * -.25 * (k * 4. + 1.) * (k * 4. - 1.) / t2;
/* L35: */
	    g += r__;
	}
	t0 = t - (integer) (t / (pi * 2.)) * 2. * pi;
	*c__ = (f * sin(t0) - g * cos(t0)) / px + .5;
	*s = .5 - (f * cos(t0) + g * sin(t0)) / px;
    }
L40:
    if (*x < 0.) {
	*c__ = -(*c__);
	*s = -(*s);
    }
    return 0;
} /* fcs_ */

/*       ********************************** */
/* Subroutine */ int rctj_(integer *n, doublereal *x, integer *nm, doublereal 
	*rj, doublereal *dj)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, m;
    static doublereal f0, f1, cs, rj0, rj1;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ======================================================== */
/*       Purpose: Compute Riccati-Bessel functions of the first */
/*                kind and their derivatives */
/*       Input:   x --- Argument of Riccati-Bessel function */
/*                n --- Order of jn(x)  ( n = 0,1,2,... ) */
/*       Output:  RJ(n) --- x·jn(x) */
/*                DJ(n) --- [x·jn(x)]' */
/*                NM --- Highest order computed */
/*       Routines called: */
/*                MSTA1 and MSTA2 for computing the starting */
/*                point for backward recurrence */
/*       ======================================================== */

    *nm = *n;
    if (abs(*x) < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    rj[k] = 0.;
/* L10: */
	    dj[k] = 0.;
	}
	dj[0] = 1.;
	return 0;
    }
    rj[0] = sin(*x);
    rj[1] = rj[0] / *x - cos(*x);
    rj0 = rj[0];
    rj1 = rj[1];
    cs = 0.;
    f = 0.;
    if (*n >= 2) {
	m = msta1_(x, &c__200);
	if (m < *n) {
	    *nm = m;
	} else {
	    m = msta2_(x, n, &c__15);
	}
	f0 = 0.;
	f1 = 1e-100;
	for (k = m; k >= 0; --k) {
	    f = (k * 2. + 3.) * f1 / *x - f0;
	    if (k <= *nm) {
		rj[k] = f;
	    }
	    f0 = f1;
/* L15: */
	    f1 = f;
	}
	if (abs(rj0) > abs(rj1)) {
	    cs = rj0 / f;
	}
	if (abs(rj0) <= abs(rj1)) {
	    cs = rj1 / f0;
	}
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
/* L20: */
	    rj[k] = cs * rj[k];
	}
    }
    dj[0] = cos(*x);
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	dj[k] = -k * rj[k] / *x + rj[k - 1];
    }
    return 0;
} /* rctj_ */

/*       ********************************** */
/* Subroutine */ int herzo_(integer *n, doublereal *x, doublereal *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal p, q, r__, z__, f0, f1, r1, r2, z0, fd, gd, hd, hf, hn;
    static integer it, nr;
    static doublereal zl, wp;


/*       ======================================================== */
/*       Purpose : Compute the zeros of Hermite polynomial Ln(x) */
/*                 in the interval [-∞,∞], and the corresponding */
/*                 weighting coefficients for Gauss-Hermite */
/*                 integration */
/*       Input :   n    --- Order of the Hermite polynomial */
/*                 X(n) --- Zeros of the Hermite polynomial */
/*                 W(n) --- Corresponding weighting coefficients */
/*       ======================================================== */

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    hn = 1. / *n;
    d__1 = (doublereal) (*n);
    zl = pow_dd(&d__1, &c_b50) * 1.46 - 1.1611;
    z__ = 0.;
    hf = 0.;
    hd = 0.;
    i__1 = *n / 2;
    for (nr = 1; nr <= i__1; ++nr) {
	if (nr == 1) {
	    z__ = zl;
	}
	if (nr != 1) {
	    z__ -= hn * (*n / 2 + 1 - nr);
	}
	it = 0;
L10:
	++it;
	z0 = z__;
	f0 = 1.;
	f1 = z__ * 2.;
	i__2 = *n;
	for (k = 2; k <= i__2; ++k) {
	    hf = z__ * 2. * f1 - (k - 1.) * 2. * f0;
	    hd = k * 2. * f1;
	    f0 = f1;
/* L15: */
	    f1 = hf;
	}
	p = 1.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    p *= z__ - x[i__];
	}
	fd = hf / p;
	q = 0.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wp = 1.;
	    i__3 = nr - 1;
	    for (j = 1; j <= i__3; ++j) {
		if (j == i__) {
		    goto L25;
		}
		wp *= z__ - x[j];
L25:
		;
	    }
/* L30: */
	    q += wp;
	}
	gd = (hd - q * fd) / p;
	z__ -= fd / gd;
	if (it <= 40 && (d__1 = (z__ - z0) / z__, abs(d__1)) > 1e-15) {
	    goto L10;
	}
	x[nr] = z__;
	x[*n + 1 - nr] = -z__;
	r__ = 1.;
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
/* L35: */
	    r__ = r__ * 2. * k;
	}
	w[nr] = r__ * 3.544907701811 / (hd * hd);
/* L40: */
	w[*n + 1 - nr] = w[nr];
    }
    if (*n != *n / 2 << 1) {
	r1 = 1.;
	r2 = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    r1 = r1 * 2. * j;
	    if (j >= (*n + 1) / 2) {
		r2 *= j;
	    }
/* L45: */
	}
	w[*n / 2 + 1] = r1 * .88622692545276 / (r2 * r2);
	x[*n / 2 + 1] = 0.;
    }
    return 0;
} /* herzo_ */

/*       ********************************** */
/* Subroutine */ int jy01b_(doublereal *x, doublereal *bj0, doublereal *dj0, 
	doublereal *bj1, doublereal *dj1, doublereal *by0, doublereal *dy0, 
	doublereal *by1, doublereal *dy1)
{
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static doublereal t, a0, p0, q0, p1, q1, t2, pi, ta0, ta1;


/*       ======================================================= */
/*       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x), */
/*                Y1(x), and their derivatives */
/*       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ≥ 0 ) */
/*       Output:  BJ0 --- J0(x) */
/*                DJ0 --- J0'(x) */
/*                BJ1 --- J1(x) */
/*                DJ1 --- J1'(x) */
/*                BY0 --- Y0(x) */
/*                DY0 --- Y0'(x) */
/*                BY1 --- Y1(x) */
/*                DY1 --- Y1'(x) */
/*       ======================================================= */

    pi = 3.141592653589793;
    if (*x == 0.) {
	*bj0 = 1.;
	*bj1 = 0.;
	*dj0 = 0.;
	*dj1 = .5;
	*by0 = -1e300;
	*by1 = -1e300;
	*dy0 = 1e300;
	*dy1 = 1e300;
	return 0;
    } else if (*x <= 4.) {
	t = *x / 4.;
	t2 = t * t;
	*bj0 = ((((((t2 * -5.014415e-4 + .0076771853) * t2 - .0709253492) * 
		t2 + .4443584263) * t2 - 1.7777560599) * t2 + 3.9999973021) * 
		t2 - 3.9999998721) * t2 + 1.;
	*bj1 = t * (((((((t2 * -1.289769e-4 + .0022069155) * t2 - .0236616773)
		 * t2 + .1777582922) * t2 - .8888839649) * t2 + 2.6666660544) 
		* t2 - 3.999999971) * t2 + 1.9999999998);
	*by0 = (((((((t2 * -5.67433e-5 + 8.59977e-4) * t2 - .0094855882) * t2 
		+ .0772975809) * t2 - .4261737419) * t2 + 1.4216421221) * t2 
		- 2.3498519931) * t2 + 1.0766115157) * t2 + .3674669052;
	*by0 = 2. / pi * log(*x / 2.) * *bj0 + *by0;
	*by1 = ((((((((t2 * 6.535773e-4 - .0108175626) * t2 + .107657606) * 
		t2 - .7268945577) * t2 + 3.1261399273) * t2 - 7.3980241381) * 
		t2 + 6.8529236342) * t2 + .3932562018) * t2 - .6366197726) / *
		x;
	*by1 = 2. / pi * log(*x / 2.) * *bj1 + *by1;
    } else {
	t = 4. / *x;
	t2 = t * t;
	a0 = sqrt(2. / (pi * *x));
	p0 = ((((t2 * -9.285e-6 + 4.3506e-5) * t2 - 1.22226e-4) * t2 + 
		4.34725e-4) * t2 - .004394275) * t2 + .999999997;
	q0 = t * (((((t2 * 8.099e-6 - 3.5614e-5) * t2 + 8.5844e-5) * t2 - 
		2.18024e-4) * t2 + .001144106) * t2 - .031249995);
	ta0 = *x - pi * .25;
	*bj0 = a0 * (p0 * cos(ta0) - q0 * sin(ta0));
	*by0 = a0 * (p0 * sin(ta0) + q0 * cos(ta0));
	p1 = ((((t2 * 1.0632e-5 - 5.0363e-5) * t2 + 1.45575e-4) * t2 - 
		5.59487e-4) * t2 + .007323931) * t2 + 1.000000004;
	q1 = t * (((((t2 * -9.173e-6 + 4.0658e-5) * t2 - 9.9941e-5) * t2 + 
		2.66891e-4) * t2 - .001601836) * t2 + .093749994);
	ta1 = *x - pi * .75;
	*bj1 = a0 * (p1 * cos(ta1) - q1 * sin(ta1));
	*by1 = a0 * (p1 * sin(ta1) + q1 * cos(ta1));
    }
    *dj0 = -(*bj1);
    *dj1 = *bj0 - *bj1 / *x;
    *dy0 = -(*by1);
    *dy1 = *by0 - *by1 / *x;
    return 0;
} /* jy01b_ */

/*       ********************************** */
/* Subroutine */ int enxb_(integer *n, doublereal *x, doublereal *en)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), log(doublereal);

    /* Local variables */
    static integer j, k, l, m;
    static doublereal r__, s, t, s0, t0, rp, ps, ens;


/*       =============================================== */
/*       Purpose: Compute exponential integral En(x) */
/*       Input :  x --- Argument of En(x) */
/*                n --- Order of En(x)  (n = 0,1,2,...) */
/*       Output:  EN(n) --- En(x) */
/*       =============================================== */

    if (*x == 0.f) {
	en[0] = 1e300;
	en[1] = 1e300;
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
/* L10: */
	    en[k] = 1. / (k - 1.f);
	}
	return 0;
    } else if (*x <= 1.f) {
	en[0] = exp(-(*x)) / *x;
	s0 = 0.;
	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    rp = 1.;
	    i__2 = l - 1;
	    for (j = 1; j <= i__2; ++j) {
/* L15: */
		rp = -rp * *x / j;
	    }
	    ps = -.5772156649015328;
	    i__2 = l - 1;
	    for (m = 1; m <= i__2; ++m) {
/* L20: */
		ps += 1. / m;
	    }
	    ens = rp * (-log(*x) + ps);
	    s = 0.;
	    for (m = 0; m <= 20; ++m) {
		if (m == l - 1) {
		    goto L30;
		}
		r__ = 1.;
		i__2 = m;
		for (j = 1; j <= i__2; ++j) {
/* L25: */
		    r__ = -r__ * *x / j;
		}
		s += r__ / (m - l + 1.);
		if ((d__1 = s - s0, abs(d__1)) < abs(s) * 1e-15) {
		    goto L35;
		}
		s0 = s;
L30:
		;
	    }
L35:
	    en[l] = ens - s;
/* L40: */
	}
    } else {
	en[0] = exp(-(*x)) / *x;
	m = (integer) (100.f / *x) + 15;
	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    t0 = 0.;
	    for (k = m; k >= 1; --k) {
/* L45: */
		t0 = (l + k - 1.) / (k / (*x + t0) + 1.);
	    }
	    t = 1. / (*x + t0);
/* L50: */
	    en[l] = exp(-(*x)) * t;
	}
    }
    return 0;
} /* enxb_ */

/*       ********************************** */
/* Subroutine */ int sphk_(integer *n, doublereal *x, integer *nm, doublereal 
	*sk, doublereal *dk)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k;
    static doublereal f0, f1, pi;


/*       ===================================================== */
/*       Purpose: Compute modified spherical Bessel functions */
/*                of the second kind, kn(x) and kn'(x) */
/*       Input :  x --- Argument of kn(x)  ( x ≥ 0 ) */
/*                n --- Order of kn(x) ( n = 0,1,2,... ) */
/*       Output:  SK(n) --- kn(x) */
/*                DK(n) --- kn'(x) */
/*                NM --- Highest order computed */
/*       ===================================================== */

    pi = 3.141592653589793;
    *nm = *n;
    if (*x < 1e-60) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    sk[k] = 1e300;
/* L10: */
	    dk[k] = -1e300;
	}
	return 0;
    }
    sk[0] = pi * .5 / *x * exp(-(*x));
    sk[1] = sk[0] * (1. / *x + 1.);
    f0 = sk[0];
    f1 = sk[1];
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	f = (k * 2. - 1.) * f1 / *x + f0;
	sk[k] = f;
	if (abs(f) > 1e300) {
	    goto L20;
	}
	f0 = f1;
/* L15: */
	f1 = f;
    }
L20:
    *nm = k - 1;
    dk[0] = -sk[1];
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	dk[k] = -sk[k - 1] - (k + 1.) / *x * sk[k];
    }
    return 0;
} /* sphk_ */

/*       ********************************** */
/* Subroutine */ int enxa_(integer *n, doublereal *x, doublereal *en)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer k;
    static doublereal e1, ek;
    extern /* Subroutine */ int e1xb_(doublereal *, doublereal *);


/*       ============================================ */
/*       Purpose: Compute exponential integral En(x) */
/*       Input :  x --- Argument of En(x) ( x ≤ 20 ) */
/*                n --- Order of En(x) */
/*       Output:  EN(n) --- En(x) */
/*       Routine called: E1XB for computing E1(x) */
/*       ============================================ */

    en[0] = exp(-(*x)) / *x;
    e1xb_(x, &e1);
    en[1] = e1;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	ek = (exp(-(*x)) - *x * e1) / (k - 1.);
	en[k] = ek;
/* L10: */
	e1 = ek;
    }
    return 0;
} /* enxa_ */

/*       ********************************** */
/* Subroutine */ int gaih_(doublereal *x, doublereal *ga)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, m1;
    static doublereal pi;


/*       ===================================================== */
/*       Purpose: Compute gamma function Г(x) */
/*       Input :  x  --- Argument of Г(x), x = n/2, n=1,2,… */
/*       Output:  GA --- Г(x) */
/*       ===================================================== */

    pi = 3.141592653589793;
    if (*x == (doublereal) ((integer) (*x)) && *x > 0.f) {
	*ga = 1.;
	m1 = (integer) (*x - 1.f);
	i__1 = m1;
	for (k = 2; k <= i__1; ++k) {
/* L10: */
	    *ga *= k;
	}
    } else if (*x + .5 == (doublereal) ((integer) (*x + .5)) && *x > 0.f) {
	m = (integer) (*x);
	*ga = sqrt(pi);
	i__1 = m;
	for (k = 1; k <= i__1; ++k) {
/* L15: */
	    *ga = *ga * .5 * (k * 2. - 1.);
	}
    }
    return 0;
} /* gaih_ */

/*       ********************************** */
/* Subroutine */ int pbvv_(doublereal *v, doublereal *x, doublereal *vv, 
	doublereal *vp, doublereal *pvf, doublereal *pvd)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), exp(doublereal), sqrt(
	    doublereal);

    /* Local variables */
    static doublereal f;
    static integer k, l, m;
    static doublereal f0, f1, s0, v0, v1, v2;
    static integer ja, na;
    static doublereal qe, pi, xa, vh;
    static integer kv, nv;
    static doublereal q2p, pv0;
    extern /* Subroutine */ int vvla_(doublereal *, doublereal *, doublereal *
	    ), vvsa_(doublereal *, doublereal *, doublereal *);


/*       =================================================== */
/*       Purpose: Compute parabolic cylinder functions Vv(x) */
/*                and their derivatives */
/*       Input:   x --- Argument of Vv(x) */
/*                v --- Order of Vv(x) */
/*       Output:  VV(na) --- Vv(x) */
/*                VP(na) --- Vv'(x) */
/*                ( na = |n|, v = n+v0, |v0| < 1 */
/*                  n = 0,±1,±2,… ) */
/*                PVF --- Vv(x) */
/*                PVD --- Vv'(x) */
/*       Routines called: */
/*             (1) VVSA for computing Vv(x) for small |x| */
/*             (2) VVLA for computing Vv(x) for large |x| */
/*       =================================================== */

    pi = 3.141592653589793;
    xa = abs(*x);
    vh = *v;
    *v += d_sign(&c_b524, v);
    nv = (integer) (*v);
    v0 = *v - nv;
    na = abs(nv);
    qe = exp(*x * .25 * *x);
    q2p = sqrt(2. / pi);
    ja = 0;
    if (na >= 1) {
	ja = 1;
    }
    f = 0.;
    if (*v <= 0.f) {
	if (v0 == 0.f) {
	    if (xa <= 7.5f) {
		vvsa_(&v0, x, &pv0);
	    }
	    if (xa > 7.5f) {
		vvla_(&v0, x, &pv0);
	    }
	    f0 = q2p * qe;
	    f1 = *x * f0;
	    vv[0] = pv0;
	    vv[1] = f0;
	    vv[2] = f1;
	} else {
	    i__1 = ja;
	    for (l = 0; l <= i__1; ++l) {
		v1 = v0 - l;
		if (xa <= 7.5f) {
		    vvsa_(&v1, x, &f1);
		}
		if (xa > 7.5f) {
		    vvla_(&v1, x, &f1);
		}
		if (l == 0) {
		    f0 = f1;
		}
/* L10: */
	    }
	    vv[0] = f0;
	    vv[1] = f1;
	}
	kv = 2;
	if (v0 == 0.f) {
	    kv = 3;
	}
	i__1 = na;
	for (k = kv; k <= i__1; ++k) {
	    f = *x * f1 + (k - v0 - 2.) * f0;
	    vv[k] = f;
	    f0 = f1;
/* L15: */
	    f1 = f;
	}
    } else {
	if (*x >= 0.f && *x <= 7.5) {
	    v2 = *v;
	    if (v2 < 1.f) {
		v2 += 1.;
	    }
	    vvsa_(&v2, x, &f1);
	    v1 = v2 - 1.;
	    kv = (integer) v2;
	    vvsa_(&v1, x, &f0);
	    vv[kv] = f1;
	    vv[kv - 1] = f0;
	    for (k = kv - 2; k >= 0; --k) {
		f = *x * f0 - (k + v0 + 2.) * f1;
		if (k <= na) {
		    vv[k] = f;
		}
		f1 = f0;
/* L20: */
		f0 = f;
	    }
	} else if (*x > 7.5) {
	    vvla_(&v0, x, &pv0);
	    m = abs(na) + 100;
	    vv[1] = pv0;
	    f1 = 0.;
	    f0 = 1e-40;
	    for (k = m; k >= 0; --k) {
		f = *x * f0 - (k + v0 + 2.) * f1;
		if (k <= na) {
		    vv[k] = f;
		}
		f1 = f0;
/* L25: */
		f0 = f;
	    }
	    s0 = pv0 / f;
	    i__1 = na;
	    for (k = 0; k <= i__1; ++k) {
/* L30: */
		vv[k] = s0 * vv[k];
	    }
	} else {
	    if (xa <= 7.5) {
		vvsa_(&v0, x, &f0);
		v1 = v0 + 1.f;
		vvsa_(&v1, x, &f1);
	    } else {
		vvla_(&v0, x, &f0);
		v1 = v0 + 1.;
		vvla_(&v1, x, &f1);
	    }
	    vv[0] = f0;
	    vv[1] = f1;
	    i__1 = na;
	    for (k = 2; k <= i__1; ++k) {
		f = (*x * f1 - f0) / (k + v0);
		vv[k] = f;
		f0 = f1;
/* L35: */
		f1 = f;
	    }
	}
    }
    i__1 = na - 1;
    for (k = 0; k <= i__1; ++k) {
	v1 = v0 + k;
	if (*v >= 0.) {
	    vp[k] = *x * .5 * vv[k] - (v1 + 1.) * vv[k + 1];
	} else {
	    vp[k] = *x * -.5 * vv[k] + vv[k + 1];
	}
/* L40: */
    }
    *pvf = vv[na - 1];
    *pvd = vp[na - 1];
    *v = vh;
    return 0;
} /* pbvv_ */

/*       ********************************** */
/* Subroutine */ int clqmn_(integer *mm, integer *m, integer *n, doublereal *
	x, doublereal *y, doublecomplex *cqm, doublecomplex *cqd)
{
    /* System generated locals */
    integer cqm_dim1, cqm_offset, cqd_dim1, cqd_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);
    void z_sqrt(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *), z_log(doublecomplex *, 
	    doublecomplex *);
    double log(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex z__;
    static integer km;
    static doublereal xc;
    static integer ls;
    static doublecomplex zq, zs, cq0, cq1, cq10, cqf, cqf0, cqf1, cqf2;


/*       ======================================================= */
/*       Purpose: Compute the associated Legendre functions of */
/*                the second kind, Qmn(z) and Qmn'(z), for a */
/*                complex argument */
/*       Input :  x  --- Real part of z */
/*                y  --- Imaginary part of z */
/*                m  --- Order of Qmn(z)  ( m = 0,1,2,… ) */
/*                n  --- Degree of Qmn(z) ( n = 0,1,2,… ) */
/*                mm --- Physical dimension of CQM and CQD */
/*       Output:  CQM(m,n) --- Qmn(z) */
/*                CQD(m,n) --- Qmn'(z) */
/*       ======================================================= */

    /* Parameter adjustments */
    cqd_dim1 = *mm - 0 + 1;
    cqd_offset = 0 + cqd_dim1 * 0;
    cqd -= cqd_offset;
    cqm_dim1 = *mm - 0 + 1;
    cqm_offset = 0 + cqm_dim1 * 0;
    cqm -= cqm_offset;

    /* Function Body */
    z__1.r = *x, z__1.i = *y;
    z__.r = z__1.r, z__.i = z__1.i;
    if (abs(*x) == 1. && *y == 0.) {
	i__1 = *m;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = 0; j <= i__2; ++j) {
		i__3 = i__ + j * cqm_dim1;
		cqm[i__3].r = 1e300, cqm[i__3].i = 0.;
		i__3 = i__ + j * cqd_dim1;
		cqd[i__3].r = 1e300, cqd[i__3].i = 0.;
/* L10: */
	    }
	}
	return 0;
    }
    xc = z_abs(&z__);
    ls = 0;
    if (d_imag(&z__) == 0. || xc < 1.) {
	ls = 1;
    }
    if (xc > 1.) {
	ls = -1;
    }
    z__4.r = z__.r * z__.r - z__.i * z__.i, z__4.i = z__.r * z__.i + z__.i * 
	    z__.r;
    z__3.r = 1. - z__4.r, z__3.i = -z__4.i;
    d__1 = (doublereal) ls;
    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
    z_sqrt(&z__1, &z__2);
    zq.r = z__1.r, zq.i = z__1.i;
    z__3.r = z__.r * z__.r - z__.i * z__.i, z__3.i = z__.r * z__.i + z__.i * 
	    z__.r;
    z__2.r = 1. - z__3.r, z__2.i = -z__3.i;
    d__1 = (doublereal) ls;
    z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
    zs.r = z__1.r, zs.i = z__1.i;
    z__5.r = z__.r + 1., z__5.i = z__.i;
    d__1 = (doublereal) ls;
    z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
    z__6.r = 1. - z__.r, z__6.i = -z__.i;
    z_div(&z__3, &z__4, &z__6);
    z_log(&z__2, &z__3);
    z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
    cq0.r = z__1.r, cq0.i = z__1.i;
    if (xc < 1.0001) {
	cqm[0].r = cq0.r, cqm[0].i = cq0.i;
	i__2 = cqm_dim1;
	z__2.r = z__.r * cq0.r - z__.i * cq0.i, z__2.i = z__.r * cq0.i + 
		z__.i * cq0.r;
	z__1.r = z__2.r - 1., z__1.i = z__2.i;
	cqm[i__2].r = z__1.r, cqm[i__2].i = z__1.i;
	z_div(&z__1, &c_b1387, &zq);
	cqm[1].r = z__1.r, cqm[1].i = z__1.i;
	i__2 = cqm_dim1 + 1;
	z__2.r = -zq.r, z__2.i = -zq.i;
	z__6.r = z__.r * z__.r - z__.i * z__.i, z__6.i = z__.r * z__.i + 
		z__.i * z__.r;
	z__5.r = 1. - z__6.r, z__5.i = -z__6.i;
	z_div(&z__4, &z__, &z__5);
	z__3.r = cq0.r + z__4.r, z__3.i = cq0.i + z__4.i;
	z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * z__3.i 
		+ z__2.i * z__3.r;
	cqm[i__2].r = z__1.r, cqm[i__2].i = z__1.i;
	for (i__ = 0; i__ <= 1; ++i__) {
	    i__2 = *n;
	    for (j = 2; j <= i__2; ++j) {
		i__1 = i__ + j * cqm_dim1;
		d__1 = j * 2. - 1.;
		z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
		i__3 = i__ + (j - 1) * cqm_dim1;
		z__3.r = z__4.r * cqm[i__3].r - z__4.i * cqm[i__3].i, z__3.i =
			 z__4.r * cqm[i__3].i + z__4.i * cqm[i__3].r;
		d__2 = j + i__ - 1.;
		i__4 = i__ + (j - 2) * cqm_dim1;
		z__5.r = d__2 * cqm[i__4].r, z__5.i = d__2 * cqm[i__4].i;
		z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
		i__5 = j - i__;
		d__3 = (doublereal) i__5;
		z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
		cqm[i__1].r = z__1.r, cqm[i__1].i = z__1.i;
/* L15: */
	    }
	}
	i__2 = *n;
	for (j = 0; j <= i__2; ++j) {
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		i__3 = i__ + j * cqm_dim1;
		d__1 = (i__ - 1.) * -2.;
		z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
		z_div(&z__3, &z__4, &zq);
		i__4 = i__ - 1 + j * cqm_dim1;
		z__2.r = z__3.r * cqm[i__4].r - z__3.i * cqm[i__4].i, z__2.i =
			 z__3.r * cqm[i__4].i + z__3.i * cqm[i__4].r;
		d__2 = ls * (j + i__ - 1.) * (j - i__ + 2.);
		i__5 = i__ - 2 + j * cqm_dim1;
		z__5.r = d__2 * cqm[i__5].r, z__5.i = d__2 * cqm[i__5].i;
		z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - z__5.i;
		cqm[i__3].r = z__1.r, cqm[i__3].i = z__1.i;
/* L20: */
	    }
	}
    } else {
	if (xc > 1.1f) {
	    km = *m + 40 + *n;
	} else {
	    km = (*m + 40 + *n) * (integer) (-1.f - log(xc - 1.f) * 1.8f);
	}
	cqf2.r = 0., cqf2.i = 0.;
	cqf1.r = 1., cqf1.i = 0.;
	for (k = km; k >= 0; --k) {
	    d__1 = (k << 1) + 3.;
	    z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	    z__3.r = z__4.r * cqf1.r - z__4.i * cqf1.i, z__3.i = z__4.r * 
		    cqf1.i + z__4.i * cqf1.r;
	    d__2 = k + 2.;
	    z__5.r = d__2 * cqf2.r, z__5.i = d__2 * cqf2.i;
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	    d__3 = k + 1.;
	    z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
	    cqf0.r = z__1.r, cqf0.i = z__1.i;
	    if (k <= *n) {
		i__1 = k * cqm_dim1;
		cqm[i__1].r = cqf0.r, cqm[i__1].i = cqf0.i;
	    }
	    cqf2.r = cqf1.r, cqf2.i = cqf1.i;
/* L25: */
	    cqf1.r = cqf0.r, cqf1.i = cqf0.i;
	}
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
/* L30: */
	    i__2 = k * cqm_dim1;
	    i__3 = k * cqm_dim1;
	    z__2.r = cq0.r * cqm[i__3].r - cq0.i * cqm[i__3].i, z__2.i = 
		    cq0.r * cqm[i__3].i + cq0.i * cqm[i__3].r;
	    z_div(&z__1, &z__2, &cqf0);
	    cqm[i__2].r = z__1.r, cqm[i__2].i = z__1.i;
	}
	cqf2.r = 0., cqf2.i = 0.;
	cqf1.r = 1., cqf1.i = 0.;
	for (k = km; k >= 0; --k) {
	    d__1 = (k << 1) + 3.;
	    z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	    z__3.r = z__4.r * cqf1.r - z__4.i * cqf1.i, z__3.i = z__4.r * 
		    cqf1.i + z__4.i * cqf1.r;
	    d__2 = k + 1.;
	    z__5.r = d__2 * cqf2.r, z__5.i = d__2 * cqf2.i;
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	    d__3 = k + 2.;
	    z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
	    cqf0.r = z__1.r, cqf0.i = z__1.i;
	    if (k <= *n) {
		i__2 = k * cqm_dim1 + 1;
		cqm[i__2].r = cqf0.r, cqm[i__2].i = cqf0.i;
	    }
	    cqf2.r = cqf1.r, cqf2.i = cqf1.i;
/* L35: */
	    cqf1.r = cqf0.r, cqf1.i = cqf0.i;
	}
	z_div(&z__1, &c_b1387, &zq);
	cq10.r = z__1.r, cq10.i = z__1.i;
	i__2 = *n;
	for (k = 0; k <= i__2; ++k) {
/* L40: */
	    i__3 = k * cqm_dim1 + 1;
	    i__1 = k * cqm_dim1 + 1;
	    z__2.r = cq10.r * cqm[i__1].r - cq10.i * cqm[i__1].i, z__2.i = 
		    cq10.r * cqm[i__1].i + cq10.i * cqm[i__1].r;
	    z_div(&z__1, &z__2, &cqf0);
	    cqm[i__3].r = z__1.r, cqm[i__3].i = z__1.i;
	}
	i__3 = *n;
	for (j = 0; j <= i__3; ++j) {
	    i__1 = j * cqm_dim1;
	    cq0.r = cqm[i__1].r, cq0.i = cqm[i__1].i;
	    i__1 = j * cqm_dim1 + 1;
	    cq1.r = cqm[i__1].r, cq1.i = cqm[i__1].i;
	    i__1 = *m - 2;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		d__1 = (i__ + 1) * -2.;
		z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
		z_div(&z__3, &z__4, &zq);
		z__2.r = z__3.r * cq1.r - z__3.i * cq1.i, z__2.i = z__3.r * 
			cq1.i + z__3.i * cq1.r;
		d__2 = (j - i__) * (j + i__ + 1.);
		z__5.r = d__2 * cq0.r, z__5.i = d__2 * cq0.i;
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
		cqf.r = z__1.r, cqf.i = z__1.i;
		i__2 = i__ + 2 + j * cqm_dim1;
		cqm[i__2].r = cqf.r, cqm[i__2].i = cqf.i;
		cq0.r = cq1.r, cq0.i = cq1.i;
		cq1.r = cqf.r, cq1.i = cqf.i;
/* L45: */
	    }
	}
    }
    z__2.r = (doublereal) ls, z__2.i = 0.;
    z_div(&z__1, &z__2, &zs);
    cqd[0].r = z__1.r, cqd[0].i = z__1.i;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L50: */
	i__3 = j * cqd_dim1;
	i__2 = ls * j;
	i__4 = (j - 1) * cqm_dim1;
	i__5 = j * cqm_dim1;
	z__4.r = z__.r * cqm[i__5].r - z__.i * cqm[i__5].i, z__4.i = z__.r * 
		cqm[i__5].i + z__.i * cqm[i__5].r;
	z__3.r = cqm[i__4].r - z__4.r, z__3.i = cqm[i__4].i - z__4.i;
	d__1 = (doublereal) i__2;
	z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
	z_div(&z__1, &z__2, &zs);
	cqd[i__3].r = z__1.r, cqd[i__3].i = z__1.i;
    }
    i__3 = *n;
    for (j = 0; j <= i__3; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__4 = i__ + j * cqd_dim1;
	    i__5 = ls * i__;
	    d__1 = (doublereal) i__5;
	    z__4.r = d__1 * z__.r, z__4.i = d__1 * z__.i;
	    z_div(&z__3, &z__4, &zs);
	    i__1 = i__ + j * cqm_dim1;
	    z__2.r = z__3.r * cqm[i__1].r - z__3.i * cqm[i__1].i, z__2.i = 
		    z__3.r * cqm[i__1].i + z__3.i * cqm[i__1].r;
	    d__2 = (i__ + j) * (j - i__ + 1.);
	    z__7.r = d__2, z__7.i = 0.;
	    z_div(&z__6, &z__7, &zq);
	    i__6 = i__ - 1 + j * cqm_dim1;
	    z__5.r = z__6.r * cqm[i__6].r - z__6.i * cqm[i__6].i, z__5.i = 
		    z__6.r * cqm[i__6].i + z__6.i * cqm[i__6].r;
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
	    cqd[i__4].r = z__1.r, cqd[i__4].i = z__1.i;
/* L55: */
	}
    }
    return 0;
} /* clqmn_ */

/*       ********************************** */
/* Subroutine */ int segv_(integer *m, integer *n, doublereal *c__, integer *
	kd, doublereal *cv, doublereal *eg)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a[300], b[100], d__[300], e[300], f[300], g[300], h__[
	    100];
    static integer i__, j, k, l;
    static doublereal s, t;
    static integer k1;
    static doublereal t1, x1, cs, xa;
    static integer nm;
    static doublereal xb, dk0, dk1, dk2, d2k, cv0[100];
    static integer nm1, icm;


/*       ========================================================= */
/*       Purpose: Compute the characteristic values of spheroidal */
/*                wave functions */
/*       Input :  m  --- Mode parameter */
/*                n  --- Mode parameter */
/*                c  --- Spheroidal parameter */
/*                KD --- Function code */
/*                       KD=1 for Prolate; KD=-1 for Oblate */
/*       Output:  CV --- Characteristic value for given m, n and c */
/*                EG(L) --- Characteristic value for mode m and n' */
/*                          ( L = n' - m + 1 ) */
/*       ========================================================= */

    /* Parameter adjustments */
    --eg;

    /* Function Body */
    if (*c__ < 1e-10) {
	i__1 = *n - *m + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	    eg[i__] = (i__ + *m) * (i__ + *m - 1.);
	}
	goto L70;
    }
    icm = (*n - *m + 2) / 2;
    nm = (integer) ((*n - *m) * .5f + *c__) + 10;
    cs = *c__ * *c__ * *kd;
    k = 0;
    for (l = 0; l <= 1; ++l) {
	i__1 = nm;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (l == 0) {
		k = i__ - 1 << 1;
	    }
	    if (l == 1) {
		k = (i__ << 1) - 1;
	    }
	    dk0 = (doublereal) (*m + k);
	    dk1 = (doublereal) (*m + k + 1);
	    dk2 = (doublereal) (*m + k << 1);
	    d2k = (doublereal) ((*m << 1) + k);
	    a[i__ - 1] = (d2k + 2.f) * (d2k + 1.f) / ((dk2 + 3.f) * (dk2 + 
		    5.f)) * cs;
	    d__[i__ - 1] = dk0 * dk1 + (dk0 * 2.f * dk1 - *m * 2.f * *m - 1.f)
		     / ((dk2 - 1.f) * (dk2 + 3.f)) * cs;
/* L10: */
	    g[i__ - 1] = k * (k - 1.f) / ((dk2 - 3.f) * (dk2 - 1.f)) * cs;
	}
	i__1 = nm;
	for (k = 2; k <= i__1; ++k) {
	    e[k - 1] = sqrt(a[k - 2] * g[k - 1]);
/* L15: */
	    f[k - 1] = e[k - 1] * e[k - 1];
	}
	f[0] = 0.;
	e[0] = 0.;
	xa = d__[nm - 1] + (d__1 = e[nm - 1], abs(d__1));
	xb = d__[nm - 1] - (d__1 = e[nm - 1], abs(d__1));
	nm1 = nm - 1;
	i__1 = nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = (d__1 = e[i__ - 1], abs(d__1)) + (d__2 = e[i__], abs(d__2));
	    t1 = d__[i__ - 1] + t;
	    if (xa < t1) {
		xa = t1;
	    }
	    t1 = d__[i__ - 1] - t;
	    if (t1 < xb) {
		xb = t1;
	    }
/* L20: */
	}
	i__1 = icm;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__ - 1] = xa;
/* L25: */
	    h__[i__ - 1] = xb;
	}
	i__1 = icm;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = icm;
	    for (k1 = k; k1 <= i__2; ++k1) {
		if (b[k1 - 1] < b[k - 1]) {
		    b[k - 1] = b[k1 - 1];
		    goto L35;
		}
/* L30: */
	    }
L35:
	    if (k != 1) {
		if (h__[k - 1] < h__[k - 2]) {
		    h__[k - 1] = h__[k - 2];
		}
	    }
L40:
	    x1 = (b[k - 1] + h__[k - 1]) / 2.;
	    cv0[k - 1] = x1;
	    if ((d__1 = (b[k - 1] - h__[k - 1]) / x1, abs(d__1)) < 1e-14) {
		goto L50;
	    }
	    j = 0;
	    s = 1.;
	    i__2 = nm;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (s == 0.) {
		    s += 1e-30;
		}
		t = f[i__ - 1] / s;
		s = d__[i__ - 1] - t - x1;
		if (s < 0.) {
		    ++j;
		}
/* L45: */
	    }
	    if (j < k) {
		h__[k - 1] = x1;
	    } else {
		b[k - 1] = x1;
		if (j >= icm) {
		    b[icm - 1] = x1;
		} else {
		    if (h__[j] < x1) {
			h__[j] = x1;
		    }
		    if (x1 < b[j - 1]) {
			b[j - 1] = x1;
		    }
		}
	    }
	    goto L40;
L50:
	    cv0[k - 1] = x1;
	    if (l == 0) {
		eg[(k << 1) - 1] = cv0[k - 1];
	    }
	    if (l == 1) {
		eg[k * 2] = cv0[k - 1];
	    }
/* L55: */
	}
/* L60: */
    }
L70:
    *cv = eg[*n - *m + 1];
    return 0;
} /* segv_ */

/*       ********************************** */
/* Subroutine */ int ciknb_(integer *n, doublecomplex *z__, integer *nm, 
	doublecomplex *cbi, doublecomplex *cdi, doublecomplex *cbk, 
	doublecomplex *cdk)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_exp(
	    doublecomplex *, doublecomplex *), z_log(doublecomplex *, 
	    doublecomplex *), z_sqrt(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer k, l, m;
    static doublereal a0;
    static integer k0;
    static doublecomplex z1, cf, cg, ci;
    static doublereal el;
    static doublecomplex cr;
    static doublereal pi, vt;
    static doublecomplex ca0, cf0, cf1, cg0, cg1, cs0;
    static doublereal fac;
    static doublecomplex cbs, csk0, cbkl;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ============================================================ */
/*       Purpose: Compute modified Bessel functions In(z) and Kn(z), */
/*                and their derivatives for a complex argument */
/*       Input:   z --- Complex argument */
/*                n --- Order of In(z) and Kn(z) */
/*       Output:  CBI(n) --- In(z) */
/*                CDI(n) --- In'(z) */
/*                CBK(n) --- Kn(z) */
/*                CDK(n) --- Kn'(z) */
/*                NM --- Highest order computed */
/*       Routones called: */
/*                MSTA1 and MSTA2 to compute the starting point for */
/*                backward recurrence */
/*       =========================================================== */

    pi = 3.141592653589793;
    el = .57721566490153;
    a0 = z_abs(z__);
    *nm = *n;
    if (a0 < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbi[i__2].r = 0., cbi[i__2].i = 0.;
	    i__2 = k;
	    cbk[i__2].r = 1e300, cbk[i__2].i = 0.;
	    i__2 = k;
	    cdi[i__2].r = 0., cdi[i__2].i = 0.;
/* L10: */
	    i__2 = k;
	    cdk[i__2].r = -1e300, cdk[i__2].i = -0.;
	}
	cbi[0].r = 1., cbi[0].i = 0.;
	cdi[1].r = .5, cdi[1].i = 0.;
	return 0;
    }
    z1.r = z__->r, z1.i = z__->i;
    ci.r = 0., ci.i = 1.;
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (*n == 0) {
	*nm = 1;
    }
    m = msta1_(&a0, &c__200);
    if (m < *nm) {
	*nm = m;
    } else {
	m = msta2_(&a0, nm, &c__15);
    }
    cbs.r = 0., cbs.i = 0.;
    csk0.r = 0., csk0.i = 0.;
    cf0.r = 0., cf0.i = 0.;
    cf1.r = 1e-100, cf1.i = 0.;
    for (k = m; k >= 0; --k) {
	d__1 = (k + 1.) * 2.;
	z__3.r = d__1 * cf1.r, z__3.i = d__1 * cf1.i;
	z_div(&z__2, &z__3, &z1);
	z__1.r = z__2.r + cf0.r, z__1.i = z__2.i + cf0.i;
	cf.r = z__1.r, cf.i = z__1.i;
	if (k <= *nm) {
	    i__2 = k;
	    cbi[i__2].r = cf.r, cbi[i__2].i = cf.i;
	}
	if (k != 0 && k == k / 2 << 1) {
	    z__3.r = cf.r * 4., z__3.i = cf.i * 4.;
	    d__1 = (doublereal) k;
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
	    z__1.r = csk0.r + z__2.r, z__1.i = csk0.i + z__2.i;
	    csk0.r = z__1.r, csk0.i = z__1.i;
	}
	z__2.r = cf.r * 2., z__2.i = cf.i * 2.;
	z__1.r = cbs.r + z__2.r, z__1.i = cbs.i + z__2.i;
	cbs.r = z__1.r, cbs.i = z__1.i;
	cf0.r = cf1.r, cf0.i = cf1.i;
/* L15: */
	cf1.r = cf.r, cf1.i = cf.i;
    }
    z_exp(&z__2, &z1);
    z__3.r = cbs.r - cf.r, z__3.i = cbs.i - cf.i;
    z_div(&z__1, &z__2, &z__3);
    cs0.r = z__1.r, cs0.i = z__1.i;
    i__2 = *nm;
    for (k = 0; k <= i__2; ++k) {
/* L20: */
	i__1 = k;
	i__3 = k;
	z__1.r = cs0.r * cbi[i__3].r - cs0.i * cbi[i__3].i, z__1.i = cs0.r * 
		cbi[i__3].i + cs0.i * cbi[i__3].r;
	cbi[i__1].r = z__1.r, cbi[i__1].i = z__1.i;
    }
    if (a0 <= 9.f) {
	z__6.r = z1.r * .5, z__6.i = z1.i * .5;
	z_log(&z__5, &z__6);
	z__4.r = z__5.r + el, z__4.i = z__5.i;
	z__3.r = -z__4.r, z__3.i = -z__4.i;
	z__2.r = z__3.r * cbi[0].r - z__3.i * cbi[0].i, z__2.i = z__3.r * cbi[
		0].i + z__3.i * cbi[0].r;
	z__7.r = cs0.r * csk0.r - cs0.i * csk0.i, z__7.i = cs0.r * csk0.i + 
		cs0.i * csk0.r;
	z__1.r = z__2.r + z__7.r, z__1.i = z__2.i + z__7.i;
	cbk[0].r = z__1.r, cbk[0].i = z__1.i;
	z_div(&z__3, &c_b6, &z1);
	z__4.r = cbi[1].r * cbk[0].r - cbi[1].i * cbk[0].i, z__4.i = cbi[1].r 
		* cbk[0].i + cbi[1].i * cbk[0].r;
	z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	z_div(&z__1, &z__2, cbi);
	cbk[1].r = z__1.r, cbk[1].i = z__1.i;
    } else {
	z__4.r = pi, z__4.i = 0.;
	z__5.r = z1.r * 2., z__5.i = z1.i * 2.;
	z_div(&z__3, &z__4, &z__5);
	z_sqrt(&z__2, &z__3);
	z__7.r = -z1.r, z__7.i = -z1.i;
	z_exp(&z__6, &z__7);
	z__1.r = z__2.r * z__6.r - z__2.i * z__6.i, z__1.i = z__2.r * z__6.i 
		+ z__2.i * z__6.r;
	ca0.r = z__1.r, ca0.i = z__1.i;
	k0 = 16;
	if (a0 >= 25.f) {
	    k0 = 10;
	}
	if (a0 >= 80.f) {
	    k0 = 8;
	}
	if (a0 >= 200.f) {
	    k0 = 6;
	}
	for (l = 0; l <= 1; ++l) {
	    cbkl.r = 1., cbkl.i = 0.;
	    vt = l * 4.;
	    cr.r = 1., cr.i = 0.;
	    i__1 = k0;
	    for (k = 1; k <= i__1; ++k) {
		z__3.r = cr.r * .125, z__3.i = cr.i * .125;
/* Computing 2nd power */
		r__1 = k * 2.f - 1.f;
		d__1 = vt - r__1 * r__1;
		z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
		d__2 = (doublereal) k;
		z__4.r = d__2 * z1.r, z__4.i = d__2 * z1.i;
		z_div(&z__1, &z__2, &z__4);
		cr.r = z__1.r, cr.i = z__1.i;
/* L25: */
		z__1.r = cbkl.r + cr.r, z__1.i = cbkl.i + cr.i;
		cbkl.r = z__1.r, cbkl.i = z__1.i;
	    }
	    i__1 = l;
	    z__1.r = ca0.r * cbkl.r - ca0.i * cbkl.i, z__1.i = ca0.r * cbkl.i 
		    + ca0.i * cbkl.r;
	    cbk[i__1].r = z__1.r, cbk[i__1].i = z__1.i;
/* L30: */
	}
    }
    cg0.r = cbk[0].r, cg0.i = cbk[0].i;
    cg1.r = cbk[1].r, cg1.i = cbk[1].i;
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
	d__1 = (k - 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, &z1);
	z__2.r = z__3.r * cg1.r - z__3.i * cg1.i, z__2.i = z__3.r * cg1.i + 
		z__3.i * cg1.r;
	z__1.r = z__2.r + cg0.r, z__1.i = z__2.i + cg0.i;
	cg.r = z__1.r, cg.i = z__1.i;
	i__3 = k;
	cbk[i__3].r = cg.r, cbk[i__3].i = cg.i;
	cg0.r = cg1.r, cg0.i = cg1.i;
/* L35: */
	cg1.r = cg.r, cg1.i = cg.i;
    }
    if (z__->r < 0.f) {
	fac = 1.;
	i__1 = *nm;
	for (k = 0; k <= i__1; ++k) {
	    if (d_imag(z__) < 0.f) {
		i__3 = k;
		i__2 = k;
		z__2.r = fac * cbk[i__2].r, z__2.i = fac * cbk[i__2].i;
		z__4.r = pi * ci.r, z__4.i = pi * ci.i;
		i__4 = k;
		z__3.r = z__4.r * cbi[i__4].r - z__4.i * cbi[i__4].i, z__3.i =
			 z__4.r * cbi[i__4].i + z__4.i * cbi[i__4].r;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		cbk[i__3].r = z__1.r, cbk[i__3].i = z__1.i;
	    } else {
		i__3 = k;
		i__2 = k;
		z__2.r = fac * cbk[i__2].r, z__2.i = fac * cbk[i__2].i;
		z__4.r = pi * ci.r, z__4.i = pi * ci.i;
		i__4 = k;
		z__3.r = z__4.r * cbi[i__4].r - z__4.i * cbi[i__4].i, z__3.i =
			 z__4.r * cbi[i__4].i + z__4.i * cbi[i__4].r;
		z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
		cbk[i__3].r = z__1.r, cbk[i__3].i = z__1.i;
	    }
	    i__3 = k;
	    i__2 = k;
	    z__1.r = fac * cbi[i__2].r, z__1.i = fac * cbi[i__2].i;
	    cbi[i__3].r = z__1.r, cbi[i__3].i = z__1.i;
	    fac = -fac;
/* L45: */
	}
    }
    cdi[0].r = cbi[1].r, cdi[0].i = cbi[1].i;
    z__1.r = -cbk[1].r, z__1.i = -cbk[1].i;
    cdk[0].r = z__1.r, cdk[0].i = z__1.i;
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
	i__3 = k;
	i__2 = k - 1;
	z__4.r = (doublereal) k, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__4 = k;
	z__2.r = z__3.r * cbi[i__4].r - z__3.i * cbi[i__4].i, z__2.i = z__3.r 
		* cbi[i__4].i + z__3.i * cbi[i__4].r;
	z__1.r = cbi[i__2].r - z__2.r, z__1.i = cbi[i__2].i - z__2.i;
	cdi[i__3].r = z__1.r, cdi[i__3].i = z__1.i;
/* L50: */
	i__3 = k;
	i__2 = k - 1;
	z__2.r = -cbk[i__2].r, z__2.i = -cbk[i__2].i;
	z__5.r = (doublereal) k, z__5.i = 0.;
	z_div(&z__4, &z__5, z__);
	i__4 = k;
	z__3.r = z__4.r * cbk[i__4].r - z__4.i * cbk[i__4].i, z__3.i = z__4.r 
		* cbk[i__4].i + z__4.i * cbk[i__4].r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	cdk[i__3].r = z__1.r, cdk[i__3].i = z__1.i;
    }
    return 0;
} /* ciknb_ */

/*       ********************************** */
/* Subroutine */ int cikna_(integer *n, doublecomplex *z__, integer *nm, 
	doublecomplex *cbi, doublecomplex *cdi, doublecomplex *cbk, 
	doublecomplex *cdk)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, m;
    static doublereal a0;
    static doublecomplex cf, cs, cf1, cf2, ckk, cbi0, cbi1, cdi0, cdi1, cbk0, 
	    cdk0, cbk1, cdk1;
    extern /* Subroutine */ int cik01_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , doublecomplex *, doublecomplex *, doublecomplex *);
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);


/*       ======================================================== */
/*       Purpose: Compute modified Bessel functions In(z), Kn(x) */
/*                and their derivatives for a complex argument */
/*       Input :  z --- Complex argument of In(z) and Kn(z) */
/*                n --- Order of In(z) and Kn(z) */
/*       Output:  CBI(n) --- In(z) */
/*                CDI(n) --- In'(z) */
/*                CBK(n) --- Kn(z) */
/*                CDK(n) --- Kn'(z) */
/*                NM --- Highest order computed */
/*       Routines called: */
/*             (1) CIK01 to compute I0(z), I1(z) K0(z) & K1(z) */
/*             (2) MSTA1 and MSTA2 to compute the starting */
/*                 point for backward recurrence */
/*       ======================================================== */

    a0 = z_abs(z__);
    *nm = *n;
    if (a0 < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    i__2 = k;
	    cbi[i__2].r = 0., cbi[i__2].i = 0.;
	    i__2 = k;
	    cdi[i__2].r = 0., cdi[i__2].i = 0.;
	    i__2 = k;
	    cbk[i__2].r = -1e300, cbk[i__2].i = -0.;
/* L10: */
	    i__2 = k;
	    cdk[i__2].r = 1e300, cdk[i__2].i = 0.;
	}
	cbi[0].r = 1., cbi[0].i = 0.;
	cdi[1].r = .5, cdi[1].i = 0.;
	return 0;
    }
    cik01_(z__, &cbi0, &cdi0, &cbi1, &cdi1, &cbk0, &cdk0, &cbk1, &cdk1);
    cbi[0].r = cbi0.r, cbi[0].i = cbi0.i;
    cbi[1].r = cbi1.r, cbi[1].i = cbi1.i;
    cbk[0].r = cbk0.r, cbk[0].i = cbk0.i;
    cbk[1].r = cbk1.r, cbk[1].i = cbk1.i;
    cdi[0].r = cdi0.r, cdi[0].i = cdi0.i;
    cdi[1].r = cdi1.r, cdi[1].i = cdi1.i;
    cdk[0].r = cdk0.r, cdk[0].i = cdk0.i;
    cdk[1].r = cdk1.r, cdk[1].i = cdk1.i;
    if (*n <= 1) {
	return 0;
    }
    m = msta1_(&a0, &c__200);
    if (m < *n) {
	*nm = m;
    } else {
	m = msta2_(&a0, n, &c__15);
    }
    cf2.r = 0., cf2.i = 0.;
    cf1.r = 1e-100, cf1.i = 0.;
    for (k = m; k >= 0; --k) {
	d__1 = (k + 1.) * 2.;
	z__4.r = d__1, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	z__2.r = z__3.r * cf1.r - z__3.i * cf1.i, z__2.i = z__3.r * cf1.i + 
		z__3.i * cf1.r;
	z__1.r = z__2.r + cf2.r, z__1.i = z__2.i + cf2.i;
	cf.r = z__1.r, cf.i = z__1.i;
	if (k <= *nm) {
	    i__2 = k;
	    cbi[i__2].r = cf.r, cbi[i__2].i = cf.i;
	}
	cf2.r = cf1.r, cf2.i = cf1.i;
/* L45: */
	cf1.r = cf.r, cf1.i = cf.i;
    }
    z_div(&z__1, &cbi0, &cf);
    cs.r = z__1.r, cs.i = z__1.i;
    i__2 = *nm;
    for (k = 0; k <= i__2; ++k) {
/* L50: */
	i__1 = k;
	i__3 = k;
	z__1.r = cs.r * cbi[i__3].r - cs.i * cbi[i__3].i, z__1.i = cs.r * cbi[
		i__3].i + cs.i * cbi[i__3].r;
	cbi[i__1].r = z__1.r, cbi[i__1].i = z__1.i;
    }
    i__1 = *nm;
    for (k = 2; k <= i__1; ++k) {
	if (z_abs(&cbi[k - 1]) > z_abs(&cbi[k - 2])) {
	    z_div(&z__3, &c_b6, z__);
	    i__3 = k;
	    i__2 = k - 1;
	    z__4.r = cbi[i__3].r * cbk[i__2].r - cbi[i__3].i * cbk[i__2].i, 
		    z__4.i = cbi[i__3].r * cbk[i__2].i + cbi[i__3].i * cbk[
		    i__2].r;
	    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
	    z_div(&z__1, &z__2, &cbi[k - 1]);
	    ckk.r = z__1.r, ckk.i = z__1.i;
	} else {
	    i__3 = k;
	    i__2 = k - 2;
	    z__3.r = cbi[i__3].r * cbk[i__2].r - cbi[i__3].i * cbk[i__2].i, 
		    z__3.i = cbi[i__3].r * cbk[i__2].i + cbi[i__3].i * cbk[
		    i__2].r;
	    d__1 = (k - 1.) * 2.;
	    z__5.r = d__1, z__5.i = 0.;
	    z__6.r = z__->r * z__->r - z__->i * z__->i, z__6.i = z__->r * 
		    z__->i + z__->i * z__->r;
	    z_div(&z__4, &z__5, &z__6);
	    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
	    z_div(&z__1, &z__2, &cbi[k - 2]);
	    ckk.r = z__1.r, ckk.i = z__1.i;
	}
/* L60: */
	i__3 = k;
	cbk[i__3].r = ckk.r, cbk[i__3].i = ckk.i;
    }
    i__3 = *nm;
    for (k = 2; k <= i__3; ++k) {
	i__1 = k;
	i__2 = k - 1;
	z__4.r = (doublereal) k, z__4.i = 0.;
	z_div(&z__3, &z__4, z__);
	i__4 = k;
	z__2.r = z__3.r * cbi[i__4].r - z__3.i * cbi[i__4].i, z__2.i = z__3.r 
		* cbi[i__4].i + z__3.i * cbi[i__4].r;
	z__1.r = cbi[i__2].r - z__2.r, z__1.i = cbi[i__2].i - z__2.i;
	cdi[i__1].r = z__1.r, cdi[i__1].i = z__1.i;
/* L70: */
	i__1 = k;
	i__2 = k - 1;
	z__2.r = -cbk[i__2].r, z__2.i = -cbk[i__2].i;
	z__5.r = (doublereal) k, z__5.i = 0.;
	z_div(&z__4, &z__5, z__);
	i__4 = k;
	z__3.r = z__4.r * cbk[i__4].r - z__4.i * cbk[i__4].i, z__3.i = z__4.r 
		* cbk[i__4].i + z__4.i * cbk[i__4].r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	cdk[i__1].r = z__1.r, cdk[i__1].i = z__1.i;
    }
    return 0;
} /* cikna_ */

/*       ********************************** */
/* Subroutine */ int mtu12_(integer *kf, integer *kc, integer *m, doublereal *
	q, doublereal *x, doublereal *f1r, doublereal *d1r, doublereal *f2r, 
	doublereal *d2r)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal a;
    static integer k;
    static doublereal c1, c2, u1, u2, w1, w2;
    static integer ic;
    static doublereal fg[251];
    static integer kd, km, nm;
    static doublereal qm, bj1[252], bj2[252], dj1[252], dj2[252], by1[252], 
	    by2[252], dy1[252], dy2[252], eps;
    extern /* Subroutine */ int cva2_(integer *, integer *, doublereal *, 
	    doublereal *);
    extern doublereal dnan_(void);
    extern /* Subroutine */ int jynb_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), fcoef_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *);


/*       ============================================================== */
/*       Purpose: Compute modified Mathieu functions of the first and */
/*                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q), */
/*                and their derivatives */
/*       Input:   KF --- Function code */
/*                       KF=1 for computing Mcm(x,q) */
/*                       KF=2 for computing Msm(x,q) */
/*                KC --- Function Code */
/*                       KC=1 for computing the first kind */
/*                       KC=2 for computing the second kind */
/*                            or Msm(2)(x,q) and Msm(2)'(x,q) */
/*                       KC=3 for computing both the first */
/*                            and second kinds */
/*                m  --- Order of Mathieu functions */
/*                q  --- Parameter of Mathieu functions ( q ≥ 0 ) */
/*                x  --- Argument of Mathieu functions */
/*       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q) */
/*                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q) */
/*                F2R --- Mcm(2)(x,q) or Msm(2)(x,q) */
/*                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q) */
/*       Routines called: */
/*            (1) CVA2 for computing the characteristic values */
/*            (2) FCOEF for computing expansion coefficients */
/*            (3) JYNB for computing Jn(x), Yn(x) and their */
/*                derivatives */
/*       ============================================================== */

    eps = 1e-14;
    if (*kf == 1 && *m == *m / 2 << 1) {
	kd = 1;
    }
    if (*kf == 1 && *m != *m / 2 << 1) {
	kd = 2;
    }
    if (*kf == 2 && *m != *m / 2 << 1) {
	kd = 3;
    }
    if (*kf == 2 && *m == *m / 2 << 1) {
	kd = 4;
    }
    cva2_(&kd, m, q, &a);
    if (*q <= 1.) {
	qm = sqrt(*q) * 56.1f + 7.5f - *q * 134.7f + sqrt(*q) * 90.7f * *q;
    } else {
	qm = sqrt(*q) * 3.1f + 17.f - *q * .126f + sqrt(*q) * .0037f * *q;
    }
    km = (integer) (qm + *m * .5f);
    if (km >= 251) {
	*f1r = dnan_();
	*d1r = dnan_();
	*f2r = dnan_();
	*d2r = dnan_();
	return 0;
    }
    fcoef_(&kd, m, q, &a, fg);
    ic = *m / 2 + 1;
    if (kd == 4) {
	ic = *m / 2;
    }
    c1 = exp(-(*x));
    c2 = exp(*x);
    u1 = sqrt(*q) * c1;
    u2 = sqrt(*q) * c2;
    i__1 = km + 1;
    jynb_(&i__1, &u1, &nm, bj1, dj1, by1, dy1);
    i__1 = km + 1;
    jynb_(&i__1, &u2, &nm, bj2, dj2, by2, dy2);
    w1 = 0.;
    w2 = 0.;
    if (*kc == 2) {
	goto L50;
    }
    *f1r = 0.;
    i__1 = km;
    for (k = 1; k <= i__1; ++k) {
	if (kd == 1) {
	    i__2 = ic + k;
	    *f1r += pow_ii(&c_n1, &i__2) * fg[k - 1] * bj1[k - 1] * bj2[k - 1]
		    ;
	} else if (kd == 2 || kd == 3) {
	    i__2 = ic + k;
	    *f1r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (bj1[k - 1] * bj2[k] + 
		    pow_ii(&c_n1, &kd) * bj1[k] * bj2[k - 1]);
	} else {
	    i__2 = ic + k;
	    *f1r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (bj1[k - 1] * bj2[k + 
		    1] - bj1[k + 1] * bj2[k - 1]);
	}
	if (k >= 5 && (d__1 = *f1r - w1, abs(d__1)) < abs(*f1r) * eps) {
	    goto L35;
	}
/* L30: */
	w1 = *f1r;
    }
L35:
    *f1r /= fg[0];
    *d1r = 0.;
    i__1 = km;
    for (k = 1; k <= i__1; ++k) {
	if (kd == 1) {
	    i__2 = ic + k;
	    *d1r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (c2 * bj1[k - 1] * dj2[
		    k - 1] - c1 * dj1[k - 1] * bj2[k - 1]);
	} else if (kd == 2 || kd == 3) {
	    i__2 = ic + k;
	    *d1r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (c2 * (bj1[k - 1] * 
		    dj2[k] + pow_ii(&c_n1, &kd) * bj1[k] * dj2[k - 1]) - c1 * 
		    (dj1[k - 1] * bj2[k] + pow_ii(&c_n1, &kd) * dj1[k] * bj2[
		    k - 1]));
	} else {
	    i__2 = ic + k;
	    *d1r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (c2 * (bj1[k - 1] * 
		    dj2[k + 1] - bj1[k + 1] * dj2[k - 1]) - c1 * (dj1[k - 1] *
		     bj2[k + 1] - dj1[k + 1] * bj2[k - 1]));
	}
	if (k >= 5 && (d__1 = *d1r - w2, abs(d__1)) < abs(*d1r) * eps) {
	    goto L45;
	}
/* L40: */
	w2 = *d1r;
    }
L45:
    *d1r = *d1r * sqrt(*q) / fg[0];
    if (*kc == 1) {
	return 0;
    }
L50:
    *f2r = 0.;
    i__1 = km;
    for (k = 1; k <= i__1; ++k) {
	if (kd == 1) {
	    i__2 = ic + k;
	    *f2r += pow_ii(&c_n1, &i__2) * fg[k - 1] * bj1[k - 1] * by2[k - 1]
		    ;
	} else if (kd == 2 || kd == 3) {
	    i__2 = ic + k;
	    *f2r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (bj1[k - 1] * by2[k] + 
		    pow_ii(&c_n1, &kd) * bj1[k] * by2[k - 1]);
	} else {
	    i__2 = ic + k;
	    *f2r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (bj1[k - 1] * by2[k + 
		    1] - bj1[k + 1] * by2[k - 1]);
	}
	if (k >= 5 && (d__1 = *f2r - w1, abs(d__1)) < abs(*f2r) * eps) {
	    goto L60;
	}
/* L55: */
	w1 = *f2r;
    }
L60:
    *f2r /= fg[0];
    *d2r = 0.;
    i__1 = km;
    for (k = 1; k <= i__1; ++k) {
	if (kd == 1) {
	    i__2 = ic + k;
	    *d2r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (c2 * bj1[k - 1] * dy2[
		    k - 1] - c1 * dj1[k - 1] * by2[k - 1]);
	} else if (kd == 2 || kd == 3) {
	    i__2 = ic + k;
	    *d2r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (c2 * (bj1[k - 1] * 
		    dy2[k] + pow_ii(&c_n1, &kd) * bj1[k] * dy2[k - 1]) - c1 * 
		    (dj1[k - 1] * by2[k] + pow_ii(&c_n1, &kd) * dj1[k] * by2[
		    k - 1]));
	} else {
	    i__2 = ic + k;
	    *d2r += pow_ii(&c_n1, &i__2) * fg[k - 1] * (c2 * (bj1[k - 1] * 
		    dy2[k + 1] - bj1[k + 1] * dy2[k - 1]) - c1 * (dj1[k - 1] *
		     by2[k + 1] - dj1[k + 1] * by2[k - 1]));
	}
	if (k >= 5 && (d__1 = *d2r - w2, abs(d__1)) < abs(*d2r) * eps) {
	    goto L70;
	}
/* L65: */
	w2 = *d2r;
    }
L70:
    *d2r = *d2r * sqrt(*q) / fg[0];
    return 0;
} /* mtu12_ */

/*       ********************************** */
/* Subroutine */ int cik01_(doublecomplex *z__, doublecomplex *cbi0, 
	doublecomplex *cdi0, doublecomplex *cbi1, doublecomplex *cdi1, 
	doublecomplex *cbk0, doublecomplex *cdk0, doublecomplex *cbk1, 
	doublecomplex *cdk1)
{
    /* Initialized data */

    static doublereal a[12] = { .125,.0703125,.0732421875,.11215209960938,
	    .22710800170898,.57250142097473,1.7277275025845,6.0740420012735,
	    24.380529699556,110.01714026925,551.33589612202,3038.0905109224 };
    static doublereal b[12] = { -.375,-.1171875,-.1025390625,-.14419555664063,
	    -.2775764465332,-.67659258842468,-1.9935317337513,
	    -6.8839142681099,-27.248827311269,-121.59789187654,
	    -603.84407670507,-3302.2722944809 };
    static doublereal a1[10] = { .125,.2109375,1.0986328125,11.775970458984,
	    214.61706161499,5951.1522710323,233476.45606175,12312234.987631,
	    840139034.6421,72031420482.627 };

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_exp(
	    doublecomplex *, doublecomplex *), z_sqrt(doublecomplex *, 
	    doublecomplex *), pow_zi(doublecomplex *, doublecomplex *, 
	    integer *), z_log(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer k;
    static doublereal a0;
    static integer k0;
    static doublereal w0;
    static doublecomplex z1, z2, ca, cb, ci, cr, cs, ct;
    static doublereal pi;
    static doublecomplex cw, zr, zr2;


/*       ========================================================== */
/*       Purpose: Compute modified Bessel functions I0(z), I1(z), */
/*                K0(z), K1(z), and their derivatives for a */
/*                complex argument */
/*       Input :  z --- Complex argument */
/*       Output:  CBI0 --- I0(z) */
/*                CDI0 --- I0'(z) */
/*                CBI1 --- I1(z) */
/*                CDI1 --- I1'(z) */
/*                CBK0 --- K0(z) */
/*                CDK0 --- K0'(z) */
/*                CBK1 --- K1(z) */
/*                CDK1 --- K1'(z) */
/*       ========================================================== */

    pi = 3.141592653589793;
    ci.r = 0., ci.i = 1.;
    a0 = z_abs(z__);
    z__1.r = z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
	    z__->i * z__->r;
    z2.r = z__1.r, z2.i = z__1.i;
    z1.r = z__->r, z1.i = z__->i;
    if (a0 == 0.) {
	cbi0->r = 1., cbi0->i = 0.;
	cbi1->r = 0., cbi1->i = 0.;
	cdi0->r = 0., cdi0->i = 0.;
	cdi1->r = .5, cdi1->i = 0.;
	cbk0->r = 1e300, cbk0->i = 0.;
	cbk1->r = 1e300, cbk1->i = 0.;
	cdk0->r = -1e300, cdk0->i = -0.;
	cdk1->r = -1e300, cdk1->i = -0.;
	return 0;
    }
    if (z__->r < 0.f) {
	z__1.r = -z__->r, z__1.i = -z__->i;
	z1.r = z__1.r, z1.i = z__1.i;
    }
    if (a0 <= 18.f) {
	cbi0->r = 1., cbi0->i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 50; ++k) {
	    z__3.r = cr.r * .25, z__3.i = cr.i * .25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    i__1 = k * k;
	    d__1 = (doublereal) i__1;
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cbi0->r + cr.r, z__1.i = cbi0->i + cr.i;
	    cbi0->r = z__1.r, cbi0->i = z__1.i;
	    z_div(&z__1, &cr, cbi0);
	    if (z_abs(&z__1) < 1e-15) {
		goto L15;
	    }
/* L10: */
	}
L15:
	cbi1->r = 1., cbi1->i = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 50; ++k) {
	    z__3.r = cr.r * .25, z__3.i = cr.i * .25;
	    z__2.r = z__3.r * z2.r - z__3.i * z2.i, z__2.i = z__3.r * z2.i + 
		    z__3.i * z2.r;
	    i__1 = k * (k + 1);
	    d__1 = (doublereal) i__1;
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__1.r = cbi1->r + cr.r, z__1.i = cbi1->i + cr.i;
	    cbi1->r = z__1.r, cbi1->i = z__1.i;
	    z_div(&z__1, &cr, cbi1);
	    if (z_abs(&z__1) < 1e-15) {
		goto L25;
	    }
/* L20: */
	}
L25:
	z__2.r = z1.r * .5, z__2.i = z1.i * .5;
	z__1.r = z__2.r * cbi1->r - z__2.i * cbi1->i, z__1.i = z__2.r * 
		cbi1->i + z__2.i * cbi1->r;
	cbi1->r = z__1.r, cbi1->i = z__1.i;
    } else {
	k0 = 12;
	if (a0 >= 35.f) {
	    k0 = 9;
	}
	if (a0 >= 50.f) {
	    k0 = 7;
	}
	z_exp(&z__2, &z1);
	d__1 = pi * 2.;
	z__4.r = d__1 * z1.r, z__4.i = d__1 * z1.i;
	z_sqrt(&z__3, &z__4);
	z_div(&z__1, &z__2, &z__3);
	ca.r = z__1.r, ca.i = z__1.i;
	cbi0->r = 1., cbi0->i = 0.;
	z_div(&z__1, &c_b6, &z1);
	zr.r = z__1.r, zr.i = z__1.i;
	i__1 = k0;
	for (k = 1; k <= i__1; ++k) {
/* L30: */
	    i__2 = k - 1;
	    pow_zi(&z__3, &zr, &k);
	    z__2.r = a[i__2] * z__3.r, z__2.i = a[i__2] * z__3.i;
	    z__1.r = cbi0->r + z__2.r, z__1.i = cbi0->i + z__2.i;
	    cbi0->r = z__1.r, cbi0->i = z__1.i;
	}
	z__1.r = ca.r * cbi0->r - ca.i * cbi0->i, z__1.i = ca.r * cbi0->i + 
		ca.i * cbi0->r;
	cbi0->r = z__1.r, cbi0->i = z__1.i;
	cbi1->r = 1., cbi1->i = 0.;
	i__2 = k0;
	for (k = 1; k <= i__2; ++k) {
/* L35: */
	    i__1 = k - 1;
	    pow_zi(&z__3, &zr, &k);
	    z__2.r = b[i__1] * z__3.r, z__2.i = b[i__1] * z__3.i;
	    z__1.r = cbi1->r + z__2.r, z__1.i = cbi1->i + z__2.i;
	    cbi1->r = z__1.r, cbi1->i = z__1.i;
	}
	z__1.r = ca.r * cbi1->r - ca.i * cbi1->i, z__1.i = ca.r * cbi1->i + 
		ca.i * cbi1->r;
	cbi1->r = z__1.r, cbi1->i = z__1.i;
    }
    if (a0 <= 9.f) {
	cs.r = 0., cs.i = 0.;
	z__4.r = z1.r * .5, z__4.i = z1.i * .5;
	z_log(&z__3, &z__4);
	z__2.r = -z__3.r, z__2.i = -z__3.i;
	z__1.r = z__2.r - .5772156649015329, z__1.i = z__2.i;
	ct.r = z__1.r, ct.i = z__1.i;
	w0 = 0.;
	cr.r = 1., cr.i = 0.;
	for (k = 1; k <= 50; ++k) {
	    w0 += 1. / k;
	    z__3.r = cr.r * .25, z__3.i = cr.i * .25;
	    i__1 = k * k;
	    d__1 = (doublereal) i__1;
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
	    z__1.r = z__2.r * z2.r - z__2.i * z2.i, z__1.i = z__2.r * z2.i + 
		    z__2.i * z2.r;
	    cr.r = z__1.r, cr.i = z__1.i;
	    z__3.r = w0 + ct.r, z__3.i = ct.i;
	    z__2.r = cr.r * z__3.r - cr.i * z__3.i, z__2.i = cr.r * z__3.i + 
		    cr.i * z__3.r;
	    z__1.r = cs.r + z__2.r, z__1.i = cs.i + z__2.i;
	    cs.r = z__1.r, cs.i = z__1.i;
	    z__2.r = cs.r - cw.r, z__2.i = cs.i - cw.i;
	    z_div(&z__1, &z__2, &cs);
	    if (z_abs(&z__1) < 1e-15) {
		goto L45;
	    }
/* L40: */
	    cw.r = cs.r, cw.i = cs.i;
	}
L45:
	z__1.r = ct.r + cs.r, z__1.i = ct.i + cs.i;
	cbk0->r = z__1.r, cbk0->i = z__1.i;
    } else {
	z_div(&z__1, &c_b340, &z1);
	cb.r = z__1.r, cb.i = z__1.i;
	z_div(&z__1, &c_b6, &z2);
	zr2.r = z__1.r, zr2.i = z__1.i;
	cbk0->r = 1., cbk0->i = 0.;
	for (k = 1; k <= 10; ++k) {
/* L50: */
	    i__1 = k - 1;
	    pow_zi(&z__3, &zr2, &k);
	    z__2.r = a1[i__1] * z__3.r, z__2.i = a1[i__1] * z__3.i;
	    z__1.r = cbk0->r + z__2.r, z__1.i = cbk0->i + z__2.i;
	    cbk0->r = z__1.r, cbk0->i = z__1.i;
	}
	z__2.r = cb.r * cbk0->r - cb.i * cbk0->i, z__2.i = cb.r * cbk0->i + 
		cb.i * cbk0->r;
	z_div(&z__1, &z__2, cbi0);
	cbk0->r = z__1.r, cbk0->i = z__1.i;
    }
    z_div(&z__3, &c_b6, &z1);
    z__4.r = cbi1->r * cbk0->r - cbi1->i * cbk0->i, z__4.i = cbi1->r * 
	    cbk0->i + cbi1->i * cbk0->r;
    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
    z_div(&z__1, &z__2, cbi0);
    cbk1->r = z__1.r, cbk1->i = z__1.i;
    if (z__->r < 0.f) {
	if (d_imag(z__) < 0.f) {
	    z__3.r = pi * ci.r, z__3.i = pi * ci.i;
	    z__2.r = z__3.r * cbi0->r - z__3.i * cbi0->i, z__2.i = z__3.r * 
		    cbi0->i + z__3.i * cbi0->r;
	    z__1.r = cbk0->r + z__2.r, z__1.i = cbk0->i + z__2.i;
	    cbk0->r = z__1.r, cbk0->i = z__1.i;
	}
	if (d_imag(z__) > 0.f) {
	    z__3.r = pi * ci.r, z__3.i = pi * ci.i;
	    z__2.r = z__3.r * cbi0->r - z__3.i * cbi0->i, z__2.i = z__3.r * 
		    cbi0->i + z__3.i * cbi0->r;
	    z__1.r = cbk0->r - z__2.r, z__1.i = cbk0->i - z__2.i;
	    cbk0->r = z__1.r, cbk0->i = z__1.i;
	}
	if (d_imag(z__) < 0.f) {
	    z__2.r = -cbk1->r, z__2.i = -cbk1->i;
	    z__4.r = pi * ci.r, z__4.i = pi * ci.i;
	    z__3.r = z__4.r * cbi1->r - z__4.i * cbi1->i, z__3.i = z__4.r * 
		    cbi1->i + z__4.i * cbi1->r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    cbk1->r = z__1.r, cbk1->i = z__1.i;
	}
	if (d_imag(z__) > 0.f) {
	    z__2.r = -cbk1->r, z__2.i = -cbk1->i;
	    z__4.r = pi * ci.r, z__4.i = pi * ci.i;
	    z__3.r = z__4.r * cbi1->r - z__4.i * cbi1->i, z__3.i = z__4.r * 
		    cbi1->i + z__4.i * cbi1->r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cbk1->r = z__1.r, cbk1->i = z__1.i;
	}
	z__1.r = -cbi1->r, z__1.i = -cbi1->i;
	cbi1->r = z__1.r, cbi1->i = z__1.i;
    }
    cdi0->r = cbi1->r, cdi0->i = cbi1->i;
    z_div(&z__3, &c_b6, z__);
    z__2.r = z__3.r * cbi1->r - z__3.i * cbi1->i, z__2.i = z__3.r * cbi1->i + 
	    z__3.i * cbi1->r;
    z__1.r = cbi0->r - z__2.r, z__1.i = cbi0->i - z__2.i;
    cdi1->r = z__1.r, cdi1->i = z__1.i;
    z__1.r = -cbk1->r, z__1.i = -cbk1->i;
    cdk0->r = z__1.r, cdk0->i = z__1.i;
    z__2.r = -cbk0->r, z__2.i = -cbk0->i;
    z_div(&z__4, &c_b6, z__);
    z__3.r = z__4.r * cbk1->r - z__4.i * cbk1->i, z__3.i = z__4.r * cbk1->i + 
	    z__4.i * cbk1->r;
    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
    cdk1->r = z__1.r, cdk1->i = z__1.i;
    return 0;
} /* cik01_ */

/*       ********************************** */
/* Subroutine */ int cpsi_(doublereal *x, doublereal *y, doublereal *psr, 
	doublereal *psi)
{
    /* Initialized data */

    static doublereal a[8] = { -.08333333333333,.0083333333333333333,
	    -.0039682539682539683,.0041666666666666667,-.0075757575757575758,
	    .021092796092796093,-.083333333333333333,.4432598039215686 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), log(doublereal), cos(
	    doublereal), pow_di(doublereal *, integer *), sin(doublereal), 
	    pow_dd(doublereal *, doublereal *), tan(doublereal), tanh(
	    doublereal);

    /* Local variables */
    static integer k, n;
    static doublereal x0, x1, y1, z0, z2, pi, ri, th, tm, tn, rr, ct2;


/*       ============================================= */
/*       Purpose: Compute the psi function for a */
/*                complex argument */
/*       Input :  x   --- Real part of z */
/*                y   --- Imaginary part of z */
/*       Output:  PSR --- Real part of psi(z) */
/*                PSI --- Imaginary part of psi(z) */
/*       ============================================= */

    pi = 3.141592653589793;
    if (*y == 0. && *x == (doublereal) ((integer) (*x)) && *x <= 0.) {
	*psr = 1e300;
	*psi = 0.;
    } else {
	x1 = *x;
	y1 = *y;
	if (*x < 0.) {
	    *x = -(*x);
	    *y = -(*y);
	}
	x0 = *x;
	n = 0;
	if (*x < 8.) {
	    n = 8 - (integer) (*x);
	    x0 = *x + n;
	}
	th = 0.;
	if (x0 == 0. && *y != 0.) {
	    th = pi * .5;
	}
	if (x0 != 0.) {
	    th = atan(*y / x0);
	}
	z2 = x0 * x0 + *y * *y;
	z0 = sqrt(z2);
	*psr = log(z0) - x0 * .5 / z2;
	*psi = th + *y * .5 / z2;
	for (k = 1; k <= 8; ++k) {
	    i__1 = -k;
	    *psr += a[k - 1] * pow_di(&z2, &i__1) * cos(k * 2. * th);
/* L10: */
	    i__1 = -k;
	    *psi -= a[k - 1] * pow_di(&z2, &i__1) * sin(k * 2. * th);
	}
	if (*x < 8.) {
	    rr = 0.;
	    ri = 0.;
	    i__1 = n;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = x0 - k;
		rr += (x0 - k) / (pow_dd(&d__1, &c_b4) + *y * *y);
/* L20: */
		d__1 = x0 - k;
		ri += *y / (pow_dd(&d__1, &c_b4) + *y * *y);
	    }
	    *psr -= rr;
	    *psi += ri;
	}
	if (x1 < 0.) {
	    tn = tan(pi * *x);
	    tm = tanh(pi * *y);
	    ct2 = tn * tn + tm * tm;
	    *psr = *psr + *x / (*x * *x + *y * *y) + pi * (tn - tn * tm * tm) 
		    / ct2;
	    *psi = *psi - *y / (*x * *x + *y * *y) - pi * tm * (tn * tn + 1.) 
		    / ct2;
	    *x = x1;
	    *y = y1;
	}
    }
    return 0;
} /* cpsi_ */

/*       ********************************** */
/* Subroutine */ int sphy_(integer *n, doublereal *x, integer *nm, doublereal 
	*sy, doublereal *dy)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k;
    static doublereal f0, f1;


/*       ====================================================== */
/*       Purpose: Compute spherical Bessel functions yn(x) and */
/*                their derivatives */
/*       Input :  x --- Argument of yn(x) ( x ≥ 0 ) */
/*                n --- Order of yn(x) ( n = 0,1,… ) */
/*       Output:  SY(n) --- yn(x) */
/*                DY(n) --- yn'(x) */
/*                NM --- Highest order computed */
/*       ====================================================== */

    *nm = *n;
    if (*x < 1e-60) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    sy[k] = -1e300;
/* L10: */
	    dy[k] = 1e300;
	}
	return 0;
    }
    sy[0] = -cos(*x) / *x;
    f0 = sy[0];
    dy[0] = (sin(*x) + cos(*x) / *x) / *x;
    if (*n < 1) {
	return 0;
    }
    sy[1] = (sy[0] - sin(*x)) / *x;
    f1 = sy[1];
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	f = (k * 2. - 1.) * f1 / *x - f0;
	sy[k] = f;
	if (abs(f) >= 1e300) {
	    goto L20;
	}
	f0 = f1;
/* L15: */
	f1 = f;
    }
L20:
    *nm = k - 1;
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	dy[k] = sy[k - 1] - (k + 1.) * sy[k] / *x;
    }
    return 0;
} /* sphy_ */

/*       ********************************** */
/* Subroutine */ int jelp_(doublereal *u, doublereal *hk, doublereal *esn, 
	doublereal *ecn, doublereal *edn, doublereal *eph)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_di(doublereal *, integer *), sin(doublereal),
	     atan(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal a, b, c__, d__;
    static integer j, n;
    static doublereal r__[40], t, a0, b0, dn, sa, pi;


/*       ======================================================== */
/*       Purpose: Compute Jacobian elliptic functions sn u, cn u */
/*                and dn u */
/*       Input  : u   --- Argument of Jacobian elliptic functions */
/*                Hk  --- Modulus k ( 0 ≤ k ≤ 1 ) */
/*       Output : ESN --- sn u */
/*                ECN --- cn u */
/*                EDN --- dn u */
/*                EPH --- phi ( in degrees ) */
/*       ======================================================== */

    pi = 3.14159265358979;
    a0 = 1.;
    b0 = sqrt(1. - *hk * *hk);
    for (n = 1; n <= 40; ++n) {
	a = (a0 + b0) / 2.;
	b = sqrt(a0 * b0);
	c__ = (a0 - b0) / 2.;
	r__[n - 1] = c__ / a;
	if (c__ < 1e-7) {
	    goto L15;
	}
	a0 = a;
/* L10: */
	b0 = b;
    }
L15:
    dn = pow_di(&c_b4, &n) * a * *u;
    d__ = 0.;
    for (j = n; j >= 1; --j) {
	t = r__[j - 1] * sin(dn);
	sa = atan(t / sqrt((d__1 = 1. - t * t, abs(d__1))));
	d__ = (dn + sa) * .5;
/* L20: */
	dn = d__;
    }
    *eph = d__ * 180. / pi;
    *esn = sin(d__);
    *ecn = cos(d__);
    *edn = sqrt(1. - *hk * *hk * *esn * *esn);
    return 0;
} /* jelp_ */

