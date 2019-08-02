/* basym.f -- translated by f2c (version 20190311).
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

doublereal basym_(doublereal *a, doublereal *b, doublereal *lambda, 
	doublereal *eps)
{
    /* Initialized data */

    static integer num = 20;
    static doublereal e0 = 1.12837916709551;
    static doublereal e1 = .353553390593274;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal c__[21], d__[21], f, h__;
    static integer i__, j, m, n;
    static doublereal r__, s, t, u, w, z__, a0[21], b0[21], h2, j0, j1, r0, 
	    r1, t0, t1, w0, z0, z2, hn, zn;
    static integer im1, mm1, np1, imj, mmj;
    static doublereal sum, znm1, bsum, dsum;
    extern doublereal erfc1_(integer *, doublereal *), rlog1_(doublereal *), 
	    bcorr_(doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B. */
/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/* ------------------------ */
/*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

/* ------------------------ */
/*     E0 = 2/SQRT(PI) */
/*     E1 = 2**(-3/2) */
/* ------------------------ */
/*     .. */
/*     .. Executable Statements .. */
/* ------------------------ */
    ret_val = 0.;
    if (*a >= *b) {
	goto L10;
    }
    h__ = *a / *b;
    r0 = 1. / (h__ + 1.);
    r1 = (*b - *a) / *b;
    w0 = 1. / sqrt(*a * (h__ + 1.));
    goto L20;
L10:
    h__ = *b / *a;
    r0 = 1. / (h__ + 1.);
    r1 = (*b - *a) / *a;
    w0 = 1. / sqrt(*b * (h__ + 1.));

L20:
    d__1 = -(*lambda) / *a;
    d__2 = *lambda / *b;
    f = *a * rlog1_(&d__1) + *b * rlog1_(&d__2);
    t = exp(-f);
    if (t == 0.) {
	return ret_val;
    }
    z0 = sqrt(f);
    z__ = z0 / e1 * .5;
    z2 = f + f;

    a0[0] = r1 * .66666666666666663;
    c__[0] = a0[0] * -.5;
    d__[0] = -c__[0];
    j0 = .5 / e0 * erfc1_(&c__1, &z0);
    j1 = e1;
    sum = j0 + d__[0] * w0 * j1;

    s = 1.;
    h2 = h__ * h__;
    hn = 1.;
    w = w0;
    znm1 = z__;
    zn = z2;
    i__1 = num;
    for (n = 2; n <= i__1; n += 2) {
	hn = h2 * hn;
	a0[n - 1] = r0 * 2. * (h__ * hn + 1.) / (n + 2.);
	np1 = n + 1;
	s += hn;
	a0[np1 - 1] = r1 * 2. * s / (n + 3.);

	i__2 = np1;
	for (i__ = n; i__ <= i__2; ++i__) {
	    r__ = (i__ + 1.) * -.5;
	    b0[0] = r__ * a0[0];
	    i__3 = i__;
	    for (m = 2; m <= i__3; ++m) {
		bsum = 0.;
		mm1 = m - 1;
		i__4 = mm1;
		for (j = 1; j <= i__4; ++j) {
		    mmj = m - j;
		    bsum += (j * r__ - mmj) * a0[j - 1] * b0[mmj - 1];
/* L30: */
		}
		b0[m - 1] = r__ * a0[m - 1] + bsum / m;
/* L40: */
	    }
	    c__[i__ - 1] = b0[i__ - 1] / (i__ + 1.);

	    dsum = 0.;
	    im1 = i__ - 1;
	    i__3 = im1;
	    for (j = 1; j <= i__3; ++j) {
		imj = i__ - j;
		dsum += d__[imj - 1] * c__[j - 1];
/* L50: */
	    }
	    d__[i__ - 1] = -(dsum + c__[i__ - 1]);
/* L60: */
	}

	j0 = e1 * znm1 + (n - 1.) * j0;
	j1 = e1 * zn + n * j1;
	znm1 = z2 * znm1;
	zn = z2 * zn;
	w = w0 * w;
	t0 = d__[n - 1] * w * j0;
	w = w0 * w;
	t1 = d__[np1 - 1] * w * j1;
	sum += t0 + t1;
	if (abs(t0) + abs(t1) <= *eps * sum) {
	    goto L80;
	}
/* L70: */
    }

L80:
    u = exp(-bcorr_(a, b));
    ret_val = e0 * t * u * sum;
    return ret_val;
} /* basym_ */

