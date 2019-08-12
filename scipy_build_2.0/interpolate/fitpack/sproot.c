/* sproot.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int sproot_(doublereal *t, integer *n, doublereal *c__, 
	doublereal *zero, integer *mest, integer *m, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal y[3], a0, a1, a2, a3, b0, b1, c1, c2, c3, c4;
    static integer j1;
    static doublereal c5, d4, d5, h1, h2;
    static integer n4;
    static doublereal t1, t2, t3, t4, t5;
    static logical z0, z1, z2, z3, z4;
    static doublereal ah, bh, zz;
    static logical nz0, nz1, nz2, nz3, nz4;
    static doublereal two, three;
    extern /* Subroutine */ int fpcuro_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);

/*  subroutine sproot finds the zeros of a cubic spline s(x),which is */
/*  given in its normalized b-spline representation. */

/*  calling sequence: */
/*     call sproot(t,n,c,zero,mest,m,ier) */

/*  input parameters: */
/*    t    : real array,length n, containing the knots of s(x). */
/*    n    : integer, containing the number of knots.  n>=8 */
/*    c    : real array,length n, containing the b-spline coefficients. */
/*    mest : integer, specifying the dimension of array zero. */

/*  output parameters: */
/*    zero : real array,length mest, containing the zeros of s(x). */
/*    m    : integer,giving the number of zeros. */
/*    ier  : error flag: */
/*      ier = 0: normal return. */
/*      ier = 1: the number of zeros exceeds mest. */
/*      ier =10: invalid input data (see restrictions). */

/*  other subroutines required: fpcuro */

/*  restrictions: */
/*    1) n>= 8. */
/*    2) t(4) < t(5) < ... < t(n-4) < t(n-3). */
/*       t(1) <= t(2) <= t(3) <= t(4) */
/*       t(n-3) <= t(n-2) <= t(n-1) <= t(n) */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/* .. */
/* ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local array.. */
/*  .. */
/*  set some constants */
    /* Parameter adjustments */
    --c__;
    --t;
    --zero;

    /* Function Body */
    two = 2.;
    three = 3.;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    n4 = *n - 4;
    *ier = 10;
    if (*n < 8) {
	goto L800;
    }
    j = *n;
    for (i__ = 1; i__ <= 3; ++i__) {
	if (t[i__] > t[i__ + 1]) {
	    goto L800;
	}
	if (t[j] < t[j - 1]) {
	    goto L800;
	}
	--j;
/* L10: */
    }
    i__1 = n4;
    for (i__ = 4; i__ <= i__1; ++i__) {
	if (t[i__] >= t[i__ + 1]) {
	    goto L800;
	}
/* L20: */
    }
/*  the problem considered reduces to finding the zeros of the cubic */
/*  polynomials pl(x) which define the cubic spline in each knot */
/*  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on */
/*  the condition that it belongs to the knot interval. */
/*  the cubic polynomial pl(x) is determined by computing s(t(l)), */
/*  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute */
/*  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of */
/*  splines and their derivatives, the value of s(t(l)) and s'(t(l)) */
/*  is already known from the foregoing knot interval. */
    *ier = 0;
/*  evaluate some constants for the first knot interval */
    h1 = t[4] - t[3];
    h2 = t[5] - t[4];
    t1 = t[4] - t[2];
    t2 = t[5] - t[3];
    t3 = t[6] - t[4];
    t4 = t[5] - t[2];
    t5 = t[6] - t[3];
/*  calculate a0 = s(t(4)) and ah = s'(t(4)). */
    c1 = c__[1];
    c2 = c__[2];
    c3 = c__[3];
    c4 = (c2 - c1) / t4;
    c5 = (c3 - c2) / t5;
    d4 = (h2 * c1 + t1 * c2) / t4;
    d5 = (t3 * c2 + h1 * c3) / t5;
    a0 = (h2 * d4 + h1 * d5) / t2;
    ah = three * (h2 * c4 + h1 * c5) / t2;
    z1 = TRUE_;
    if (ah < 0.) {
	z1 = FALSE_;
    }
    nz1 = ! z1;
    *m = 0;
/*  main loop for the different knot intervals. */
    i__1 = n4;
    for (l = 4; l <= i__1; ++l) {
/*  evaluate some constants for the knot interval t(l) <= x <= t(l+1). */
	h1 = h2;
	h2 = t[l + 2] - t[l + 1];
	t1 = t2;
	t2 = t3;
	t3 = t[l + 3] - t[l + 1];
	t4 = t5;
	t5 = t[l + 3] - t[l];
/*  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)). */
	c1 = c2;
	c2 = c3;
	c3 = c__[l];
	c4 = c5;
	c5 = (c3 - c2) / t5;
	d4 = (h2 * c1 + t1 * c2) / t4;
	d5 = (h1 * c3 + t3 * c2) / t5;
	b0 = (h2 * d4 + h1 * d5) / t2;
	bh = three * (h2 * c4 + h1 * c5) / t2;
/*  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial */
/*  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)). */
	a1 = ah * h1;
	b1 = bh * h1;
	a2 = three * (b0 - a0) - b1 - two * a1;
	a3 = two * (a0 - b0) + b1 + a1;
/*  test whether or not pl(x) could have a zero in the range */
/*  t(l) <= x <= t(l+1). */
	z3 = TRUE_;
	if (b1 < 0.) {
	    z3 = FALSE_;
	}
	nz3 = ! z3;
	if (a0 * b0 <= 0.) {
	    goto L100;
	}
	z0 = TRUE_;
	if (a0 < 0.) {
	    z0 = FALSE_;
	}
	nz0 = ! z0;
	z2 = TRUE_;
	if (a2 < 0.f) {
	    z2 = FALSE_;
	}
	nz2 = ! z2;
	z4 = TRUE_;
	if (a3 * 3. + a2 < 0.) {
	    z4 = FALSE_;
	}
	nz4 = ! z4;
	if (! (z0 && (nz1 && (z3 || z2 && nz4) || nz2 && z3 && z4) || nz0 && (
		z1 && (nz3 || nz2 && z4) || z2 && nz3 && nz4))) {
	    goto L200;
	}
/*  find the zeros of ql(y). */
L100:
	fpcuro_(&a3, &a2, &a1, &a0, y, &j);
	if (j == 0) {
	    goto L200;
	}
/*  find which zeros of pl(x) are zeros of s(x). */
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (y[i__ - 1] < 0. || y[i__ - 1] > 1.) {
		goto L150;
	    }
/*  test whether the number of zeros of s(x) exceeds mest. */
	    if (*m >= *mest) {
		goto L700;
	    }
	    ++(*m);
	    zero[*m] = t[l] + h1 * y[i__ - 1];
L150:
	    ;
	}
L200:
	a0 = b0;
	ah = bh;
	z1 = z3;
	nz1 = nz3;
/* L300: */
    }
/*  the zeros of s(x) are arranged in increasing order. */
    if (*m < 2) {
	goto L800;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = i__;
L350:
	j1 = j - 1;
	if (j1 == 0) {
	    goto L400;
	}
	if (zero[j] >= zero[j1]) {
	    goto L400;
	}
	zz = zero[j];
	zero[j] = zero[j1];
	zero[j1] = zz;
	j = j1;
	goto L350;
L400:
	;
    }
    j = *m;
    *m = 1;
    i__1 = j;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (zero[i__] == zero[*m]) {
	    goto L500;
	}
	++(*m);
	zero[*m] = zero[i__];
L500:
	;
    }
    goto L800;
L700:
    *ier = 1;
L800:
    return 0;
} /* sproot_ */

