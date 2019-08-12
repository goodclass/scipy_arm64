/* fpcoco.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpcoco_(integer *iopt, integer *m, doublereal *x, 
	doublereal *y, doublereal *w, doublereal *v, doublereal *s, integer *
	nest, integer *maxtr, integer *maxbin, integer *n, doublereal *t, 
	doublereal *c__, doublereal *sq, doublereal *sx, logical *bind, 
	doublereal *e, doublereal *wrk, integer *lwrk, integer *iwrk, integer 
	*kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, i1, l1, m1, n4, n6, n8, ia, ib, ic, mb, ji, 
	    jl, iq, nm, jr, iu, ju, nr;
    static doublereal tj, xi;
    static integer iz, jib, jjb;
    static doublereal sql;
    static integer izz;
    static doublereal half;
    static integer nmax;
    static doublereal term, sqmax;
    extern /* Subroutine */ int fpcosp_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, logical *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *)
	    ;

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..subroutine references.. */
/*    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno */
/*  .. */
/*  set constant */
    /* Parameter adjustments */
    --sx;
    --v;
    --w;
    --y;
    --x;
    --e;
    --bind;
    --c__;
    --t;
    --wrk;
    --iwrk;

    /* Function Body */
    half = .5f;
/*  determine the maximal admissible number of knots. */
    nmax = *m + 4;
/*  the initial choice of knots depends on the value of iopt. */
/*    if iopt=0 the program starts with the minimal number of knots */
/*    so that can be guarantied that the concavity/convexity constraints */
/*    will be satisfied. */
/*    if iopt = 1 the program will continue from the point on where she */
/*    left at the foregoing call. */
    if (*iopt > 0) {
	goto L80;
    }
/*  find the minimal number of knots. */
/*  a knot is located at the data point x(i), i=2,3,...m-1 if */
/*    1) v(i) ^= 0    and */
/*    2) v(i)*v(i-1) <= 0  or  v(i)*v(i+1) <= 0. */
    m1 = *m - 1;
    *n = 4;
    i__1 = m1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (v[i__] == 0.f || v[i__] * v[i__ - 1] > 0.f && v[i__] * v[i__ + 1] 
		> 0.f) {
	    goto L20;
	}
	++(*n);
/*  test whether the required storage space exceeds the available one. */
	if (*n + 4 > *nest) {
	    goto L200;
	}
	t[*n] = x[i__];
L20:
	;
    }
/*  find the position of the knots t(1),...t(4) and t(n-3),...t(n) which */
/*  are needed for the b-spline representation of s(x). */
    for (i__ = 1; i__ <= 4; ++i__) {
	t[i__] = x[1];
	++(*n);
	t[*n] = x[*m];
/* L30: */
    }
/*  test whether the minimum number of knots exceeds the maximum number. */
    if (*n > nmax) {
	goto L210;
    }
/*  main loop for the different sets of knots. */
/*  find corresponding values e(j) to the knots t(j+3),j=1,2,...n-6 */
/*    e(j) will take the value -1,1, or 0 according to the requirement */
/*    that s(x) must be locally convex or concave at t(j+3) or that the */
/*    sign of s''(x) is unrestricted at that point. */
L40:
    i__ = 1;
    xi = x[1];
    j = 4;
    tj = t[4];
    n6 = *n - 6;
    i__1 = n6;
    for (l = 1; l <= i__1; ++l) {
L50:
	if (xi == tj) {
	    goto L60;
	}
	++i__;
	xi = x[i__];
	goto L50;
L60:
	e[l] = v[i__];
	++j;
	tj = t[j];
/* L70: */
    }
/*  we partition the working space */
    nm = *n + *maxbin;
    mb = *maxbin + 1;
    ia = 1;
    ib = ia + (*n << 2);
    ic = ib + nm * *maxbin;
    iz = ic + *n;
    izz = iz + *n;
    iu = izz + *n;
    iq = iu + *maxbin;
    ji = 1;
    ju = ji + *maxtr;
    jl = ju + *maxtr;
    jr = jl + *maxtr;
    jjb = jr + *maxtr;
    jib = jjb + mb;
/*  given the set of knots t(j),j=1,2,...n, find the least-squares cubic */
/*  spline which satisfies the imposed concavity/convexity constraints. */
    fpcosp_(m, &x[1], &y[1], &w[1], n, &t[1], &e[1], maxtr, maxbin, &c__[1], 
	    sq, &sx[1], &bind[1], &nm, &mb, &wrk[ia], &wrk[ib], &wrk[ic], &
	    wrk[iz], &wrk[izz], &wrk[iu], &wrk[iq], &iwrk[ji], &iwrk[ju], &
	    iwrk[jl], &iwrk[jr], &iwrk[jjb], &iwrk[jib], ier);
/*  if sq <= s or in case of abnormal exit from fpcosp, control is */
/*  repassed to the driver program. */
    if (*sq <= *s || *ier > 0) {
	goto L300;
    }
/*  calculate for each knot interval t(l-1) <= xi <= t(l) the */
/*  sum((wi*(yi-s(xi)))**2). */
/*  find the interval t(k-1) <= x <= t(k) for which this sum is maximal */
/*  on the condition that this interval contains at least one interior */
/*  data point x(nr) and that s(x) is not given there by a straight line. */
L80:
    sqmax = 0.f;
    sql = 0.f;
    l = 5;
    nr = 0;
    i1 = 1;
    n4 = *n - 4;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = w[i__] * (sx[i__] - y[i__]);
	term = d__1 * d__1;
	if (x[i__] < t[l] || l > n4) {
	    goto L100;
	}
	term *= half;
	sql += term;
	if (i__ - i1 <= 1 || bind[l - 4] && bind[l - 3]) {
	    goto L90;
	}
	if (sql <= sqmax) {
	    goto L90;
	}
	k = l;
	sqmax = sql;
	nr = i1 + (i__ - i1) / 2;
L90:
	++l;
	i1 = i__;
	sql = 0.f;
L100:
	sql += term;
/* L110: */
    }
    if (*m - i1 <= 1 || bind[l - 4] && bind[l - 3]) {
	goto L120;
    }
    if (sql <= sqmax) {
	goto L120;
    }
    k = l;
    nr = i1 + (*m - i1) / 2;
/*  if no such interval is found, control is repassed to the driver */
/*  program (ier = -1). */
L120:
    if (nr == 0) {
	goto L190;
    }
/*  if s(x) is given by the same straight line in two succeeding knot */
/*  intervals t(l-1) <= x <= t(l) and t(l) <= x <= t(l+1),delete t(l) */
    n8 = *n - 8;
    l1 = 0;
    if (n8 <= 0) {
	goto L150;
    }
    i__1 = n8;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (! (bind[i__] && bind[i__ + 1] && bind[i__ + 2])) {
	    goto L140;
	}
	l = i__ + 4 - l1;
	if (k > l) {
	    --k;
	}
	--(*n);
	++l1;
	i__2 = *n;
	for (j = l; j <= i__2; ++j) {
	    t[j] = t[j + 1];
/* L130: */
	}
L140:
	;
    }
/*  test whether we cannot further increase the number of knots. */
L150:
    if (*n == nmax) {
	goto L180;
    }
    if (*n == *nest) {
	goto L170;
    }
/*  locate an additional knot at the point x(nr). */
    j = *n;
    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) {
	t[j + 1] = t[j];
	--j;
/* L160: */
    }
    t[k] = x[nr];
    ++(*n);
/*  restart the computations with the new set of knots. */
    goto L40;
/*  error codes and messages. */
L170:
    *ier = -3;
    goto L300;
L180:
    *ier = -2;
    goto L300;
L190:
    *ier = -1;
    goto L300;
L200:
    *ier = 4;
    goto L300;
L210:
    *ier = 5;
L300:
    return 0;
} /* fpcoco_ */

