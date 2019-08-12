/* fpchep.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpchep_(doublereal *x, integer *m, doublereal *t, 
	integer *n, integer *k, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, i1, i2, j1, k1, k2, l1, l2, m1, mm;
    static doublereal tj, tl, xi;
    static integer nk1, nk2;
    static doublereal per;

/*  subroutine fpchep verifies the number and the position of the knots */
/*  t(j),j=1,2,...,n of a periodic spline of degree k, in relation to */
/*  the number and the position of the data points x(i),i=1,2,...,m. */
/*  if all of the following conditions are fulfilled, ier is set */
/*  to zero. if one of the conditions is violated ier is set to ten. */
/*      1) k+1 <= n-k-1 <= m+k-1 */
/*      2) t(1) <= t(2) <= ... <= t(k+1) */
/*         t(n-k) <= t(n-k+1) <= ... <= t(n) */
/*      3) t(k+1) < t(k+2) < ... < t(n-k) */
/*      4) t(k+1) <= x(i) <= t(n-k) */
/*      5) the conditions specified by schoenberg and whitney must hold */
/*         for at least one subset of data points, i.e. there must be a */
/*         subset of data points y(j) such that */
/*             t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1 */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --x;
    --t;

    /* Function Body */
    k1 = *k + 1;
    k2 = k1 + 1;
    nk1 = *n - k1;
    nk2 = nk1 + 1;
    m1 = *m - 1;
    *ier = 10;
/*  check condition no 1 */
    if (nk1 < k1 || *n > *m + (*k << 1)) {
	goto L130;
    }
/*  check condition no 2 */
    j = *n;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (t[i__] > t[i__ + 1]) {
	    goto L130;
	}
	if (t[j] < t[j - 1]) {
	    goto L130;
	}
	--j;
/* L20: */
    }
/*  check condition no 3 */
    i__1 = nk2;
    for (i__ = k2; i__ <= i__1; ++i__) {
	if (t[i__] <= t[i__ - 1]) {
	    goto L130;
	}
/* L30: */
    }
/*  check condition no 4 */
    if (x[1] < t[k1] || x[*m] > t[nk2]) {
	goto L130;
    }
/*  check condition no 5 */
    l1 = k1;
    l2 = 1;
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
	xi = x[l];
L40:
	if (xi < t[l1 + 1] || l == nk1) {
	    goto L50;
	}
	++l1;
	++l2;
	if (l2 > k1) {
	    goto L60;
	}
	goto L40;
L50:
	;
    }
    l = *m;
L60:
    per = t[nk2] - t[k1];
    i__1 = l;
    for (i1 = 2; i1 <= i__1; ++i1) {
	i__ = i1 - 1;
	mm = i__ + m1;
	i__2 = nk1;
	for (j = k1; j <= i__2; ++j) {
	    tj = t[j];
	    j1 = j + k1;
	    tl = t[j1];
L70:
	    ++i__;
	    if (i__ > mm) {
		goto L120;
	    }
	    i2 = i__ - m1;
	    if (i2 <= 0) {
		goto L80;
	    }
	    goto L90;
L80:
	    xi = x[i__];
	    goto L100;
L90:
	    xi = x[i2] + per;
L100:
	    if (xi <= tj) {
		goto L70;
	    }
	    if (xi >= tl) {
		goto L120;
	    }
/* L110: */
	}
	*ier = 0;
	goto L130;
L120:
	;
    }
L130:
    return 0;
} /* fpchep_ */

