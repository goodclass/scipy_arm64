/* fppocu.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fppocu_(integer *idim, integer *k, doublereal *a, 
	doublereal *b, integer *ib, doublereal *db, integer *nb, integer *ie, 
	doublereal *de, integer *ne, doublereal *cp, integer *np)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, k1, k2;
    static doublereal ab;
    static integer id, jj, ll;
    static doublereal aki, work[36]	/* was [6][6] */;

/*  subroutine fppocu finds a idim-dimensional polynomial curve p(u) = */
/*  (p1(u),p2(u),...,pidim(u)) of degree k, satisfying certain derivative */
/*  constraints at the end points a and b, i.e. */
/*                  (l) */
/*    if ib > 0 : pj   (a) = db(idim*l+j), l=0,1,...,ib-1 */
/*                  (l) */
/*    if ie > 0 : pj   (b) = de(idim*l+j), l=0,1,...,ie-1 */

/*  the polynomial curve is returned in its b-spline representation */
/*  ( coefficients cp(j), j=1,2,...,np ) */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local array.. */
/*  .. */
    /* Parameter adjustments */
    --db;
    --de;
    --cp;

    /* Function Body */
    k1 = *k + 1;
    k2 = k1 << 1;
    ab = *b - *a;
    i__1 = *idim;
    for (id = 1; id <= i__1; ++id) {
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
	    work[j - 1] = 0.f;
/* L10: */
	}
	if (*ib == 0) {
	    goto L50;
	}
	l = id;
	i__2 = *ib;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__ * 6 - 6] = db[l];
	    l += *idim;
/* L20: */
	}
	if (*ib == 1) {
	    goto L50;
	}
	ll = *ib;
	i__2 = *ib;
	for (j = 2; j <= i__2; ++j) {
	    --ll;
	    i__3 = ll;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		aki = (doublereal) (k1 - i__);
		work[j + i__ * 6 - 7] = ab * work[j - 1 + (i__ + 1) * 6 - 7] /
			 aki + work[j - 1 + i__ * 6 - 7];
/* L30: */
	    }
/* L40: */
	}
L50:
	if (*ie == 0) {
	    goto L90;
	}
	l = id;
	j = k1;
	i__2 = *ie;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[j + i__ * 6 - 7] = de[l];
	    l += *idim;
	    --j;
/* L60: */
	}
	if (*ie == 1) {
	    goto L90;
	}
	ll = *ie;
	i__2 = *ie;
	for (jj = 2; jj <= i__2; ++jj) {
	    --ll;
	    j = k1 + 1 - jj;
	    i__3 = ll;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		aki = (doublereal) (k1 - i__);
		work[j + i__ * 6 - 7] = work[j + 1 + i__ * 6 - 7] - ab * work[
			j + (i__ + 1) * 6 - 7] / aki;
		--j;
/* L70: */
	    }
/* L80: */
	}
L90:
	l = (id - 1) * k2;
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
	    ++l;
	    cp[l] = work[j - 1];
/* L100: */
	}
/* L110: */
    }
    return 0;
} /* fppocu_ */

