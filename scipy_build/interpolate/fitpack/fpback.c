/* fpback.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpback_(doublereal *a, doublereal *z__, integer *n, 
	integer *k, doublereal *c__, integer *nest)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, m, i1, k1;
    static doublereal store;

/*  subroutine fpback calculates the solution of the system of */
/*  equations a*c = z with a a n x n upper triangular matrix */
/*  of bandwidth k. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --c__;
    --z__;
    a_dim1 = *nest;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    k1 = *k - 1;
    c__[*n] = z__[*n] / a[*n + a_dim1];
    i__ = *n - 1;
    if (i__ == 0) {
	goto L30;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	store = z__[i__];
	i1 = k1;
	if (j <= k1) {
	    i1 = j - 1;
	}
	m = i__;
	i__2 = i1;
	for (l = 1; l <= i__2; ++l) {
	    ++m;
	    store -= c__[m] * a[i__ + (l + 1) * a_dim1];
/* L10: */
	}
	c__[i__] = store / a[i__ + a_dim1];
	--i__;
/* L20: */
    }
L30:
    return 0;
} /* fpback_ */

