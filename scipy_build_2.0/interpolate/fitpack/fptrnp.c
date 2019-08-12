/* fptrnp.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fptrnp_(integer *m, integer *mm, integer *idim, integer *
	n, integer *nr, doublereal *sp, doublereal *p, doublereal *b, 
	doublereal *z__, doublereal *a, doublereal *q, doublereal *right)
{
    /* System generated locals */
    integer sp_dim1, sp_offset, b_dim1, b_offset, a_dim1, a_offset, i__1, 
	    i__2, i__3, i__4;

    /* Local variables */
    static doublereal h__[7];
    static integer i__, j, l, i2, i3, m2, m3, n1, n4, ii, jj, it, mid, nmd;
    static doublereal one, cos__, sin__, piv, pinv;
    static integer irot, iband, nrold, number;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*  subroutine fptrnp reduces the (m+n-7) x (n-4) matrix a to upper */
/*  triangular form and applies the same givens transformations to */
/*  the (m) x (mm) x (idim) matrix z to obtain the (n-4) x (mm) x */
/*  (idim) matrix q */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..subroutine references.. */
/*    fpgivs,fprota */
/*  .. */
    /* Parameter adjustments */
    sp_dim1 = *m;
    sp_offset = 1 + sp_dim1;
    sp -= sp_offset;
    --nr;
    --right;
    --z__;
    --q;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    one = 1.;
    if (*p > 0.f) {
	pinv = one / *p;
    }
    n4 = *n - 4;
    mid = *mm * *idim;
    m2 = *m * *mm;
    m3 = n4 * *mm;
/*  reduce the matrix (a) to upper triangular form (r) using givens */
/*  rotations. apply the same transformations to the rows of matrix z */
/*  to obtain the mm x (n-4) matrix g. */
/*  store matrix (r) into (a) and g into q. */
/*  initialization. */
    nmd = n4 * mid;
    i__1 = nmd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = 0.f;
/* L50: */
    }
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 5; ++j) {
	    a[i__ + j * a_dim1] = 0.f;
/* L100: */
	}
    }
    nrold = 0;
/*  iband denotes the bandwidth of the matrices (a) and (r). */
    iband = 4;
    i__1 = *m;
    for (it = 1; it <= i__1; ++it) {
	number = nr[it];
L150:
	if (nrold == number) {
	    goto L300;
	}
	if (*p <= 0.f) {
	    goto L700;
	}
	iband = 5;
/*  fetch a new row of matrix (b). */
	n1 = nrold + 1;
	for (j = 1; j <= 5; ++j) {
	    h__[j - 1] = b[n1 + j * b_dim1] * pinv;
/* L200: */
	}
/*  find the appropriate column of q. */
	i__2 = mid;
	for (j = 1; j <= i__2; ++j) {
	    right[j] = 0.f;
/* L250: */
	}
	irot = nrold;
	goto L450;
/*  fetch a new row of matrix (sp). */
L300:
	h__[iband - 1] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    h__[j - 1] = sp[it + j * sp_dim1];
/* L350: */
	}
/*  find the appropriate column of q. */
	j = 0;
	i__2 = *idim;
	for (ii = 1; ii <= i__2; ++ii) {
	    l = (ii - 1) * m2 + (it - 1) * *mm;
	    i__3 = *mm;
	    for (jj = 1; jj <= i__3; ++jj) {
		++j;
		++l;
		right[j] = z__[l];
/* L400: */
	    }
	}
	irot = number;
/*  rotate the new row of matrix (a) into triangle. */
L450:
	i__3 = iband;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ++irot;
	    piv = h__[i__ - 1];
	    if (piv == 0.f) {
		goto L600;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &a[irot + a_dim1], &cos__, &sin__);
/*  apply that transformation to the rows of matrix q. */
	    j = 0;
	    i__2 = *idim;
	    for (ii = 1; ii <= i__2; ++ii) {
		l = (ii - 1) * m3 + irot;
		i__4 = *mm;
		for (jj = 1; jj <= i__4; ++jj) {
		    ++j;
		    fprota_(&cos__, &sin__, &right[j], &q[l]);
		    l += n4;
/* L500: */
		}
	    }
/*  apply that transformation to the columns of (a). */
	    if (i__ == iband) {
		goto L650;
	    }
	    i2 = 1;
	    i3 = i__ + 1;
	    i__4 = iband;
	    for (j = i3; j <= i__4; ++j) {
		++i2;
		fprota_(&cos__, &sin__, &h__[j - 1], &a[irot + i2 * a_dim1]);
/* L550: */
	    }
L600:
	    ;
	}
L650:
	if (nrold == number) {
	    goto L750;
	}
L700:
	++nrold;
	goto L150;
L750:
	;
    }
    return 0;
} /* fptrnp_ */

