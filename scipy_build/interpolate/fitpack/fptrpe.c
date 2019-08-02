/* fptrpe.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fptrpe_(integer *m, integer *mm, integer *idim, integer *
	n, integer *nr, doublereal *sp, doublereal *p, doublereal *b, 
	doublereal *z__, doublereal *a, doublereal *aa, doublereal *q, 
	doublereal *right)
{
    /* System generated locals */
    integer sp_dim1, sp_offset, b_dim1, b_offset, a_dim1, a_offset, aa_dim1, 
	    aa_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal h__[5];
    static integer i__, j, l;
    static doublereal h1[5], h2[4];
    static integer i2, i3, l0, m1, m2, m3, n1, n4, l1, i1, n7, j1, n11;
    static doublereal co;
    static integer ii, jj, jk, ik, ij;
    static doublereal si;
    static integer it, mid, nmd;
    static doublereal one, piv;
    static integer jper;
    static doublereal pinv;
    static integer irot, nrold, number;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*  subroutine fptrpe reduces the (m+n-7) x (n-7) cyclic bandmatrix a */
/*  to upper triangular form and applies the same givens transformations */
/*  to the (m) x (mm) x (idim) matrix z to obtain the (n-7) x (mm) x */
/*  (idim) matrix q. */
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
    aa_dim1 = *n;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
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
    n7 = *n - 7;
    n11 = *n - 11;
    mid = *mm * *idim;
    m2 = *m * *mm;
    m3 = n7 * *mm;
    m1 = *m - 1;
/*  we determine the matrix (a) and then we reduce her to */
/*  upper triangular form (r) using givens rotations. */
/*  we apply the same transformations to the rows of matrix */
/*  z to obtain the (mm) x (n-7) matrix g. */
/*  we store matrix (r) into a and aa, g into q. */
/*  the n7 x n7 upper triangular matrix (r) has the form */
/*             | a1 '     | */
/*       (r) = |    ' a2  | */
/*             |  0 '     | */
/*  with (a2) a n7 x 4 matrix and (a1) a n11 x n11 upper */
/*  triangular matrix of bandwidth 5. */
/*  initialization. */
    nmd = n7 * mid;
    i__1 = nmd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = 0.f;
/* L50: */
    }
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + a_dim1 * 5] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    a[i__ + j * a_dim1] = 0.f;
	    aa[i__ + j * aa_dim1] = 0.f;
/* L100: */
	}
    }
    jper = 0;
    nrold = 0;
    i__1 = m1;
    for (it = 1; it <= i__1; ++it) {
	number = nr[it];
L120:
	if (nrold == number) {
	    goto L180;
	}
	if (*p <= 0.f) {
	    goto L740;
	}
/*  fetch a new row of matrix (b). */
	n1 = nrold + 1;
	for (j = 1; j <= 5; ++j) {
	    h__[j - 1] = b[n1 + j * b_dim1] * pinv;
/* L140: */
	}
/*  find the appropriate row of q. */
	i__2 = mid;
	for (j = 1; j <= i__2; ++j) {
	    right[j] = 0.f;
/* L160: */
	}
	goto L240;
/*  fetch a new row of matrix (sp) */
L180:
	h__[4] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    h__[j - 1] = sp[it + j * sp_dim1];
/* L200: */
	}
/*  find the appropriate row of q. */
	j = 0;
	i__2 = *idim;
	for (ii = 1; ii <= i__2; ++ii) {
	    l = (ii - 1) * m2 + (it - 1) * *mm;
	    i__3 = *mm;
	    for (jj = 1; jj <= i__3; ++jj) {
		++j;
		++l;
		right[j] = z__[l];
/* L220: */
	    }
	}
/*  test whether there are non-zero values in the new row of (a) */
/*  corresponding to the b-splines n(j,*),j=n7+1,...,n4. */
L240:
	if (nrold < n11) {
	    goto L640;
	}
	if (jper != 0) {
	    goto L320;
	}
/*  initialize the matrix (aa). */
	jk = n11 + 1;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ik = jk;
	    for (j = 1; j <= 5; ++j) {
		if (ik <= 0) {
		    goto L280;
		}
		aa[ik + i__ * aa_dim1] = a[ik + j * a_dim1];
		--ik;
/* L260: */
	    }
L280:
	    ++jk;
/* L300: */
	}
	jper = 1;
/*  if one of the non-zero elements of the new row corresponds to one of */
/*  the b-splines n(j;*),j=n7+1,...,n4,we take account of the periodicity */
/*  conditions for setting up this row of (a). */
L320:
	for (i__ = 1; i__ <= 4; ++i__) {
	    h1[i__ - 1] = 0.f;
	    h2[i__ - 1] = 0.f;
/* L340: */
	}
	h1[4] = 0.f;
	j = nrold - n11;
	for (i__ = 1; i__ <= 5; ++i__) {
	    ++j;
	    l0 = j;
L360:
	    l1 = l0 - 4;
	    if (l1 <= 0) {
		goto L400;
	    }
	    if (l1 <= n11) {
		goto L380;
	    }
	    l0 = l1 - n11;
	    goto L360;
L380:
	    h1[l1 - 1] = h__[i__ - 1];
	    goto L420;
L400:
	    h2[l0 - 1] += h__[i__ - 1];
L420:
	    ;
	}
/*  rotate the new row of (a) into triangle. */
	if (n11 <= 0) {
	    goto L560;
	}
/*  rotations with the rows 1,2,...,n11 of (a). */
	i__3 = n11;
	for (irot = 1; irot <= i__3; ++irot) {
	    piv = h1[0];
/* Computing MIN */
	    i__2 = n11 - irot;
	    i2 = min(i__2,4);
	    if (piv == 0.f) {
		goto L500;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the columns of matrix q. */
	    j = 0;
	    i__2 = *idim;
	    for (ii = 1; ii <= i__2; ++ii) {
		l = (ii - 1) * m3 + irot;
		i__4 = *mm;
		for (jj = 1; jj <= i__4; ++jj) {
		    ++j;
		    fprota_(&co, &si, &right[j], &q[l]);
		    l += n7;
/* L440: */
		}
	    }
/*  apply that transformation to the rows of (a) with respect to aa. */
	    for (i__ = 1; i__ <= 4; ++i__) {
		fprota_(&co, &si, &h2[i__ - 1], &aa[irot + i__ * aa_dim1]);
/* L460: */
	    }
/*  apply that transformation to the rows of (a) with respect to a. */
	    if (i2 == 0) {
		goto L560;
	    }
	    i__4 = i2;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i1 = i__ + 1;
		fprota_(&co, &si, &h1[i1 - 1], &a[irot + i1 * a_dim1]);
/* L480: */
	    }
L500:
	    i__4 = i2;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		h1[i__ - 1] = h1[i__];
/* L520: */
	    }
	    h1[i2] = 0.f;
/* L540: */
	}
/*  rotations with the rows n11+1,...,n7 of a. */
L560:
	for (irot = 1; irot <= 4; ++irot) {
	    ij = n11 + irot;
	    if (ij <= 0) {
		goto L620;
	    }
	    piv = h2[irot - 1];
	    if (piv == 0.f) {
		goto L620;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &aa[ij + irot * aa_dim1], &co, &si);
/*  apply that transformation to the columns of matrix q. */
	    j = 0;
	    i__3 = *idim;
	    for (ii = 1; ii <= i__3; ++ii) {
		l = (ii - 1) * m3 + ij;
		i__4 = *mm;
		for (jj = 1; jj <= i__4; ++jj) {
		    ++j;
		    fprota_(&co, &si, &right[j], &q[l]);
		    l += n7;
/* L580: */
		}
	    }
	    if (irot == 4) {
		goto L620;
	    }
/*  apply that transformation to the rows of (a) with respect to aa. */
	    j1 = irot + 1;
	    for (i__ = j1; i__ <= 4; ++i__) {
		fprota_(&co, &si, &h2[i__ - 1], &aa[ij + i__ * aa_dim1]);
/* L600: */
	    }
L620:
	    ;
	}
	goto L720;
/*  rotation into triangle of the new row of (a), in case the elements */
/*  corresponding to the b-splines n(j;*),j=n7+1,...,n4 are all zero. */
L640:
	irot = nrold;
	for (i__ = 1; i__ <= 5; ++i__) {
	    ++irot;
	    piv = h__[i__ - 1];
	    if (piv == 0.f) {
		goto L700;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
	    j = 0;
	    i__4 = *idim;
	    for (ii = 1; ii <= i__4; ++ii) {
		l = (ii - 1) * m3 + irot;
		i__3 = *mm;
		for (jj = 1; jj <= i__3; ++jj) {
		    ++j;
		    fprota_(&co, &si, &right[j], &q[l]);
		    l += n7;
/* L660: */
		}
	    }
/*  apply that transformation to the rows of (a). */
	    if (i__ == 5) {
		goto L700;
	    }
	    i2 = 1;
	    i3 = i__ + 1;
	    for (j = i3; j <= 5; ++j) {
		++i2;
		fprota_(&co, &si, &h__[j - 1], &a[irot + i2 * a_dim1]);
/* L680: */
	    }
L700:
	    ;
	}
L720:
	if (nrold == number) {
	    goto L760;
	}
L740:
	++nrold;
	goto L120;
L760:
	;
    }
    return 0;
} /* fptrpe_ */

