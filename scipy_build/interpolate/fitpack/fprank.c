/* fprank.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fprank_(doublereal *a, doublereal *f, integer *n, 
	integer *m, integer *na, doublereal *tol, doublereal *c__, doublereal 
	*sq, integer *rank, doublereal *aa, doublereal *ff, doublereal *h__)
{
    /* System generated locals */
    integer a_dim1, a_offset, aa_dim1, aa_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i1, i2, j1, j2, j3, m1, ii, ij, jj, kk, nl;
    static doublereal yi, fac, cos__, sin__, piv, stor1, stor2, stor3, store;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*  subroutine fprank finds the minimum norm solution of a least- */
/*  squares problem in case of rank deficiency. */

/*  input parameters: */
/*    a : array, which contains the non-zero elements of the observation */
/*        matrix after triangularization by givens transformations. */
/*    f : array, which contains the transformed right hand side. */
/*    n : integer,which contains the dimension of a. */
/*    m : integer, which denotes the bandwidth of a. */
/*  tol : real value, giving a threshold to determine the rank of a. */

/*  output parameters: */
/*    c : array, which contains the minimum norm solution. */
/*   sq : real value, giving the contribution of reducing the rank */
/*        to the sum of squared residuals. */
/* rank : integer, which contains the rank of matrix a. */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpgivs,fprota */
/*  .. */
    /* Parameter adjustments */
    --ff;
    --c__;
    --f;
    --h__;
    aa_dim1 = *n;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    a_dim1 = *na;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    m1 = *m - 1;
/*  the rank deficiency nl is considered to be the number of sufficient */
/*  small diagonal elements of a. */
    nl = 0;
    *sq = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[i__ + a_dim1] > *tol) {
	    goto L90;
	}
/*  if a sufficient small diagonal element is found, we put it to */
/*  zero. the remainder of the row corresponding to that zero diagonal */
/*  element is then rotated into triangle by givens rotations . */
/*  the rank deficiency is increased by one. */
	++nl;
	if (i__ == *n) {
	    goto L90;
	}
	yi = f[i__];
	i__2 = m1;
	for (j = 1; j <= i__2; ++j) {
	    h__[j] = a[i__ + (j + 1) * a_dim1];
/* L10: */
	}
	h__[*m] = 0.f;
	i1 = i__ + 1;
	i__2 = *n;
	for (ii = i1; ii <= i__2; ++ii) {
/* Computing MIN */
	    i__3 = *n - ii;
	    i2 = min(i__3,m1);
	    piv = h__[1];
	    if (piv == 0.f) {
		goto L30;
	    }
	    fpgivs_(&piv, &a[ii + a_dim1], &cos__, &sin__);
	    fprota_(&cos__, &sin__, &yi, &f[ii]);
	    if (i2 == 0) {
		goto L70;
	    }
	    i__3 = i2;
	    for (j = 1; j <= i__3; ++j) {
		j1 = j + 1;
		fprota_(&cos__, &sin__, &h__[j1], &a[ii + j1 * a_dim1]);
		h__[j] = h__[j1];
/* L20: */
	    }
	    goto L50;
L30:
	    if (i2 == 0) {
		goto L70;
	    }
	    i__3 = i2;
	    for (j = 1; j <= i__3; ++j) {
		h__[j] = h__[j + 1];
/* L40: */
	    }
L50:
	    h__[i2 + 1] = 0.f;
/* L60: */
	}
/*  add to the sum of squared residuals the contribution of deleting */
/*  the row with small diagonal element. */
L70:
/* Computing 2nd power */
	d__1 = yi;
	*sq += d__1 * d__1;
L90:
	;
    }
/*  rank denotes the rank of a. */
    *rank = *n - nl;
/*  let b denote the (rank*n) upper trapezoidal matrix which can be */
/*  obtained from the (n*n) upper triangular matrix a by deleting */
/*  the rows and interchanging the columns corresponding to a zero */
/*  diagonal element. if this matrix is factorized using givens */
/*  transformations as  b = (r) (u)  where */
/*    r is a (rank*rank) upper triangular matrix, */
/*    u is a (rank*n) orthonormal matrix */
/*  then the minimal least-squares solution c is given by c = b' v, */
/*  where v is the solution of the system  (r) (r)' v = g  and */
/*  g denotes the vector obtained from the old right hand side f, by */
/*  removing the elements corresponding to a zero diagonal element of a. */
/*  initialization. */
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    aa[i__ + j * aa_dim1] = 0.f;
/* L100: */
	}
    }
/*  form in aa the upper triangular matrix obtained from a by */
/*  removing rows and columns with zero diagonal elements. form in ff */
/*  the new right hand side by removing the elements of the old right */
/*  hand side corresponding to a deleted row. */
    ii = 0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (a[i__ + a_dim1] <= *tol) {
	    goto L120;
	}
	++ii;
	ff[ii] = f[i__];
	aa[ii + aa_dim1] = a[i__ + a_dim1];
	jj = ii;
	kk = 1;
	j = i__;
/* Computing MIN */
	i__1 = j - 1;
	j1 = min(i__1,m1);
	if (j1 == 0) {
	    goto L120;
	}
	i__1 = j1;
	for (k = 1; k <= i__1; ++k) {
	    --j;
	    if (a[j + a_dim1] <= *tol) {
		goto L110;
	    }
	    ++kk;
	    --jj;
	    aa[jj + kk * aa_dim1] = a[j + (k + 1) * a_dim1];
L110:
	    ;
	}
L120:
	;
    }
/*  form successively in h the columns of a with a zero diagonal element. */
    ii = 0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	++ii;
	if (a[i__ + a_dim1] > *tol) {
	    goto L200;
	}
	--ii;
	if (ii == 0) {
	    goto L200;
	}
	jj = 1;
	j = i__;
/* Computing MIN */
	i__1 = j - 1;
	j1 = min(i__1,m1);
	i__1 = j1;
	for (k = 1; k <= i__1; ++k) {
	    --j;
	    if (a[j + a_dim1] <= *tol) {
		goto L130;
	    }
	    h__[jj] = a[j + (k + 1) * a_dim1];
	    ++jj;
L130:
	    ;
	}
	i__1 = *m;
	for (kk = jj; kk <= i__1; ++kk) {
	    h__[kk] = 0.f;
/* L140: */
	}
/*  rotate this column into aa by givens transformations. */
	jj = ii;
	i__1 = ii;
	for (i1 = 1; i1 <= i__1; ++i1) {
/* Computing MIN */
	    i__3 = jj - 1;
	    j1 = min(i__3,m1);
	    piv = h__[1];
	    if (piv != 0.f) {
		goto L160;
	    }
	    if (j1 == 0) {
		goto L200;
	    }
	    i__3 = j1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
		j3 = j2 + 1;
		h__[j2] = h__[j3];
/* L150: */
	    }
	    goto L180;
L160:
	    fpgivs_(&piv, &aa[jj + aa_dim1], &cos__, &sin__);
	    if (j1 == 0) {
		goto L200;
	    }
	    kk = jj;
	    i__3 = j1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
		j3 = j2 + 1;
		--kk;
		fprota_(&cos__, &sin__, &h__[j3], &aa[kk + j3 * aa_dim1]);
		h__[j2] = h__[j3];
/* L170: */
	    }
L180:
	    --jj;
	    h__[j3] = 0.f;
/* L190: */
	}
L200:
	;
    }
/*  solve the system (aa) (f1) = ff */
    ff[*rank] /= aa[*rank + aa_dim1];
    i__ = *rank - 1;
    if (i__ == 0) {
	goto L230;
    }
    i__2 = *rank;
    for (j = 2; j <= i__2; ++j) {
	store = ff[i__];
/* Computing MIN */
	i__1 = j - 1;
	i1 = min(i__1,m1);
	k = i__;
	i__1 = i1;
	for (ii = 1; ii <= i__1; ++ii) {
	    ++k;
	    stor1 = ff[k];
	    stor2 = aa[i__ + (ii + 1) * aa_dim1];
	    store -= stor1 * stor2;
/* L210: */
	}
	stor1 = aa[i__ + aa_dim1];
	ff[i__] = store / stor1;
	--i__;
/* L220: */
    }
/*  solve the system  (aa)' (f2) = f1 */
L230:
    ff[1] /= aa[aa_dim1 + 1];
    if (*rank == 1) {
	goto L260;
    }
    i__2 = *rank;
    for (j = 2; j <= i__2; ++j) {
	store = ff[j];
/* Computing MIN */
	i__1 = j - 1;
	i1 = min(i__1,m1);
	k = j;
	i__1 = i1;
	for (ii = 1; ii <= i__1; ++ii) {
	    --k;
	    stor1 = ff[k];
	    stor2 = aa[k + (ii + 1) * aa_dim1];
	    store -= stor1 * stor2;
/* L240: */
	}
	stor1 = aa[j + aa_dim1];
	ff[j] = store / stor1;
/* L250: */
    }
/*  premultiply f2 by the transpoze of a. */
L260:
    k = 0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	store = 0.f;
	if (a[i__ + a_dim1] > *tol) {
	    ++k;
	}
	j1 = min(i__,*m);
	kk = k;
	ij = i__ + 1;
	i__1 = j1;
	for (j = 1; j <= i__1; ++j) {
	    --ij;
	    if (a[ij + a_dim1] <= *tol) {
		goto L270;
	    }
	    stor1 = a[ij + j * a_dim1];
	    stor2 = ff[kk];
	    store += stor1 * stor2;
	    --kk;
L270:
	    ;
	}
	c__[i__] = store;
/* L280: */
    }
/*  add to the sum of squared residuals the contribution of putting */
/*  to zero the small diagonal elements of matrix (a). */
    stor3 = 0.f;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (a[i__ + a_dim1] > *tol) {
	    goto L310;
	}
	store = f[i__];
/* Computing MIN */
	i__1 = *n - i__;
	i1 = min(i__1,m1);
	if (i1 == 0) {
	    goto L300;
	}
	i__1 = i1;
	for (j = 1; j <= i__1; ++j) {
	    ij = i__ + j;
	    stor1 = c__[ij];
	    stor2 = a[i__ + (j + 1) * a_dim1];
	    store -= stor1 * stor2;
/* L290: */
	}
L300:
	fac = a[i__ + a_dim1] * c__[i__];
	stor1 = a[i__ + a_dim1];
	stor2 = c__[i__];
	stor1 *= stor2;
	stor3 += stor1 * (stor1 - store - store);
L310:
	;
    }
    fac = stor3;
    *sq += fac;
    return 0;
} /* fprank_ */

