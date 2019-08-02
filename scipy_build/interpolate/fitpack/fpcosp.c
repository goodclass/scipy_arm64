/* fpcosp.f -- translated by f2c (version 20190311).
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

static integer c__3 = 3;

/* Subroutine */ int fpcosp_(integer *m, doublereal *x, doublereal *y, 
	doublereal *w, integer *n, doublereal *t, doublereal *e, integer *
	maxtr, integer *maxbin, doublereal *c__, doublereal *sq, doublereal *
	sx, logical *bind, integer *nm, integer *mb, doublereal *a, 
	doublereal *b, doublereal *const__, doublereal *z__, doublereal *zz, 
	doublereal *u, doublereal *q, integer *info, integer *up, integer *
	left, integer *right, integer *jbind, integer *ibind, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal f, h__[4];
    static integer i__, j, k, l, i1, j1, j2, j3, k1, k2, k3, k4, k5, k6, l1, 
	    l2, l3, n1, n4, n6;
    static doublereal wi, xi;
    static integer lp1, kdim, merk, nbind, count;
    extern /* Subroutine */ int fpadno_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), fpdeno_(integer *, integer *, integer *, integer *, 
	    integer *, integer *), fpbspl_(doublereal *, integer *, integer *,
	     doublereal *, integer *, doublereal *);
    static integer number;
    extern /* Subroutine */ int fpseno_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local array.. */
/*  ..subroutine references.. */
/*    fpbspl,fpadno,fpdeno,fpfrno,fpseno */
/*  .. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  if we use the b-spline representation of s(x) our approximation     c */
/*  problem results in a quadratic programming problem:                 c */
/*    find the b-spline coefficients c(j),j=1,2,...n-4 such that        c */
/*        (1) sumi((wi*(yi-sumj(cj*nj(xi))))**2),i=1,2,...m is minimal  c */
/*        (2) sumj(cj*n''j(t(l+3)))*e(l) <= 0, l=1,2,...n-6.            c */
/*  to solve this problem we use the theil-van de panne procedure.      c */
/*  if the inequality constraints (2) are numbered from 1 to n-6,       c */
/*  this algorithm finds a subset of constraints ibind(1)..ibind(nbind) c */
/*  such that the solution of the minimization problem (1) with these   c */
/*  constraints in equality form, satisfies all constraints. such a     c */
/*  feasible solution is optimal if the lagrange parameters associated  c */
/*  with that problem with equality constraints, are all positive.      c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  determine n6, the number of inequality constraints. */
    /* Parameter adjustments */
    q_dim1 = *m;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --sx;
    --w;
    --y;
    --x;
    --zz;
    --z__;
    --const__;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --bind;
    --c__;
    --e;
    --t;
    --right;
    --left;
    --up;
    --info;
    --u;
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --ibind;
    --jbind;

    /* Function Body */
    n6 = *n - 6;
/*  fix the parameters which determine these constraints. */
    i__1 = n6;
    for (i__ = 1; i__ <= i__1; ++i__) {
	const__[i__] = e[i__] * (t[i__ + 4] - t[i__ + 1]) / (t[i__ + 5] - t[
		i__ + 2]);
/* L10: */
    }
/*  initialize the triply linked tree which is used to find the subset */
/*  of constraints ibind(1),...ibind(nbind). */
    count = 1;
    info[1] = 0;
    left[1] = 0;
    right[1] = 0;
    up[1] = 1;
    merk = 1;
/*  set up the normal equations  n'nc=n'y  where n denotes the m x (n-4) */
/*  observation matrix with elements ni,j = wi*nj(xi)  and y is the */
/*  column vector with elements yi*wi. */
/*  from the properties of the b-splines nj(x),j=1,2,...n-4, it follows */
/*  that  n'n  is a (n-4) x (n-4)  positive definite bandmatrix of */
/*  bandwidth 7. the matrices n'n and n'y are built up in a and z. */
    n4 = *n - 4;
/*  initialization */
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    a[i__ + j * a_dim1] = 0.f;
/* L20: */
	}
    }
    l = 4;
    lp1 = l + 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch the current row of the observation matrix. */
	xi = x[i__];
/* Computing 2nd power */
	d__1 = w[i__];
	wi = d__1 * d__1;
/*  search for knot interval  t(l) <= xi < t(l+1) */
L30:
	if (xi < t[lp1] || l == n4) {
	    goto L40;
	}
	l = lp1;
	lp1 = l + 1;
	goto L30;
/*  evaluate the four non-zero cubic b-splines nj(xi),j=l-3,...l. */
L40:
	fpbspl_(&t[1], n, &c__3, &xi, &l, h__);
/*  store in q these values h(1),h(2),...h(4). */
	for (j = 1; j <= 4; ++j) {
	    q[i__ + j * q_dim1] = h__[j - 1];
/* L50: */
	}
/*  add the contribution of the current row of the observation matrix */
/*  n to the normal equations. */
	l3 = l - 3;
	k1 = 0;
	i__2 = l;
	for (j1 = l3; j1 <= i__2; ++j1) {
	    ++k1;
	    f = h__[k1 - 1];
	    z__[j1] += f * wi * y[i__];
	    k2 = k1;
	    j2 = 4;
	    i__3 = l;
	    for (j3 = j1; j3 <= i__3; ++j3) {
		a[j3 + j2 * a_dim1] += f * wi * h__[k2 - 1];
		++k2;
		--j2;
/* L60: */
	    }
	}
/* L70: */
    }
/*  since n'n is a symmetric matrix it can be factorized as */
/*        (3)  n'n = (r1)'(d1)(r1) */
/*  with d1 a diagonal matrix and r1 an (n-4) x (n-4)  unit upper */
/*  triangular matrix of bandwidth 4. the matrices r1 and d1 are built */
/*  up in a. at the same time we solve the systems of equations */
/*        (4)  (r1)'(z2) = n'y */
/*        (5)  (d1) (z1) = (z2) */
/*  the vectors z2 and z1 are kept in zz and z. */
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k1 = 1;
	if (i__ < 4) {
	    k1 = 5 - i__;
	}
	k2 = i__ - 4 + k1;
	k3 = k2;
	for (j = k1; j <= 4; ++j) {
	    k4 = j - 1;
	    k5 = 4 - j + k1;
	    f = a[i__ + j * a_dim1];
	    if (k1 > k4) {
		goto L90;
	    }
	    k6 = k2;
	    i__3 = k4;
	    for (k = k1; k <= i__3; ++k) {
		f -= a[i__ + k * a_dim1] * a[k3 + k5 * a_dim1] * a[k6 + (
			a_dim1 << 2)];
		++k5;
		++k6;
/* L80: */
	    }
L90:
	    if (j == 4) {
		goto L110;
	    }
	    a[i__ + j * a_dim1] = f / a[k3 + (a_dim1 << 2)];
	    ++k3;
/* L100: */
	}
L110:
	a[i__ + (a_dim1 << 2)] = f;
	f = z__[i__];
	if (i__ == 1) {
	    goto L130;
	}
	k4 = i__;
	for (j = k1; j <= 3; ++j) {
	    k = k1 + 3 - j;
	    --k4;
	    f -= a[i__ + k * a_dim1] * z__[k4] * a[k4 + (a_dim1 << 2)];
/* L120: */
	}
L130:
	z__[i__] = f / a[i__ + (a_dim1 << 2)];
	zz[i__] = f;
/* L140: */
    }
/*  start computing the least-squares cubic spline without taking account */
/*  of any constraint. */
    nbind = 0;
    n1 = 1;
    ibind[1] = 0;
/*  main loop for the least-squares problems with different subsets of */
/*  the constraints (2) in equality form. the resulting b-spline coeff. */
/*  c and lagrange parameters u are the solution of the system */
/*            ! n'n  b' ! ! c !   ! n'y ! */
/*        (6) !         ! !   ! = !     ! */
/*            !  b   0  ! ! u !   !  0  ! */
/*  z1 is stored into array c. */
L150:
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = z__[i__];
/* L160: */
    }
/*  if there are no equality constraints, compute the coeff. c directly. */
    if (nbind == 0) {
	goto L370;
    }
/*  initialization */
    kdim = n4 + nbind;
    i__1 = nbind;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = kdim;
	for (j = 1; j <= i__3; ++j) {
	    b[j + i__ * b_dim1] = 0.f;
/* L170: */
	}
    }
/*  matrix b is built up,expressing that the constraints nrs ibind(1),... */
/*  ibind(nbind) must be satisfied in equality form. */
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
	l = ibind[i__];
	b[l + i__ * b_dim1] = e[l];
	b[l + 1 + i__ * b_dim1] = -(e[l] + const__[l]);
	b[l + 2 + i__ * b_dim1] = const__[l];
/* L180: */
    }
/*  find the matrix (b1) as the solution of the system of equations */
/*        (7)  (r1)'(d1)(b1) = b' */
/*  (b1) is built up in the upper part of the array b(rows 1,...n-4). */
    i__3 = nbind;
    for (k1 = 1; k1 <= i__3; ++k1) {
	l = ibind[k1];
	i__1 = n4;
	for (i__ = l; i__ <= i__1; ++i__) {
	    f = b[i__ + k1 * b_dim1];
	    if (i__ == 1) {
		goto L200;
	    }
	    k2 = 3;
	    if (i__ < 4) {
		k2 = i__ - 1;
	    }
	    i__2 = k2;
	    for (k3 = 1; k3 <= i__2; ++k3) {
		l1 = i__ - k3;
		l2 = 4 - k3;
		f -= b[l1 + k1 * b_dim1] * a[i__ + l2 * a_dim1] * a[l1 + (
			a_dim1 << 2)];
/* L190: */
	    }
L200:
	    b[i__ + k1 * b_dim1] = f / a[i__ + (a_dim1 << 2)];
/* L210: */
	}
/* L220: */
    }
/*  factorization of the symmetric matrix  -(b1)'(d1)(b1) */
/*        (8)  -(b1)'(d1)(b1) = (r2)'(d2)(r2) */
/*  with (d2) a diagonal matrix and (r2) an nbind x nbind unit upper */
/*  triangular matrix. the matrices r2 and d2 are built up in the lower */
/*  part of the array b (rows n-3,n-2,...n-4+nbind). */
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i1 = i__ - 1;
	i__1 = nbind;
	for (j = i__; j <= i__1; ++j) {
	    f = 0.f;
	    i__2 = n4;
	    for (k = 1; k <= i__2; ++k) {
		f += b[k + i__ * b_dim1] * b[k + j * b_dim1] * a[k + (a_dim1 
			<< 2)];
/* L230: */
	    }
	    k1 = n4 + 1;
	    if (i1 == 0) {
		goto L250;
	    }
	    i__2 = i1;
	    for (k = 1; k <= i__2; ++k) {
		f += b[k1 + i__ * b_dim1] * b[k1 + j * b_dim1] * b[k1 + k * 
			b_dim1];
		++k1;
/* L240: */
	    }
L250:
	    b[k1 + j * b_dim1] = -f;
	    if (j == i__) {
		goto L260;
	    }
	    b[k1 + j * b_dim1] /= b[k1 + i__ * b_dim1];
L260:
	    ;
	}
/* L270: */
    }
/*  according to (3),(7) and (8) the system of equations (6) becomes */
/*         ! (r1)'    0  ! ! (d1)    0  ! ! (r1)  (b1) ! ! c !   ! n'y ! */
/*    (9)  !             ! !            ! !            ! !   ! = !     ! */
/*         ! (b1)'  (r2)'! !   0   (d2) ! !   0   (r2) ! ! u !   !  0  ! */
/*  backward substitution to obtain the b-spline coefficients c(j),j=1,.. */
/*  n-4 and the lagrange parameters u(j),j=1,2,...nbind. */
/*  first step of the backward substitution: solve the system */
/*             ! (r1)'(d1)      0     ! ! (c1) !   ! n'y ! */
/*        (10) !                      ! !      ! = !     ! */
/*             ! (b1)'(d1)  (r2)'(d2) ! ! (u1) !   !  0  ! */
/*  from (4) and (5) we know that this is equivalent to */
/*        (11)  (c1) = (z1) */
/*        (12)  (r2)'(d2)(u1) = -(b1)'(z2) */
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
	f = 0.f;
	i__1 = n4;
	for (j = 1; j <= i__1; ++j) {
	    f += b[j + i__ * b_dim1] * zz[j];
/* L280: */
	}
	i1 = i__ - 1;
	k1 = n4 + 1;
	if (i1 == 0) {
	    goto L300;
	}
	i__1 = i1;
	for (j = 1; j <= i__1; ++j) {
	    f += u[j] * b[k1 + i__ * b_dim1] * b[k1 + j * b_dim1];
	    ++k1;
/* L290: */
	}
L300:
	u[i__] = -f / b[k1 + i__ * b_dim1];
/* L310: */
    }
/*  second step of the backward substitution: solve the system */
/*             ! (r1)  (b1) ! ! c !   ! c1 ! */
/*        (13) !            ! !   ! = !    ! */
/*             !   0   (r2) ! ! u !   ! u1 ! */
    k1 = nbind;
    k2 = kdim;
/*  find the lagrange parameters u. */
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
	f = u[k1];
	if (i__ == 1) {
	    goto L330;
	}
	k3 = k1 + 1;
	i__1 = nbind;
	for (j = k3; j <= i__1; ++j) {
	    f -= u[j] * b[k2 + j * b_dim1];
/* L320: */
	}
L330:
	u[k1] = f;
	--k1;
	--k2;
/* L340: */
    }
/*  find the b-spline coefficients c. */
    i__3 = n4;
    for (i__ = 1; i__ <= i__3; ++i__) {
	f = c__[i__];
	i__1 = nbind;
	for (j = 1; j <= i__1; ++j) {
	    f -= u[j] * b[i__ + j * b_dim1];
/* L350: */
	}
	c__[i__] = f;
/* L360: */
    }
L370:
    k1 = n4;
    i__3 = n4;
    for (i__ = 2; i__ <= i__3; ++i__) {
	--k1;
	f = c__[k1];
	k2 = 1;
	if (i__ < 5) {
	    k2 = 5 - i__;
	}
	k3 = k1;
	l = 3;
	for (j = k2; j <= 3; ++j) {
	    ++k3;
	    f -= a[k3 + l * a_dim1] * c__[k3];
	    --l;
/* L380: */
	}
	c__[k1] = f;
/* L390: */
    }
/*  test whether the solution of the least-squares problem with the */
/*  constraints ibind(1),...ibind(nbind) in equality form, satisfies */
/*  all of the constraints (2). */
    k = 1;
/*  number counts the number of violated inequality constraints. */
    number = 0;
    i__3 = n6;
    for (j = 1; j <= i__3; ++j) {
	l = ibind[k];
	++k;
	if (j == l) {
	    goto L440;
	}
	--k;
/*  test whether constraint j is satisfied */
	f = e[j] * (c__[j] - c__[j + 1]) + const__[j] * (c__[j + 2] - c__[j + 
		1]);
	if (f <= 0.f) {
	    goto L440;
	}
/*  if constraint j is not satisfied, add a branch of length nbind+1 */
/*  to the tree. the nodes of this branch contain in their information */
/*  field the number of the constraints ibind(1),...ibind(nbind) and j, */
/*  arranged in increasing order. */
	++number;
	k1 = k - 1;
	if (k1 == 0) {
	    goto L410;
	}
	i__1 = k1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jbind[i__] = ibind[i__];
/* L400: */
	}
L410:
	jbind[k] = j;
	if (l == 0) {
	    goto L430;
	}
	i__1 = nbind;
	for (i__ = k; i__ <= i__1; ++i__) {
	    jbind[i__ + 1] = ibind[i__];
/* L420: */
	}
L430:
	fpadno_(maxtr, &up[1], &left[1], &right[1], &info[1], &count, &merk, &
		jbind[1], &n1, ier);
/*  test whether the storage space which is required for the tree,exceeds */
/*  the available storage space. */
	if (*ier != 0) {
	    goto L560;
	}
L440:
	;
    }
/*  test whether the solution of the least-squares problem with equality */
/*  constraints is a feasible solution. */
    if (number == 0) {
	goto L470;
    }
/*  test whether there are still cases with nbind constraints in */
/*  equality form to be considered. */
L450:
    if (merk > 1) {
	goto L460;
    }
    nbind = n1;
/*  test whether the number of knots where s''(x)=0 exceeds maxbin. */
    if (nbind > *maxbin) {
	goto L550;
    }
    ++n1;
    ibind[n1] = 0;
/*  search which cases with nbind constraints in equality form */
/*  are going to be considered. */
    fpdeno_(maxtr, &up[1], &left[1], &right[1], &nbind, &merk);
/*  test whether the quadratic programming problem has a solution. */
    if (merk == 1) {
	goto L570;
    }
/*  find a new case with nbind constraints in equality form. */
L460:
    fpseno_(maxtr, &up[1], &left[1], &right[1], &info[1], &merk, &ibind[1], &
	    nbind);
    goto L150;
/*  test whether the feasible solution is optimal. */
L470:
    *ier = 0;
    i__3 = n6;
    for (i__ = 1; i__ <= i__3; ++i__) {
	bind[i__] = FALSE_;
/* L480: */
    }
    if (nbind == 0) {
	goto L500;
    }
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
	if (u[i__] <= 0.f) {
	    goto L450;
	}
	j = ibind[i__];
	bind[j] = TRUE_;
/* L490: */
    }
/*  evaluate s(x) at the data points x(i) and calculate the weighted */
/*  sum of squared residual right hand sides sq. */
L500:
    *sq = 0.f;
    l = 4;
    lp1 = 5;
    i__3 = *m;
    for (i__ = 1; i__ <= i__3; ++i__) {
L510:
	if (x[i__] < t[lp1] || l == n4) {
	    goto L520;
	}
	l = lp1;
	lp1 = l + 1;
	goto L510;
L520:
	sx[i__] = c__[l - 3] * q[i__ + q_dim1] + c__[l - 2] * q[i__ + (q_dim1 
		<< 1)] + c__[l - 1] * q[i__ + q_dim1 * 3] + c__[l] * q[i__ + (
		q_dim1 << 2)];
/* Computing 2nd power */
	d__1 = w[i__] * (y[i__] - sx[i__]);
	*sq += d__1 * d__1;
/* L530: */
    }
    goto L600;
/*  error codes and messages. */
L550:
    *ier = 1;
    goto L600;
L560:
    *ier = 2;
    goto L600;
L570:
    *ier = 3;
L600:
    return 0;
} /* fpcosp_ */

