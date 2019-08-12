/* fpclos.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpclos_(integer *iopt, integer *idim, integer *m, 
	doublereal *u, integer *mx, doublereal *x, doublereal *w, integer *k, 
	doublereal *s, integer *nest, doublereal *tol, integer *maxit, 
	integer *k1, integer *k2, integer *n, doublereal *t, integer *nc, 
	doublereal *c__, doublereal *fp, doublereal *fpint, doublereal *z__, 
	doublereal *a1, doublereal *a2, doublereal *b, doublereal *g1, 
	doublereal *g2, doublereal *q, integer *nrdata, integer *ier)
{
    /* System generated locals */
    integer a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, b_offset, g1_dim1,
	     g1_offset, g2_dim1, g2_offset, q_dim1, q_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l;
    static doublereal p, d1, f1, f2, f3;
    static integer i1, i2, i3;
    static doublereal p1, p2, p3;
    static integer j1, j2, k3, l0, l1, l5, m1, n7, n8;
    static doublereal h1[7], h2[6];
    static integer n10, n11, ij, ik, jj, jk, kk, mm, it;
    static doublereal ui, wi, rn, xi[10], fp0;
    static integer kk1, nk1, nk2;
    static doublereal acc, fac, one, cos__, per, sin__;
    static integer new__;
    static doublereal piv;
    static integer ich1, ich3;
    static doublereal con1, con4, con9;
    static integer npl1;
    static doublereal half;
    static integer jper, nmin, iter, nmax;
    static doublereal fpms, term, pinv, fpold, fpart;
    static integer nrint;
    static doublereal store;
    static integer nplus;
    extern /* Subroutine */ int fpbacp_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), fpdisc_(doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), fprota_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpgivs_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpknot_(doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota */
/*  .. */
/*  set constants */
    /* Parameter adjustments */
    --w;
    --u;
    --x;
    --nrdata;
    a2_dim1 = *nest;
    a2_offset = 1 + a2_dim1;
    a2 -= a2_offset;
    --fpint;
    --t;
    q_dim1 = *m;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    g2_dim1 = *nest;
    g2_offset = 1 + g2_dim1;
    g2 -= g2_offset;
    a1_dim1 = *nest;
    a1_offset = 1 + a1_dim1;
    a1 -= a1_offset;
    g1_dim1 = *nest;
    g1_offset = 1 + g1_dim1;
    g1 -= g1_offset;
    b_dim1 = *nest;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --z__;
    --c__;

    /* Function Body */
    one = 1.f;
    con1 = .1f;
    con9 = .9f;
    con4 = .04f;
    half = .5f;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 1: determination of the number of knots and their position     c */
/*  **************************************************************      c */
/*  given a set of knots we compute the least-squares closed curve      c */
/*  sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    c */
/*  if iopt=-1 sinf(u) is the requested curve                           c */
/*  if iopt=0 or iopt=1 we check whether we can accept the knots:       c */
/*    if fp <=s we will continue with the current set of knots.         c */
/*    if fp > s we will increase the number of knots and compute the    c */
/*       corresponding least-squares curve until finally fp<=s.         c */
/*  the initial choice of knots depends on the value of s and iopt.     c */
/*    if s=0 we have spline interpolation; in that case the number of   c */
/*    knots equals nmax = m+2*k.                                        c */
/*    if s > 0 and                                                      c */
/*      iopt=0 we first compute the least-squares polynomial curve of   c */
/*      degree k; n = nmin = 2*k+2. since s(u) must be periodic we      c */
/*      find that s(u) reduces to a fixed point.                        c */
/*      iopt=1 we start with the set of knots found at the last         c */
/*      call of the routine, except for the case that s > fp0; then     c */
/*      we compute directly the least-squares polynomial curve.         c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    m1 = *m - 1;
    kk = *k;
    kk1 = *k1;
    k3 = *k * 3 + 1;
    nmin = *k1 << 1;
/*  determine the length of the period of the splines. */
    per = u[*m] - u[1];
    if (*iopt < 0) {
	goto L50;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = *tol * *s;
/*  determine nmax, the number of knots for periodic spline interpolation */
    nmax = *m + (*k << 1);
    if (*s > 0.f || nmax == nmin) {
	goto L30;
    }
/*  if s=0, s(u) is an interpolating curve. */
    *n = nmax;
/*  test whether the required storage space exceeds the available one. */
    if (*n > *nest) {
	goto L620;
    }
/*  find the position of the interior knots in case of interpolation. */
L5:
    if (*k / 2 << 1 == *k) {
	goto L20;
    }
    i__1 = m1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = i__ + *k;
	t[j] = u[i__];
/* L10: */
    }
    if (*s > 0.f) {
	goto L50;
    }
    kk = *k - 1;
    kk1 = *k;
    if (kk > 0) {
	goto L50;
    }
    t[1] = t[*m] - per;
    t[2] = u[1];
    t[*m + 1] = u[*m];
    t[*m + 2] = t[3] + per;
    jj = 0;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__;
	i__2 = *idim;
	for (j1 = 1; j1 <= i__2; ++j1) {
	    ++jj;
	    c__[j] = x[jj];
	    j += *n;
/* L12: */
	}
/* L15: */
    }
    jj = 1;
    j = *m;
    i__1 = *idim;
    for (j1 = 1; j1 <= i__1; ++j1) {
	c__[j] = c__[jj];
	j += *n;
	jj += *n;
/* L17: */
    }
    *fp = 0.f;
    fpint[*n] = fp0;
    fpint[*n - 1] = 0.f;
    nrdata[*n] = 0;
    goto L630;
L20:
    i__1 = m1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = i__ + *k;
	t[j] = (u[i__] + u[i__ - 1]) * half;
/* L25: */
    }
    goto L50;
/*  if s > 0 our initial choice depends on the value of iopt. */
/*  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares */
/*  polynomial curve. (i.e. a constant point). */
/*  if iopt=1 and fp0>s we start computing the least-squares closed */
/*  curve according the set of knots found at the last call of the */
/*  routine. */
L30:
    if (*iopt == 0) {
	goto L35;
    }
    if (*n == nmin) {
	goto L35;
    }
    fp0 = fpint[*n];
    fpold = fpint[*n - 1];
    nplus = nrdata[*n];
    if (fp0 > *s) {
	goto L50;
    }
/*  the case that s(u) is a fixed point is treated separetely. */
/*  fp0 denotes the corresponding sum of squared residuals. */
L35:
    fp0 = 0.f;
    d1 = 0.f;
    i__1 = *idim;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L37: */
    }
    jj = 0;
    i__1 = m1;
    for (it = 1; it <= i__1; ++it) {
	wi = w[it];
	fpgivs_(&wi, &d1, &cos__, &sin__);
	i__2 = *idim;
	for (j = 1; j <= i__2; ++j) {
	    ++jj;
	    fac = wi * x[jj];
	    fprota_(&cos__, &sin__, &fac, &z__[j]);
/* Computing 2nd power */
	    d__1 = fac;
	    fp0 += d__1 * d__1;
/* L40: */
	}
/* L45: */
    }
    i__1 = *idim;
    for (j = 1; j <= i__1; ++j) {
	z__[j] /= d1;
/* L47: */
    }
/*  test whether that fixed point is a solution of our problem. */
    fpms = fp0 - *s;
    if (fpms < acc || nmax == nmin) {
	goto L640;
    }
    fpold = fp0;
/*  test whether the required storage space exceeds the available one. */
    if (*n >= *nest) {
	goto L620;
    }
/*  start computing the least-squares closed curve with one */
/*  interior knot. */
    nplus = 1;
    *n = nmin + 1;
    mm = (*m + 1) / 2;
    t[*k2] = u[mm];
    nrdata[1] = mm - 2;
    nrdata[2] = m1 - mm;
/*  main loop for the different sets of knots. m is a save upper */
/*  bound for the number of trials. */
L50:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find nrint, the number of knot intervals. */
	nrint = *n - nmin + 1;
/*  find the position of the additional knots which are needed for */
/*  the b-spline representation of s(u). if we take */
/*      t(k+1) = u(1), t(n-k) = u(m) */
/*      t(k+1-j) = t(n-k-j) - per, j=1,2,...k */
/*      t(n-k+j) = t(k+1+j) + per, j=1,2,...k */
/*  then s(u) will be a smooth closed curve if the b-spline */
/*  coefficients satisfy the following conditions */
/*      c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**) */
/*  with n7=n-2*k-1. */
	t[*k1] = u[1];
	nk1 = *n - *k1;
	nk2 = nk1 + 1;
	t[nk2] = u[*m];
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
	    i1 = nk2 + j;
	    i2 = nk2 - j;
	    j1 = *k1 + j;
	    j2 = *k1 - j;
	    t[i1] = t[j1] + per;
	    t[j2] = t[i2] - per;
/* L60: */
	}
/*  compute the b-spline coefficients of the least-squares closed curve */
/*  sinf(u). the observation matrix a is built up row by row while */
/*  taking into account condition (**) and is reduced to triangular */
/*  form by givens transformations . */
/*  at the same time fp=f(p=inf) is computed. */
/*  the n7 x n7 triangularised upper matrix a has the form */
/*            ! a1 '    ! */
/*        a = !    ' a2 ! */
/*            ! 0  '    ! */
/*  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular */
/*  matrix of bandwidth k+1 ( n10 = n7-k). */
/*  initialization. */
	i__2 = *nc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] = 0.f;
/* L65: */
	}
	i__2 = nk1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = kk1;
	    for (j = 1; j <= i__3; ++j) {
		a1[i__ + j * a1_dim1] = 0.f;
/* L70: */
	    }
	}
	n7 = nk1 - *k;
	n10 = n7 - kk;
	jper = 0;
	*fp = 0.f;
	l = *k1;
	jj = 0;
	i__3 = m1;
	for (it = 1; it <= i__3; ++it) {
/*  fetch the current data point u(it),x(it) */
	    ui = u[it];
	    wi = w[it];
	    i__2 = *idim;
	    for (j = 1; j <= i__2; ++j) {
		++jj;
		xi[j - 1] = x[jj] * wi;
/* L75: */
	    }
/*  search for knot interval t(l) <= ui < t(l+1). */
L80:
	    if (ui < t[l + 1]) {
		goto L85;
	    }
	    ++l;
	    goto L80;
/*  evaluate the (k+1) non-zero b-splines at ui and store them in q. */
L85:
	    fpbspl_(&t[1], n, k, &ui, &l, h__);
	    i__2 = *k1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		q[it + i__ * q_dim1] = h__[i__ - 1];
		h__[i__ - 1] *= wi;
/* L90: */
	    }
	    l5 = l - *k1;
/*  test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui */
	    if (l5 < n10) {
		goto L285;
	    }
	    if (jper != 0) {
		goto L160;
	    }
/*  initialize the matrix a2. */
	    i__2 = n7;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__4 = kk;
		for (j = 1; j <= i__4; ++j) {
		    a2[i__ + j * a2_dim1] = 0.f;
/* L95: */
		}
	    }
	    jk = n10 + 1;
	    i__4 = kk;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		ik = jk;
		i__2 = kk1;
		for (j = 1; j <= i__2; ++j) {
		    if (ik <= 0) {
			goto L105;
		    }
		    a2[ik + i__ * a2_dim1] = a1[ik + j * a1_dim1];
		    --ik;
/* L100: */
		}
L105:
		++jk;
/* L110: */
	    }
	    jper = 1;
/*  if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui */
/*  we take account of condition (**) for setting up the new row */
/*  of the observation matrix a. this row is stored in the arrays h1 */
/*  (the part with respect to a1) and h2 (the part with */
/*  respect to a2). */
L160:
	    i__4 = kk;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		h1[i__ - 1] = 0.f;
		h2[i__ - 1] = 0.f;
/* L170: */
	    }
	    h1[kk1 - 1] = 0.f;
	    j = l5 - n10;
	    i__4 = kk1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		++j;
		l0 = j;
L180:
		l1 = l0 - kk;
		if (l1 <= 0) {
		    goto L200;
		}
		if (l1 <= n10) {
		    goto L190;
		}
		l0 = l1 - n10;
		goto L180;
L190:
		h1[l1 - 1] = h__[i__ - 1];
		goto L210;
L200:
		h2[l0 - 1] += h__[i__ - 1];
L210:
		;
	    }
/*  rotate the new row of the observation matrix into triangle */
/*  by givens transformations. */
	    if (n10 <= 0) {
		goto L250;
	    }
/*  rotation with the rows 1,2,...n10 of matrix a. */
	    i__4 = n10;
	    for (j = 1; j <= i__4; ++j) {
		piv = h1[0];
		if (piv != 0.f) {
		    goto L214;
		}
		i__2 = kk;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    h1[i__ - 1] = h1[i__];
/* L212: */
		}
		h1[kk1 - 1] = 0.f;
		goto L240;
/*  calculate the parameters of the givens transformation. */
L214:
		fpgivs_(&piv, &a1[j + a1_dim1], &cos__, &sin__);
/*  transformation to the right hand side. */
		j1 = j;
		i__2 = *idim;
		for (j2 = 1; j2 <= i__2; ++j2) {
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &z__[j1]);
		    j1 += *n;
/* L217: */
		}
/*  transformations to the left hand side with respect to a2. */
		i__2 = kk;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &a2[j + i__ * 
			    a2_dim1]);
/* L220: */
		}
		if (j == n10) {
		    goto L250;
		}
/* Computing MIN */
		i__2 = n10 - j;
		i2 = min(i__2,kk);
/*  transformations to the left hand side with respect to a1. */
		i__2 = i2;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i1 = i__ + 1;
		    fprota_(&cos__, &sin__, &h1[i1 - 1], &a1[j + i1 * a1_dim1]
			    );
		    h1[i__ - 1] = h1[i1 - 1];
/* L230: */
		}
		h1[i1 - 1] = 0.f;
L240:
		;
	    }
/*  rotation with the rows n10+1,...n7 of matrix a. */
L250:
	    i__4 = kk;
	    for (j = 1; j <= i__4; ++j) {
		ij = n10 + j;
		if (ij <= 0) {
		    goto L270;
		}
		piv = h2[j - 1];
		if (piv == 0.f) {
		    goto L270;
		}
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &a2[ij + j * a2_dim1], &cos__, &sin__);
/*  transformations to right hand side. */
		j1 = ij;
		i__2 = *idim;
		for (j2 = 1; j2 <= i__2; ++j2) {
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &z__[j1]);
		    j1 += *n;
/* L255: */
		}
		if (j == kk) {
		    goto L280;
		}
		j1 = j + 1;
/*  transformations to left hand side. */
		i__2 = kk;
		for (i__ = j1; i__ <= i__2; ++i__) {
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &a2[ij + i__ * 
			    a2_dim1]);
/* L260: */
		}
L270:
		;
	    }
/*  add contribution of this row to the sum of squares of residual */
/*  right hand sides. */
L280:
	    i__4 = *idim;
	    for (j2 = 1; j2 <= i__4; ++j2) {
/* Computing 2nd power */
		d__1 = xi[j2 - 1];
		*fp += d__1 * d__1;
/* L282: */
	    }
	    goto L290;
/*  rotation of the new row of the observation matrix into */
/*  triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero */
/*  at ui. */
L285:
	    j = l5;
	    i__4 = kk1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		++j;
		piv = h__[i__ - 1];
		if (piv == 0.f) {
		    goto L140;
		}
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &a1[j + a1_dim1], &cos__, &sin__);
/*  transformations to right hand side. */
		j1 = j;
		i__2 = *idim;
		for (j2 = 1; j2 <= i__2; ++j2) {
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &z__[j1]);
		    j1 += *n;
/* L125: */
		}
		if (i__ == kk1) {
		    goto L150;
		}
		i2 = 1;
		i3 = i__ + 1;
/*  transformations to left hand side. */
		i__2 = kk1;
		for (i1 = i3; i1 <= i__2; ++i1) {
		    ++i2;
		    fprota_(&cos__, &sin__, &h__[i1 - 1], &a1[j + i2 * 
			    a1_dim1]);
/* L130: */
		}
L140:
		;
	    }
/*  add contribution of this row to the sum of squares of residual */
/*  right hand sides. */
L150:
	    i__4 = *idim;
	    for (j2 = 1; j2 <= i__4; ++j2) {
/* Computing 2nd power */
		d__1 = xi[j2 - 1];
		*fp += d__1 * d__1;
/* L155: */
	    }
L290:
	    ;
	}
	fpint[*n] = fp0;
	fpint[*n - 1] = fpold;
	nrdata[*n] = nplus;
/*  backward substitution to obtain the b-spline coefficients . */
	j1 = 1;
	i__3 = *idim;
	for (j2 = 1; j2 <= i__3; ++j2) {
	    fpbacp_(&a1[a1_offset], &a2[a2_offset], &z__[j1], &n7, &kk, &c__[
		    j1], &kk1, nest);
	    j1 += *n;
/* L292: */
	}
/*  calculate from condition (**) the remaining coefficients. */
	i__3 = *k;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    j1 = i__;
	    i__4 = *idim;
	    for (j = 1; j <= i__4; ++j) {
		j2 = j1 + n7;
		c__[j2] = c__[j1];
		j1 += *n;
/* L295: */
	    }
/* L297: */
	}
	if (*iopt < 0) {
	    goto L660;
	}
/*  test whether the approximation sinf(u) is an acceptable solution. */
	fpms = *fp - *s;
	if (abs(fpms) < acc) {
	    goto L660;
	}
/*  if f(p=inf) < s accept the choice of knots. */
	if (fpms < 0.f) {
	    goto L350;
	}
/*  if n=nmax, sinf(u) is an interpolating curve. */
	if (*n == nmax) {
	    goto L630;
	}
/*  increase the number of knots. */
/*  if n=nest we cannot increase the number of knots because of the */
/*  storage capacity limitation. */
	if (*n == *nest) {
	    goto L620;
	}
/*  determine the number of knots nplus we are going to add. */
	npl1 = nplus << 1;
	rn = (doublereal) nplus;
	if (fpold - *fp > acc) {
	    npl1 = (integer) (rn * fpms / (fpold - *fp));
	}
/* Computing MIN */
/* Computing MAX */
	i__2 = npl1, i__5 = nplus / 2, i__2 = max(i__2,i__5);
	i__3 = nplus << 1, i__4 = max(i__2,1);
	nplus = min(i__3,i__4);
	fpold = *fp;
/*  compute the sum of squared residuals for each knot interval */
/*  t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
	fpart = 0.f;
	i__ = 1;
	l = *k1;
	jj = 0;
	i__3 = m1;
	for (it = 1; it <= i__3; ++it) {
	    if (u[it] < t[l]) {
		goto L300;
	    }
	    new__ = 1;
	    ++l;
L300:
	    term = 0.f;
	    l0 = l - *k2;
	    i__4 = *idim;
	    for (j2 = 1; j2 <= i__4; ++j2) {
		fac = 0.f;
		j1 = l0;
		i__2 = *k1;
		for (j = 1; j <= i__2; ++j) {
		    ++j1;
		    fac += c__[j1] * q[it + j * q_dim1];
/* L305: */
		}
		++jj;
/* Computing 2nd power */
		d__1 = w[it] * (fac - x[jj]);
		term += d__1 * d__1;
		l0 += *n;
/* L310: */
	    }
	    fpart += term;
	    if (new__ == 0) {
		goto L320;
	    }
	    if (l > *k2) {
		goto L315;
	    }
	    fpint[nrint] = term;
	    new__ = 0;
	    goto L320;
L315:
	    store = term * half;
	    fpint[i__] = fpart - store;
	    ++i__;
	    fpart = store;
	    new__ = 0;
L320:
	    ;
	}
	fpint[nrint] += fpart;
	i__3 = nplus;
	for (l = 1; l <= i__3; ++l) {
/*  add a new knot */
	    fpknot_(&u[1], m, &t[1], n, &fpint[1], &nrdata[1], &nrint, nest, &
		    c__1);
/*  if n=nmax we locate the knots as for interpolation */
	    if (*n == nmax) {
		goto L5;
	    }
/*  test whether we cannot further increase the number of knots. */
	    if (*n == *nest) {
		goto L340;
	    }
/* L330: */
	}
/*  restart the computations with the new set of knots. */
L340:
	;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 2: determination of the smoothing closed curve sp(u).          c */
/*  **********************************************************          c */
/*  we have determined the number of knots and their position.          c */
/*  we now compute the b-spline coefficients of the smoothing curve     c */
/*  sp(u). the observation matrix a is extended by the rows of matrix   c */
/*  b expressing that the kth derivative discontinuities of sp(u) at    c */
/*  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c */
/*  ponding weights of these additional rows are set to 1/p.            c */
/*  iteratively we then have to determine the value of p such that f(p),c */
/*  the sum of squared residuals be = s. we already know that the least-c */
/*  squares polynomial curve corresponds to p=0, and that the least-    c */
/*  squares periodic spline curve corresponds to p=infinity. the        c */
/*  iteration process which is proposed here, makes use of rational     c */
/*  interpolation. since f(p) is a convex and strictly decreasing       c */
/*  function of p, it can be approximated by a rational function        c */
/*  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c */
/*  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c */
/*  to calculate the new value of p such that r(p)=s. convergence is    c */
/*  guaranteed by taking f1>0 and f3<0.                                 c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  evaluate the discontinuity jump of the kth derivative of the */
/*  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b. */
L350:
    fpdisc_(&t[1], n, k2, &b[b_offset], nest);
/*  initial value for p. */
    p1 = 0.f;
    f1 = fp0 - *s;
    p3 = -one;
    f3 = fpms;
    n11 = n10 - 1;
    n8 = n7 - 1;
    p = 0.f;
    l = n7;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *k + 1 - i__;
	p += a2[l + j * a2_dim1];
	--l;
	if (l == 0) {
	    goto L356;
	}
/* L352: */
    }
    i__1 = n10;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p += a1[i__ + a1_dim1];
/* L354: */
    }
L356:
    rn = (doublereal) n7;
    p = rn / p;
    ich1 = 0;
    ich3 = 0;
/*  iteration process to find the root of f(p) = s. */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*  form the matrix g  as the matrix a extended by the rows of matrix b. */
/*  the rows of matrix b with weight 1/p are rotated into */
/*  the triangularised observation matrix a. */
/*  after triangularisation our n7 x n7 matrix g takes the form */
/*            ! g1 '    ! */
/*        g = !    ' g2 ! */
/*            ! 0  '    ! */
/*  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular */
/*  matrix of bandwidth k+2. ( n11 = n7-k-1) */
	pinv = one / p;
/*  store matrix a into g */
	i__3 = *nc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    c__[i__] = z__[i__];
/* L358: */
	}
	i__3 = n7;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    g1[i__ + *k1 * g1_dim1] = a1[i__ + *k1 * a1_dim1];
	    g1[i__ + *k2 * g1_dim1] = 0.f;
	    g2[i__ + g2_dim1] = 0.f;
	    i__4 = *k;
	    for (j = 1; j <= i__4; ++j) {
		g1[i__ + j * g1_dim1] = a1[i__ + j * a1_dim1];
		g2[i__ + (j + 1) * g2_dim1] = a2[i__ + j * a2_dim1];
/* L360: */
	    }
	}
	l = n10;
	i__4 = *k1;
	for (j = 1; j <= i__4; ++j) {
	    if (l <= 0) {
		goto L375;
	    }
	    g2[l + g2_dim1] = a1[l + j * a1_dim1];
	    --l;
/* L370: */
	}
L375:
	i__4 = n8;
	for (it = 1; it <= i__4; ++it) {
/*  fetch a new row of matrix b and store it in the arrays h1 (the part */
/*  with respect to g1) and h2 (the part with respect to g2). */
	    i__3 = *idim;
	    for (j = 1; j <= i__3; ++j) {
		xi[j - 1] = 0.f;
/* L380: */
	    }
	    i__3 = *k1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		h1[i__ - 1] = 0.f;
		h2[i__ - 1] = 0.f;
/* L385: */
	    }
	    h1[*k2 - 1] = 0.f;
	    if (it > n11) {
		goto L420;
	    }
	    l = it;
	    l0 = it;
	    i__3 = *k2;
	    for (j = 1; j <= i__3; ++j) {
		if (l0 == n10) {
		    goto L400;
		}
		h1[j - 1] = b[it + j * b_dim1] * pinv;
		++l0;
/* L390: */
	    }
	    goto L470;
L400:
	    l0 = 1;
	    i__3 = *k2;
	    for (l1 = j; l1 <= i__3; ++l1) {
		h2[l0 - 1] = b[it + l1 * b_dim1] * pinv;
		++l0;
/* L410: */
	    }
	    goto L470;
L420:
	    l = 1;
	    i__ = it - n10;
	    i__3 = *k2;
	    for (j = 1; j <= i__3; ++j) {
		++i__;
		l0 = i__;
L430:
		l1 = l0 - *k1;
		if (l1 <= 0) {
		    goto L450;
		}
		if (l1 <= n11) {
		    goto L440;
		}
		l0 = l1 - n11;
		goto L430;
L440:
		h1[l1 - 1] = b[it + j * b_dim1] * pinv;
		goto L460;
L450:
		h2[l0 - 1] += b[it + j * b_dim1] * pinv;
L460:
		;
	    }
	    if (n11 <= 0) {
		goto L510;
	    }
/*  rotate this row into triangle by givens transformations */
/*  rotation with the rows l,l+1,...n11. */
L470:
	    i__3 = n11;
	    for (j = l; j <= i__3; ++j) {
		piv = h1[0];
/*  calculate the parameters of the givens transformation. */
		fpgivs_(&piv, &g1[j + g1_dim1], &cos__, &sin__);
/*  transformation to right hand side. */
		j1 = j;
		i__2 = *idim;
		for (j2 = 1; j2 <= i__2; ++j2) {
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &c__[j1]);
		    j1 += *n;
/* L475: */
		}
/*  transformation to the left hand side with respect to g2. */
		i__2 = *k1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &g2[j + i__ * 
			    g2_dim1]);
/* L480: */
		}
		if (j == n11) {
		    goto L510;
		}
/* Computing MIN */
		i__2 = n11 - j;
		i2 = min(i__2,*k1);
/*  transformation to the left hand side with respect to g1. */
		i__2 = i2;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i1 = i__ + 1;
		    fprota_(&cos__, &sin__, &h1[i1 - 1], &g1[j + i1 * g1_dim1]
			    );
		    h1[i__ - 1] = h1[i1 - 1];
/* L490: */
		}
		h1[i1 - 1] = 0.f;
/* L500: */
	    }
/*  rotation with the rows n11+1,...n7 */
L510:
	    i__3 = *k1;
	    for (j = 1; j <= i__3; ++j) {
		ij = n11 + j;
		if (ij <= 0) {
		    goto L530;
		}
		piv = h2[j - 1];
/*  calculate the parameters of the givens transformation */
		fpgivs_(&piv, &g2[ij + j * g2_dim1], &cos__, &sin__);
/*  transformation to the right hand side. */
		j1 = ij;
		i__2 = *idim;
		for (j2 = 1; j2 <= i__2; ++j2) {
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &c__[j1]);
		    j1 += *n;
/* L515: */
		}
		if (j == *k1) {
		    goto L540;
		}
		j1 = j + 1;
/*  transformation to the left hand side. */
		i__2 = *k1;
		for (i__ = j1; i__ <= i__2; ++i__) {
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &g2[ij + i__ * 
			    g2_dim1]);
/* L520: */
		}
L530:
		;
	    }
L540:
	    ;
	}
/*  backward substitution to obtain the b-spline coefficients */
	j1 = 1;
	i__4 = *idim;
	for (j2 = 1; j2 <= i__4; ++j2) {
	    fpbacp_(&g1[g1_offset], &g2[g2_offset], &c__[j1], &n7, k1, &c__[
		    j1], k2, nest);
	    j1 += *n;
/* L542: */
	}
/*  calculate from condition (**) the remaining b-spline coefficients. */
	i__4 = *k;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    j1 = i__;
	    i__3 = *idim;
	    for (j = 1; j <= i__3; ++j) {
		j2 = j1 + n7;
		c__[j2] = c__[j1];
		j1 += *n;
/* L545: */
	    }
/* L547: */
	}
/*  computation of f(p). */
	*fp = 0.f;
	l = *k1;
	jj = 0;
	i__4 = m1;
	for (it = 1; it <= i__4; ++it) {
	    if (u[it] < t[l]) {
		goto L550;
	    }
	    ++l;
L550:
	    l0 = l - *k2;
	    term = 0.f;
	    i__3 = *idim;
	    for (j2 = 1; j2 <= i__3; ++j2) {
		fac = 0.f;
		j1 = l0;
		i__2 = *k1;
		for (j = 1; j <= i__2; ++j) {
		    ++j1;
		    fac += c__[j1] * q[it + j * q_dim1];
/* L560: */
		}
		++jj;
/* Computing 2nd power */
		d__1 = fac - x[jj];
		term += d__1 * d__1;
		l0 += *n;
/* L565: */
	    }
/* Computing 2nd power */
	    d__1 = w[it];
	    *fp += term * (d__1 * d__1);
/* L570: */
	}
/*  test whether the approximation sp(u) is an acceptable solution. */
	fpms = *fp - *s;
	if (abs(fpms) < acc) {
	    goto L660;
	}
/*  test whether the maximal number of iterations is reached. */
	if (iter == *maxit) {
	    goto L600;
	}
/*  carry out one more step of the iteration process. */
	p2 = p;
	f2 = fpms;
	if (ich3 != 0) {
	    goto L580;
	}
	if (f2 - f3 > acc) {
	    goto L575;
	}
/*  our initial choice of p is too large. */
	p3 = p2;
	f3 = f2;
	p *= con4;
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
	goto L595;
L575:
	if (f2 < 0.f) {
	    ich3 = 1;
	}
L580:
	if (ich1 != 0) {
	    goto L590;
	}
	if (f1 - f2 > acc) {
	    goto L585;
	}
/*  our initial choice of p is too small */
	p1 = p2;
	f1 = f2;
	p /= con4;
	if (p3 < 0.f) {
	    goto L595;
	}
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
	goto L595;
L585:
	if (f2 > 0.f) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
L590:
	if (f2 >= f1 || f2 <= f3) {
	    goto L610;
	}
/*  find the new value for p. */
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
L595:
	;
    }
/*  error codes and messages. */
L600:
    *ier = 3;
    goto L660;
L610:
    *ier = 2;
    goto L660;
L620:
    *ier = 1;
    goto L660;
L630:
    *ier = -1;
    goto L660;
L640:
    *ier = -2;
/*  the point (z(1),z(2),...,z(idim)) is a solution of our problem. */
/*  a constant function is a spline of degree k with all b-spline */
/*  coefficients equal to that constant. */
    i__1 = *k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rn = (doublereal) (*k1 - i__);
	t[i__] = u[1] - rn * per;
	j = i__ + *k1;
	rn = (doublereal) (i__ - 1);
	t[j] = u[*m] + rn * per;
/* L650: */
    }
    *n = nmin;
    j1 = 0;
    i__1 = *idim;
    for (j = 1; j <= i__1; ++j) {
	fac = z__[j];
	j2 = j1;
	i__4 = *k1;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    ++j2;
	    c__[j2] = fac;
/* L654: */
	}
	j1 += *n;
/* L658: */
    }
    *fp = fp0;
    fpint[*n] = fp0;
    fpint[*n - 1] = 0.f;
    nrdata[*n] = 0;
L660:
    return 0;
} /* fpclos_ */

