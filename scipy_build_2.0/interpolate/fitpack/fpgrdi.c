/* fpgrdi.f -- translated by f2c (version 20190311).
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
static integer c__5 = 5;
static integer c__4 = 4;

/* Subroutine */ int fpgrdi_(integer *ifsu, integer *ifsv, integer *ifbu, 
	integer *ifbv, integer *iback, doublereal *u, integer *mu, doublereal 
	*v, integer *mv, doublereal *z__, integer *mz, doublereal *dz, 
	integer *iop0, integer *iop1, doublereal *tu, integer *nu, doublereal 
	*tv, integer *nv, doublereal *p, doublereal *c__, integer *nc, 
	doublereal *sq, doublereal *fp, doublereal *fpu, doublereal *fpv, 
	integer *mm, integer *mvnu, doublereal *spu, doublereal *spv, 
	doublereal *right, doublereal *q, doublereal *au, doublereal *av1, 
	doublereal *av2, doublereal *bu, doublereal *bv, doublereal *aa, 
	doublereal *bb, doublereal *cc, doublereal *cosi, integer *nru, 
	integer *nrv)
{
    /* System generated locals */
    integer spu_dim1, spu_offset, spv_dim1, spv_offset, au_dim1, au_offset, 
	    av1_dim1, av1_offset, av2_dim1, av2_offset, bu_dim1, bu_offset, 
	    bv_dim1, bv_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal h__[5];
    static integer i__, j, k, l, i0, i1, i2, i3, j0, j1, k1, k2, l0, l1, l2, 
	    n1;
    static doublereal h1[5], h2[4];
    static integer ic;
    static doublereal co;
    static integer ii, ij, ik, iq;
    static doublereal si;
    static integer it, jj, jk, iz;
    static doublereal dz1, dz2, dz3;
    static integer nu4, nv4, nu7, nu8, nu9, nv7, nv8;
    static doublereal fac, arg, one;
    static integer nv11;
    static doublereal piv, fac0;
    static integer mvv, nuu;
    static doublereal half;
    static integer ncof, jper;
    static doublereal term, pinv;
    static integer irot, numu, numv, numu1, numv1;
    static doublereal three;
    static integer nrold;
    extern /* Subroutine */ int fpcyt1_(doublereal *, integer *, integer *), 
	    fpcyt2_(doublereal *, integer *, doublereal *, doublereal *, 
	    integer *), fpback_(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *), fpbacp_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), fpdisc_(doublereal *, integer *, integer *,
	     doublereal *, integer *), fpbspl_(doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer number;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer nroldu, nroldv;

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota */
/*  .. */
/*  let */
/*               |   (spu)    |            |   (spv)    | */
/*        (au) = | ---------- |     (av) = | ---------- | */
/*               | (1/p) (bu) |            | (1/p) (bv) | */

/*                                | z  ' 0 | */
/*                            q = | ------ | */
/*                                | 0  ' 0 | */

/*  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline */
/*                coefficients. */
/*       z      : the mu x mv matrix which contains the function values. */
/*       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices */
/*                according to the least-squares problems in the u-,resp. */
/*                v-direction. */
/*       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices */
/*                containing the discontinuity jumps of the derivatives */
/*                of the b-splines in the u-,resp.v-variable at the knots */
/*  the b-spline coefficients of the smoothing spline are then calculated */
/*  as the least-squares solution of the following over-determined linear */
/*  system of equations */

/*    (1)  (av) c (au)' = q */

/*  subject to the constraints */

/*    (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4 */

/*    (3)  if iop0 = 0  c(1,j) = dz(1) */
/*            iop0 = 1  c(1,j) = dz(1) */
/*                      c(2,j) = dz(1)+(dz(2)*cosi(1,j)+dz(3)*cosi(2,j))* */
/*                               tu(5)/3. = cc(j) , j=1,2,...nv-4 */

/*    (4)  if iop1 = 1  c(nu-4,j) = 0, j=1,2,...,nv-4. */

/*  set constants */
    /* Parameter adjustments */
    --nru;
    spu_dim1 = *mu;
    spu_offset = 1 + spu_dim1;
    spu -= spu_offset;
    --u;
    --nrv;
    aa -= 3;
    spv_dim1 = *mv;
    spv_offset = 1 + spv_dim1;
    spv -= spv_offset;
    --v;
    --z__;
    --dz;
    bu_dim1 = *nu;
    bu_offset = 1 + bu_dim1;
    bu -= bu_offset;
    au_dim1 = *nu;
    au_offset = 1 + au_dim1;
    au -= au_offset;
    --fpu;
    --tu;
    cosi -= 3;
    --cc;
    bb -= 3;
    bv_dim1 = *nv;
    bv_offset = 1 + bv_dim1;
    bv -= bv_offset;
    av2_dim1 = *nv;
    av2_offset = 1 + av2_dim1;
    av2 -= av2_offset;
    av1_dim1 = *nv;
    av1_offset = 1 + av1_dim1;
    av1 -= av1_offset;
    --fpv;
    --tv;
    --c__;
    --right;
    --q;

    /* Function Body */
    one = 1.;
    three = 3.;
    half = .5f;
/*  initialization */
    nu4 = *nu - 4;
    nu7 = *nu - 7;
    nu8 = *nu - 8;
    nu9 = *nu - 9;
    nv4 = *nv - 4;
    nv7 = *nv - 7;
    nv8 = *nv - 8;
    nv11 = *nv - 11;
    nuu = nu4 - *iop0 - *iop1 - 1;
    if (*p > 0.f) {
	pinv = one / *p;
    }
/*  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv and iop0 and */
/*  on the value of p whether the matrices (spu), (spv), (bu), (bv) and */
/*  (cosi) still must be determined. */
    if (*ifsu != 0) {
	goto L30;
    }
/*  calculate the non-zero elements of the matrix (spu) which is the ob- */
/*  servation matrix according to the least-squares spline approximation */
/*  problem in the u-direction. */
    l = 4;
    l1 = 5;
    number = 0;
    i__1 = *mu;
    for (it = 1; it <= i__1; ++it) {
	arg = u[it];
L10:
	if (arg < tu[l1] || l == nu4) {
	    goto L15;
	}
	l = l1;
	l1 = l + 1;
	++number;
	goto L10;
L15:
	fpbspl_(&tu[1], nu, &c__3, &arg, &l, h__);
	for (i__ = 1; i__ <= 4; ++i__) {
	    spu[it + i__ * spu_dim1] = h__[i__ - 1];
/* L20: */
	}
	nru[it] = number;
/* L25: */
    }
    *ifsu = 1;
/*  calculate the non-zero elements of the matrix (spv) which is the ob- */
/*  servation matrix according to the least-squares spline approximation */
/*  problem in the v-direction. */
L30:
    if (*ifsv != 0) {
	goto L85;
    }
    l = 4;
    l1 = 5;
    number = 0;
    i__1 = *mv;
    for (it = 1; it <= i__1; ++it) {
	arg = v[it];
L35:
	if (arg < tv[l1] || l == nv4) {
	    goto L40;
	}
	l = l1;
	l1 = l + 1;
	++number;
	goto L35;
L40:
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
	for (i__ = 1; i__ <= 4; ++i__) {
	    spv[it + i__ * spv_dim1] = h__[i__ - 1];
/* L45: */
	}
	nrv[it] = number;
/* L50: */
    }
    *ifsv = 1;
    if (*iop0 == 0) {
	goto L85;
    }
/*  calculate the coefficients of the interpolating splines for cos(v) */
/*  and sin(v). */
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cosi[(i__ << 1) + 1] = 0.f;
	cosi[(i__ << 1) + 2] = 0.f;
/* L55: */
    }
    if (nv7 < 4) {
	goto L85;
    }
    i__1 = nv7;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = i__ + 3;
	arg = tv[l];
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
	for (j = 1; j <= 3; ++j) {
	    av1[i__ + j * av1_dim1] = h__[j - 1];
/* L60: */
	}
	cosi[(i__ << 1) + 1] = cos(arg);
	cosi[(i__ << 1) + 2] = sin(arg);
/* L65: */
    }
    fpcyt1_(&av1[av1_offset], &nv7, nv);
    for (j = 1; j <= 2; ++j) {
	i__1 = nv7;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    right[i__] = cosi[j + (i__ << 1)];
/* L70: */
	}
	fpcyt2_(&av1[av1_offset], &nv7, &right[1], &right[1], nv);
	i__1 = nv7;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cosi[j + (i__ + 1 << 1)] = right[i__];
/* L75: */
	}
	cosi[j + 2] = cosi[j + (nv7 + 1 << 1)];
	cosi[j + (nv7 + 2 << 1)] = cosi[j + 4];
	cosi[j + (nv4 << 1)] = cosi[j + 6];
/* L80: */
    }
L85:
    if (*p <= 0.f) {
	goto L150;
    }
/*  calculate the non-zero elements of the matrix (bu). */
    if (*ifbu != 0 || nu8 == 0) {
	goto L90;
    }
    fpdisc_(&tu[1], nu, &c__5, &bu[bu_offset], nu);
    *ifbu = 1;
/*  calculate the non-zero elements of the matrix (bv). */
L90:
    if (*ifbv != 0 || nv8 == 0) {
	goto L150;
    }
    fpdisc_(&tv[1], nv, &c__5, &bv[bv_offset], nv);
    *ifbv = 1;
/*  substituting (2),(3) and (4) into (1), we obtain the overdetermined */
/*  system */
/*         (5)  (avv) (cr) (auu)' = (qq) */
/*  from which the nuu*nv7 remaining coefficients */
/*         c(i,j) , i=2+iop0,3+iop0,...,nu-4-iop1 ; j=1,2,...,nv-7 , */
/*  the elements of (cr), are then determined in the least-squares sense. */
/*  simultaneously, we compute the resulting sum of squared residuals sq. */
L150:
    dz1 = dz[1];
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aa[(i__ << 1) + 1] = dz1;
/* L155: */
    }
    if (nv8 == 0 || *p <= 0.f) {
	goto L165;
    }
    i__1 = nv8;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bb[(i__ << 1) + 1] = 0.f;
/* L160: */
    }
L165:
    mvv = *mv;
    if (*iop0 == 0) {
	goto L220;
    }
    fac = tu[5] / three;
    dz2 = dz[2] * fac;
    dz3 = dz[3] * fac;
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cc[i__] = dz1 + dz2 * cosi[(i__ << 1) + 1] + dz3 * cosi[(i__ << 1) + 
		2];
/* L170: */
    }
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	number = nrv[i__];
	fac = 0.f;
	for (j = 1; j <= 4; ++j) {
	    ++number;
	    fac += cc[number] * spv[i__ + j * spv_dim1];
/* L180: */
	}
	aa[(i__ << 1) + 2] = fac;
/* L190: */
    }
    if (nv8 == 0 || *p <= 0.f) {
	goto L220;
    }
    i__1 = nv8;
    for (i__ = 1; i__ <= i__1; ++i__) {
	number = i__;
	fac = 0.f;
	for (j = 1; j <= 5; ++j) {
	    fac += cc[number] * bv[i__ + j * bv_dim1];
	    ++number;
/* L200: */
	}
	bb[(i__ << 1) + 2] = fac * pinv;
/* L210: */
    }
    mvv += nv8;
/*  we first determine the matrices (auu) and (qq). then we reduce the */
/*  matrix (auu) to upper triangular form (ru) using givens rotations. */
/*  we apply the same transformations to the rows of matrix qq to obtain */
/*  the (mv+nv8) x nuu matrix g. */
/*  we store matrix (ru) into au and g into q. */
L220:
    l = mvv * nuu;
/*  initialization. */
    *sq = 0.f;
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = 0.f;
/* L230: */
    }
    i__1 = nuu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 5; ++j) {
	    au[i__ + j * au_dim1] = 0.f;
/* L240: */
	}
    }
    l = 0;
    nrold = 0;
    n1 = nrold + 1;
    i__1 = *mu;
    for (it = 1; it <= i__1; ++it) {
	number = nru[it];
/*  find the appropriate column of q. */
L250:
	i__2 = mvv;
	for (j = 1; j <= i__2; ++j) {
	    right[j] = 0.f;
/* L260: */
	}
	if (nrold == number) {
	    goto L280;
	}
	if (*p <= 0.f) {
	    goto L410;
	}
/*  fetch a new row of matrix (bu). */
	for (j = 1; j <= 5; ++j) {
	    h__[j - 1] = bu[n1 + j * bu_dim1] * pinv;
/* L270: */
	}
	i0 = 1;
	i1 = 5;
	goto L310;
/*  fetch a new row of matrix (spu). */
L280:
	for (j = 1; j <= 4; ++j) {
	    h__[j - 1] = spu[it + j * spu_dim1];
/* L290: */
	}
/*  find the appropriate column of q. */
	i__2 = *mv;
	for (j = 1; j <= i__2; ++j) {
	    ++l;
	    right[j] = z__[l];
/* L300: */
	}
	i0 = 1;
	i1 = 4;
L310:
	if (nu7 - number == *iop1) {
	    --i1;
	}
	j0 = n1;
/*  take into account that we eliminate the constraints (3) */
L320:
	if (j0 - 1 > *iop0) {
	    goto L360;
	}
	fac0 = h__[i0 - 1];
	i__2 = *mv;
	for (j = 1; j <= i__2; ++j) {
	    right[j] -= fac0 * aa[j0 + (j << 1)];
/* L330: */
	}
	if (*mv == mvv) {
	    goto L350;
	}
	j = *mv;
	i__2 = nv8;
	for (jj = 1; jj <= i__2; ++jj) {
	    ++j;
	    right[j] -= fac0 * bb[j0 + (jj << 1)];
/* L340: */
	}
L350:
	++j0;
	++i0;
	goto L320;
L360:
	irot = nrold - *iop0 - 1;
	if (irot < 0) {
	    irot = 0;
	}
/*  rotate the new row of matrix (auu) into triangle. */
	i__2 = i1;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    ++irot;
	    piv = h__[i__ - 1];
	    if (piv == 0.f) {
		goto L390;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &au[irot + au_dim1], &co, &si);
/*  apply that transformation to the rows of matrix (qq). */
	    iq = (irot - 1) * mvv;
	    i__3 = mvv;
	    for (j = 1; j <= i__3; ++j) {
		++iq;
		fprota_(&co, &si, &right[j], &q[iq]);
/* L370: */
	    }
/*  apply that transformation to the columns of (auu). */
	    if (i__ == i1) {
		goto L390;
	    }
	    i2 = 1;
	    i3 = i__ + 1;
	    i__3 = i1;
	    for (j = i3; j <= i__3; ++j) {
		++i2;
		fprota_(&co, &si, &h__[j - 1], &au[irot + i2 * au_dim1]);
/* L380: */
	    }
L390:
	    ;
	}
/* we update the sum of squared residuals */
	i__2 = mvv;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = right[j];
	    *sq += d__1 * d__1;
/* L395: */
	}
	if (nrold == number) {
	    goto L420;
	}
L410:
	nrold = n1;
	++n1;
	goto L250;
L420:
	;
    }
/*  we determine the matrix (avv) and then we reduce her to */
/*  upper triangular form (rv) using givens rotations. */
/*  we apply the same transformations to the columns of matrix */
/*  g to obtain the (nv-7) x (nu-5-iop0-iop1) matrix h. */
/*  we store matrix (rv) into av1 and av2, h into c. */
/*  the nv7 x nv7 upper triangular matrix (rv) has the form */
/*              | av1 '     | */
/*       (rv) = |     ' av2 | */
/*              |  0  '     | */
/*  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 upper */
/*  triangular matrix of bandwidth 5. */
    ncof = nuu * nv7;
/*  initialization. */
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = 0.f;
/* L430: */
    }
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	av1[i__ + av1_dim1 * 5] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    av1[i__ + j * av1_dim1] = 0.f;
	    av2[i__ + j * av2_dim1] = 0.f;
/* L440: */
	}
    }
    jper = 0;
    nrold = 0;
    i__1 = *mv;
    for (it = 1; it <= i__1; ++it) {
	number = nrv[it];
L450:
	if (nrold == number) {
	    goto L480;
	}
	if (*p <= 0.f) {
	    goto L760;
	}
/*  fetch a new row of matrix (bv). */
	n1 = nrold + 1;
	for (j = 1; j <= 5; ++j) {
	    h__[j - 1] = bv[n1 + j * bv_dim1] * pinv;
/* L460: */
	}
/*  find the appropriate row of g. */
	i__2 = nuu;
	for (j = 1; j <= i__2; ++j) {
	    right[j] = 0.f;
/* L465: */
	}
	if (*mv == mvv) {
	    goto L510;
	}
	l = *mv + n1;
	i__2 = nuu;
	for (j = 1; j <= i__2; ++j) {
	    right[j] = q[l];
	    l += mvv;
/* L470: */
	}
	goto L510;
/*  fetch a new row of matrix (spv) */
L480:
	h__[4] = 0.f;
	for (j = 1; j <= 4; ++j) {
	    h__[j - 1] = spv[it + j * spv_dim1];
/* L490: */
	}
/*  find the appropriate row of g. */
	l = it;
	i__2 = nuu;
	for (j = 1; j <= i__2; ++j) {
	    right[j] = q[l];
	    l += mvv;
/* L500: */
	}
/*  test whether there are non-zero values in the new row of (avv) */
/*  corresponding to the b-splines n(j,v),j=nv7+1,...,nv4. */
L510:
	if (nrold < nv11) {
	    goto L710;
	}
	if (jper != 0) {
	    goto L550;
	}
/*  initialize the matrix (av2). */
	jk = nv11 + 1;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ik = jk;
	    for (j = 1; j <= 5; ++j) {
		if (ik <= 0) {
		    goto L530;
		}
		av2[ik + i__ * av2_dim1] = av1[ik + j * av1_dim1];
		--ik;
/* L520: */
	    }
L530:
	    ++jk;
/* L540: */
	}
	jper = 1;
/*  if one of the non-zero elements of the new row corresponds to one of */
/*  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition */
/*  (2) for setting up this row of (avv). the row is stored in h1( the */
/*  part with respect to av1) and h2 (the part with respect to av2). */
L550:
	for (i__ = 1; i__ <= 4; ++i__) {
	    h1[i__ - 1] = 0.f;
	    h2[i__ - 1] = 0.f;
/* L560: */
	}
	h1[4] = 0.f;
	j = nrold - nv11;
	for (i__ = 1; i__ <= 5; ++i__) {
	    ++j;
	    l0 = j;
L570:
	    l1 = l0 - 4;
	    if (l1 <= 0) {
		goto L590;
	    }
	    if (l1 <= nv11) {
		goto L580;
	    }
	    l0 = l1 - nv11;
	    goto L570;
L580:
	    h1[l1 - 1] = h__[i__ - 1];
	    goto L600;
L590:
	    h2[l0 - 1] += h__[i__ - 1];
L600:
	    ;
	}
/*  rotate the new row of (avv) into triangle. */
	if (nv11 <= 0) {
	    goto L670;
	}
/*  rotations with the rows 1,2,...,nv11 of (avv). */
	i__2 = nv11;
	for (j = 1; j <= i__2; ++j) {
	    piv = h1[0];
/* Computing MIN */
	    i__3 = nv11 - j;
	    i2 = min(i__3,4);
	    if (piv == 0.f) {
		goto L640;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &av1[j + av1_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
	    ic = j;
	    i__3 = nuu;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		fprota_(&co, &si, &right[i__], &c__[ic]);
		ic += nv7;
/* L610: */
	    }
/*  apply that transformation to the rows of (avv) with respect to av2. */
	    for (i__ = 1; i__ <= 4; ++i__) {
		fprota_(&co, &si, &h2[i__ - 1], &av2[j + i__ * av2_dim1]);
/* L620: */
	    }
/*  apply that transformation to the rows of (avv) with respect to av1. */
	    if (i2 == 0) {
		goto L670;
	    }
	    i__3 = i2;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i1 = i__ + 1;
		fprota_(&co, &si, &h1[i1 - 1], &av1[j + i1 * av1_dim1]);
/* L630: */
	    }
L640:
	    i__3 = i2;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		h1[i__ - 1] = h1[i__];
/* L650: */
	    }
	    h1[i2] = 0.f;
/* L660: */
	}
/*  rotations with the rows nv11+1,...,nv7 of avv. */
L670:
	for (j = 1; j <= 4; ++j) {
	    ij = nv11 + j;
	    if (ij <= 0) {
		goto L700;
	    }
	    piv = h2[j - 1];
	    if (piv == 0.f) {
		goto L700;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &av2[ij + j * av2_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
	    ic = ij;
	    i__2 = nuu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		fprota_(&co, &si, &right[i__], &c__[ic]);
		ic += nv7;
/* L680: */
	    }
	    if (j == 4) {
		goto L700;
	    }
/*  apply that transformation to the rows of (avv) with respect to av2. */
	    j1 = j + 1;
	    for (i__ = j1; i__ <= 4; ++i__) {
		fprota_(&co, &si, &h2[i__ - 1], &av2[ij + i__ * av2_dim1]);
/* L690: */
	    }
L700:
	    ;
	}
/* we update the sum of squared residuals */
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = right[i__];
	    *sq += d__1 * d__1;
/* L705: */
	}
	goto L750;
/*  rotation into triangle of the new row of (avv), in case the elements */
/*  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero. */
L710:
	irot = nrold;
	for (i__ = 1; i__ <= 5; ++i__) {
	    ++irot;
	    piv = h__[i__ - 1];
	    if (piv == 0.f) {
		goto L740;
	    }
/*  calculate the parameters of the givens transformation. */
	    fpgivs_(&piv, &av1[irot + av1_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
	    ic = irot;
	    i__2 = nuu;
	    for (j = 1; j <= i__2; ++j) {
		fprota_(&co, &si, &right[j], &c__[ic]);
		ic += nv7;
/* L720: */
	    }
/*  apply that transformation to the rows of (avv). */
	    if (i__ == 5) {
		goto L740;
	    }
	    i2 = 1;
	    i3 = i__ + 1;
	    for (j = i3; j <= 5; ++j) {
		++i2;
		fprota_(&co, &si, &h__[j - 1], &av1[irot + i2 * av1_dim1]);
/* L730: */
	    }
L740:
	    ;
	}
/* we update the sum of squared residuals */
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = right[i__];
	    *sq += d__1 * d__1;
/* L745: */
	}
L750:
	if (nrold == number) {
	    goto L770;
	}
L760:
	++nrold;
	goto L450;
L770:
	;
    }
/*  test whether the b-spline coefficients must be determined. */
    if (*iback != 0) {
	return 0;
    }
/*  backward substitution to obtain the b-spline coefficients as the */
/*  solution of the linear system    (rv) (cr) (ru)' = h. */
/*  first step: solve the system  (rv) (c1) = h. */
    k = 1;
    i__1 = nuu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fpbacp_(&av1[av1_offset], &av2[av2_offset], &c__[k], &nv7, &c__4, &
		c__[k], &c__5, nv);
	k += nv7;
/* L780: */
    }
/*  second step: solve the system  (cr) (ru)' = (c1). */
    k = 0;
    i__1 = nv7;
    for (j = 1; j <= i__1; ++j) {
	++k;
	l = k;
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    right[i__] = c__[l];
	    l += nv7;
/* L790: */
	}
	fpback_(&au[au_offset], &right[1], &nuu, &c__5, &right[1], nu);
	l = k;
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[l] = right[i__];
	    l += nv7;
/* L795: */
	}
/* L800: */
    }
/*  calculate from the conditions (2)-(3)-(4), the remaining b-spline */
/*  coefficients. */
    ncof = nu4 * nv4;
    i__ = nv4;
    j = 0;
    i__1 = nv4;
    for (l = 1; l <= i__1; ++l) {
	q[l] = dz1;
/* L805: */
    }
    if (*iop0 == 0) {
	goto L815;
    }
    i__1 = nv4;
    for (l = 1; l <= i__1; ++l) {
	++i__;
	q[i__] = cc[l];
/* L810: */
    }
L815:
    if (nuu == 0) {
	goto L850;
    }
    i__1 = nuu;
    for (l = 1; l <= i__1; ++l) {
	ii = i__;
	i__2 = nv7;
	for (k = 1; k <= i__2; ++k) {
	    ++i__;
	    ++j;
	    q[i__] = c__[j];
/* L820: */
	}
	for (k = 1; k <= 3; ++k) {
	    ++ii;
	    ++i__;
	    q[i__] = q[ii];
/* L830: */
	}
/* L840: */
    }
L850:
    if (*iop1 == 0) {
	goto L870;
    }
    i__1 = nv4;
    for (l = 1; l <= i__1; ++l) {
	++i__;
	q[i__] = 0.f;
/* L860: */
    }
L870:
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = q[i__];
/* L880: */
    }
/*  calculate the quantities */
/*    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv */
/*    fp = sumi=1,mu(sumj=1,mv(res(i,j))) */
/*    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7 */
/*                  tu(r+3) <= u(i) <= tu(r+4) */
/*    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7 */
/*                  tv(r+3) <= v(j) <= tv(r+4) */
    *fp = 0.f;
    i__1 = *nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fpu[i__] = 0.f;
/* L890: */
    }
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fpv[i__] = 0.f;
/* L900: */
    }
    iz = 0;
    nroldu = 0;
/*  main loop for the different grid points. */
    i__1 = *mu;
    for (i1 = 1; i1 <= i__1; ++i1) {
	numu = nru[i1];
	numu1 = numu + 1;
	nroldv = 0;
	i__2 = *mv;
	for (i2 = 1; i2 <= i__2; ++i2) {
	    numv = nrv[i2];
	    numv1 = numv + 1;
	    ++iz;
/*  evaluate s(u,v) at the current grid point by making the sum of the */
/*  cross products of the non-zero b-splines at (u,v), multiplied with */
/*  the appropriate b-spline coefficients. */
	    term = 0.f;
	    k1 = numu * nv4 + numv;
	    for (l1 = 1; l1 <= 4; ++l1) {
		k2 = k1;
		fac = spu[i1 + l1 * spu_dim1];
		for (l2 = 1; l2 <= 4; ++l2) {
		    ++k2;
		    term += fac * spv[i2 + l2 * spv_dim1] * c__[k2];
/* L910: */
		}
		k1 += nv4;
/* L920: */
	    }
/*  calculate the squared residual at the current grid point. */
/* Computing 2nd power */
	    d__1 = z__[iz] - term;
	    term = d__1 * d__1;
/*  adjust the different parameters. */
	    *fp += term;
	    fpu[numu1] += term;
	    fpv[numv1] += term;
	    fac = term * half;
	    if (numv == nroldv) {
		goto L930;
	    }
	    fpv[numv1] -= fac;
	    fpv[numv] += fac;
L930:
	    nroldv = numv;
	    if (numu == nroldu) {
		goto L940;
	    }
	    fpu[numu1] -= fac;
	    fpu[numu] += fac;
L940:
	    ;
	}
	nroldu = numu;
/* L950: */
    }
    return 0;
} /* fpgrdi_ */

