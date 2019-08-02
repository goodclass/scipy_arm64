/* fpbfout.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpbfou_(doublereal *t, integer *n, doublereal *par, 
	doublereal *ress, doublereal *resc)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal c1, c2, f1, f2, f3, s1, s2, hc[5];
    static integer ic;
    static doublereal ak, co[5];
    static integer jj, li, lj;
    static doublereal rc[3];
    static integer ll;
    static doublereal hs[5];
    static integer is;
    static doublereal si[5], rs[3];
    static integer jp1, jp4, nm3, nm7;
    static doublereal fac, one;
    static integer ipj, nmj;
    static doublereal eps, six, con1, con2, beta, sign, term, delta, quart;
    extern /* Subroutine */ int fpcsin_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/*  subroutine fpbfou calculates the integrals */
/*                    /t(n-3) */
/*    ress(j) =      !        nj,4(x)*sin(par*x) dx    and */
/*              t(4)/ */
/*                    /t(n-3) */
/*    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4 */
/*              t(4)/ */
/*  where nj,4(x) denotes the cubic b-spline defined on the knots */
/*  t(j),t(j+1),...,t(j+4). */

/*  calling sequence: */
/*     call fpbfou(t,n,par,ress,resc) */

/*  input parameters: */
/*    t    : real array,length n, containing the knots. */
/*    n    : integer, containing the number of knots. */
/*    par  : real, containing the value of the parameter par. */

/*  output parameters: */
/*    ress  : real array,length n, containing the integrals ress(j). */
/*    resc  : real array,length n, containing the integrals resc(j). */

/*  restrictions: */
/*    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3). */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  .. */
/*  initialization. */
    /* Parameter adjustments */
    --resc;
    --ress;
    --t;

    /* Function Body */
    one = 1.f;
    six = 6.f;
    eps = 1e-8f;
    quart = .25f;
    con1 = .05f;
    con2 = 120.f;
    nm3 = *n - 3;
    nm7 = *n - 7;
    if (*par != 0.f) {
	term = six / *par;
    }
    beta = *par * t[4];
    co[0] = cos(beta);
    si[0] = sin(beta);
/*  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up */
/*  a divided difference table. */
    for (j = 1; j <= 3; ++j) {
	jp1 = j + 1;
	jp4 = j + 4;
	beta = *par * t[jp4];
	co[jp1 - 1] = cos(beta);
	si[jp1 - 1] = sin(beta);
	fpcsin_(&t[4], &t[jp4], par, si, co, &si[jp1 - 1], &co[jp1 - 1], &rs[
		j - 1], &rc[j - 1]);
	i__ = 5 - j;
	hs[i__ - 1] = 0.f;
	hc[i__ - 1] = 0.f;
	i__1 = j;
	for (jj = 1; jj <= i__1; ++jj) {
	    ipj = i__ + jj;
	    hs[ipj - 1] = rs[jj - 1];
	    hc[ipj - 1] = rc[jj - 1];
/* L10: */
	}
	for (jj = 1; jj <= 3; ++jj) {
	    if (i__ < jj) {
		i__ = jj;
	    }
	    k = 5;
	    li = jp4;
	    for (ll = i__; ll <= 4; ++ll) {
		lj = li - jj;
		fac = t[li] - t[lj];
		hs[k - 1] = (hs[k - 1] - hs[k - 2]) / fac;
		hc[k - 1] = (hc[k - 1] - hc[k - 2]) / fac;
		--k;
		--li;
/* L20: */
	    }
	}
	ress[j] = hs[4] - hs[3];
	resc[j] = hc[4] - hc[3];
/* L30: */
    }
    if (nm7 < 4) {
	goto L160;
    }
/*  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7. */
    i__1 = nm7;
    for (j = 4; j <= i__1; ++j) {
	jp4 = j + 4;
	beta = *par * t[jp4];
	co[4] = cos(beta);
	si[4] = sin(beta);
	delta = t[jp4] - t[j];
/*  the way of computing ress(j) and resc(j) depends on the value of */
/*  beta = par*(t(j+4)-t(j)). */
	beta = delta * *par;
	if (abs(beta) <= one) {
	    goto L60;
	}
/*  if !beta! > 1 the integrals are calculated by setting up a divided */
/*  difference table. */
	for (k = 1; k <= 5; ++k) {
	    hs[k - 1] = si[k - 1];
	    hc[k - 1] = co[k - 1];
/* L40: */
	}
	for (jj = 1; jj <= 3; ++jj) {
	    k = 5;
	    li = jp4;
	    for (ll = jj; ll <= 4; ++ll) {
		lj = li - jj;
		fac = *par * (t[li] - t[lj]);
		hs[k - 1] = (hs[k - 1] - hs[k - 2]) / fac;
		hc[k - 1] = (hc[k - 1] - hc[k - 2]) / fac;
		--k;
		--li;
/* L50: */
	    }
	}
	s2 = (hs[4] - hs[3]) * term;
	c2 = (hc[4] - hc[3]) * term;
	goto L130;
/*  if !beta! <= 1 the integrals are calculated by evaluating a series */
/*  expansion. */
L60:
	f3 = 0.f;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipj = i__ + j;
	    hs[i__ - 1] = *par * (t[ipj] - t[j]);
	    hc[i__ - 1] = hs[i__ - 1];
	    f3 += hs[i__ - 1];
/* L70: */
	}
	f3 *= con1;
	c1 = quart;
	s1 = f3;
	if (abs(f3) <= eps) {
	    goto L120;
	}
	sign = one;
	fac = con2;
	k = 5;
	is = 0;
	for (ic = 1; ic <= 20; ++ic) {
	    ++k;
	    ak = (doublereal) k;
	    fac *= ak;
	    f1 = 0.f;
	    f3 = 0.f;
	    for (i__ = 1; i__ <= 4; ++i__) {
		f1 += hc[i__ - 1];
		f2 = f1 * hs[i__ - 1];
		hc[i__ - 1] = f2;
		f3 += f2;
/* L80: */
	    }
	    f3 = f3 * six / fac;
	    if (is == 0) {
		goto L90;
	    }
	    is = 0;
	    s1 += f3 * sign;
	    goto L100;
L90:
	    sign = -sign;
	    is = 1;
	    c1 += f3 * sign;
L100:
	    if (abs(f3) <= eps) {
		goto L120;
	    }
/* L110: */
	}
L120:
	s2 = delta * (co[0] * s1 + si[0] * c1);
	c2 = delta * (co[0] * c1 - si[0] * s1);
L130:
	ress[j] = s2;
	resc[j] = c2;
	for (i__ = 1; i__ <= 4; ++i__) {
	    co[i__ - 1] = co[i__];
	    si[i__ - 1] = si[i__];
/* L140: */
	}
/* L150: */
    }
/*  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting */
/*  up a divided difference table. */
L160:
    for (j = 1; j <= 3; ++j) {
	nmj = nm3 - j;
	i__ = 5 - j;
	fpcsin_(&t[nm3], &t[nmj], par, &si[3], &co[3], &si[i__ - 2], &co[i__ 
		- 2], &rs[j - 1], &rc[j - 1]);
	hs[i__ - 1] = 0.f;
	hc[i__ - 1] = 0.f;
	i__1 = j;
	for (jj = 1; jj <= i__1; ++jj) {
	    ipj = i__ + jj;
	    hc[ipj - 1] = rc[jj - 1];
	    hs[ipj - 1] = rs[jj - 1];
/* L170: */
	}
	for (jj = 1; jj <= 3; ++jj) {
	    if (i__ < jj) {
		i__ = jj;
	    }
	    k = 5;
	    li = nmj;
	    for (ll = i__; ll <= 4; ++ll) {
		lj = li + jj;
		fac = t[lj] - t[li];
		hs[k - 1] = (hs[k - 2] - hs[k - 1]) / fac;
		hc[k - 1] = (hc[k - 2] - hc[k - 1]) / fac;
		--k;
		++li;
/* L180: */
	    }
	}
	ress[nmj] = hs[3] - hs[4];
	resc[nmj] = hc[3] - hc[4];
/* L190: */
    }
    return 0;
} /* fpbfou_ */

