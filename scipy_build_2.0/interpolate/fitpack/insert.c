/* insert.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int insert_(integer *iopt, doublereal *t, integer *n, 
	doublereal *c__, integer *k, doublereal *x, doublereal *tt, integer *
	nn, doublereal *cc, integer *nest, integer *ier)
{
    static integer l, k1, kk, nk;
    extern /* Subroutine */ int fpinst_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *);

/*  subroutine insert inserts a new knot x into a spline function s(x) */
/*  of degree k and calculates the b-spline representation of s(x) with */
/*  respect to the new set of knots. in addition, if iopt.ne.0, s(x) */
/*  will be considered as a periodic spline with period per=t(n-k)-t(k+1) */
/*  satisfying the boundary constraints */
/*       t(i+n-2*k-1) = t(i)+per  ,i=1,2,...,2*k+1 */
/*       c(i+n-2*k-1) = c(i)      ,i=1,2,...,k */
/*  in that case, the knots and b-spline coefficients returned will also */
/*  satisfy these boundary constraints, i.e. */
/*       tt(i+nn-2*k-1) = tt(i)+per  ,i=1,2,...,2*k+1 */
/*       cc(i+nn-2*k-1) = cc(i)      ,i=1,2,...,k */

/*  calling sequence: */
/*     call insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier) */

/*  input parameters: */
/*    iopt : integer flag, specifying whether (iopt.ne.0) or not (iopt=0) */
/*           the given spline must be considered as being periodic. */
/*    t    : array,length nest, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length nest, which contains the b-spline coefficients. */
/*    k    : integer, giving the degree of s(x). */
/*    x    : real, which gives the location of the knot to be inserted. */
/*    nest : integer specifying the dimension of the arrays t,c,tt and cc */
/*           nest > n. */

/*  output parameters: */
/*    tt   : array,length nest, which contains the position of the knots */
/*           after insertion. */
/*    nn   : integer, giving the total number of knots after insertion */
/*    cc   : array,length nest, which contains the b-spline coefficients */
/*           of s(x) with respect to the new set of knots. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    nest > n */
/*    t(k+1) <= x <= t(n-k) */
/*    in case of a periodic spline (iopt.ne.0) there must be */
/*       either at least k interior knots t(j) satisfying t(k+1)<t(j)<=x */
/*       or at least k interior knots t(j) satisfying x<=t(j)<t(n-k) */

/*  other subroutines required: fpinst. */

/*  further comments: */
/*   subroutine insert may be called as follows */
/*        call insert(iopt,t,n,c,k,x,t,n,c,nest,ier) */
/*   in which case the new representation will simply replace the old one */

/*  references : */
/*    boehm w : inserting new knots into b-spline curves. computer aided */
/*              design 12 (1980) 199-201. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : february 2007 (second interval search added) */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --cc;
    --tt;
    --c__;
    --t;

    /* Function Body */
    *ier = 10;
    if (*nest <= *n) {
	goto L40;
    }
    k1 = *k + 1;
    nk = *n - *k;
    if (*x < t[k1] || *x > t[nk]) {
	goto L40;
    }
/*  search for knot interval t(l) <= x < t(l+1). */
    l = k1;
L10:
    if (*x < t[l + 1]) {
	goto L20;
    }
    ++l;
    if (l == nk) {
	goto L14;
    }
    goto L10;
/*  if no interval found above, then reverse the search and */
/*  look for knot interval t(l) < x <= t(l+1). */
L14:
    l = nk - 1;
L16:
    if (*x > t[l]) {
	goto L20;
    }
    --l;
    if (l == *k) {
	goto L40;
    }
    goto L16;
L20:
    if (t[l] >= t[l + 1]) {
	goto L40;
    }
    if (*iopt == 0) {
	goto L30;
    }
    kk = *k << 1;
    if (l <= kk && l >= *n - kk) {
	goto L40;
    }
L30:
    *ier = 0;
/*  insert the new knot. */
    fpinst_(iopt, &t[1], n, &c__[1], k, x, &l, &tt[1], nn, &cc[1], nest);
L40:
    return 0;
} /* insert_ */

