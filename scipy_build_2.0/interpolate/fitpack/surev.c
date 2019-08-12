/* surev.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int surev_(integer *idim, doublereal *tu, integer *nu, 
	doublereal *tv, integer *nv, doublereal *c__, doublereal *u, integer *
	mu, doublereal *v, integer *mv, doublereal *f, integer *mf, 
	doublereal *wrk, integer *lwrk, integer *iwrk, integer *kwrk, integer 
	*ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, muv;
    extern /* Subroutine */ int fpsuev_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *);

/*  subroutine surev evaluates on a grid (u(i),v(j)),i=1,...,mu; j=1,... */
/*  ,mv a bicubic spline surface of dimension idim, given in the */
/*  b-spline representation. */

/*  calling sequence: */
/*     call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk, */
/*    * iwrk,kwrk,ier) */

/*  input parameters: */
/*   idim  : integer, specifying the dimension of the spline surface. */
/*   tu    : real array, length nu, which contains the position of the */
/*           knots in the u-direction. */
/*   nu    : integer, giving the total number of knots in the u-direction */
/*   tv    : real array, length nv, which contains the position of the */
/*           knots in the v-direction. */
/*   nv    : integer, giving the total number of knots in the v-direction */
/*   c     : real array, length (nu-4)*(nv-4)*idim, which contains the */
/*           b-spline coefficients. */
/*   u     : real array of dimension (mu). */
/*           before entry u(i) must be set to the u co-ordinate of the */
/*           i-th grid point along the u-axis. */
/*           tu(4)<=u(i-1)<=u(i)<=tu(nu-3), i=2,...,mu. */
/*   mu    : on entry mu must specify the number of grid points along */
/*           the u-axis. mu >=1. */
/*   v     : real array of dimension (mv). */
/*           before entry v(j) must be set to the v co-ordinate of the */
/*           j-th grid point along the v-axis. */
/*           tv(4)<=v(j-1)<=v(j)<=tv(nv-3), j=2,...,mv. */
/*   mv    : on entry mv must specify the number of grid points along */
/*           the v-axis. mv >=1. */
/*   mf    : on entry, mf must specify the dimension of the array f. */
/*           mf >= mu*mv*idim */
/*   wrk   : real array of dimension lwrk. used as workspace. */
/*   lwrk  : integer, specifying the dimension of wrk. */
/*           lwrk >= 4*(mu+mv) */
/*   iwrk  : integer array of dimension kwrk. used as workspace. */
/*   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mu+mv. */

/*  output parameters: */
/*   f     : real array of dimension (mf). */
/*           on successful exit f(mu*mv*(l-1)+mv*(i-1)+j) contains the */
/*           l-th co-ordinate of the bicubic spline surface at the */
/*           point (u(i),v(j)),l=1,...,idim,i=1,...,mu;j=1,...,mv. */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   mu >=1, mv >=1, lwrk>=4*(mu+mv), kwrk>=mu+mv , mf>=mu*mv*idim */
/*   tu(4) <= u(i-1) <= u(i) <= tu(nu-3), i=2,...,mu */
/*   tv(4) <= v(j-1) <= v(j) <= tv(nv-3), j=2,...,mv */

/*  other subroutines required: */
/*    fpsuev,fpbspl */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths */
/*                applics 10 (1972) 134-149. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --tu;
    --c__;
    --tv;
    --u;
    --v;
    --f;
    --wrk;
    --iwrk;

    /* Function Body */
    *ier = 10;
    if (*mf < *mu * *mv * *idim) {
	goto L100;
    }
    muv = *mu + *mv;
    if (*lwrk < muv << 2) {
	goto L100;
    }
    if (*kwrk < muv) {
	goto L100;
    }
    if (*mu < 1) {
	goto L100;
    }
    if (*mu == 1) {
	goto L30;
    }
    goto L10;
L10:
    i__1 = *mu;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (u[i__] < u[i__ - 1]) {
	    goto L100;
	}
/* L20: */
    }
L30:
    if (*mv < 1) {
	goto L100;
    }
    if (*mv == 1) {
	goto L60;
    }
    goto L40;
L40:
    i__1 = *mv;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (v[i__] < v[i__ - 1]) {
	    goto L100;
	}
/* L50: */
    }
L60:
    *ier = 0;
    fpsuev_(idim, &tu[1], nu, &tv[1], nv, &c__[1], &u[1], mu, &v[1], mv, &f[1]
	    , &wrk[1], &wrk[(*mu << 2) + 1], &iwrk[1], &iwrk[*mu + 1]);
L100:
    return 0;
} /* surev_ */

