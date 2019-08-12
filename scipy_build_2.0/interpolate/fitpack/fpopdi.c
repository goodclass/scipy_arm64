/* fpopdi.f -- translated by f2c (version 20190311).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int fpopdi_(integer *ifsu, integer *ifsv, integer *ifbu, 
	integer *ifbv, doublereal *u, integer *mu, doublereal *v, integer *mv,
	 doublereal *z__, integer *mz, doublereal *z0, doublereal *dz, 
	integer *iopt, integer *ider, doublereal *tu, integer *nu, doublereal 
	*tv, integer *nv, integer *nuest, integer *nvest, doublereal *p, 
	doublereal *step, doublereal *c__, integer *nc, doublereal *fp, 
	doublereal *fpu, doublereal *fpv, integer *nru, integer *nrv, 
	doublereal *wrk, integer *lwrk)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal a[36]	/* was [6][6] */, g[6];
    static integer i__, j, l, i1, l1, l2, mm, lq, nr[3];
    static doublereal sq;
    static integer id0, laa, lbb, lcc, lau, lbu, lbv, lcs, lri;
    static doublereal res, sqq;
    static integer lsu, lsv;
    static doublereal dzz[3], sum[3];
    static integer lav1, lav2, iop0, iop1, mvnu;
    static doublereal step1, step2, delta[3], three;
    extern /* Subroutine */ int fpgrdi_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    static integer number;
    extern /* Subroutine */ int fpsysy_(doublereal *, integer *, doublereal *)
	    ;

/*  given the set of function values z(i,j) defined on the rectangular */
/*  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopdi determines a */
/*  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu */
/*  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this */
/*  spline sp(u,v) will be periodic in the variable v and will satisfy */
/*  the following constraints */

/*     s(tu(1),v) = dz(1) , tv(4) <=v<= tv(nv-3) */

/*  and (if iopt(2) = 1) */

/*     d s(tu(1),v) */
/*     ------------ =  dz(2)*cos(v)+dz(3)*sin(v) , tv(4) <=v<= tv(nv-3) */
/*     d u */

/*  and (if iopt(3) = 1) */

/*     s(tu(nu),v)  =  0   tv(4) <=v<= tv(nv-3) */

/*  where the parameters dz(i) correspond to the derivative values g(i,j) */
/*  as defined in subroutine pogrid. */

/*  the b-spline coefficients of sp(u,v) are determined as the least- */
/*  squares solution  of an overdetermined linear system which depends */
/*  on the value of p and on the values dz(i),i=1,2,3. the correspond- */
/*  ing sum of squared residuals sq is a simple quadratic function in */
/*  the variables dz(i). these may or may not be provided. the values */
/*  dz(i) which are not given will be determined so as to minimize the */
/*  resulting sum of squared residuals sq. in that case the user must */
/*  provide some initial guess dz(i) and some estimate (dz(i)-step, */
/*  dz(i)+step) of the range of possible values for these latter. */

/*  sp(u,v) also depends on the parameter p (p>0) in such a way that */
/*    - if p tends to infinity, sp(u,v) becomes the least-squares spline */
/*      with given knots, satisfying the constraints. */
/*    - if p tends to zero, sp(u,v) becomes the least-squares polynomial, */
/*      satisfying the constraints. */
/*    - the function  f(p)=sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2) */
/*      is continuous and strictly decreasing for p>0. */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..function references.. */
/*  ..subroutine references.. */
/*    fpgrdi,fpsysy */
/*  .. */
/*  set constant */
    /* Parameter adjustments */
    --nru;
    --u;
    --nrv;
    --v;
    --z__;
    --dz;
    --iopt;
    --ider;
    --fpu;
    --tu;
    --fpv;
    --tv;
    --c__;
    --wrk;

    /* Function Body */
    three = 3.;
/*  we partition the working space */
    lsu = 1;
    lsv = lsu + (*mu << 2);
    lri = lsv + (*mv << 2);
/* Computing MAX */
    i__1 = *nuest, i__2 = *mv + *nvest;
    mm = max(i__1,i__2);
    lq = lri + mm;
    mvnu = *nuest * (*mv + *nvest - 8);
    lau = lq + mvnu;
    lav1 = lau + *nuest * 5;
    lav2 = lav1 + *nvest * 6;
    lbu = lav2 + (*nvest << 2);
    lbv = lbu + *nuest * 5;
    laa = lbv + *nvest * 5;
    lbb = laa + (*mv << 1);
    lcc = lbb + (*nvest << 1);
    lcs = lcc + *nvest;
/*  we calculate the smoothing spline sp(u,v) according to the input */
/*  values dz(i),i=1,2,3. */
    iop0 = iopt[2];
    iop1 = iopt[3];
    fpgrdi_(ifsu, ifsv, ifbu, ifbv, &c__0, &u[1], mu, &v[1], mv, &z__[1], mz, 
	    &dz[1], &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc, &sq,
	     fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv], &wrk[lri]
	    , &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[lav2], &wrk[lbu], &wrk[
	    lbv], &wrk[laa], &wrk[lbb], &wrk[lcc], &wrk[lcs], &nru[1], &nrv[1]
	    );
    id0 = ider[1];
    if (id0 != 0) {
	goto L5;
    }
/* Computing 2nd power */
    d__1 = *z0 - dz[1];
    res = d__1 * d__1;
    *fp += res;
    sq += res;
/* in case all derivative values dz(i) are given (step<=0) or in case */
/* we have spline interpolation, we accept this spline as a solution. */
L5:
    if (*step <= 0.f || sq <= 0.f) {
	return 0;
    }
    dzz[0] = dz[1];
    dzz[1] = dz[2];
    dzz[2] = dz[3];
/* number denotes the number of derivative values dz(i) that still must */
/* be optimized. let us denote these parameters by g(j),j=1,...,number. */
    number = 0;
    if (id0 > 0) {
	goto L10;
    }
    number = 1;
    nr[0] = 1;
    delta[0] = *step;
L10:
    if (iop0 == 0) {
	goto L20;
    }
    if (ider[2] != 0) {
	goto L20;
    }
    step2 = *step * three / tu[5];
    nr[number] = 2;
    nr[number + 1] = 3;
    delta[number] = step2;
    delta[number + 1] = step2;
    number += 2;
L20:
    if (number == 0) {
	return 0;
    }
/* the sum of squared residuals sq is a quadratic polynomial in the */
/* parameters g(j). we determine the unknown coefficients of this */
/* polymomial by calculating (number+1)*(number+2)/2 different splines */
/* according to specific values for g(j). */
    i__1 = number;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = nr[i__ - 1];
	step1 = delta[i__ - 1];
	dzz[l - 1] = dz[l] + step1;
	fpgrdi_(ifsu, ifsv, ifbu, ifbv, &c__1, &u[1], mu, &v[1], mv, &z__[1], 
		mz, dzz, &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc,
		 &sum[i__ - 1], fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &
		wrk[lsv], &wrk[lri], &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[
		lav2], &wrk[lbu], &wrk[lbv], &wrk[laa], &wrk[lbb], &wrk[lcc], 
		&wrk[lcs], &nru[1], &nrv[1]);
	if (id0 == 0) {
/* Computing 2nd power */
	    d__1 = *z0 - dzz[0];
	    sum[i__ - 1] += d__1 * d__1;
	}
	dzz[l - 1] = dz[l] - step1;
	fpgrdi_(ifsu, ifsv, ifbu, ifbv, &c__1, &u[1], mu, &v[1], mv, &z__[1], 
		mz, dzz, &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc,
		 &sqq, fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv],
		 &wrk[lri], &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[lav2], &wrk[
		lbu], &wrk[lbv], &wrk[laa], &wrk[lbb], &wrk[lcc], &wrk[lcs], &
		nru[1], &nrv[1]);
	if (id0 == 0) {
/* Computing 2nd power */
	    d__1 = *z0 - dzz[0];
	    sqq += d__1 * d__1;
	}
/* Computing 2nd power */
	d__1 = step1;
	a[i__ + i__ * 6 - 7] = (sum[i__ - 1] + sqq - sq - sq) / (d__1 * d__1);
	if (a[i__ + i__ * 6 - 7] <= 0.f) {
	    goto L80;
	}
	g[i__ - 1] = (sqq - sum[i__ - 1]) / (step1 + step1);
	dzz[l - 1] = dz[l];
/* L30: */
    }
    if (number == 1) {
	goto L60;
    }
    i__1 = number;
    for (i__ = 2; i__ <= i__1; ++i__) {
	l1 = nr[i__ - 1];
	step1 = delta[i__ - 1];
	dzz[l1 - 1] = dz[l1] + step1;
	i1 = i__ - 1;
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
	    l2 = nr[j - 1];
	    step2 = delta[j - 1];
	    dzz[l2 - 1] = dz[l2] + step2;
	    fpgrdi_(ifsu, ifsv, ifbu, ifbv, &c__1, &u[1], mu, &v[1], mv, &z__[
		    1], mz, dzz, &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &
		    c__[1], nc, &sqq, fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[
		    lsu], &wrk[lsv], &wrk[lri], &wrk[lq], &wrk[lau], &wrk[
		    lav1], &wrk[lav2], &wrk[lbu], &wrk[lbv], &wrk[laa], &wrk[
		    lbb], &wrk[lcc], &wrk[lcs], &nru[1], &nrv[1]);
	    if (id0 == 0) {
/* Computing 2nd power */
		d__1 = *z0 - dzz[0];
		sqq += d__1 * d__1;
	    }
	    a[i__ + j * 6 - 7] = (sq + sqq - sum[i__ - 1] - sum[j - 1]) / (
		    step1 * step2);
	    dzz[l2 - 1] = dz[l2];
/* L40: */
	}
	dzz[l1 - 1] = dz[l1];
/* L50: */
    }
/* the optimal values g(j) are found as the solution of the system */
/* d (sq) / d (g(j)) = 0 , j=1,...,number. */
L60:
    fpsysy_(a, &number, g);
    i__1 = number;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = nr[i__ - 1];
	dz[l] += g[i__ - 1];
/* L70: */
    }
/* we determine the spline sp(u,v) according to the optimal values g(j). */
L80:
    fpgrdi_(ifsu, ifsv, ifbu, ifbv, &c__0, &u[1], mu, &v[1], mv, &z__[1], mz, 
	    &dz[1], &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc, &sq,
	     fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv], &wrk[lri]
	    , &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[lav2], &wrk[lbu], &wrk[
	    lbv], &wrk[laa], &wrk[lbb], &wrk[lcc], &wrk[lcs], &nru[1], &nrv[1]
	    );
    if (id0 == 0) {
/* Computing 2nd power */
	d__1 = *z0 - dz[1];
	*fp += d__1 * d__1;
    }
    return 0;
} /* fpopdi_ */

