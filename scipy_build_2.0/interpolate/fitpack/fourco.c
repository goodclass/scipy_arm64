/* fourco.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fourco_(doublereal *t, integer *n, doublereal *c__, 
	doublereal *alfa, integer *m, doublereal *ress, doublereal *resc, 
	doublereal *wrk1, doublereal *wrk2, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, n4;
    static doublereal rc, rs;
    extern /* Subroutine */ int fpbfou_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *);

/*  subroutine fourco calculates the integrals */
/*                    /t(n-3) */
/*    ress(i) =      !        s(x)*sin(alfa(i)*x) dx    and */
/*              t(4)/ */
/*                    /t(n-3) */
/*    resc(i) =      !        s(x)*cos(alfa(i)*x) dx, i=1,...,m, */
/*              t(4)/ */
/*  where s(x) denotes a cubic spline which is given in its */
/*  b-spline representation. */

/*  calling sequence: */
/*     call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier) */

/*  input parameters: */
/*    t    : real array,length n, containing the knots of s(x). */
/*    n    : integer, containing the total number of knots. n>=10. */
/*    c    : real array,length n, containing the b-spline coefficients. */
/*    alfa : real array,length m, containing the parameters alfa(i). */
/*    m    : integer, specifying the number of integrals to be computed. */
/*    wrk1 : real array,length n. used as working space */
/*    wrk2 : real array,length n. used as working space */

/*  output parameters: */
/*    ress : real array,length m, containing the integrals ress(i). */
/*    resc : real array,length m, containing the integrals resc(i). */
/*    ier  : error flag: */
/*      ier=0 : normal return. */
/*      ier=10: invalid input data (see restrictions). */

/*  restrictions: */
/*    n >= 10 */
/*    t(4) < t(5) < ... < t(n-4) < t(n-3). */
/*    t(1) <= t(2) <= t(3) <= t(4). */
/*    t(n-3) <= t(n-2) <= t(n-1) <= t(n). */

/*  other subroutines required: fpbfou,fpcsin */

/*  references : */
/*    dierckx p. : calculation of fouriercoefficients of discrete */
/*                 functions using cubic splines. j. computational */
/*                 and applied mathematics 3 (1977) 207-209. */
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
    /* Parameter adjustments */
    --wrk2;
    --wrk1;
    --c__;
    --t;
    --resc;
    --ress;
    --alfa;

    /* Function Body */
    n4 = *n - 4;
/*  before starting computations a data check is made. in the input data */
/*  are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (*n < 10) {
	goto L50;
    }
    j = *n;
    for (i__ = 1; i__ <= 3; ++i__) {
	if (t[i__] > t[i__ + 1]) {
	    goto L50;
	}
	if (t[j] < t[j - 1]) {
	    goto L50;
	}
	--j;
/* L10: */
    }
    i__1 = n4;
    for (i__ = 4; i__ <= i__1; ++i__) {
	if (t[i__] >= t[i__ + 1]) {
	    goto L50;
	}
/* L20: */
    }
    *ier = 0;
/*  main loop for the different alfa(i). */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  calculate the integrals */
/*    wrk1(j) = integral(nj,4(x)*sin(alfa*x))    and */
/*    wrk2(j) = integral(nj,4(x)*cos(alfa*x)),  j=1,2,...,n-4, */
/*  where nj,4(x) denotes the normalised cubic b-spline defined on the */
/*  knots t(j),t(j+1),...,t(j+4). */
	fpbfou_(&t[1], n, &alfa[i__], &wrk1[1], &wrk2[1]);
/*  calculate the integrals ress(i) and resc(i). */
	rs = 0.f;
	rc = 0.f;
	i__2 = n4;
	for (j = 1; j <= i__2; ++j) {
	    rs += c__[j] * wrk1[j];
	    rc += c__[j] * wrk2[j];
/* L30: */
	}
	ress[i__] = rs;
	resc[i__] = rc;
/* L40: */
    }
L50:
    return 0;
} /* fourco_ */

