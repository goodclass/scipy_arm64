/* dqelg.f -- translated by f2c (version 20190311).
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

static integer c__4 = 4;
static integer c__2 = 2;

/* Subroutine */ int dqelg_(integer *n, doublereal *epstab, doublereal *
	result, doublereal *abserr, doublereal *res3la, integer *nres)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal e0, e1, e2, e3;
    static integer k1, k2, k3, ib, ie;
    static doublereal ss;
    static integer ib2;
    static doublereal res;
    static integer num;
    static doublereal err1, err2, err3, tol1, tol2, tol3;
    static integer indx;
    static doublereal e1abs, oflow, error;
    extern doublereal d1mach_(integer *);
    static doublereal delta1, delta2, delta3, epmach, epsinf;
    static integer newelm, limexp;

/* ***begin prologue  dqelg */
/* ***refer to  dqagie,dqagoe,dqagpe,dqagse */
/* ***routines called  d1mach */
/* ***revision date  830518   (yymmdd) */
/* ***keywords  epsilon algorithm, convergence acceleration, */
/*             extrapolation */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
/* ***purpose  the routine determines the limit of a given sequence of */
/*            approximations, by means of the epsilon algorithm of */
/*            p.wynn. an estimate of the absolute error is also given. */
/*            the condensed epsilon table is computed. only those */
/*            elements needed for the computation of the next diagonal */
/*            are preserved. */
/* ***description */

/*           epsilon algorithm */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters */
/*              n      - integer */
/*                       epstab(n) contains the new element in the */
/*                       first column of the epsilon table. */

/*              epstab - double precision */
/*                       vector of dimension 52 containing the elements */
/*                       of the two lower diagonals of the triangular */
/*                       epsilon table. the elements are numbered */
/*                       starting at the right-hand corner of the */
/*                       triangle. */

/*              result - double precision */
/*                       resulting approximation to the integral */

/*              abserr - double precision */
/*                       estimate of the absolute error computed from */
/*                       result and the 3 previous results */

/*              res3la - double precision */
/*                       vector of dimension 3 containing the last 3 */
/*                       results */

/*              nres   - integer */
/*                       number of calls to the routine */
/*                       (should be zero at first call) */

/* ***end prologue  dqelg */


/*           list of major variables */
/*           ----------------------- */

/*           e0     - the 4 elements on which the computation of a new */
/*           e1       element in the epsilon table is based */
/*           e2 */
/*           e3                 e0 */
/*                        e3    e1    new */
/*                              e2 */
/*           newelm - number of elements to be computed in the new */
/*                    diagonal */
/*           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2) */
/*           result - the element in the new diagonal with least value */
/*                    of error */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           oflow is the largest positive magnitude. */
/*           limexp is the maximum number of elements the epsilon */
/*           table can contain. if this number is reached, the upper */
/*           diagonal of the epsilon table is deleted. */

/* ***first executable statement  dqelg */
    /* Parameter adjustments */
    --res3la;
    --epstab;

    /* Function Body */
    epmach = d1mach_(&c__4);
    oflow = d1mach_(&c__2);
    ++(*nres);
    *abserr = oflow;
    *result = epstab[*n];
    if (*n < 3) {
	goto L100;
    }
    limexp = 50;
    epstab[*n + 2] = epstab[*n];
    newelm = (*n - 1) / 2;
    epstab[*n] = oflow;
    num = *n;
    k1 = *n;
    i__1 = newelm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k2 = k1 - 1;
	k3 = k1 - 2;
	res = epstab[k1 + 2];
	e0 = epstab[k3];
	e1 = epstab[k2];
	e2 = res;
	e1abs = abs(e1);
	delta2 = e2 - e1;
	err2 = abs(delta2);
/* Computing MAX */
	d__1 = abs(e2);
	tol2 = max(d__1,e1abs) * epmach;
	delta3 = e1 - e0;
	err3 = abs(delta3);
/* Computing MAX */
	d__1 = e1abs, d__2 = abs(e0);
	tol3 = max(d__1,d__2) * epmach;
	if (err2 > tol2 || err3 > tol3) {
	    goto L10;
	}

/*           if e0, e1 and e2 are equal to within machine */
/*           accuracy, convergence is assumed. */
/*           result = e2 */
/*           abserr = abs(e1-e0)+abs(e2-e1) */

	*result = res;
	*abserr = err2 + err3;
/* ***jump out of do-loop */
	goto L100;
L10:
	e3 = epstab[k1];
	epstab[k1] = e1;
	delta1 = e1 - e3;
	err1 = abs(delta1);
/* Computing MAX */
	d__1 = e1abs, d__2 = abs(e3);
	tol1 = max(d__1,d__2) * epmach;

/*           if two elements are very close to each other, omit */
/*           a part of the table by adjusting the value of n */

	if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
	    goto L20;
	}
	ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
	epsinf = (d__1 = ss * e1, abs(d__1));

/*           test to detect irregular behaviour in the table, and */
/*           eventually omit a part of the table adjusting the value */
/*           of n. */

	if (epsinf > 1e-4) {
	    goto L30;
	}
L20:
	*n = i__ + i__ - 1;
/* ***jump out of do-loop */
	goto L50;

/*           compute a new element and eventually adjust */
/*           the value of result. */

L30:
	res = e1 + 1. / ss;
	epstab[k1] = res;
	k1 += -2;
	error = err2 + (d__1 = res - e2, abs(d__1)) + err3;
	if (error > *abserr) {
	    goto L40;
	}
	*abserr = error;
	*result = res;
L40:
	;
    }

/*           shift the table. */

L50:
    if (*n == limexp) {
	*n = (limexp / 2 << 1) - 1;
    }
    ib = 1;
    if (num / 2 << 1 == num) {
	ib = 2;
    }
    ie = newelm + 1;
    i__1 = ie;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib2 = ib + 2;
	epstab[ib] = epstab[ib2];
	ib = ib2;
/* L60: */
    }
    if (num == *n) {
	goto L80;
    }
    indx = num - *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epstab[i__] = epstab[indx];
	++indx;
/* L70: */
    }
L80:
    if (*nres >= 4) {
	goto L90;
    }
    res3la[*nres] = *result;
    *abserr = oflow;
    goto L100;

/*           compute error estimate */

L90:
    *abserr = (d__1 = *result - res3la[3], abs(d__1)) + (d__2 = *result - 
	    res3la[2], abs(d__2)) + (d__3 = *result - res3la[1], abs(d__3));
    res3la[1] = res3la[2];
    res3la[2] = res3la[3];
    res3la[3] = *result;
L100:
/* Computing MAX */
    d__1 = *abserr, d__2 = epmach * 5. * abs(*result);
    *abserr = max(d__1,d__2);
    return 0;
} /* dqelg_ */

