/* cfode.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cfode_(integer *meth, doublereal *elco, doublereal *
	tesco)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ib;
    static doublereal pc[12];
    static integer nq;
    static doublereal fnq;
    static integer nqm1, nqp1;
    static doublereal ragq, pint, xpin, fnqm1, agamq, rqfac, tsign, rq1fac;

/* lll. optimize */
/* ----------------------------------------------------------------------- */
/* cfode is called by the integrator routine to set coefficients */
/* needed there.  the coefficients for the current method, as */
/* given by the value of meth, are set for all orders and saved. */
/* the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2. */
/* (a smaller value of the maximum order is also allowed.) */
/* cfode is called once at the beginning of the problem, */
/* and is not called again unless and until meth is changed. */

/* the elco array contains the basic method coefficients. */
/* the coefficients el(i), 1 .le. i .le. nq+1, for the method of */
/* order nq are stored in elco(i,nq).  they are given by a genetrating */
/* polynomial, i.e., */
/*     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq. */
/* for the implicit adams methods, l(x) is given by */
/*     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0. */
/* for the bdf methods, l(x) is given by */
/*     l(x) = (x+1)*(x+2)* ... *(x+nq)/k, */
/* where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq). */

/* the tesco array contains test constants used for the */
/* local error test and the selection of step size and/or order. */
/* at order nq, tesco(k,nq) is used for the selection of step */
/* size at order nq - 1 if k = 1, at order nq if k = 2, and at order */
/* nq + 1 if k = 3. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    tesco -= 4;
    elco -= 14;

    /* Function Body */
    switch (*meth) {
	case 1:  goto L100;
	case 2:  goto L200;
    }

L100:
    elco[14] = 1.;
    elco[15] = 1.;
    tesco[4] = 0.;
    tesco[5] = 2.;
    tesco[7] = 1.;
    tesco[39] = 0.;
    pc[0] = 1.;
    rqfac = 1.;
    for (nq = 2; nq <= 12; ++nq) {
/* ----------------------------------------------------------------------- */
/* the pc array will contain the coefficients of the polynomial */
/*     p(x) = (x+1)*(x+2)*...*(x+nq-1). */
/* initially, p(x) = 1. */
/* ----------------------------------------------------------------------- */
	rq1fac = rqfac;
	rqfac /= (doublereal) nq;
	nqm1 = nq - 1;
	fnqm1 = (doublereal) nqm1;
	nqp1 = nq + 1;
/* form coefficients of p(x)*(x+nq-1). ---------------------------------- */
	pc[nq - 1] = 0.;
	i__1 = nqm1;
	for (ib = 1; ib <= i__1; ++ib) {
	    i__ = nqp1 - ib;
/* L110: */
	    pc[i__ - 1] = pc[i__ - 2] + fnqm1 * pc[i__ - 1];
	}
	pc[0] = fnqm1 * pc[0];
/* compute integral, -1 to 0, of p(x) and x*p(x). ----------------------- */
	pint = pc[0];
	xpin = pc[0] / 2.;
	tsign = 1.;
	i__1 = nq;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    tsign = -tsign;
	    pint += tsign * pc[i__ - 1] / (doublereal) i__;
/* L120: */
	    xpin += tsign * pc[i__ - 1] / (doublereal) (i__ + 1);
	}
/* store coefficients in elco and tesco. -------------------------------- */
	elco[nq * 13 + 1] = pint * rq1fac;
	elco[nq * 13 + 2] = 1.;
	i__1 = nq;
	for (i__ = 2; i__ <= i__1; ++i__) {
/* L130: */
	    elco[i__ + 1 + nq * 13] = rq1fac * pc[i__ - 1] / (doublereal) i__;
	}
	agamq = rqfac * xpin;
	ragq = 1. / agamq;
	tesco[nq * 3 + 2] = ragq;
	if (nq < 12) {
	    tesco[nqp1 * 3 + 1] = ragq * rqfac / (doublereal) nqp1;
	}
	tesco[nqm1 * 3 + 3] = ragq;
/* L140: */
    }
    return 0;

L200:
    pc[0] = 1.;
    rq1fac = 1.;
    for (nq = 1; nq <= 5; ++nq) {
/* ----------------------------------------------------------------------- */
/* the pc array will contain the coefficients of the polynomial */
/*     p(x) = (x+1)*(x+2)*...*(x+nq). */
/* initially, p(x) = 1. */
/* ----------------------------------------------------------------------- */
	fnq = (doublereal) nq;
	nqp1 = nq + 1;
/* form coefficients of p(x)*(x+nq). ------------------------------------ */
	pc[nqp1 - 1] = 0.;
	i__1 = nq;
	for (ib = 1; ib <= i__1; ++ib) {
	    i__ = nq + 2 - ib;
/* L210: */
	    pc[i__ - 1] = pc[i__ - 2] + fnq * pc[i__ - 1];
	}
	pc[0] = fnq * pc[0];
/* store coefficients in elco and tesco. -------------------------------- */
	i__1 = nqp1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L220: */
	    elco[i__ + nq * 13] = pc[i__ - 1] / pc[1];
	}
	elco[nq * 13 + 2] = 1.;
	tesco[nq * 3 + 1] = rq1fac;
	tesco[nq * 3 + 2] = (doublereal) nqp1 / elco[nq * 13 + 1];
	tesco[nq * 3 + 3] = (doublereal) (nq + 2) / elco[nq * 13 + 1];
	rq1fac /= fnq;
/* L230: */
    }
    return 0;
/* ----------------------- end of subroutine cfode ----------------------- */
} /* cfode_ */

