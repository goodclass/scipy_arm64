/* dqc25c.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dqc25c_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*c__, doublereal *result, doublereal *abserr, integer *krul, integer *
	neval)
{
    /* Initialized data */

    static doublereal x[11] = { .991444861373810411144557526928563,
	    .965925826289068286749743199728897,
	    .923879532511286756128183189396788,
	    .866025403784438646763723170752936,
	    .793353340291235164579776961501299,
	    .707106781186547524400844362104849,
	    .608761429008720639416097542898164,.5,
	    .382683432365089771728459984030399,
	    .258819045102520762348898837624048,
	    .130526192220051591548406227895489 };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal u, p2, p3, p4, cc;
    static integer kp;
    static doublereal ak22, fval[25], res12, res24;
    static integer isym;
    static doublereal amom0, amom1, amom2, cheb12[13], cheb24[25], hlgth, 
	    centr;
    extern /* Subroutine */ int dqk15w_(D_fp, D_fp, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    dqcheb_(doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal resabs, resasc;
    extern doublereal dqwgtc_();

/* ***begin prologue  dqc25c */
/* ***date written   810101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a2,j4 */
/* ***keywords  25-point clenshaw-curtis integration */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f*w over (a,b) with */
/*            error estimate, where w(x) = 1/(x-c) */
/* ***description */

/*        integration rules for the computation of cauchy */
/*        principal value integrals */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*           f      - double precision */
/*                    function subprogram defining the integrand function */
/*                    f(x). the actual name for f needs to be declared */
/*                    e x t e r n a l  in the driver program. */

/*           a      - double precision */
/*                    left end point of the integration interval */

/*           b      - double precision */
/*                    right end point of the integration interval, b.gt.a */

/*           c      - double precision */
/*                    parameter in the weight function */

/*           result - double precision */
/*                    approximation to the integral */
/*                    result is computed by using a generalized */
/*                    clenshaw-curtis method if c lies within ten percent */
/*                    of the integration interval. in the other case the */
/*                    15-point kronrod rule obtained by optimal addition */
/*                    of abscissae to the 7-point gauss rule, is applied. */

/*           abserr - double precision */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed abs(i-result) */

/*           krul   - integer */
/*                    key which is decreased by 1 if the 15-point */
/*                    gauss-kronrod scheme has been used */

/*           neval  - integer */
/*                    number of integrand evaluations */

/* ....................................................................... */
/* ***references  (none) */
/* ***routines called  dqcheb,dqk15w,dqwgtc */
/* ***end prologue  dqc25c */




/*           the vector x contains the values cos(k*pi/24), */
/*           k = 1, ..., 11, to be used for the chebyshev series */
/*           expansion of f */


/*           list of major variables */
/*           ---------------------- */
/*           fval   - value of the function f at the points */
/*                    cos(k*pi/24),  k = 0, ..., 24 */
/*           cheb12 - chebyshev series expansion coefficients, */
/*                    for the function f, of degree 12 */
/*           cheb24 - chebyshev series expansion coefficients, */
/*                    for the function f, of degree 24 */
/*           res12  - approximation to the integral corresponding */
/*                    to the use of cheb12 */
/*           res24  - approximation to the integral corresponding */
/*                    to the use of cheb24 */
/*           dqwgtc - external function subprogram defining */
/*                    the weight function */
/*           hlgth  - half-length of the interval */
/*           centr  - mid point of the interval */


/*           check the position of c. */

/* ***first executable statement  dqc25c */
    cc = (2. * *c__ - *b - *a) / (*b - *a);
    if (abs(cc) < 1.1) {
	goto L10;
    }

/*           apply the 15-point gauss-kronrod scheme. */

    --(*krul);
    dqk15w_((D_fp)f, (D_fp)dqwgtc_, c__, &p2, &p3, &p4, &kp, a, b, result, 
	    abserr, &resabs, &resasc);
    *neval = 15;
    if (resasc == *abserr) {
	++(*krul);
    }
    goto L50;

/*           use the generalized clenshaw-curtis method. */

L10:
    hlgth = (*b - *a) * .5;
    centr = (*b + *a) * .5;
    *neval = 25;
    d__1 = hlgth + centr;
    fval[0] = (*f)(&d__1) * .5;
    fval[12] = (*f)(&centr);
    d__1 = centr - hlgth;
    fval[24] = (*f)(&d__1) * .5;
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	fval[i__ - 1] = (*f)(&d__1);
	d__1 = centr - u;
	fval[isym - 1] = (*f)(&d__1);
/* L20: */
    }

/*           compute the chebyshev series expansion. */

    dqcheb_(x, fval, cheb12, cheb24);

/*           the modified chebyshev moments are computed by forward */
/*           recursion, using amom0 and amom1 as starting values. */

    amom0 = log((d__1 = (1. - cc) / (cc + 1.), abs(d__1)));
    amom1 = cc * amom0 + 2.;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 3; k <= 13; ++k) {
	amom2 = cc * 2. * amom1 - amom0;
	ak22 = (doublereal) ((k - 2) * (k - 2));
	if (k / 2 << 1 == k) {
	    amom2 -= 4. / (ak22 - 1.);
	}
	res12 += cheb12[k - 1] * amom2;
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
/* L30: */
    }
    for (k = 14; k <= 25; ++k) {
	amom2 = cc * 2. * amom1 - amom0;
	ak22 = (doublereal) ((k - 2) * (k - 2));
	if (k / 2 << 1 == k) {
	    amom2 -= 4. / (ak22 - 1.);
	}
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
/* L40: */
    }
    *result = res24;
    *abserr = (d__1 = res24 - res12, abs(d__1));
L50:
    return 0;
} /* dqc25c_ */

