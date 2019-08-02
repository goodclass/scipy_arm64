/* dqwgts.f -- translated by f2c (version 20190311).
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

doublereal dqwgts_(doublereal *x, doublereal *a, doublereal *b, doublereal *
	alfa, doublereal *beta, integer *integr)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    static doublereal xma, bmx;

/* ***begin prologue  dqwgts */
/* ***refer to dqk15w */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  weight function, algebraico-logarithmic */
/*             end-point singularities */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this function subprogram is used together with the */
/*            routine dqaws and defines the weight function. */
/* ***end prologue  dqwgts */

/* ***first executable statement  dqwgts */
    xma = *x - *a;
    bmx = *b - *x;
    ret_val = pow_dd(&xma, alfa) * pow_dd(&bmx, beta);
    switch (*integr) {
	case 1:  goto L40;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L30;
    }
L10:
    ret_val *= log(xma);
    goto L40;
L20:
    ret_val *= log(bmx);
    goto L40;
L30:
    ret_val = ret_val * log(xma) * log(bmx);
L40:
    return ret_val;
} /* dqwgts_ */

