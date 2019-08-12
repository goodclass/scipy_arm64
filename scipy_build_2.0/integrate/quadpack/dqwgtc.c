/* dqwgtc.f -- translated by f2c (version 20190311).
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

doublereal dqwgtc_(doublereal *x, doublereal *c__, doublereal *p2, doublereal 
	*p3, doublereal *p4, integer *kp)
{
    /* System generated locals */
    doublereal ret_val;

/* ***begin prologue  dqwgtc */
/* ***refer to dqk15w */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  weight function, cauchy principal value */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this function subprogram is used together with the */
/*            routine qawc and defines the weight function. */
/* ***end prologue  dqwgtc */

/* ***first executable statement  dqwgtc */
    ret_val = 1. / (*x - *c__);
    return ret_val;
} /* dqwgtc_ */

