/* fprati.f -- translated by f2c (version 20190311).
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

doublereal fprati_(doublereal *p1, doublereal *f1, doublereal *p2, doublereal 
	*f2, doublereal *p3, doublereal *f3)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal p, h1, h2, h3;

/*  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati */
/*  gives the value of p such that the rational interpolating function */
/*  of the form r(p) = (u*p+v)/(p+w) equals zero at p. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..local scalars.. */
/*  .. */
    if (*p3 > 0.f) {
	goto L10;
    }
/*  value of p in case p3 = infinity. */
    p = (*p1 * (*f1 - *f3) * *f2 - *p2 * (*f2 - *f3) * *f1) / ((*f1 - *f2) * *
	    f3);
    goto L20;
/*  value of p in case p3 ^= infinity. */
L10:
    h1 = *f1 * (*f2 - *f3);
    h2 = *f2 * (*f3 - *f1);
    h3 = *f3 * (*f1 - *f2);
    p = -(*p1 * *p2 * h3 + *p2 * *p3 * h1 + *p3 * *p1 * h2) / (*p1 * h1 + *p2 
	    * h2 + *p3 * h3);
/*  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0. */
L20:
    if (*f2 < 0.f) {
	goto L30;
    }
    *p1 = *p2;
    *f1 = *f2;
    goto L40;
L30:
    *p3 = *p2;
    *f3 = *f2;
L40:
    ret_val = p;
    return ret_val;
} /* fprati_ */

