/* gsumln.f -- translated by f2c (version 20190311).
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

doublereal gsumln_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal x;
    extern doublereal gamln1_(doublereal *), alnrel_(doublereal *);

/* ----------------------------------------------------------------------- */
/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B)) */
/*          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2 */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
    x = *a + *b - 2.;
    if (x > .25) {
	goto L10;
    }
    d__1 = x + 1.;
    ret_val = gamln1_(&d__1);
    return ret_val;
L10:
    if (x > 1.25) {
	goto L20;
    }
    ret_val = gamln1_(&x) + alnrel_(&x);
    return ret_val;
L20:
    d__1 = x - 1.;
    ret_val = gamln1_(&d__1) + log(x * (x + 1.));
    return ret_val;
} /* gsumln_ */

