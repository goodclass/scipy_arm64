/* vmnorm.f -- translated by f2c (version 20190311).
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

doublereal vmnorm_(integer *n, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal vm;

/* lll. optimize */
/* ----------------------------------------------------------------------- */
/* this function routine computes the weighted max-norm */
/* of the vector of length n contained in the array v, with weights */
/* contained in the array w of length n.. */
/*   vmnorm = max(i=1,...,n) abs(v(i))*w(i) */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    vm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing MAX */
	d__2 = vm, d__3 = (d__1 = v[i__], abs(d__1)) * w[i__];
	vm = max(d__2,d__3);
    }
    ret_val = vm;
    return ret_val;
/* ----------------------- end of function vmnorm ------------------------ */
} /* vmnorm_ */

