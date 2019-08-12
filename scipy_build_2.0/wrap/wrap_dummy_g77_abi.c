/* wrap_dummy_g77_abi.f -- translated by f2c (version 20190311).
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

/* Complex */ VOID wcdotc_(complex * ret_val, integer *n, complex *cx, 
	integer *incx, complex *cy, integer *incy)
{
    /* System generated locals */
    complex q__1;

    /* Local variables */
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    cdotc_(&q__1, n, &cx[1], incx, &cy[1], incy);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
} /* wcdotc_ */

/* Complex */ VOID wcdotu_(complex * ret_val, integer *n, complex *cx, 
	integer *incx, complex *cy, integer *incy)
{
    /* System generated locals */
    complex q__1;

    /* Local variables */
    extern /* Complex */ VOID cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    cdotu_(&q__1, n, &cx[1], incx, &cy[1], incy);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
} /* wcdotu_ */

/* Double Complex */ VOID wzdotc_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    doublecomplex z__1;

    /* Local variables */
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    zdotc_(&z__1, n, &cx[1], incx, &cy[1], incy);
     ret_val->r = z__1.r,  ret_val->i = z__1.i;
} /* wzdotc_ */

/* Double Complex */ VOID wzdotu_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    doublecomplex z__1;

    /* Local variables */
    extern /* Double Complex */ VOID zdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    zdotu_(&z__1, n, &cx[1], incx, &cy[1], incy);
     ret_val->r = z__1.r,  ret_val->i = z__1.i;
} /* wzdotu_ */

/* Complex */ VOID wcladiv_(complex * ret_val, complex *x, complex *y)
{
    /* System generated locals */
    complex q__1;

    /* Local variables */
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);

    cladiv_(&q__1, x, y);
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
} /* wcladiv_ */

/* Double Complex */ VOID wzladiv_(doublecomplex * ret_val, doublecomplex *x, 
	doublecomplex *y)
{
    /* System generated locals */
    doublecomplex z__1;

    /* Local variables */
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);

    zladiv_(&z__1, x, y);
     ret_val->r = z__1.r,  ret_val->i = z__1.i;
} /* wzladiv_ */

