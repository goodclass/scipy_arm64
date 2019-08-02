/* fpknot.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpknot_(doublereal *x, integer *m, doublereal *t, 
	integer *n, doublereal *fpint, integer *nrdata, integer *nrint, 
	integer *nest, integer *istart)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static doublereal am, an;
    static integer jj, jk, nrx, next, ihalf;
    static doublereal fpmax;
    static integer maxpt, jbegin, maxbeg, number, jpoint;

/*  subroutine fpknot locates an additional knot for a spline of degree */
/*  k and adjusts the corresponding parameters,i.e. */
/*    t     : the position of the knots. */
/*    n     : the number of knots. */
/*    nrint : the number of knotintervals. */
/*    fpint : the sum of squares of residual right hand sides */
/*            for each knot interval. */
/*    nrdata: the number of data points inside each knot interval. */
/*  istart indicates that the smallest data point at which the new knot */
/*  may be added is x(istart+1) */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --x;
    --nrdata;
    --fpint;
    --t;

    /* Function Body */
    k = (*n - *nrint - 1) / 2;
/*  search for knot interval t(number+k) <= x <= t(number+k+1) where */
/*  fpint(number) is maximal on the condition that nrdata(number) */
/*  not equals zero. */
    fpmax = 0.f;
    jbegin = *istart;
    i__1 = *nrint;
    for (j = 1; j <= i__1; ++j) {
	jpoint = nrdata[j];
	if (fpmax >= fpint[j] || jpoint == 0) {
	    goto L10;
	}
	fpmax = fpint[j];
	number = j;
	maxpt = jpoint;
	maxbeg = jbegin;
L10:
	jbegin = jbegin + jpoint + 1;
/* L20: */
    }
/*  let coincide the new knot t(number+k+1) with a data point x(nrx) */
/*  inside the old knot interval t(number+k) <= x <= t(number+k+1). */
    ihalf = maxpt / 2 + 1;
    nrx = maxbeg + ihalf;
    next = number + 1;
    if (next > *nrint) {
	goto L40;
    }
/*  adjust the different parameters. */
    i__1 = *nrint;
    for (j = next; j <= i__1; ++j) {
	jj = next + *nrint - j;
	fpint[jj + 1] = fpint[jj];
	nrdata[jj + 1] = nrdata[jj];
	jk = jj + k;
	t[jk + 1] = t[jk];
/* L30: */
    }
L40:
    nrdata[number] = ihalf - 1;
    nrdata[next] = maxpt - ihalf;
    am = (doublereal) maxpt;
    an = (doublereal) nrdata[number];
    fpint[number] = fpmax * an / am;
    an = (doublereal) nrdata[next];
    fpint[next] = fpmax * an / am;
    jk = next + k;
    t[jk] = x[nrx];
    ++(*n);
    ++(*nrint);
    return 0;
} /* fpknot_ */

