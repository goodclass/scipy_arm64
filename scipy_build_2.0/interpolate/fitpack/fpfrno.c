/* fpfrno.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpfrno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *info, integer *point, integer *merk, integer 
	*n1, integer *count, integer *ier)
{
    static integer i__, j, k, l, n, niveau;

/*  subroutine fpfrno collects the free nodes (up field zero) of the */
/*  triply linked tree the information of which is kept in the arrays */
/*  up,left,right and info. the maximal length of the branches of the */
/*  tree is given by n1. if no free nodes are found, the error flag */
/*  ier is set to 1. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars */
/*  .. */
    /* Parameter adjustments */
    --info;
    --right;
    --left;
    --up;

    /* Function Body */
    *ier = 1;
    if (*n1 == 2) {
	goto L140;
    }
    niveau = 1;
    *count = 2;
L10:
    j = 0;
    i__ = 1;
L20:
    if (j == niveau) {
	goto L30;
    }
    k = 0;
    l = left[i__];
    if (l == 0) {
	goto L110;
    }
    i__ = l;
    ++j;
    goto L20;
L30:
    if (i__ < *count) {
	goto L110;
    }
    if (i__ == *count) {
	goto L100;
    }
    goto L40;
L40:
    if (up[*count] == 0) {
	goto L50;
    }
    ++(*count);
    goto L30;
L50:
    up[*count] = up[i__];
    left[*count] = left[i__];
    right[*count] = right[i__];
    info[*count] = info[i__];
    if (*merk == i__) {
	*merk = *count;
    }
    if (*point == i__) {
	*point = *count;
    }
    if (k == 0) {
	goto L60;
    }
    right[k] = *count;
    goto L70;
L60:
    n = up[i__];
    left[n] = *count;
L70:
    l = left[i__];
L80:
    if (l == 0) {
	goto L90;
    }
    up[l] = *count;
    l = right[l];
    goto L80;
L90:
    up[i__] = 0;
    i__ = *count;
L100:
    ++(*count);
L110:
    l = right[i__];
    k = i__;
    if (l == 0) {
	goto L120;
    }
    i__ = l;
    goto L20;
L120:
    l = up[i__];
    --j;
    if (j == 0) {
	goto L130;
    }
    i__ = l;
    goto L110;
L130:
    ++niveau;
    if (niveau <= *n1) {
	goto L10;
    }
    if (*count > *maxtr) {
	goto L140;
    }
    *ier = 0;
L140:
    return 0;
} /* fpfrno_ */

