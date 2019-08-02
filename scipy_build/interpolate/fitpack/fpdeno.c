/* fpdeno.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpdeno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *nbind, integer *merk)
{
    static integer i__, j, k, l, point, niveau;

/*  subroutine fpdeno frees the nodes of all branches of a triply linked */
/*  tree with length < nbind by putting to zero their up field. */
/*  on exit the parameter merk points to the terminal node of the */
/*  most left branch of length nbind or takes the value 1 if there */
/*  is no such branch. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars .. */
/*  .. */
    /* Parameter adjustments */
    --right;
    --left;
    --up;

    /* Function Body */
    i__ = 1;
    niveau = 0;
L10:
    point = i__;
    i__ = left[point];
    if (i__ == 0) {
	goto L20;
    }
    ++niveau;
    goto L10;
L20:
    if (niveau == *nbind) {
	goto L70;
    }
L30:
    i__ = right[point];
    j = up[point];
    up[point] = 0;
    k = left[j];
    if (point != k) {
	goto L50;
    }
    if (i__ != 0) {
	goto L40;
    }
    --niveau;
    if (niveau == 0) {
	goto L80;
    }
    point = j;
    goto L30;
L40:
    left[j] = i__;
    goto L10;
L50:
    l = right[k];
    if (point == l) {
	goto L60;
    }
    k = l;
    goto L50;
L60:
    right[k] = i__;
    point = k;
L70:
    i__ = right[point];
    if (i__ != 0) {
	goto L10;
    }
    i__ = up[point];
    --niveau;
    if (niveau == 0) {
	goto L80;
    }
    point = i__;
    goto L70;
L80:
    k = 1;
    l = left[k];
    if (up[l] == 0) {
	return 0;
    }
L90:
    *merk = k;
    k = left[k];
    if (k != 0) {
	goto L90;
    }
    return 0;
} /* fpdeno_ */

