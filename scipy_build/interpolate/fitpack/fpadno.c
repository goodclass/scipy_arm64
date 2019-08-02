/* fpadno.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpadno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *info, integer *count, integer *merk, integer 
	*jbind, integer *n1, integer *ier)
{
    static integer k;
    static logical bool;
    static integer point, niveau;
    extern /* Subroutine */ int fpfrno_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);

/*  subroutine fpadno adds a branch of length n1 to the triply linked */
/*  tree,the information of which is kept in the arrays up,left,right */
/*  and info. the information field of the nodes of this new branch is */
/*  given in the array jbind. in linking the new branch fpadno takes */
/*  account of the property of the tree that */
/*    info(k) < info(right(k)) ; info(k) < info(left(k)) */
/*  if necessary the subroutine calls subroutine fpfrno to collect the */
/*  free nodes of the tree. if no computer words are available at that */
/*  moment, the error parameter ier is set to 1. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..subroutine references.. */
/*    fpfrno */
/*  .. */
    /* Parameter adjustments */
    --info;
    --right;
    --left;
    --up;
    --jbind;

    /* Function Body */
    point = 1;
    niveau = 1;
L10:
    k = left[point];
    bool = TRUE_;
L20:
    if (k == 0) {
	goto L50;
    }
    if (info[k] - jbind[niveau] < 0) {
	goto L30;
    }
    if (info[k] - jbind[niveau] == 0) {
	goto L40;
    }
    goto L50;
L30:
    point = k;
    k = right[point];
    bool = FALSE_;
    goto L20;
L40:
    point = k;
    ++niveau;
    goto L10;
L50:
    if (niveau > *n1) {
	goto L90;
    }
    ++(*count);
    if (*count <= *maxtr) {
	goto L60;
    }
    fpfrno_(maxtr, &up[1], &left[1], &right[1], &info[1], &point, merk, n1, 
	    count, ier);
    if (*ier != 0) {
	goto L100;
    }
L60:
    info[*count] = jbind[niveau];
    left[*count] = 0;
    right[*count] = k;
    if (bool) {
	goto L70;
    }
    bool = TRUE_;
    right[point] = *count;
    up[*count] = up[point];
    goto L80;
L70:
    up[*count] = point;
    left[point] = *count;
L80:
    point = *count;
    ++niveau;
    k = 0;
    goto L50;
L90:
    *ier = 0;
L100:
    return 0;
} /* fpadno_ */

