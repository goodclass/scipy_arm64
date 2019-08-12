/* fpseno.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpseno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *info, integer *merk, integer *ibind, integer 
	*nbind)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;

/*  subroutine fpseno fetches a branch of a triply linked tree the */
/*  information of which is kept in the arrays up,left,right and info. */
/*  the branch has a specified length nbind and is determined by the */
/*  parameter merk which points to its terminal node. the information */
/*  field of the nodes of this branch is stored in the array ibind. on */
/*  exit merk points to a new branch of length nbind or takes the value */
/*  1 if no such branch was found. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..scalar arguments.. */
/*  .. */
    /* Parameter adjustments */
    --info;
    --right;
    --left;
    --up;
    --ibind;

    /* Function Body */
    k = *merk;
    j = *nbind;
    i__1 = *nbind;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ibind[j] = info[k];
	k = up[k];
	--j;
/* L10: */
    }
L20:
    k = right[*merk];
    if (k != 0) {
	goto L30;
    }
    *merk = up[*merk];
    if (*merk <= 1) {
	goto L40;
    }
    goto L20;
L30:
    *merk = k;
    k = left[*merk];
    if (k != 0) {
	goto L30;
    }
L40:
    return 0;
} /* fpseno_ */

