/* exparg.f -- translated by f2c (version 20190311).
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

/* Table of constant values */

static integer c__4 = 4;
static integer c__9 = 9;
static integer c__10 = 10;

doublereal exparg_(integer *l)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer b, m;
    static doublereal lnb;
    extern integer ipmpar_(integer *);

/* -------------------------------------------------------------------- */
/*     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH */
/*     EXP(W) CAN BE COMPUTED. */

/*     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR */
/*     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO. */

/*     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED. */
/* -------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    b = ipmpar_(&c__4);
    if (b != 2) {
	goto L10;
    }
    lnb = .69314718055995;
    goto L40;
L10:
    if (b != 8) {
	goto L20;
    }
    lnb = 2.0794415416798;
    goto L40;
L20:
    if (b != 16) {
	goto L30;
    }
    lnb = 2.7725887222398;
    goto L40;
L30:
    lnb = log((doublereal) b);

L40:
    if (*l == 0) {
	goto L50;
    }
    m = ipmpar_(&c__9) - 1;
    ret_val = m * lnb * .99999;
    return ret_val;
L50:
    m = ipmpar_(&c__10);
    ret_val = m * lnb * .99999;
    return ret_val;
} /* exparg_ */

