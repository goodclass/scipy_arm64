/* spmpar.f -- translated by f2c (version 20190311).
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
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;

doublereal spmpar_(integer *i__)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal b;
    static integer m;
    static doublereal w, z__, bm1, one;
    static integer emin, emax;
    static doublereal binv;
    static integer ibeta;
    extern integer ipmpar_(integer *);

/* ----------------------------------------------------------------------- */

/*     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR */
/*     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT */
/*     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE */
/*     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND */
/*     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN */

/*        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION, */

/*        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE, */

/*        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE. */

/* ----------------------------------------------------------------------- */
/*     WRITTEN BY */
/*        ALFRED H. MORRIS, JR. */
/*        NAVAL SURFACE WARFARE CENTER */
/*        DAHLGREN VIRGINIA */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE */
/*     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS */
/*     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION */
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

    if (*i__ > 1) {
	goto L10;
    }
    b = (doublereal) ipmpar_(&c__4);
    m = ipmpar_(&c__8);
    i__1 = 1 - m;
    ret_val = pow_di(&b, &i__1);
    return ret_val;

L10:
    if (*i__ > 2) {
	goto L20;
    }
    b = (doublereal) ipmpar_(&c__4);
    emin = ipmpar_(&c__9);
    one = 1.;
    binv = one / b;
    i__1 = emin + 2;
    w = pow_di(&b, &i__1);
    ret_val = w * binv * binv * binv;
    return ret_val;

L20:
    ibeta = ipmpar_(&c__4);
    m = ipmpar_(&c__8);
    emax = ipmpar_(&c__10);

    b = (doublereal) ibeta;
    bm1 = (doublereal) (ibeta - 1);
    one = 1.;
    i__1 = m - 1;
    z__ = pow_di(&b, &i__1);
    w = ((z__ - one) * b + bm1) / (b * z__);

    i__1 = emax - 2;
    z__ = pow_di(&b, &i__1);
    ret_val = w * z__ * b * b;
    return ret_val;
} /* spmpar_ */

