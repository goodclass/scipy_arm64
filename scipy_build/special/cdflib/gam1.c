/* gam1.f -- translated by f2c (version 20190311).
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

doublereal gam1_(doublereal *a)
{
    /* Initialized data */

    static doublereal p[7] = { .577215664901533,-.409078193005776,
	    -.230975380857675,.0597275330452234,.0076696818164949,
	    -.00514889771323592,5.89597428611429e-4 };
    static doublereal q[5] = { 1.,.427569613095214,.158451672430138,
	    .0261132021441447,.00423244297896961 };
    static doublereal r__[9] = { -.422784335098468,-.771330383816272,
	    -.244757765222226,.118378989872749,9.30357293360349e-4,
	    -.0118290993445146,.00223047661158249,2.66505979058923e-4,
	    -1.32674909766242e-4 };
    static doublereal s1 = .273076135303957;
    static doublereal s2 = .0559398236957378;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal d__, t, w, bot, top;

/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5 */
/*     ------------------------------------------------------------------ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Data statements .. */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     .. */
/*     .. Executable Statements .. */
/*     ------------------- */
    t = *a;
    d__ = *a - .5;
    if (d__ > 0.) {
	t = d__ - .5;
    }
    if (t < 0.) {
	goto L40;
    }
    if (t == 0.) {
	goto L10;
    }
    goto L20;

L10:
    ret_val = 0.;
    return ret_val;

L20:
    top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]
	    ) * t + p[0];
    bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
    w = top / bot;
    if (d__ > 0.) {
	goto L30;
    }
    ret_val = *a * w;
    return ret_val;
L30:
    ret_val = t / *a * (w - .5 - .5);
    return ret_val;

L40:
    top = (((((((r__[8] * t + r__[7]) * t + r__[6]) * t + r__[5]) * t + r__[4]
	    ) * t + r__[3]) * t + r__[2]) * t + r__[1]) * t + r__[0];
    bot = (s2 * t + s1) * t + 1.;
    w = top / bot;
    if (d__ > 0.) {
	goto L50;
    }
    ret_val = *a * (w + .5 + .5);
    return ret_val;
L50:
    ret_val = t * w / *a;
    return ret_val;
} /* gam1_ */

