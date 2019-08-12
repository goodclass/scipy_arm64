/* gamln1.f -- translated by f2c (version 20190311).
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

doublereal gamln1_(doublereal *a)
{
    /* Initialized data */

    static doublereal p0 = .577215664901533;
    static doublereal q3 = 1.56875193295039;
    static doublereal q4 = .361951990101499;
    static doublereal q5 = .0325038868253937;
    static doublereal q6 = 6.67465618796164e-4;
    static doublereal r0 = .422784335098467;
    static doublereal r1 = .848044614534529;
    static doublereal r2 = .565221050691933;
    static doublereal r3 = .156513060486551;
    static doublereal r4 = .017050248402265;
    static doublereal r5 = 4.97958207639485e-4;
    static doublereal p1 = .844203922187225;
    static doublereal s1 = 1.24313399877507;
    static doublereal s2 = .548042109832463;
    static doublereal s3 = .10155218743983;
    static doublereal s4 = .00713309612391;
    static doublereal s5 = 1.16165475989616e-4;
    static doublereal p2 = -.168860593646662;
    static doublereal p3 = -.780427615533591;
    static doublereal p4 = -.402055799310489;
    static doublereal p5 = -.0673562214325671;
    static doublereal p6 = -.00271935708322958;
    static doublereal q1 = 2.88743195473681;
    static doublereal q2 = 3.12755088914843;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal w, x;

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25 */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Data statements .. */
/* ---------------------- */
/*     .. */
/*     .. Executable Statements .. */
/* ---------------------- */
    if (*a >= .6) {
	goto L10;
    }
    w = ((((((p6 * *a + p5) * *a + p4) * *a + p3) * *a + p2) * *a + p1) * *a 
	    + p0) / ((((((q6 * *a + q5) * *a + q4) * *a + q3) * *a + q2) * *a 
	    + q1) * *a + 1.);
    ret_val = -(*a) * w;
    return ret_val;

L10:
    x = *a - .5 - .5;
    w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) / (((((s5 * 
	    x + s4) * x + s3) * x + s2) * x + s1) * x + 1.);
    ret_val = x * w;
    return ret_val;
} /* gamln1_ */

