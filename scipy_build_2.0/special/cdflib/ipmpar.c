/* ipmpar.f -- translated by f2c (version 20190311).
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

integer ipmpar_(integer *i__)
{
    /* Initialized data */

    static integer imach[10] = { 2,31,2147483647,2,24,-125,128,53,-1021,1024 }
	    ;

    /* System generated locals */
    integer ret_val;

/* ----------------------------------------------------------------------- */

/*     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER */
/*     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER */
/*     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ... */

/*  INTEGERS. */

/*     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM */

/*               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) ) */

/*               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1. */

/*     IPMPAR(1) = A, THE BASE. */

/*     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS. */

/*     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE. */

/*  FLOATING-POINT NUMBERS. */

/*     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING */
/*     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE */
/*     NONZERO NUMBERS ARE REPRESENTED IN THE FORM */

/*               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M) */

/*               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M, */
/*               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX. */

/*     IPMPAR(4) = B, THE BASE. */

/*  SINGLE-PRECISION */

/*     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS. */

/*     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E. */

/*     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E. */

/*  DOUBLE-PRECISION */

/*     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS. */

/*     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E. */

/*     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E. */

/* ----------------------------------------------------------------------- */

/*     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE */
/*     THE DATA STATEMENTS FOR THE COMPUTER BY REMOVING THE C FROM */
/*     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN */
/*     COLUMN 1.) */

/* ----------------------------------------------------------------------- */

/*     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY */
/*     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES). */
/*     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE */
/*     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES. */

/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Data statements .. */

/*     MACHINE CONSTANTS FOR AMDAHL MACHINES. */

/*     DATA IMACH( 1) /   2 / */
/*     DATA IMACH( 2) /  31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /  16 / */
/*     DATA IMACH( 5) /   6 / */
/*     DATA IMACH( 6) / -64 / */
/*     DATA IMACH( 7) /  63 / */
/*     DATA IMACH( 8) /  14 / */
/*     DATA IMACH( 9) / -64 / */
/*     DATA IMACH(10) /  63 / */

/*     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T */
/*     PC 7300, AND AT&T 6300. */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    24 / */
/*     DATA IMACH( 6) /  -125 / */
/*     DATA IMACH( 7) /   128 / */
/*     DATA IMACH( 8) /    53 / */
/*     DATA IMACH( 9) / -1021 / */
/*     DATA IMACH(10) /  1024 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   33 / */
/*     DATA IMACH( 3) / 8589934591 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   24 / */
/*     DATA IMACH( 6) / -256 / */
/*     DATA IMACH( 7) /  255 / */
/*     DATA IMACH( 8) /   60 / */
/*     DATA IMACH( 9) / -256 / */
/*     DATA IMACH(10) /  255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   39 / */
/*     DATA IMACH( 3) / 549755813887 / */
/*     DATA IMACH( 4) /    8 / */
/*     DATA IMACH( 5) /   13 / */
/*     DATA IMACH( 6) /  -50 / */
/*     DATA IMACH( 7) /   76 / */
/*     DATA IMACH( 8) /   26 / */
/*     DATA IMACH( 9) /  -50 / */
/*     DATA IMACH(10) /   76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS. */

/*     DATA IMACH( 1) /      2 / */
/*     DATA IMACH( 2) /     39 / */
/*     DATA IMACH( 3) / 549755813887 / */
/*     DATA IMACH( 4) /      8 / */
/*     DATA IMACH( 5) /     13 / */
/*     DATA IMACH( 6) /    -50 / */
/*     DATA IMACH( 7) /     76 / */
/*     DATA IMACH( 8) /     26 / */
/*     DATA IMACH( 9) / -32754 / */
/*     DATA IMACH(10) /  32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */
/*     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT */
/*     ARITHMETIC (NOS OPERATING SYSTEM). */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   48 / */
/*     DATA IMACH( 3) / 281474976710655 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   48 / */
/*     DATA IMACH( 6) / -974 / */
/*     DATA IMACH( 7) / 1070 / */
/*     DATA IMACH( 8) /   95 / */
/*     DATA IMACH( 9) / -926 / */
/*     DATA IMACH(10) / 1070 / */

/*     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT */
/*     ARITHMETIC (NOS/VE OPERATING SYSTEM). */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    63 / */
/*     DATA IMACH( 3) / 9223372036854775807 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    48 / */
/*     DATA IMACH( 6) / -4096 / */
/*     DATA IMACH( 7) /  4095 / */
/*     DATA IMACH( 8) /    96 / */
/*     DATA IMACH( 9) / -4096 / */
/*     DATA IMACH(10) /  4095 / */

/*     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3. */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    63 / */
/*     DATA IMACH( 3) / 9223372036854775807 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    47 / */
/*     DATA IMACH( 6) / -8189 / */
/*     DATA IMACH( 7) /  8190 / */
/*     DATA IMACH( 8) /    94 / */
/*     DATA IMACH( 9) / -8099 / */
/*     DATA IMACH(10) /  8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   15 / */
/*     DATA IMACH( 3) / 32767 / */
/*     DATA IMACH( 4) /   16 / */
/*     DATA IMACH( 5) /    6 / */
/*     DATA IMACH( 6) /  -64 / */
/*     DATA IMACH( 7) /   63 / */
/*     DATA IMACH( 8) /   14 / */
/*     DATA IMACH( 9) /  -64 / */
/*     DATA IMACH(10) /   63 / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   23 / */
/*     DATA IMACH( 3) / 8388607 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   23 / */
/*     DATA IMACH( 6) / -127 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   38 / */
/*     DATA IMACH( 9) / -127 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 */
/*     AND DPS 8/70 SERIES. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   35 / */
/*     DATA IMACH( 3) / 34359738367 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   27 / */
/*     DATA IMACH( 6) / -127 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   63 / */
/*     DATA IMACH( 9) / -127 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   15 / */
/*     DATA IMACH( 3) / 32767 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   23 / */
/*     DATA IMACH( 6) / -128 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   39 / */
/*     DATA IMACH( 9) / -128 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   15 / */
/*     DATA IMACH( 3) / 32767 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   23 / */
/*     DATA IMACH( 6) / -128 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   55 / */
/*     DATA IMACH( 9) / -128 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE HP 9000. */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    24 / */
/*     DATA IMACH( 6) /  -126 / */
/*     DATA IMACH( 7) /   128 / */
/*     DATA IMACH( 8) /    53 / */
/*     DATA IMACH( 9) / -1021 / */
/*     DATA IMACH(10) /  1024 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA */
/*     5/7/9 AND THE SEL SYSTEMS 85/86. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /   16 / */
/*     DATA IMACH( 5) /    6 / */
/*     DATA IMACH( 6) /  -64 / */
/*     DATA IMACH( 7) /   63 / */
/*     DATA IMACH( 8) /   14 / */
/*     DATA IMACH( 9) /  -64 / */
/*     DATA IMACH(10) /   63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC. */

/*      DATA imach(1)/2/ */
/*      DATA imach(2)/31/ */
/*      DATA imach(3)/2147483647/ */
/*      DATA imach(4)/2/ */
/*      DATA imach(5)/24/ */
/*      DATA imach(6)/-125/ */
/*      DATA imach(7)/128/ */
/*      DATA imach(8)/53/ */
/*      DATA imach(9)/-1021/ */
/*      DATA imach(10)/1024/ */

/*     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT */
/*     MACFORTRAN II. */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    24 / */
/*     DATA IMACH( 6) /  -125 / */
/*     DATA IMACH( 7) /   128 / */
/*     DATA IMACH( 8) /    53 / */
/*     DATA IMACH( 9) / -1021 / */
/*     DATA IMACH(10) /  1024 / */

/*     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   24 / */
/*     DATA IMACH( 6) / -127 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   56 / */
/*     DATA IMACH( 9) / -127 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   35 / */
/*     DATA IMACH( 3) / 34359738367 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   27 / */
/*     DATA IMACH( 6) / -128 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   54 / */
/*     DATA IMACH( 9) / -101 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   35 / */
/*     DATA IMACH( 3) / 34359738367 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   27 / */
/*     DATA IMACH( 6) / -128 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   62 / */
/*     DATA IMACH( 9) / -128 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   24 / */
/*     DATA IMACH( 6) / -127 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   56 / */
/*     DATA IMACH( 9) / -127 / */
/*     DATA IMACH(10) /  127 / */

/*     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000. */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    24 / */
/*     DATA IMACH( 6) /  -125 / */
/*     DATA IMACH( 7) /   128 / */
/*     DATA IMACH( 8) /    53 / */
/*     DATA IMACH( 9) / -1021 / */
/*     DATA IMACH(10) /  1024 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D */
/*     SERIES (MIPS R3000 PROCESSOR). */

/*     DATA IMACH( 1) /     2 / */
/*     DATA IMACH( 2) /    31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /     2 / */
/*     DATA IMACH( 5) /    24 / */
/*     DATA IMACH( 6) /  -125 / */
/*     DATA IMACH( 7) /   128 / */
/*     DATA IMACH( 8) /    53 / */
/*     DATA IMACH( 9) / -1021 / */
/*     DATA IMACH(10) /  1024 / */

/*     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T */
/*     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T */
/*     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300). */


/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   35 / */
/*     DATA IMACH( 3) / 34359738367 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   27 / */
/*     DATA IMACH( 6) / -128 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   60 / */
/*     DATA IMACH( 9) /-1024 / */
/*     DATA IMACH(10) / 1023 / */

/*     MACHINE CONSTANTS FOR THE VAX 11/780. */

/*     DATA IMACH( 1) /    2 / */
/*     DATA IMACH( 2) /   31 / */
/*     DATA IMACH( 3) / 2147483647 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   24 / */
/*     DATA IMACH( 6) / -127 / */
/*     DATA IMACH( 7) /  127 / */
/*     DATA IMACH( 8) /   56 / */
/*     DATA IMACH( 9) / -127 / */
/*     DATA IMACH(10) /  127 / */

    ret_val = imach[(0 + (0 + (*i__ - 1 << 2))) / 4];
    return ret_val;
} /* ipmpar_ */

