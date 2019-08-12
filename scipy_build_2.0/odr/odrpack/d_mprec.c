/* d_mprec.f -- translated by f2c (version 20190311).
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

/* DMPREC */
doublereal dmprec_(void)
{
    /* Initialized data */

    static doublereal b = 2.;
    static integer td = 53;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

/* ***BEGIN PROLOGUE  DPREC */
/* ***REFER TO  DODR,DODRC */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   860529   (YYMMDD) */
/* ***REVISION DATE  920304   (YYMMDD) */
/* ***PURPOSE  DETERMINE MACHINE PRECISION FOR TARGET MACHINE AND COMPILER */
/*            ASSUMING FLOATING-POINT NUMBERS ARE REPRESENTED IN THE */
/*            T-DIGIT, BASE-B FORM */
/*                  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */
/*            WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T, AND */
/*                  0 .LT. X(1). */
/*            TO ALTER THIS FUNCTION FOR A PARTICULAR TARGET MACHINE, */
/*            EITHER */
/*                  ACTIVATE THE DESIRED SET OF DATA STATEMENTS BY */
/*                  REMOVING THE C FROM COLUMN 1 */
/*            OR */
/*                  SET B, TD AND TS USING I1MACH BY ACTIVATING */
/*                  THE DECLARATION STATEMENTS FOR I1MACH */
/*                  AND THE STATEMENTS PRECEDING THE FIRST */
/*                  EXECUTABLE STATEMENT BELOW. */
/* ***END PROLOGUE  DPREC */
/* ...LOCAL SCALARS */
/* ...EXTERNAL FUNCTIONS */
/*     INTEGER */
/*    +   I1MACH */
/*     EXTERNAL */
/*    +   I1MACH */
/* ...VARIABLE DEFINITIONS (ALPHABETICALLY) */
/*     DOUBLE PRECISION B */
/*        THE BASE OF THE TARGET MACHINE. */
/*        (MAY BE DEFINED USING I1MACH(10).) */
/*     INTEGER TD */
/*        THE NUMBER OF BASE-B DIGITS IN DOUBLE PRECISION. */
/*        (MAY BE DEFINED USING I1MACH(14).) */
/*     INTEGER TS */
/*        THE NUMBER OF BASE-B DIGITS IN SINGLE PRECISION. */
/*        (MAY BE DEFINED USING I1MACH(11).) */
/*   MACHINE CONSTANTS FOR COMPUTERS FOLLOWING IEEE ARITHMETIC STANDARD */
/*   (E.G., MOTOROLA 68000 BASED MACHINES SUCH AS SUN AND SPARC */
/*   WORKSTATIONS, AND AT&T PC 7300; AND 8087 BASED MICROS SUCH AS THE */
/*   IBM PC AND THE AT&T 6300). */
/*   MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM. */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  60 / */
/*   MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */
/*                         THE BURROUGHS 6700/7700 SYSTEMS */
/*     DATA B  /   8 / */
/*     DATA TS /  13 / */
/*     DATA TD /  26 / */
/*   MACHINE CONSTANTS FOR THE CDC 6000/7000 (FTN5 COMPILER) */
/*                         THE CYBER 170/180 SERIES UNDER NOS */
/*     DATA B  /   2 / */
/*     DATA TS /  48 / */
/*     DATA TD /  96 / */
/*   MACHINE CONSTANTS FOR THE CDC 6000/7000 (FTN COMPILER) */
/*                         THE CYBER 170/180 SERIES UNDER NOS/VE */
/*                         THE CYBER 200 SERIES */
/*     DATA B  /   2 / */
/*     DATA TS /  47 / */
/*     DATA TD /  94 / */
/*   MACHINE CONSTANTS FOR THE CRAY */
/*     DATA B  /   2 / */
/*     DATA TS /  47 / */
/*     DATA TD /  94 / */
/*   MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     DATA B  /  16 / */
/*     DATA TS /   6 / */
/*     DATA TD /  14 / */
/*   MACHINE CONSTANTS FOR THE HARRIS COMPUTER */
/*     DATA B  /   2 / */
/*     DATA TS /  23 / */
/*     DATA TD /  38 / */
/*   MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 */
/*                         THE HONEYWELL 600/6000 SERIES */
/*     DATA B  /   2 / */
/*     DATA TS /  27 / */
/*     DATA TD /  63 / */
/*   MACHINE CONSTANTS FOR THE HP 2100 */
/*      (3 WORD DOUBLE PRECISION OPTION WITH FTN4) */
/*     DATA B  /   2 / */
/*     DATA TS /  23 / */
/*     DATA TD /  39 / */
/*   MACHINE CONSTANTS FOR THE HP 2100 */
/*      (4 WORD DOUBLE PRECISION OPTION WITH FTN4) */
/*     DATA B  /   2 / */
/*     DATA TS /  23 / */
/*     DATA TD /  55 / */
/*   MACHINE CONSTANTS FOR THE IBM 360/370 SERIES */
/*     DATA B  /  16 / */
/*     DATA TS /   6 / */
/*     DATA TD /  14 / */
/*   MACHINE CONSTANTS FOR THE IBM PC */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  53 / */
/*   MACHINE CONSTANTS FOR THE INTERDATA (PERKIN ELMER) 7/32 */
/*                             INTERDATA (PERKIN ELMER) 8/32 */
/*     DATA B  /  16 / */
/*     DATA TS /   6 / */
/*     DATA TD /  14 / */
/*   MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */
/*     DATA B  /   2 / */
/*     DATA TS /  27 / */
/*     DATA TD /  54 / */
/*   MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */
/*     DATA B  /   2 / */
/*     DATA TS /  27 / */
/*     DATA TD /  62 / */
/*   MACHINE CONSTANTS FOR THE PDP-11 SYSTEM */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  56 / */
/*   MACHINE CONSTANTS FOR THE PERKIN-ELMER 3230 */
/*     DATA B  /  16 / */
/*     DATA TS /   6 / */
/*     DATA TD /  14 / */
/*   MACHINE CONSTANTS FOR THE PRIME 850 AND PRIME 4050 */
/*     DATA B  /   2 / */
/*     DATA TS /  23 / */
/*     DATA TD /  47 / */
/*   MACHINE CONSTANTS FOR THE SEL SYSTEMS 85/86 */
/*     DATA B  /  16 / */
/*     DATA TS /   6 / */
/*     DATA TD /  14 / */
/*   MACHINE CONSTANTS FOR SUN AND SPARC WORKSTATIONS */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  53 / */
/*   MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES */
/*     DATA B  /   2 / */
/*     DATA TS /  27 / */
/*     DATA TD /  60 / */
/*   MACHINE CONSTANTS FOR THE VAX-11 WITH FORTRAN IV-PLUS COMPILER */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  56 / */
/*   MACHINE CONSTANTS FOR THE VAX/VMS SYSTEM WITHOUT  G_FLOATING */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  56 / */
/*   MACHINE CONSTANTS FOR THE VAX/VMS SYSTEM WITH G_FLOATING */
/*     DATA B  /   2 / */
/*     DATA TS /  24 / */
/*     DATA TD /  53 / */
/*   MACHINE CONSTANTS FOR THE XEROX SIGMA 5/7/9 */
/*     DATA B  /  16 / */
/*     DATA TS /   6 / */
/*     DATA TD /  14 / */
/* ***FIRST EXECUTABLE STATEMENT  DMPREC */
/*     B = I1MACH(10) */
/*     TS = I1MACH(11) */
/*     TD = I1MACH(14) */
    i__1 = 1 - td;
    ret_val = pow_di(&b, &i__1);
    return ret_val;
} /* dmprec_ */

