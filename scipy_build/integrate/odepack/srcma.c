/* srcma.f -- translated by f2c (version 20190311).
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

/* Common Block Declarations */

struct {
    doublereal rls[218];
    integer ils[39];
} ls0001_;

#define ls0001_1 ls0001_

struct {
    doublereal rlsa[22];
    integer ilsa[9];
} lsa001_;

#define lsa001_1 lsa001_

struct {
    integer ieh[2];
} eh0001_;

#define eh0001_1 eh0001_

/* Subroutine */ int srcma_(doublereal *rsav, integer *isav, integer *job)
{
    /* Initialized data */

    static integer lenrls = 218;
    static integer lenils = 39;
    static integer lenrla = 22;
    static integer lenila = 9;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ----------------------------------------------------------------------- */
/* this routine saves or restores (depending on job) the contents of */
/* the common blocks ls0001, lsa001, and eh0001, which are used */
/* internally by one or more odepack solvers. */

/* rsav = real array of length 240 or more. */
/* isav = integer array of length 50 or more. */
/* job  = flag indicating to save or restore the common blocks.. */
/*        job  = 1 if common is to be saved (written to rsav/isav) */
/*        job  = 2 if common is to be restored (read from rsav/isav) */
/*        a call with job = 2 presumes a prior call with job = 1. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

    if (*job == 2) {
	goto L100;
    }
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	rsav[i__] = ls0001_1.rls[i__ - 1];
    }
    i__1 = lenrla;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	rsav[lenrls + i__] = lsa001_1.rlsa[i__ - 1];
    }

    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	isav[i__] = ls0001_1.ils[i__ - 1];
    }
    i__1 = lenila;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	isav[lenils + i__] = lsa001_1.ilsa[i__ - 1];
    }

    isav[lenils + lenila + 1] = eh0001_1.ieh[0];
    isav[lenils + lenila + 2] = eh0001_1.ieh[1];
    return 0;

L100:
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	ls0001_1.rls[i__ - 1] = rsav[i__];
    }
    i__1 = lenrla;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L115: */
	lsa001_1.rlsa[i__ - 1] = rsav[lenrls + i__];
    }

    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	ls0001_1.ils[i__ - 1] = isav[i__];
    }
    i__1 = lenila;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	lsa001_1.ilsa[i__ - 1] = isav[lenils + i__];
    }

    eh0001_1.ieh[0] = isav[lenils + lenila + 1];
    eh0001_1.ieh[1] = isav[lenils + lenila + 2];
    return 0;
/* ----------------------- end of subroutine srcma ----------------------- */
} /* srcma_ */

