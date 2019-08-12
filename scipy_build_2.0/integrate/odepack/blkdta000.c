/* blkdta000.f -- translated by f2c (version 20190311).
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

struct ls0001_1_ {
    doublereal rowns[209], rcomm[9];
    integer illin, iduma[10], ntrep, idumb[2], iowns[6], icomm[19];
};

#define ls0001_1 (*(struct ls0001_1_ *) &ls0001_)

struct eh0001_1_ {
    integer mesflg, lunit;
};

#define eh0001_1 (*(struct eh0001_1_ *) &eh0001_)

/* Initialized data */

struct {
    doublereal fill_1[218];
    integer e_2;
    integer fill_3[10];
    integer e_4;
    integer fill_5[27];
    } ls0001_ = { {0}, 0, {0}, 0 };

struct {
    integer e_1[2];
    } eh0001_ = { 1, 6 };


/* ----------------------------------------------------------------------- */
/* this data subprogram loads variables into the internal common */
/* blocks used by the odepack solvers.  the variables are */
/* defined as follows.. */
/*   illin  = counter for the number of consecutive times the package */
/*            was called with illegal input.  the run is stopped when */
/*            illin reaches 5. */
/*   ntrep  = counter for the number of consecutive times the package */
/*            was called with istate = 1 and tout = t.  the run is */
/*            stopped when ntrep reaches 5. */
/*   mesflg = flag to control printing of error messages.  1 means print, */
/*            0 means no printing. */
/*   lunit  = default value of logical unit number for printing of error */
/*            messages. */
/* ----------------------------------------------------------------------- */

/* ----------------------- end of block data ----------------------------- */

