/* xsetun.f -- translated by f2c (version 20190311).
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
    integer mesflg, lunit;
} eh0001_;

#define eh0001_1 eh0001_

/* Subroutine */ int xsetun_(integer *lun)
{

/* this routine resets the logical unit number for messages. */


    if (*lun > 0) {
	eh0001_1.lunit = *lun;
    }
    return 0;
/* ----------------------- end of subroutine xsetun ---------------------- */
} /* xsetun_ */

