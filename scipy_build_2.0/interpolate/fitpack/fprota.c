/* fprota.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fprota_(doublereal *cos__, doublereal *sin__, doublereal 
	*a, doublereal *b)
{
    static doublereal stor1, stor2;

/*  subroutine fprota applies a givens rotation to a and b. */
/*  .. */
/*  ..scalar arguments.. */
/* ..local scalars.. */
/*  .. */
    stor1 = *a;
    stor2 = *b;
    *b = *cos__ * stor2 + *sin__ * stor1;
    *a = *cos__ * stor1 - *sin__ * stor2;
    return 0;
} /* fprota_ */

