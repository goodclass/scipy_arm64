/* xerror.f -- translated by f2c (version 20190311).
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

static integer c__9 = 9;
static integer c__1 = 1;

/* Subroutine */ int xerror_(char *mess, integer *nmess, integer *l1, integer 
	*l2, ftnlen mess_len)
{
    /* Format strings */
    static char fmt_900[] = "(/)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), do_lio(integer *
	    , integer *, char *, ftnlen), e_wsle(void);

    /* Local variables */
    static integer i__, k, nn, nr, kmin;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, fmt_900, 0 };



/*     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS */
/*     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL */
/*     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77 */
/*     ROUTINE. */

    nn = *nmess / 70;
    nr = *nmess - nn * 70;
    if (nr != 0) {
	++nn;
    }
    k = 1;
    s_wsfe(&io___4);
    e_wsfe();
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = k + 69;
	kmin = min(i__2,*nmess);
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, mess + (k - 1), kmin - (k - 1));
	e_wsle();
	k += 70;
/* L10: */
    }
    s_wsfe(&io___8);
    e_wsfe();
    return 0;
} /* xerror_ */

