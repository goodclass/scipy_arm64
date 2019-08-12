/* dlunoc.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dluno_(integer *lun, char *fn, ftnlen fn_len)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);

    o__1.oerr = 0;
    o__1.ounit = *lun;
    o__1.ofnmlen = fn_len;
    o__1.ofnm = fn;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    return 0;
} /* dluno_ */

/* Subroutine */ int dlunc_(integer *lun)
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    integer f_clos(cllist *);

    cl__1.cerr = 0;
    cl__1.cunit = *lun;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* dlunc_ */

