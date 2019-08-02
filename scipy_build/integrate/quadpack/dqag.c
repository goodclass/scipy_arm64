/* dqag.f -- translated by f2c (version 20190311).
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

static integer c__26 = 26;

/* Subroutine */ int dqag_(D_fp f, doublereal *a, doublereal *b, doublereal *
	epsabs, doublereal *epsrel, integer *key, doublereal *result, 
	doublereal *abserr, integer *neval, integer *ier, integer *limit, 
	integer *lenw, integer *last, integer *iwork, doublereal *work)
{
    static integer l1, l2, l3, lvl;
    extern /* Subroutine */ int dqage_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), xerror_(char *,
	     integer *, integer *, integer *, ftnlen);

/* ***begin prologue  dqag */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a1 */
/* ***keywords  automatic integrator, general-purpose, */
/*             integrand examinator, globally adaptive, */
/*             gauss-kronrod */
/* ***author  piessens,robert,appl. math. & progr. div - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral i = integral of f over (a,b), */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result)le.max(epsabs,epsrel*abs(i)). */
/* ***description */

/*        computation of a definite integral */
/*        standard fortran subroutine */
/*        double precision version */

/*            f      - double precision */
/*                     function subprogam defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            a      - double precision */
/*                     lower limit of integration */

/*            b      - double precision */
/*                     upper limit of integration */

/*            epsabs - double precision */
/*                     absolute accoracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            key    - integer */
/*                     key for choice of local integration rule */
/*                     a gauss-kronrod pair is used with */
/*                       7 - 15 points if key.lt.2, */
/*                      10 - 21 points if key = 2, */
/*                      15 - 31 points if key = 3, */
/*                      20 - 41 points if key = 4, */
/*                      25 - 51 points if key = 5, */
/*                      30 - 61 points if key.gt.5. */

/*         on return */
/*            result - double precision */
/*                     approximation to the integral */

/*            abserr - double precision */
/*                     estimate of the modulus of the absolute error, */
/*                     which should equal or exceed abs(i-result) */

/*            neval  - integer */
/*                     number of integrand evaluations */

/*            ier    - integer */
/*                     ier = 0 normal and reliable termination of the */
/*                             routine. it is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     ier.gt.0 abnormal termination of the routine */
/*                             the estimates for result and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*                      error messages */
/*                     ier = 1 maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking the according dimension */
/*                             adjustments into account). however, if */
/*                             this yield no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulaties. */
/*                             if the position of a local difficulty can */
/*                             be determined (i.e.singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. if possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 the occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 the input is invalid, because */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
/*                             or limit.lt.1 or lenw.lt.limit*4. */
/*                             result, abserr, neval, last are set */
/*                             to zero. */
/*                             except when lenw is invalid, iwork(1), */
/*                             work(limit*2+1) and work(limit*3+1) are */
/*                             set to zero, work(1) is set to a and */
/*                             work(limit+1) to b. */

/*         dimensioning parameters */
/*            limit - integer */
/*                    dimensioning parameter for iwork */
/*                    limit determines the maximum number of subintervals */
/*                    in the partition of the given integration interval */
/*                    (a,b), limit.ge.1. */
/*                    if limit.lt.1, the routine will end with ier = 6. */

/*            lenw  - integer */
/*                    dimensioning parameter for work */
/*                    lenw must be at least limit*4. */
/*                    if lenw.lt.limit*4, the routine will end with */
/*                    ier = 6. */

/*            last  - integer */
/*                    on return, last equals the number of subintervals */
/*                    produced in the subdiviosion process, which */
/*                    determines the number of significant elements */
/*                    actually in the work arrays. */

/*         work arrays */
/*            iwork - integer */
/*                    vector of dimension at least limit, the first k */
/*                    elements of which contain pointers to the error */
/*                    estimates over the subintervals, such that */
/*                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k)) */
/*                    form a decreasing sequence with k = last if */
/*                    last.le.(limit/2+2), and k = limit+1-last otherwise */

/*            work  - double precision */
/*                    vector of dimension at least lenw */
/*                    on return */
/*                    work(1), ..., work(last) contain the left end */
/*                    points of the subintervals in the partition of */
/*                     (a,b), */
/*                    work(limit+1), ..., work(limit+last) contain the */
/*                     right end points, */
/*                    work(limit*2+1), ..., work(limit*2+last) contain */
/*                     the integral approximations over the subintervals, */
/*                    work(limit*3+1), ..., work(limit*3+last) contain */
/*                     the error estimates. */

/* ***references  (none) */
/* ***routines called  dqage,xerror */
/* ***end prologue  dqag */



/*         check validity of lenw. */

/* ***first executable statement  dqag */
    /* Parameter adjustments */
    --iwork;
    --work;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*limit < 1 || *lenw < *limit << 2) {
	goto L10;
    }

/*         prepare call for dqage. */

    l1 = *limit + 1;
    l2 = *limit + l1;
    l3 = *limit + l2;

    dqage_((D_fp)f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, 
	    ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

/*         call error handler if necessary. */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xerror_("abnormal return from dqag", &c__26, ier, &lvl, (ftnlen)25);
    }
    return 0;
} /* dqag_ */

