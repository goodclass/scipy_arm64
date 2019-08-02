/* dqagi.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dqagi_(D_fp f, doublereal *bound, integer *inf, 
	doublereal *epsabs, doublereal *epsrel, doublereal *result, 
	doublereal *abserr, integer *neval, integer *ier, integer *limit, 
	integer *lenw, integer *last, integer *iwork, doublereal *work)
{
    static integer l1, l2, l3, lvl;
    extern /* Subroutine */ int dqagie_(D_fp, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), xerror_(char *, integer *, 
	    integer *, integer *, ftnlen);

/* ***begin prologue  dqagi */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a3a1,h2a4a1 */
/* ***keywords  automatic integrator, infinite intervals, */
/*             general-purpose, transformation, extrapolation, */
/*             globally adaptive */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. -k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            integral   i = integral of f over (bound,+infinity) */
/*            or i = integral of f over (-infinity,bound) */
/*            or i = integral of f over (-infinity,+infinity) */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result).le.max(epsabs,epsrel*abs(i)). */
/* ***description */

/*        integration over infinite intervals */
/*        standard fortran subroutine */

/*        parameters */
/*         on entry */
/*            f      - double precision */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            bound  - double precision */
/*                     finite bound of integration range */
/*                     (has no meaning if interval is doubly-infinite) */

/*            inf    - integer */
/*                     indicating the kind of integration range involved */
/*                     inf = 1 corresponds to  (bound,+infinity), */
/*                     inf = -1            to  (-infinity,bound), */
/*                     inf = 2             to (-infinity,+infinity). */

/*            epsabs - double precision */
/*                     absolute accuracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */


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
/*                   - ier.gt.0 abnormal termination of the routine. the */
/*                             estimates for result and error are less */
/*                             reliable. it is assumed that the requested */
/*                             accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking the according dimension */
/*                             adjustments into account). however, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. if */
/*                             the position of a local difficulty can be */
/*                             determined (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. if possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used, which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 the occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             the error may be under-estimated. */
/*                         = 3 extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 4 the algorithm does not converge. */
/*                             roundoff error is detected in the */
/*                             extrapolation table. */
/*                             it is assumed that the requested tolerance */
/*                             cannot be achieved, and that the returned */
/*                             result is the best which can be obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier. */
/*                         = 6 the input is invalid, because */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
/*                              or limit.lt.1 or leniw.lt.limit*4. */
/*                             result, abserr, neval, last are set to */
/*                             zero. except when limit or leniw is */
/*                             invalid, iwork(1), work(limit*2+1) and */
/*                             work(limit*3+1) are set to zero, work(1) */
/*                             is set to a and work(limit+1) to b. */

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
/*                    if lenw.lt.limit*4, the routine will end */
/*                    with ier = 6. */

/*            last  - integer */
/*                    on return, last equals the number of subintervals */
/*                    produced in the subdivision process, which */
/*                    determines the number of significant elements */
/*                    actually in the work arrays. */

/*         work arrays */
/*            iwork - integer */
/*                    vector of dimension at least limit, the first */
/*                    k elements of which contain pointers */
/*                    to the error estimates over the subintervals, */
/*                    such that work(limit*3+iwork(1)),... , */
/*                    work(limit*3+iwork(k)) form a decreasing */
/*                    sequence, with k = last if last.le.(limit/2+2), and */
/*                    k = limit+1-last otherwise */

/*            work  - double precision */
/*                    vector of dimension at least lenw */
/*                    on return */
/*                    work(1), ..., work(last) contain the left */
/*                     end points of the subintervals in the */
/*                     partition of (a,b), */
/*                    work(limit+1), ..., work(limit+last) contain */
/*                     the right end points, */
/*                    work(limit*2+1), ...,work(limit*2+last) contain the */
/*                     integral approximations over the subintervals, */
/*                    work(limit*3+1), ..., work(limit*3+last) */
/*                     contain the error estimates. */
/* ***references  (none) */
/* ***routines called  dqagie,xerror */
/* ***end prologue  dqagi */




/*         check validity of limit and lenw. */

/* ***first executable statement  dqagi */
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

/*         prepare call for dqagie. */

    l1 = *limit + 1;
    l2 = *limit + l1;
    l3 = *limit + l2;

    dqagie_((D_fp)f, bound, inf, epsabs, epsrel, limit, result, abserr, neval,
	     ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

/*         call error handler if necessary. */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xerror_("abnormal return from dqagi", &c__26, ier, &lvl, (ftnlen)26);
    }
    return 0;
} /* dqagi_ */

