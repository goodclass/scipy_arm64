/* dqaws.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dqaws_(D_fp f, doublereal *a, doublereal *b, doublereal *
	alfa, doublereal *beta, integer *integr, doublereal *epsabs, 
	doublereal *epsrel, doublereal *result, doublereal *abserr, integer *
	neval, integer *ier, integer *limit, integer *lenw, integer *last, 
	integer *iwork, doublereal *work)
{
    static integer l1, l2, l3, lvl;
    extern /* Subroutine */ int dqawse_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), xerror_(char *, integer *, integer *, integer *, 
	    ftnlen);

/* ***begin prologue  dqaws */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1 */
/* ***keywords  automatic integrator, special-purpose, */
/*             algebraico-logarithmic end-point singularities, */
/*             clenshaw-curtis, globally adaptive */
/* ***author  piessens,robert,appl. math. & progr. div. -k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral i = integral of f*w over (a,b), */
/*            (where w shows a singular behaviour at the end points */
/*            see parameter integr). */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result).le.max(epsabs,epsrel*abs(i)). */
/* ***description */

/*        integration of functions having algebraico-logarithmic */
/*        end point singularities */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*         on entry */
/*            f      - double precision */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            a      - double precision */
/*                     lower limit of integration */

/*            b      - double precision */
/*                     upper limit of integration, b.gt.a */
/*                     if b.le.a, the routine will end with ier = 6. */

/*            alfa   - double precision */
/*                     parameter in the integrand function, alfa.gt.(-1) */
/*                     if alfa.le.(-1), the routine will end with */
/*                     ier = 6. */

/*            beta   - double precision */
/*                     parameter in the integrand function, beta.gt.(-1) */
/*                     if beta.le.(-1), the routine will end with */
/*                     ier = 6. */

/*            integr - integer */
/*                     indicates which weight function is to be used */
/*                     = 1  (x-a)**alfa*(b-x)**beta */
/*                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a) */
/*                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x) */
/*                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x) */
/*                     if integr.lt.1 or integr.gt.4, the routine */
/*                     will end with ier = 6. */

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
/*                     ier.gt.0 abnormal termination of the routine */
/*                             the estimates for the integral and error */
/*                             are less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking the according dimension */
/*                             adjustments into account). however, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand, in order to */
/*                             determine the integration difficulties */
/*                             which prevent the requested tolerance from */
/*                             being achieved. in case of a jump */
/*                             discontinuity or a local singularity */
/*                             of algebraico-logarithmic type at one or */
/*                             more interior points of the integration */
/*                             range, one should proceed by splitting up */
/*                             the interval at these points and calling */
/*                             the integrator on the subranges. */
/*                         = 2 the occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 6 the input is invalid, because */
/*                             b.le.a or alfa.le.(-1) or beta.le.(-1) or */
/*                             or integr.lt.1 or integr.gt.4 or */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
/*                             or limit.lt.2 or lenw.lt.limit*4. */
/*                             result, abserr, neval, last are set to */
/*                             zero. except when lenw or limit is invalid */
/*                             iwork(1), work(limit*2+1) and */
/*                             work(limit*3+1) are set to zero, work(1) */
/*                             is set to a and work(limit+1) to b. */

/*         dimensioning parameters */
/*            limit  - integer */
/*                     dimensioning parameter for iwork */
/*                     limit determines the maximum number of */
/*                     subintervals in the partition of the given */
/*                     integration interval (a,b), limit.ge.2. */
/*                     if limit.lt.2, the routine will end with ier = 6. */

/*            lenw   - integer */
/*                     dimensioning parameter for work */
/*                     lenw must be at least limit*4. */
/*                     if lenw.lt.limit*4, the routine will end */
/*                     with ier = 6. */

/*            last   - integer */
/*                     on return, last equals the number of */
/*                     subintervals produced in the subdivision process, */
/*                     which determines the significant number of */
/*                     elements actually in the work arrays. */

/*         work arrays */
/*            iwork  - integer */
/*                     vector of dimension limit, the first k */
/*                     elements of which contain pointers */
/*                     to the error estimates over the subintervals, */
/*                     such that work(limit*3+iwork(1)), ..., */
/*                     work(limit*3+iwork(k)) form a decreasing */
/*                     sequence with k = last if last.le.(limit/2+2), */
/*                     and k = limit+1-last otherwise */

/*            work   - double precision */
/*                     vector of dimension lenw */
/*                     on return */
/*                     work(1), ..., work(last) contain the left */
/*                      end points of the subintervals in the */
/*                      partition of (a,b), */
/*                     work(limit+1), ..., work(limit+last) contain */
/*                      the right end points, */
/*                     work(limit*2+1), ..., work(limit*2+last) */
/*                      contain the integral approximations over */
/*                      the subintervals, */
/*                     work(limit*3+1), ..., work(limit*3+last) */
/*                      contain the error estimates. */

/* ***references  (none) */
/* ***routines called  dqawse,xerror */
/* ***end prologue  dqaws */




/*         check validity of limit and lenw. */

/* ***first executable statement  dqaws */
    /* Parameter adjustments */
    --iwork;
    --work;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*limit < 2 || *lenw < *limit << 2) {
	goto L10;
    }

/*         prepare call for dqawse. */

    l1 = *limit + 1;
    l2 = *limit + l1;
    l3 = *limit + l2;

    dqawse_((D_fp)f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, 
	    abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
	    iwork[1], last);

/*         call error handler if necessary. */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xerror_("abnormal return from dqaws", &c__26, ier, &lvl, (ftnlen)26);
    }
    return 0;
} /* dqaws_ */

