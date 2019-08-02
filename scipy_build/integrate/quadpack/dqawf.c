/* dqawf.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dqawf_(D_fp f, doublereal *a, doublereal *omega, integer 
	*integr, doublereal *epsabs, doublereal *result, doublereal *abserr, 
	integer *neval, integer *ier, integer *limlst, integer *lst, integer *
	leniw, integer *maxp1, integer *lenw, integer *iwork, doublereal *
	work)
{
    static integer l1, l2, l3, l4, l5, l6, ll2, lvl, last, limit;
    extern /* Subroutine */ int dqawfe_(D_fp, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *), 
	    xerror_(char *, integer *, integer *, integer *, ftnlen);

/* ***begin prologue  dqawf */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a3a1 */
/* ***keywords  automatic integrator, special-purpose,fourier */
/*             integral, integration between zeros with dqawoe, */
/*             convergence acceleration with dqelg */
/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            fourier integral i=integral of f(x)*w(x) over (a,infinity) */
/*            where w(x) = cos(omega*x) or w(x) = sin(omega*x). */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result).le.epsabs. */
/* ***description */

/*        computation of fourier integrals */
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

/*            omega  - double precision */
/*                     parameter in the integrand weight function */

/*            integr - integer */
/*                     indicates which of the weight functions is used */
/*                     integr = 1      w(x) = cos(omega*x) */
/*                     integr = 2      w(x) = sin(omega*x) */
/*                     if integr.ne.1.and.integr.ne.2, the routine */
/*                     will end with ier = 6. */

/*            epsabs - double precision */
/*                     absolute accuracy requested, epsabs.gt.0. */
/*                     if epsabs.le.0, the routine will end with ier = 6. */

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
/*                     ier.gt.0 abnormal termination of the routine. */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                    if omega.ne.0 */
/*                     ier = 1 maximum number of cycles allowed */
/*                             has been achieved, i.e. of subintervals */
/*                             (a+(k-1)c,a+kc) where */
/*                             c = (2*int(abs(omega))+1)*pi/abs(omega), */
/*                             for k = 1, 2, ..., lst. */
/*                             one can allow more cycles by increasing */
/*                             the value of limlst (and taking the */
/*                             according dimension adjustments into */
/*                             account). examine the array iwork which */
/*                             contains the error flags on the cycles, in */
/*                             order to look for eventual local */
/*                             integration difficulties. */
/*                             if the position of a local difficulty */
/*                             can be determined (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling */
/*                             appropriate integrators on the subranges. */
/*                         = 4 the extrapolation table constructed for */
/*                             convergence accelaration of the series */
/*                             formed by the integral contributions over */
/*                             the cycles, does not converge to within */
/*                             the requested accuracy. */
/*                             as in the case of ier = 1, it is advised */
/*                             to examine the array iwork which contains */
/*                             the error flags on the cycles. */
/*                         = 6 the input is invalid because */
/*                             (integr.ne.1 and integr.ne.2) or */
/*                              epsabs.le.0 or limlst.lt.1 or */
/*                              leniw.lt.(limlst+2) or maxp1.lt.1 or */
/*                              lenw.lt.(leniw*2+maxp1*25). */
/*                              result, abserr, neval, lst are set to */
/*                              zero. */
/*                         = 7 bad integrand behaviour occurs within */
/*                             one or more of the cycles. location and */
/*                             type of the difficulty involved can be */
/*                             determined from the first lst elements of */
/*                             vector iwork.  here lst is the number of */
/*                             cycles actually needed (see below). */
/*                             iwork(k) = 1 the maximum number of */
/*                                          subdivisions (=(leniw-limlst) */
/*                                          /2) has been achieved on the */
/*                                          k th cycle. */
/*                                      = 2 occurrence of roundoff error */
/*                                          is detected and prevents the */
/*                                          tolerance imposed on the k th */
/*                                          cycle, from being achieved */
/*                                          on this cycle. */
/*                                      = 3 extremely bad integrand */
/*                                          behaviour occurs at some */
/*                                          points of the k th cycle. */
/*                                      = 4 the integration procedure */
/*                                          over the k th cycle does */
/*                                          not converge (to within the */
/*                                          required accuracy) due to */
/*                                          roundoff in the extrapolation */
/*                                          procedure invoked on this */
/*                                          cycle. it is assumed that the */
/*                                          result on this interval is */
/*                                          the best which can be */
/*                                          obtained. */
/*                                      = 5 the integral over the k th */
/*                                          cycle is probably divergent */
/*                                          or slowly convergent. it must */
/*                                          be noted that divergence can */
/*                                          occur with any other value of */
/*                                          iwork(k). */
/*                    if omega = 0 and integr = 1, */
/*                    the integral is calculated by means of dqagie, */
/*                    and ier = iwork(1) (with meaning as described */
/*                    for iwork(k),k = 1). */

/*         dimensioning parameters */
/*            limlst - integer */
/*                     limlst gives an upper bound on the number of */
/*                     cycles, limlst.ge.3. */
/*                     if limlst.lt.3, the routine will end with ier = 6. */

/*            lst    - integer */
/*                     on return, lst indicates the number of cycles */
/*                     actually needed for the integration. */
/*                     if omega = 0, then lst is set to 1. */

/*            leniw  - integer */
/*                     dimensioning parameter for iwork. on entry, */
/*                     (leniw-limlst)/2 equals the maximum number of */
/*                     subintervals allowed in the partition of each */
/*                     cycle, leniw.ge.(limlst+2). */
/*                     if leniw.lt.(limlst+2), the routine will end with */
/*                     ier = 6. */

/*            maxp1  - integer */
/*                     maxp1 gives an upper bound on the number of */
/*                     chebyshev moments which can be stored, i.e. for */
/*                     the intervals of lengths abs(b-a)*2**(-l), */
/*                     l = 0,1, ..., maxp1-2, maxp1.ge.1. */
/*                     if maxp1.lt.1, the routine will end with ier = 6. */
/*            lenw   - integer */
/*                     dimensioning parameter for work */
/*                     lenw must be at least leniw*2+maxp1*25. */
/*                     if lenw.lt.(leniw*2+maxp1*25), the routine will */
/*                     end with ier = 6. */

/*         work arrays */
/*            iwork  - integer */
/*                     vector of dimension at least leniw */
/*                     on return, iwork(k) for k = 1, 2, ..., lst */
/*                     contain the error flags on the cycles. */

/*            work   - double precision */
/*                     vector of dimension at least */
/*                     on return, */
/*                     work(1), ..., work(lst) contain the integral */
/*                      approximations over the cycles, */
/*                     work(limlst+1), ..., work(limlst+lst) contain */
/*                      the error extimates over the cycles. */
/*                     further elements of work have no specific */
/*                     meaning for the user. */

/* ***references  (none) */
/* ***routines called  dqawfe,xerror */
/* ***end prologue  dqawf */




/*         check validity of limlst, leniw, maxp1 and lenw. */

/* ***first executable statement  dqawf */
    /* Parameter adjustments */
    --iwork;
    --work;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*limlst < 3 || *leniw < *limlst + 2 || *maxp1 < 1 || *lenw < (*leniw 
	    << 1) + *maxp1 * 25) {
	goto L10;
    }

/*         prepare call for dqawfe */

    limit = (*leniw - *limlst) / 2;
    l1 = *limlst + 1;
    l2 = *limlst + l1;
    l3 = limit + l2;
    l4 = limit + l3;
    l5 = limit + l4;
    l6 = limit + l5;
    ll2 = limit + l1;
    dqawfe_((D_fp)f, a, omega, integr, epsabs, limlst, &limit, maxp1, result, 
	    abserr, neval, ier, &work[1], &work[l1], &iwork[1], lst, &work[l2]
	    , &work[l3], &work[l4], &work[l5], &iwork[l1], &iwork[ll2], &work[
	    l6]);

/*         call error handler if necessary */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 1;
    }
    if (*ier != 0) {
	xerror_("abnormal return from dqawf", &c__26, ier, &lvl, (ftnlen)26);
    }
    return 0;
} /* dqawf_ */

