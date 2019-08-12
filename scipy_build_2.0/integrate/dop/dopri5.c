/* dopri5.f -- translated by f2c (version 20190311).
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

union {
    struct {
	doublereal xold, hout;
    } _1;
    struct {
	doublereal xold, h__;
    } _2;
} condo5_;

#define condo5_1 (condo5_._1)
#define condo5_2 (condo5_._2)

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b47 = 1.;

/* Subroutine */ int dopri5_(integer *n, U_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *rtol, doublereal *atol, integer *
	itol, U_fp solout, integer *iout, doublereal *work, integer *lwork, 
	integer *iwork, integer *liwork, doublereal *rpar, integer *ipar, 
	integer *idid)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal fac1, fac2;
    static integer iek1, iek2, iek3, iek4, iek5, iek6, iey1;
    static doublereal beta, safe;
    static integer ieco, nfcn, meth;
    static doublereal hmax;
    static integer nmax, ieys, icomp;
    static logical arret;
    static integer nstep, naccpt, nrejct;
    extern /* Subroutine */ int dopcor_(integer *, U_fp, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, U_fp, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer nstiff, nrdens, iprint, istore;
    static doublereal uround;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 0, 0, 0, 0 };
    static cilist io___10 = { 0, 0, 0, 0, 0 };
    static cilist io___13 = { 0, 0, 0, 0, 0 };
    static cilist io___14 = { 0, 0, 0, 0, 0 };
    static cilist io___17 = { 0, 0, 0, 0, 0 };
    static cilist io___19 = { 0, 0, 0, 0, 0 };
    static cilist io___23 = { 0, 0, 0, 0, 0 };
    static cilist io___36 = { 0, 0, 0, 0, 0 };
    static cilist io___38 = { 0, 0, 0, 0, 0 };


/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER */
/*     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y). */
/*     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER (4)5 */
/*     DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND */
/*     DENSE OUTPUT). */

/*     AUTHORS: E. HAIRER AND G. WANNER */
/*              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES */
/*              CH-1211 GENEVE 24, SWITZERLAND */
/*              E-MAIL:  Ernst.Hairer@math.unige.ch */
/*                       Gerhard.Wanner@math.unige.ch */

/*     THIS CODE IS DESCRIBED IN: */
/*         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY */
/*         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION. */
/*         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS, */
/*         SPRINGER-VERLAG (1993) */

/*     VERSION OF APRIL 25, 1996 */
/*     (latest correction of a small bug: August 8, 2005) */

/*     INPUT PARAMETERS */
/*     ---------------- */
/*     N           DIMENSION OF THE SYSTEM */

/*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
/*                 VALUE OF F(X,Y): */
/*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),F(N) */
/*                    F(1)=...   ETC. */

/*     X           INITIAL X-VALUE */

/*     Y(N)        INITIAL VALUES FOR Y */

/*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

/*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY */
/*                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. */

/*     ITOL        SWITCH FOR RTOL AND ATOL: */
/*                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. */
/*                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF */
/*                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL */
/*                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. */
/*                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW */
/*                     RTOL(I)*ABS(Y(I))+ATOL(I). */

/*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
/*                 NUMERICAL SOLUTION DURING INTEGRATION. */
/*                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
/*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
/*                 IT MUST HAVE THE FORM */
/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND, */
/*                                       RPAR,IPAR,IRTRN) */
/*                    DIMENSION Y(N),CON(5*ND),ICOMP(ND) */
/*                    .... */
/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
/*                    THE FIRST GRID-POINT). */
/*                 "XOLD" IS THE PRECEDING GRID-POINT. */
/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
/*                    IS SET <0, DOPRI5 WILL RETURN TO THE CALLING PROGRAM. */
/*                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT, */
/*                    SET  IRTRN = 2 */

/*          -----  CONTINUOUS OUTPUT: ----- */
/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
/*                 THE FUNCTION */
/*                        >>>   CONTD5(I,S,CON,ICOMP,ND)   <<< */
/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */

/*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT. */
/*                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT */
/*                            (IN THIS CASE WORK(5) MUST BE SPECIFIED) */

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
/*                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE. */
/*                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING. */
/*                 "LWORK" MUST BE AT LEAST  8*N+5*NRDENS+21 */
/*                 WHERE  NRDENS = IWORK(5) */

/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

/*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
/*                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS FOR THE CODE. */
/*                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING. */
/*                 "LIWORK" MUST BE AT LEAST NRDENS+21 . */

/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
/*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
/*                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. */

/* ----------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*     ----------------------------------- */
/*              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW */
/*              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF */
/*              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES. */

/*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16. */

/*    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, */
/*              DEFAULT 0.9D0. */

/*    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION */
/*              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION */
/*                 WORK(3) <= HNEW/HOLD <= WORK(4) */
/*              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=10.D0 */

/*    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL */
/*              (SEE SECTION IV.2). LARGER VALUES OF BETA ( <= 0.1 ) */
/*              MAKE THE STEP SIZE CONTROL MORE STABLE. DOPRI5 NEEDS */
/*              A LARGER BETA THAN HIGHAM & HALL. NEGATIVE WORK(5) */
/*              PROVOKE BETA=0. */
/*              DEFAULT 0.04D0. */

/*    WORK(6)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

/*    WORK(7)   INITIAL STEP SIZE, FOR WORK(7)=0.D0 AN INITIAL GUESS */
/*              IS COMPUTED WITH HELP OF THE FUNCTION HINIT */

/*    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
/*              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000. */

/*    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS */
/*              IF IWORK(2).EQ.1  METHOD DOPRI5 OF DORMAND AND PRINCE */
/*              (TABLE 5.2 OF SECTION II.5). */
/*              AT THE MOMENT THIS IS THE ONLY POSSIBLE CHOICE. */
/*              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1. */

/*    IWORK(3)  SWITCH FOR PRINTING ERROR MESSAGES */
/*              IF IWORK(3).LT.0 NO MESSAGES ARE BEING PRINTED */
/*              IF IWORK(3).GT.0 MESSAGES ARE PRINTED WITH */
/*              WRITE (IWORK(3),*) ... */
/*              DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=6 */

/*    IWORK(4)  TEST FOR STIFFNESS IS ACTIVATED AFTER STEP NUMBER */
/*              J*IWORK(4) (J INTEGER), PROVIDED IWORK(4).GT.0. */
/*              FOR NEGATIVE IWORK(4) THE STIFFNESS TEST IS */
/*              NEVER ACTIVATED; DEFAULT VALUE IS IWORK(4)=1000 */

/*    IWORK(5)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT */
/*              IS REQUIRED; DEFAULT VALUE IS IWORK(5)=0; */
/*              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE */
/*              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN */
/*              IWORK(21),...,IWORK(NRDENS+20); */
/*              FOR  NRDENS=N  THIS IS DONE BY THE CODE. */

/* ---------------------------------------------------------------------- */

/*     OUTPUT PARAMETERS */
/*     ----------------- */
/*     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED */
/*                 (AFTER SUCCESSFUL RETURN X=XEND). */

/*     Y(N)        NUMERICAL SOLUTION AT X */

/*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

/*     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: */
/*                   IDID= 1  COMPUTATION SUCCESSFUL, */
/*                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) */
/*                   IDID=-1  INPUT IS NOT CONSISTENT, */
/*                   IDID=-2  LARGER NMAX IS NEEDED, */
/*                   IDID=-3  STEP SIZE BECOMES TOO SMALL. */
/*                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED). */

/*   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS */
/*   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS */
/*   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS */
/*   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
/*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
/* ----------------------------------------------------------------------- */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*          DECLARATIONS */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/* *** *** *** *** *** *** *** */
/*        SETTING THE PARAMETERS */
/* *** *** *** *** *** *** *** */
    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --work;
    --iwork;
    --rpar;
    --ipar;

    /* Function Body */
    nfcn = 0;
    nstep = 0;
    naccpt = 0;
    nrejct = 0;
    arret = FALSE_;
/* -------- IPRINT FOR MONITORING THE PRINTING */
    if (iwork[3] == 0) {
	iprint = 6;
    } else {
	iprint = iwork[3];
    }
/* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
    if (iwork[1] == 0) {
	nmax = 100000;
    } else {
	nmax = iwork[1];
	if (nmax <= 0) {
	    if (iprint > 0) {
		io___8.ciunit = iprint;
		s_wsle(&io___8);
		do_lio(&c__9, &c__1, " WRONG INPUT IWORK(1)=", (ftnlen)22);
		do_lio(&c__3, &c__1, (char *)&iwork[1], (ftnlen)sizeof(
			integer));
		e_wsle();
	    }
	    arret = TRUE_;
	}
    }
/* -------- METH   COEFFICIENTS OF THE METHOD */
    if (iwork[2] == 0) {
	meth = 1;
    } else {
	meth = iwork[2];
	if (meth <= 0 || meth >= 4) {
	    if (iprint > 0) {
		io___10.ciunit = iprint;
		s_wsle(&io___10);
		do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(2)=", (ftnlen)24);
		do_lio(&c__3, &c__1, (char *)&iwork[2], (ftnlen)sizeof(
			integer));
		e_wsle();
	    }
	    arret = TRUE_;
	}
    }
/* -------- NSTIFF   PARAMETER FOR STIFFNESS DETECTION */
    nstiff = iwork[4];
    if (nstiff == 0) {
	nstiff = 1000;
    }
    if (nstiff < 0) {
	nstiff = nmax + 10;
    }
/* -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS */
    nrdens = iwork[5];
    if (nrdens < 0 || nrdens > *n) {
	if (iprint > 0) {
	    io___13.ciunit = iprint;
	    s_wsle(&io___13);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(5)=", (ftnlen)24);
	    do_lio(&c__3, &c__1, (char *)&iwork[5], (ftnlen)sizeof(integer));
	    e_wsle();
	}
	arret = TRUE_;
    } else {
	if (nrdens > 0 && *iout < 2) {
	    if (iprint > 0) {
		io___14.ciunit = iprint;
		s_wsle(&io___14);
		do_lio(&c__9, &c__1, " WARNING: PUT IOUT=2 FOR DENSE OUTPUT ",
			 (ftnlen)38);
		e_wsle();
	    }
	}
	if (nrdens == *n) {
	    i__1 = nrdens;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L16: */
		iwork[i__ + 20] = i__;
	    }
	}
    }
/* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
    if (work[1] == 0.) {
	uround = 2.3e-16;
    } else {
	uround = work[1];
	if (uround <= 1e-35 || uround >= 1.) {
	    if (iprint > 0) {
		io___17.ciunit = iprint;
		s_wsle(&io___17);
		do_lio(&c__9, &c__1, " WHICH MACHINE DO YOU HAVE? YOUR UROUN"
			"D WAS:", (ftnlen)44);
		do_lio(&c__5, &c__1, (char *)&work[1], (ftnlen)sizeof(
			doublereal));
		e_wsle();
	    }
	    arret = TRUE_;
	}
    }
/* -------  SAFETY FACTOR ------------- */
    if (work[2] == 0.) {
	safe = .9;
    } else {
	safe = work[2];
	if (safe >= 1. || safe <= 1e-4) {
	    if (iprint > 0) {
		io___19.ciunit = iprint;
		s_wsle(&io___19);
		do_lio(&c__9, &c__1, " CURIOUS INPUT FOR SAFETY FACTOR WORK("
			"2)=", (ftnlen)41);
		do_lio(&c__5, &c__1, (char *)&work[2], (ftnlen)sizeof(
			doublereal));
		e_wsle();
	    }
	    arret = TRUE_;
	}
    }
/* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
    if (work[3] == 0.) {
	fac1 = .2;
    } else {
	fac1 = work[3];
    }
    if (work[4] == 0.) {
	fac2 = 10.;
    } else {
	fac2 = work[4];
    }
/* --------- BETA FOR STEP CONTROL STABILIZATION ----------- */
    if (work[5] == 0.) {
	beta = .04;
    } else {
	if (work[5] < 0.) {
	    beta = 0.;
	} else {
	    beta = work[5];
	    if (beta > .2) {
		if (iprint > 0) {
		    io___23.ciunit = iprint;
		    s_wsle(&io___23);
		    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR BETA: WORK(5)=", 
			    (ftnlen)33);
		    do_lio(&c__5, &c__1, (char *)&work[5], (ftnlen)sizeof(
			    doublereal));
		    e_wsle();
		}
		arret = TRUE_;
	    }
	}
    }
/* -------- MAXIMAL STEP SIZE */
    if (work[6] == 0.) {
	hmax = *xend - *x;
    } else {
	hmax = work[6];
    }
/* -------- INITIAL STEP SIZE */
    h__ = work[7];
/* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
    iey1 = 21;
    iek1 = iey1 + *n;
    iek2 = iek1 + *n;
    iek3 = iek2 + *n;
    iek4 = iek3 + *n;
    iek5 = iek4 + *n;
    iek6 = iek5 + *n;
    ieys = iek6 + *n;
    ieco = ieys + *n;
/* ------ TOTAL STORAGE REQUIREMENT ----------- */
    istore = ieys + nrdens * 5 - 1;
    if (istore > *lwork) {
	if (iprint > 0) {
	    io___36.ciunit = iprint;
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR WORK, MIN. LWORK="
		    , (ftnlen)43);
	    do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	    e_wsle();
	}
	arret = TRUE_;
    }
    icomp = 21;
    istore = icomp + nrdens - 1;
    if (istore > *liwork) {
	if (iprint > 0) {
	    io___38.ciunit = iprint;
	    s_wsle(&io___38);
	    do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR IWORK, MIN. LIWO"
		    "RK=", (ftnlen)45);
	    do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	    e_wsle();
	}
	arret = TRUE_;
    }
/* ------ WHEN A FAIL HAS OCCURRED, WE RETURN WITH IDID=-1 */
    if (arret) {
	*idid = -1;
	return 0;
    }
/* -------- CALL TO CORE INTEGRATOR ------------ */
    dopcor_(n, (U_fp)fcn, x, &y[1], xend, &hmax, &h__, &rtol[1], &atol[1], 
	    itol, &iprint, (U_fp)solout, iout, idid, &nmax, &uround, &meth, &
	    nstiff, &safe, &beta, &fac1, &fac2, &work[iey1], &work[iek1], &
	    work[iek2], &work[iek3], &work[iek4], &work[iek5], &work[iek6], &
	    work[ieys], &work[ieco], &iwork[icomp], &nrdens, &rpar[1], &ipar[
	    1], &nfcn, &nstep, &naccpt, &nrejct);
    work[7] = h__;
    iwork[17] = nfcn;
    iwork[18] = nstep;
    iwork[19] = naccpt;
    iwork[20] = nrejct;
/* ----------- RETURN ----------- */
    return 0;
} /* dopri5_ */




/*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

/* Subroutine */ int dopcor_(integer *n, S_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, integer *iprint, S_fp solout, 
	integer *iout, integer *idid, integer *nmax, doublereal *uround, 
	integer *meth, integer *nstiff, doublereal *safe, doublereal *beta, 
	doublereal *fac1, doublereal *fac2, doublereal *y1, doublereal *k1, 
	doublereal *k2, doublereal *k3, doublereal *k4, doublereal *k5, 
	doublereal *k6, doublereal *ysti, doublereal *cont, integer *icomp, 
	integer *nrd, doublereal *rpar, integer *ipar, integer *nfcn, integer 
	*nstep, integer *naccpt, integer *nrejct)
{
    /* Format strings */
    static char fmt_979[] = "(\002 EXIT OF DOPRI5 AT X=\002,e18.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal), pow_dd(
	    doublereal *, doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal c2, c3, c4, c5, e1, e3, e4, e5, e6, e7, d1, d3, d4, d5, 
	    d6, d7, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, 
	    a62, a63, a64, a65, a71, a73, a74, a75, a76, sk, yd0, fac, err, 
	    xph, fac11;
    static integer iord;
    static logical last;
    static doublereal hnew, bspl, facc1, facc2, expo1, hlamb, ydiff, atoli;
    static integer iasti;
    extern doublereal hinit_(integer *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    static doublereal stden, rtoli;
    static integer irtrn;
    static doublereal stnum, facold;
    static logical reject;
    extern /* Subroutine */ int cdopri_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal posneg;
    static integer nonsti;

    /* Fortran I/O blocks */
    static cilist io___99 = { 0, 0, 0, 0, 0 };
    static cilist io___103 = { 0, 0, 0, fmt_979, 0 };
    static cilist io___104 = { 0, 0, 0, 0, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_979, 0 };
    static cilist io___106 = { 0, 0, 0, 0, 0 };
    static cilist io___107 = { 0, 0, 0, fmt_979, 0 };


/* ---------------------------------------------------------- */
/*     CORE INTEGRATOR FOR DOPRI5 */
/*     PARAMETERS SAME AS IN DOPRI5 WITH WORKSPACE ADDED */
/* ---------------------------------------------------------- */
/*         DECLARATIONS */
/* ---------------------------------------------------------- */
/* *** *** *** *** *** *** *** */
/*  INITIALISATIONS */
/* *** *** *** *** *** *** *** */
    /* Parameter adjustments */
    --ysti;
    --k6;
    --k5;
    --k4;
    --k3;
    --k2;
    --k1;
    --y1;
    --y;
    --rtol;
    --atol;
    --icomp;
    --cont;
    --rpar;
    --ipar;

    /* Function Body */
    if (*meth == 1) {
	cdopri_(&c2, &c3, &c4, &c5, &e1, &e3, &e4, &e5, &e6, &e7, &a21, &a31, 
		&a32, &a41, &a42, &a43, &a51, &a52, &a53, &a54, &a61, &a62, &
		a63, &a64, &a65, &a71, &a73, &a74, &a75, &a76, &d1, &d3, &d4, 
		&d5, &d6, &d7);
    }
    facold = 1e-4;
    expo1 = .2 - *beta * .75;
    facc1 = 1. / *fac1;
    facc2 = 1. / *fac2;
    d__1 = *xend - *x;
    posneg = d_sign(&c_b47, &d__1);
/* --- INITIAL PREPARATIONS */
    atoli = atol[1];
    rtoli = rtol[1];
    last = FALSE_;
    hlamb = 0.;
    iasti = 0;
    (*fcn)(n, x, &y[1], &k1[1], &rpar[1], &ipar[1]);
    *hmax = abs(*hmax);
    iord = 5;
    if (*h__ == 0.) {
	*h__ = hinit_(n, (S_fp)fcn, x, &y[1], xend, &posneg, &k1[1], &k2[1], &
		k3[1], &iord, hmax, &atol[1], &rtol[1], itol, &rpar[1], &ipar[
		1]);
    }
    *nfcn += 2;
    reject = FALSE_;
    condo5_1.xold = *x;
    if (*iout != 0) {
	irtrn = 1;
	condo5_1.hout = *h__;
	i__1 = *naccpt + 1;
	(*solout)(&i__1, &condo5_1.xold, x, &y[1], n, &cont[1], &icomp[1], 
		nrd, &rpar[1], &ipar[1], &irtrn);
	if (irtrn < 0) {
	    goto L79;
	}
    } else {
	irtrn = 0;
    }
/* --- BASIC INTEGRATION STEP */
L1:
    if (*nstep > *nmax) {
	goto L78;
    }
    if (abs(*h__) * .1 <= abs(*x) * *uround) {
	goto L77;
    }
    if ((*x + *h__ * 1.01 - *xend) * posneg > 0.) {
	*h__ = *xend - *x;
	last = TRUE_;
    }
    ++(*nstep);
/* --- THE FIRST 6 STAGES */
    if (irtrn >= 2) {
	(*fcn)(n, x, &y[1], &k1[1], &rpar[1], &ipar[1]);
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L22: */
	y1[i__] = y[i__] + *h__ * a21 * k1[i__];
    }
    d__1 = *x + c2 * *h__;
    (*fcn)(n, &d__1, &y1[1], &k2[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L23: */
	y1[i__] = y[i__] + *h__ * (a31 * k1[i__] + a32 * k2[i__]);
    }
    d__1 = *x + c3 * *h__;
    (*fcn)(n, &d__1, &y1[1], &k3[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L24: */
	y1[i__] = y[i__] + *h__ * (a41 * k1[i__] + a42 * k2[i__] + a43 * k3[
		i__]);
    }
    d__1 = *x + c4 * *h__;
    (*fcn)(n, &d__1, &y1[1], &k4[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	y1[i__] = y[i__] + *h__ * (a51 * k1[i__] + a52 * k2[i__] + a53 * k3[
		i__] + a54 * k4[i__]);
    }
    d__1 = *x + c5 * *h__;
    (*fcn)(n, &d__1, &y1[1], &k5[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L26: */
	ysti[i__] = y[i__] + *h__ * (a61 * k1[i__] + a62 * k2[i__] + a63 * k3[
		i__] + a64 * k4[i__] + a65 * k5[i__]);
    }
    xph = *x + *h__;
    (*fcn)(n, &xph, &ysti[1], &k6[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L27: */
	y1[i__] = y[i__] + *h__ * (a71 * k1[i__] + a73 * k3[i__] + a74 * k4[
		i__] + a75 * k5[i__] + a76 * k6[i__]);
    }
    (*fcn)(n, &xph, &y1[1], &k2[1], &rpar[1], &ipar[1]);
    if (*iout >= 2) {
	i__1 = *nrd;
	for (j = 1; j <= i__1; ++j) {
	    i__ = icomp[j];
	    cont[(*nrd << 2) + j] = *h__ * (d1 * k1[i__] + d3 * k3[i__] + d4 *
		     k4[i__] + d5 * k5[i__] + d6 * k6[i__] + d7 * k2[i__]);
/* L40: */
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L28: */
	k4[i__] = (e1 * k1[i__] + e3 * k3[i__] + e4 * k4[i__] + e5 * k5[i__] 
		+ e6 * k6[i__] + e7 * k2[i__]) * *h__;
    }
    *nfcn += 6;
/* --- ERROR ESTIMATION */
    err = 0.;
    if (*itol == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = y1[i__], abs(
		    d__2));
	    sk = atoli + rtoli * max(d__3,d__4);
/* L41: */
/* Computing 2nd power */
	    d__1 = k4[i__] / sk;
	    err += d__1 * d__1;
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = y1[i__], abs(
		    d__2));
	    sk = atol[i__] + rtol[i__] * max(d__3,d__4);
/* L42: */
/* Computing 2nd power */
	    d__1 = k4[i__] / sk;
	    err += d__1 * d__1;
	}
    }
    err = sqrt(err / *n);
/* --- COMPUTATION OF HNEW */
    fac11 = pow_dd(&err, &expo1);
/* --- LUND-STABILIZATION */
    fac = fac11 / pow_dd(&facold, beta);
/* --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2 */
/* Computing MAX */
/* Computing MIN */
    d__3 = facc1, d__4 = fac / *safe;
    d__1 = facc2, d__2 = min(d__3,d__4);
    fac = max(d__1,d__2);
    hnew = *h__ / fac;
    if (err <= 1.) {
/* --- STEP IS ACCEPTED */
	facold = max(err,1e-4);
	++(*naccpt);
/* ------- STIFFNESS DETECTION */
	if (*naccpt % *nstiff == 0 || iasti > 0) {
	    stnum = 0.;
	    stden = 0.;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = k2[i__] - k6[i__];
		stnum += d__1 * d__1;
/* Computing 2nd power */
		d__1 = y1[i__] - ysti[i__];
		stden += d__1 * d__1;
/* L64: */
	    }
	    if (stden > 0.) {
		hlamb = *h__ * sqrt(stnum / stden);
	    }
	    if (hlamb > 3.25) {
		nonsti = 0;
		++iasti;
		if (iasti == 15) {
		    if (*iprint > 0) {
			io___99.ciunit = *iprint;
			s_wsle(&io___99);
			do_lio(&c__9, &c__1, " THE PROBLEM SEEMS TO BECOME S"
				"TIFF AT X = ", (ftnlen)42);
			do_lio(&c__5, &c__1, (char *)&(*x), (ftnlen)sizeof(
				doublereal));
			e_wsle();
		    }
		    if (*iprint <= 0) {
			goto L76;
		    }
		}
	    } else {
		++nonsti;
		if (nonsti == 6) {
		    iasti = 0;
		}
	    }
	}
	if (*iout >= 2) {
	    i__1 = *nrd;
	    for (j = 1; j <= i__1; ++j) {
		i__ = icomp[j];
		yd0 = y[i__];
		ydiff = y1[i__] - yd0;
		bspl = *h__ * k1[i__] - ydiff;
		cont[j] = y[i__];
		cont[*nrd + j] = ydiff;
		cont[(*nrd << 1) + j] = bspl;
		cont[*nrd * 3 + j] = -(*h__) * k2[i__] + ydiff - bspl;
/* L43: */
	    }
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k1[i__] = k2[i__];
/* L44: */
	    y[i__] = y1[i__];
	}
	condo5_1.xold = *x;
	*x = xph;
	if (*iout != 0) {
	    condo5_1.hout = *h__;
	    i__1 = *naccpt + 1;
	    (*solout)(&i__1, &condo5_1.xold, x, &y[1], n, &cont[1], &icomp[1],
		     nrd, &rpar[1], &ipar[1], &irtrn);
	    if (irtrn < 0) {
		goto L79;
	    }
	}
/* ------- NORMAL EXIT */
	if (last) {
	    *h__ = hnew;
	    *idid = 1;
	    return 0;
	}
	if (abs(hnew) > *hmax) {
	    hnew = posneg * *hmax;
	}
	if (reject) {
/* Computing MIN */
	    d__1 = abs(hnew), d__2 = abs(*h__);
	    hnew = posneg * min(d__1,d__2);
	}
	reject = FALSE_;
    } else {
/* --- STEP IS REJECTED */
/* Computing MIN */
	d__1 = facc1, d__2 = fac11 / *safe;
	hnew = *h__ / min(d__1,d__2);
	reject = TRUE_;
	if (*naccpt >= 1) {
	    ++(*nrejct);
	}
	last = FALSE_;
    }
    *h__ = hnew;
    goto L1;
/* --- FAIL EXIT */
L76:
    *idid = -4;
    return 0;
L77:
    if (*iprint > 0) {
	io___103.ciunit = *iprint;
	s_wsfe(&io___103);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*iprint > 0) {
	io___104.ciunit = *iprint;
	s_wsle(&io___104);
	do_lio(&c__9, &c__1, " STEP SIZE T0O SMALL, H=", (ftnlen)24);
	do_lio(&c__5, &c__1, (char *)&(*h__), (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    *idid = -3;
    return 0;
L78:
    if (*iprint > 0) {
	io___105.ciunit = *iprint;
	s_wsfe(&io___105);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*iprint > 0) {
	io___106.ciunit = *iprint;
	s_wsle(&io___106);
	do_lio(&c__9, &c__1, " MORE THAN NMAX =", (ftnlen)17);
	do_lio(&c__3, &c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, "STEPS ARE NEEDED", (ftnlen)16);
	e_wsle();
    }
    *idid = -2;
    return 0;
L79:
    if (*iprint > 0) {
	io___107.ciunit = *iprint;
	s_wsfe(&io___107);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    *idid = 2;
    return 0;
} /* dopcor_ */


doublereal hinit_(integer *n, S_fp fcn, doublereal *x, doublereal *y, 
	doublereal *xend, doublereal *posneg, doublereal *f0, doublereal *f1, 
	doublereal *y1, integer *iord, doublereal *hmax, doublereal *atol, 
	doublereal *rtol, integer *itol, doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *), pow_dd(
	    doublereal *, doublereal *);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal h1, sk, dnf, dny, der2, der12, atoli, rtoli;

/* ---------------------------------------------------------- */
/* ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS */
/* ---------------------------------------------------------- */
/* ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS */
/* ----   H = 0.01 * NORM (Y0) / NORM (F0) */
/* ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL */
/* ---- COMPARED TO THE SOLUTION */
    /* Parameter adjustments */
    --y1;
    --f1;
    --f0;
    --y;
    --atol;
    --rtol;
    --rpar;
    --ipar;

    /* Function Body */
    dnf = 0.;
    dny = 0.;
    atoli = atol[1];
    rtoli = rtol[1];
    if (*itol == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sk = atoli + rtoli * (d__1 = y[i__], abs(d__1));
/* Computing 2nd power */
	    d__1 = f0[i__] / sk;
	    dnf += d__1 * d__1;
/* L10: */
/* Computing 2nd power */
	    d__1 = y[i__] / sk;
	    dny += d__1 * d__1;
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sk = atol[i__] + rtol[i__] * (d__1 = y[i__], abs(d__1));
/* Computing 2nd power */
	    d__1 = f0[i__] / sk;
	    dnf += d__1 * d__1;
/* L11: */
/* Computing 2nd power */
	    d__1 = y[i__] / sk;
	    dny += d__1 * d__1;
	}
    }
    if (dnf <= 1e-10 || dny <= 1e-10) {
	h__ = 1e-6;
    } else {
	h__ = sqrt(dny / dnf) * .01;
    }
    h__ = min(h__,*hmax);
    h__ = d_sign(&h__, posneg);
/* ---- PERFORM AN EXPLICIT EULER STEP */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L12: */
	y1[i__] = y[i__] + h__ * f0[i__];
    }
    d__1 = *x + h__;
    (*fcn)(n, &d__1, &y1[1], &f1[1], &rpar[1], &ipar[1]);
/* ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION */
    der2 = 0.;
    if (*itol == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sk = atoli + rtoli * (d__1 = y[i__], abs(d__1));
/* L15: */
/* Computing 2nd power */
	    d__1 = (f1[i__] - f0[i__]) / sk;
	    der2 += d__1 * d__1;
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sk = atol[i__] + rtol[i__] * (d__1 = y[i__], abs(d__1));
/* L16: */
/* Computing 2nd power */
	    d__1 = (f1[i__] - f0[i__]) / sk;
	    der2 += d__1 * d__1;
	}
    }
    der2 = sqrt(der2) / h__;
/* ---- STEP SIZE IS COMPUTED SUCH THAT */
/* ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01 */
/* Computing MAX */
    d__1 = abs(der2), d__2 = sqrt(dnf);
    der12 = max(d__1,d__2);
    if (der12 <= 1e-15) {
/* Computing MAX */
	d__1 = 1e-6, d__2 = abs(h__) * .001;
	h1 = max(d__1,d__2);
    } else {
	d__1 = .01 / der12;
	d__2 = 1. / *iord;
	h1 = pow_dd(&d__1, &d__2);
    }
/* Computing MIN */
    d__1 = abs(h__) * 100, d__1 = min(d__1,h1);
    h__ = min(d__1,*hmax);
    ret_val = d_sign(&h__, posneg);
    return ret_val;
} /* hinit_ */


doublereal contd5_(integer *ii, doublereal *x, doublereal *con, integer *
	icomp, integer *nd)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, j;
    static doublereal theta, theta1;

    /* Fortran I/O blocks */
    static cilist io___120 = { 0, 6, 0, 0, 0 };


/* ---------------------------------------------------------- */
/*     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION */
/*     WITH THE OUTPUT-SUBROUTINE FOR DOPRI5. IT PROVIDES AN */
/*     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X. */
/* ---------------------------------------------------------- */
/* ----- COMPUTE PLACE OF II-TH COMPONENT */
    /* Parameter adjustments */
    --icomp;
    --con;

    /* Function Body */
    i__ = 0;
    i__1 = *nd;
    for (j = 1; j <= i__1; ++j) {
	if (icomp[j] == *ii) {
	    i__ = j;
	}
/* L5: */
    }
    if (i__ == 0) {
	s_wsle(&io___120);
	do_lio(&c__9, &c__1, " NO DENSE OUTPUT AVAILABLE FOR COMP.", (ftnlen)
		36);
	do_lio(&c__3, &c__1, (char *)&(*ii), (ftnlen)sizeof(integer));
	e_wsle();
	ret_val = -1.;
	return ret_val;
    }
    theta = (*x - condo5_2.xold) / condo5_2.h__;
    theta1 = 1. - theta;
    ret_val = con[i__] + theta * (con[*nd + i__] + theta1 * (con[(*nd << 1) + 
	    i__] + theta * (con[*nd * 3 + i__] + theta1 * con[(*nd << 2) + 
	    i__])));
    return ret_val;
} /* contd5_ */


/* Subroutine */ int cdopri_(doublereal *c2, doublereal *c3, doublereal *c4, 
	doublereal *c5, doublereal *e1, doublereal *e3, doublereal *e4, 
	doublereal *e5, doublereal *e6, doublereal *e7, doublereal *a21, 
	doublereal *a31, doublereal *a32, doublereal *a41, doublereal *a42, 
	doublereal *a43, doublereal *a51, doublereal *a52, doublereal *a53, 
	doublereal *a54, doublereal *a61, doublereal *a62, doublereal *a63, 
	doublereal *a64, doublereal *a65, doublereal *a71, doublereal *a73, 
	doublereal *a74, doublereal *a75, doublereal *a76, doublereal *d1, 
	doublereal *d3, doublereal *d4, doublereal *d5, doublereal *d6, 
	doublereal *d7)
{
/* ---------------------------------------------------------- */
/*     RUNGE-KUTTA COEFFICIENTS OF DORMAND AND PRINCE (1980) */
/* ---------------------------------------------------------- */
    *c2 = .2;
    *c3 = .3;
    *c4 = .8;
    *c5 = .88888888888888884;
    *a21 = .2;
    *a31 = .074999999999999997;
    *a32 = .22500000000000001;
    *a41 = .97777777777777775;
    *a42 = -3.7333333333333334;
    *a43 = 3.5555555555555554;
    *a51 = 2.9525986892242035;
    *a52 = -11.595793324188385;
    *a53 = 9.8228928516994358;
    *a54 = -.29080932784636487;
    *a61 = 2.8462752525252526;
    *a62 = -10.757575757575758;
    *a63 = 8.9064227177434727;
    *a64 = .27840909090909088;
    *a65 = -.2735313036020583;
    *a71 = .091145833333333329;
    *a73 = .44923629829290207;
    *a74 = .65104166666666663;
    *a75 = -.322376179245283;
    *a76 = .13095238095238096;
    *e1 = .0012326388888888888;
    *e3 = -.0042527702905061394;
    *e4 = .036979166666666667;
    *e5 = -.05086379716981132;
    *e6 = .041904761904761903;
    *e7 = -.025000000000000001;
/* ---- DENSE OUTPUT OF SHAMPINE (1986) */
    *d1 = -1.1270175653862835;
    *d3 = 2.675424484351598;
    *d4 = -5.6855269615885042;
    *d5 = 3.5219323679207912;
    *d6 = -1.7672812570757455;
    *d7 = 2.3824689317781438;
    return 0;
} /* cdopri_ */

