/* dop853.f -- translated by f2c (version 20190311).
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
} condo8_;

#define condo8_1 (condo8_._1)
#define condo8_2 (condo8_._2)

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b46 = 1.;

/* Subroutine */ int dop853_(integer *n, U_fp fcn, doublereal *x, doublereal *
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
    static integer iek1, iek2, iek3, iek4, iek5, iek6, iek7, iek8, iek9, iey1,
	     iek10;
    static doublereal beta, safe;
    static integer ieco, nfcn, meth;
    static doublereal hmax;
    static integer nmax;
    extern /* Subroutine */ int dp86co_(integer *, U_fp, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, U_fp, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer icomp;
    static logical arret;
    static integer nstep, naccpt, nrejct, nstiff, nrdens, iprint, istore;
    static doublereal uround;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 0, 0, 0, 0 };
    static cilist io___10 = { 0, 0, 0, 0, 0 };
    static cilist io___13 = { 0, 0, 0, 0, 0 };
    static cilist io___14 = { 0, 0, 0, 0, 0 };
    static cilist io___17 = { 0, 0, 0, 0, 0 };
    static cilist io___19 = { 0, 0, 0, 0, 0 };
    static cilist io___23 = { 0, 0, 0, 0, 0 };
    static cilist io___39 = { 0, 0, 0, 0, 0 };
    static cilist io___41 = { 0, 0, 0, 0, 0 };


/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER */
/*     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y). */
/*     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER 8(5,3) */
/*     DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND */
/*     DENSE OUTPUT) */

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

/*     Edited (22 Feb 2009) by J.C. Travers: */
/*       renamed HINIT->HINIT853 to avoid name collision with dopri5 */

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
/*                 ATOL SHOULD BE STRICTLY POSITIVE (POSSIBLY VERY SMALL) */

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
/*                    DIMENSION Y(N),CON(8*ND),ICOMP(ND) */
/*                    .... */
/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
/*                    THE FIRST GRID-POINT). */
/*                 "XOLD" IS THE PRECEDING GRID-POINT. */
/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
/*                    IS SET <0, DOP853 WILL RETURN TO THE CALLING PROGRAM. */
/*                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT, */
/*                    SET  IRTRN = 2 */

/*          -----  CONTINUOUS OUTPUT: ----- */
/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
/*                 THE FUNCTION */
/*                        >>>   CONTD8(I,S,CON,ICOMP,ND)   <<< */
/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */

/*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT */
/*                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT */
/*                            (IN THIS CASE WORK(5) MUST BE SPECIFIED) */

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
/*                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE. */
/*                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING. */
/*                 "LWORK" MUST BE AT LEAST  11*N+8*NRDENS+21 */
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
/*              DEFAULT VALUES: WORK(3)=0.333D0, WORK(4)=6.D0 */

/*    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL */
/*              (SEE SECTION IV.2). POSITIVE VALUES OF BETA ( <= 0.04 ) */
/*              MAKE THE STEP SIZE CONTROL MORE STABLE. */
/*              NEGATIVE WORK(5) PROVOKE BETA=0. */
/*              DEFAULT 0.0D0. */

/*    WORK(6)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

/*    WORK(7)   INITIAL STEP SIZE, FOR WORK(7)=0.D0 AN INITIAL GUESS */
/*              IS COMPUTED WITH HELP OF THE FUNCTION HINIT */

/*    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
/*              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000. */

/*    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS */
/*              IF IWORK(2).EQ.1  METHOD DOP853 OF DORMAND AND PRINCE */
/*              (SECTION II.6). */
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
	fac1 = .333;
    } else {
	fac1 = work[3];
    }
    if (work[4] == 0.) {
	fac2 = 6.;
    } else {
	fac2 = work[4];
    }
/* --------- BETA FOR STEP CONTROL STABILIZATION ----------- */
    if (work[5] == 0.) {
	beta = 0.;
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
    iek1 = 21;
    iek2 = iek1 + *n;
    iek3 = iek2 + *n;
    iek4 = iek3 + *n;
    iek5 = iek4 + *n;
    iek6 = iek5 + *n;
    iek7 = iek6 + *n;
    iek8 = iek7 + *n;
    iek9 = iek8 + *n;
    iek10 = iek9 + *n;
    iey1 = iek10 + *n;
    ieco = iey1 + *n;
/* ------ TOTAL STORAGE REQUIREMENT ----------- */
    istore = ieco + (nrdens << 3) - 1;
    if (istore > *lwork) {
	if (iprint > 0) {
	    io___39.ciunit = iprint;
	    s_wsle(&io___39);
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
	    io___41.ciunit = iprint;
	    s_wsle(&io___41);
	    do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR IWORK, MIN. LIWO"
		    "RK=", (ftnlen)45);
	    do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	    e_wsle();
	}
	arret = TRUE_;
    }
/* -------- WHEN A FAIL HAS OCCURRED, WE RETURN WITH IDID=-1 */
    if (arret) {
	*idid = -1;
	return 0;
    }
/* -------- CALL TO CORE INTEGRATOR ------------ */
    dp86co_(n, (U_fp)fcn, x, &y[1], xend, &hmax, &h__, &rtol[1], &atol[1], 
	    itol, &iprint, (U_fp)solout, iout, idid, &nmax, &uround, &meth, &
	    nstiff, &safe, &beta, &fac1, &fac2, &work[iek1], &work[iek2], &
	    work[iek3], &work[iek4], &work[iek5], &work[iek6], &work[iek7], &
	    work[iek8], &work[iek9], &work[iek10], &work[iey1], &work[ieco], &
	    iwork[icomp], &nrdens, &rpar[1], &ipar[1], &nfcn, &nstep, &naccpt,
	     &nrejct);
    work[7] = h__;
    iwork[17] = nfcn;
    iwork[18] = nstep;
    iwork[19] = naccpt;
    iwork[20] = nrejct;
/* ----------- RETURN ----------- */
    return 0;
} /* dop853_ */




/*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

/* Subroutine */ int dp86co_(integer *n, S_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, integer *iprint, S_fp solout, 
	integer *iout, integer *idid, integer *nmax, doublereal *uround, 
	integer *meth, integer *nstiff, doublereal *safe, doublereal *beta, 
	doublereal *fac1, doublereal *fac2, doublereal *k1, doublereal *k2, 
	doublereal *k3, doublereal *k4, doublereal *k5, doublereal *k6, 
	doublereal *k7, doublereal *k8, doublereal *k9, doublereal *k10, 
	doublereal *y1, doublereal *cont, integer *icomp, integer *nrd, 
	doublereal *rpar, integer *ipar, integer *nfcn, integer *nstep, 
	integer *naccpt, integer *nrejct)
{
    /* Format strings */
    static char fmt_979[] = "(\002 EXIT OF DOP853 AT X=\002,e18.4)";

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
    static doublereal xph, err, sk, fac, err2, fac11, deno;
    static integer iord;
    static logical last;
    static doublereal erri, hnew, bspl, facc1, facc2, expo1, hlamb, ydiff, 
	    atoli;
    static integer iasti;
    static doublereal stden, rtoli;
    static integer irtrn;
    static doublereal stnum, facold;
    static logical reject;
    static doublereal posneg;
    static integer nonsti;
    extern doublereal hinit853_(integer *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___68 = { 0, 0, 0, 0, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_979, 0 };
    static cilist io___73 = { 0, 0, 0, 0, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_979, 0 };
    static cilist io___75 = { 0, 0, 0, 0, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_979, 0 };


/* ---------------------------------------------------------- */
/*     CORE INTEGRATOR FOR DOP853 */
/*     PARAMETERS SAME AS IN DOP853 WITH WORKSPACE ADDED */
/* ---------------------------------------------------------- */
/*         DECLARATIONS */
/* ---------------------------------------------------------- */
/* *** *** *** *** *** *** *** */
/*  INITIALISATIONS */
/* *** *** *** *** *** *** *** */
    /* Parameter adjustments */
    --y1;
    --k10;
    --k9;
    --k8;
    --k7;
    --k6;
    --k5;
    --k4;
    --k3;
    --k2;
    --k1;
    --y;
    --rtol;
    --atol;
    --icomp;
    --cont;
    --rpar;
    --ipar;

    /* Function Body */
    facold = 1e-4;
    expo1 = .125 - *beta * .2;
    facc1 = 1. / *fac1;
    facc2 = 1. / *fac2;
    d__1 = *xend - *x;
    posneg = d_sign(&c_b46, &d__1);
/* --- INITIAL PREPARATIONS */
    atoli = atol[1];
    rtoli = rtol[1];
    last = FALSE_;
    hlamb = 0.;
    iasti = 0;
    (*fcn)(n, x, &y[1], &k1[1], &rpar[1], &ipar[1]);
    *hmax = abs(*hmax);
    iord = 8;
    if (*h__ == 0.) {
	*h__ = hinit853_(n, (S_fp)fcn, x, &y[1], xend, &posneg, &k1[1], &k2[1]
		, &k3[1], &iord, hmax, &atol[1], &rtol[1], itol, &rpar[1], &
		ipar[1]);
    }
    *nfcn += 2;
    reject = FALSE_;
    condo8_1.xold = *x;
    if (*iout >= 1) {
	irtrn = 1;
	condo8_1.hout = 1.;
	i__1 = *naccpt + 1;
	(*solout)(&i__1, &condo8_1.xold, x, &y[1], n, &cont[1], &icomp[1], 
		nrd, &rpar[1], &ipar[1], &irtrn);
	if (irtrn < 0) {
	    goto L79;
	}
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
/* --- THE TWELVE STAGES */
    if (irtrn >= 2) {
	(*fcn)(n, x, &y[1], &k1[1], &rpar[1], &ipar[1]);
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L22: */
	y1[i__] = y[i__] + *h__ * .0526001519587677318785587544488 * k1[i__];
    }
    d__1 = *x + *h__ * .0526001519587677318785587544488;
    (*fcn)(n, &d__1, &y1[1], &k2[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L23: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .0197250569845378994544595329183 
		+ k2[i__] * .0591751709536136983633785987549);
    }
    d__1 = *x + *h__ * .0789002279381515978178381316732;
    (*fcn)(n, &d__1, &y1[1], &k3[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L24: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .0295875854768068491816892993775 
		+ k3[i__] * .0887627564304205475450678981324);
    }
    d__1 = *x + *h__ * .11835034190722739672675719751;
    (*fcn)(n, &d__1, &y1[1], &k4[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .241365134159266685502369798665 
		+ k3[i__] * -.884549479328286085344864962717 + k4[i__] * 
		.924834003261792003115737966543);
    }
    d__1 = *x + *h__ * .28164965809277260327324280249;
    (*fcn)(n, &d__1, &y1[1], &k5[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L26: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .037037037037037037037037037037 
		+ k4[i__] * .170828608729473871279604482173 + k5[i__] * 
		.125467687566822425016691814123);
    }
    d__1 = *x + *h__ * .333333333333333333333333333333;
    (*fcn)(n, &d__1, &y1[1], &k6[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L27: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .037109375 + k4[i__] * 
		.170252211019544039314978060272 + k5[i__] * 
		.0602165389804559606850219397283 + k6[i__] * -.017578125);
    }
    d__1 = *x + *h__ * .25;
    (*fcn)(n, &d__1, &y1[1], &k7[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L28: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .0370920001185047927108779319836 
		+ k4[i__] * .170383925712239993810214054705 + k5[i__] * 
		.107262030446373284651809199168 + k6[i__] * 
		-.0153194377486244017527936158236 + k7[i__] * 
		.00827378916381402288758473766002);
    }
    d__1 = *x + *h__ * .307692307692307692307692307692;
    (*fcn)(n, &d__1, &y1[1], &k8[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L29: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .624110958716075717114429577812 
		+ k4[i__] * -3.36089262944694129406857109825 + k5[i__] * 
		-.868219346841726006818189891453 + k6[i__] * 
		27.5920996994467083049415600797 + k7[i__] * 
		20.1540675504778934086186788979 + k8[i__] * 
		-43.4898841810699588477366255144);
    }
    d__1 = *x + *h__ * .651282051282051282051282051282;
    (*fcn)(n, &d__1, &y1[1], &k9[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * .477662536438264365890433908527 
		+ k4[i__] * -2.48811461997166764192642586468 + k5[i__] * 
		-.590290826836842996371446475743 + k6[i__] * 
		21.2300514481811942347288949897 + k7[i__] * 
		15.2792336328824235832596922938 + k8[i__] * 
		-33.2882109689848629194453265587 + k9[i__] * 
		-.0203312017085086261358222928593);
    }
    d__1 = *x + *h__ * .6;
    (*fcn)(n, &d__1, &y1[1], &k10[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L31: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * -.93714243008598732571704021658 
		+ k4[i__] * 5.18637242884406370830023853209 + k5[i__] * 
		1.09143734899672957818500254654 + k6[i__] * 
		-8.14978701074692612513997267357 + k7[i__] * 
		-18.5200656599969598641566180701 + k8[i__] * 
		22.7394870993505042818970056734 + k9[i__] * 
		2.49360555267965238987089396762 + k10[i__] * 
		-3.0467644718982195003823669022);
    }
    d__1 = *x + *h__ * .857142857142857142857142857142;
    (*fcn)(n, &d__1, &y1[1], &k2[1], &rpar[1], &ipar[1]);
    xph = *x + *h__;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L32: */
	y1[i__] = y[i__] + *h__ * (k1[i__] * 2.27331014751653820792359768449 
		+ k4[i__] * -10.5344954667372501984066689879 + k5[i__] * 
		-2.00087205822486249909675718444 + k6[i__] * 
		-17.9589318631187989172765950534 + k7[i__] * 
		27.9488845294199600508499808837 + k8[i__] * 
		-2.85899827713502369474065508674 + k9[i__] * 
		-8.87285693353062954433549289258 + k10[i__] * 
		12.3605671757943030647266201528 + k2[i__] * 
		.643392746015763530355970484046);
    }
    (*fcn)(n, &xph, &y1[1], &k3[1], &rpar[1], &ipar[1]);
    *nfcn += 11;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k4[i__] = k1[i__] * .0542937341165687622380535766363 + k6[i__] * 
		4.45031289275240888144113950566 + k7[i__] * 
		1.89151789931450038304281599044 + k8[i__] * 
		-5.8012039600105847814672114227 + k9[i__] * 
		.31116436695781989440891606237 + k10[i__] * 
		-.152160949662516078556178806805 + k2[i__] * 
		.201365400804030348374776537501 + k3[i__] * 
		.0447106157277725905176885569043;
/* L35: */
	k5[i__] = y[i__] + *h__ * k4[i__];
    }
/* --- ERROR ESTIMATION */
    err = 0.;
    err2 = 0.;
    if (*itol == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = k5[i__], abs(
		    d__2));
	    sk = atoli + rtoli * max(d__3,d__4);
	    erri = k4[i__] - k1[i__] * .244094488188976377952755905512 - k9[
		    i__] * .733846688281611857341361741547 - k3[i__] * 
		    .0220588235294117647058823529412;
/* Computing 2nd power */
	    d__1 = erri / sk;
	    err2 += d__1 * d__1;
	    erri = k1[i__] * .01312004499419488073250102996 + k6[i__] * 
		    -1.225156446376204440720569753 + k7[i__] * 
		    -.4957589496572501915214079952 + k8[i__] * 
		    1.664377182454986536961530415 + k9[i__] * 
		    -.350328848749973681688648729 + k10[i__] * 
		    .3341791187130174790297318841 + k2[i__] * 
		    .08192320648511571246570742613 + k3[i__] * 
		    -.02235530786388629525884427845;
/* L41: */
/* Computing 2nd power */
	    d__1 = erri / sk;
	    err += d__1 * d__1;
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = k5[i__], abs(
		    d__2));
	    sk = atol[i__] + rtol[i__] * max(d__3,d__4);
	    erri = k4[i__] - k1[i__] * .244094488188976377952755905512 - k9[
		    i__] * .733846688281611857341361741547 - k3[i__] * 
		    .0220588235294117647058823529412;
/* Computing 2nd power */
	    d__1 = erri / sk;
	    err2 += d__1 * d__1;
	    erri = k1[i__] * .01312004499419488073250102996 + k6[i__] * 
		    -1.225156446376204440720569753 + k7[i__] * 
		    -.4957589496572501915214079952 + k8[i__] * 
		    1.664377182454986536961530415 + k9[i__] * 
		    -.350328848749973681688648729 + k10[i__] * 
		    .3341791187130174790297318841 + k2[i__] * 
		    .08192320648511571246570742613 + k3[i__] * 
		    -.02235530786388629525884427845;
/* L42: */
/* Computing 2nd power */
	    d__1 = erri / sk;
	    err += d__1 * d__1;
	}
    }
    deno = err + err2 * .01;
    if (deno <= 0.) {
	deno = 1.;
    }
    err = abs(*h__) * err * sqrt(1. / (*n * deno));
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
	(*fcn)(n, &xph, &k5[1], &k4[1], &rpar[1], &ipar[1]);
	++(*nfcn);
/* ------- STIFFNESS DETECTION */
	if (*naccpt % *nstiff == 0 || iasti > 0) {
	    stnum = 0.;
	    stden = 0.;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = k4[i__] - k3[i__];
		stnum += d__1 * d__1;
/* Computing 2nd power */
		d__1 = k5[i__] - y1[i__];
		stden += d__1 * d__1;
/* L64: */
	    }
	    if (stden > 0.) {
		hlamb = abs(*h__) * sqrt(stnum / stden);
	    }
	    if (hlamb > 6.1) {
		nonsti = 0;
		++iasti;
		if (iasti == 15) {
		    if (*iprint > 0) {
			io___68.ciunit = *iprint;
			s_wsle(&io___68);
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
/* ------- FINAL PREPARATION FOR DENSE OUTPUT */
	if (*iout >= 2) {
/* ----    SAVE THE FIRST FUNCTION EVALUATIONS */
	    i__1 = *nrd;
	    for (j = 1; j <= i__1; ++j) {
		i__ = icomp[j];
		cont[j] = y[i__];
		ydiff = k5[i__] - y[i__];
		cont[j + *nrd] = ydiff;
		bspl = *h__ * k1[i__] - ydiff;
		cont[j + (*nrd << 1)] = bspl;
		cont[j + *nrd * 3] = ydiff - *h__ * k4[i__] - bspl;
		cont[j + (*nrd << 2)] = k1[i__] * 
			-8.4289382761090128651353491142 + k6[i__] * 
			.5667149535193777696253178359 + k7[i__] * 
			-3.0689499459498916912797304727 + k8[i__] * 
			2.384667656512069828772814968 + k9[i__] * 
			2.1170345824450282767155149946 + k10[i__] * 
			-.8713915837779729920678990749 + k2[i__] * 
			2.240437430260788275854177165 + k3[i__] * 
			.6315787787694688181557024929;
		cont[j + *nrd * 5] = k1[i__] * 10.427508642579134603413151009 
			+ k6[i__] * 242.28349177525818288430175319 + k7[i__] *
			 165.20045171727028198505394887 + k8[i__] * 
			-374.54675472269020279518312152 + k9[i__] * 
			-22.113666853125306036270938578 + k10[i__] * 
			7.7334326684722638389603898808 + k2[i__] * 
			-30.674084731089398182061213626 + k3[i__] * 
			-9.3321305264302278729567221706;
		cont[j + *nrd * 6] = k1[i__] * 19.985053242002433820987653617 
			+ k6[i__] * -387.03730874935176555105901742 + k7[i__] 
			* -189.17813819516756882830838328 + k8[i__] * 
			527.80815920542364900561016686 + k9[i__] * 
			-11.573902539959630126141871134 + k10[i__] * 
			6.8812326946963000169666922661 + k2[i__] * 
			-1.000605096691083840318386098 + k3[i__] * 
			.7777137798053443209286926574;
		cont[j + *nrd * 7] = k1[i__] * 
			-25.693933462703749003312586129 + k6[i__] * 
			-154.18974869023643374053993627 + k7[i__] * 
			-231.52937917604549567536039109 + k8[i__] * 
			357.6391179106141237828534991 + k9[i__] * 
			93.405324183624310003907691704 + k10[i__] * 
			-37.458323136451633156875139351 + k2[i__] * 
			104.09964950896230045147246184 + k3[i__] * 
			29.840293426660503123344363579;
/* L62: */
	    }
/* ---     THE NEXT THREE FUNCTION EVALUATIONS */
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L51: */
		y1[i__] = y[i__] + *h__ * (k1[i__] * 
			.0561675022830479523392909219681 + k7[i__] * 
			.253500210216624811088794765333 + k8[i__] * 
			-.246239037470802489917441475441 + k9[i__] * 
			-.124191423263816360469010140626 + k10[i__] * 
			.15329179827876569731206322685 + k2[i__] * 
			.00820105229563468988491666602057 + k3[i__] * 
			.00756789766054569976138603589584 + k4[i__] * 
			-.008298);
	    }
	    d__1 = *x + *h__ * .1;
	    (*fcn)(n, &d__1, &y1[1], &k10[1], &rpar[1], &ipar[1]);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L52: */
		y1[i__] = y[i__] + *h__ * (k1[i__] * 
			.0318346481635021405060768473261 + k6[i__] * 
			.0283009096723667755288322961402 + k7[i__] * 
			.0535419883074385676223797384372 + k8[i__] * 
			-.0549237485713909884646569340306 + k2[i__] * 
			-1.08347328697249322858509316994e-4 + k3[i__] * 
			3.82571090835658412954920192323e-4 + k4[i__] * 
			-3.40465008687404560802977114492e-4 + k10[i__] * 
			.141312443674632500278074618366);
	    }
	    d__1 = *x + *h__ * .2;
	    (*fcn)(n, &d__1, &y1[1], &k2[1], &rpar[1], &ipar[1]);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L53: */
		y1[i__] = y[i__] + *h__ * (k1[i__] * 
			-.428896301583791923408573538692 + k6[i__] * 
			-4.69762141536116384314449447206 + k7[i__] * 
			7.68342119606259904184240953878 + k8[i__] * 
			4.06898981839711007970213554331 + k9[i__] * 
			.356727187455281109270669543021 + k4[i__] * 
			-.00139902416515901462129418009734 + k10[i__] * 
			2.9475147891527723389556272149 + k2[i__] * 
			-9.15095847217987001081870187138);
	    }
	    d__1 = *x + *h__ * .777777777777777777777777777778;
	    (*fcn)(n, &d__1, &y1[1], &k3[1], &rpar[1], &ipar[1]);
	    *nfcn += 3;
/* ---     FINAL PREPARATION */
	    i__1 = *nrd;
	    for (j = 1; j <= i__1; ++j) {
		i__ = icomp[j];
		cont[j + (*nrd << 2)] = *h__ * (cont[j + (*nrd << 2)] + k4[
			i__] * -.0889903364513333108206981174 + k10[i__] * 
			18.148505520854727256656404962 + k2[i__] * 
			-9.1946323924783554000451984436 + k3[i__] * 
			-4.4360363875948939664310572);
		cont[j + *nrd * 5] = *h__ * (cont[j + *nrd * 5] + k4[i__] * 
			15.697238121770843886131091075 + k10[i__] * 
			-31.139403219565177677282850411 + k2[i__] * 
			-9.3529243588444783865713862664 + k3[i__] * 
			35.81684148639408375246589854);
		cont[j + *nrd * 6] = *h__ * (cont[j + *nrd * 6] + k4[i__] * 
			-2.7782057523535084065932004339 + k10[i__] * 
			-60.196695231264120758267380846 + k2[i__] * 
			84.320405506677161018159903784 + k3[i__] * 
			11.99229113618278932803513003);
		cont[j + *nrd * 7] = *h__ * (cont[j + *nrd * 7] + k4[i__] * 
			-43.533456590011143754432175058 + k10[i__] * 
			96.3245539591882829483949506 + k2[i__] * 
			-39.177261675615439165231486172 + k3[i__] * 
			-149.72683625798562581422125276);
/* L63: */
	    }
	    condo8_1.hout = *h__;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k1[i__] = k4[i__];
/* L67: */
	    y[i__] = k5[i__];
	}
	condo8_1.xold = *x;
	*x = xph;
	if (*iout >= 1) {
	    i__1 = *naccpt + 1;
	    (*solout)(&i__1, &condo8_1.xold, x, &y[1], n, &cont[1], &icomp[1],
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
	io___72.ciunit = *iprint;
	s_wsfe(&io___72);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*iprint > 0) {
	io___73.ciunit = *iprint;
	s_wsle(&io___73);
	do_lio(&c__9, &c__1, " STEP SIZE TOO SMALL, H=", (ftnlen)24);
	do_lio(&c__5, &c__1, (char *)&(*h__), (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    *idid = -3;
    return 0;
L78:
    if (*iprint > 0) {
	io___74.ciunit = *iprint;
	s_wsfe(&io___74);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*iprint > 0) {
	io___75.ciunit = *iprint;
	s_wsle(&io___75);
	do_lio(&c__9, &c__1, " MORE THAN NMAX =", (ftnlen)17);
	do_lio(&c__3, &c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, "STEPS ARE NEEDED", (ftnlen)16);
	e_wsle();
    }
    *idid = -2;
    return 0;
L79:
    if (*iprint > 0) {
	io___76.ciunit = *iprint;
	s_wsfe(&io___76);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    *idid = 2;
    return 0;
} /* dp86co_ */


doublereal hinit853_(integer *n, S_fp fcn, doublereal *x, doublereal *y, 
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
} /* hinit853_ */


doublereal contd8_(integer *ii, doublereal *x, doublereal *con, integer *
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
    static doublereal s, s1, conpar;

    /* Fortran I/O blocks */
    static cilist io___89 = { 0, 6, 0, 0, 0 };


/* ---------------------------------------------------------- */
/*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONNECTION */
/*     WITH THE OUTPUT-SUBROUTINE FOR DOP853. IT PROVIDES AN */
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
	s_wsle(&io___89);
	do_lio(&c__9, &c__1, " NO DENSE OUTPUT AVAILABLE FOR COMP.", (ftnlen)
		36);
	do_lio(&c__3, &c__1, (char *)&(*ii), (ftnlen)sizeof(integer));
	e_wsle();
	ret_val = -1.;
	return ret_val;
    }
    s = (*x - condo8_2.xold) / condo8_2.h__;
    s1 = 1. - s;
    conpar = con[i__ + (*nd << 2)] + s * (con[i__ + *nd * 5] + s1 * (con[i__ 
	    + *nd * 6] + s * con[i__ + *nd * 7]));
    ret_val = con[i__] + s * (con[i__ + *nd] + s1 * (con[i__ + (*nd << 1)] + 
	    s * (con[i__ + *nd * 3] + s1 * conpar)));
    return ret_val;
} /* contd8_ */

