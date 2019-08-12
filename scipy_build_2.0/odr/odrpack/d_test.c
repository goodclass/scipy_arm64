/* d_test.f -- translated by f2c (version 20190311).
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
    integer ntest;
} tstset_;

#define tstset_1 tstset_

struct {
    integer setno;
} setid_;

#define setid_1 setid_

/* Table of constant values */

static integer c__1 = 1;
static integer c__4380 = 4380;

/* DTEST */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);

    /* Local variables */
    extern /* Subroutine */ int dodrx_(doublereal *, logical *, integer *);
    static logical passed;
    static doublereal tstfac;
    static integer lunerr, lunsum, lunrpt;

/* ***BEGIN PROLOGUE  TEST */
/* ***REFER TO DODR,DODRC */
/* ***ROUTINES CALLED  DODRX */
/* ***DATE WRITTEN   861229   (YYMMDD) */
/* ***REVISION DATE  920619   (YYMMDD) */
/* ***PURPOSE  EXERCISE FEATURES OF ODRPACK SOFTWARE */
/* ***END PROLOGUE  ODRPACK */
/* ...SCALARS IN COMMON */
/* ...LOCAL SCALARS */
/* ...EXTERNAL SUBROUTINES */
/* ...COMMON BLOCKS */
/* ***VARIABLE DECLARATIONS (ALPHABETICALLY) */
/*   LUNERR:  THE LOGICAL UNIT NUMBER USED FOR ERROR MESSAGES. */
/*   LUNRPT:  THE LOGICAL UNIT NUMBER USED FOR COMPUTATION REPORTS. */
/*   LUNSUM:  THE LOGICAL UNIT NUMBER USED FOR A SUMMARY REPORT LISTING */
/*            ONLY THE TEST COMPARISONS AND NOT THE ODRPACK GENERATED */
/*            REPORTS. */
/*   NTEST:   THE NUMBER OF TESTS TO BE RUN. */
/*   PASSED:  THE VARIABLE DESIGNATING WHETHER THE RESULTS OF ALL OF THE */
/*            TESTS AGREE WITH THOSE FROM THE CRAY YMP USING DOUBLE */
/*            PRECISION (PASSED=TRUE), OR WHETHER SOME OF THE RESULTS */
/*            DISAGREED (PASSED=FALSE). */
/*   TSTFAC:  THE USER-SUPPLIED FACTOR FOR SCALING THE TEST TOLERANCES */
/*            USED TO CHECK FOR AGREEMENT BETWEEN COMPUTED RESULTS AND */
/*            RESULTS OBTAINED USING DOUBLE PRECISION VERSION ON CRAY */
/*            YMP.  VALUES OF TSTFAC GREATER THAN ONE INCREASE THE */
/*            TEST TOLERANCES, MAKING THE TESTS EASIER TO PASS AND */
/*            ALLOWING SMALL DISCREPANCIES BETWEEN THE COMPUTED AND */
/*            EXPECTED RESULTS TO BE AUTOMATICALLY DISCOUNTED. */
/* ***FIRST EXECUTABLE STATEMENT  TEST */
/*  SET UP NECESSARY FILES */
/*  NOTE:  ODRPACK GENERATES COMPUTATION AND ERROR REPORTS ON */
/*         LOGICAL UNIT 6 BY DEFAULT; */
/*         LOGICAL UNIT 'LUNSUM' USED TO SUMMARIZE RESULTS OF COMPARISONS */
/*         FROM EXERCISE ROUTINE DODRX. */
    lunrpt = 18;
    lunerr = 18;
    lunsum = 19;
    o__1.oerr = 0;
    o__1.ounit = lunrpt;
    o__1.ofnmlen = 6;
    o__1.ofnm = "REPORT";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = lunerr;
    o__1.ofnmlen = 6;
    o__1.ofnm = "REPORT";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = lunsum;
    o__1.ofnmlen = 7;
    o__1.ofnm = "SUMMARY";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*  EXERCISE DOUBLE PRECISION VERSION OF ODRPACK */
/*  (TEST REPORTS GENERATED ON FILE 'RESULTS' AND */
/*   SUMMARIZED IN FILE 'SUMMARY') */
    tstset_1.ntest = 12;
    tstfac = 1.;
    dodrx_(&tstfac, &passed, &lunsum);
    return 0;
} /* MAIN__ */

/* DODRX */
/* Subroutine */ int dodrx_(doublereal *tstfac, logical *passed, integer *
	lunsum)
{
    /* Initialized data */

    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal hundrd = 100.;
    static doublereal dpymp[24]	/* was [2][12] */ = { 
	    27627.33195780256808978449342964,
	    7.532639569022918943695104672512e-4,
	    27627.32630143673024399942947263,
	    7.53846772268713150687427931494e-4,
	    1069944100.000000027940905194068,
	    1.212808593256056359629660672046e-5,
	    1069944100.000000026623461142867,
	    5.452084633790606017572015067556e-7,
	    1.426988156377258617521571734503,1.084728687127432219753903919409,
	    4.261321829513978871872508874025,.0147796721039842073356542432928,
	    4.261272307142888464011486769858,
	    .01477966125465374336804138554559,
	    43.71487317909745009110272283622,
	    .00114441947440828606711223359255,
	    3.099048849376848610380977303924,
	    .08824708863783850023783338218501,
	    9.469917836739932584221023234527,.420538921558810465119853680988,
	    39.50949253027682207109233363651,66.51838750834910819636881506915,
	    39.50949253027682207109233363651,66.51838750834910819636881506915 
	    };
    static integer idpymp[12] = { 1,1,3,1,1,4,1,1,2,1,1023,40100 };
    static doublereal zero = 0.;
    static doublereal p01 = .01;
    static doublereal p2 = .2;
    static doublereal one = 1.;

    /* Format strings */
    static char fmt_1000[] = "(\0021\002)";
    static char fmt_1001[] = "(\002 EXAMPLE \002,i2/)";
    static char fmt_1010[] = "(\002 TEST SIMPLE ODR PROBLEM\002/\002 WITH AN"
	    "ALYTIC DERIVATIVES\002,\002 USING DODR.\002)";
    static char fmt_1020[] = "(\002 TEST SIMPLE OLS PROBLEM\002/\002 WITH FI"
	    "NITE DIFFERENCE DERIVATIVES\002,\002 USING DODR.\002)";
    static char fmt_1030[] = "(\002 TEST PARAMETER FIXING CAPABILITIES\002"
	    ",\002 FOR POORLY SCALED OLS PROBLEM\002/\002 WITH ANALYTIC DERIV"
	    "ATIVES\002,\002 USING DODRC.\002)";
    static char fmt_1040[] = "(\002 TEST WEIGHTING CAPABILITIES\002,\002 FOR"
	    " ODR PROBLEM\002/\002 WITH ANALYTIC DERIVATIVES\002,\002 USING D"
	    "ODRC. \002/\002 ALSO SHOWS SOLUTION OF POORLY SCALED\002,\002 OD"
	    "R PROBLEM.\002/\002 (DERIVATIVE CHECKING TURNED OFF.)\002)";
    static char fmt_1050[] = "(\002 TEST DELTA INITIALIZATION CAPABILITIE"
	    "S\002/\002 AND USE OF ISTOP TO RESTRICT PARAMETER VALUES\002,"
	    "\002 FOR ODR PROBLEM\002/\002 WITH ANALYTIC DERIVATIVES\002,\002"
	    " USING DODRC.\002)";
    static char fmt_1060[] = "(\002 TEST STIFF STOPPING CONDITIONS\002,\002 "
	    "FOR UNSCALED ODR PROBLEM\002/\002 WITH ANALYTIC DERIVATIVES\002"
	    ",\002 USING DODRC.\002)";
    static char fmt_1070[] = "(\002 TEST RESTART\002,\002 FOR UNSCALED ODR P"
	    "ROBLEM\002/\002 WITH ANALYTIC DERIVATIVES\002,\002 USING DODRC"
	    ".\002)";
    static char fmt_1080[] = "(\002 TEST USE OF TAUFAC TO RESTRICT FIRST S"
	    "TEP\002,\002 FOR ODR PROBLEM\002/\002 WITH FINITE DIFFERENCE DER"
	    "IVATIVES\002,\002 USING DODRC.\002)";
    static char fmt_1090[] = "(\002 TEST IMPLICIT MODEL\002,\002 FOR OLS PRO"
	    "BLEM\002/\002 USING DODRC.\002)";
    static char fmt_1100[] = "(\002 TEST MULTIRESPONSE MODEL\002,\002 FOR OD"
	    "R PROBLEM\002/\002 WITH FINITE DIFFERENCE DERIVATIVES\002,\002 U"
	    "SING DODRC.\002)";
    static char fmt_1110[] = "(\002 TEST DETECTION OF QUESTIONABLE ANALYTIC "
	    "DERIVATIVES\002,\002 FOR OLS PROBLEM\002/\002 USING DODRC.\002)";
    static char fmt_1120[] = "(\002 TEST DETECTION OF INCORRECT ANALYTIC DER"
	    "IVATIVES\002,\002 FOR ODR PROBLEM\002/\002 WITH ANALYTIC DERIVAT"
	    "IVES\002,\002 USING DODRC.\002)";
    static char fmt_2200[] = "(\002 DATA SET REFERENCE: \002,a80)";
    static char fmt_3100[] = "(/\002 COMPARISON OF NEW RESULTS WITH\002,\002"
	    " DOUBLE PRECISION CRAY YMP RESULT:\002//\002                    "
	    "     NORM OF BETA\002,\002        SUM OF SQUARED WTD OBS ERRORS "
	    " INFO\002)";
    static char fmt_3210[] = "(/a25/1p,2d37.30,i6)";
    static char fmt_3220[] = "(/a25,1p,d12.5,25x,d12.5,i6)";
    static char fmt_3310[] = "(/\002 *** STOPPING CONDITIONS\002,\002 SHOW C"
	    "ONVERGENCE NOT ATTAINED. ***\002/\002        NO FURTHER COMPARIS"
	    "ONS MADE BETWEEN RESULTS.\002//)";
    static char fmt_3320[] = "(//\002 *** WARNING ***\002,\002 RESULTS DO NO"
	    "T AGREE TO WITHIN STOPPING TOLERANCE. ***\002//)";
    static char fmt_3330[] = "(//\002 *** RESULTS AGREE TO WITHIN STOPPING T"
	    "OLERANCE. ***\002//)";
    static char fmt_3340[] = "(//\002 *** WARNING ***\002,\002 STOPPING COND"
	    "ITIONS DO NOT AGREE. ***\002//)";
    static char fmt_3350[] = "(//\002 *** WARNING ***\002,\002 UNEXPECTED ST"
	    "OPPING CONDITION.\002,\002  PLEASE CONTACT PACKAGE AUTHORS. **"
	    "*\002//)";
    static char fmt_4100[] = "(///\002 *** SUMMARY:\002,\002 ONE OR MORE TES"
	    "TS DO NOT AGREE WITH EXPECTED RESULTS. ***\002)";
    static char fmt_4200[] = "(///\002 *** SUMMARY:\002,\002 ALL TESTS AGREE"
	    " WITH EXPECTED RESULTS. ***\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, l, m, n;
    static doublereal x[150]	/* was [50][3] */, y[100]	/* was [50][2]
	     */, wd[450]	/* was [50][3][3] */, we[200]	/* was [50][2]
	    [2] */;
    static integer np, nq, job, msg, ldx, ldy, lun;
    static doublereal wrk[250], wss, beta[10], sclb[10], scld[150]	/* 
	    was [50][3] */;
    extern /* Subroutine */ int dodr_(U_fp, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static integer info;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal bnrm, stpb[10], stpd[150]	/* was [50][3] */, work[4380];
    static integer ldwd1, ldwe1;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dodrc_(U_fp, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *);
    static integer ld2wd1, ld2we1;
    static logical fails;
    static integer ifixb[10], ldifx;
    extern /* Subroutine */ int dwght_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical isodr;
    static integer maxit;
    static char title[80];
    static integer lwmin, ifixx[150]	/* was [50][3] */, itest, iwork[56];
    extern /* Subroutine */ int dzero_(integer *, integer *, doublereal *, 
	    integer *);
    static logical short__;
    static doublereal sstol;
    static logical failed;
    static doublereal taufac;
    static integer ldscld;
    static doublereal epsmac;
    extern doublereal dmprec_(void);
    static integer ndigit;
    extern /* Subroutine */ int dodrxd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dodrxf_();
    static integer ldstpd, liwmin;
    static doublereal partol, wssdel;
    static integer iprint, lunerr;
    extern /* Subroutine */ int dodrxw_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, logical *, integer *, integer *);
    static integer lunrpt;
    static doublereal wsseps, tsttol;

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_1030, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_1040, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_1050, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_1060, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_1070, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_1080, 0 };
    static cilist io___82 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_1090, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_1110, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_1120, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_2200, 0 };
    static cilist io___97 = { 0, 0, 0, fmt_2200, 0 };
    static cilist io___108 = { 0, 0, 0, fmt_3100, 0 };
    static cilist io___109 = { 0, 0, 0, fmt_3210, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_3210, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_3220, 0 };
    static cilist io___112 = { 0, 0, 0, fmt_3220, 0 };
    static cilist io___113 = { 0, 0, 0, fmt_3310, 0 };
    static cilist io___114 = { 0, 0, 0, fmt_3320, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_3330, 0 };
    static cilist io___116 = { 0, 0, 0, fmt_3340, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_3350, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___119 = { 0, 0, 0, fmt_4100, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_4100, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_4200, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_4200, 0 };


/* ***BEGIN PROLOGUE  DODRX */
/* ***REFER TO DODR,DODRC */
/* ***ROUTINES CALLED  DDOT,DMPREC,DNRM2,DODR,DODRC,DODRXD, */
/*                    DODRXF,DODRXW,DWGHT,DZERO */
/* ***DATE WRITTEN   860529   (YYMMDD) */
/* ***REVISION DATE  920619   (YYMMDD) */
/* ***PURPOSE  EXERCISE FEATURES OF ODRPACK SOFTWARE */
/* ***END PROLOGUE  DODRX */
/* ...PARAMETERS */
/* ...SCALAR ARGUMENTS */
/* ...SCALARS IN COMMON */
/* ...LOCAL SCALARS */
/* ...LOCAL ARRAYS */
/* ...EXTERNAL FUNCTIONS */
/* ...EXTERNAL SUBROUTINES */
/* ...INTRINSIC FUNCTIONS */
/* ...COMMON BLOCKS */
/* ...DATA STATEMENTS */
/* ...ROUTINE NAMES USED AS SUBPROGRAM ARGUMENTS */
/*   DODRXF:  THE USER-SUPPLIED ROUTINE FOR EVALUATING THE MODEL. */
/* ...VARIABLE DEFINITIONS (ALPHABETICALLY) */
/*   BETA:    THE FUNCTION PARAMETERS. */
/*   BNRM:    THE NORM OF BETA. */
/*   DPYMP:   THE FLOATING POINT RESULTS FROM A CRAY YMP USING */
/*            DOUBLE PRECISION. */
/*   EPSMAC:  THE VALUE OF MACHINE PRECISION. */
/*   FAILED:  THE VARIABLE DESIGNATING WHETHER THE RESULTS OF ALL OF THE */
/*            DEMONSTRATION RUNS AGREED WITH THOSE FROM THE CRAY YMP */
/*            USING DOUBLE PRECISION (FAILED=FALSE) OR WHETHER SOME OF */
/*            THE TESTS DISAGREED (FAILED=TRUE). */
/*   FAILS:   THE VARIABLE DESIGNATING WHETHER THE RESULTS OF AN */
/*            INDIVIDUAL DEMONSTRATION RUN AGREED WITH THOSE FROM THE */
/*            CRAY YMP USING DOUBLE PRECISION (FAILS=FALSE) OR */
/*            DISAGREE (FAILS=TRUE). */
/*   HUNDRD:  THE VALUE 100.0D0. */
/*   I:       AN INDEX VARIABLE. */
/*   IDPYMP:  THE INTEGER RESULTS FROM A CRAY YMP USING */
/*            DOUBLE PRECISION. */
/*   IFIXB:   THE VALUES DESIGNATING WHETHER THE ELEMENTS OF BETA ARE */
/*            FIXED AT THEIR INPUT VALUES OR NOT. */
/*   IFIXX:   THE VALUES DESIGNATING WHETHER THE ELEMENTS OF DELTA ARE */
/*            FIXED AT THEIR INPUT VALUES OR NOT. */
/*   INFO:    THE VARIABLE DESIGNATING WHY THE COMPUTATIONS STOPPED. */
/*   IPRINT:  THE PRINT CONTROL VARIABLE. */
/*   ISODR:   THE VARIABLE DESIGNATING WHETHER THE SOLUTION IS BY ODR */
/*            (ISODR=TRUE) OR BY OLS (ISODR=FALSE). */
/*   ITEST:   THE NUMBER OF THE CURRENT TEST BEING RUN. */
/*   IWORK:   THE INTEGER WORK SPACE. */
/*   J:       AN INDEX VARIABLE. */
/*   JOB:     THE VARIABLE CONTROLLING PROBLEM INITIALIZATION AND */
/*            COMPUTATIONAL METHOD. */
/*   LDIFX:   THE LEADING DIMENSION OF ARRAY IFIXX. */
/*   LDSCLD:  THE LEADING DIMENSION OF ARRAY SCLD. */
/*   LDWD:    THE LEADING DIMENSION OF ARRAY WD. */
/*   LDWD1:   THE LEADING DIMENSION OF ARRAY WD AS PASSED TO ODRPACK. */
/*   LDWE:    THE LEADING DIMENSION OF ARRAY WE. */
/*   LDWE1:   THE LEADING DIMENSION OF ARRAY WE AS PASSED TO ODRPACK. */
/*   LDX:     THE LEADING DIMENSION OF ARRAY X. */
/*   LDY:     THE LEADING DIMENSION OF ARRAY Y. */
/*   LD2WD:   THE SECOND DIMENSION OF ARRAY WD. */
/*   LD2WD1:  THE SECOND DIMENSION OF ARRAY WD AS PASSED TO ODRPACK. */
/*   LD2WE:   THE SECOND DIMENSION OF ARRAY WE. */
/*   LD2WE1:  THE SECOND DIMENSION OF ARRAY WE AS PASSED TO ODRPACK. */
/*   LIWKMN:  THE MINIMUM ACCEPTABLE LENGTH OF ARRAY IWORK. */
/*   LIWMIN:  THE MINIMUM LENGTH OF VECTOR IWORK FOR A GIVEN PROBLEM. */
/*   LIWORK:  THE LENGTH OF VECTOR IWORK. */
/*   LUN:     THE LOGICAL UNIT NUMBER CURRENTLY BEING USED. */
/*   LUNERR:  THE LOGICAL UNIT NUMBER USED FOR ERROR MESSAGES. */
/*   LUNRPT:  THE LOGICAL UNIT NUMBER USED FOR COMPUTATION REPORTS. */
/*   LUNSUM:  THE LOGICAL UNIT NUMBER USED FOR A SUMMARY REPORT. */
/*   LWKMN:   THE MINIMUM ACCEPTABLE LENGTH OF ARRAY WORK. */
/*   LWMIN:   THE MINIMUM LENGTH OF VECTOR WORK FOR A GIVEN PROBLEM. */
/*   LWORK:   THE LENGTH OF VECTOR WORK. */
/*   M:       THE NUMBER OF COLUMNS OF DATA IN THE EXPLANATORY VARIABLE. */
/*   MAXIT:   THE MAXIMUM NUMBER OF ITERATIONS ALLOWED. */
/*   MSG:     THE VARIABLE DESIGNATING WHICH MESSAGE IS TO BE PRINTED AS */
/*            A RESULT OF THE COMPARISON WITH THE CRAY YMP RESULTS. */
/*   N:       THE NUMBER OF OBSERVATIONS. */
/*   NDIGIT:  THE NUMBER OF ACCURATE DIGITS IN THE FUNCTION RESULTS, AS */
/*            SUPPLIED BY THE USER. */
/*   NP:      THE NUMBER OF FUNCTION PARAMETERS. */
/*   NTEST:   THE NUMBER OF TESTS TO BE RUN. */
/*   NTESTS:  THE NUMBER OF DIFFERENT TESTS AVAILABLE. */
/*   ONE:     THE VALUE 1.0D0. */
/*   PASSED:  THE VARIABLE DESIGNATING WHETHER THE RESULTS OF ALL OF THE */
/*            DEMONSTRATION RUNS AGREED WITH THOSE FROM THE CRAY YMP */
/*            USING DOUBLE PRECISION (PASSED=TRUE), OR WHETHER SOME OF */
/*            THE RESULTS DISAGREED (PASSED=FALSE). */
/*   P01:     THE VALUE 0.01D0. */
/*   P2:      THE VALUE 0.2D0. */
/*   PARTOL:  THE PARAMETER CONVERGENCE STOPPING CRITERIA. */
/*   SCLB:    THE SCALING VALUES FOR BETA. */
/*   SCLD:    THE SCALING VALUES FOR DELTA. */
/*   SETNO:   THE NUMBER OF THE DATA SET BEING ANALYZED. */
/*   SHORT:   THE VARIABLE DESIGNATING WHETHER ODRPACK IS INVOKED BY THE */
/*            SHORT-CALL (SHORT=.TRUE.) OR THE LONG-CALL (SHORT=.FALSE.). */
/*   SSTOL:   THE SUM-OF-SQUARES CONVERGENCE STOPPING TOLERANCE. */
/*   TAUFAC:  THE FACTOR USED TO COMPUTE THE INITIAL TRUST REGION */
/*            DIAMETER. */
/*   THREE:   THE VALUE 3.0D0. */
/*   TITLE:   THE REFERENCE FOR THE DATA SET BEING ANALYZED. */
/*   TSTFAC:  THE USER-SUPPLIED FACTOR FOR SCALING THE TEST TOLERANCES */
/*            USED TO CHECK FOR AGREEMENT BETWEEN COMPUTED RESULTS AND */
/*            RESULTS OBTAINED USING DOUBLE PRECISION VERSION ON CRAY */
/*            YMP. */
/*   TSTTOL:  THE TEST TOLERANCE USED IN CHECKING COMPUTED VALUES FOR */
/*            PURPOSES OF DETERMINING PROPER INSTALLATION. */
/*   TWO:     THE VALUE 2.0D0. */
/*   WD:      THE DELTA WEIGHTS. */
/*   WE:      THE EPSILON WEIGHTS. */
/*   WORK:    THE DOUBLE PRECISION WORK SPACE. */
/*   WRK:     THE DOUBLE PRECISION WORK SPACE FOR COMPUTING TEST RESULTS. */
/*   WSS:     THE SUM OF THE SQUARED WEIGHTED ERRORS. */
/*   WSSDEL:  THE SUM OF THE SQUARED WEIGHTED ERRORS IN X. */
/*   WSSEPS:  THE SUM OF THE SQUARED WEIGHTED ERRORS IN Y. */
/*   X:       THE EXPLANATORY VARIABLE. */
/*   Y:       THE RESPONSE VARIABLE. */
/*   ZERO:    THE VALUE 0.0D0. */
/* ***FIRST EXECUTABLE STATEMENT  DODRX */
/*  SET LOGICAL UNITS FOR ERROR AND COMPUTATION REPORTS */
    lunerr = 18;
    lunrpt = 18;
/*  INITIALIZE TEST TOLERANCE */
    if (*tstfac > one) {
	tsttol = *tstfac;
    } else {
	tsttol = one;
    }
/*  INITIALIZE MACHINE PRECISION */
    epsmac = dmprec_();
/*  INITIALIZE LEADING DIMENSION OF X */
    ldx = 50;
    ldy = 50;
/*  INITIALIZE MISCELLANEOUS VARIABLES USED IN THE EXERCISE PROCEDURE */
    failed = FALSE_;
    short__ = TRUE_;
    isodr = TRUE_;
    n = 0;
/*  BEGIN EXERCISING ODRPACK */
    i__1 = tstset_1.ntest;
    for (itest = 1; itest <= i__1; ++itest) {
/*  SET CONTROL VALUES TO INVOKE DEFAULT VALUES */
	we[0] = -one;
	ldwe1 = 50;
	ld2we1 = 2;
	wd[0] = -one;
	ldwd1 = 50;
	ld2wd1 = 3;
	ifixb[0] = -1;
	ifixx[0] = -1;
	ldifx = 50;
	ndigit = -1;
	taufac = -one;
	sstol = -one;
	partol = -one;
	maxit = -1;
	iprint = 2112;
	stpb[0] = -one;
	stpd[0] = -one;
	ldstpd = 1;
	sclb[0] = -one;
	scld[0] = -one;
	ldscld = 1;
	if (itest == 1) {
/*  TEST SIMPLE ODR PROBLEM */
/*  WITH ANALYTIC DERIVATIVES. */
	    lun = lunrpt;
	    io___48.ciunit = lun;
	    s_wsfe(&io___48);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___50.ciunit = lun;
		s_wsfe(&io___50);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___51.ciunit = lun;
		s_wsfe(&io___51);
		e_wsfe();
		lun = *lunsum;
/* L10: */
	    }
	    setid_1.setno = 5;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 20;
	    short__ = TRUE_;
	    isodr = TRUE_;
	} else if (itest == 2) {
/*  TEST SIMPLE OLS PROBLEM */
/*  WITH FORWARD DIFFERENCE DERIVATIVES. */
	    lun = lunrpt;
	    io___61.ciunit = lun;
	    s_wsfe(&io___61);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___62.ciunit = lun;
		s_wsfe(&io___62);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___63.ciunit = lun;
		s_wsfe(&io___63);
		e_wsfe();
		lun = *lunsum;
/* L20: */
	    }
	    setid_1.setno = 5;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 2;
	    short__ = TRUE_;
	    isodr = FALSE_;
	} else if (itest == 3) {
/*  TEST PARAMETER FIXING CAPABILITIES FOR POORLY SCALED OLS PROBLEM */
/*  WITH ANALYTIC DERIVATIVES. */
/*  (DERIVATIVE CHECKING TURNED OFF.) */
	    lun = lunrpt;
	    io___64.ciunit = lun;
	    s_wsfe(&io___64);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___65.ciunit = lun;
		s_wsfe(&io___65);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___66.ciunit = lun;
		s_wsfe(&io___66);
		e_wsfe();
		lun = *lunsum;
/* L30: */
	    }
	    setid_1.setno = 3;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    ifixb[0] = 1;
	    ifixb[1] = 1;
	    ifixb[2] = 1;
	    ifixb[3] = 0;
	    ifixb[4] = 1;
	    ifixb[5] = 0;
	    ifixb[6] = 0;
	    ifixb[7] = 0;
	    ifixb[8] = 0;
	    job = 42;
	    short__ = FALSE_;
	    isodr = FALSE_;
	} else if (itest == 4) {
/*  TEST WEIGHTING CAPABILITIES FOR ODR PROBLEM WITH */
/*  ANALYTIC DERIVATIVES. */
/*  ALSO SHOWS SOLUTION OF POORLY SCALED ODR PROBLEM. */
/*  (DERIVATIVE CHECKING TURNED OFF.) */
/*  N.B., THIS RUN CONTINUES FROM WHERE TEST 3 LEFT OFF. */
	    lun = lunrpt;
	    io___67.ciunit = lun;
	    s_wsfe(&io___67);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___68.ciunit = lun;
		s_wsfe(&io___68);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___69.ciunit = lun;
		s_wsfe(&io___69);
		e_wsfe();
		lun = *lunsum;
/* L40: */
	    }
	    setid_1.setno = 3;
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    ldwd1 = 50;
	    ldwe1 = 50;
	    ld2wd1 = 3;
	    ld2we1 = 2;
	    i__2 = n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
		d__2 = p01 / (d__1 = x[i__ - 1], abs(d__1));
		wd[i__ - 1] = d__2 * d__2;
		we[i__ - 1] = one;
/* L45: */
	    }
	    we[27] = zero;
	    ifixb[0] = 1;
	    ifixb[1] = 1;
	    ifixb[2] = 1;
	    ifixb[3] = 0;
	    ifixb[4] = 1;
	    ifixb[5] = 1;
	    ifixb[6] = 1;
	    ifixb[7] = 0;
	    ifixb[8] = 0;
	    job = 30;
	    iprint = 2232;
	    short__ = FALSE_;
	    isodr = TRUE_;
	} else if (itest == 5) {
/*  TEST DELTA INITIALIZATION CAPABILITIES AND USER-SUPPLIED SCALING */
/*  AND USE OF ISTOP TO RESTRICT PARAMETER VALUES */
/*  FOR ODR PROBLEM WITH ANALYTIC DERIVATIVES. */
	    lun = lunrpt;
	    io___70.ciunit = lun;
	    s_wsfe(&io___70);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___71.ciunit = lun;
		s_wsfe(&io___71);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___72.ciunit = lun;
		s_wsfe(&io___72);
		e_wsfe();
		lun = *lunsum;
/* L50: */
	    }
	    setid_1.setno = 1;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 1020;
	    ldscld = 1;
	    scld[0] = two;
	    sclb[0] = p2;
	    sclb[1] = one;
	    ldwe1 = 1;
	    ld2we1 = 1;
	    we[0] = -one;
	    ldwd1 = 1;
	    ld2wd1 = 1;
	    wd[0] = -one;
	    for (i__ = 20; i__ <= 21; ++i__) {
		work[i__ - 1] = beta[0] / y[i__ - 1] + beta[1] - x[i__ - 1];
/* L55: */
	    }
	    short__ = FALSE_;
	    isodr = TRUE_;
	} else if (itest == 6) {
/*  TEST STIFF STOPPING CONDITIONS FOR UNSCALED ODR PROBLEM */
/*  WITH ANALYTIC DERIVATIVES. */
	    lun = lunrpt;
	    io___73.ciunit = lun;
	    s_wsfe(&io___73);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___74.ciunit = lun;
		s_wsfe(&io___74);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___75.ciunit = lun;
		s_wsfe(&io___75);
		e_wsfe();
		lun = *lunsum;
/* L60: */
	    }
	    setid_1.setno = 4;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 20;
	    sstol = hundrd * epsmac;
	    partol = epsmac;
	    maxit = 2;
	    short__ = FALSE_;
	    isodr = TRUE_;
	} else if (itest == 7) {
/*  TEST RESTART FOR UNSCALED ODR PROBLEM */
/*  WITH ANALYTIC DERIVATIVES. */
	    lun = lunrpt;
	    io___76.ciunit = lun;
	    s_wsfe(&io___76);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___77.ciunit = lun;
		s_wsfe(&io___77);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___78.ciunit = lun;
		s_wsfe(&io___78);
		e_wsfe();
		lun = *lunsum;
/* L70: */
	    }
	    setid_1.setno = 4;
	    job = 20220;
	    sstol = hundrd * epsmac;
	    partol = epsmac;
	    maxit = 50;
	    short__ = FALSE_;
	    isodr = TRUE_;
	} else if (itest == 8) {
/*  TEST USE OF TAUFAC TO RESTRICT FIRST STEP */
/*  FOR ODR PROBLEM WITH CENTRAL DIFFERENCE DERIVATIVES. */
	    lun = lunrpt;
	    io___79.ciunit = lun;
	    s_wsfe(&io___79);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___80.ciunit = lun;
		s_wsfe(&io___80);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___81.ciunit = lun;
		s_wsfe(&io___81);
		e_wsfe();
		lun = *lunsum;
/* L80: */
	    }
	    setid_1.setno = 6;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 210;
	    taufac = p01;
	    short__ = FALSE_;
	    isodr = TRUE_;
	} else if (itest == 9) {
/*  TEST IMPLICIT ODR PROBLEM */
/*  WITH FORWARD FINITE DIFFERENCE DERIVATIVES */
/*  AND COVARIANCE MATRIX CONSTRUCTED WITH RECOMPUTED DERIVATIVES. */
	    lun = lunrpt;
	    io___82.ciunit = lun;
	    s_wsfe(&io___82);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___83.ciunit = lun;
		s_wsfe(&io___83);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___84.ciunit = lun;
		s_wsfe(&io___84);
		e_wsfe();
		lun = *lunsum;
/* L90: */
	    }
	    setid_1.setno = 7;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 1;
	    d__1 = one / three;
	    partol = pow_dd(&epsmac, &d__1);
	    short__ = TRUE_;
	    isodr = TRUE_;
	} else if (itest == 10) {
/*  TEST MULTIRESPONSE ODR PROBLEM */
/*  WITH CENTRAL DIFFERENCE DERIVATIVES , */
/*  DELTA INITIALIZED TO NONZERO VALUES, */
/*  VARIABLE FIXING,  AND WEIGHTING. */
	    lun = lunrpt;
	    io___85.ciunit = lun;
	    s_wsfe(&io___85);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___86.ciunit = lun;
		s_wsfe(&io___86);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___87.ciunit = lun;
		s_wsfe(&io___87);
		e_wsfe();
		lun = *lunsum;
/* L100: */
	    }
	    setid_1.setno = 8;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    ldwd1 = 50;
	    ldwe1 = 50;
	    ld2wd1 = 3;
	    ld2we1 = 2;
	    i__2 = n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*  INITIALIZE DELTA, AND SPECIFY FIRST DECADE OF FREQUENCIES AS FIXED */
		if (x[i__ - 1] < 100.) {
		    work[i__ - 1] = 0.;
		    ifixx[i__ - 1] = 0;
		} else if (x[i__ - 1] <= 150.) {
		    work[i__ - 1] = 0.;
		    ifixx[i__ - 1] = 1;
		} else if (x[i__ - 1] <= 1e3) {
		    work[i__ - 1] = 25.;
		    ifixx[i__ - 1] = 1;
		} else if (x[i__ - 1] <= 1e4) {
		    work[i__ - 1] = 560.;
		    ifixx[i__ - 1] = 1;
		} else if (x[i__ - 1] <= 1e5) {
		    work[i__ - 1] = 9500.;
		    ifixx[i__ - 1] = 1;
		} else {
		    work[i__ - 1] = 1.44e5;
		    ifixx[i__ - 1] = 1;
		}
/*  SET WEIGHTS */
		if (x[i__ - 1] == 100. || x[i__ - 1] == 150.) {
		    we[i__ - 1] = 0.;
		    we[i__ + 99] = 0.;
		    we[i__ + 49] = 0.;
		    we[i__ + 149] = 0.;
		} else {
		    we[i__ - 1] = 559.6;
		    we[i__ + 99] = -1634.;
		    we[i__ + 49] = -1634.;
		    we[i__ + 149] = 8397.;
		}
/* Computing 2nd power */
		d__1 = x[i__ - 1];
		wd[i__ - 1] = 1e-4 / (d__1 * d__1);
/* L105: */
	    }
	    job = 210;
	    short__ = FALSE_;
	    isodr = TRUE_;
	} else if (itest == 11) {
/*  TEST DETECTION OF INCORRECT DERIVATIVES */
	    lun = lunrpt;
	    io___88.ciunit = lun;
	    s_wsfe(&io___88);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___89.ciunit = lun;
		s_wsfe(&io___89);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___90.ciunit = lun;
		s_wsfe(&io___90);
		e_wsfe();
		lun = *lunsum;
/* L110: */
	    }
	    setid_1.setno = 6;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 22;
	    short__ = FALSE_;
	    isodr = FALSE_;
	} else if (itest == 12) {
/*  TEST DETECTION OF INCORRECT DERIVATIVES */
	    lun = lunrpt;
	    io___91.ciunit = lun;
	    s_wsfe(&io___91);
	    e_wsfe();
	    for (i__ = 1; i__ <= 2; ++i__) {
		io___92.ciunit = lun;
		s_wsfe(&io___92);
		do_fio(&c__1, (char *)&itest, (ftnlen)sizeof(integer));
		e_wsfe();
		io___93.ciunit = lun;
		s_wsfe(&io___93);
		e_wsfe();
		lun = *lunsum;
/* L120: */
	    }
	    setid_1.setno = 6;
	    dodrxd_(title, &n, &m, &np, &nq, &ldx, x, &ldy, y, beta, (ftnlen)
		    80);
	    dzero_(&c__4380, &c__1, work, &c__4380);
	    job = 20;
	    short__ = FALSE_;
	    isodr = TRUE_;
	}
	dodrxw_(&n, &m, &np, &nq, &ldwe1, &ld2we1, &isodr, &liwmin, &lwmin);
/*  COMPUTE SOLUTION */
	io___96.ciunit = lunrpt;
	s_wsfe(&io___96);
	do_fio(&c__1, title, (ftnlen)80);
	e_wsfe();
	io___97.ciunit = *lunsum;
	s_wsfe(&io___97);
	do_fio(&c__1, title, (ftnlen)80);
	e_wsfe();
	if (short__) {
	    dodr_((U_fp)dodrxf_, &n, &m, &np, &nq, beta, y, &ldy, x, &ldx, we,
		     &ldwe1, &ld2we1, wd, &ldwd1, &ld2wd1, &job, &iprint, &
		    lunerr, &lunrpt, work, &lwmin, iwork, &liwmin, &info);
	} else {
	    dodrc_((U_fp)dodrxf_, &n, &m, &np, &nq, beta, y, &ldy, x, &ldx, 
		    we, &ldwe1, &ld2we1, wd, &ldwd1, &ld2wd1, ifixb, ifixx, &
		    ldifx, &job, &ndigit, &taufac, &sstol, &partol, &maxit, &
		    iprint, &lunerr, &lunrpt, stpb, stpd, &ldstpd, sclb, scld,
		     &ldscld, work, &lwmin, iwork, &liwmin, &info);
	}
/*  COMPARE RESULTS WITH THOSE OBTAINED ON THE CRAY YMP */
/*  USING DOUBLE PRECISION VERSION OF ODRPACK */
	bnrm = dnrm2_(&np, beta, &c__1);
	dwght_(&n, &m, wd, &ldwd1, &ld2wd1, work, &n, wrk, &n);
	i__2 = n * m;
	wssdel = ddot_(&i__2, work, &c__1, wrk, &c__1);
	dwght_(&n, &nq, we, &ldwe1, &ld2we1, &work[n * m], &n, &wrk[n * m], &
		n);
	i__2 = n * nq;
	wsseps = ddot_(&i__2, &work[n * m], &c__1, &wrk[n * m], &c__1);
	wss = wsseps + wssdel;
	if (sstol < zero) {
	    sstol = sqrt(epsmac);
	} else {
	    sstol = min(sstol,one);
	}
	if (partol < zero) {
	    d__1 = two / three;
	    partol = pow_dd(&epsmac, &d__1);
	} else {
	    partol = min(partol,one);
	}
	if (info >= 10000) {
	    if (idpymp[itest - 1] == info) {
		fails = FALSE_;
		msg = 1;
	    } else {
		fails = TRUE_;
		msg = 3;
	    }
	} else if (info % 10 == 1) {
	    fails = (d__1 = wss - dpymp[(itest << 1) - 1], abs(d__1)) > dpymp[
		    (itest << 1) - 1] * sstol * tsttol;
	    msg = 2;
	} else if (info % 10 == 2) {
	    fails = (d__1 = bnrm - dpymp[(itest << 1) - 2], abs(d__1)) > 
		    dpymp[(itest << 1) - 2] * partol * tsttol;
	    msg = 2;
	} else if (info % 10 == 3) {
	    fails = (d__1 = wss - dpymp[(itest << 1) - 1], abs(d__1)) > dpymp[
		    (itest << 1) - 1] * sstol * tsttol && (d__2 = bnrm - 
		    dpymp[(itest << 1) - 2], abs(d__2)) > dpymp[(itest << 1) 
		    - 2] * partol * tsttol;
	    msg = 2;
	} else if (info % 10 == 4 && idpymp[itest - 1] == 4) {
	    fails = FALSE_;
	    msg = 1;
	} else if (info == idpymp[itest - 1]) {
	    fails = TRUE_;
	    msg = 4;
	} else {
	    fails = TRUE_;
	    msg = 3;
	}
	failed = failed || fails;
	lun = lunrpt;
	for (l = 1; l <= 2; ++l) {
	    io___108.ciunit = lun;
	    s_wsfe(&io___108);
	    e_wsfe();
	    io___109.ciunit = lun;
	    s_wsfe(&io___109);
	    do_fio(&c__1, " CRAY YMP RESULT = ", (ftnlen)19);
	    do_fio(&c__1, (char *)&dpymp[(itest << 1) - 2], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&dpymp[(itest << 1) - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&idpymp[itest - 1], (ftnlen)sizeof(integer))
		    ;
	    e_wsfe();
	    io___110.ciunit = lun;
	    s_wsfe(&io___110);
	    do_fio(&c__1, " NEW TEST RESULT      = ", (ftnlen)24);
	    do_fio(&c__1, (char *)&bnrm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&wss, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&info, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___111.ciunit = lun;
	    s_wsfe(&io___111);
	    do_fio(&c__1, " DIFFERENCE           = ", (ftnlen)24);
	    d__2 = (d__1 = dpymp[(itest << 1) - 2] - bnrm, abs(d__1));
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    d__4 = (d__3 = dpymp[(itest << 1) - 1] - wss, abs(d__3));
	    do_fio(&c__1, (char *)&d__4, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___112.ciunit = lun;
	    s_wsfe(&io___112);
	    do_fio(&c__1, " RELATIVE ERROR       = ", (ftnlen)24);
	    d__5 = (d__1 = dpymp[(itest << 1) - 2] - bnrm, abs(d__1)) / (d__2 
		    = dpymp[(itest << 1) - 2], abs(d__2));
	    do_fio(&c__1, (char *)&d__5, (ftnlen)sizeof(doublereal));
	    d__6 = (d__3 = dpymp[(itest << 1) - 1] - wss, abs(d__3)) / (d__4 =
		     dpymp[(itest << 1) - 1], abs(d__4));
	    do_fio(&c__1, (char *)&d__6, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    if (msg == 1) {
		io___113.ciunit = lun;
		s_wsfe(&io___113);
		e_wsfe();
	    } else if (msg == 2) {
		if (fails) {
		    io___114.ciunit = lun;
		    s_wsfe(&io___114);
		    e_wsfe();
		} else {
		    io___115.ciunit = lun;
		    s_wsfe(&io___115);
		    e_wsfe();
		}
	    } else if (msg == 3) {
		io___116.ciunit = lun;
		s_wsfe(&io___116);
		e_wsfe();
	    } else if (msg == 4) {
		io___117.ciunit = lun;
		s_wsfe(&io___117);
		e_wsfe();
	    }
	    lun = *lunsum;
/* L300: */
	}
/* L400: */
    }
    io___118.ciunit = lunrpt;
    s_wsfe(&io___118);
    e_wsfe();
    if (failed) {
	io___119.ciunit = lunrpt;
	s_wsfe(&io___119);
	e_wsfe();
	io___120.ciunit = *lunsum;
	s_wsfe(&io___120);
	e_wsfe();
	*passed = FALSE_;
    } else {
	io___121.ciunit = lunrpt;
	s_wsfe(&io___121);
	e_wsfe();
	io___122.ciunit = *lunsum;
	s_wsfe(&io___122);
	e_wsfe();
	*passed = TRUE_;
    }
/*  FORMAT STATEMENTS */
    return 0;
} /* dodrx_ */

/* DODRXD */
/* Subroutine */ int dodrxd_(char *title, integer *n, integer *m, integer *np,
	 integer *nq, integer *ldx, doublereal *x, integer *ldy, doublereal *
	y, doublereal *beta, ftnlen title_len)
{
    /* Initialized data */

    static char tdata[80*10] = " BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 1  "
	    "                                    " " BOGGS, BYRD AND SCHNABEL"
	    ", 1985, EXAMPLE 2                                      " " BOGGS"
	    ", BYRD AND SCHNABEL, 1985, EXAMPLE 3                            "
	    "          " " HIMMELBLAU, 1970, EXAMPLE 6.2-4, PAGE 188         "
	    "                             " " DRAPER AND SMITH, 1981, EXERCIS"
	    "E I, PAGE 521-522                               " " POWELL AND M"
	    "ACDONALD, 1972, TABLES 7 AND 8, PAGES 153-154                   "
	    "   " " FULLER, 1987, TABLE 3.2.10, PAGES 244-245                "
	    "                      " " BATES AND WATTS, 1988, TABLE A1.13, PA"
	    "GES 280-281                              ";
    static integer ndata[10] = { 40,50,44,13,8,14,20,23 };
    static integer mdata[10] = { 1,2,1,2,2,1,2,1 };
    static integer npdata[10] = { 2,3,9,3,2,3,5,5 };
    static integer nqdata[10] = { 1,1,1,1,1,1,1,2 };
    static struct {
	doublereal e_1[2];
	doublereal fill_2[8];
	doublereal e_3[3];
	doublereal fill_4[7];
	doublereal e_5[9];
	doublereal fill_6[1];
	doublereal e_7[3];
	doublereal fill_8[7];
	doublereal e_9[2];
	doublereal fill_10[8];
	doublereal e_11[3];
	doublereal fill_12[7];
	doublereal e_13[5];
	doublereal fill_14[5];
	doublereal e_15[5];
	doublereal fill_16[25];
	} equiv_135 = { 1., 1., {0}, -1., 1., 1., {0}, 2.81887509408440189e-6,
		 -.00231290549212363845, 5.83035555572801965, 0., 
		40691077.6203121026, .00138001105225, .0596038513209999999, 
		6.70582099359999998, 1069944100., {0}, 3., 3., -.5, {0}, 
		.01155, 5e3, {0}, 25., 30., 6., {0}, -1., -3., .09, .02, .08, 
		{0}, 4., 2., 7., .4, .5 };

#define bdata ((doublereal *)&equiv_135)

    static struct {
	doublereal e_1[40];
	doublereal fill_2[110];
	doublereal e_3[50];
	doublereal fill_4[100];
	doublereal e_5[44];
	doublereal fill_6[106];
	doublereal e_7[13];
	doublereal fill_8[137];
	doublereal e_9[8];
	doublereal fill_10[142];
	doublereal e_11[14];
	doublereal fill_12[136];
	doublereal e_13[20];
	doublereal fill_14[130];
	doublereal e_15[23];
	doublereal fill_16[27];
	doublereal e_17[23];
	doublereal fill_18[377];
	} equiv_136 = { -1.19569795672791172, -1.28023349509594288, 
		-1.25270693343174591, -.996698267935287383, 
		-1.04681033065801934, -1.46724952092847308, 
		-1.23366891873487528, -1.65665097907185554, 
		-1.68476460930907119, -1.98571971169224491, 
		-1.95691696638051344, -2.11871342665769836, 
		-2.6864293255867102, -2.81123260058024347, 
		-3.2870448658178592, -4.23062993461887032, 
		-5.12043906552226903, -7.31032616379005535, 
		-10.9002759485608993, -25.1810238510370206, 
		100.123028650879944, 16.8225085871915048, 8.94830510866913009,
		 6.45853815227747004, 4.98218564760117328, 
		3.82971664718710476, 3.44116492497344184, 2.76840496973858949,
		 2.59521665196956666, 2.05996022794557661, 
		1.97939614345337836, 1.56739340562905589, 1.59032057073028366,
		 1.73102268158937949, 1.55512561664824758, 1.4963599494413326,
		 1.47487601463073568, 1.17244575233306998, .91093133606917258,
		 1.26172980914513272, {0}, .6808327772179429, 
		1.221835945953022, 1.189586787346082, 1.469826237640946, 
		1.677753381893553, 2.024857219060262, 2.589128519359388, 
		3.668942032541548, 5.746095833513473, 12.76764240264893, 
		1.234730796936231, 1.422561208640828, 1.698895340130247, 
		1.734855779012044, 2.777612639728346, 3.391633246626173, 
		5.896151373121475, 12.44156252145768, -49.84917391538616, 
		-8.327955090006186, 1.849346177742399, 1.751929791768392, 
		2.539493812385358, 3.735007749285017, 5.48408128950331, 
		12.52568805217743, -49.35877971649166, -8.011589749654127, 
		-4.373994870619341, -2.978001034259796, 2.718110574546613, 
		3.770358656133924, 5.601110539171431, 12.81523761749268, 
		-49.87091777324672, -8.157976969083143, -4.402404911951586, 
		-2.767239570617675, -2.232036672887348, -1.69728270310622, 
		5.51015652153227, 12.80361804962158, -49.8257683396339, 
		-8.773345502217619, -4.538201921568676, -2.974993157386779, 
		-2.127432559785389, -2.09703205365401, -1.552872920420862, 
		-1.613566737704807, {0}, .988227696721327788, 
		.988268083998559958, .988341022958438831, .988380557606306446,
		 .988275062411751338, .988326680176446987, 
		.988306058860433439, .988292880079125555, .988305279259496905,
		 .988278142019574202, .988224953369819946, 
		.988111989169778223, .988045627103840613, .987913715667047655,
		 .987841994238525678, .98763845043243427, .987587364331771395,
		 .987576264149633684, .987539209110983643, 
		.987621143807705698, .988023229785526217, .988558376710994197,
		 .989304775352439885, .990210452265710472, .9910959505922639, 
		.991475677297119272, .991901306250746771, .992619222425303263,
		 .993617037631973475, .994727321698030676, 
		.996523114720326189, .99803690956376402, .999151968626971372, 
		1.00017083706131769, 1.00110046382923523, 1.00059103180404652,
		 .999211829791257561, .994711451526761862, 
		.989844132928847109, .987234104554490439, .980928240178404887,
		 .970888680366055576, .960043769857327398, 
		.947277159259551068, {0}, 2.93, 1.95, .81, .58, 5.9, 4.74, 
		4.18, 4.05, 9.03, 7.85, 7.22, 8.5, 9.81, {0}, .912, .382, 
		.397, .376, .342, .358, .348, .376, {0}, 26.38, 25.79, 25.29, 
		24.86, 24.46, 24.1, 23.78, 23.5, 23.24, 23., 22.78, 22.58, 
		22.39, 22.22, {0}, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
		 0., 0., 0., 0., 0., 0., 0., 0., 0., {0}, 4.22, 4.167, 4.132, 
		4.038, 4.019, 3.956, 3.884, 3.784, 3.713, 3.633, 3.54, 3.433, 
		3.358, 3.258, 3.193, 3.128, 3.059, 2.984, 2.934, 2.876, 2.838,
		 2.798, 2.759, {0}, .136, .167, .188, .212, .236, .257, .276, 
		.297, .309, .311, .314, .311, .305, .289, .277, .255, .24, 
		.218, .202, .182, .168, .153, .139 };

#define ydata ((doublereal *)&equiv_136)

    static struct {
	doublereal e_1[40];
	doublereal fill_2[110];
	doublereal e_3[100];
	doublereal fill_4[50];
	doublereal e_5[44];
	doublereal fill_6[106];
	doublereal e_7[13];
	doublereal fill_8[37];
	doublereal e_9[13];
	doublereal fill_10[87];
	doublereal e_11[8];
	doublereal fill_12[42];
	doublereal e_13[8];
	doublereal fill_14[92];
	doublereal e_15[14];
	doublereal fill_16[136];
	doublereal e_17[20];
	doublereal fill_18[30];
	doublereal e_19[20];
	doublereal fill_20[80];
	doublereal e_21[23];
	doublereal fill_22[427];
	} equiv_137 = { -.0213701920211315155, .0494813247025012969, 
		.127889194935560226, .128615394085645676, .232544285655021667,
		 .268151108026504516, .309041029810905456, 
		.405991539210081099, .376611424833536147, .475875890851020811,
		 .49924693539738655, .536615037024021147, .58183076590299606, 
		.684512710422277446, .660219819694757458, .766990323960781092,
		 .808270426690578456, .897410020083189004, 
		.959199774116277687, .914675474762916558, .997759691476821892,
		 1.07136870384216308, 1.08033321037888526, 
		1.16064198672771453, 1.19080889359116553, 1.2941887518763542, 
		1.35594148099422453, 1.35302808716893195, 1.37994666010141371,
		 1.47630019545555113, 1.5345070807635784, 1.52805351451039313,
		 1.57147316247224806, 1.66649596005678175, 
		1.66505665838718412, 1.75214128553867338, 1.80567992463707922,
		 1.84624404296278952, 1.95568727388978002, 
		1.99326394036412237, {0}, .06254745988339948, 
		.2025003436206424, .1649437385998765, .3048741376105061, 
		.5327274455806651, .5088237075989102, .704227041878554, 
		.592077736111512, 1.049409456464216, .9793825175586192, 
		.06378704531655387, .1761233129060257, .310965082300263, 
		.3113942691167821, .4470761261906125, .3847862309982111, 
		.6490931764507805, .6856120053725255, .9687471394250884, 
		.8697893679895329, -.004653099303327366, .00604753397196646, 
		.239418809621756, .4566624689116998, .3711153205220795, 
		.586442107042503, .579796274973298, .8050080949038999, 
		.63724234083571, .9821328179361187, -.02235156571212627, 
		.1360814275450336, .1453670530198706, .3082219195764355, 
		.4326587691335283, .4777855010799803, .727986827616619, 
		.7459503855882651, .7325375035271135, .9673523614338463, 
		.01297617843108911, .1701632439506297, .162768461906274, 
		.2229148079461658, .4029100956046249, .2337708125934432, 
		.6465286934869147, .8028116585689694, .8371378598912229, 
		1.031659807565266, .1101790642097831, -.01961408628913276, 
		.1665148747509966, .006129086880414905, .09382487875524446, 
		.004996057750205054, .08193548490923262, .01271139606723891, 
		.02580952436583161, .1242807551810279, .3048564011371964, 
		.2623870280788969, .2264307654747588, .2713758404102818, 
		.2550008589026183, .154958003178364, .2583016854637732, 
		.1073912606032286, .1519325261357407, .06255075005864, 
		.5467956625953752, .2309057494739227, .1907520696811707, 
		.3288706151709844, .4399785566406605, .4906890437522867, 
		.5218609982033831, .2922835389553916, .402261740352486, 
		.392546836419047, .6504790197089788, .7530201018976618, 
		.6111535320030931, .4552172832904239, .678607663414113, 
		.536178207572157, .6684979205734939, .7860775890072637, 
		.5826251640468284, .4607793960168328, .70000953793186, 
		.8531318307643487, .865315129048175, .7975117585020945, 
		.7614929587270231, .8960000958442235, .9685743337007557, 
		.9048664504767116, .8356844249900219, .7939021919123461, {0}, 
		2.5e-9, 6.4e-9, 1e-8, 9e-8, 1e-6, 4e-6, 9e-6, 1.6e-5, 3.6e-5, 
		6.4e-5, 1e-4, 1.44e-4, 2.25e-4, 4e-4, 6.25e-4, 9e-4, .001225, 
		.0016, .002025, .0025, .0036, .0049, .0064, .0081, .01, 
		.011025, .0121, .0144, .0169, .0196, .0256, .0324, .04, 
		.050625, .075625, .1225, .16, .25, .3364, .3844, .49, .64, 
		.81, 1., {0}, 0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.5,
		 2.9, {0}, 0., 1., 2., 3., 0., 1., 2., 2., 0., 1., 2., 2., 
		1.8, {0}, 109., 65., 1180., 66., 1270., 69., 1230., 68., {0}, 
		600., 640., 600., 640., 600., 640., 600., 640., {0}, 1., 2., 
		3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., {0}, .5, 
		1.2, 1.6, 1.86, 2.12, 2.36, 2.44, 2.36, 2.06, 1.74, 1.34, .9, 
		-.28, -.78, -1.36, -1.9, -2.5, -2.88, -3.18, -3.44, {0}, -.12,
		 -.6, -1., -1.4, -2.54, -3.36, -4., -4.75, -5.25, -5.64, 
		-5.97, -6.32, -6.44, -6.44, -6.41, -6.25, -5.88, -5.5, -5.24, 
		-4.86, {0}, 30., 50., 70., 100., 150., 200., 300., 500., 700.,
		 1e3, 1500., 2e3, 3e3, 5e3, 7e3, 1e4, 1.5e4, 2e4, 3e4, 5e4, 
		7e4, 1e5, 1.5e5 };

#define xdata ((doublereal *)&equiv_137)


    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, l;

/* ***BEGIN PROLOGUE  DODRXD */
/* ***REFER TO  DODR,DODRC */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   860529   (YYMMDD) */
/* ***REVISION DATE  920619   (YYMMDD) */
/* ***PURPOSE  SET UP DATA FOR ODRPACK EXERCISER */
/* ***END PROLOGUE  DODRXD */
/* ...PARAMETERS */
/* ...SCALAR ARGUMENTS */
/* ...ARRAY ARGUMENTS */
/* ...SCALARS IN COMMON */
/* ...LOCAL SCALARS */
/* ...LOCAL ARRAYS */
/* ...COMMON BLOCKS */
/* ...DATA STATEMENTS */
    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --beta;

    /* Function Body */
/* ...VARIABLE DEFINITIONS (ALPHABETICALLY) */
/*   BDATA:   THE FUNCTION PARAMETER FOR EACH DATA SET. */
/*   BETA:    THE FUNCTION PARAMETERS. */
/*   I:       AN INDEXING VARIABLE. */
/*   J:       AN INDEXING VARIABLE. */
/*   L:       AN INDEXING VARIABLE. */
/*   LDX:     THE LEADING DIMENSION OF ARRAY X. */
/*   M:       THE NUMBER OF COLUMNS OF DATA IN THE EXPLANATORY VARIABLE. */
/*   MDATA:   THE NUMBER OF COLUMNS OF DATA IN THE EXPLANATORY VARIABLE */
/*            IN EACH DATA SET. */
/*   N:       THE NUMBER OF OBSERVATIONS. */
/*   NDATA:   THE NUMBER OF OBSERVATIONS PER DATA SET. */
/*   NP:      THE NUMBER OF FUNCTION PARAMETERS. */
/*   NPDATA:  THE NUMBER OF FUNCTION PARAMETERS IN EACH DATA SET. */
/*   NQDATA:  THE NUMBER OF RESPONSES PER OBSERVATION IN EACH DATA SET. */
/*   SETNO:   THE NUMBER OF THE DATA SET BEING ANALYZED. */
/*   TDATA:   THE REFERENCE FOR THE EACH OF THE DATA SETS. */
/*   TITLE:   THE REFERENCE FOR THE DATA SET BEING ANALYZED. */
/*   X:       THE EXPLANATORY VARIABLES. */
/*   XDATA:   THE EXPLANATORY VARIABLES FOR EACH DATA SET. */
/*   Y:       THE RESPONSE VARIABLE. */
/*   YDATA:   THE RESPONSE VARIABLES FOR EACH DATA SET. */
/* ***FIRST EXECUTABLE STATEMENT  DODRXD */
    s_copy(title, tdata + (0 + (0 + (setid_1.setno - 1) * 80)), (ftnlen)80, (
	    ftnlen)80);
    *n = ndata[setid_1.setno - 1];
    *m = mdata[setid_1.setno - 1];
    *np = npdata[setid_1.setno - 1];
    *nq = nqdata[setid_1.setno - 1];
    i__1 = *nq;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + l * y_dim1] = ydata[i__ + (l + setid_1.setno * 3) * 50 - 
		    201];
/* L10: */
	}
/* L20: */
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ + j * x_dim1] = xdata[i__ + (j + setid_1.setno * 3) * 50 - 
		    201];
/* L30: */
	}
/* L40: */
    }
    i__1 = *np;
    for (k = 1; k <= i__1; ++k) {
	beta[k] = bdata[k + setid_1.setno * 10 - 11];
/* L50: */
    }
    return 0;
} /* dodrxd_ */

#undef xdata
#undef ydata
#undef bdata


/* DODRXF */
/* Subroutine */ int dodrxf_(integer *n, integer *m, integer *np, integer *nq,
	 integer *ldn, integer *ldm, integer *ldnp, doublereal *beta, 
	doublereal *xplusd, integer *ifixb, integer *ifixx, integer *ldifx, 
	integer *ideval, doublereal *f, doublereal *fjacb, doublereal *fjacd, 
	integer *istop)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer f_dim1, f_offset, fjacb_dim1, fjacb_dim2, fjacb_offset, 
	    fjacd_dim1, fjacd_dim2, fjacd_offset, xplusd_dim1, xplusd_offset, 
	    ifixx_dim1, ifixx_offset, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double exp(doublereal), pow_dd(doublereal *, doublereal *), cos(
	    doublereal), sin(doublereal), atan2(doublereal, doublereal), sqrt(
	    doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, pi, phi, fac1, fac2, fac3, fac4, freq, omega, 
	    theta, ctheta, stheta;

/* ***BEGIN PROLOGUE  DODRXF */
/* ***REFER TO  DODR,DODRC */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   860529   (YYMMDD) */
/* ***REVISION DATE  920619   (YYMMDD) */
/* ***PURPOSE  COMPUTE JACOBIAN MATRICIES FOR ODRPACK EXERCISER */
/* ***END PROLOGUE  DODRXF */
/* ...SCALAR ARGUMENTS */
/* ...ARRAY ARGUMENTS */
/* ...SCALARS IN COMMON */
/* ...LOCAL SCALARS */
/* ...INTRINSIC FUNCTIONS */
/* ...COMMON BLOCKS */
/* ...DATA STATEMENTS */
    /* Parameter adjustments */
    --ifixb;
    --beta;
    f_dim1 = *ldn;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    xplusd_dim1 = *ldn;
    xplusd_offset = 1 + xplusd_dim1;
    xplusd -= xplusd_offset;
    fjacd_dim1 = *ldn;
    fjacd_dim2 = *ldm;
    fjacd_offset = 1 + fjacd_dim1 * (1 + fjacd_dim2);
    fjacd -= fjacd_offset;
    fjacb_dim1 = *ldn;
    fjacb_dim2 = *ldnp;
    fjacb_offset = 1 + fjacb_dim1 * (1 + fjacb_dim2);
    fjacb -= fjacb_offset;
    ifixx_dim1 = *ldifx;
    ifixx_offset = 1 + ifixx_dim1;
    ifixx -= ifixx_offset;

    /* Function Body */
/* ...VARIABLE DEFINITIONS (ALPHABETICALLY) */
/*   BETA:    CURRENT VALUES OF PARAMETERS */
/*   F:       PREDICTED FUNCTION VALUES */
/*   FAC1:    A FACTORS OR TERMS USED IN COMPUTING THE JACOBIANS. */
/*   FAC2:    A FACTORS OR TERMS USED IN COMPUTING THE JACOBIANS. */
/*   FAC3:    A FACTORS OR TERMS USED IN COMPUTING THE JACOBIANS. */
/*   FAC4:    A FACTORS OR TERMS USED IN COMPUTING THE JACOBIANS. */
/*   FJACB:   JACOBIAN WITH RESPECT TO BETA */
/*   FJACD:   JACOBIAN WITH RESPECT TO ERRORS DELTA */
/*   IDEVAL:  INDICATOR FOR SELECTING COMPUTATION TO BE PERFORMED */
/*   IFIXB:   INDICATORS FOR "FIXING" PARAMETERS (BETA) */
/*   IFIXX:   INDICATORS FOR "FIXING" EXPLANATORY VARIABLE (X) */
/*   LDIFX:   LEADING DIMENSION OF ARRAY IFIXX */
/*   ISTOP:   STOPPING CONDITION, WHERE */
/*                     0 MEANS CURRENT BETA AND X+DELTA WERE */
/*                       ACCEPTABLE AND VALUES WERE COMPUTED SUCCESSFULLY */
/*                     1 MEANS CURRENT BETA AND X+DELTA ARE */
/*                       NOT ACCEPTABLE;  ODRPACK SHOULD SELECT */
/*                       VALUES CLOSER TO MOST RECENTLY USED VALUES */
/*                    -1 MEANS CURRENT BETA AND X+DELTA ARE */
/*                       NOT ACCEPTABLE;  ODRPACK SHOULD STOP */
/*   LDN:     LEADING DIMENSION DECLARATOR EQUAL OR EXCEEDING N */
/*   LDM:     LEADING DIMENSION DECLARATOR EQUAL OR EXCEEDING M */
/*   LDNP:    LEADING DIMENSION DECLARATOR EQUAL OR EXCEEDING NP */
/*   M:       THE NUMBER OF COLUMNS OF DATA IN THE EXPLANATORY VARIABLE. */
/*   N:       THE NUMBER OF OBSERVATIONS. */
/*   NP:      THE NUMBER OF FUNCTION PARAMETERS. */
/*   NQ:      THE NUMBER OF RESPONSES PER OBSERVATION. */
/*   ONE:     THE VALUE 1.0D0. */
/*   SETNO:   THE NUMBER OF THE DATA SET BEING ANALYZED. */
/*   XPLUSD:  CURRENT VALUE OF EXPLANATORY VARIABLE, I.E., X + DELTA */
/*   ZERO:    THE VALUE 0.0D0. */
/* ***FIRST EXECUTABLE STATEMENT  DODRXF */
    if (setid_1.setno == 1) {
/*  SETNO. 1:  BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 1 */
	if (beta[1] <= 1.01) {
	    *istop = 0;
	    if (*ideval % 10 != 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    f[i__ + f_dim1] = beta[1] / (xplusd[i__ + xplusd_dim1] - 
			    beta[2]);
/* L100: */
		}
	    }
	    if (*ideval / 10 % 10 != 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = one / (
			    xplusd[i__ + xplusd_dim1] - beta[2]);
/* Computing 2nd power */
		    d__1 = 1 / (xplusd[i__ + xplusd_dim1] - beta[2]);
		    fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = beta[1] * (
			    d__1 * d__1);
/* L110: */
		}
	    }
	    if (*ideval / 100 % 10 != 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		    d__1 = 1 / (xplusd[i__ + xplusd_dim1] - beta[2]);
		    fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = -beta[1] * (
			    d__1 * d__1);
/* L120: */
		}
	    }
	} else {
	    *istop = 1;
	}
    } else if (setid_1.setno == 2) {
/*  SETNO. 2:  BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 2 */
	*istop = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fac1 = beta[2] * xplusd[i__ + xplusd_dim1] + beta[3] * xplusd[i__ 
		    + (xplusd_dim1 << 1)] - one;
	    if (*ideval % 10 != 0) {
		f[i__ + f_dim1] = beta[1] / fac1;
	    }
	    if (*ideval / 10 % 10 != 0) {
		fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = one / fac1;
/* Computing 2nd power */
		d__1 = 1 / fac1;
		fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = -beta[1] * (d__1 
			* d__1) * xplusd[i__ + xplusd_dim1];
/* Computing 2nd power */
		d__1 = 1 / fac1;
		fjacb[i__ + (fjacb_dim2 + 3) * fjacb_dim1] = -beta[1] * (d__1 
			* d__1) * xplusd[i__ + (xplusd_dim1 << 1)];
	    }
	    if (*ideval / 100 % 10 != 0) {
/* Computing 2nd power */
		d__1 = 1 / fac1;
		fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = -beta[1] * (d__1 
			* d__1) * beta[2];
/* Computing 2nd power */
		d__1 = 1 / fac1;
		fjacd[i__ + (fjacd_dim2 + 2) * fjacd_dim1] = -beta[1] * (d__1 
			* d__1) * beta[3];
	    }
/* L200: */
	}
    } else if (setid_1.setno == 3) {
/*  SETNO. 3:  BOGGS, BYRD AND SCHNABEL, 1985, EXAMPLE 3 */
	*istop = 0;
	if (*ideval % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f[i__ + f_dim1] = zero;
		for (j = 1; j <= 4; ++j) {
		    f[i__ + f_dim1] += beta[j] / (xplusd[i__ + xplusd_dim1] + 
			    beta[j + 5]);
/* L300: */
		}
		f[i__ + f_dim1] += beta[5];
/* L310: */
	    }
	}
	if (*ideval / 10 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacb[i__ + (fjacb_dim2 + 5) * fjacb_dim1] = one;
		for (k = 1; k <= 4; ++k) {
		    fjacb[i__ + (k + fjacb_dim2) * fjacb_dim1] = one / (
			    xplusd[i__ + xplusd_dim1] + beta[k + 5]);
/* Computing 2nd power */
		    d__1 = 1 / (xplusd[i__ + xplusd_dim1] + beta[k + 5]);
		    fjacb[i__ + (k + 5 + fjacb_dim2) * fjacb_dim1] = -beta[k] 
			    * (d__1 * d__1);
/* L320: */
		}
/* L330: */
	    }
	}
	if (*ideval / 100 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = zero;
		for (k = 4; k >= 1; --k) {
/* Computing 2nd power */
		    d__1 = 1 / (xplusd[i__ + xplusd_dim1] + beta[k + 5]);
		    fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] -= beta[k] * (
			    d__1 * d__1);
/* L340: */
		}
/* L350: */
	    }
	}
    } else if (setid_1.setno == 4) {
/*  SETNO. 4:  HIMMELBLAU, 1970, EXAMPLE 6.2-4, PAGE 188 */
	*istop = 0;
	if (*ideval % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f[i__ + f_dim1] = beta[1] * xplusd[i__ + xplusd_dim1] + beta[
			2] * exp(beta[3] * xplusd[i__ + (xplusd_dim1 << 1)]);
/* L400: */
	    }
	}
	if (*ideval / 10 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = xplusd[i__ + 
			xplusd_dim1];
		fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = exp(beta[3] * 
			xplusd[i__ + (xplusd_dim1 << 1)]);
		fjacb[i__ + (fjacb_dim2 + 3) * fjacb_dim1] = beta[2] * exp(
			beta[3] * xplusd[i__ + (xplusd_dim1 << 1)]) * xplusd[
			i__ + (xplusd_dim1 << 1)];
/* L410: */
	    }
	}
	if (*ideval / 100 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = beta[1];
		fjacd[i__ + (fjacd_dim2 + 2) * fjacd_dim1] = beta[2] * exp(
			beta[3] * xplusd[i__ + (xplusd_dim1 << 1)]) * beta[3];
/* L420: */
	    }
	}
    } else if (setid_1.setno == 5) {
/*  SETNO. 5:  DRAPER AND SMITH, 1981, EXERCISE I, PAGE 521-522 */
	*istop = 0;
	if (*ideval % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f[i__ + f_dim1] = exp(-beta[1] * xplusd[i__ + xplusd_dim1] * 
			exp(-beta[2] * (one / xplusd[i__ + (xplusd_dim1 << 1)]
			 - one / 620.)));
/* L500: */
	    }
	}
	if (*ideval / 10 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fac1 = one / xplusd[i__ + (xplusd_dim1 << 1)] - one / 620.;
		fac2 = exp(-beta[2] * fac1);
		fac3 = beta[1] * xplusd[i__ + xplusd_dim1];
		fac4 = exp(-fac3 * fac2);
		fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = -fac4 * xplusd[
			i__ + xplusd_dim1] * fac2;
		fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = fac4 * fac3 * 
			fac2 * fac1;
		if (*ideval / 100 % 10 != 0) {
		    fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = -fac4 * beta[
			    1] * fac2;
/* Computing 2nd power */
		    d__1 = xplusd[i__ + (xplusd_dim1 << 1)];
		    fjacd[i__ + (fjacd_dim2 + 2) * fjacd_dim1] = -fac4 * fac3 
			    * fac2 * beta[2] / (d__1 * d__1);
		}
/* L510: */
	    }
	}
    } else if (setid_1.setno == 6) {
/*  SETNO. 6:  POWELL AND MACDONALD, 1972, TABLES 7 AND 8, PAGE 153-154 */
/*             N.B.  THIS DERIVATIVE IS INTENTIONALLY CODED INCORRECTLY */
	*istop = 0;
	if (*ideval % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = one + beta[3] * xplusd[i__ + xplusd_dim1] / beta[2];
		d__2 = -one / beta[3];
		f[i__ + f_dim1] = beta[1] * pow_dd(&d__1, &d__2);
/* L600: */
	    }
	}
	if (*ideval / 10 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 3) * fjacb_dim1] = zero;
		if (*ideval / 100 % 10 != 0) {
		    fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = xplusd[i__ + 
			    xplusd_dim1];
		}
/* L610: */
	    }
	}
    } else if (setid_1.setno == 7) {
/*  SETNO. 7:  FULLER, 1987, TABLE 3.2.10, PAGES 244-245 */
/*             N.B.  THIS DERIVATIVE IS INTENTIONALLY CODED INCORRECTLY */
	*istop = 0;
	if (*ideval % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = xplusd[i__ + xplusd_dim1] - beta[1];
/* Computing 2nd power */
		d__2 = xplusd[i__ + (xplusd_dim1 << 1)] - beta[2];
		f[i__ + f_dim1] = beta[3] * (d__1 * d__1) + beta[4] * 2 * (
			xplusd[i__ + xplusd_dim1] - beta[1]) * (xplusd[i__ + (
			xplusd_dim1 << 1)] - beta[2]) + beta[5] * (d__2 * 
			d__2) - 1.;
/* L700: */
	    }
	}
	if (*ideval / 10 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 3) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 4) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 5) * fjacb_dim1] = zero;
		if (*ideval / 100 % 10 != 0) {
		    fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = zero;
		    fjacd[i__ + (fjacd_dim2 + 2) * fjacd_dim1] = zero;
		}
/* L710: */
	    }
	}
    } else if (setid_1.setno == 8) {
/*  SETNO. 8:  BATES AND WATTS, 1988, TABLE A1.13, PAGES 280-281 */
/*             N.B.  THIS DERIVATIVE IS INTENTIONALLY CODED INCORRECTLY */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (xplusd[i__ + xplusd_dim1] < 0.) {
		*istop = 1;
		return 0;
	    }
/* L800: */
	}
	*istop = 0;
	if (*ideval % 10 != 0) {
	    pi = 3.141592653589793238462643383279;
	    theta = pi * beta[4] * .5;
	    ctheta = cos(theta);
	    stheta = sin(theta);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		freq = xplusd[i__ + xplusd_dim1];
		d__1 = pi * 2. * freq * exp(-beta[3]);
		omega = pow_dd(&d__1, &beta[4]);
		phi = atan2(omega * stheta, omega * ctheta + 1);
/* Computing 2nd power */
		d__2 = omega * ctheta + 1;
/* Computing 2nd power */
		d__3 = omega * stheta;
		d__1 = sqrt(d__2 * d__2 + d__3 * d__3);
		d__4 = -beta[5];
		r__ = (beta[1] - beta[2]) * pow_dd(&d__1, &d__4);
		f[i__ + f_dim1] = beta[2] + r__ * cos(beta[5] * phi);
		f[i__ + (f_dim1 << 1)] = r__ * sin(beta[5] * phi);
/* L810: */
	    }
	}
	if (*ideval / 10 % 10 != 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fjacb[i__ + (fjacb_dim2 + 1) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 2) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 3) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 4) * fjacb_dim1] = zero;
		fjacb[i__ + (fjacb_dim2 + 5) * fjacb_dim1] = zero;
		fjacb[i__ + ((fjacb_dim2 << 1) + 1) * fjacb_dim1] = zero;
		fjacb[i__ + ((fjacb_dim2 << 1) + 2) * fjacb_dim1] = zero;
		fjacb[i__ + ((fjacb_dim2 << 1) + 3) * fjacb_dim1] = zero;
		fjacb[i__ + ((fjacb_dim2 << 1) + 4) * fjacb_dim1] = zero;
		fjacb[i__ + ((fjacb_dim2 << 1) + 5) * fjacb_dim1] = zero;
		if (*ideval / 100 % 10 != 0) {
		    fjacd[i__ + (fjacd_dim2 + 1) * fjacd_dim1] = zero;
		    fjacd[i__ + ((fjacd_dim2 << 1) + 1) * fjacd_dim1] = zero;
		}
/* L820: */
	    }
	}
    }
    return 0;
} /* dodrxf_ */

/* DODRXW */
/* Subroutine */ int dodrxw_(integer *maxn, integer *maxm, integer *maxnp, 
	integer *maxnq, integer *ldwe, integer *ld2we, logical *isodr, 
	integer *liwmin, integer *lwmin)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

/* ***BEGIN PROLOGUE  DODRXW */
/* ***REFER TO  DODR,DODRC */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   890205   (YYMMDD) */
/* ***REVISION DATE  920619   (YYMMDD) */
/* ***PURPOSE  COMPUTE MINIMUM LENGTHS FOR WORK VECTORS */
/* ***ROUTINES CALLED  NONE */
/* ***END PROLOGUE  DODRXW */
/* ...SCALAR ARGUMENTS */
/* ...VARIABLE DEFINITIONS (ALPHABETICALLY) */
/*   ISODR:   THE VARIABLE DESIGNATING WHETHER THE SOLUTION IS BY ODR */
/*            (ISODR=TRUE) OR BY OLS (ISODR=FALSE). */
/*   LDWE:    THE LEADING DIMENSION OF ARRAY WE. */
/*   LD2WE:   THE SECOND DIMENSION OF ARRAY WE. */
/*   LIWMIN:  THE MINIMUM LENGTH OF VECTOR IWORK FOR A GIVEN PROBLEM. */
/*   LWMIN:   THE MINIMUM LENGTH OF VECTOR WORK FOR A GIVEN PROBLEM. */
/*   MAXM:    THE NUMBER OF COLUMNS IN THE EXPLANATORY VARIABLE. */
/*   MAXN:    THE NUMBER OF OBSERVATIONS. */
/*   MAXNP:   THE NUMBER OF FUNCTION PARAMETERS. */
/*   MAXNQ:   THE NUMBER OF RESPONSES PER OBSERVATION. */
/* ***FIRST EXECUTABLE STATEMENT  DODRXW */
    *liwmin = *maxnp + 20 + *maxnq * (*maxnp + *maxm);
    if (*isodr) {
/* Computing 2nd power */
	i__1 = *maxnp;
/* Computing 2nd power */
	i__2 = *maxm;
/* Computing 2nd power */
	i__3 = *maxnq;
	*lwmin = *maxnp * 11 + 18 + i__1 * i__1 + *maxm + i__2 * i__2 + (*
		maxn << 2) * *maxnq + *maxn * 6 * *maxm + (*maxn << 1) * *
		maxnq * *maxnp + (*maxn << 1) * *maxnq * *maxm + i__3 * i__3 
		+ *maxnq * 5 + *maxnq * (*maxnp + *maxm) + *ldwe * *ld2we * *
		maxnq;
    } else {
/* Computing 2nd power */
	i__1 = *maxnp;
/* Computing 2nd power */
	i__2 = *maxm;
	*lwmin = *maxnp * 11 + 18 + i__1 * i__1 + *maxm + i__2 * i__2 + (*
		maxn << 2) * *maxnq + (*maxn << 1) * *maxm + (*maxn << 1) * *
		maxnq * *maxnp + *maxnq * 5 + *maxnq * (*maxnp + *maxm) + *
		ldwe * *ld2we * *maxnq;
    }
    return 0;
} /* dodrxw_ */

/* Main program alias */ int dtest_ () { MAIN__ (); return 0; }
