/* dgamln.f -- translated by f2c (version 20190311).
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

static integer c__4 = 4;
static integer c__14 = 14;
static integer c__5 = 5;

doublereal dgamln_(doublereal *z__, integer *ierr)
{
    /* Initialized data */

    static doublereal gln[100] = { 0.,0.,.693147180559945309,
	    1.791759469228055,3.17805383034794562,4.78749174278204599,
	    6.579251212010101,8.5251613610654143,10.6046029027452502,
	    12.8018274800814696,15.1044125730755153,17.5023078458738858,
	    19.9872144956618861,22.5521638531234229,25.1912211827386815,
	    27.8992713838408916,30.6718601060806728,33.5050734501368889,
	    36.3954452080330536,39.339884187199494,42.335616460753485,
	    45.380138898476908,48.4711813518352239,51.6066755677643736,
	    54.7847293981123192,58.0036052229805199,61.261701761002002,
	    64.5575386270063311,67.889743137181535,71.257038967168009,
	    74.6582363488301644,78.0922235533153106,81.5579594561150372,
	    85.0544670175815174,88.5808275421976788,92.1361756036870925,
	    95.7196945421432025,99.3306124547874269,102.968198614513813,
	    106.631760260643459,110.320639714757395,114.034211781461703,
	    117.771881399745072,121.533081515438634,125.317271149356895,
	    129.123933639127215,132.95257503561631,136.802722637326368,
	    140.673923648234259,144.565743946344886,148.477766951773032,
	    152.409592584497358,156.360836303078785,160.331128216630907,
	    164.320112263195181,168.327445448427652,172.352797139162802,
	    176.395848406997352,180.456291417543771,184.533828861449491,
	    188.628173423671591,192.739047287844902,196.866181672889994,
	    201.009316399281527,205.168199482641199,209.342586752536836,
	    213.532241494563261,217.736934113954227,221.956441819130334,
	    226.190548323727593,230.439043565776952,234.701723442818268,
	    238.978389561834323,243.268849002982714,247.572914096186884,
	    251.890402209723194,256.221135550009525,260.564940971863209,
	    264.921649798552801,269.291097651019823,273.673124285693704,
	    278.067573440366143,282.474292687630396,286.893133295426994,
	    291.323950094270308,295.766601350760624,300.220948647014132,
	    304.686856765668715,309.164193580146922,313.652829949879062,
	    318.152639620209327,322.663499126726177,327.185287703775217,
	    331.717887196928473,336.261181979198477,340.815058870799018,
	    345.379407062266854,349.954118040770237,354.539085519440809,
	    359.134205369575399 };
    static doublereal cf[22] = { .0833333333333333333,-.00277777777777777778,
	    7.93650793650793651e-4,-5.95238095238095238e-4,
	    8.41750841750841751e-4,-.00191752691752691753,
	    .00641025641025641026,-.0295506535947712418,.179644372368830573,
	    -1.39243221690590112,13.402864044168392,-156.848284626002017,
	    2193.10333333333333,-36108.7712537249894,691472.268851313067,
	    -15238221.5394074162,382900751.391414141,-10882266035.7843911,
	    347320283765.002252,-12369602142269.2745,488788064793079.335,
	    -21320333960919373.9 };
    static doublereal con = 1.83787706640934548;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal s, t1, fz, zm;
    static integer mz, nz;
    static doublereal zp;
    static integer i1m;
    static doublereal fln, tlg, rln, trm, tst, zsq, zinc, zmin, zdmy, wdtol;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  DGAMLN */
/* ***DATE WRITTEN   830501   (YYMMDD) */
/* ***REVISION DATE  830501   (YYMMDD) */
/* ***CATEGORY NO.  B5F */
/* ***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION */
/* ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES */
/* ***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION */
/* ***DESCRIPTION */

/*               **** A DOUBLE PRECISION ROUTINE **** */
/*         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR */
/*         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES */
/*         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION */
/*         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS */
/*         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE */
/*         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18) */
/*         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY. */

/*         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100 */
/*         VALUES IS USED FOR SPEED OF EXECUTION. */

/*     DESCRIPTION OF ARGUMENTS */

/*         INPUT      Z IS D0UBLE PRECISION */
/*           Z      - ARGUMENT, Z.GT.0.0D0 */

/*         OUTPUT      DGAMLN IS DOUBLE PRECISION */
/*           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0 */
/*           IERR    - ERROR FLAG */
/*                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED */
/*                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION */


/* ***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT */
/*                 BY D. E. AMOS, SAND83-0083, MAY, 1983. */
/* ***ROUTINES CALLED  I1MACH,D1MACH */
/* ***END PROLOGUE  DGAMLN */
/*           LNGAMMA(N), N=1,100 */
/*             COEFFICIENTS OF ASYMPTOTIC EXPANSION */

/*             LN(2*PI) */

/* ***FIRST EXECUTABLE STATEMENT  DGAMLN */
    *ierr = 0;
    if (*z__ <= 0.) {
	goto L70;
    }
    if (*z__ > 101.) {
	goto L10;
    }
    nz = (integer) ((real) (*z__));
    fz = *z__ - (real) nz;
    if (fz > 0.) {
	goto L10;
    }
    if (nz > 100) {
	goto L10;
    }
    ret_val = gln[nz - 1];
    return ret_val;
L10:
    wdtol = d1mach_(&c__4);
    wdtol = max(wdtol,5e-19);
    i1m = i1mach_(&c__14);
    rln = d1mach_(&c__5) * (real) i1m;
    fln = min(rln,20.);
    fln = max(fln,3.);
    fln += -3.;
    zm = fln * .3875 + 1.8;
    mz = (integer) ((real) zm) + 1;
    zmin = (real) mz;
    zdmy = *z__;
    zinc = 0.;
    if (*z__ >= zmin) {
	goto L20;
    }
    zinc = zmin - (real) nz;
    zdmy = *z__ + zinc;
L20:
    zp = 1. / zdmy;
    t1 = cf[0] * zp;
    s = t1;
    if (zp < wdtol) {
	goto L40;
    }
    zsq = zp * zp;
    tst = t1 * wdtol;
    for (k = 2; k <= 22; ++k) {
	zp *= zsq;
	trm = cf[k - 1] * zp;
	if (abs(trm) < tst) {
	    goto L40;
	}
	s += trm;
/* L30: */
    }
L40:
    if (zinc != 0.) {
	goto L50;
    }
    tlg = log(*z__);
    ret_val = *z__ * (tlg - 1.) + (con - tlg) * .5 + s;
    return ret_val;
L50:
    zp = 1.;
    nz = (integer) ((real) zinc);
    i__1 = nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zp *= *z__ + (real) (i__ - 1);
/* L60: */
    }
    tlg = log(zdmy);
    ret_val = zdmy * (tlg - 1.) - log(zp) + (con - tlg) * .5 + s;
    return ret_val;


L70:
    *ierr = 1;
    return ret_val;
} /* dgamln_ */

