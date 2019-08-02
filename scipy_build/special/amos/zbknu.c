/* zbknu.f -- translated by f2c (version 20190311).
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

static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b10 = .5;
static doublereal c_b11 = 0.;
static integer c__14 = 14;
static integer c__5 = 5;

/* Subroutine */ int zbknu_(doublereal *zr, doublereal *zi, doublereal *fnu, 
	integer *kode, integer *n, doublereal *yr, doublereal *yi, integer *
	nz, doublereal *tol, doublereal *elim, doublereal *alim)
{
    /* Initialized data */

    static integer kmax = 30;
    static doublereal spi = 1.90985931710274403;
    static doublereal hpi = 1.57079632679489662;
    static doublereal fpi = 1.89769999331517738;
    static doublereal tth = .666666666666666666;
    static doublereal cc[8] = { .577215664901532861,-.0420026350340952355,
	    -.0421977345555443367,.00721894324666309954,
	    -2.15241674114950973e-4,-2.01348547807882387e-5,
	    1.13302723198169588e-6,6.11609510448141582e-9 };
    static doublereal czeror = 0.;
    static doublereal czeroi = 0.;
    static doublereal coner = 1.;
    static doublereal conei = 0.;
    static doublereal ctwor = 2.;
    static doublereal r1 = 2.;
    static doublereal dpi = 3.14159265358979324;
    static doublereal rthpi = 1.25331413731550025;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), exp(doublereal), cos(doublereal), atan(doublereal)
	    , sqrt(doublereal), log(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, a1, a2, g1, g2, t1, t2, aa, bb, fc, ak, bk;
    static integer ic;
    static doublereal fi, fk, as;
    static integer kk;
    static doublereal fr, pi, qi, tm, pr, qr;
    static integer nw;
    static doublereal p1i, p2i, s1i, s2i, p2m, p1r, p2r, s1r, s2r, cbi, cbr, 
	    cki, caz, csi, ckr, fhs, fks, rak, czi, dnu, csr, elm, zdi, bry[3]
	    , pti, czr, sti, zdr, cyr[2], rzi, ptr, cyi[2];
    static integer inu;
    static doublereal str, rzr, dnu2, cchi, cchr, alas, cshi;
    static integer inub, idum;
    static doublereal cshr, fmui, rcaz, csrr[3], cssr[3], fmur;
    extern /* Subroutine */ int zdiv_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static doublereal smui, smur;
    extern /* Subroutine */ int zmlt_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static integer iflag, kflag;
    static doublereal coefi;
    static integer koded;
    static doublereal ascle, coefr, helim;
    extern doublereal azabs_(doublereal *, doublereal *);
    static doublereal celmr, csclr, crscr;
    extern /* Subroutine */ int azlog_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *), zshch_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal etest;
    extern /* Subroutine */ int zuchk_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *), azexp_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), zkscl_(doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    extern doublereal dgamln_(doublereal *, integer *);
    extern /* Subroutine */ int azsqrt_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  ZBKNU */
/* ***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH */

/*     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE. */

/* ***ROUTINES CALLED  DGAMLN,I1MACH,D1MACH,ZKSCL,ZSHCH,ZUCHK,AZABS,ZDIV, */
/*                    AZEXP,AZLOG,ZMLT,AZSQRT */
/* ***END PROLOGUE  ZBKNU */

/*     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH */
/*     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK */

    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    caz = azabs_(zr, zi);
    csclr = 1. / *tol;
    crscr = *tol;
    cssr[0] = csclr;
    cssr[1] = 1.;
    cssr[2] = crscr;
    csrr[0] = crscr;
    csrr[1] = 1.;
    csrr[2] = csclr;
    bry[0] = d1mach_(&c__1) * 1e3 / *tol;
    bry[1] = 1. / bry[0];
    bry[2] = d1mach_(&c__2);
    *nz = 0;
    iflag = 0;
    koded = *kode;
    rcaz = 1. / caz;
    str = *zr * rcaz;
    sti = -(*zi) * rcaz;
    rzr = (str + str) * rcaz;
    rzi = (sti + sti) * rcaz;
    inu = (integer) ((real) (*fnu + .5));
    dnu = *fnu - (doublereal) ((real) inu);
    if (abs(dnu) == .5) {
	goto L110;
    }
    dnu2 = 0.;
    if (abs(dnu) > *tol) {
	dnu2 = dnu * dnu;
    }
    if (caz > r1) {
	goto L110;
    }
/* ----------------------------------------------------------------------- */
/*     SERIES FOR CABS(Z).LE.R1 */
/* ----------------------------------------------------------------------- */
    fc = 1.;
    azlog_(&rzr, &rzi, &smur, &smui, &idum);
    fmur = smur * dnu;
    fmui = smui * dnu;
    zshch_(&fmur, &fmui, &cshr, &cshi, &cchr, &cchi);
    if (dnu == 0.) {
	goto L10;
    }
    fc = dnu * dpi;
    fc /= sin(fc);
    smur = cshr / dnu;
    smui = cshi / dnu;
L10:
    a2 = dnu + 1.;
/* ----------------------------------------------------------------------- */
/*     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU) */
/* ----------------------------------------------------------------------- */
    t2 = exp(-dgamln_(&a2, &idum));
    t1 = 1. / (t2 * fc);
    if (abs(dnu) > .1) {
	goto L40;
    }
/* ----------------------------------------------------------------------- */
/*     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU) */
/* ----------------------------------------------------------------------- */
    ak = 1.;
    s = cc[0];
    for (k = 2; k <= 8; ++k) {
	ak *= dnu2;
	tm = cc[k - 1] * ak;
	s += tm;
	if (abs(tm) < *tol) {
	    goto L30;
	}
/* L20: */
    }
L30:
    g1 = -s;
    goto L50;
L40:
    g1 = (t1 - t2) / (dnu + dnu);
L50:
    g2 = (t1 + t2) * .5;
    fr = fc * (cchr * g1 + smur * g2);
    fi = fc * (cchi * g1 + smui * g2);
    azexp_(&fmur, &fmui, &str, &sti);
    pr = str * .5 / t2;
    pi = sti * .5 / t2;
    zdiv_(&c_b10, &c_b11, &str, &sti, &ptr, &pti);
    qr = ptr / t1;
    qi = pti / t1;
    s1r = fr;
    s1i = fi;
    s2r = pr;
    s2i = pi;
    ak = 1.;
    a1 = 1.;
    ckr = coner;
    cki = conei;
    bk = 1. - dnu2;
    if (inu > 0 || *n > 1) {
	goto L80;
    }
/* ----------------------------------------------------------------------- */
/*     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1 */
/* ----------------------------------------------------------------------- */
    if (caz < *tol) {
	goto L70;
    }
    zmlt_(zr, zi, zr, zi, &czr, &czi);
    czr *= .25;
    czi *= .25;
    t1 = caz * .25 * caz;
L60:
    fr = (fr * ak + pr + qr) / bk;
    fi = (fi * ak + pi + qi) / bk;
    str = 1. / (ak - dnu);
    pr *= str;
    pi *= str;
    str = 1. / (ak + dnu);
    qr *= str;
    qi *= str;
    str = ckr * czr - cki * czi;
    rak = 1. / ak;
    cki = (ckr * czi + cki * czr) * rak;
    ckr = str * rak;
    s1r = ckr * fr - cki * fi + s1r;
    s1i = ckr * fi + cki * fr + s1i;
    a1 = a1 * t1 * rak;
    bk = bk + ak + ak + 1.;
    ak += 1.;
    if (a1 > *tol) {
	goto L60;
    }
L70:
    yr[1] = s1r;
    yi[1] = s1i;
    if (koded == 1) {
	return 0;
    }
    azexp_(zr, zi, &str, &sti);
    zmlt_(&s1r, &s1i, &str, &sti, &yr[1], &yi[1]);
    return 0;
/* ----------------------------------------------------------------------- */
/*     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE */
/* ----------------------------------------------------------------------- */
L80:
    if (caz < *tol) {
	goto L100;
    }
    zmlt_(zr, zi, zr, zi, &czr, &czi);
    czr *= .25;
    czi *= .25;
    t1 = caz * .25 * caz;
L90:
    fr = (fr * ak + pr + qr) / bk;
    fi = (fi * ak + pi + qi) / bk;
    str = 1. / (ak - dnu);
    pr *= str;
    pi *= str;
    str = 1. / (ak + dnu);
    qr *= str;
    qi *= str;
    str = ckr * czr - cki * czi;
    rak = 1. / ak;
    cki = (ckr * czi + cki * czr) * rak;
    ckr = str * rak;
    s1r = ckr * fr - cki * fi + s1r;
    s1i = ckr * fi + cki * fr + s1i;
    str = pr - fr * ak;
    sti = pi - fi * ak;
    s2r = ckr * str - cki * sti + s2r;
    s2i = ckr * sti + cki * str + s2i;
    a1 = a1 * t1 * rak;
    bk = bk + ak + ak + 1.;
    ak += 1.;
    if (a1 > *tol) {
	goto L90;
    }
L100:
    kflag = 2;
    a1 = *fnu + 1.;
    ak = a1 * abs(smur);
    if (ak > *alim) {
	kflag = 3;
    }
    str = cssr[kflag - 1];
    p2r = s2r * str;
    p2i = s2i * str;
    zmlt_(&p2r, &p2i, &rzr, &rzi, &s2r, &s2i);
    s1r *= str;
    s1i *= str;
    if (koded == 1) {
	goto L210;
    }
    azexp_(zr, zi, &fr, &fi);
    zmlt_(&s1r, &s1i, &fr, &fi, &s1r, &s1i);
    zmlt_(&s2r, &s2i, &fr, &fi, &s2r, &s2i);
    goto L210;
/* ----------------------------------------------------------------------- */
/*     IFLAG=0 MEANS NO UNDERFLOW OCCURRED */
/*     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH */
/*     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD */
/*     RECURSION */
/* ----------------------------------------------------------------------- */
L110:
    azsqrt_(zr, zi, &str, &sti);
    zdiv_(&rthpi, &czeroi, &str, &sti, &coefr, &coefi);
    kflag = 2;
    if (koded == 2) {
	goto L120;
    }
    if (*zr > *alim) {
	goto L290;
    }
/*     BLANK LINE */
    str = exp(-(*zr)) * cssr[kflag - 1];
    sti = -str * sin(*zi);
    str *= cos(*zi);
    zmlt_(&coefr, &coefi, &str, &sti, &coefr, &coefi);
L120:
    if (abs(dnu) == .5) {
	goto L300;
    }
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM FOR CABS(Z).GT.R1 */
/* ----------------------------------------------------------------------- */
    ak = cos(dpi * dnu);
    ak = abs(ak);
    if (ak == czeror) {
	goto L300;
    }
    fhs = (d__1 = .25 - dnu2, abs(d__1));
    if (fhs == czeror) {
	goto L300;
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO */
/*     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON */
/*     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))= */
/*     TOL WHERE B IS THE BASE OF THE ARITHMETIC. */
/* ----------------------------------------------------------------------- */
    t1 = (doublereal) ((real) (i1mach_(&c__14) - 1));
    t1 = t1 * d1mach_(&c__5) * 3.321928094;
    t1 = max(t1,12.);
    t1 = min(t1,60.);
    t2 = tth * t1 - 6.;
    if (*zr != 0.) {
	goto L130;
    }
    t1 = hpi;
    goto L140;
L130:
    t1 = atan(*zi / *zr);
    t1 = abs(t1);
L140:
    if (t2 > caz) {
	goto L170;
    }
/* ----------------------------------------------------------------------- */
/*     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2 */
/* ----------------------------------------------------------------------- */
    etest = ak / (dpi * caz * *tol);
    fk = coner;
    if (etest < coner) {
	goto L180;
    }
    fks = ctwor;
    ckr = caz + caz + ctwor;
    p1r = czeror;
    p2r = coner;
    i__1 = kmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ak = fhs / fks;
	cbr = ckr / (fk + coner);
	ptr = p2r;
	p2r = cbr * p2r - p1r * ak;
	p1r = ptr;
	ckr += ctwor;
	fks = fks + fk + fk + ctwor;
	fhs = fhs + fk + fk;
	fk += coner;
	str = abs(p2r) * fk;
	if (etest < str) {
	    goto L160;
	}
/* L150: */
    }
    goto L310;
L160:
    fk += spi * t1 * sqrt(t2 / caz);
    fhs = (d__1 = .25 - dnu2, abs(d__1));
    goto L180;
L170:
/* ----------------------------------------------------------------------- */
/*     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2 */
/* ----------------------------------------------------------------------- */
    a2 = sqrt(caz);
    ak = fpi * ak / (*tol * sqrt(a2));
    aa = t1 * 3. / (caz + 1.);
    bb = t1 * 14.7 / (caz + 28.);
    ak = (log(ak) + caz * cos(aa) / (caz * .008 + 1.)) / cos(bb);
    fk = ak * .12125 * ak / caz + 1.5;
L180:
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM */
/* ----------------------------------------------------------------------- */
    k = (integer) ((real) fk);
    fk = (doublereal) ((real) k);
    fks = fk * fk;
    p1r = czeror;
    p1i = czeroi;
    p2r = *tol;
    p2i = czeroi;
    csr = p2r;
    csi = p2i;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a1 = fks - fk;
	ak = (fks + fk) / (a1 + fhs);
	rak = 2. / (fk + coner);
	cbr = (fk + *zr) * rak;
	cbi = *zi * rak;
	ptr = p2r;
	pti = p2i;
	p2r = (ptr * cbr - pti * cbi - p1r) * ak;
	p2i = (pti * cbr + ptr * cbi - p1i) * ak;
	p1r = ptr;
	p1i = pti;
	csr += p2r;
	csi += p2i;
	fks = a1 - fk + coner;
	fk -= coner;
/* L190: */
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER */
/*     SCALING */
/* ----------------------------------------------------------------------- */
    tm = azabs_(&csr, &csi);
    ptr = 1. / tm;
    s1r = p2r * ptr;
    s1i = p2i * ptr;
    csr *= ptr;
    csi = -csi * ptr;
    zmlt_(&coefr, &coefi, &s1r, &s1i, &str, &sti);
    zmlt_(&str, &sti, &csr, &csi, &s1r, &s1i);
    if (inu > 0 || *n > 1) {
	goto L200;
    }
    zdr = *zr;
    zdi = *zi;
    if (iflag == 1) {
	goto L270;
    }
    goto L240;
L200:
/* ----------------------------------------------------------------------- */
/*     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING */
/* ----------------------------------------------------------------------- */
    tm = azabs_(&p2r, &p2i);
    ptr = 1. / tm;
    p1r *= ptr;
    p1i *= ptr;
    p2r *= ptr;
    p2i = -p2i * ptr;
    zmlt_(&p1r, &p1i, &p2r, &p2i, &ptr, &pti);
    str = dnu + .5 - ptr;
    sti = -pti;
    zdiv_(&str, &sti, zr, zi, &str, &sti);
    str += 1.;
    zmlt_(&str, &sti, &s1r, &s1i, &s2r, &s2i);
/* ----------------------------------------------------------------------- */
/*     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH */
/*     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3 */
/* ----------------------------------------------------------------------- */
L210:
    str = dnu + 1.;
    ckr = str * rzr;
    cki = str * rzi;
    if (*n == 1) {
	--inu;
    }
    if (inu > 0) {
	goto L220;
    }
    if (*n > 1) {
	goto L215;
    }
    s1r = s2r;
    s1i = s2i;
L215:
    zdr = *zr;
    zdi = *zi;
    if (iflag == 1) {
	goto L270;
    }
    goto L240;
L220:
    inub = 1;
    if (iflag == 1) {
	goto L261;
    }
L225:
    p1r = csrr[kflag - 1];
    ascle = bry[kflag - 1];
    i__1 = inu;
    for (i__ = inub; i__ <= i__1; ++i__) {
	str = s2r;
	sti = s2i;
	s2r = ckr * str - cki * sti + s1r;
	s2i = ckr * sti + cki * str + s1i;
	s1r = str;
	s1i = sti;
	ckr += rzr;
	cki += rzi;
	if (kflag >= 3) {
	    goto L230;
	}
	p2r = s2r * p1r;
	p2i = s2i * p1r;
	str = abs(p2r);
	sti = abs(p2i);
	p2m = max(str,sti);
	if (p2m <= ascle) {
	    goto L230;
	}
	++kflag;
	ascle = bry[kflag - 1];
	s1r *= p1r;
	s1i *= p1r;
	s2r = p2r;
	s2i = p2i;
	str = cssr[kflag - 1];
	s1r *= str;
	s1i *= str;
	s2r *= str;
	s2i *= str;
	p1r = csrr[kflag - 1];
L230:
	;
    }
    if (*n != 1) {
	goto L240;
    }
    s1r = s2r;
    s1i = s2i;
L240:
    str = csrr[kflag - 1];
    yr[1] = s1r * str;
    yi[1] = s1i * str;
    if (*n == 1) {
	return 0;
    }
    yr[2] = s2r * str;
    yi[2] = s2i * str;
    if (*n == 2) {
	return 0;
    }
    kk = 2;
L250:
    ++kk;
    if (kk > *n) {
	return 0;
    }
    p1r = csrr[kflag - 1];
    ascle = bry[kflag - 1];
    i__1 = *n;
    for (i__ = kk; i__ <= i__1; ++i__) {
	p2r = s2r;
	p2i = s2i;
	s2r = ckr * p2r - cki * p2i + s1r;
	s2i = cki * p2r + ckr * p2i + s1i;
	s1r = p2r;
	s1i = p2i;
	ckr += rzr;
	cki += rzi;
	p2r = s2r * p1r;
	p2i = s2i * p1r;
	yr[i__] = p2r;
	yi[i__] = p2i;
	if (kflag >= 3) {
	    goto L260;
	}
	str = abs(p2r);
	sti = abs(p2i);
	p2m = max(str,sti);
	if (p2m <= ascle) {
	    goto L260;
	}
	++kflag;
	ascle = bry[kflag - 1];
	s1r *= p1r;
	s1i *= p1r;
	s2r = p2r;
	s2i = p2i;
	str = cssr[kflag - 1];
	s1r *= str;
	s1i *= str;
	s2r *= str;
	s2i *= str;
	p1r = csrr[kflag - 1];
L260:
	;
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW */
/* ----------------------------------------------------------------------- */
L261:
    helim = *elim * .5;
    elm = exp(-(*elim));
    celmr = elm;
    ascle = bry[0];
    zdr = *zr;
    zdi = *zi;
    ic = -1;
    j = 2;
    i__1 = inu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	str = s2r;
	sti = s2i;
	s2r = str * ckr - sti * cki + s1r;
	s2i = sti * ckr + str * cki + s1i;
	s1r = str;
	s1i = sti;
	ckr += rzr;
	cki += rzi;
	as = azabs_(&s2r, &s2i);
	alas = log(as);
	p2r = -zdr + alas;
	if (p2r < -(*elim)) {
	    goto L263;
	}
	azlog_(&s2r, &s2i, &str, &sti, &idum);
	p2r = -zdr + str;
	p2i = -zdi + sti;
	p2m = exp(p2r) / *tol;
	p1r = p2m * cos(p2i);
	p1i = p2m * sin(p2i);
	zuchk_(&p1r, &p1i, &nw, &ascle, tol);
	if (nw != 0) {
	    goto L263;
	}
	j = 3 - j;
	cyr[j - 1] = p1r;
	cyi[j - 1] = p1i;
	if (ic == i__ - 1) {
	    goto L264;
	}
	ic = i__;
	goto L262;
L263:
	if (alas < helim) {
	    goto L262;
	}
	zdr -= *elim;
	s1r *= celmr;
	s1i *= celmr;
	s2r *= celmr;
	s2i *= celmr;
L262:
	;
    }
    if (*n != 1) {
	goto L270;
    }
    s1r = s2r;
    s1i = s2i;
    goto L270;
L264:
    kflag = 1;
    inub = i__ + 1;
    s2r = cyr[j - 1];
    s2i = cyi[j - 1];
    j = 3 - j;
    s1r = cyr[j - 1];
    s1i = cyi[j - 1];
    if (inub <= inu) {
	goto L225;
    }
    if (*n != 1) {
	goto L240;
    }
    s1r = s2r;
    s1i = s2i;
    goto L240;
L270:
    yr[1] = s1r;
    yi[1] = s1i;
    if (*n == 1) {
	goto L280;
    }
    yr[2] = s2r;
    yi[2] = s2i;
L280:
    ascle = bry[0];
    zkscl_(&zdr, &zdi, fnu, n, &yr[1], &yi[1], nz, &rzr, &rzi, &ascle, tol, 
	    elim);
    inu = *n - *nz;
    if (inu <= 0) {
	return 0;
    }
    kk = *nz + 1;
    s1r = yr[kk];
    s1i = yi[kk];
    yr[kk] = s1r * csrr[0];
    yi[kk] = s1i * csrr[0];
    if (inu == 1) {
	return 0;
    }
    kk = *nz + 2;
    s2r = yr[kk];
    s2i = yi[kk];
    yr[kk] = s2r * csrr[0];
    yi[kk] = s2i * csrr[0];
    if (inu == 2) {
	return 0;
    }
    t2 = *fnu + (doublereal) ((real) (kk - 1));
    ckr = t2 * rzr;
    cki = t2 * rzi;
    kflag = 1;
    goto L250;
L290:
/* ----------------------------------------------------------------------- */
/*     SCALE BY DEXP(Z), IFLAG = 1 CASES */
/* ----------------------------------------------------------------------- */
    koded = 2;
    iflag = 1;
    kflag = 2;
    goto L120;
/* ----------------------------------------------------------------------- */
/*     FNU=HALF ODD INTEGER CASE, DNU=-0.5 */
/* ----------------------------------------------------------------------- */
L300:
    s1r = coefr;
    s1i = coefi;
    s2r = coefr;
    s2i = coefi;
    goto L210;


L310:
    *nz = -2;
    return 0;
} /* zbknu_ */

