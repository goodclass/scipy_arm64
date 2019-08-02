/* idzp_rsvd_subr_0.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int idzp_rsvd_(integer *lw, doublereal *eps, integer *m, 
	integer *n, U_fp matveca, doublecomplex *p1t, doublecomplex *p2t, 
	doublecomplex *p3t, doublecomplex *p4t, U_fp matvec, doublecomplex *
	p1, doublecomplex *p2, doublecomplex *p3, doublecomplex *p4, integer *
	krank, integer *iu, integer *iv, integer *is, doublecomplex *w, 
	integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int idz_reco_(integer *, doublecomplex *, 
	    doublecomplex *), idzp_rid_(integer *, doublereal *, integer *, 
	    integer *, U_fp, doublecomplex *, doublecomplex *, doublecomplex *
	    , doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *);
    static integer k, lp, ls, lu, lv;
    extern /* Subroutine */ int idzp_rsvd0_(integer *, integer *, U_fp, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , U_fp, doublecomplex *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    static integer lw2, isi, iui, ivi, icol, lcol, iproj, ilist, lproj, llist,
	     iwork, lwork;


/*       constructs a rank-krank SVD  U Sigma V^*  approximating a */
/*       to precision eps, where matveca is a routine which applies a^* */
/*       to an arbitrary vector, and matvec is a routine */
/*       which applies a to an arbitrary vector; U is an m x krank */
/*       matrix whose columns are orthonormal, V is an n x krank */
/*       matrix whose columns are orthonormal, and Sigma is a diagonal */
/*       krank x krank matrix whose entries are all nonnegative. */
/*       The entries of U are stored in w, starting at w(iu); */
/*       the entries of V are stored in w, starting at w(iv). */
/*       The diagonal entries of Sigma are stored in w, */
/*       starting at w(is). This routine uses a randomized algorithm. */

/*       input: */
/*       lw -- maximum usable length (in complex*16 elements) */
/*             of the array w */
/*       eps -- precision of the desired approximation */
/*       m -- number of rows in a */
/*       n -- number of columns in a */
/*       matveca -- routine which applies the adjoint */
/*                  of the matrix to be SVD'd */
/*                  to an arbitrary vector; this routine must have */
/*                  a calling sequence of the form */

/*                  matveca(m,x,n,y,p1t,p2t,p3t,p4t), */

/*                  where m is the length of x, */
/*                  x is the vector to which the adjoint */
/*                  of the matrix is to be applied, */
/*                  n is the length of y, */
/*                  y is the product of the adjoint of the matrix and x, */
/*                  and p1t, p2t, p3t, and p4t are user-specified */
/*                  parameters */
/*       p1t -- parameter to be passed to routine matveca */
/*       p2t -- parameter to be passed to routine matveca */
/*       p3t -- parameter to be passed to routine matveca */
/*       p4t -- parameter to be passed to routine matveca */
/*       matvec -- routine which applies the matrix to be SVD'd */
/*                 to an arbitrary vector; this routine must have */
/*                 a calling sequence of the form */

/*                 matvec(n,x,m,y,p1,p2,p3,p4), */

/*                 where n is the length of x, */
/*                 x is the vector to which the matrix is to be applied, */
/*                 m is the length of y, */
/*                 y is the product of the matrix and x, */
/*                 and p1, p2, p3, and p4 are user-specified parameters */
/*       p1 -- parameter to be passed to routine matvec */
/*       p2 -- parameter to be passed to routine matvec */
/*       p3 -- parameter to be passed to routine matvec */
/*       p4 -- parameter to be passed to routine matvec */

/*       output: */
/*       krank -- rank of the SVD constructed */
/*       iu -- index in w of the first entry of the matrix */
/*             of orthonormal left singular vectors of a */
/*       iv -- index in w of the first entry of the matrix */
/*             of orthonormal right singular vectors of a */
/*       is -- index in w of the first entry of the array */
/*             of singular values of a; the singular values are stored */
/*             as complex*16 numbers whose imaginary parts are zeros */
/*       w -- array containing the singular values and singular vectors */
/*            of a; w doubles as a work array, and so must be at least */
/*            (krank+1)*(3*m+5*n+11)+8*krank^2 complex*16 elements long, */
/*            where krank is the rank returned by the present routine */
/*       ier -- 0 when the routine terminates successfully; */
/*              -1000 when lw is too small; */
/*              other nonzero values when idz_id2svd bombs */

/*       _N.B._: w must be at least (krank+1)*(3*m+5*n+11)+8*krank**2 */
/*               complex*16 elements long, where krank is the rank */
/*               returned by the present routine. Also, the algorithm */
/*               used by the present routine is randomized. */



/*       Allocate some memory. */

    /* Parameter adjustments */
    --w;

    /* Function Body */
    lw2 = 0;

    ilist = lw2 + 1;
    llist = *n;
    lw2 += llist;

    iproj = lw2 + 1;


/*       ID a. */

    lp = *lw - lw2;
    idzp_rid_(&lp, eps, m, n, (U_fp)matveca, p1t, p2t, p3t, p4t, krank, &w[
	    ilist], &w[iproj], ier);
    if (*ier != 0) {
	return 0;
    }


    if (*krank > 0) {


/*         Allocate more memory. */

	lproj = *krank * (*n - *krank);
	lw2 += lproj;

	icol = lw2 + 1;
	lcol = *m * *krank;
	lw2 += lcol;

	iui = lw2 + 1;
	lu = *m * *krank;
	lw2 += lu;

	ivi = lw2 + 1;
	lv = *n * *krank;
	lw2 += lv;

	isi = lw2 + 1;
	ls = *krank;
	lw2 += ls;

	iwork = lw2 + 1;
/* Computing 2nd power */
	i__1 = *krank;
	lwork = (*krank + 1) * (*m + *n * 3 + 10) + i__1 * i__1 * 9;
	lw2 += lwork;


	if (*lw < lw2) {
	    *ier = -1000;
	    return 0;
	}


	idzp_rsvd0_(m, n, (U_fp)matveca, p1t, p2t, p3t, p4t, (U_fp)matvec, 
		p1, p2, p3, p4, krank, &w[iui], &w[ivi], &w[isi], ier, &w[
		ilist], &w[iproj], &w[icol], &w[iwork]);
	if (*ier != 0) {
	    return 0;
	}


	*iu = 1;
	*iv = *iu + lu;
	*is = *iv + lv;


/*         Copy the singular values and singular vectors */
/*         into their proper locations. */

	i__1 = lu;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *iu + k - 1;
	    i__3 = iui + k - 1;
	    w[i__2].r = w[i__3].r, w[i__2].i = w[i__3].i;
	}

/* k */
	i__1 = lv;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *iv + k - 1;
	    i__3 = ivi + k - 1;
	    w[i__2].r = w[i__3].r, w[i__2].i = w[i__3].i;
	}

/* k */
	idz_reco_(&ls, &w[isi], &w[*is]);


    }


/* krank .gt. 0 */
    return 0;
} /* idzp_rsvd__ */

