


//_lib
extern PyMODINIT_FUNC PyInit__ccallback_c(void);
extern PyMODINIT_FUNC PyInit__fpumode(void);
extern PyMODINIT_FUNC PyInit_messagestream(void);

//cluster
extern PyMODINIT_FUNC PyInit__vq(void);
extern PyMODINIT_FUNC PyInit__hierarchy(void);
extern PyMODINIT_FUNC PyInit__optimal_leaf_ordering(void);

//fftpack
extern PyMODINIT_FUNC PyInit__fftpack(void);
extern PyMODINIT_FUNC PyInit_convolve(void);

//integrate
extern PyMODINIT_FUNC PyInit__dop(void);
extern PyMODINIT_FUNC PyInit__odepack(void);
extern PyMODINIT_FUNC PyInit__quadpack(void);
extern PyMODINIT_FUNC PyInit_lsoda(void);
extern PyMODINIT_FUNC PyInit_vode(void);

//interpolate
extern PyMODINIT_FUNC PyInit__bspl(void);
extern PyMODINIT_FUNC PyInit__fitpack(void);
extern PyMODINIT_FUNC PyInit__interpolate(void);
extern PyMODINIT_FUNC PyInit__ppoly(void);
extern PyMODINIT_FUNC PyInit_dfitpack(void);
extern PyMODINIT_FUNC PyInit_interpnd(void);

//io_matlab
extern PyMODINIT_FUNC PyInit_mio5_utils(void);
extern PyMODINIT_FUNC PyInit_mio_utils(void);
extern PyMODINIT_FUNC PyInit_streams(void);

//linalg
extern PyMODINIT_FUNC PyInit__decomp_update(void);
extern PyMODINIT_FUNC PyInit__fblas(void);
extern PyMODINIT_FUNC PyInit__flapack(void);
extern PyMODINIT_FUNC PyInit__flinalg(void);
extern PyMODINIT_FUNC PyInit__interpolative(void);
extern PyMODINIT_FUNC PyInit__solve_toeplitz(void);
extern PyMODINIT_FUNC PyInit_cython_blas(void);
extern PyMODINIT_FUNC PyInit_cython_lapack(void);

//ndimage
extern PyMODINIT_FUNC PyInit__nd_image(void);
extern PyMODINIT_FUNC PyInit__ni_label(void);
extern PyMODINIT_FUNC PyInit___odrpack(void);

//optimize
extern PyMODINIT_FUNC PyInit__cobyla(void);
extern PyMODINIT_FUNC PyInit__group_columns(void);
extern PyMODINIT_FUNC PyInit__lbfgsb(void);
extern PyMODINIT_FUNC PyInit__minpack(void);
extern PyMODINIT_FUNC PyInit__nnls(void);
extern PyMODINIT_FUNC PyInit__slsqp(void);
extern PyMODINIT_FUNC PyInit__zeros(void);
extern PyMODINIT_FUNC PyInit_minpack2(void);
extern PyMODINIT_FUNC PyInit_moduleTNC(void);

extern PyMODINIT_FUNC PyInit_givens_elimination(void);
extern PyMODINIT_FUNC PyInit__trlib(void);

//signal
extern PyMODINIT_FUNC PyInit__max_len_seq_inner(void);
extern PyMODINIT_FUNC PyInit__peak_finding_utils(void);
extern PyMODINIT_FUNC PyInit__spectral(void);
extern PyMODINIT_FUNC PyInit__upfirdn_apply(void);
extern PyMODINIT_FUNC PyInit_sigtools(void);
extern PyMODINIT_FUNC PyInit_spline(void);

//sparse_csgraph
extern PyMODINIT_FUNC PyInit__min_spanning_tree(void);
extern PyMODINIT_FUNC PyInit__reordering(void);
extern PyMODINIT_FUNC PyInit__shortest_path(void);
extern PyMODINIT_FUNC PyInit__tools(void);
extern PyMODINIT_FUNC PyInit__traversal(void);
//sparse_linalg
extern PyMODINIT_FUNC PyInit__superlu(void);
extern PyMODINIT_FUNC PyInit__arpack(void);
extern PyMODINIT_FUNC PyInit__iterative(void);
//sparse
extern PyMODINIT_FUNC PyInit__csparsetools(void);
extern PyMODINIT_FUNC PyInit__sparsetools(void);

//spatial
extern PyMODINIT_FUNC PyInit__distance_wrap(void);
extern PyMODINIT_FUNC PyInit__hausdorff(void);
extern PyMODINIT_FUNC PyInit__voronoi(void);
extern PyMODINIT_FUNC PyInit_ckdtree(void);
extern PyMODINIT_FUNC PyInit_qhull(void);

//special
extern PyMODINIT_FUNC PyInit__comb(void);
extern PyMODINIT_FUNC PyInit__ellip_harm_2(void);
extern PyMODINIT_FUNC PyInit__ufuncs(void);
extern PyMODINIT_FUNC PyInit__ufuncs_cxx(void);
extern PyMODINIT_FUNC PyInit_cython_special(void);
extern PyMODINIT_FUNC PyInit_specfun(void);

//stats
extern PyMODINIT_FUNC PyInit__stats(void);
extern PyMODINIT_FUNC PyInit_mvn(void);
extern PyMODINIT_FUNC PyInit_statlib(void);

void init_scipy(){
    
    PyImport_AppendInittab("__scipy__lib__ccallback_c", &PyInit__ccallback_c);
    PyImport_AppendInittab("__scipy__lib__fpumode", &PyInit__fpumode);
    PyImport_AppendInittab("__scipy__lib_messagestream", &PyInit_messagestream);
    
    PyImport_AppendInittab("__scipy_cluster__vq", &PyInit__vq);
    PyImport_AppendInittab("__scipy_cluster__hierarchy", &PyInit__hierarchy);
    PyImport_AppendInittab("__scipy_cluster__optimal_leaf_ordering", &PyInit__optimal_leaf_ordering);
    
    PyImport_AppendInittab("__scipy_fftpack__fftpack", &PyInit__fftpack);
    PyImport_AppendInittab("__scipy_fftpack_convolve", &PyInit_convolve);
    
    PyImport_AppendInittab("__scipy_integrate__dop", &PyInit__dop);
    PyImport_AppendInittab("__scipy_integrate__odepack", &PyInit__odepack);
    PyImport_AppendInittab("__scipy_integrate__quadpack", &PyInit__quadpack);
    PyImport_AppendInittab("__scipy_integrate_lsoda", &PyInit_lsoda);
    PyImport_AppendInittab("__scipy_integrate_vode", &PyInit_vode);
    
    PyImport_AppendInittab("__scipy_interpolate__bspl", &PyInit__bspl);
    PyImport_AppendInittab("__scipy_interpolate__fitpack", &PyInit__fitpack);
    PyImport_AppendInittab("__scipy_interpolate__interpolate", &PyInit__interpolate);
    PyImport_AppendInittab("__scipy_interpolate__ppoly", &PyInit__ppoly);
    PyImport_AppendInittab("__scipy_interpolate_dfitpack", &PyInit_dfitpack);
    PyImport_AppendInittab("__scipy_interpolate_interpnd", &PyInit_interpnd);
    
    PyImport_AppendInittab("__scipy_io_matlab_mio5_utils", &PyInit_mio5_utils);
    PyImport_AppendInittab("__scipy_io_matlab_mio_utils", &PyInit_mio_utils);
    PyImport_AppendInittab("__scipy_io_matlab_streams", &PyInit_streams);
    
    PyImport_AppendInittab("__scipy_linalg__decomp_update", &PyInit__decomp_update);
    PyImport_AppendInittab("__scipy_linalg__fblas", &PyInit__fblas);
    PyImport_AppendInittab("__scipy_linalg__flapack", &PyInit__flapack);
    PyImport_AppendInittab("__scipy_linalg__flinalg", &PyInit__flinalg);
    PyImport_AppendInittab("__scipy_linalg__interpolative", &PyInit__interpolative);
    PyImport_AppendInittab("__scipy_linalg__solve_toeplitz", &PyInit__solve_toeplitz);
    PyImport_AppendInittab("__scipy_linalg_cython_blas", &PyInit_cython_blas);
    PyImport_AppendInittab("__scipy_linalg_cython_lapack", &PyInit_cython_lapack);
    
    PyImport_AppendInittab("__scipy_ndimage__nd_image", &PyInit__nd_image);
    PyImport_AppendInittab("__scipy_ndimage__ni_label", &PyInit__ni_label);
    PyImport_AppendInittab("__scipy_odr___odrpack", &PyInit___odrpack);
    
    PyImport_AppendInittab("__scipy_optimize__cobyla", &PyInit__cobyla);
    PyImport_AppendInittab("__scipy_optimize__group_columns", &PyInit__group_columns);
    PyImport_AppendInittab("__scipy_optimize__lbfgsb", &PyInit__lbfgsb);
    PyImport_AppendInittab("__scipy_optimize__minpack", &PyInit__minpack);
    PyImport_AppendInittab("__scipy_optimize__nnls", &PyInit__nnls);
    PyImport_AppendInittab("__scipy_optimize__slsqp", &PyInit__slsqp);
    PyImport_AppendInittab("__scipy_optimize__zeros", &PyInit__zeros);
    PyImport_AppendInittab("__scipy_optimize_minpack2", &PyInit_minpack2);
    PyImport_AppendInittab("__scipy_optimize_moduleTNC", &PyInit_moduleTNC);
    PyImport_AppendInittab("__scipy_optimize__lsq_givens_elimination", &PyInit_givens_elimination);
    PyImport_AppendInittab("__scipy_optimize__trlib__trlib", &PyInit__trlib);
    
    PyImport_AppendInittab("__scipy_signal__max_len_seq_inner", &PyInit__max_len_seq_inner);
    PyImport_AppendInittab("__scipy_signal__peak_finding_utils", &PyInit__peak_finding_utils);
    PyImport_AppendInittab("__scipy_signal__spectral", &PyInit__spectral);
    PyImport_AppendInittab("__scipy_signal__upfirdn_apply",  &PyInit__upfirdn_apply);
    PyImport_AppendInittab("__scipy_signal_sigtools", &PyInit_sigtools);
    PyImport_AppendInittab("__scipy_signal_spline", &PyInit_spline);
    
    PyImport_AppendInittab("__scipy_sparse_csgraph__min_spanning_tree", &PyInit__min_spanning_tree);
    PyImport_AppendInittab("__scipy_sparse_csgraph__reordering", &PyInit__reordering);
    PyImport_AppendInittab("__scipy_sparse_csgraph__shortest_path", &PyInit__shortest_path);
    PyImport_AppendInittab("__scipy_sparse_csgraph__tools", &PyInit__tools);
    PyImport_AppendInittab("__scipy_sparse_csgraph__traversal", &PyInit__traversal);
    PyImport_AppendInittab("__scipy_sparse_linalg_dsolve__superlu", &PyInit__superlu);
    PyImport_AppendInittab("__scipy_sparse_linalg_eigen_arpack__arpack", &PyInit__arpack);
    PyImport_AppendInittab("__scipy_sparse_linalg_isolve__iterative", &PyInit__iterative);
    PyImport_AppendInittab("__scipy_sparse__csparsetools", &PyInit__csparsetools);
    PyImport_AppendInittab("__scipy_sparse__sparsetools", &PyInit__sparsetools);
    
    PyImport_AppendInittab("__scipy_spatial__distance_wrap", &PyInit__distance_wrap);
    PyImport_AppendInittab("__scipy_spatial__hausdorff", &PyInit__hausdorff);
    PyImport_AppendInittab("__scipy_spatial__voronoi", &PyInit__voronoi);
    PyImport_AppendInittab("__scipy_spatial_ckdtree", &PyInit_ckdtree);
    PyImport_AppendInittab("__scipy_spatial_qhull", &PyInit_qhull);
    
    PyImport_AppendInittab("__scipy_special__comb", &PyInit__comb);
    PyImport_AppendInittab("__scipy_special__ellip_harm_2", &PyInit__ellip_harm_2);
    PyImport_AppendInittab("__scipy_special__ufuncs_cxx", &PyInit__ufuncs_cxx);
    PyImport_AppendInittab("__scipy_special__ufuncs", &PyInit__ufuncs);
    PyImport_AppendInittab("__scipy_special_cython_special", &PyInit_cython_special);
    PyImport_AppendInittab("__scipy_special_specfun", &PyInit_specfun);
    
    PyImport_AppendInittab("__scipy_stats__stats", &PyInit__stats);
    PyImport_AppendInittab("__scipy_stats_mvn", &PyInit_mvn);
    PyImport_AppendInittab("__scipy_stats_statlib", &PyInit_statlib);
}
