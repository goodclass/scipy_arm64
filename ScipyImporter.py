
# MARK: - scipy
class ScipyImporter(object):
    
    __is_importing__ = False
    
    def find_module(self, fullname, mpath=None):
        if fullname in ('scipy.odr.__odrpack', 'scipy.cluster._optimal_leaf_ordering', 'scipy.cluster._vq', 'scipy.cluster._hierarchy', 'scipy.ndimage._ni_label', 'scipy.ndimage._nd_image', 'scipy.ndimage._ctest_oldapi', 'scipy.ndimage._cytest', 'scipy.ndimage._ctest', 'scipy.linalg._solve_toeplitz', 'scipy.linalg._flinalg', 'scipy.linalg._decomp_update', 'scipy.linalg._interpolative', 'scipy.linalg.cython_blas', 'scipy.linalg._flapack', 'scipy.linalg._fblas', 'scipy.linalg.cython_lapack', 'scipy.optimize._zeros', 'scipy.optimize._minpack', 'scipy.optimize._trlib._trlib', 'scipy.optimize._slsqp', 'scipy.optimize._group_columns', 'scipy.optimize._cobyla', 'scipy.optimize._lsq.givens_elimination', 'scipy.optimize.minpack2', 'scipy.optimize._lbfgsb', 'scipy.optimize._nnls', 'scipy.optimize.moduleTNC', 'scipy.integrate._odepack', 'scipy.integrate._test_multivariate', 'scipy.integrate._test_odeint_banded', 'scipy.integrate.lsoda', 'scipy.integrate.vode', 'scipy.integrate._quadpack', 'scipy.integrate._dop', 'scipy.io.matlab.streams', 'scipy.io.matlab.mio5_utils', 'scipy.io.matlab.mio_utils', 'scipy.io._test_fortran', 'scipy._lib._fpumode', 'scipy._lib._ccallback_c', 'scipy._lib.messagestream', 'scipy._lib._test_ccallback', 'scipy.special._comb', 'scipy.special.cython_special', 'scipy.special._ufuncs', 'scipy.special._test_round', 'scipy.special.specfun', 'scipy.special._ufuncs_cxx', 'scipy.special._ellip_harm_2', 'scipy.fftpack._fftpack', 'scipy.fftpack.convolve', 'scipy.interpolate.dfitpack', 'scipy.interpolate._bspl', 'scipy.interpolate._ppoly', 'scipy.interpolate.interpnd', 'scipy.interpolate._fitpack', 'scipy.interpolate._interpolate', 'scipy.sparse.linalg.isolve._iterative', 'scipy.sparse.linalg.eigen.arpack._arpack', 'scipy.sparse.linalg.dsolve._superlu', 'scipy.sparse._sparsetools', 'scipy.sparse.csgraph._reordering', 'scipy.sparse.csgraph._min_spanning_tree', 'scipy.sparse.csgraph._tools', 'scipy.sparse.csgraph._traversal', 'scipy.sparse.csgraph._shortest_path', 'scipy.sparse._csparsetools', 'scipy.spatial.qhull', 'scipy.spatial._voronoi', 'scipy.spatial._hausdorff', 'scipy.spatial.ckdtree', 'scipy.spatial._distance_wrap', 'scipy.signal._upfirdn_apply', 'scipy.signal.sigtools', 'scipy.signal._peak_finding_utils', 'scipy.signal._spectral', 'scipy.signal.spline', 'scipy.signal._max_len_seq_inner', 'scipy.stats.statlib', 'scipy.stats.mvn', 'scipy.stats._stats'):
            return self
        
        if fullname == 'scipy' and not self.__is_importing__:
            return self
        return
    
    def load_module(self, fullname):
        f = fullname
        if f != 'scipy':
            f = '__' + fullname.replace('.', '_')
        mod = sys.modules.get(f)
        
        if mod is None:
            def importMod():
                mod = importlib.__import__(f)
                sys.modules[fullname] = mod
            
            if fullname != 'scipy':
                importMod()
            else:
                try:
                    self.__is_importing__ = True
                    importMod()
                    self.__is_importing__ = False
                except KeyboardInterrupt:
                    pass
                except SystemExit:
                    pass
                except Exception as e:
                    print(e)
                    report_error('scipy', traceback.format_exc())
                    raise
                finally:
                    self.__is_importing__ = False
            
            return mod
        return mod
