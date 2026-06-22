from .resources import BaseSTIS
from functools import partial
import pytest
import numpy as np
from stistools.splice import splice, merge_overlap

approx = partial(pytest.approx, nan_ok=True, rel=1e-14, abs=1e-16)

# DQ flags for testing:
X = 2 | 512  # in SDQFLAGS
y = 1024     # not in SDQFLAGS
z = 4        # in SDQFLAGS by default, but we will exclude it from SDQFLAGS for testing purposes


@pytest.mark.bigdata
@pytest.mark.slow
class TestSplice(BaseSTIS):
    # Test the entire splicing pipeline
    def test_pipeline(self, precision_threshold=1E-6):
        dataset = 'oblh01040_x1d.fits'
        truth = 'spliced_spectrum_truth.dat'
        path_input = self.get_data("splice/input", dataset)
        path_truth = self.get_data("splice/truth", truth)

        spectrum_table = splice(path_input, weight='sensitivity')

        truth = np.loadtxt(path_truth, skiprows=1)
        wl_truth = truth[:, 0]
        f_truth = truth[:, 1]
        u_truth = truth[:, 2]

        f_scale = 1E-12
        wl_diff = abs(np.sum(spectrum_table['WAVELENGTH'].data - wl_truth))
        f_diff = abs(np.sum(spectrum_table['FLUX'].data - f_truth)) / f_scale
        u_diff = abs(np.sum(spectrum_table['ERROR'].data - u_truth)) / f_scale
        diff = wl_diff + f_diff + u_diff

        assert(diff < precision_threshold)


class TestMergeAnalytic:
    '''Test merge_overlap().
    '''
    def setup_method(self, method):
        section1 = {
            'wavelength'  : np.array([ 7., 21., 35., 49., 63., 77.,], dtype=float),
            'net'         : np.array([20., 20., 20., 20., 20., 20.,], dtype=float),
            'flux'        : np.array([ 1.,  1.,  1.,  1.,  1.,  1.,], dtype=float),
            'uncertainty' : np.array([ 5.,  5.,  5.,  5.,  5.,  5.,], dtype=float),
            'data_quality': np.array([ 0,   0,   0,   0,   0,   0, ], dtype=np.uint16),}
        section2 = {
            'wavelength'  : np.array([ 6., 18., 30., 42., 54., 66.,], dtype=float),
            'net'         : np.array([90., 90., 90., 90., 90., 90.,], dtype=float),
            'flux'        : np.array([ 3.,  3.,  3.,  3.,  3.,  3.,], dtype=float),
            'uncertainty' : np.array([10., 10., 10., 10., 10., 10.,], dtype=float),
            'data_quality': np.array([ 0,   0,   0,   0,   0,   0, ], dtype=np.uint16),}
        self.overlap_pair_section = [section1, section2,]

    @pytest.mark.skip
    def test_all_good(self, weight='equal'):
        '''
        Wl 7      21     35     49     63
        |  1   |  1   |  1   |  1   |  1   |
        |      |      |      |      |      |
        
        wl 6     18    30    42    54    66
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |     |     |     |     |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |      |      |      |      |      | 
        
        Combined:
        |  2   |  2   |  2   |  2   |  2   | 
        |      |      |      |      |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,0,0,0,0,], dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality'][:-1] = np.array([0,0,0,0,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['wavelength'], self.overlap_pair_section[0]['wavelength']), 'λ grid'  # higher SNR λ grid
        assert np.array_equal(merged['flux'][:-1],         [2., 2., 2., 2., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  0,  0,  0,  0, ]), 'DQ'
        raise NotImplementedError('uncertainty test')
        #assert np.array_equal(merged['uncertainty'][:-1],  [__________________ ]), 'σ'  # TBD

    @pytest.mark.skip
    def test_all_good_sensitivityswap(self, weight='equal'):
        '''
        Wl 7      21     35     49     63
        |  1   |  1   |  1   |  1   |  1   |
        |      |      |      |      |      |
        
        wl 6     18    30    42    54    66
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |     |     |     |     |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |      |      |      |      |      | 
        
        Combined:
        |  2   |  2   |  2   |  2   |  2   | 
        |      |      |      |      |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,0,0,0,0,], dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality'][:-1] = np.array([0,0,0,0,0,], dtype=np.uint16)

        raise NotImplementedError('Evaluation of sensitivity arrays')

        # Swap the sensitivity values:
        #self.overlap_pair_section[0][''] = np.array([], dtype=float)
        #self.overlap_pair_section[1][''] = np.array([], dtype=float)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['wavelength'], self.overlap_pair_section[1]['wavelength']), 'λ grid'  # higher SNR λ grid
        assert np.array_equal(merged['flux'][:-1],         [2., 2., 2., 2., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  0,  0,  0,  0, ]), 'DQ'

    def test_serious1(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |     |     |     |     |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |      |      |      |      |      | 
        
        Combined:
        |  2   |  3   |  2   |  3   |  2   | 
        |      |      |      |      |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,0,0,0,0,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [2., 3., 2., 3., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  0,  0,  0,  0, ]), 'DQ'

    def test_serious2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |      |      |      |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |     |     |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  X   |      |  X   |  X   | 
        
        Combined:
        |  1   |  1   |  2   |  1   |  1   | 
        |      |      |      |      |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,0,0,0,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,0,0,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., 1., 2., 1., 1.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  0,  0,  0,  0, ]), 'DQ'

    def test_minor1(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |   y  |      |   y  |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |     |     |     |     |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |      |      |      |      |      | 
        
        Combined:
        |  2   |  2   |  2   |  2   |  2   | 
        |      |   y  |      |   y  |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,y,0,y,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,0,0,0,0,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [2., 2., 2., 2., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  y,  0,  y,  0, ]), 'DQ'

    def test_minor2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |      |      |      |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |   y |     |     |   y |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |   y  |   y  |      |   y  |   y  | 
        
        Combined:
        |  2   |  2   |  2   |  2   |  2   | 
        |   y  |   y  |      |   y  |   y  | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,0,0,0,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,y,0,0,y,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [2., 2., 2., 2., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [y,  y,  0,  y,  y, ]), 'DQ'

    def test_serious1_minor1(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  Xy  |      |  Xy  |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |     |     |     |     |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |      |      |      |      |      | 
        
        Combined:
        |  2   |  3   |  2   |  3   |  2   | 
        |      |      |      |      |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X|y,0,X|y,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,0,0,0,0,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [2., 3., 2., 3., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  0,  0,  0,  0, ]), 'DQ'

    def test_serious2_minor1(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |   y  |      |      |   y  |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |     |     |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  X   |      |  X   |  X   | 
        
        Combined:
        |  1   |  1   |  2   |  1   |  1   | 
        |      |   y  |      |      |   y  | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,y,0,0,y,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,0,0,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., 1., 2., 1., 1.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  y,  0,  0,  y, ]), 'DQ'

    def test_serious1_minor2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |   y |     |     |   y |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |   y  |   y  |      |   y  |   y  | 
        
        Combined:
        |  2   |  3   |  2   |  3   |  2   | 
        |   y  |   y  |      |   y  |   y  | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,y,0,0,y,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [2., 3., 2., 3., 2.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [y,  y,  0,  y,  y, ]), 'DQ'

    def test_serious2_minor2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |      |      |      |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  Xy |     |     |  Xy |     |
        
        Interpolate second array to reference
        |  3   |  3   |  3   |  3   |  3   | 
        |  Xy  |  Xy  |      |  Xy  |  Xy  | 
        
        Combined:
        |  1   |  1   |  2   |  1   |  1   | 
        |      |      |      |      |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,0,0,0,0,],       dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X|y,0,0,X|y,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., 1., 2., 1., 1.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0,  0,  0,  0,  0, ]), 'DQ'

    def test_serious12(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |     |     |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  X   |      |  X   |  X   | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  X   |      |  X   |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,0,0,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X, 0, X, 0,]), 'DQ'

    def test_serious12_minor1(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  Xy  |      |  Xy  |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |     |     |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  X   |      |  X   |  X   | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  Xy  |      |  Xy  |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X|y,0,X|y,0,], dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,0,0,X,0,],   dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X|y, 0, X|y, 0,]), 'DQ'

    def test_serious12_minor2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  Xy |     |     |  Xy |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  Xy  |  Xy  |      |  Xy  |  Xy  | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  Xy  |      |  Xy  |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],       dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X|y,0,0,X|y,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X|y, 0, X|y, 0,]), 'DQ'

    def test_serious12_minor12(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  Xy  |      |  Xy  |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  Xz |     |     |  Xz |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  Xz  |  Xz  |      |  Xz  |  Xz  | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  Xyz |      |  Xyz |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X|y,0,X|y,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X|z,0,0,X|z,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight, sdqflags=int(31743 & ~np.uint16(z)))

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X|y|z, 0, X|y|z, 0,]), 'DQ'

    def test_serious12_minor1neighbor(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |   y  |  X   |      |  X   |   y  |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |     |     |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  X   |      |  X   |  X   | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |   y  |  X   |      |  X   |   y  | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([y,X,0,X,y,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,0,0,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [y, X, 0, X, y,]), 'DQ'

    def test_serious12_minor1neighbor2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |   y  |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |     |     |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  X   |      |  X   |  X   | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  X   |   y  |  X   |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,y,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,0,0,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X, y, X, 0,]), 'DQ'

    def test_serious12_minor2neighbor(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |   y |  X  |     |     |  X  |   y |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  Xy  |  X   |      |  X   |  Xy  | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  X   |      |  X   |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([y,X,0,0,X,y,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X, 0, X, 0,]), 'DQ'

    def test_serious12_minor2neighbor2(self, weight='equal'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |   y |   y |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  Xy  |   y  |  Xy  |  X   | 
        
        Combined:
        |  1   | NaN  |  2   | NaN  |  1   | 
        |      |  Xy  |   y  |  Xy  |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,y,y,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert np.array_equal(merged['flux'][:-1],         [1., np.nan, 2., np.nan, 1.,], equal_nan=True), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X|y, y, X|y, 0,]), 'DQ'

    def test_snr(self, weight='snr'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |   y |   y |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  Xy  |   y  |  Xy  |  X   | 
        
        Combined:
        |  1   | NaN  |  1.4 | NaN  |  1   |   <-- 1.4 comes from weighted mean
        |      |  Xy  |   y  |  Xy  |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,y,y,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert merged['flux'][:-1] == approx([1., np.nan, 1.4, np.nan, 1.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X|y, y, X|y, 0,]), 'DQ'

    def test_binary(self, weight='binary'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  Xz  |  y   |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |   y |   y |  X  |     |
        
        Interpolate second array to reference grid (IGNORED):
        |  3   |  3   |  3   |  3   |  3   |
        |  X   |  Xy  |   y  |  Xy  |  X   |
        
        Combined:
        |  1   | NaN  |  1   | NaN  |  1   |
        |      |  X   |      |  Xz  |  y   |
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X|z,y,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,y,y,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert merged['flux'][:-1] == approx([1., np.nan, 1., np.nan, 1.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X, 0, X|z, y,]), 'DQ'

    def test_sensitivity(self, weight='sensitivity-dataset'):
        '''
        |  1   |  1   |  1   |  1   |  1   |
        |      |  X   |      |  X   |      |
        
        |  3  |  3  |  3  |  3  |  3  |  3  |
        |     |  X  |   y |   y |  X  |     |
        
        Interpolate second array to reference grid:
        |  3   |  3   |  3   |  3   |  3   | 
        |  X   |  Xy  |   y  |  Xy  |  X   | 
        
        Combined:
        |  1   | NaN  |  2.2 | NaN  |  1   |   <-- 2.2 comes from NET/FLUX weighting
        |      |  Xy  |   y  |  Xy  |      | 
        '''
        self.overlap_pair_section[0]['data_quality'][:-1] = np.array([0,X,0,X,0,],   dtype=np.uint16)
        self.overlap_pair_section[1]['data_quality']      = np.array([0,X,y,y,X,0,], dtype=np.uint16)

        merged = merge_overlap(self.overlap_pair_section, weight=weight)

        assert merged['flux'][:-1] == approx([1., np.nan, 2.2, np.nan, 1.,]), 'flux'
        assert np.array_equal(merged['data_quality'][:-1], [0, X|y, y, X|y, 0,]), 'DQ'

        raise NotImplementedError('Uncertainties for weight="sensitivity"')





## PROBLEM:  sections[0]['wavelength'][N] and sections[1]['wavelength'][N] are different!
##           Must identify a different N index for the first array to test multiple flag inputs! *** <--- WORKING HERE
#
#class TestMerge:
#    '''Test merge_overlap().
#    '''
#    def setup_method(self, method):
#        self.overlap_pair_sections = get_real_test_data()
#
#    def test_realdata_dq(self, i=-6, weight='snr'):
#        overlap_merged = merge_overlap(self.overlap_pair_sections[i], weight=weight)
#        assert (overlap_merged['data_quality'] == 0).all()
#
#    def test_dq1(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] =   0
#        sections[1]['data_quality'][11] = 514  # In SDQFLAGS; higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] == 0
#        assert overlap_merged['data_quality'][11] == 0
#        assert overlap_merged['data_quality'][12] == 0
#
#        # Also check that output wavelengths are from the second (higher SNR) spectrum:
#        assert not np.array_equal(overlap_merged['wavelength'], sections[0]['wavelength'])  # low SNR
#        assert     np.array_equal(overlap_merged['wavelength'], sections[1]['wavelength'])  # high SNR; target λ
#
#    def test_dq1_serious_and_unserious(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] =   0
#        sections[1]['data_quality'][11] = 514 | 1024  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] == 0
#        #assert overlap_merged['data_quality'][11] == 0 # WRONGLY RETURNING 1024 ***
#        assert overlap_merged['data_quality'][12] == 0
#
#        #assert overlap_merged['flux'][11] == sections[0]['flux'][11]
#
#    def test_dq2(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] = 514  # In SDQFLAGS
#        sections[1]['data_quality'][11] =   0  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] == 0
#        assert overlap_merged['data_quality'][11] == 0
#        assert overlap_merged['data_quality'][12] == 0
#
#    def test_dq2_serious_and_unserious(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] = 514 | 1024
#        sections[1]['data_quality'][11] = 0  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] == 0
#        assert overlap_merged['data_quality'][11] == 0
#        assert overlap_merged['data_quality'][12] == 0
#
#        #assert overlap_merged['flux'][11] == sections[1]['flux'][11]
#
#    def test_dq12_serious_and_unserious(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] = 514
#        sections[1]['data_quality'][11] = 1024  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] ==    0
#        assert overlap_merged['data_quality'][11] == 1024
#        assert overlap_merged['data_quality'][12] ==    0
#
#    def test_dq12_v2_serious_and_unserious(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] = 1024
#        sections[1]['data_quality'][11] =  514  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][ 9] ==    0
#        assert overlap_merged['data_quality'][10] == 1024
#        assert overlap_merged['data_quality'][11] == 1024
#        assert overlap_merged['data_quality'][12] ==    0
#
#    def test_dq12_serious(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] = 514
#        sections[1]['data_quality'][11] = 514  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        #print('wavelength: ', sections[0]['wavelength'][11], sections[1]['wavelength'][11], overlap_merged['wavelength'][11])
#        #assert overlap_merged['data_quality'][10] == 
#        #assert overlap_merged['data_quality'][11] == 
#        #assert overlap_merged['data_quality'][12] == 
#
#        assert np.isnan(overlap_merged['flux'][11])
#
#    def test_dq1_minordq(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] =    0
#        sections[1]['data_quality'][11] = 1024  # Not in SDQFLAGS; higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] ==    0
#        assert overlap_merged['data_quality'][11] == 1024
#        assert overlap_merged['data_quality'][12] ==    0
#
#    def test_dq2_minordq(self, i=-6, weight='snr'):
#        sections = deepcopy(self.overlap_pair_sections[i])
#        sections[0]['data_quality'][11] = 1024  # Not in SDQFLAGS
#        sections[1]['data_quality'][11] =    0  # higher SNR section, so target λ
#
#        overlap_merged = merge_overlap(sections, weight=weight)
#        assert overlap_merged['data_quality'][10] == 1024
#        assert overlap_merged['data_quality'][11] == 1024
#        assert overlap_merged['data_quality'][12] ==    0

