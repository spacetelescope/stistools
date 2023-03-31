from .resources import BaseSTIS
import pytest
import numpy as np
from stistools.splice import splice


@pytest.mark.bigdata
@pytest.mark.slow
class TestSplice(BaseSTIS):
    # Test the entire splicing pipeline
    def test_pipeline(self, precision_threshold=1E-6):
        dataset = 'oblh01040_x1d.fits'
        truth = 'spliced_spectrum_truth.dat'
        path_input = self.get_data("splice/input", dataset)
        path_truth = self.get_data("splice/truth", truth)

        spectrum_table = splice(path_input, weight='snr')

        truth = np.loadtxt(path_truth)
        wl_truth = truth[:, 0]
        f_truth = truth[:, 1]
        u_truth = truth[:, 2]

        f_scale = 1E-12
        wl_diff = abs(np.sum(spectrum_table['WAVELENGTH'].data - wl_truth))
        f_diff = abs(np.sum(spectrum_table['FLUX'].data - f_truth)) / f_scale
        u_diff = abs(np.sum(spectrum_table['ERROR'].data - u_truth)) / f_scale
        diff = wl_diff + f_diff + u_diff

        assert(diff < precision_threshold)
