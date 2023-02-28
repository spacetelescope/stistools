from .resources import BaseSTIS
import pytest
from stistools.splice import splice_pipeline, nearest_index


@pytest.mark.bigdata
@pytest.mark.slow
class TestSplice(BaseSTIS):
    # Test the entire splicing pipeline
    def test_pipeline(self, precision_threshold=1E-6):
        dataset = 'oblh01040_x1d.fits'
        path = self.get_data("splice/input", dataset)

        spectrum_table = splice_pipeline(path)
        i0 = nearest_index(spectrum_table['WAVELENGTH'].data, 2310)
        test = spectrum_table['FLUX'][i0]
        assert(abs(test - 2.4648798E-12) / test < precision_threshold)
