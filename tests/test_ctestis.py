import numpy as np
from stistools.ctestis import ctestis
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestCtestis(BaseSTIS):
    """
    We need to add more tests for this
    """

    input_loc = 'ctestis'
    ref_loc = 'ctestis'

    def test_single_value_ctestis(self, capsys):
        """
        Test ctestis with a single input set.
        """

        # Prepare input files.
        self.get_data("input", "o4qp9g010_crj.fits")
        self.get_data("input", "o4qp9g010_spt.fits")

        capsys.readouterr()

        # Run ctestis
        fluxc, dmagc, dyc = ctestis(182., 5000., 150.,
                                    stisimage='o4qp9g010_crj.fits')

        captured = capsys.readouterr()
        assert captured.out == "\n" \
                               "mjd: 50893.30\n" \
                               "nread: 2     \n" \
                               "ybin: 1     \n" \
                               "gain:  1.0\n" \
                               "amp: D\n" \
                               "\n" \
                               "tt0: -2.3865942\n" \
                               "lcts: -0.67595399\n" \
                               "bck: 75.0\n" \
                               "lbck: 2.317577\n" \
                               "cti: 1.7314006e-05\n" \
                               "fluxc: 2536.7133\n" \
                               "dmagc: -0.015828427\n" \
                               "cti10000: 0.17314006\n" \
                               "dy512: 0.0043051192\n" \
                               "dyc: 0.007079903\n" \
                               "\n" \
                               "net: 5000.0\n" \
                               "sky: 150.0\n" \
                               "ycol: 182.0\n" \
                               "fluxc: 2536.7133\n" \
                               "dmagc: -0.015828427\n" \
                               "dyc: 0.007079903\n" \
                               "\n"

        assert np.allclose([2536.7133], [fluxc])
        assert np.allclose([-0.015828427], [dmagc])
        assert np.allclose([0.007079903], [dyc])

    def test_list_values_ctestis(self, capsys):
        """
        Test ctestis with a list of values.
        """

        # Prepare input files.
        self.get_data("input", "o4qp9g010_crj.fits")
        self.get_data("input", "o4qp9g010_spt.fits")

        capsys.readouterr()

        # Run ctestis
        fluxc, dmagc, dyc = ctestis([182., 182.], [5000., 1000.], [150., 150.],
                                    stisimage='o4qp9g010_crj.fits')

        captured = capsys.readouterr()
        assert captured.out == "\n" \
                               "mjd: 50893.30\n" \
                               "nread: 2     \n" \
                               "ybin: 1     \n" \
                               "gain:  1.0\n" \
                               "amp: D\n" \
                               "\n" \
                               "tt0: -2.3865942\n" \
                               "lcts: [-0.67595399 -2.2853919 ]\n" \
                               "bck: [75. 75.]\n" \
                               "lbck: [2.31757699 2.31757699]\n" \
                               "cti: [1.73140064e-05 2.15163301e-05]\n" \
                               "fluxc: [2536.71326136  509.14102614]\n" \
                               "dmagc: [-0.01582843 -0.01967022]\n" \
                               "cti10000: [0.17314006 0.2151633 ]\n" \
                               "dy512: [0.00430512 0.00534297]\n" \
                               "dyc: [0.0070799  0.00878668]\n" \
                               "\n" \
                               "net: [5000. 1000.]\n" \
                               "sky: [150. 150.]\n" \
                               "ycol: [182. 182.]\n" \
                               "fluxc: [2536.71326136  509.14102614]\n" \
                               "dmagc: [-0.01582843 -0.01967022]\n" \
                               "dyc: [0.0070799  0.00878668]\n" \
                               "\n"

        assert np.allclose([2536.71326136, 509.14102614], [fluxc])
        assert np.allclose([-0.01582843, -0.01967022], [dmagc])
        assert np.allclose([0.0070799, 0.00878668], [dyc])
