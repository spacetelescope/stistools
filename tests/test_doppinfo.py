from stistools.doppinfo import Doppinfo
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestDoppinfo(BaseSTIS):

    input_loc = 'doppinfo'
    ref_loc = 'doppinfo'

    def test_doppinfo_basic(self, capsys):
        """
        Test stis.doppnfo.Doppinfo with defaults, no updating

        Data info: CENWAVE -- 2263.000
                   FILTER -- E230H
                   ocb6o2020 # IMAGE
        """

        # Prepare input files.
        self.get_data("input", "ocb6o2020_raw.fits")
        self.get_data("input", "ocb6o2020_spt.fits")

        capsys.readouterr()

        # Run doppinfo
        Doppinfo("ocb6o2020_raw.fits", dt=100, spt="ocb6o2020_spt.fits")

        # Compare results
        captured = capsys.readouterr()
        assert captured.out == "# orbitper  doppzero      doppmag      doppmag_v    file\n" \
                               "  5728.67  56752.114170  11.68643135   7.40391177   ocb6o2020_raw.fits[sci,1]\n" \
                               "# time (MJD)   shift   radvel\n" \
                               "56752.165175  -11.59   -7.345\n" \
                               "56752.166333  -11.37   -7.203\n" \
                               "56752.167490  -11.01   -6.975\n" \
                               "56752.168647  -10.52   -6.663\n" \
                               "56752.169805   -9.90   -6.272\n" \
                               "56752.170962   -9.16   -5.805\n" \
                               "56752.172120   -8.32   -5.269\n" \
                               "56752.173277   -7.37   -4.669\n" \
                               "56752.174434   -6.34   -4.014\n" \
                               "56752.175592   -5.23   -3.311\n" \
                               "56752.176749   -4.05   -2.568\n" \
                               "56752.177907   -2.83   -1.794\n" \
                               "56752.179064   -1.58   -0.998\n" \
                               "56752.180222   -0.30   -0.190\n" \
                               "\n" \
                               "# orbitper  doppzero      doppmag      doppmag_v    file\n" \
                               "  5728.67  56752.180505  11.68734454   7.40449032   ocb6o2020_raw.fits[sci,2]\n" \
                               "# time (MJD)   shift   radvel\n" \
                               "56752.181784    1.42    0.902\n" \
                               "56752.182941    2.68    1.700\n" \
                               "56752.184099    3.91    2.477\n" \
                               "56752.185256    5.09    3.225\n" \
                               "56752.186413    6.21    3.935\n" \
                               "56752.187571    7.26    4.598\n" \
                               "56752.188728    8.22    5.205\n" \
                               "56752.189886    9.08    5.750\n" \
                               "56752.191043    9.83    6.227\n" \
                               "56752.192200   10.46    6.628\n" \
                               "56752.193358   10.97    6.950\n" \
                               "56752.194515   11.35    7.189\n" \
                               "56752.195673   11.59    7.342\n" \
                               "56752.196830   11.69    7.406\n" \
                               "\n"

    def test_doppinfo_update(self, capsys):
            """
            This tests stis.doppinfo.Doppinfo with updating.
            Data info: CENWAVE -- 7283.000
                       FILTER -- G750M
            """

            # Prepare input files.
            self.get_data("input", "oac6010a0_raw.fits")
            self.get_data("input", "oac6010a0_spt.fits")

            capsys.readouterr()

            # Run doppinfo
            Doppinfo("oac6010a0_raw.fits", dt=2, spt="oac6010a0_spt.fits",
                     update=True)

            # Compare results
            captured = capsys.readouterr()
            assert captured.out == "# orbitper  doppzero      doppmag      doppmag_v    file\n" \
                                   "  5749.288  55011.984418  0.14207078   3.23985005   oac6010a0_raw.fits[sci,1]\n" \
                                   "# time (MJD)   shift   radvel\n" \
                                   "55012.034045   -0.14   -3.241\n" \
                                   "55012.034068   -0.14   -3.241\n" \
                                   "55012.034092   -0.14   -3.242\n" \
                                   "55012.034115   -0.14   -3.242\n" \
                                   "55012.034138   -0.14   -3.242\n" \
                                   "55012.034161   -0.14   -3.242\n" \
                                   "55012.034184   -0.14   -3.242\n" \
                                   "55012.034207   -0.14   -3.242\n" \
                                   "55012.034230   -0.14   -3.242\n" \
                                   "55012.034254   -0.14   -3.242\n" \
                                   "\n" \
                                   "oac6010a0_raw.fits[sci,1] has been updated as follows:\n" \
                                   "orbitper:  5749.2879 (added)\n" \
                                   "doppzero:  55011.9844183 (added)\n" \
                                   "doppmag:   0.142071 (added)\n" \
                                   "doppmagv:  3.239850 (added)\n" \
                                   "\n" \
                                   "# orbitper  doppzero      doppmag      doppmag_v    file\n" \
                                   "  5749.288  55011.984418  0.14208142   3.24009286   oac6010a0_raw.fits[sci,2]\n" \
                                   "# time (MJD)   shift   radvel\n" \
                                   "55012.034393   -0.14   -3.242\n" \
                                   "55012.034416   -0.14   -3.242\n" \
                                   "55012.034439   -0.14   -3.242\n" \
                                   "55012.034462   -0.14   -3.242\n" \
                                   "55012.034485   -0.14   -3.242\n" \
                                   "55012.034508   -0.14   -3.242\n" \
                                   "55012.034532   -0.14   -3.242\n" \
                                   "55012.034555   -0.14   -3.242\n" \
                                   "55012.034578   -0.14   -3.242\n" \
                                   "55012.034601   -0.14   -3.241\n" \
                                   "\n" \
                                   "oac6010a0_raw.fits[sci,2] has been updated as follows:\n" \
                                   "orbitper:  5749.2879 (added)\n" \
                                   "doppzero:  55011.9844183 (added)\n" \
                                   "doppmag:   0.142081 (added)\n" \
                                   "doppmagv:  3.240093 (added)\n" \
                                   "\n"

            outputs = [("oac6010a0_raw.fits", "oac6010a0_raw_ref.fits")]
            self.compare_outputs(outputs)
