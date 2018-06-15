import pytest
import shutil
import os

from astropy.io import fits

from ..tastis import tastis

filelist = ['data/tastis/oc7w11viq_raw.fits',
            'data/tastis/ocmv0lw6q_raw.fits',
            'data/tastis/octka7jeq_raw.fits',
            'data/tastis/octr11h4q_raw.fits',
            'data/tastis/octr11hrq_raw.fits',
            'data/tastis/ocu252cmq_raw.fits',
            'data/tastis/ocui04xeq_raw.fits',
            'data/tastis/ocyw05afq_raw.fits']


# Need to test keyword ACQSTAT update
def setup_module():
    for filename in filelist:
        shutil.copy(filename, filename.replace("raw.fits", "copy_raw.fits"))
        shutil.copy(filename.replace("raw.fits", "spt.fits"),
                    filename.replace("raw.fits", "copy_spt.fits"))


def teardown_module():
    for filename in filelist:
        os.remove(filename.replace("raw.fits", "copy_raw.fits"))
        os.remove(filename.replace("raw.fits", "copy_spt.fits"))


def test_header_update1():
    tastis('data/tastis/oc7w11viq_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/oc7w11viq_copy_raw.fits')
    assert head['ACQP_END'] == 'OK_END'
    assert head['ACQP_FLX'] == 'LO_FLUX'
    assert head['ACQP_RAT'] == 'HIRATIO'
    assert head['ACQP_SAT'] == 'UNSAT'
    assert head['ACQSTAT'] == 'FAILED'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update2():
    tastis('data/tastis/octka7jeq_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/octka7jeq_copy_raw.fits')
    assert head['ACQP_END'] == 'OK_END'
    assert head['ACQP_FLX'] == 'OK_FLUX'
    assert head['ACQP_RAT'] == 'LORATIO'
    assert head['ACQP_SAT'] == 'UNSAT'
    assert head['ACQSTAT'] == 'FAILED'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update3():
    tastis('data/tastis/octr11hrq_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/octr11hrq_copy_raw.fits')
    assert head['ACQP_END'] == 'HI_END'
    assert head['ACQP_FLX'] == 'LO_FLUX'
    assert head['ACQP_RAT'] == 'HIRATIO'
    assert head['ACQP_SAT'] == 'UNSAT'
    assert head['ACQSTAT'] == 'FAILED'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update4():
    tastis('data/tastis/ocui04xeq_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/ocui04xeq_copy_raw.fits')
    assert head['ACQP_END'] == 'HI_END'
    assert head['ACQP_FLX'] == 'OK_FLUX'
    assert head['ACQP_RAT'] == 'LORATIO'
    assert head['ACQP_SAT'] == 'UNSAT'
    assert head['ACQSTAT'] == 'FAILED'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update5():
    tastis('data/tastis/ocyw05afq_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/ocyw05afq_copy_raw.fits')
    assert head['ACQP_END'] == 'OK_END'
    assert head['ACQP_FLX'] == 'OK_FLUX'
    assert head['ACQP_RAT'] == 'OKRATIO'
    assert head['ACQP_SAT'] == 'UNSAT'
    assert head['ACQSTAT'] == 'OK'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update6():
    tastis('data/tastis/ocmv0lw6q_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/ocmv0lw6q_copy_raw.fits')
    assert head['ACQSTAT'] == 'OK'
    assert head['ACQ_LAMP'] == 'OK_LAMP'
    assert head['ACQ_RAT'] == 'OKRATIO'
    assert head['ACQ_SAT'] == 'UNSAT'
    assert head['ACQ_SLEW'] == 'OK_SLEW'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update7():
    tastis('data/tastis/octr11h4q_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/octr11h4q_copy_raw.fits')
    assert head['ACQSTAT'] == 'FAILED'
    assert head['ACQ_LAMP'] == 'OK_LAMP'
    assert head['ACQ_RAT'] == 'OKRATIO'
    assert head['ACQ_SAT'] == 'UNSAT'
    assert head['ACQ_SLEW'] == 'BIGSLEW'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_header_update8():
    tastis('data/tastis/ocu252cmq_copy_raw.fits', update=True)
    head = fits.getheader('data/tastis/ocu252cmq_copy_raw.fits')
    assert head['ACQSTAT'] == 'FAILED'
    assert head['ACQ_LAMP'] == 'OK_LAMP'
    assert head['ACQ_RAT'] == 'LORATIO'
    assert head['ACQ_SAT'] == 'UNSAT'
    assert head['ACQ_SLEW'] == 'OK_SLEW'
    assert head['DATAFLAG'] == 'TDF_Up'
    assert head['GROUPS'] is False


def test_tastis_ouput_f25nd3(capsys):
    tastis('data/tastis/ocmv0lw6q_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
       "ocmv0lw6q       HST/STIS    MIRVIS      F25ND3             ACQ/POINT\n" \
       "prop: 13760      visit: 0L    line: 1   target: CD-59D3300\n" \
       "obs date, time: 2016-09-29    23:43:50   exposure time:  1.10\n" \
       "dom GS/FGS: S4B0000993F2    sub-dom GS/FGS: S4B0000953F1\n" \
       "ACQ params:     bias sub: 1510   checkbox: 3      method: FLUX CENTROID\n" \
       "subarray (axis1,axis2):   size=(100,100)          corner=(487,466)\n" \
       "-------------------------------------------------------------------------------\n" \
       "Coarse locate phase:           Target flux in max checkbox (DN): 1560\n" \
       "\n" \
       "                       global          local\n" \
       "                    axis1 axis2     axis1 axis2\n" \
       "Target location:    534.2  507.0    48.2  42.0\n" \
       "\n" \
       "                    axis1 axis2     axis1  axis2         V2      V3\n" \
       "                      (pixels)        (arcsec)            (arcsec)\n" \
       "Estimated slew:     -1.5  -9.0    -0.079 -0.457       -0.379  0.268\n" \
       "-------------------------------------------------------------------------------\n" \
       "Fine locate phase:            Target flux in max checkbox (DN): 1559\n" \
       "\n" \
       "                       global            local\n" \
       "                    axis1 axis2     axis1 axis2\n" \
       "Target location:    534.2  516.8    48.2  51.8\n" \
       "Ref ap location:    537.5  517.0    19.5  17.0\n" \
       "\n" \
       "                    axis1 axis2     axis1  axis2         V2      V3\n" \
       "                      (pixels)        (arcsec)           (arcsec)\n" \
       "Estimated slew:     -2.1  -0.2     -0.104 -0.010      -0.081 -0.067\n" \
       "-------------------------------------------------------------------------------\n" \
       "Total est. slew:    -3.6  -9.2    -0.183 -0.467        -0.460  0.201\n" \
       "-------------------------------------------------------------------------------\n" \
       "Your ACQ appears to have succeeded, as the fluxes in the coarse\n" \
       "and fine stages agree within 25% and the fine slews were less than\n" \
       "4 pixels as expected\n" \
       "\n" \
       "===============================================================================\n"


def test_tastis_zero_divide():
    with pytest.raises(ZeroDivisionError):
        tastis('data/tastis/o4er06llq_raw.fits')


def test_tastis_output_toobright(capsys):
    tastis('data/tastis/oc7w11viq_raw.fits')
    captured = capsys.readouterr()
    print(captured.out)
    assert captured.out == "===============================================================================\n" \
       "oc7w11viq       HST/STIS    G430L        0.3X0.05ND             ACQ/PEAK-UP\n" \
       "prop: 13465      visit: 11    line: 3   target: HD128621-2\n" \
       "obs date, time: 2014-07-24    22:05:06   exposure time:  0.10\n" \
       "dom GS/FGS: S7QX000330F1    sub-dom GS/FGS: S7QX000694F2\n" \
       "ACQ params:     bias sub: 1510                     method: RETURN-TO-BRIGHTEST\n" \
       "subarray (axis1,axis2):   size=(1022,32)          corner=(25,500)\n" \
       "-------------------------------------------------------------------------------\n" \
       "Scan type: LINEARAXIS2                  Step size (mas): 250\n" \
       "\n" \
       " [210 753   0]\n" \
       "\n" \
       "                    axis1 axis2     axis1  axis2         V2   V3\n" \
       "                      (pixels)        (arcsec)           (arcsec)\n" \
       "Estimated slew:      0.0  -0.1     0.000 -0.005     -0.004  0.004\n" \
       "Flux in post-slew confirmation image (751752) - Pedestal (748587) = 3165 DN\n" \
       "-------------------------------------------------------------------------------\n" \
       "The flux in the confirmation image is 320% greater than the maximum flux\n" \
       "in the ACQ/PEAK scan.  An excess greater than 100% indicates\n" \
       "problems in the ACQ/PEAK.\n" \
       "\n" \
       "The flux in the confirmation image is 16% of the recommended minimum\n" \
       "of 20000 DN for a dispersed-light ACQ/PEAK.  The signal-to-noise in\n" \
       "the ACQ/PEAK may be inadequate for an accurate centering.\n" \
       "\n" \
       "===============================================================================\n"


def test_tastis_output_geometric_center(capsys):
    tastis('data/tastis/ocoa03q2q_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "ocoa03q2q       HST/STIS    MIRVIS     F28X50LP             ACQ/DIFFUSE\n" \
          "prop: 13693      visit: 03    line: 1   target: CERES-2\n" \
          "obs date, time: 2015-08-16    16:34:18   exposure time:  1.10\n" \
          "dom GS/FGS: SCHN000911F2    sub-dom GS/FGS: SCHK000675F1\n" \
          "ACQ params:     bias sub: 1510   checkbox: 3      method: GEOMETRIC-CENTER\n" \
          "subarray (axis1,axis2):   size=(104,104)          corner=(485,464)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Coarse locate phase:           Target flux in max checkbox (DN): 87956\n" \
          "\n" \
          "                       global          local\n" \
          "                    axis1 axis2     axis1 axis2\n" \
          "Target location:    528.0  515.0    44.0  52.0\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2      V3\n" \
          "                      (pixels)        (arcsec)            (arcsec)\n" \
          "Estimated slew:     -7.7  -1.0    -0.393 -0.051       -0.314 -0.242\n" \
          "-------------------------------------------------------------------------------\n" \
          "Fine locate phase:            Target flux in max checkbox (DN): 87849\n" \
          "\n" \
          "                       global            local\n" \
          "                    axis1 axis2     axis1 axis2\n" \
          "Target location:    534.0  517.0    50.0  54.0\n" \
          "Ref ap location:    537.3  515.4    21.3  17.4\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2      V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:     -2.0   1.6     -0.103  0.081      -0.015 -0.130\n" \
          "-------------------------------------------------------------------------------\n" \
          "Total est. slew:    -9.8   0.6    -0.496  0.030        -0.329 -0.372\n" \
          "-------------------------------------------------------------------------------\n" \
          "Your ACQ appears to have succeeded, as the fluxes in the coarse\n" \
          "and fine stages agree within 25% and the fine slews were less than\n" \
          "4 pixels as expected\n" \
          "\n" \
          "===============================================================================\n"


def test_tastis_output_lowflux(capsys):
    tastis('data/tastis/octka7jeq_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "octka7jeq       HST/STIS    G430L        0.2X0.09             ACQ/PEAK-UP\n" \
          "prop: 14161      visit: A7    line: 2   target: HD-84937\n" \
          "obs date, time: 2016-05-09    23:15:29   exposure time:  0.20\n" \
          "dom GS/FGS: N6U6000023F2    sub-dom GS/FGS: N6U7000178F1\n" \
          "ACQ params:     bias sub: 1510                     method: MAX-FLUX-CENTROID\n" \
          "subarray (axis1,axis2):   size=(1022,32)          corner=(26,500)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Scan type: LINEARAXIS1                  Step size (mas): 69\n" \
          "\n" \
          " [    0 16309 83580 21884  8029]\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2   V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:      0.2   0.0     0.010  0.000      0.007  0.007\n" \
          "Flux in post-slew confirmation image (852814) - Pedestal (791686) = 61128 DN\n" \
          "-------------------------------------------------------------------------------\n" \
          "The flux in the confirmation image is only 73% of the maximum flux\n" \
          "in the ACQ/PEAK scan.  Percentages below 80% often indicate problems\n" \
          "in the ACQ/PEAK.\n" \
          "\n" \
          "===============================================================================\n"


def test_tastis_output_slewerr(capsys):
    tastis('data/tastis/octr11h4q_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "octr11h4q       HST/STIS    MIRVIS      F25ND5             ACQ/POINT\n" \
          "prop: 14341      visit: 11    line: 1   target: HD128620\n" \
          "obs date, time: 2016-08-28    19:57:49   exposure time:  0.30\n" \
          "dom GS/FGS: S7QX000303F1    sub-dom GS/FGS: S7QX000751F2\n" \
          "ACQ params:     bias sub: 1510   checkbox: 3      method: FLUX CENTROID\n" \
          "subarray (axis1,axis2):   size=(100,100)          corner=(487,466)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Coarse locate phase:           Target flux in max checkbox (DN): 278\n" \
          "\n" \
          "                       global          local\n" \
          "                    axis1 axis2     axis1 axis2\n" \
          "Target location:    557.0  473.0    71.0   8.0\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2      V3\n" \
          "                      (pixels)        (arcsec)            (arcsec)\n" \
          "Estimated slew:     21.3  -43.0     1.080 -2.184       -0.781  2.308\n" \
          "-------------------------------------------------------------------------------\n" \
          "Fine locate phase:            Target flux in max checkbox (DN): 280\n" \
          "\n" \
          "                       global            local\n" \
          "                    axis1 axis2     axis1 axis2\n" \
          "Target location:    547.0  564.0    61.0  99.0\n" \
          "Ref ap location:    537.6  517.3    19.6  17.3\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2      V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:     10.6  46.7      0.541  2.372       2.060 -1.295\n" \
          "-------------------------------------------------------------------------------\n" \
          "Total est. slew:    31.9   3.7     1.621  0.188         1.279  1.013\n" \
          "-------------------------------------------------------------------------------\n" \
          "The fine slew (to center the target in the reference aperture) is larger\n" \
          "than 4 pixels.  This may indicate a problem with your acquisition.\n" \
          "\n" \
          "===============================================================================\n"


def test_tastis_output_highflux(capsys):
    tastis('data/tastis/octr11hrq_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "octr11hrq       HST/STIS    G430M        31X0.05NDA             ACQ/PEAK-UP\n" \
          "prop: 14341      visit: 11    line: 9   target: HD128621-2\n" \
          "obs date, time: 2016-08-28    22:33:14   exposure time:  0.10\n" \
          "dom GS/FGS: S7QX000303F1    sub-dom GS/FGS: S7QX000751F2\n" \
          "ACQ params:     bias sub: 1510                     method: MAX-FLUX-CENTROID\n" \
          "subarray (axis1,axis2):   size=(1022,32)          corner=(25,500)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Scan type: LINEARAXIS1                  Step size (mas): 39\n" \
          "\n" \
          " [5478    0  798 3264 4796 1923 4876]\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2   V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:      0.2   0.0     0.010  0.000      0.007  0.007\n" \
          "Flux in post-slew confirmation image (882661) - Pedestal (871184) = 11477 DN\n" \
          "-------------------------------------------------------------------------------\n" \
          "The flux in the confirmation image is 110% greater than the maximum flux\n" \
          "in the ACQ/PEAK scan.  An excess greater than 100% indicates\n" \
          "problems in the ACQ/PEAK.\n" \
          "\n" \
          "The flux in the confirmation image is 57% of the recommended minimum\n" \
          "of 20000 DN for a dispersed-light ACQ/PEAK.  The signal-to-noise in\n" \
          "the ACQ/PEAK may be inadequate for an accurate centering.\n" \
          "\n" \
          "The maximum flux in the sequence occurred at one end.\n" \
          "This may indicate that the target was beyond that end\n" \
          "or that a neighboring object affected the acquisition.\n" \
          "===============================================================================\n"


def test_tastis_output_diffflux(capsys):
    tastis('data/tastis/ocu252cmq_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "ocu252cmq       HST/STIS    MIRVIS     F28X50OII             ACQ/POINT\n" \
          "prop: 14143      visit: 52    line: 1   target: BD+41-3306\n" \
          "obs date, time: 2016-06-06    08:30:05   exposure time:  2.10\n" \
          "dom GS/FGS: N2JU001340F2    sub-dom GS/FGS: N2K1001229F1\n" \
          "ACQ params:     bias sub: 1510   checkbox: 3      method: FLUX CENTROID\n" \
          "subarray (axis1,axis2):   size=(100,100)          corner=(487,466)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Coarse locate phase:           Target flux in max checkbox (DN): 1442\n" \
          "\n" \
          "                       global          local\n" \
          "                    axis1 axis2     axis1 axis2\n" \
          "Target location:    527.8  513.1    41.8  48.1\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2      V3\n" \
          "                      (pixels)        (arcsec)            (arcsec)\n" \
          "Estimated slew:     -7.9  -2.9    -0.400 -0.147       -0.387 -0.179\n" \
          "-------------------------------------------------------------------------------\n" \
          "Fine locate phase:            Target flux in max checkbox (DN): 611\n" \
          "\n" \
          "                       global            local\n" \
          "                    axis1 axis2     axis1 axis2\n" \
          "Target location:    534.1  516.1    48.1  51.1\n" \
          "Ref ap location:    537.5  516.5    19.5  16.5\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2      V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:     -2.1  -0.4     -0.106 -0.020      -0.089 -0.061\n" \
          "-------------------------------------------------------------------------------\n" \
          "Total est. slew:    -10.0  -3.3    -0.506 -0.168        -0.477 -0.239\n" \
          "-------------------------------------------------------------------------------\n" \
          "The fluxes in the maximum checkbox in the fine and coarse stages differ\n" \
          "by more than 25%.  This may indicate a problem with your acquisition.\n" \
          "\n" \
          "===============================================================================\n"


def test_tastis_output_unbalanced_flux(capsys):
    tastis('data/tastis/ocui04xeq_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "ocui04xeq       HST/STIS    MIRVIS        52X0.1E1             ACQ/PEAK-UP\n" \
          "prop: 14086      visit: 04    line: 2   target: M62-VLA1\n" \
          "obs date, time: 2016-07-22    06:10:30   exposure time: 20.00\n" \
          "dom GS/FGS: S8ES000684F2    sub-dom GS/FGS: S8ES000207F1\n" \
          "ACQ params:     bias sub: 1510                     method: MAX-FLUX-CENTROID\n" \
          "subarray (axis1,axis2):   size=(32,32)          corner=(524,883)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Scan type: LINEARAXIS1                  Step size (mas): 75\n" \
          "\n" \
          " [17007  5446  1717   993     0]\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2   V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:     -2.6   0.0    -0.132  0.000     -0.093 -0.093\n" \
          "Flux in post-slew confirmation image (56705) - Pedestal (43530) = 13175 DN\n" \
          "-------------------------------------------------------------------------------\n" \
          "The flux in the confirmation image is only 77% of the maximum flux\n" \
          "in the ACQ/PEAK scan.  Percentages below 80% often indicate problems\n" \
          "in the ACQ/PEAK.\n" \
          "\n" \
          "The maximum flux in the sequence occurred at one end.\n" \
          "This may indicate that the target was beyond that end\n" \
          "or that a neighboring object affected the acquisition.\n" \
          "===============================================================================\n"



def test_tastis_output_linearaxis2(capsys):
    tastis('data/tastis/ocyw05afq_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
          "ocyw05afq       HST/STIS    G430L        0.2X0.09             ACQ/PEAK-UP\n" \
          "prop: 14084      visit: 05    line: 2   target: BD-11D916\n" \
          "obs date, time: 2016-09-22    08:33:17   exposure time:  1.80\n" \
          "dom GS/FGS: S2AE000156F1    sub-dom GS/FGS: S2AE000086F2\n" \
          "ACQ params:     bias sub: 1510                     method: MAX-FLUX-CENTROID\n" \
          "subarray (axis1,axis2):   size=(1022,32)          corner=(26,500)\n" \
          "-------------------------------------------------------------------------------\n" \
          "Scan type: LINEARAXIS2                  Step size (mas): 150\n" \
          "\n" \
          " [ 5139 67252     0]\n" \
          "\n" \
          "                    axis1 axis2     axis1  axis2         V2   V3\n" \
          "                      (pixels)        (arcsec)           (arcsec)\n" \
          "Estimated slew:      0.0   0.0     0.000  0.000      0.000  0.000\n" \
          "Flux in post-slew confirmation image (907707) - Pedestal (838752) = 68955 DN\n" \
          "-------------------------------------------------------------------------------\n" \
          "The confirmation image has a flux between 0.8 and 2.0 times the\n" \
          "maximum flux in the peakup, which is typical of a successful ACQ/PEAK.\n" \
          "===============================================================================\n"


def test_tastis_output_linearaxis1(capsys):
    tastis('data/tastis/od3v01bfq_raw.fits')
    captured = capsys.readouterr()
    assert captured.out == "===============================================================================\n" \
       "od3v01bfq       HST/STIS    MIRVIS        52X0.05             ACQ/PEAK-UP\n" \
       "prop: 14493      visit: 01    line: 5   target: 2MASS-J23062928-0502285\n" \
       "obs date, time: 2016-09-26    04:17:57   exposure time: 10.00\n" \
       "dom GS/FGS: SB5F000135F1    sub-dom GS/FGS: SB5F000156F2\n" \
       "ACQ params:     bias sub: 1510                     method: MAX-FLUX-CENTROID\n" \
       "subarray (axis1,axis2):   size=(32,32)          corner=(521,500)\n" \
       "-------------------------------------------------------------------------------\n" \
       "Scan type: LINEARAXIS1                  Step size (mas): 26\n" \
       "\n" \
       " [   0  848 3069 3432 1912  555  228]\n" \
       "\n" \
       "                    axis1 axis2     axis1  axis2         V2   V3\n" \
       "                      (pixels)        (arcsec)           (arcsec)\n" \
       "Estimated slew:     -0.2   0.0    -0.010  0.000     -0.007 -0.007\n" \
       "Flux in post-slew confirmation image (40210) - Pedestal (35982) = 4228 DN\n" \
       "-------------------------------------------------------------------------------\n" \
       "The confirmation image has a flux between 0.8 and 2.0 times the\n" \
       "maximum flux in the peakup, which is typical of a successful ACQ/PEAK.\n" \
       "===============================================================================\n"
