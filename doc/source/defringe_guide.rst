.. code:: ipython3

    import stistools
    from astropy.io import fits

Defringing STIS CCD Spectra With The Stistools Defringing Tool Suite
====================================================================

STIS CCD spectra taken at wavelengths >7000\ :math:`\mathring A` (the
G750M and G750L modes) are impacted by “fringing”, a phenomenon caused
by interference of multiple reflections between the two surfaces of the
CCD in cases where the wavelength of the incident light is a small
integer multiple of the distance between the two surfaces of the CCD
(Goudfrooij et al. 1998). Stistools contains four tools used for
correcting the fringing effect, ``normspflat``, ``prepspec``,
``mkfringeflat``, and ``defringe``, which are located in the
``stistools.defringe`` package. These tools are python ports of the
original IRAF/PyRAF ``stis`` defringing tools. The following serves as a
guide for using these tools to defringe STIS CCD data, and presents
several practical examples for common use cases.

This guide assumes that the ``stistools`` package has been installed.
For instructions on how to install ``stistools``, consult the Getting
Started Guide:
https://stistools.readthedocs.io/en/latest/gettingstarted.html.
Additionally, it assumes you’ve setup CRDS which is also walked through
in the Getting Started Guide.

Data Setup
----------

To defringe a STIS CCD G750M/G750L observation, the fits files for the
fringed science image in addition to the associated contemporaneous
fringe flat for the observation are required. It’s standard practice to
take contemporaneous fringe flats for G750M/G750L spectra, the filename
for the fringe flat is stored in the primary header of the science
image, specifically in the ``'FRNGFLAT'`` keyword. Additionally, the
``'OPT_ELEM'`` keyword in the primary header contains the mode used for
the observation, which will be important as there are a few differences
in the defringing process depending on the mode. These differences will
be highlighted in each step.

.. code:: ipython3

    #setup data paths
    sci_file = "odvkl3050" # the science file
    flat_file = "odvkl3080" # the associated fringe flat
    
    #Confirm that the flat file is indeed the associated fringe flat for odvkl3050
    print("Associated Fringe Flat: "+fits.getheader(f"{sci_file}_raw.fits",0)['FRNGFLAT'])
    print("Observing Mode: "+fits.getheader(f"{sci_file}_raw.fits",0)['OPT_ELEM'])


.. parsed-literal::

    Associated Fringe Flat: ODVKL3080
    Observing Mode: G750L


1. Normalize the Fringe Flat
----------------------------

The ``normspflat`` tool is used to normalize the contemporaneous flat
field image. The ``do_cal`` parameter tells ``normspflat`` whether to
calibrate the flat file using ``calstis`` or not. If using a raw input
image, ``do_cal`` should be set to true, if using a calibrated image
(``crj`` or ``sx2``) ``do_cal`` can be set to false. The output file,
which has the nsp identifier in our example, is the normalized flat
field for use in the next steps.

.. code:: ipython3

    stistools.defringe.normspflat(f"{flat_file}_raw.fits", 
                                  f"{flat_file}_nsp.fits", do_cal=True, 
                                  wavecal=f"{sci_file}_wav.fits")


.. parsed-literal::

    File written:  /Users/stisuser/data/path/odvkl3080_crj.fits


2. Prepare the Science File for the Defringing Correction (Optional)
--------------------------------------------------------------------

The ``prepspec`` tool is used to calibrate the raw science image. It
essentially runs the science through calstis with a specific set of
calibration flags. ``prepspec`` will not overwrite image products, so be
sure to remove any higher level science products (``flt``, ``crj``,
``sx1``, ``sx2``) from the working directory before you run it.

**Note:** Running ``prepspec`` is an optional step, if you already have
calibrated science data (``crj``/``sx2``) then running ``prepspec`` is
not essential. However, you may wish to delete these data products and
rebuild from the raw science file with ``prepspec`` to ensure that the
correct calibration was done on the files.

.. code:: ipython3

    %%capture 
    #capture the long calstis output
    
    stistools.defringe.prepspec(f"{sci_file}_raw.fits")

3. Match Fringes in the Fringe Flat Field and the Science Spectra
-----------------------------------------------------------------

The ``mkfringeflat`` tool is used to calculate the appropriate shifts
and scale factors needed to match the fringes in the fringe flat and the
science spectra. The best shift and scale factors are obtained by
finding the values that minimize the RMS within a user-specified search
range. The parameters that control the range and step size for the shift
and scale have default values (shown explicitly below) that should serve
most use cases well. ``mkfringeflat`` will warn the user if the best
shift and scale values were found at the edge of the range, suggesting
the range may need to be expanded further to find the best values.

**G750M/G750L Point of Difference:** The appropriate file type to use as
the input science file depends on the observation mode. For G750L,
``crj`` files should be used. For G750M, geometric correction is
required before defringing can take place, so ``sx2`` products should be
used.

.. code:: ipython3

    # choose the correct science product type based on the mode
    mode = fits.getheader(f"{sci_file}_raw.fits",0)['OPT_ELEM']
    if mode == "G750L":
        prod_type = "crj"
    elif mode == "G750M":
        prod_type = "sx2"
    
    stistools.defringe.mkfringeflat(f"{sci_file}_{prod_type}.fits", f"{flat_file}_nsp.fits", 
                                    f"{flat_file}_frr.fits", beg_shift=-0.5, end_shift=0.5, shift_step=0.1, 
                                    beg_scale=0.8, end_scale=1.2, scale_step=0.04)


.. parsed-literal::

    mkfringeflat.py version 0.1
     - matching fringes in a flatfield to those in science data
     Extraction center: row 583
       Extraction size: 11.0 pixels  [Aperture: 52X2]
    Range to be normalized: [578:589,4:1020]
    
    Determining best shift for fringe flat
    
    shift = -0.5, rms = 8.86830499443427
    shift = -0.4, rms = 9.537416205896829
    shift = -0.3, rms = 10.313446231610731
    shift = -0.19999999999999996, rms = 10.859321697580842
    shift = -0.09999999999999998, rms = 12.815115948375244
    shift = 0.0, rms = 2.865668927113401
    shift = 0.10000000000000009, rms = 2.930035418866309
    shift = 0.20000000000000007, rms = 2.9325634070997313
    shift = 0.30000000000000004, rms = 3.0000779821345067
    shift = 0.4, rms = 3.048857113598449
    shift = 0.5, rms = 3.099803105862569
     
     Best shift : 0.08832022657826488 pixels
     Shifted flat : odvkl3080_nsp_sh.fits
                    (Can be used as input flat for next iteration)
    
    Determining best scaling of amplitude of fringes in flat
    
    Fringes scaled  0.8: RMS = 3.465684888742133
    Fringes scaled  0.8400000000000001: RMS = 3.293885574776519
    Fringes scaled  0.88: RMS = 3.173002670273641
    Fringes scaled  0.92: RMS = 3.127387585002071
    Fringes scaled  0.9600000000000001: RMS = 2.970040264399544
    Fringes scaled  1.0: RMS = 2.924781030708238
    Fringes scaled  1.04: RMS = 12.653574547469535
    Fringes scaled  1.08: RMS = 10.852646024830669
    Fringes scaled  1.12: RMS = 9.53614375927818
    Fringes scaled  1.1600000000000001: RMS = 8.531799668011566
    Fringes scaled  1.2000000000000002: RMS = 7.740201743078294
     
     Best scale : 0.9660612672342681
    Output flat : odvkl3080_frr.fits
      (to be used as input to task 'defringe.py')


4. Defringe the Science Spectra
-------------------------------

The final step is to use the ``defringe`` tool to divide the scaled and
shifted fringe flat off of the calibrated science spectra, removing the
fringing pattern.

**G750M/G750L Point of Difference:** As in the previous step, the input
science product type is dependent on mode. (G750L: ``crj``, G750M:
``sx2``)

.. code:: ipython3

    stistools.defringe.defringe(f"{sci_file}_{prod_type}.fits", f"{flat_file}_frr.fits", overwrite=True)


.. parsed-literal::

    Fringe flat data were read from the primary HDU
    66 pixels in the fringe flat were less than or equal to 0
    Imset 1 done
    Removing and recreating odvkl3050_drj.fits
    Defringed science saved to odvkl3050_drj.fits




.. parsed-literal::

    'odvkl3050_drj.fits'



We now have a ``drj`` file that is the fully defringed calibrated
science spectra. This file functionally behaves as the ``crj`` file used
to produce it, and may be worked with in the same manner.

**G750M/G750L Point of Difference**: If working with a G750M
observation, the output product by default will have the ``s2d``
identifier.
