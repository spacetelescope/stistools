#!/usr/bin/env python
import os
import pkgutil
import sys
from glob import glob
from setuptools import setup, find_packages, Extension
from subprocess import check_call, CalledProcessError


if not pkgutil.find_loader('relic'):
    relic_local = os.path.exists('relic')
    relic_submodule = (relic_local and
                       os.path.exists('.gitmodules') and
                       not os.listdir('relic'))
    try:
        if relic_submodule:
            check_call(['git', 'submodule', 'update', '--init', '--recursive'])
        elif not relic_local:
            check_call(['git', 'clone', 'https://github.com/spacetelescope/relic.git'])

        sys.path.insert(1, 'relic')
    except CalledProcessError as e:
        print(e)
        exit(1)

import relic.release

version = relic.release.get_info()
relic.release.write_template(version, 'stistools')

setup(
    name = 'stistools',
    version = version.pep386,
    author = 'Paul Barrett, Phil Hodge',
    author_email = 'help@stsci.edu',
    description = 'Tools for STIS (Space Telescope Imaging Spectrograph)',
    url = 'https://github.com/spacetelescope/stistools',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires = [
        'astropy>=4.1',
        'numpy',
        'scipy',
        'stsci.tools',
        'pysiaf',
        'astroquery',
    ],
    tests_require = [
        'pytest',
        'six',
        #'requests',
        #'crds',
    ],
    packages = find_packages(),
    package_data = {
        'stistools': [
            'pars/*',
            'LICENSE.txt',
        ]
    },
    entry_points = {
        'console_scripts': [
            'add_stis_s_region = stistools.add_stis_s_region:call_main',
        ],
    },
)
