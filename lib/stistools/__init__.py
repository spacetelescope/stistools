from __future__ import division         # confidence high
__version__ = '1.0.0'

try:
    from stistools.svninfo import(__svn_version__, __full_svn_info__,
                                  __setup_datetime__)
except ImportError:
    pass

import calstis
import basic2d
import mktrace
import sshift
import stisnoise
import wx2d

# These lines allow TEAL to print out the names of TEAL-enabled tasks
# upon importing this package.
import os
from stsci.tools import teal
teal.print_tasknames(__name__, os.path.dirname(__file__))
