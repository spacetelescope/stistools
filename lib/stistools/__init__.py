try:
    from stistools.svninfo import(__version__, __svn_revision__,
                                  __svn_full_info__, __setup_datetime__)
except ImportError:
    __version__ = '1.0.0'
    __svn_revision__ = ''
    __svn_full_info__ = ''
    __setup_datetime__ = None

import calstis
import basic2d
import ocrreject
import x1d
import mktrace
import sshift
import stisnoise
import wx2d

# These lines allow TEAL to print out the names of TEAL-enabled tasks
# upon importing this package.
import os
from stsci.tools import teal
teal.print_tasknames(__name__, os.path.dirname(__file__))
