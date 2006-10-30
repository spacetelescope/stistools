from distutils.core import setup
import sys, os.path

if not hasattr(sys, 'version_info') or sys.version_info < (2,3,0,'alpha',0):
    raise SystemExit, "Python 2.3 or later required to build stistools."

args = sys.argv[:]

for a in args:
    if a.startswith('--local='):
        dir = os.path.abspath(a.split("=")[1])
        sys.argv.extend([
                "--install-lib="+dir,
                "--install-scripts=%s" % dir])
        args.remove(a)
        sys.argv.remove(a)

setup(name = "stistools",
      version = "1.0.0",
      description = "Python Tools for STIS Data",
      author = "Paul Barrett",
      author_email = "help@stsci.edu",
      license = "http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE",
      platforms = ["Linux","Solaris","Mac OS X","Win"],
      py_modules = ['sshift', 'stisnoise', 'gettable', 'r_util', 'wavelen', 'evaldisp', 'radialvel', 'mktrace', 'wx2d'],
      package_dir={'':'lib'},
      )



