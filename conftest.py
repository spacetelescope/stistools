"""Project default for pytest"""
import os
import pytest
import re
import crds

from astropy.tests.helper import enable_deprecations_as_exceptions


# Uncomment the following line to treat all DeprecationWarnings as exceptions
enable_deprecations_as_exceptions()

def pytest_addoption(parser):
    # Add option to run slow tests
    parser.addoption(
        "--runslow",
        action="store_true",
        help="run slow tests"
    )

    parser.addoption(
        "--slow",
        action="store_true",
        help="run slow tests"
    )

    # Add option to use big data sets
    parser.addoption(
        "--bigdata",
        action="store_true",
        help="use big data sets (intranet)"
    )

    parser.addoption(
        "--env",
        choices=['dev', 'stable', ''],
        default='',
        help="specify what environment to test"
    )


@pytest.fixture(scope='function', autouse=True)
def _jail(tmpdir):
    """ Perform test in a pristine temporary working directory
    """
    os.chdir(tmpdir.strpath)
    yield

@pytest.fixture
def envopt(request):
    return request.config.getoption("env")


def require_crds_context(required_context):
    """Ensure CRDS context is a certain level
    Parameters
    ----------
    level: int
        The minimal level required
    Returns
    -------
    pytest.mark.skipif decorator
    """
    current_context_string = crds.get_context_name('jwst')
    match = re.match('jwst_(\d\d\d\d)\.pmap', current_context_string)
    current_context = int(match.group(1))
    return pytest.mark.skipif(
        current_context < required_context,
        reason='CRDS context {} less than required context {}'.format(
            current_context_string, required_context
        )
    )

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "bigdata: use big data sets (intranet)")
    config.addinivalue_line("markers", "not_under_travis: mark as test to skip if running under a TravisCI")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    skip_bigdata = pytest.mark.skip(reason="need --bigdata option to run")
    skip_travis = pytest.mark.skip(reason="temporarily disabled due to performance issues")
    for item in items:
        if "slow" not in item.keywords:
            item.add_marker(skip_slow)
        if "bigdata" not in item.keywords:
            item.add_marker(skip_bigdata)
        if "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true":
            item.add_marker(skip_travis)
    
