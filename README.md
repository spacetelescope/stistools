# stistools

[![Jenkins CI Status](https://ssbjenkins.stsci.edu/job/STScI/job/stistools/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/stistools/job/master/)

Tools for HST/STIS.

Code Contribution Guide:

- For new additions to stistools, a new branch off the main repository is encouraged.  Use initials for the beginning of the branch title. It is also acceptable to put a PR in from a fork of stistools (necessary for an external contributor).

- Each PR should have at least one approved review by at least one STIS team member AND one DATB/SCSB member (this could be either Sara or Robert).

- After approved reviews, test your new content with readthedocs.  You can do this by pushing to the doc_updates_rtd branch.  This branch is setup to re-build the https://stistools.readthedocs.io/en/doc_updates_rtd/ page after any new commits.

## Documentation
Minimum requirement to have at least some inline numpy style API docstrings.  PR author should also review narrative docs to make sure those are appropriately updated. Any new tasks will need a new rst file for sphinx/rtd to pick up the new docs.

## Testing
New functions and or new functionality should have appropriate unit tests.  Tests that use any input and or truth files, will need to use artifactory for hosting test files.

## Pep 8
Try to adhere to pep 8 standards when reasonable.  Code comments are heartily encouraged!
