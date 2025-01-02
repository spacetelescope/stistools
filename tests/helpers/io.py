import copy
import json
import os
import shutil

UPLOAD_SCHEMA = {"files": [
                    {"pattern": "",
                     "target": "",
                     "props": None,
                     "recursive": "false",
                     "flat": "true",
                     "regexp": "false",
                     "explode": "false",
                     "excludePatterns": []
                    }
                  ]
                }

__all__ = ['upload_results',]


def upload_results(**kwargs):
    """Write out JSON file to upload results from test to storage area.
    This function relies on the JFROG JSON schema for uploading data into
    artifactory using the Jenkins plugin.  Docs can be found at::
        https://www.jfrog.com/confluence/display/RTF/Using+File+Specs
    Parameters
    ----------
    pattern : str or list of strings
        Specifies the local file system path to test results which should be
        uploaded to Artifactory. You can specify multiple artifacts by using
        wildcards or a regular expression as designated by the regexp property.
    target : str
        Specifies the target path in Artifactory in the following format:
            [repository_name]/[repository_path]
    testname : str
        Name of test that generate the results.  This will be used to create the
        name of the JSON file to enable these results to be uploaded to Artifactory.
    recursive : bool, optional
        Specify whether or not to identify files listed in sub-directories
        for uploading.  Default: False
    """
    # Interpret mandatory inputs
    pattern = kwargs.get("pattern")
    target = kwargs.get("target")
    testname = kwargs.get("testname")

    # Finish interpreting inputs
    jsonfile = "{}_results.json".format(testname)
    recursive = repr(kwargs.get("recursive", False)).lower()

    if isinstance(pattern, list):
        # Populate schema for this test's data
        upload_schema = {"files": []}

        for p in pattern:
            temp_schema = copy.deepcopy(UPLOAD_SCHEMA["files"][0])
            temp_schema.update({"pattern": p, "target": target, "recursive": recursive})
            upload_schema["files"].append(temp_schema)

    else:
        # Populate schema for this test's data
        upload_schema = copy.deepcopy(UPLOAD_SCHEMA)
        upload_schema["files"][0].update({"pattern": pattern, "target": target, "recursive": recursive})

    # Write out JSON file with description of test results
    with open(jsonfile, 'w') as outfile:
        json.dump(upload_schema, outfile)
