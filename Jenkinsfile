if (utils.scm_checkout()) return


data_config = new DataConfig()
data_config.server_id = "bytesalad"
data_config.root = 'tests_output'
data_config.match_prefix = '(.*)_result' // .json is appended automatically


bc0 = new BuildConfig()
bc0.nodetype = "linux"
bc0.name = "Install"
bc0.conda_channels = ['https://ssb.stsci.edu/astroconda']
bc0.conda_packages = ['python=3.6',
                     'astropy']
bc0.test_configs = [data_config]
bc0.build_cmds = ["conda env update --file=doc/environment.yml -q",
                  "with_env -n stistools python setup.py install",
                  "with_env -n stistools conda install -q -y -c http://ssb.stsci.edu/astroconda hstcal",
                  "with_env -n stistools conda install -q -y -c http://ssb.stsci.edu/astroconda crds"]
bc0.test_cmds = ["with_env -n stistools conda install -q -y pytest=4.2.1",
                 "with_env -n stistools conda install -q -y pytest-remotedata",
                 "with_env -n stistools pytest --basetemp=tests_output --junit-xml=results.xml --bigdata"]

utils.run([bc0])
