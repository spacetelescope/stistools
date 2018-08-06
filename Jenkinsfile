if (utils.scm_checkout()) return

bc0 = new BuildConfig()
bc0.nodetype = "linux"
bc0.name = "Install"
bc0.env_vars = ['PREFIX=/tmp/testbed',
                'PATH=/tmp/testbed/bin:$PATH']
bc0.build_cmds = ["make",
                  "make install PREFIX=$PREFIX"]
bc0.test_cmds = ["make check",
                 "my_program --help",
                 "my_program --version"]

utils.run([bc0])