import os
from os.path import abspath, join, dirname
from shutil import copyfile

curdir = abspath(dirname(__file__))
rootdir = dirname(curdir)
apidoc_cmd = 'sphinx-apidoc'
apidoc_exclusions = ''
apidoc_params = '-f -o {} {} {}'.format(curdir, rootdir, apidoc_exclusions)


def main():
    copyfile(join(rootdir, 'README.rst'), join(curdir, 'README.rst'))
    os.system('{} {}'.format(apidoc_cmd, apidoc_params))


if __name__ == '__main__':
    main()
