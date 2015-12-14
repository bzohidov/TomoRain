#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from distutils.core import setup

BIN_FOLDER = 'bin'


def readme():
    with open('README.md') as f:
        return f.read()


def apply_folder_join(item):
    return os.path.join(BIN_FOLDER, item)


setup(
    name = "TomoRain",
    version = "0.1.0",
    description = 'Software package for rainfall monitoring based on cellular networks',
    url = 'https://github.com/bzohidov/TomoRain',
    author = "Bahtiyor Zohidov",
    author_email = "bahtiyor.zohidov@eleves.ec-nantes.fr",
    packages = ["TomoRain"],
    scripts = ["bin/signal_generator.py","bin/rain_retrieval.py"],
    license = "LICENSE.txt",
    long_description = open('README.txt').read(),
    install_requires = [
        "Scipy >= 1.1.1",
        "Numpy >= 0.1.4",
        "Pandas >= 0.11.0",
        "pyproj >= 1.9.4",
        "xlrd >= 0.9.2"
        "PIL >= 1.1.7"],
    )
