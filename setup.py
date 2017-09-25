#!/usr/bin/env python

from setuptools import setup

version = '0.1.0'

required = open('requirements.txt').read().split('\n')

setup(
    name='maca',
    version=version,
    description='Command line utilities for RNA-sequencing data',
    author='Olga Botvinnik',
    author_email='olga.botvinnik@gmail.com',
    url='https://github.com/czbiohub/maca',
    packages=['maca'],
    install_requires=required,
    long_description='See ' + 'https://github.com/czbiohub/maca',
    license='MIT',
    entry_points={"console_scripts": ['maca = maca.cli:cli']}
)
