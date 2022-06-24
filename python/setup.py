#!/usr/bin/env python3

from setuptools import setup, find_packages
import numpy as np

with open('../README.md') as f:
    long_description = f.read()

setup(
    name='IGA_quadrature',
    version='1.0.1',
    description='Generation of IGA quadrature',
    long_description_content_type='text/markdown',
    long_description=long_description,
    keywords=['Bspline', 'Splines', 'NURBS', 'Quadrature', 'Integration', 'Differentiation'],
    url='https://github.com/vikingscientist/IGA-quadrature',
    maintainer='Kjetil Andre Johannessen',
    maintainer_email='kjetijo@gmail.com',
    license='GNU public license v3',
    packages=find_packages(),
    package_data={
        'iga_quadrature': ['templates/*.bpt'],
    },
    install_requires=[
        'splipy         >= 1.5',
        'numpy          >= 1.20',
        'scipy          >= 1.2',
        'scikit-umfpack >= 0.3'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
    ],
)
