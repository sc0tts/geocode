from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='geogrid',
    version='0.1.0',
    description='classes for provideing GeoGrid and GeoPoint',
    packages=find_packages(exclude='contrib', 'docs', 'tests']),
    python_requires='>=3.6',
)

