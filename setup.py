from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='georef',
    version='0.1.0',
    author='J Scott Stewart',
    author_email='james.stewart@colorado.edu',
    description='classes for provideing GeoGrid and GeoPoint',
    long_description=long_description,
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires=('scipy', 'gdal', 'pytest'),
    python_requires='>=3.6',
)

