from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from pangolin import __version__, _program

setup(name='quokka',
      version=__version__,
      packages=find_packages(),
      scripts=[],
      package_data={"pangolin":["data/*"]},
      install_requires=[
            "biopython>=1.70",
            'pandas>=1.0.1',
            "wheel>=0.34",
            'joblib>=0.11',
            'pysam>=0.16.0',
            'scikit-learn>=0.23.1',
            "PuLP>=2"
        ],
      description='quality of known klusters assessment',
      url='https://github.com/cov-lineages/quokka',
      author='Pango Network',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = pangolin.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)