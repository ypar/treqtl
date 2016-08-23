
from setuptools import setup

setup(name='treqtl',
      version='0.0.5',
      description='scripts to run mendelian randomization analysis on gwas and eqtl summary statistics',
      url='http://github.com/ypar/treqtl',
      author='YoSon Park',
      author_email='ypar@upenn.edu',
      license='MIT',
      packages=['src'],
      zip_safe=False,
      install_requires=[
        'python >= 3.4',
        'pandas >= 1.6',
        'numpy >= 1.10'
        ]
      )


