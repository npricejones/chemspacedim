from setuptools import setup
import warnings

setup(name='chemspacedim',
      version='1.',
      description='Thesis work on chemical dimensionality',
      author='Natalie Price-Jones',
      author_email='price-jones@astro.utoronto.ca',
      url='https://github.com/npricejones/chemspacedim',
      package_dir = {'chemspacedim/': ''},
      packages=['chemspacedim/','chemspacedim/tools'],
      package_data={},
      dependency_links = [],
      install_requires=[]
      )
