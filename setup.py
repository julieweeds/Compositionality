__author__ = 'juliewe'
from setuptools import setup
#from disttest import test

# Import this to prevent spurious error: info('process shutting down')
from multiprocessing import util

setup(name='compositionality',
      version='0.1',
      description='Experiments on compositionality of compounds',
      #url='http://github.com/storborg/funniest',
      author='Julie Weeds',
      #author_email='flyingcircus@example.com',
      license='MIT',
      packages=['compositionality'],
      #cmdclass = {'test': test},
      #options = {'test' : {'test_dir':['test']}}
      #zip_safe=False
      )
