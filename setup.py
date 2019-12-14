from setuptools import setup 
  
# reading long description from file 
with open('DESCRIPTION.txt') as file: 
    long_description = file.read() 
  
  
# specify requirements of your package here 
REQUIREMENTS = ['numpy', 'matplotlib'] 
  
# some more details 
CLASSIFIERS = [ 
    'Development Status :: 4 - Beta', 
    'Intended Audience :: Science/Research',
    'Topic :: File-import', 
    #'License :: OSI Approved :: MIT License', 
    'Programming Language :: Python', 
    'Programming Language :: Python :: 2', 
    'Programming Language :: Python :: 2.7', 
    'Programming Language :: Python :: 3', 
    ] 
  
# calling the setup function  
setup(name='DOPpy', 
      version='2.12-beta', 
      description='Import binary DOP data (.bdd) of the DOP measurement systems', 
      long_description=long_description, 
      url='https://github.com/pyZerrenner/DOPpy', 
      author='pyZerrenner', 
      #author_email='', 
      #license='MIT', 
      packages=['geo'], 
      classifiers=CLASSIFIERS, 
      install_requires=REQUIREMENTS, 
      keywords='file-import'
      ) 
