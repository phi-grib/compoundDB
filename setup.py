from setuptools import setup

setup(name='phidb',
      version='0.1',
      description='Tools for creating, querying and managing the compounds and the annotations PostgreSQL database',
      url='http://github.com/phi-grib/compoundDB',
      author='Elisabet Gregori',
      author_email='elisabet.gregori@upf.edu',
      license='GNU',
      packages=['compoundDB', 'annotationDB'],
      package_data={
        'annotationDB': ['data/*'],
        'compoundDB': ['data/*'],
        },
      zip_safe=False)
