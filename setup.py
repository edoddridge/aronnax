from setuptools import setup

setup(name='MIMutils',
      version='0.1',
      description='An idealized isopycnal model with n layers and variable bathymetry.',
      url='http://github.com/edoddridge/MIM',
      author='Ed Doddridge',
      author_email='blank',
      license='MIT licence',
      packages=['MIMutils'],
      package_dir={"MIMutils": "src"},
      install_requires=[
          'numpy',],
      zip_safe=False)
