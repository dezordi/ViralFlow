from setuptools import setup

setup(name='ViralFlow',
      version='1.0.0',
      description='''
      Workflows for viral genome Assembly at FioCruz/IAM
      ''',
      url='https://github.com/dezordi/ViralFlow/',
      author='Antonio Marinho & Filipe Z. Dezordi',
      author_email='amarinhosn@gmail.com & zimmer.filipe@gmail.com',
      #packages=['viralflow'],
      py_modules = [],
      classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
          ],
      scripts=['wrapper/viralflow_dev'],
      zip_safe=False
)
