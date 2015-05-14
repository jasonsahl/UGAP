#!/usr/bin/env python
from distutils.core import setup

__author__ = "Jason Sahl"
__credits__ = ["Jason Sahl"]
__license__ = "GPL v3"
__version__ = "1.0"
__maintainer__ = "Jason Sahl"
__email__ = "jsahl@tgen.org"
__status__ = "Development"
 
long_description = """Unnamed Genome
Assembly Pipeline
"""

setup(name='UGAP',
      version=__version__,
      description='UGAP',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      packages=['ugap'],
      long_description=long_description
)
