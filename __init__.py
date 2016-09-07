"""
A module to produce and analyze output files from gPhoton v1.27.x - see https://github.com/cmillion/gPhoton

Alexander de la Vega, Luciana Bianchi & Halley Cromley (JHU)

Includes scripts for stand-alone file creation and overall analysis and 
functions for individual analysis.

last modified: 5 September 2016
"""

__version__ = '0.1'

import make_files

try:
	import astropy
except ImportError:
	print 'gPhoton Analysis requires Astropy'

try:
	import gPhoton
except ImportError:
	print 'gPhoton Analysis requires gPhoton'

try:
	import pyastronomy
except ImportError:
	print 'gPhoton Analysis requires PyAstronomy'