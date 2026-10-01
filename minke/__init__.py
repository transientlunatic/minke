# -*- coding: utf-8 -*-

__author__ = 'Daniel Williams'
__email__ = 'daniel.williams@ligo.org'

from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    pass
