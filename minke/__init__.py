# -*- coding: utf-8 -*-

__author__ = 'Daniel Williams'
__email__ = 'daniel.williams@ligo.org'

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
