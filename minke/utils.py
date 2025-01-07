"""
Various utilities and non-GPR based stuff.
"""

import yaml


def load_yaml(filename):
    with open(filename, "r") as f:
        data = yaml.safe_load(f)

    return data
