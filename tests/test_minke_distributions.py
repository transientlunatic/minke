import unittest
from minke import distribution
import numpy as np

class TestMinkeDistributions(unittest.TestCase):
    def setUp(self):
        pass
        
    def test_uniform_interval(self):
        # The uniform interval should return values which are between 0 
        # and 10
        distro = distribution.uniform_interval([0,10], 1000)
        for sample in distro:
            self.assertGreater(sample, 0)
            self.assertLess(sample, 10)

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
