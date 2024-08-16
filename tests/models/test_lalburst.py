import unittest
from minke.models.lalburst import SineGaussian

class TestSineGaussian(unittest.TestCase):

    def setUp(self):
        self.model = SineGaussian()

    def test_time_domain_simple(self):

        parameters = {"centre_frequency": 20,
                      "phase": 0,
                      "eccentricity": 0,
                      "q": 1.,
                      "hrss": 1e-22,
                      "duration": 2}
        
        data = self.model.time_domain(parameters)
        f = data['plus'].plot()
        f.savefig("sinegaussian.png")
        self.assertEqual(len(data.keys()), 2)

    
