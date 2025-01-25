import numpy as np
from lal import antenna, MSUN_SI
from astropy import units as u


class WaveformModel:

    def _convert(self, args):
        if "mass_ratio" in args and "chirp_mass" in args:
            args["total_mass"] = (
                args["chirp_mass"]
                * (1 + args["mass_ratio"]) ** 1.2
                / args["mass_ratio"] ** 0.6
            )
        if "total_mass" in args and "mass_ratio" in args:
            args = self._convert_mass_ratio_total_mass(args)
        if "luminosity_distance" in args:
            args = self._convert_luminosity_distance(args)
        if "geocent_time" in args:
            args["gpstime"] = args.pop("geocent_time")

        return args

    def _convert_luminosity_distance(self, args):
        args["distance"] = args.pop("luminosity_distance")
        return args

    def _convert_mass_ratio_total_mass(self, args):        
        """
        Convert a mass ratio and a total mass into individual component masses.
        If the masses have no units they are assumed to be in SI units.

        Parameters
        ----------
        args['total_mass'] : float, `astropy.units.Quantity`
           The total mass for the system.
        args['mass_ratio'] : float
           The mass ratio of the system using the convention m2/m1
        """
        args["m1"] = (args["total_mass"] / (1 + args["mass_ratio"]))
        args["m2"] = (args["total_mass"] / (1 + (1 / args["mass_ratio"])))
        # Do these have units?
        # If not then we can skip some relatively expensive operations and apply a heuristic.
        if isinstance(args["m1"], u.Quantity):
            args["m1"] = args["m1"].to_value(u.kilogram)
            args["m2"] = args["m2"].to_value(u.kilogram)
        if (not isinstance(args["m1"], u.Quantity)) and (args["m1"] < 1000):
            # This appears to be in solar masses
            args["m1"] *= MSUN_SI
        if (not isinstance(args["m2"], u.Quantity)) and (args["m2"] < 1000):
            # This appears to be in solar masses
            args["m2"] *= MSUN_SI
        
        args.pop("total_mass")
        args.pop("mass_ratio")

        return args


class WaveformApproximant(WaveformModel):
    """
    This class handles a waveform approximant model.
    """

    pass


class WaveformSurrogate(WaveformModel):
    """
    This class handles a waveform surrogate model.
    """

    pass


class PSDModel:

    def to_file(self, filename, *args, **kwargs):
        data = self.twocolumn(*args, **kwargs)
        np.savetxt(filename, data, fmt="%.2f %.5e")

class PSDApproximant(PSDModel):
    pass
