"""Tests for generalised LALSimulation approximant loading in minke.

Verifies that any LALSimulation approximant can be instantiated by name,
without requiring a dedicated subclass, while existing named subclasses
continue to work unchanged (backwards compatibility).

Run with: python -m pytest tests/test_approximants.py
"""

from __future__ import annotations

import pytest

lalsimulation = pytest.importorskip("lalsimulation")


class TestLALSimulationApproximantByName:
    """Tests for LALSimulationApproximant(approximant_name) factory path."""

    def test_imrphenomxphm_by_name(self):
        """LALSimulationApproximant('IMRPhenomXPHM') creates a usable object."""
        from minke.models.lalsimulation import LALSimulationApproximant
        model = LALSimulationApproximant("IMRPhenomXPHM")
        assert model is not None

    def test_imrphenompv2_by_name(self):
        """LALSimulationApproximant('IMRPhenomPv2') creates a usable object."""
        from minke.models.lalsimulation import LALSimulationApproximant
        model = LALSimulationApproximant("IMRPhenomPv2")
        assert model is not None

    def test_seobrnv2_by_name(self):
        """LALSimulationApproximant('SEOBNRv2') creates a usable object."""
        from minke.models.lalsimulation import LALSimulationApproximant
        model = LALSimulationApproximant("SEOBNRv2")
        assert model is not None

    def test_approximant_stored_correctly(self):
        """The internal approximant integer matches GetApproximantFromString."""
        from minke.models.lalsimulation import LALSimulationApproximant
        model = LALSimulationApproximant("IMRPhenomXPHM")
        expected = lalsimulation.GetApproximantFromString("IMRPhenomXPHM")
        assert model._args["approximant"] == expected

    def test_unknown_approximant_raises(self):
        """An unrecognised approximant name raises a RuntimeError from lalsim."""
        from minke.models.lalsimulation import LALSimulationApproximant
        with pytest.raises(RuntimeError):
            LALSimulationApproximant("NotARealApproximant_xyz")

    def test_get_approximant_factory(self):
        """get_approximant(name) returns an LALSimulationApproximant instance."""
        from minke.models.lalsimulation import get_approximant
        model = get_approximant("IMRPhenomXPHM")
        from minke.models.lalsimulation import LALSimulationApproximant
        assert isinstance(model, LALSimulationApproximant)


class TestNamedSubclassBackwardsCompat:
    """Named subclasses (IMRPhenomXPHM etc.) must still work unchanged."""

    def test_imrphenomxphm_class_still_works(self):
        """IMRPhenomXPHM() constructs without error."""
        from minke.models.lalsimulation import IMRPhenomXPHM
        model = IMRPhenomXPHM()
        assert model is not None

    def test_imrphenompv2_class_still_works(self):
        """IMRPhenomPv2() constructs without error."""
        from minke.models.lalsimulation import IMRPhenomPv2
        model = IMRPhenomPv2()
        assert model is not None

    def test_seobrnv2_class_still_works(self):
        """SEOBNRv2() constructs without error."""
        from minke.models.lalsimulation import SEOBNRv2
        model = SEOBNRv2()
        assert model is not None

    def test_seobrnv3_class_still_works(self):
        """SEOBNRv3() constructs without error."""
        from minke.models.lalsimulation import SEOBNRv3
        model = SEOBNRv3()
        assert model is not None

    def test_named_class_same_approximant_as_by_name(self):
        """IMRPhenomXPHM() and LALSimulationApproximant('IMRPhenomXPHM')
        store the same LALSim approximant integer."""
        from minke.models.lalsimulation import IMRPhenomXPHM, LALSimulationApproximant
        named  = IMRPhenomXPHM()
        by_name = LALSimulationApproximant("IMRPhenomXPHM")
        assert named._args["approximant"] == by_name._args["approximant"]
