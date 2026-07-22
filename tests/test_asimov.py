"""
Tests for the minke.asimov pipeline integration.
"""

import os
import sys
import tempfile
import unittest
import unittest.mock

# Stub out unavailable C extensions before minke.asimov is imported.
for _mod in ("htcondor2", "classad2", "htcondor", "classad"):
    sys.modules[_mod] = unittest.mock.MagicMock()

# Stub asimov so the Asimov class definition works without a full install.
# _FakePipeline must be a real class so that Asimov can inherit from it and
# __new__ / isinstance checks work correctly.
class _FakePipeline:
    pass

_pipeline_mod = unittest.mock.MagicMock()
_pipeline_mod.Pipeline = _FakePipeline

_asimov_mod = unittest.mock.MagicMock()
_asimov_mod.pipeline = _pipeline_mod  # asimov.pipeline.Pipeline resolves here

sys.modules["asimov"] = _asimov_mod
sys.modules["asimov.pipeline"] = _pipeline_mod
sys.modules["asimov.utils"] = unittest.mock.MagicMock()
sys.modules["asimov.config"] = unittest.mock.MagicMock()

from minke.asimov import Asimov  # noqa: E402 — must come after sys.modules setup


def _make_mock_production(rundir):
    production = unittest.mock.MagicMock()
    production.rundir = rundir
    production.event.meta = {"data": {}, "psds": {}}
    return production


class TestCollectAssets(unittest.TestCase):
    """Tests for Asimov.collect_assets()."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.pipeline = Asimov.__new__(Asimov)
        self.pipeline.production = _make_mock_production(self.tmpdir.name)
        self.pipeline.logger = unittest.mock.MagicMock()

    def tearDown(self):
        self.tmpdir.cleanup()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _touch(self, *parts):
        path = os.path.join(self.tmpdir.name, *parts)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").close()
        return path

    def _write_cache(self, ifo, content=""):
        path = os.path.join(self.tmpdir.name, "cache", f"{ifo}.cache")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(content)
        return path

    # ------------------------------------------------------------------
    # frame tests
    # ------------------------------------------------------------------

    def test_frames_collected(self):
        self._touch("H1_Injection.gwf")
        self._touch("L1_Injection.gwf")

        assets = self.pipeline.collect_assets()

        self.assertIn("frames", assets)
        self.assertIn("H1", assets["frames"])
        self.assertIn("L1", assets["frames"])

    def test_frames_advertised_to_event_meta(self):
        self._touch("H1_Injection.gwf")

        self.pipeline.collect_assets()

        self.assertIn("H1", self.pipeline.production.event.meta["data"]["data files"])

    # ------------------------------------------------------------------
    # cache tests
    # ------------------------------------------------------------------

    def test_cache_collected_from_subdirectory(self):
        self._write_cache("H1")
        self._write_cache("L1")

        assets = self.pipeline.collect_assets()

        self.assertIn("cache", assets)
        self.assertIn("H1", assets["cache"])
        self.assertIn("L1", assets["cache"])

    def test_cache_path_points_to_file(self):
        expected = self._write_cache(
            "H1",
            content="H1\tInjection\t1000000000\t4\tfile://localhost/tmp/H1_Injection.gwf\n",
        )

        assets = self.pipeline.collect_assets()

        self.assertEqual(assets["cache"]["H1"], expected)

    def test_cache_advertised_to_event_meta(self):
        self._write_cache("H1")
        self._write_cache("L1")

        self.pipeline.collect_assets()

        meta_cache = self.pipeline.production.event.meta["data"]["cache files"]
        self.assertIn("H1", meta_cache)
        self.assertIn("L1", meta_cache)

    def test_no_cache_key_without_cache_directory(self):
        """collect_assets should not include 'cache' when no cache/ dir exists."""
        assets = self.pipeline.collect_assets()

        self.assertNotIn("cache", assets)

    # ------------------------------------------------------------------
    # PSD tests
    # ------------------------------------------------------------------

    def test_psds_collected(self):
        self._touch("H1_psd.dat")
        self._touch("L1_psd.dat")

        assets = self.pipeline.collect_assets()

        self.assertIn("psds", assets)
        self.assertIn("H1", assets["psds"])
        self.assertIn("L1", assets["psds"])

    def test_psds_advertised_to_event_meta(self):
        self._touch("H1_psd.dat")

        self.pipeline.collect_assets()

        self.assertIn("H1", self.pipeline.production.event.meta["psds"])


if __name__ == "__main__":
    unittest.main()
