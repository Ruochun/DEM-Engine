"""Execute the workflow's actual validator against small wheel metadata fixtures."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest
from zipfile import ZipFile


class WheelWorkflowTests(unittest.TestCase):
    """Catch missing variables and accidental relaxation of wheel validation."""

    @classmethod
    def setUpClass(cls):
        # Extract the heredoc itself so tests exercise the code GitHub executes.
        default = Path(__file__).resolve().parents[2] / ".github/workflows/python-wheels.yml"
        workflow = Path(os.environ.get("DEME_TEST_WHEEL_WORKFLOW", default)).read_text()
        step = workflow.split("- name: Validate wheel metadata and native platform tag", 1)[1]
        block = step.split("python - <<'PY'\n", 1)[1].split("\n          PY\n", 1)[0]
        cls.validator = textwrap.dedent(block)

    def validate(self, branch, abi="cp313", metadata_name=None, extra=True, member=None, filename_name=None):
        """Render the two Actions expressions and execute against a temporary wheel."""
        distribution = "deme3" if branch == "Mesh_Particles_Py" else "deme"
        source = self.validator.replace("${{ matrix.python }}", abi).replace(
            "${{ github.ref_name == 'Mesh_Particles_Py' && 'deme3' || 'deme' }}", distribution)
        self.assertNotIn("${{", source, "Update the test renderer for new workflow expressions")
        with tempfile.TemporaryDirectory() as directory:
            wheelhouse = Path(directory) / "wheelhouse"
            wheelhouse.mkdir()
            filename = f"{filename_name or distribution}-3.0.12-{abi}-{abi}-manylinux_2_28_x86_64.whl"
            with ZipFile(wheelhouse / filename, "w") as archive:
                metadata = f"Metadata-Version: 2.1\nName: {metadata_name or distribution}\nVersion: 3.0.12\n"
                if extra:
                    metadata += "Provides-Extra: cuda12\n"
                archive.writestr(f"{distribution}-3.0.12.dist-info/METADATA", metadata)
                if member:
                    archive.writestr(member, b"fixture")
            return subprocess.run([sys.executable, "-c", source], cwd=directory, capture_output=True, text=True)

    def test_valid_wheels_for_all_branches_and_python_versions(self):
        for branch in ("Mesh_Particles_Py", "Mesh_Particles", "main"):
            for abi in ("cp39", "cp310", "cp311", "cp312", "cp313", "cp314"):
                with self.subTest(branch=branch, abi=abi):
                    result = self.validate(branch, abi)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("Validated ", result.stdout)

    def test_rejects_invalid_metadata_and_bundled_libraries(self):
        for branch in ("Mesh_Particles_Py", "Mesh_Particles"):
            for invalid in (
                {"metadata_name": "wrong-package"},
                {"filename_name": "wrong-package"},
                {"extra": False},
                {"member": "deme.libs/libnvrtc.so.12"},
                {"member": "lib/libDEMEVisualizer.so"},
            ):
                with self.subTest(branch=branch, invalid=invalid):
                    result = self.validate(branch, **invalid)
                    self.assertEqual(result.returncode, 1, result.stderr)
                    self.assertNotIn("Traceback", result.stderr)
                    self.assertTrue(result.stderr.strip())


if __name__ == "__main__":
    unittest.main()
