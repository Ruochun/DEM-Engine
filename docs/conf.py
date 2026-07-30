"""Sphinx configuration for the DEM-Engine documentation."""

from pathlib import Path


DOCS_DIR = Path(__file__).resolve().parent

project = "DEM-Engine"
author = "DEM-Engine contributors"
copyright = "2021–2026, DEM-Engine contributors"
version = "3"
release = "3"

extensions = [
    "breathe",
    "sphinx.ext.autosectionlabel",
]

breathe_projects = {
    "DEME": str(DOCS_DIR / "_build" / "doxygen" / "xml"),
}
breathe_default_project = "DEME"
breathe_default_members = ("members", "undoc-members")

autosectionlabel_prefix_document = True
# Breathe sees implementation types that are deliberately outside the public
# Doxygen input set. Keep authored directive and document warnings fatal while
# allowing those signature types to render as plain identifiers.
nitpicky = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "furo"
html_title = f"DEM-Engine {release}"
