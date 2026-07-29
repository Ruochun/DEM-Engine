SPHINXBUILD ?= sphinx-build
DOXYGEN ?= doxygen
SOURCEDIR = .
BUILDDIR = _build

.PHONY: help doxygen html clean

help:
	@echo "make html     Build Doxygen XML and the Sphinx HTML site"
	@echo "make doxygen  Build only the Doxygen XML"
	@echo "make clean    Remove generated documentation"

doxygen:
	@mkdir -p $(BUILDDIR)/doxygen
	$(DOXYGEN) Doxyfile

html: doxygen
	$(SPHINXBUILD) -W --keep-going -b html $(SOURCEDIR) $(BUILDDIR)/html
	@echo "Documentation: $(CURDIR)/$(BUILDDIR)/html/index.html"

clean:
	$(SPHINXBUILD) -M clean $(SOURCEDIR) $(BUILDDIR)
