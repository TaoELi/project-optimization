.PHONY: doc html clean

SPHINXBUILD ?= sphinx-build
DOCS_SOURCE = source
DOCS_BUILD = build/html

doc:
	LC_ALL=C LANG=C $(SPHINXBUILD) -b html $(DOCS_SOURCE) $(DOCS_BUILD)

html:
	python -c 'import pathlib, webbrowser, sys; index = pathlib.Path("build/html/index.html").resolve(); sys.exit("Docs not built; run `make doc` first.") if not index.exists() else webbrowser.open(index.as_uri())'

clean:
	rm -rf build
