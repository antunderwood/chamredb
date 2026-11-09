# Feedback loop for the builder (ticket #9). The fixture build needs MMseqs2 13.45111 on Linux amd64.
.PHONY: test test-docker test-invariants

test:
	python -m pytest -p no:cacheprovider -q builder/tests

test-docker:
	docker build -q -f builder/tests/Dockerfile -t chamredb-builder-test builder
	docker run --rm -e PYTHONDONTWRITEBYTECODE=1 -v "$(CURDIR)":/work chamredb-builder-test make test

# Graph invariants on the shipped graph (no MMseqs2 needed; included in `test` too).
test-invariants:
	python -m pytest -p no:cacheprovider -q builder/tests/test_graph_invariants.py
