# Feedback loop for the builder (ticket #9). The fixture build needs MMseqs2 18.8cc5c on Linux amd64 (or macOS arm64).
.PHONY: test test-docker test-invariants test-package lint-workflows

test:
	python -m pytest -p no:cacheprovider -q builder/tests

test-docker:
	docker build -q -f builder/tests/Dockerfile -t chamredb-builder-test builder
	docker run --rm -e PYTHONDONTWRITEBYTECODE=1 -v "$(CURDIR)":/work chamredb-builder-test make test

# Graph invariants on the shipped graph (no MMseqs2 needed; included in `test` too).
test-invariants:
	python -m pytest -p no:cacheprovider -q builder/tests/test_graph_invariants.py

# Package unit tests and CLI smoke test (needs `pip install -e .[test]`; CI tests the installed wheel).
test-package:
	python -m pytest -p no:cacheprovider -q chamredb/tests

# Lint the GitHub Actions workflows (CI runs the same check).
lint-workflows:
	uvx --from actionlint-py actionlint .github/workflows/*.yml
