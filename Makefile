.PHONY : install test check build docs clean push_release

test:
	nosetests --verbose --with-coverage --cover-package fqc

check:
	flake8 fqc && echo OK
	yapf -r --diff fqc && echo OK

build:
	python setup.py sdist bdist_wheel

docs:
	sphinx-build -a docs docs/_build

clean:
	rm -rf build
	rm -rf dist
	rm -rf fqc.egg-info
	rm -rf docs/_build
	rm -rf docs/api

bump_patch:
	bumpversion patch

bump_minor:
	bumpversion minor

bump_major:
	bumpversion major

push_release:
	git push && git push --tags
