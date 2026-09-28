# Makefile for Sphinx documentation

SHELL = /bin/bash

PACKAGE_NAME = ribbon
MODULE_PATH = ../src/$(PACKAGE_NAME)
RST_SOURCE_DIR  = $(CURDIR)/source
BUILD_DIR = $(CURDIR)/build

apidoc: export SPHINX_APIDOC_OPTIONS = members,inherited-members,show-inheritance
#available options = members,undoc-members,inherited-members,show-inheritance

.PHONY: doc
doc:
	sphinx-build -b html "$(RST_SOURCE_DIR)" "$(BUILD_DIR)"/html

.PHONY: rst   
rst:
	sphinx-apidoc -M -e -f -d 2 -o "$(RST_SOURCE_DIR)" "$(MODULE_PATH)"

.PHONY: clean
clean :
	rm -rf "$(BUILD_DIR)"/*

.PHONY: clobber
clobber :
	rm -rf "$(RST_SOURCE_DIR)"/$(PACKAGE_NAME)*.rst \
		"$(RST_SOURCE_DIR)"/modules.rst 	\
		"$(BUILD_DIR)"/*
