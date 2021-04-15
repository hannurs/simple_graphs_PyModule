#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from distutils.core import setup, Extension
setup( name = "simple_graphs", \
ext_modules = [Extension( "simple_graphs", \
sources = ["graphsmodule.cpp"] )] )