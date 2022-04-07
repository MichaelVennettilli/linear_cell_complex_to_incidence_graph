#!/usr/bin/env python
# Created by michael at 4/6/22
# Filename: test_plot.py
""" Description goes here.
"""

from plot_from_incidence_graph import *


incidence_csv = 'test_inc.csv'
point_csv = 'test_points.csv'
plot_triangulated_complex_random_color(incidence_csv, point_csv)
