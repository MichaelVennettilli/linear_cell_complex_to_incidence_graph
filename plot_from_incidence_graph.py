#!/usr/bin/env python
# Created by michael at 4/6/22
# Filename: plot_from_incidence_graph.py
""" Description goes here.
MIT License

Copyright (c) 2022 Michael Vennettilli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import csv
import mpl_toolkits.mplot3d as a3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def generate_adjacency_structures(csv_with_adjacency):
    # Read in the csv, convert to a list and extract the numbers of geometric elements, then remove the first row.
    with open(csv_with_adjacency) as f:
        incidence_reader = csv.reader(f)
        total_adjacency_list = list(incidence_reader)

    num_vertices = int(total_adjacency_list[0][0])
    num_edges = int(total_adjacency_list[0][1])
    num_faces = int(total_adjacency_list[0][2])
    num_volume = int(total_adjacency_list[0][3])
    del total_adjacency_list[0]

    # Create the lists of sets for storing incidence
    volume_face_incidences = []
    face_edge_incidences = []
    face_vertex_incidences = []
    edge_vertex_incidences = []

    for volumes in np.arange(num_volume):
        volume_face_incidences.append([])

    for faces in np.arange(num_faces):
        face_edge_incidences.append([])
        face_vertex_incidences.append(set())

    for edges in np.arange(num_edges):
        edge_vertex_incidences.append([])

    # Populate the adjacency entries for all but the face-vertex incidences
    for row in total_adjacency_list:
        entry_info = row[0].split('_')
        which_table_type = int(entry_info[0])
        which_high_d_entry = int(entry_info[1])
        which_low_d_entry = int(row[1].split('_')[1])
        if which_table_type == 1:
            edge_vertex_incidences[which_high_d_entry].append(which_low_d_entry)
        elif which_table_type == 2:
            face_edge_incidences[which_high_d_entry].append(which_low_d_entry)
        else:
            volume_face_incidences[which_high_d_entry].append(which_low_d_entry)

    # Populate the face-vertex incidence list.
    for face in np.arange(len(face_vertex_incidences)):
        for edge in face_edge_incidences[face]:
            for vertex in edge_vertex_incidences[edge]:
                face_vertex_incidences[face].add(vertex)
        face_vertex_incidences[face] = list(face_vertex_incidences[face])

    return [volume_face_incidences, face_vertex_incidences]


def generate_points(csv_with_points):
    with open(csv_with_points) as f:
        reader = csv.reader(f)
        point_list = list(reader)

    for i in np.arange(len(point_list)):
        point_list[i] = tuple(map(float, point_list[i]))

    return point_list


def plot_triangulated_complex_random_color(csv_with_adjacency, csv_with_points):
    [volume_to_face_list, face_to_vertex_list] = generate_adjacency_structures(csv_with_adjacency)
    vertex_positions = generate_points(csv_with_points)
    fig = plt.figure()
    ax = a3.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.dist = 30
    ax.azim = +40
    ax.elev = 20

    max_coordinates = list(map(max, zip(*vertex_positions)))
    min_coordinates = list(map(min, zip(*vertex_positions)))

    ax.set_xlim([min_coordinates[0], max_coordinates[0]])
    ax.set_ylim([min_coordinates[1], max_coordinates[1]])
    ax.set_zlim([min_coordinates[2], max_coordinates[2]])

    for volume in np.arange(len(volume_to_face_list)):
        color_vec = np.random.rand(3)
        faces_for_one_vol = volume_to_face_list[volume]
        for face in faces_for_one_vol:
            triangle = [vertex_positions[face_to_vertex_list[face][0]], vertex_positions[face_to_vertex_list[face][1]],
                        vertex_positions[face_to_vertex_list[face][2]]]
            drawn_face = a3.art3d.Poly3DCollection([triangle])
            drawn_face.set_color(colors.rgb2hex(color_vec))
            # face.set_edgecolor()
            drawn_face.set_alpha(0.5)
            ax.add_collection3d(drawn_face)

    # ax = plt.gca(projection='3d')
    ax.set_axis_off()
    plt.show()

