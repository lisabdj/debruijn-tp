#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import statistics
import random
from random import randint
from operator import itemgetter
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
random.seed(9001)


__author__ = "BOUARROUDJ Lisa"
__copyright__ = "Boudjii Corporation"
__credits__ = ["BOUARROUDJ Lisa"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "BOUARROUDJ Lisa"
__email__ = "lisa.bdj.95@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """
      :Parameters:
          xxx:
      Returns:

    """
    with open(fastq_file, "r") as filin:
        ligne=filin.readline()
        while ligne != "":
            while ligne[0]== "@":
                ligne = filin.readline()
                yield ligne[:-1]
                break
            ligne = filin.readline()


def cut_kmer(read, kmer_size):
    """
      :Parameters:
          xxx:
      Returns:

    """
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]



def build_kmer_dict(fastq_file, kmer_size):
    """
      :Parameters:
          xxx:
      Returns:

    """
    sequences=read_fastq(fastq_file)
    list_kmer = []
    for seq in sequences:
        kmers = cut_kmer(seq, kmer_size)
        for kmer in kmers:
            list_kmer.append(kmer)
    dico_kmer = dict((x, list_kmer.count(x)) for x in set(list_kmer))
    return dico_kmer



def build_graph(kmer_dict):
    """
      :Parameters:
          xxx:
      Returns:

    """
    graph = nx.DiGraph()
    for kmer,weight in kmer_dict.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=weight)
    return graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
      :Parameters:
          xxx:
      Returns:

    """
    for path in path_list:
        if delete_entry_node is True and delete_sink_node is True:
            graph.remove_nodes_from(path)
        elif delete_entry_node is True:
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node is True:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph


def std(data):
    """
      :Parameters:
          xxx:
      Returns:

    """
    std = statistics.stdev(data)
    return std


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """
      :Parameters:
          xxx:
      Returns:

    """
    if std(weight_avg_list) > 0:
        max_poids = max(weight_avg_list)
        i_max = weight_avg_list.index(max_poids)
    elif std(path_length) > 0:
        max_len = max(path_length)
        i_max = path_length.index(max_len)
    elif std(path_length) == 0:
        i_max = randint(0, len(path_length))

    path_list.pop(i_max)
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    """for i in range(len(path_list)):
        for j in range(i, len(path_list)-1):
            poids_1 = weight_avg_list[i]
            poids_2 = weight_avg_list[j]
            if std([poids_1, poids_2]) > 0:
                if poids_1 > poids_2:
                    graph.remove_paths(graph, path_list[j], delete_entry_node, delete_sink_node)
                else:
                    graph.remove_paths(graph, path_list[i], delete_entry_node, delete_sink_node)
            if std([poids_1, poids_2]) == 0:
                len_1 = path_length[i]
                len_2 = path_length[j]
                if std([len_1, len_2]) > 0:
                    if len_1 > len_2:
                        graph.remove_paths(graph, path_list[j], delete_entry_node, delete_sink_node)
                    else:
                        graph.remove_paths(graph, path_list[i], delete_entry_node, delete_sink_node)
                elif std([len_1, len_2]) == 0:
                    n = randint(i, j)
                    graph.remove_paths(graph, path_list[n], delete_entry_node, delete_sink_node)"""
    return graph

def path_average_weight(graph, path):
    """
      :Parameters:
          xxx:
      Returns:

    """
    poids = 0
    for noeud_1, noeud_2 in zip(path[:-1], path[1:]):
        poids = poids + graph[noeud_1][noeud_2]["weight"]
    poids_moyen = poids / (len(path)-1)
    return poids_moyen

def solve_bubble(graph, ancestor_node, descendant_node):
    """
      :Parameters:
          xxx:
      Returns:

    """
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = []
    weight_avg_list = []
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
    return select_best_path(graph, path_list, path_length, weight_avg_list)

def simplify_bubbles(graph):
    """
      :Parameters:
          xxx:
      Returns:

    """
    bubble = False
    for noeud in graph.nodes():
        if noeud in graph.nodes():
            liste_prede = list(graph.predecessors(noeud))
            if len(liste_prede) > 1:
                for i in range(len(liste_prede)):
                    for j in range(i+1, len(liste_prede)):
                        noeud_ancetre = nx.lowest_common_ancestor(graph, liste_prede[i], liste_prede[j])
                        if noeud_ancetre is not None:
                            bubble = True
                            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, noeud_ancetre, noeud))
    return graph


def solve_entry_tips(graph, starting_nodes):
    """
      :Parameters:
          xxx:
      Returns:

    """
    return graph

def solve_out_tips(graph, ending_nodes):
    """
      :Parameters:
          xxx:
      Returns:

    """
    return graph

def get_starting_nodes(graph):
    """
      :Parameters:
          xxx:
      Returns:

    """
    noeuds_entree = []
    for noeud in graph.nodes():
        if not list(graph.predecessors(noeud)):
            noeuds_entree.append(noeud)
    return noeuds_entree

def get_sink_nodes(graph):
    """
      :Parameters:
          xxx:
      Returns:

    """
    noeuds_sortie = []
    for noeud in graph.nodes():
        if not list(graph.successors(noeud)):
            noeuds_sortie.append(noeud)
    return noeuds_sortie

def get_contigs(graph, starting_nodes, ending_nodes):
    """
      :Parameters:
          xxx:
      Returns:

    """
    contigs = []
    for n_entree in starting_nodes:
        for n_sortie in ending_nodes:
            if nx.has_path(graph, n_entree, n_sortie) is True:
                paths=nx.all_simple_paths(graph, n_entree, n_sortie)
                for path in paths:
                    contig = "".join([noeud[0] for noeud in path[:-1]]+ [path[-1]])
                    contigs.append((contig, len(contig)))
    return contigs


def save_contigs(contigs_list, output_file):
    """
      :Parameters:
          xxx:
      Returns:

    """
    with open(output_file, 'w+') as filout:
        for i, contig in enumerate(contigs_list):
            filout.write('>contig_{} len={}\n'.format(i, contig[1]))
            filout.write(fill(contig[0]) + '\n')
    filout.close()


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    fastq_file = "/home/sdv/m2bi/lbouarroudj/Documents/Assemblage/debruijn-tp/data/eva71_two_reads.fq"
    #sequence=read_fastq(fastq_file)
    #for seq in sequence:
    #    print(seq)
    dico = build_kmer_dict(fastq_file, 2)
    graph = build_graph(dico)
    draw_graph(graph, "graph.png")
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
