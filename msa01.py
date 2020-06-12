# -*- coding: utf-8 -*-
# @Time    : 6/12/20 00:26
# @Author  : Sizhen Li
# @Email   : annlisizhen@gmail.com
# @FileName: msa01.py

from util import  alignment, vertex_cover_solver, graph_to_multi_seqs, multi_seqs_to_vertex_cover
import sys

def main(seqs):
    seq_left_zeros = []
    for seq in seqs:
        num_one = 0
        for i, v in enumerate(seq):
            if v == "1":
                num_one += 1
            else:
                if num_one == 2:
                    seq_left_zeros.append(len(seq) - i)
                elif num_one > 2:
                    print("%s has more than two ones!" % seq)

    nodes, edges = multi_seqs_to_vertex_cover(seqs)
    num_nodes = len(nodes)
    graph = (nodes, edges)

    for k in range(num_nodes + 1):
        cover, formula, size = vertex_cover_solver(graph, k)
        if cover != None:
            break

    template = graph_to_multi_seqs(graph, cover)
    template += '0' * max(seq_left_zeros)
    new_template, aligned_seqs = alignment(seqs, template)

    print(new_template)
    for i in range(len(seqs)):
        print(aligned_seqs[i])


if __name__ == '__main__':
    N_dash = 20
    # seqs = ['0100000100', "000010001", "0001000100"]
    seqs = sys.argv[1:]
    main(seqs)
