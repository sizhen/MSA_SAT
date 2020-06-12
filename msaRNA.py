# -*- coding: utf-8 -*-
# @Time    : 6/5/20 15:17
# @Author  : Sizhen Li
# @Email   : annlisizhen@gmail.com
# @FileName: msaRNA.py

from util import multi_seqs_to_vertex_cover, alignment, graph_to_multi_seqs, vertex_cover_solver
import sys

print_info = 0
# 0: don't print any middle info
# 1: only print info in level 1
# 2: print info in level 1 and 2

def multi_seqs_with_two_ones(seqs):
    seqs_pos_ones = []
    for seq in seqs:
        pos_ones = []
        for i, item in enumerate(seq):
            if item == '1':
                pos_ones.append(i)
        seqs_pos_ones.append(pos_ones)

    seq_left_zeros = []
    pruned_seqs = []
    for i in range(len(seqs)):
        if len(seqs_pos_ones[i]) > 0:
            seq_left_zeros.append(len(seqs[i][seqs_pos_ones[i][-1]: ]) - 1)
            pruned_seqs.append(seqs[i][:seqs_pos_ones[i][-1] + 1])
        else:
            seq_left_zeros.append(len(seqs[i]))
            pruned_seqs.append("")

    num_ones = [len(pos_ones) for pos_ones in seqs_pos_ones]
    max_num_ones = max(num_ones)
    num_gene = int((max_num_ones + 1) / 2)

    for i_gene in range(num_gene):
        # decide boundary
        # decide two ones
        new_seqs = []
        for i, seq in enumerate(pruned_seqs):
            new_seq = ''
            pos_ones = seqs_pos_ones[i]
            if (i_gene * 2 + 2) <= len(pos_ones):   # have two more ones
                if i_gene > 0: lnode = pos_ones[i_gene*2 - 1] + 1
                else: lnode = 0
                rnode = seqs_pos_ones[i][i_gene*2 + 1]
                for item in seq[lnode: rnode + 1]:
                    new_seq += item

            elif (i_gene * 2) <= len(pos_ones):
                # only single one or no one left
                if i_gene > 0: lnode = pos_ones[i_gene*2 - 1] + 1
                else: lnode = 0
                for item in seq[lnode: ]:
                    new_seq += item

            new_seqs.append(new_seq)

        yield new_seqs, seq_left_zeros

    yield None, seq_left_zeros

def align(target, seqs):
    # seqs = ["101001", "0010101", "01010"]  # seqs have more than two ones
    new_seqs = []
    for seq in seqs:
        new_seq = ""
        for nuc in seq.upper():
            if nuc == target:
                new_seq += '1'
            else:
                new_seq += '0'
        new_seqs.append(new_seq)

    if print_info >= 1: print("0/1 seqs: ", new_seqs)
    seqs_gen = multi_seqs_with_two_ones(new_seqs)
    template = ''
    for i, (sub_seqs, seq_left_zeros) in enumerate(seqs_gen):
        if sub_seqs == None:
            break

        if print_info >= 2: print("  sub seqs: ", sub_seqs)
        nodes, edges = multi_seqs_to_vertex_cover(sub_seqs)
        num_nodes = len(nodes)
        graph = (nodes, edges)
        if print_info >= 2: print("  graph: ", graph)

        for k in range(num_nodes+1):
            cover, formula, size = vertex_cover_solver(graph, k)
            if cover != None:
                if print_info >= 2: print("  formula: ", formula)
                break

        if print_info >= 2: print("  vertex cover:", cover)
        if print_info >= 2: print("  with minimum k: ", k)
        sub_template = graph_to_multi_seqs(graph, cover)
        if print_info >= 2: print("  sub template: ", sub_template)
        template += sub_template
        if print_info >= 2: print("  " + "--" * 10)

    template += '0' * max(seq_left_zeros)
    new_template, aligned_seqs = alignment(new_seqs, template)

    if print_info >= 1: print("template and aligned seqs:")
    if print_info >= 1: print(new_template)
    if print_info >= 1:
        for i in range(len(seqs)):
            print(aligned_seqs[i])
    if print_info >= 1: print("")
    return new_template, aligned_seqs

def sub_align(target, seqs, template, aligned_seqs):
    # align Gs
    sub_seqs = [""] * len(seqs)
    seqs_pos = [0] * len(seqs)
    new_aligned_seqs = [""] * len(seqs)
    new_template = ""
    for i, item in enumerate(template):
        if item == "1":
            # print(sub_seqs)
            test = sum([len(sub_seq) for sub_seq in sub_seqs])  # all sub_seqs are empty
            if test > 0:
                sub_template, sub_aligned_seqs = align(target, sub_seqs)  # align As
                new_template += sub_template
                for ii in range(len(seqs)):
                    new_aligned_seqs[ii] += sub_aligned_seqs[ii]

            sub_seqs = [""] * len(seqs)
            new_template += "1"
            for j in range(len(seqs)):
                if aligned_seqs[j][i] == "1":
                    seqs_pos[j] += 1
                    new_aligned_seqs[j] += "1"
                else:
                    new_aligned_seqs[j] += "-"
        else:
            for j in range(len(seqs)):
                if aligned_seqs[j][i] == "0":
                    assert seqs[j][seqs_pos[j]] != 'A'
                    sub_seqs[j] += seqs[j][seqs_pos[j]]
                    seqs_pos[j] += 1

    test = sum([len(sub_seq) for sub_seq in sub_seqs])  # all sub_seqs are empty
    if test > 0:
        sub_template, sub_aligned_seqs = align(target, sub_seqs)  # align As
        new_template += sub_template
        for i in range(len(seqs)):
            new_aligned_seqs[i] += sub_aligned_seqs[i]

    if print_info >= 1: print("full template and aligned seqs:")
    if print_info >= 1: print(new_template)
    if print_info >= 1:
        for i in range(len(seqs)):
            print(new_aligned_seqs[i])

    return new_template, new_aligned_seqs

def main(seqs):
    # seqs = ['GAUCAGUC', "AUCAGC", "GACGAUAAAA"]

    if print_info >= 1: print("{}Align A{}".format('-' * N_dash, '-' * N_dash))
    template, aligned_seqs = align("A", seqs)  # align A

    if print_info >= 1: print("{}Align G{}".format('-' * N_dash, '-' * N_dash))
    template, aligned_seqs = sub_align("G", seqs, template, aligned_seqs)  # align G

    if print_info >= 1: print("{}Align U{}".format('-' * N_dash, '-' * N_dash))
    template, aligned_seqs = sub_align("U", seqs, template, aligned_seqs)  # align U

    if print_info >= 1: print("{}FINISH{}".format('-' * N_dash, '-' * N_dash))

    aligned_nucs = []

    for i in range(len(seqs)):
        nuc_j = 0
        nucs = ""
        for j in range(len(aligned_seqs[i])):
            if aligned_seqs[i][j] != "-":
                nucs += seqs[i][nuc_j]
                nuc_j += 1
            else:
                nucs += "-"
        print(nucs)
        aligned_nucs.append(nucs)


if __name__ == '__main__':
    N_dash = 20
    # seqs = ['GAUCAGUC', "AUCAGC", "GACGAUAAAA"]
    seqs = sys.argv[1:]
    main(seqs)
