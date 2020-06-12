# -*- coding: utf-8 -*-
# @Time    : 6/12/20 00:35
# @Author  : Sizhen Li
# @Email   : annlisizhen@gmail.com
# @FileName: util.py

from pysmt.shortcuts import Symbol, And, Not, is_sat, Or
from pysmt.shortcuts import GE, LE, LT, Plus, Equals, Int, get_model, Solver, get_formula_size
from pysmt.typing import INT
import time



def vertex_cover_solver(graph, k):
    vertexes, edges = graph

    # vertex in vertex cover 1 otherwise 0
    vertex_vars = [Symbol(v, INT) for v in vertexes]
    vertex_range = And([And(GE(v, Int(0)), LE(v, Int(1))) for v in vertex_vars])

    # size of vertex cover <= k
    sum_vertex = Plus(vertex_vars)
    vertex_constraint = LE(sum_vertex, Int(k))

    # edge constraints
    # Plus(vi, vj) >= 1 for edge = (vi, vj)
    # cannot use Or(vi, vj), because variables are integers not boolean type
    edge_constraint = And([GE(Plus([vertex_vars[lv], vertex_vars[rv]]), Int(1)) for (lv, rv) in edges])

    # combined formula
    formula = And(vertex_range, vertex_constraint, edge_constraint)

    model = get_model(formula)
    size = get_formula_size(formula)
    cover = set()
    if model:
        for i, node in enumerate(vertexes):
            if model.get_py_value(vertex_vars[i]):
                cover.add(node)
        return cover, formula, size
    else:
        return None, formula, size

def graph_to_multi_seqs(graph, cover):
    nodes, edges = graph
    # generate template
    template = ""
    if cover:
        # cover = [ord(c) - ord('a') for c in cover]
        cover = [int(c) for c in cover]
        cover.sort()
        for i in range(len(cover)):
            if i == 0: template += '01' * (cover[i]) + '1'
            else: template += '01' * (cover[i] - cover[i-1]) + '1'
        template += '0' * (len(nodes) - cover[-1]) + '1'
    else:
        template += "1" + '0' * len(nodes) + '11'

    return template


def multi_seqs_to_vertex_cover(seqs):
    # seqs: a list of strings
    nodes = set()
    edges = []
    for seq in seqs:
        if not seq: continue
        left_node  = 0
        right_node = 0
        pos = 0
        for item in seq:
            if pos == 0:   # left vertex
                if item == '0': left_node += 1
                else: pos += 1
            elif pos == 1: # right vertex
                if item == '0': right_node += 1
                else: pos += 1
            else:
                assert item == '0'

        nodes.add(left_node)
        nodes.add(left_node + right_node + 1)
        # if right_node > 0:
        edges.append((left_node, left_node + right_node + 1))

    num_nodes = max(nodes) + 1

    char_nodes = []
    for i in range(num_nodes):
        # char_nodes.append(chr(i + ord('a')))
        char_nodes.append("%d" % i)

    return char_nodes, edges

def alignment(seqs, template):
    # threading
    seqs_pos = [0] * len(seqs)
    align_seqs = [""] * len(seqs)
    num_seqs = len(seqs)
    pruned_template = ""
    for item in template:
        if_has_thread = 0
        for i in range(num_seqs):
            if seqs_pos[i] < len(seqs[i]) and item == seqs[i][seqs_pos[i]]:
                if_has_thread = 1
                break

        if if_has_thread:
            for i in range(num_seqs):
                if seqs_pos[i] < len(seqs[i]) and item == seqs[i][seqs_pos[i]]:
                    seqs_pos[i] += 1
                    align_seqs[i] += item
                else: align_seqs[i] += '-'
            pruned_template += item

    for i in range(num_seqs):
        assert seqs_pos[i] == len(seqs[i]), print(i, seqs_pos[i], len(seqs[i]))

    return pruned_template, align_seqs