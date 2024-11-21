#!/usr/bin/env python3

'''Script for computing GC-rich and GC-poor intervals in a given sequence using Bellman-Ford 
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states
    -out: file to output intervals to (1 interval per line)

Outputs:
    File with list of intervals (a_i, b_i) such that bases a_i to b_i are
    classified as GC-rich.

Example Usage:
    python 2a.py -f hmm-sequence.fasta -mu 0.01 -out -intervals.txt
'''

import argparse
import numpy as np


'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	s: string with relevant sequence
'''


def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s


''' Outputs the  decoding of a given observation.
Arguments:
	obs: observed sequence of emitted states (list of emissions)
	trans_probs: transition log-probabilities (dictionary of dictionaries)
	emiss_probs: emission log-probabilities (dictionary of dictionaries)
	init_probs: initial log-probabilities for each hidden state (dictionary)
Returns:
	l: list of most likely hidden states at each position
        (list of hidden states)
	p: log-probability of the returned hidden state sequence
'''


def bellman_ford(obs, trans_probs, emiss_probs, init_probs):
    N = len(obs)  # Number of observations
    states = ['h', 'l']  # Hidden states

    # Step 1: Represent the graph as a list of edges
    edges = []
    for t in range(1, N):
        for prev_state in states:
            for current_state in states:
                # Each edge connects (prev_state, t-1) -> (current_state, t)
                weight = -trans_probs[prev_state][current_state] - emiss_probs[current_state][obs[t]]
                edges.append(((prev_state, t - 1), (current_state, t), weight))

    # Step 2: Initialize distances and backtracking
    dist = {}  # Store the shortest path cost to each node
    backtrack = {}  # Store backtracking pointers
    for state in states:
        dist[(state, 0)] = -init_probs[state] - emiss_probs[state][obs[0]]
    for t in range(1, N):
        for state in states:
            dist[(state, t)] = float('inf')  # Initialize to infinity

    # Step 3: Relax edges |V| - 1 times
    for _ in range(N - 1):
        for (u, v, weight) in edges:
            if dist[u] + weight < dist[v]:
                dist[v] = dist[u] + weight
                backtrack[v] = u

    # Step 4: Find the best final state
    final_state = min(states, key=lambda state: dist[(state, N - 1)])
    best_prob = -dist[(final_state, N - 1)]  # Convert back to positive log-probability

    # Step 5: Reconstruct the path
    path = [final_state]
    current_node = (final_state, N - 1)
    while current_node in backtrack:
        current_node = backtrack[current_node]
        path.append(current_node[0])
    path.reverse()

    return path, best_prob

''' Returns a list of non-overlapping intervals describing the GC rich regions.
Arguments:
	sequence: list of hidden states
Returns:
	intervals: list of tuples (i, j), 1 <= i <= j <= len(sequence), that
                describe GC rich regions in the input list of hidden states.
'''


def find_intervals(sequence):
    intervals = []
    in_gc_rich_region = False
    beginning = None
    N = len(sequence)

    for i, state in enumerate(sequence):
        if state == 'h':  # rich state
            if not in_gc_rich_region:
                beginning = i + 1  # start of the region (1-based index)
                in_gc_rich_region = True
        else:  # poor state
            if in_gc_rich_region:
                intervals.append((beginning, i))  # end  current region
                in_gc_rich_region = False

    # if the sequence ends while in a G+C-rich region
    if in_gc_rich_region:
        intervals.append((beginning, N))

    return intervals


def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into GC-rich and GC-poor regions using Bellman Ford')
    parser.add_argument('-f', action="store", dest="f",
                        type=str, required=True)
    parser.add_argument('-mu', action="store", dest="mu",
                        type=float, required=True)
    parser.add_argument('-out', action="store", dest="out",
                        type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    mu = args.mu
    intervals_file = args.out

    obs_sequence = read_fasta(fasta_file)
    transition_probabilities = {
        'h': {'h': np.log(1 - mu), 'l': np.log(mu)},
        'l': {'h': np.log(mu), 'l': np.log(1 - mu)}
    }
    emission_probabilities = {
        'h': {'A': np.log(0.13), 'C': np.log(0.37), 'G': np.log(0.37), 'T': np.log(0.13)},
        'l': {'A': np.log(0.32), 'C': np.log(0.18), 'G': np.log(0.18), 'T': np.log(0.32)}
    }
    initial_probabilities = {'h': np.log(0.5), 'l': np.log(0.5)}

    sequence, p = bellman_ford(obs_sequence, transition_probabilities,
                               emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end))
                for (start, end) in intervals]))
        f.write("\n")
    print("Bellman-Ford probability in log scale: {:.2f}".format(p))


if __name__ == "__main__":
    main()
