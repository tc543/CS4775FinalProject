#!/usr/bin/env python3

import argparse
import numpy as np
import heapq


def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s


def dijkstra(obs, trans_probs, emiss_probs, init_probs):
    N = len(obs)  # Length of the observed sequence
    states = ['h', 'l']  # Hidden states

    # Priority queue for dijkstra
    pq = []  # (cost, (state, time))
    dist = {}  # Shortest distance to each node
    backtrack = {}  # Backtracking pointer for path reconstruction
    visited = set()  # Keep track of visited nodes

    # Step 1: Initialize the priority queue with initial probabilities
    for state in states:
        cost = -init_probs[state] - emiss_probs[state][obs[0]]
        heapq.heappush(pq, (cost, (state, 0)))
        dist[(state, 0)] = cost

    # Step 2: Process the priority queue
    while pq:
        current_cost, (current_state, t) = heapq.heappop(pq)

        # Skip if already visited
        if (current_state, t) in visited:
            continue
        visited.add((current_state, t))

        # Stop processing at the last time step
        if t == N - 1:
            continue

        # Get the next observation
        next_obs = obs[t + 1]

        # Step 3: Relax edges
        for next_state in states:
            weight = -trans_probs[current_state][next_state] - emiss_probs[next_state][next_obs]
            next_cost = current_cost + weight

            if (next_state, t + 1) not in dist or next_cost < dist[(next_state, t + 1)]:
                dist[(next_state, t + 1)] = next_cost
                backtrack[(next_state, t + 1)] = (current_state, t)
                heapq.heappush(pq, (next_cost, (next_state, t + 1)))

    # Step 4: Find the best final state
    final_state = min(states, key=lambda state: dist.get((state, N - 1), float('inf')))
    best_prob = -dist[(final_state, N - 1)]  # Convert back to positive log-probability

    # Step 5: Reconstruct the path
    path = [final_state]
    current_node = (final_state, N - 1)
    while current_node in backtrack:
        current_node = backtrack[current_node]
        path.append(current_node[0])
    path.reverse()

    return path, best_prob




def find_intervals(sequence):
    intervals = []
    in_gc_rich_region = False
    beginning = None
    N = len(sequence)

    for i, state in enumerate(sequence):
        if state == 'h':  # GC-rich state
            if not in_gc_rich_region:
                beginning = i + 1  # Start of the region (1-based index)
                in_gc_rich_region = True
        else:  # GC-poor state
            if in_gc_rich_region:
                intervals.append((beginning, i))  # End the current region
                in_gc_rich_region = False

    # If the sequence ends in a GC-rich region
    if in_gc_rich_region:
        intervals.append((beginning, N))

    return intervals


def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into GC-rich and GC-poor regions using Dijkstra.')
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

    sequence, p = dijkstra(obs_sequence, transition_probabilities,
                           emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end))
                for (start, end) in intervals]))
        f.write("\n")
    print("Djkstra probability in log scale: {:.2f}".format(p))


if __name__ == "__main__":
    main()
