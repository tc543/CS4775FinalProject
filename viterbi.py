#!/usr/bin/env python3

'''Script for computing GC-rich and GC-poor intervals in a given sequence.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states
    -out: file to output intervals to (1 interval per line)

Outputs:
    File with list of intervals (a_i, b_i) such that bases a_i to b_i are
    classified as GC-rich.
    
Example Usage:
    python 2a.py -f hmm-sequence.fasta -mu 0.01 -out viterbi-intervals.txt
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


''' Outputs the Viterbi decoding of a given observation.
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


def viterbi(obs, trans_probs, emiss_probs, init_probs):
    N = len(obs)  # length of the observed sequence
    states = ['h', 'l']  # hidden states

    # First, we initialize DP tables for log-probabilities and backtracking pointers
    dp = {state: [-np.inf] * len(obs) for state in states}
    backpointer = {state: [None] * len(obs) for state in states}

    # Next, we initialize base cases (t = 0)
    for state in states:
        dp[state][0] = init_probs[state] + emiss_probs[state][obs[0]]

    # Then, we fill the table by iterating through the sequence
    for t in range(1, len(obs)):
        for current_state in states:
            max_prob = -np.inf
            best_prev_state = None

            # check the probabilities of transitioning
            for prev_state in states:
                prob = dp[prev_state][t-1] + trans_probs[prev_state][current_state] + emiss_probs[current_state][obs[t]]
                if prob > max_prob:
                    max_prob = prob
                    best_prev_state = prev_state

            # store the best probability and  backpointer
            dp[current_state][t] = max_prob
            backpointer[current_state][t] = best_prev_state

    # finally, backtrack to get the most likely hidden state sequence
    final_state = max(states, key=lambda state: dp[state][len(obs)-1])
    most_likely_sequence = [final_state]

    for t in range(len(obs)-1, 0, -1):
        final_state = backpointer[final_state][t]
        most_likely_sequence.append(final_state)

    most_likely_sequence.reverse()

    # return the most likely sequence and the final log-probability
    final_log_prob = max(dp['h'][-1], dp['l'][-1])
    return most_likely_sequence, final_log_prob


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
        description='Parse a sequence into GC-rich and GC-poor regions using Viterbi.')
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
    sequence, p = viterbi(obs_sequence, transition_probabilities,
                          emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end))
                for (start, end) in intervals]))
        f.write("\n")
    print("Viterbi probability in log scale: {:.2f}".format(p))


if __name__ == "__main__":
    main()
