import argparse
import heapq
import numpy as np

def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s

def bidirectional_dijkstra(obs, trans_probs, emiss_probs, init_probs):
    N = len(obs)
    states = ['h', 'l']

    # Priority queues for forward and backward searches
    forward_pq = []  
    backward_pq = []  

    # Distances and backtracking pointers
    forward_dist = {}
    backward_dist = {}
    forward_backtrack = {}
    backward_backtrack = {}

    # Visited nodes
    forward_visited = set()
    backward_visited = set()

    # Initialize priority queues with initial probabilities
    for state in states:
        forward_cost = -init_probs[state] - emiss_probs[state][obs[0]]
        backward_cost = -init_probs[state] - emiss_probs[state][obs[-1]]

        heapq.heappush(forward_pq, (forward_cost, (state, 0)))
        heapq.heappush(backward_pq, (backward_cost, (state, N - 1)))

        forward_dist[(state, 0)] = forward_cost
        backward_dist[(state, N - 1)] = backward_cost

    best_distance = float('inf')
    meeting_node = None

    # Bidirectional search
    while forward_pq or backward_pq:
        # Forward step
        if forward_pq:
            current_cost, (current_state, t) = heapq.heappop(forward_pq)

            if (current_state, t) in forward_visited:
                continue
            forward_visited.add((current_state, t))

            # Check for meeting point
            if (current_state, t) in backward_dist:
                total_cost = current_cost + backward_dist[(current_state, t)]
                if total_cost < best_distance:
                    best_distance = total_cost
                    meeting_node = (current_state, t)

            if t < N - 1:
                next_obs = obs[t + 1]
                for next_state in states:
                    weight = -trans_probs[current_state][next_state] - emiss_probs[next_state][next_obs]
                    next_cost = current_cost + weight

                    if (next_state, t + 1) not in forward_dist or next_cost < forward_dist[(next_state, t + 1)]:
                        forward_dist[(next_state, t + 1)] = next_cost
                        forward_backtrack[(next_state, t + 1)] = (current_state, t)
                        heapq.heappush(forward_pq, (next_cost, (next_state, t + 1)))

        # Backward step
        if backward_pq:
            current_cost, (current_state, t) = heapq.heappop(backward_pq)

            if (current_state, t) in backward_visited:
                continue
            backward_visited.add((current_state, t))

            # Check for meeting point
            if (current_state, t) in forward_dist:
                total_cost = current_cost + forward_dist[(current_state, t)]
                if total_cost < best_distance:
                    best_distance = total_cost
                    meeting_node = (current_state, t)

            if t > 0:
                prev_obs = obs[t - 1]
                for prev_state in states:
                    weight = -trans_probs[prev_state][current_state] - emiss_probs[prev_state][prev_obs]
                    next_cost = current_cost + weight

                    if (prev_state, t - 1) not in backward_dist or next_cost < backward_dist[(prev_state, t - 1)]:
                        backward_dist[(prev_state, t - 1)] = next_cost
                        backward_backtrack[(prev_state, t - 1)] = (current_state, t)
                        heapq.heappush(backward_pq, (next_cost, (prev_state, t - 1)))

    # Reconstruct the path
    if meeting_node:
        forward_path = []
        current = meeting_node
        while current in forward_backtrack:
            forward_path.append(current[0])
            current = forward_backtrack[current]
        forward_path.reverse()

        backward_path = []
        current = meeting_node
        while current in backward_backtrack:
            current = backward_backtrack[current]
            backward_path.append(current[0])

        path = forward_path + [meeting_node[0]] + backward_path
    else:
        path = []

    return path, -best_distance

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

    if in_gc_rich_region:
        intervals.append((beginning, N))

    return intervals

def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into GC-rich and GC-poor regions using Bidirectional Dijkstra.'
    )
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-mu', action="store", dest="mu", type=float, required=True)
    parser.add_argument('-out', action="store", dest="out", type=str, required=True)

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

    sequence, p = bidirectional_dijkstra(obs_sequence, transition_probabilities, emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)

    with open(intervals_file, "w") as f:
        f.write("\n".join(["%d,%d" % (start, end) for (start, end) in intervals]))
        f.write("\n")

    print("Bidirectional Dijkstra probability in log scale: {:.2f}".format(p))

if __name__ == "__main__":
    main()
