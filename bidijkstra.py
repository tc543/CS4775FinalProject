import heapq
from collections import defaultdict
import math

class BiDirectionalDijkstra:
    def __init__(self, graph):
        self.graph = graph

    def dijkstra(self, source, target):
        forward_dist = defaultdict(lambda: math.inf)
        backward_dist = defaultdict(lambda: math.inf)
        forward_dist[source] = 0
        backward_dist[target] = 0

        forward_pq = [(0, source)] 
        backward_pq = [(0, target)]  

        forward_visited = set()
        backward_visited = set()

        meeting_node = None
        best_distance = math.inf

        while forward_pq or backward_pq:
            if forward_pq:
                forward_cost, forward_node = heapq.heappop(forward_pq)
                if forward_node in forward_visited:
                    continue
                forward_visited.add(forward_node)

                # If the node is visited in both directions
                if forward_node in backward_visited:
                    total_distance = forward_dist[forward_node] + backward_dist[forward_node]
                    if total_distance < best_distance:
                        best_distance = total_distance
                        meeting_node = forward_node

                for neighbor, weight in self.graph[forward_node]:
                    if neighbor not in forward_visited:
                        new_dist = forward_cost + weight
                        if new_dist < forward_dist[neighbor]:
                            forward_dist[neighbor] = new_dist
                            heapq.heappush(forward_pq, (new_dist, neighbor))

            if backward_pq:
                backward_cost, backward_node = heapq.heappop(backward_pq)
                if backward_node in backward_visited:
                    continue
                backward_visited.add(backward_node)

                # If the node is visited in both directions
                if backward_node in forward_visited:
                    total_distance = forward_dist[backward_node] + backward_dist[backward_node]
                    if total_distance < best_distance:
                        best_distance = total_distance
                        meeting_node = backward_node

                for neighbor, weight in self.graph[backward_node]:
                    if neighbor not in backward_visited:
                        new_dist = backward_cost + weight
                        if new_dist < backward_dist[neighbor]:
                            backward_dist[neighbor] = new_dist
                            heapq.heappush(backward_pq, (new_dist, neighbor))

        return best_distance, meeting_node


# TESTING 
graph1 = {
    0: [(1, 4), (2, 1)],
    1: [(3, 1)],
    2: [(1, 2), (3, 5)],
    3: []
}
source1 = 0
target1 = 3
# Expected Output: Shortest distance = 4 (Path: 0 -> 2 -> 1 -> 3)

graph2 = {
    0: [(1, 2)],
    1: [(2, 3)],
    2: [],
    3: [(4, 1)],
    4: []
}
source2 = 0
target2 = 4
# Expected Output: Shortest distance = inf (No path exists)

graph3 = {
    0: [(1, 2), (3, 10)],
    1: [(2, 4)],
    2: [(4, 3)],
    3: [(4, 6)],
    4: []
}
source3 = 0
target3 = 4
# Expected Output: Shortest distance = 9 (Path: 0 -> 1 -> 2 -> 4)

graph4 = {
    0: [(1, 1), (2, 1)],
    1: [(3, 1)],
    2: [(3, 1)],
    3: []
}
source4 = 0
target4 = 3
# Expected Output: Shortest distance = 2 (Path: 0 -> 1 -> 3 or 0 -> 2 -> 3)

graph5 = {
    0: [(1, 1)],
    1: [(2, 1)],
    2: [(0, 1), (3, 1)],
    3: []
}
source5 = 0
target5 = 3
# Expected Output: Shortest distance = 3 (Path: 0 -> 1 -> 2 -> 3)

graph6 = {
    0: []
}
source6 = 0
target6 = 0
# Expected Output: Shortest distance = 0 (Path: 0)

graph7 = {
    0: [(1, 1), (2, 4)],
    1: [(3, 3)],
    2: [(3, 2)],
    3: []
}
source7 = 0
target7 = 3
# Expected Output: Shortest distance = 4 (Path: 0 -> 1 -> 3)

graph8 = {
    0: [(1, 2), (2, -3)],
    1: [(3, 2)],
    2: [(3, 1)],
    3: []
}
source8 = 0
target8 = 3
# Expected Output: Algorithm should not support negative weights 


graph9 = {
    0: [(1, 2), (2, 4), (3, 6)],
    1: [(2, 1), (3, 3)],
    2: [(3, 2)],
    3: []
}
source9 = 0
target9 = 3
# Expected Output: Shortest distance = 5 (Path: 0 -> 1 -> 3)

graph10 = {
    0: [(1, 2), (2, 5)],
    1: [(3, 1)],
    2: [(3, 1)],
    3: [(4, 3)],
    4: []
}
source10 = 0
target10 = 4
# Expected Output: Shortest distance = 6 (Path: 0 -> 1 -> 3 -> 4)

test_cases = [
    (graph1, source1, target1),
    (graph2, source2, target2),
    (graph3, source3, target3),
    (graph4, source4, target4),
    (graph5, source5, target5),
    (graph6, source6, target6),
    (graph7, source7, target7),
    (graph8, source8, target8),
    (graph9, source9, target9),
    (graph10, source10, target10),
]

for i, (graph, source, target) in enumerate(test_cases, 1):
    bidirectional_dijkstra = BiDirectionalDijkstra(graph)
    try:
        distance, meeting_node = bidirectional_dijkstra.dijkstra(source, target)
        print(f"Test Case {i}: Shortest distance from {source} to {target} is {distance}, Meeting Node: {meeting_node}")
    except Exception as e:
        print(f"Test Case {i}: Failed with error {e}")
