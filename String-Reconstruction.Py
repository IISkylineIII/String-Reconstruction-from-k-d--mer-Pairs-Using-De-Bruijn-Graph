def de_bruijn_graph_from_pairs(paired_reads):
    graph = {}
    for pair in paired_reads:
        kmer1, kmer2 = pair.split('|')
        prefix = kmer1[:-1] + kmer2[:-1]
        suffix = kmer1[1:] + kmer2[1:]
        if prefix not in graph:
            graph[prefix] = [suffix]
        else:
            graph[prefix].append(suffix)
    return graph

def eulerian_path(graph):
    in_degree = {}
    out_degree = {}

    for node, neighbors in graph.items():
        out_degree[node] = len(neighbors)
        if node not in in_degree:
            in_degree[node] = 0
        for neighbor in neighbors:
            if neighbor not in in_degree:
                in_degree[neighbor] = 1
            else:
                in_degree[neighbor] += 1

    start_node = None
    end_node = None
    for node in in_degree:
        if node not in out_degree or in_degree[node] < out_degree[node]:
            start_node = node
        elif in_degree[node] > out_degree[node]:
            end_node = node

    stack = [start_node]
    path = []

    while stack:
        current_node = stack[-1]
        if current_node in graph and graph[current_node]:
            next_node = graph[current_node].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())

    return path[::-1]

def reconstruct_string(k, d, paired_reads):
    graph = de_bruijn_graph_from_pairs(paired_reads)
    path = eulerian_path(graph)

    reconstructed_string = path[0]
    for node in path[1:]:
        reconstructed_string += node[-1]

    return reconstructed_string

# Input
k = 3
d = 1
paired_reads = ["TAA|GCC","AAT|CCA","ATG|CAT","TGC|ATG","GCC|TGG","CCA|GGG","CAT|GGA","ATG|GAT","TGG|ATG","GGG|TGT","GGA|GTT"]

# Output
result = reconstruct_string(k, d, paired_reads)
print(result)
