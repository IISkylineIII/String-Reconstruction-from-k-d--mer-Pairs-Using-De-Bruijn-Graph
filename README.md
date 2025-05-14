# String-Reconstruction-from-k-d--mer-Pairs-Using-De-Bruijn-Graph


# Description
This Python script reconstructs a DNA string from a list of (k,d)-mer paired reads. It builds a De Bruijn graph from the pairs, finds an Eulerian path through the graph, and assembles the original sequence. This approach is commonly used in bioinformatics for genome assembly from paired-end sequencing data.

# Usage
Example
```
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
    in_degree = {}* 
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

# Sample Input
k = 3
d = 1
paired_reads = [
    "TAA|GCC", "AAT|CCA", "ATG|CAT", "TGC|ATG", "GCC|TGG",
    "CCA|GGG", "CAT|GGA", "ATG|GAT", "TGG|ATG", "GGG|TGT", "GGA|GTT"
]

# Sample Output
result = reconstruct_string(k, d, paired_reads)
print(result)
```
Output
GATT

# Function Descriptions
de_bruijn_graph_from_pairs(paired_reads)
* Constructs a De Bruijn graph from paired k-mers.
eulerian_path(graph)
* Finds an Eulerian path through the De Bruijn graph.
reconstruct_string(k, d, paired_reads)
* Reconstructs the DNA sequence by traversing the Eulerian path and merging overlaps.

# Applications
* Genome assembly from paired-end sequencing reads
* Demonstration of graph-based algorithms in bioinformatics
* Teaching material for Eulerian path and De Bruijn graph theory

# License
This project is licensed under the MIT License.



