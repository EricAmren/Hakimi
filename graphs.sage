from itertools import combinations

class Vertex:
    def __init__(self):
        self.degree = 0
        self.index = None
        self.freeBonds = 0

    def __lt__(self, other):
        return self.freeBonds < other.free_bonds

    def __repr__(self):
        return str(self.index)


class Node:
    def __init__(self):
        self.vertices = []
        self.parent = None
        self.children = []
        self.minimumGraphs = []

    def generate_minimum_graphs(self):
        """
        Generate all possible minimum graphs of a node.
        A minimum graph is a graph that is connected with
        the minimum number of edges.
        """
        # minimumGraphs = set()
        sequence = []
        for vertex in self.vertices:
            sequence.append(vertex.degree)
        sequence = sorted(sequence, reverse=True)
        if not self.children:
            all_graphs = generate_all_graphs(len(self.vertices))
            minimumGraphs = filter_graphs(sequence, all_graphs)
            self.minimumGraphs = list(minimumGraphs)
            # self.vertices.sort(key=operator.attrgetter('degree'), reverse=True)
            # degreeSequence = []
            # hydrogenes = []
            # for vertex in self.vertices:
            #     if vertex.degree == 1:
            #         hydrogenes.append(1)
            #     else:
            #         degreeSequence.append(vertex.degree)
            # permut = list(set(permutations(degreeSequence, len(degreeSequence))))
            # for seq in permut:
            #     connected_vertices = []
            #     graph = Graph(len(self.vertices), loops=False, multiedges=True)
            #     seq = list(seq)
            #     full_seq = seq + hydrogenes
            #     for i in range(len(full_seq) - 1):
            #         # if len(connected_vertices) == len(self.vertices):
            #         #     minimumGraphs.add(graph2hash(graph))
            #         # else:
            #         for j in range(i + 1, len(full_seq)):
            #             if j not in connected_vertices:
            #                 graph.add_edge(i,j)
            #                 if i not in connected_vertices:
            #                     connected_vertices.append(i)
            #                 if j not in connected_vertices:
            #                     connected_vertices.append(j)
            #                 if len(connected_vertices) == len(self.vertices):
            #                     minimumGraphs.add(graph2hash(graph))
            # return list(minimumGraphs)


def graph2hash(graph):
    hash = graph.canonical_label().sparse6_string()
    return hash

def hash2graph(hash):
    graph = Graph(hash, loops=False, multiedges=True, data_structure="sparse")
    return graph

def generate_all_graphs(numberOfVertices):
    G = Graph(numberOfVertices, multiedges=False, loops=False)
    sol = set()
    numberOfEdge = numberOfVertices - 1
    permut = list(set(combinations(range(numberOfVertices), 2)))
    combinationOfEdges = list(combinations(permut, numberOfEdge))
    for i in combinationOfEdges:
        new_graph = copy(G)
        for edge in list(i):
            new_graph.add_edge(edge[0], edge[1])

        if new_graph.is_connected():
            sol.add(graph2hash(new_graph))
    return sol

def write_to_file(solution, filename):
    with open(filename, 'a') as f:
        for i in solution:
            f.write(i)
            f.write('\n')

# def get_number_of_hydrogen(sequence):
#     numberOfHydrogen = 0
#     for degree in sequence:
#         if degree == 1:
#             numberOfHydrogen += 1
#     return numberOfHydrogen


def filter_graphs(sequence, found_graphs):
    filtered_results = copy(found_graphs)
    for hash in found_graphs:
        graph = hash2graph(hash)
        sorted_seq = sorted(sequence)
        for i in range(len(sequence)):
            if graph.degree()[i] > sorted_seq[i]:
            # if sorted_seq[i] > graph.degree()[i]:
                filtered_results.remove(hash)
                break

    return filtered_results


def set_vertex(index, degree):
    vertex = Vertex()
    vertex.degree = degree
    vertex.index = index
    vertex.freeBonds = degree
    return vertex

def set_node(index, vertices_list):
    node = Node()
    node.index = index
    node.vertices = vertices_list
    return node

def show_graph(hash):
    show(hash2graph(hash))


C1 = set_vertex("C1", 4)
C2 = set_vertex("C2", 4)
C3 = set_vertex("C3", 4)
N1 = set_vertex("N1", 3)
N2 = set_vertex("N2", 3)
O1 = set_vertex("O1", 2)
O2 = set_vertex("O2", 2)
H1 = set_vertex("H1", 1)
H2 = set_vertex("H2", 1)

N011 = set_node("CH", [C1,H1])
N012 = set_node("NO", [N1,O1])
N021 = set_node("CO", [C2,O2])
N022 = set_node("CNH", [C3,N2,H2])
patates = [N011,N012,N021,N022]
for patate in patates:
    patate.generate_minimum_graphs()
# N02 = set_node("NNOHHHH", [N1,N1,O1,H1,H1,H1,H1])
# N01.generate_minimum_graphs()






