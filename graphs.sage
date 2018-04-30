from itertools import permutations

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
        minimumGraphs = []
        if not self.children:
            self.vertices.sort(key=operator.attrgetter('degree'), reverse=True)
            degreeSequence = []
            hydrogenes = []
            for vertex in self.vertices:
                if vertex.degree == 1:
                    hydrogenes.append(1)
                else:
                    degreeSequence.append(vertex.degree)
            permut = permutations(degreeSequence, len(degreeSequence))
            for seq in permut:
                graph = Graph(len(self.vertices), loops=False, multiedges=True)
                seq = list(seq)
                full_seq = seq + hydrogenes
                # print(full_seq)
                # numberOfMinBond = len(full_seq - 1)
                # while numberOfMinBond > 0:
                #     for i in range(full_seq[0]):
                #         graph.add_edge(full_seq[0])
                #     full_seq[0]



                #     for degree in full_seq:
                #         for i in range(min(len(full_seq) - 1, degree)):
                #             if


            # graph = Graph(len(self.vertices), loops=False, multiedges=True)
            # self.vertices.sort(key=operator.attrgetter('freeBonds'), reverse=True)
            # for vertex in self.vertices:
            #     graph.relabel()
            # for vertex_ID in range(len(self.vertices)):
            #     graph.relabel({vertex_ID:self.vertices[vertex_ID].index})
            # self.minimumGraphs.append(graph)

def graph2hash(graph):
    hash = graph.canonical_label().sparse6_string()
    return hash

def hash2graph(hash):
    graph = Graph(hash, loops=False, multiedges=True, data_structure="sparse")
    return graph

def generate_all(sequence):
    G = Graph(len(sequence), multiedges=False, loops=False)
    sol = set()
    numberOfEdge = len(sequence) - 1
    permut = list(set(permutations(sequence, 2)))
    combinationOfEdges = list(permutations(permut, 3))
    for i in combinationOfEdges:
        new_graph = copy(G)
        for edge in list(i):
            new_graph.add_edge(edge[0], edge[1])

        if new_graph.is_connected():
            sol.add(graph2hash(new_graph))
    return sol


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


C1 = set_vertex("C1", 4)
H1 = set_vertex("H1", 1)
N1 = set_vertex("N1", 3)

N01 = set_node("CH", [H1, N1, C1])
N01.generate_minimum_graphs()






