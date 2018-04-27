# import bisect as *

def reverse_insort(a, x, lo=0, hi=None):
    """
    Insert item x in list a, and keep it reverse-sorted assuming a
    is reverse-sorted.

    If x is already in a, insert it to the right of the rightmost x.
    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if x > a[mid]:
            hi = mid
        else:
            lo = mid+1
    a.insert(lo, x)


class Node:
    def __init__(self, degreeSequence = None):
        self.parent = None
        self.children = []
        self.degreeSequence = [] if degreeSequence is None else degreeSequence
        self.graph = Graph(0, loops=False, multiedges=True)
        self.freeBonds = 0
        self.leaves = self.get_leaves()
        self.index = None

    def __lt__(self, other):
        return self.freeBonds < other.freeBonds

    def __repr__(self):
        return str(self.freeBonds)


    def add_children(self, list_of_children):
        for child in list_of_children:
            self.add_child(child)


    def add_child(self, child):
        if child == self:
            raise ValueError("Cannot be its own child")
        index = len(self.children)
        if child not in self.children:
            index += 1
            reverse_insort(self.children, child)
            # self.children.append(child)
            child.parent = self
            self.degreeSequence += child.degreeSequence
            child.index = index - 1
        self.degreeSequence.sort(reverse=True)

    def set_vertex(self):
        if self.is_leaf():
            graph = Graph(1, loops=False, multiedges=True)
            self.graph = graph2hash(graph)
        else:
            raise ValueError("This node is not an atom.")

    def set_graph(self):
        if self.children:
            return self.update_graph()
        else :
            return

    def is_leaf(self):
        if self.children:
            return False
        else:
            return True

    def get_graph(self):
        if self.graph != Node:
            return hash2graph(self.graph)
        else:
            raise ValueError("No graph")

    # def show_graph(self):
    #     show(self.get_graph())


    def update_graph(self):
        graph = Graph(0, loops=False, multiedges=True)
        self.graph = self.rec_update_graph(graph)

    def rec_update_graph(self, graph):
        new_graph = Graph(0, loops=False, multiedges=True)
        if self.children:
            # if self.children[0].is_leaf():
            #     for child in self.children:
            #         new_graph += child.get_graph()
            #         self.graph = graph2hash(new_graph)
            #     return new_graph
            for child in self.children:
                new_graph += child.rec_update_graph(graph)
            self.graph = graph2hash(new_graph)
            return new_graph
        else:
            graph += self.get_graph()
            return graph

    def add_virtual_edge(self, node1, node2):
        atom1 = node1.leaves[0]
        index1 = atom1.index
        atom2 = node2.leaves[0]
        index2 = atom2.index
        new_graph = copy(self.get_graph())
        new_graph.add_edge(index1, index2)

        self.graph = graph2hash(new_graph)
        # self.update_graph()
        # self.graph = graph2hash(new_graph.add_edge(vertex1, vertex2))

    def get_leaves(self):
        leaves = []
        leaves = self.rec_get_leaves(leaves)
        leaves.sort(key=operator.attrgetter('freeBonds'), reverse=True)
        return leaves

    def rec_get_leaves(self, leaves_list):
        if self.children: # If it's not a leaf
            for child in self.children:
                leaves_list.extend(child.get_leaves())
        else:
            leaves_list.append(self)
        return leaves_list







def show_graph(node):
    show(node.get_graph())

def build_atom(letter):
    letter = letter.upper()
    degree = None
    if letter == 'H':
        degree = 1
    elif letter == 'O':
        degree = 2
    elif letter == 'N':
        degree = 3
    elif letter == 'C':
        degree = 4
    if degree == None:
        raise ValueError("Character not implemented.")
    atom = Node([degree])
    atom.set_vertex()
    atom.freeBonds = degree
    return atom


def graph2hash(graph):
    hash = graph.canonical_label().sparse6_string()
    return hash

def hash2graph(hash):
    graph = Graph(hash, loops=False, multiedges=True, data_structure="sparse")
    return graph

n011 = Node()
n012 = Node()
n021 = Node()
n022 = Node()
n01 = Node()
n02 = Node()
n0 = Node()
C1 = build_atom("C")
C2 = build_atom("C")
C3 = build_atom("C")
N1 = build_atom("N")
N2 = build_atom("N")
O1 = build_atom("O")
O2 = build_atom("O")
H1 = build_atom("H")
H2 = build_atom("H")
n011.add_children([C1,H1])
n012.add_children([N1,O1])

n021.add_children([C2,O2])
n022.add_children([C3,N2,H2])
n01.add_children([n011,n012])
n02.add_children([n021,n022])
n0.add_children([n01,n02])



