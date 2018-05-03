from itertools import combinations



class Atom:
    """
    """
    def __init__(self, degree):
        self.degree = None
        self.freeBonds = None
        # self.index = None

class Vertex:
    """
    A vertex contains a vector of Atoms.
    Each vector must have the same length in a Sgraph.
    A vertex is defined by the index of a Sgraph.
    """
    def __init__(self):
        self.atoms = []
        self.graphIndex = None

class Sgraph(Graph):
    """
    Subclass of Sage.Graph :
    A Sgraph is a link between vertices of a graph and instances of class Vertex.
    Several Sgraphs can be contained by a node.
    """
    def __init__(self, vertices):
        super().__init__(self)
        self.index = None
        self.vertices = vertices

class Node:
    """
    A Node contains all conformations of a degree sequence.
    These conformations are defined by several graphs.
    """
    def __init__(self):
        self.parent = None
        self.children = []
        self.graphs = []


v0 = Vertex()
v0.atoms.append(Atom(4), Atom(3), Atom(1), Atom(1))

v1 = Vertex()
v1.atoms.append(Atom(3), Atom(4), Atom(4), Atom(3))

v2 = Vertex()
v2.atoms.append(Atom(1), Atom(1), Atom(3), Atom(4))

G = Sgraph(3,[v0,v1,v2])








