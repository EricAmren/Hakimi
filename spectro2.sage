from operator import attrgetter
from itertools import combinations

def TODO():
    raise NotImplementedError("To be implemented")



class Node:
    def __init__(self, atoms = []):
        self.parent = None
        self.children = []
        self.atoms = atoms

    def add_children(self, children):
        for child in children:
            if child not in self.children:
                self.children.append(child)
                child.parent = self
    def is_root(self):
        return self.parent == None

    def is_connected(self, G):
        vertices = [atom.vertex for atom in self.atoms]
        subgraph = G.subgraph(vertices)
        return subgraph.is_connected()

class Atom:
    def __init__(self):
        self.vertex = None  # Index of this atom in general graph
        self.valence = None
        self.freeBonds = 0

    def __gt__(self, other):
        return self.freeBonds > other.freeBonds

def build_atom(valence):
    atom = Atom()
    atom.valence = valence
    atom.freeBonds = valence
    return atom

def connect_2_atoms(G, atom1, atom2):
    assert atom1.freeBonds > 0, "atom1 cannot have more bonds"
    assert atom2.freeBonds > 0, "atom2 cannot have more bonds"
    G.add_edge(atom1.vertex, atom2.vertex)
    atom1.freeBonds -= 1
    atom2.freeBonds -= 1

def connect_leaves2(list_of_leaves, G):
    for leaf in list_of_leaves:
        connected_atoms = set()
        while not leaf.is_connected(G):
            sorted_atoms = sorted(leaf.atoms, key=attrgetter('freeBonds'), reverse=True)
            all_pairs = combinations(sorted_atoms, 2)
            for pair in all_pairs:
                atom1 = pair[0]
                atom2 = pair[1]
                if atom1 in connected_atoms and atom2 in connected_atoms:
                    pass
                # elif not atom1 in connected_atoms and not atom2 in connected_atoms:
                #     pass
                else:
                    connect_2_atoms(G, atom1, atom2)
                    connected_atoms.update(pair)
                    break


def connect_leaves(list_of_leaves, G):
    for leaf in list_of_leaves:
        free_atoms = sorted(leaf.atoms, key=attrgetter('freeBonds'),reverse=True)
        max_atom = free_atoms.pop(0)
        connected_atoms = [max_atom]

        while free_atoms:
            max_connected_atom = max(connected_atoms, key=attrgetter('freeBonds'))
            if max_connected_atom.freeBonds == 0:
                raise ValueError("Not enough free atoms.")
            max_free_atom = free_atoms.pop(0)
            connect_2_atoms(G, max_connected_atom, max_free_atom)
            connected_atoms.append(max_free_atom)

def connect_node(node, G):
    free_children = sorted(node.children, key=attrgetter('atoms'), reverse=True)
    max_child = free_children.pop(0)
    connected_children = [max_child]

    while free_children:
        max_connected_child = max(connected_children, key=attrgetter('atoms'))
        if max(max_connected_child.atoms).freeBonds == 0:
            raise ValueError("Not enough free atoms.")
        max_free_child = free_children.pop(0)
        connect_2_atoms(G, max(max_connected_child.atoms), max(max_free_child.atoms))
        connected_children.append(max_free_child)

def connect_node2(node, G):
    connected_children = set()
    while not node.is_connected(G):
        max_atoms_of_children = []
        for child in node.children:
            max_atom = max(child.atoms, key=attrgetter('freeBonds'))
            if max_atom.freeBonds > 0:
                max_atoms_of_children.append((max_atom, child))
            else:
                pass
        max_atoms_of_children.sort(reverse=True)
        all_pairs = combinations(max_atoms_of_children, 2)
        for pair in all_pairs:
            atom1 = pair[0][0]
            atom2 = pair[1][0]
            child1 = pair[0][1]
            child2 = pair[1][1]
            if len(connected_children) == 0:
                connected_children.add(child1)
            if child1 in connected_children and child2 in connected_children:
                pass
            elif not child1 in connected_children and not child2 in connected_children:
                pass
            else:
                connect_2_atoms(G, atom1, atom2)
                connected_children.add(child1)
                connected_children.add(child2)
                break

def fill_remaining_bonds(root,G):
    unfilled_atoms = []
    for atom in root.atoms:
        if atom.freeBonds != 0:
            unfilled_atoms.append(atom)
    while unfilled_atoms != []:
        sorted_atoms = sorted(unfilled_atoms, key=attrgetter('freeBonds'), reverse=True)
        atom1 = sorted_atoms[0]
        atom2 = sorted_atoms[1]
        connect_2_atoms(G, atom1, atom2)
        unfilled_atoms = []
        for atom in root.atoms:
            if atom.freeBonds != 0:
                unfilled_atoms.append(atom)

def generate_first_graph(leaves, upper_nodes):
    G = Graph(0, loops=False, multiedges=True)
    index = 0
    for leaf in leaves:
        for atom in leaf.atoms:
            G.add_vertex(index)
            atom.vertex = index
            # G.set_vertex(index, atom.valence)
            index += 1
    connect_leaves(leaves, G)
    for node in upper_nodes:
        connect_node(node, G)
    fill_remaining_bonds(root, G)
    return G


class Fragmentation_Tree:
    def __init__(self, sequence = None):
        self.sequence = sequence
        self.leaves = []
        self.upper_nodes = []
        self.root = Node()



    def set_nodes(self):
        root = self.root
        branches = self.sequence
        self.rec_set_nodes(root, branches)
        # return ()

    def rec_set_nodes(self, node, branches):
        for branch in branches:
            if type(branch) == sage.rings.integer.Integer:
                degree = branch
                new_atom = Atom()
                new_atom.valence = degree
                new_atom.freeBonds = degree
                node.atoms.append(new_atom)
                if node not in self.leaves:
                    self.leaves.append(node)
            else:
                if not node in self.upper_nodes:
                    self.upper_nodes.append(node)
                child_node = Node()
                self.rec_set_nodes(child_node, branch)
                for atom in child_node.atoms:
                    if atom not in node.atoms:
                        node.atoms.append(atom)
                node.add_children([child_node])

def build_fragmentation_tree(sequence):
    tree = Fragmentation_Tree(sequence)
    tree.set_nodes()
    return tree

T = [[[4],[3],[2,1]],[[3,2],[3,1]],[[4,1]]]

###
# def build_fragmentation_tree(tree):
#     # set_atoms(tree)
#     set_nodes(tree)
#     # leaves, upper_nodes = set_nodes(tree)

# def set_nodes(tree):
#     root = Node()
#     leaves = []
#     upper_nodes = []
#     rec_set_nodes(root, tree, leaves, upper_nodes)
#     return ()

# def rec_set_nodes(root, tree, leaves, upper_nodes):
#     for branch in tree:
#         if type(branch) == sage.rings.integer.Integer:
#             root.atoms.append(build_atom(branch))
#             leaves.append(root)
#         else:
#             upper_nodes.append(root)
#             new_node = Node()
#             rec_set_nodes(new_node, branch, leaves, upper_nodes)
#             for atom in new_node.atoms:
#                 if atom not in root.atoms:
#                     root.atoms.append(atom)
#             root.add_children([new_node])



# def set_atoms(tree):
#     atoms = []
#     rec_set_atoms(tree, atoms)
#     return atoms

# def rec_set_atoms(node1, atoms):
#     for node in node1:
#         if type(node) == sage.rings.integer.Integer:
#             atoms.append(build_atom(node))
#         elif type(node) is list:
#             rec_set_atoms(node, atoms)



def relabel_graph(graph):
    keys = graph.vertices()
    values = graph.degree()
    dict1 = dict(zip(keys, values))
    for key, value in dict1.items():
        dict1[key] = str(key) + '-' + str(value)
    graph.relabel(dict1)


################################ TESTING SETS ###################################
## Set1
# C1 = build_atom(4)
# C2 = build_atom(4)
# C3 = build_atom(4)
# N1 = build_atom(3)
# N2 = build_atom(3)
# O1 = build_atom(2)
# O2 = build_atom(2)
# H1 = build_atom(1)
# H2 = build_atom(1)

# n1 = Node([C1,H1])
# n2 = Node([N1,O1])
# n3 = Node([C2,O2])
# n4 = Node([C3,N2,H2])
# G1 = Node(n1.atoms + n2.atoms)
# G1.add_children([n1,n2])
# G2 = Node(n3.atoms + n4.atoms)
# G2.add_children([n3,n4])
# root = Node(G1.atoms + G2.atoms)
# root.add_children([G1,G2])
# leaves = [n1,n2,n3,n4]
# upper_nodes = [G1,G2,root]

## Set2
C1 = build_atom(4)
C2 = build_atom(4)
N1 = build_atom(3)
N2 = build_atom(3)
N3 = build_atom(3)
O1 = build_atom(2)
O2 = build_atom(2)
H1 = build_atom(1)
H2 = build_atom(1)
H3 = build_atom(1)

n1 = Node([C1])
n2 = Node([N1])
n3 = Node([O1,H1])
n4 = Node([N2,O2])
n5 = Node([N3,H2])
n6 = Node([C2,H3])

G1 = Node(n1.atoms + n2.atoms + n3.atoms)
G1.add_children([n1,n2,n3])
G2 = Node(n4.atoms + n5.atoms)
G2.add_children([n4,n5])
G3 = Node(n6.atoms)
G3.add_children([n6])
root = Node(G1.atoms + G2.atoms + G3.atoms)
root.add_children([G1,G2,G3])
leaves = [n1,n2,n3,n4,n5,n6]
upper_nodes = [G1,G2,G3,root]


G = generate_first_graph(leaves, upper_nodes)
# show(G)

## Second step : Generate all graphs

def generate_all_graphs(G):
    valid_hashes = {graph2hash(G)}
    valid_graphs = [G]
    stack = [G]
    while len(stack) > 0:
        new_G = stack.pop()
        rec_generate_all_graphs(new_G, valid_hashes, stack, valid_graphs)
        print(len(valid_hashes))
    # return valid_hashes
    return valid_graphs

def rec_generate_all_graphs(G, valid_hashes, stack, valid_graphs):
    all_pairs_of_edges = combinations(G.edges(), 2)
    for pair in all_pairs_of_edges:
        edge1 = pair[0]
        edge2 = pair[1]
        if edge1[0] in edge2 or edge1[1] in edge2:
            pass
        else :
            new_G1 = copy(G)
            new_G2 = copy(G)
            straight_switch(new_G1, pair[0], pair[1])
            crossed_switch(new_G2, pair[0], pair[1])
            all_nodes = leaves + upper_nodes
            for new_G in [new_G1, new_G2]:
                new_G_is_valid = True
                for node in all_nodes:
                    if not node.is_connected(new_G):
                        new_G_is_valid = False
                if new_G_is_valid:
                    if not graph2hash(new_G) in valid_hashes:
                        stack.append(new_G)
                        valid_hashes.add(graph2hash(new_G))
                        valid_graphs.append(new_G)
    return (valid_hashes, stack, valid_graphs)

def switch_all_edges(G):
    all_pairs_of_edges = combinations(G.edges(), 2)
    valid_graphs = [G]
    # valid_hashes = {graph2hash(G)}
    # stack = {graph2hash(G)}
    for pair in all_pairs_of_edges:
        edge1 = pair[0]
        edge2 = pair[1]
        if edge1[0] in edge2 or edge1[1] in edge2:
            pass
        else :
            new_G1 = copy(G)
            new_G2 = copy(G)
            straight_switch(new_G1, pair[0], pair[1])
            crossed_switch(new_G2, pair[0], pair[1])
            all_nodes = leaves + upper_nodes
            new_G1_is_valid = True
            new_G2_is_valid = True
            for node in all_nodes:
                if not node.is_connected(new_G1):
                    new_G1_is_valid = False
                if not node.is_connected(new_G2):
                    new_G2_is_valid = False
            if new_G1_is_valid:
                valid_graphs.append(new_G1)
            if new_G2_is_valid:
                valid_graphs.append(new_G2)
    return valid_graphs

def straight_switch(G, edge1, edge2):
    assert edge1 in G.edges(), "edge1 is not in G"
    assert edge2 in G.edges(), "edge2 is not in G"
    v1 = edge1[0]
    v2 = edge1[1]
    v3 = edge2[0]
    v4 = edge2[1]
    assert v1 != v2 and v1 != v3 and v1 != v4 and v2 != v3 and v2 != v4 and v3 != v4, "edges can't share a vertex."
    G.delete_edge(edge1)
    G.delete_edge(edge2)
    G.add_edge(v1,v3)
    G.add_edge(v2,v4)

def crossed_switch(G, edge1, edge2):
    assert edge1 in G.edges(), "edge1 is not in G"
    assert edge2 in G.edges(), "edge2 is not in G"
    v1 = edge1[0]
    v2 = edge1[1]
    v3 = edge2[0]
    v4 = edge2[1]
    assert v1 != v2 and v1 != v3 and v1 != v4 and v2 != v3 and v2 != v4 and v3 != v4, "edges can't share a vertex."
    G.delete_edge(edge1)
    G.delete_edge(edge2)
    G.add_edge(v1,v4)
    G.add_edge(v2,v3)

def graph2hash(graph):
    hash = graph.canonical_label().sparse6_string()
    return hash

def hash2graph(hash):
    graph = Graph(hash, loops=False, multiedges=True, data_structure="sparse")
    return graph



# sol = generate_all_graphs(G)
