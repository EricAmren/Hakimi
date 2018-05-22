import re
from operator import attrgetter
from itertools import combinations

valence_dict = dict()
valence_dict['H'] = 1
valence_dict['O'] = 2
valence_dict['N'] = 3
valence_dict['C'] = 4
valence_dict['P'] = 5
valence_dict['S'] = 6

def draw_mol(G):
    dic_colors={}
    dic_colors["black"]=[i for i in G if G.degree(i)==4]
    dic_colors["red"]=[i for i in G if G.degree(i)==2]
    dic_colors["white"]=[i for i in G if G.degree(i)==1]
    dic_colors["blue"]=[i for i in G if G.degree(i)==3]
    if G.has_multiple_edges():
        G.plot(vertex_colors=dic_colors).show()
    else:


        G.plot3d(vertex_colors=dic_colors).show()

def formula2degrees(formula):
    """
    From a chemical formula, return the corresponding degree sequence as a list.
    """
    splitted_formula = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    # print(splitted_formula)
    degree_sequence = []
    for atom in splitted_formula:
        if atom[1] == '':
            degree = valence_dict[atom[0]]
            degree_sequence.append(degree)
        else:
            for i in range(int(atom[1])):
                degree = valence_dict[atom[0]]
                degree_sequence.append(degree)
    degree_sequence.sort(reverse=True)
    return degree_sequence

    # atom = None
    # number = ''
    # for i in formula:
    #     if formula[i].isalpha():
    #         dict1[atom] = number
    #         atom = formula[i]
    #     if formula[i].isdigit():
    #         number += formula[i]



class Node:
    def __init__(self, atoms = {}):
        self.parent = None
        self.children = []
        self.atoms = list(atoms)

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

# def connect_leaves2(list_of_leaves, G):
#     for leaf in list_of_leaves:
#         connected_atoms = set()
#         while not leaf.is_connected(G):
#             sorted_atoms = sorted(leaf.atoms, key=attrgetter('freeBonds'), reverse=True)
#             all_pairs = combinations(sorted_atoms, 2)
#             for pair in all_pairs:
#                 atom1 = pair[0]
#                 atom2 = pair[1]
#                 if atom1 in connected_atoms and atom2 in connected_atoms:
#                     pass
#                 elif not atom1 in connected_atoms and not atom2 in connected_atoms:
#                     pass
#                 else:
#                     connect_2_atoms(G, atom1, atom2)
#                     connected_atoms.update(pair)
#                     break


def connect_leaves(list_of_leaves, G):
    for leaf in list_of_leaves:
        free_atoms = sorted(leaf.atoms, key=attrgetter('freeBonds'),reverse=True)
        max_atom = free_atoms.pop(0)
        connected_atoms = [max_atom]

        while free_atoms:
            max_connected_atom = max(connected_atoms, key=attrgetter('freeBonds'))
            if max_connected_atom.freeBonds == 0:
                raise ValueError("Not enough free atoms.")
            free_atoms = sorted(free_atoms, key=attrgetter('freeBonds'), reverse=True)
            max_free_atom = free_atoms.pop(0)
            connect_2_atoms(G, max_connected_atom, max_free_atom)
            connected_atoms.append(max_free_atom)
        assert leaf.is_connected(G), "Leaves are not connected!"

def connect_node(node, G):  ####TODO#####"
    free_children = node.children
    max_child = get_max_atom_of_nodes(node.children)
    # max_child is a tuple containing the atom with the most of free bonds
    # and the node child containing this atom.
    connected_children = [max_child]
    free_children.remove(max_child[1])
    while free_children:
        max_connected_child = max(connected_children)
        if max_connected_child[0].freeBonds == 0:
            raise ValueError("Not enough free atoms.")
        max_free_child = get_max_atom_of_nodes(free_children)
        if max_free_child[0].freeBonds == 0:
            raise ValueError("Not enough free atoms!")
        connect_2_atoms(G, max_connected_child[0], max_free_child[0])
        connected_children.append(max_free_child)
        free_children.remove(max_free_child[1])




def get_max_atom_of_nodes(node_list):
    max_atoms = []
    for node in node_list:
        max_atom = max(node.atoms, key=attrgetter('freeBonds'))
        max_atoms.append((max_atom, node))
    return max(max_atoms)

def connect_node2(node, G):
    free_children = sorted(node.children, key=attrgetter('atoms'), reverse=True)
    max_child = free_children.pop(0)
    connected_children = [max_child]

    while free_children:
        max_connected_child = max(connected_children, key=attrgetter('atoms'))
        if max(max_connected_child.atoms).freeBonds == 0:
            raise ValueError("Not enough free atoms.")
        max_free_child = free_children.pop(0)
        try:
            connect_2_atoms(G, max(max_connected_child.atoms), max(max_free_child.atoms))
        except AssertionError:
            print("Error in connect_node")
            print(max(max_connected_child.atoms).valence, max(max_free_child.atoms).valence)
            print(len(node.children[0].children[0].children[0].children[0].children))
            exit(0)
        connected_children.append(max_free_child)

# def connect_node2(node, G):
#     connected_children = set()
#     while not node.is_connected(G):
#         max_atoms_of_children = []
#         for child in node.children:
#             max_atom = max(child.atoms, key=attrgetter('freeBonds'))
#             if max_atom.freeBonds > 0:
#                 max_atoms_of_children.append((max_atom, child))
#             else:
#                 pass
#         max_atoms_of_children.sort(reverse=True)
#         all_pairs = combinations(max_atoms_of_children, 2)
#         for pair in all_pairs:
#             atom1 = pair[0][0]
#             atom2 = pair[1][0]
#             child1 = pair[0][1]
#             child2 = pair[1][1]
#             if len(connected_children) == 0:
#                 connected_children.add(child1)
#             if child1 in connected_children and child2 in connected_children:
#                 pass
#             elif not child1 in connected_children and not child2 in connected_children:
#                 pass
#             else:
#                 connect_2_atoms(G, atom1, atom2)
#                 connected_children.add(child1)
#                 connected_children.add(child2)
#                 break

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

def generate_first_graph(tree):
    G = Graph(0, loops=False, multiedges=True)
    index = 0
    for leaf in tree.leaves:
        for atom in leaf.atoms:
            G.add_vertex(index)
            atom.vertex = index
            # G.set_vertex(index, atom.valence)
            index += 1
    connect_leaves(tree.leaves, G)
    for node in tree.upper_nodes:
        connect_node(node, G)
    fill_remaining_bonds(tree.root, G)
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
        self.upper_nodes = list(reversed(self.upper_nodes))
        # return ()

    def rec_set_nodes(self, node, branches):
        for branch in branches:
            if type(branch) == sage.rings.integer.Integer:
                new_atom = build_atom(branch)
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
# C1 = build_atom(4)
# C2 = build_atom(4)
# N1 = build_atom(3)
# N2 = build_atom(3)
# N3 = build_atom(3)
# O1 = build_atom(2)
# O2 = build_atom(2)
# H1 = build_atom(1)
# H2 = build_atom(1)
# H3 = build_atom(1)

# n1 = Node({C1})
# n2 = Node({N1})
# n3 = Node({O1,H1})
# n4 = Node({N2,O2})
# n5 = Node({N3,H2})
# n6 = Node({C2,H3})

# G1 = Node(n1.atoms + n2.atoms + n3.atoms)
# G1.add_children([n1,n2,n3])
# G2 = Node(n4.atoms + n5.atoms)
# G2.add_children([n4,n5])
# G3 = Node(n6.atoms)
# G3.add_children([n6])
# root = Node(G1.atoms + G2.atoms + G3.atoms)
# root.add_children([G1,G2,G3])
# leaves = [n1,n2,n3,n4,n5,n6]
# upper_nodes = [G1,G2,G3,root]


# G = generate_first_graph(leaves, upper_nodes)
# show(G)

## Second step : Generate all graphs

def generate_all_graphs(tree):
    G = generate_first_graph(tree)
    valid_hashes = {graph2hash(G)}
    valid_graphs = [G]
    stack = [G]
    while len(stack) > 0:
        new_G = stack.pop()
        rec_generate_all_graphs(tree, new_G, valid_hashes, stack, valid_graphs)
        print(len(valid_hashes))
    # return valid_hashes
    return valid_graphs

def rec_generate_all_graphs(tree, G, valid_hashes, stack, valid_graphs):
    all_pairs_of_edges = combinations(G.edges(), 2)
    for pair in all_pairs_of_edges:
        edge1 = pair[0]
        edge2 = pair[1]
        if edge1[0] in edge2 or edge1[1] in edge2:
            pass
        # edges cannot share one vertex.
        # elif edge1[0].valence == 1 or edge1[1].valence == 1:
        #     if edge2[0].valence == 1 or edge2[1].valence == 1:
        #         pass
        # # If one edge has a hydrogene, there is no need to swap with an edge
        # # also containing a hydrogene.
        # elif edge2[0].valence == 1 or edge2[1].valence == 1:
        #     if edge1[0].valence == 1 or edge1[1].valence == 1:
        #         pass
        else :
            new_G1 = copy(G)
            new_G2 = copy(G)
            straight_switch(new_G1, pair[0], pair[1])
            crossed_switch(new_G2, pair[0], pair[1])
            all_nodes = tree.leaves + tree.upper_nodes
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

# def switch_all_edges(G):
#     all_pairs_of_edges = combinations(G.edges(), 2)
#     valid_graphs = [G]
#     # valid_hashes = {graph2hash(G)}
#     # stack = {graph2hash(G)}
#     for pair in all_pairs_of_edges:
#         edge1 = pair[0]
#         edge2 = pair[1]
#         if edge1[0] in edge2 or edge1[1] in edge2:
#             pass
#         else :
#             new_G1 = copy(G)
#             new_G2 = copy(G)
#             straight_switch(new_G1, pair[0], pair[1])
#             crossed_switch(new_G2, pair[0], pair[1])
#             all_nodes = leaves + upper_nodes
#             new_G1_is_valid = True
#             new_G2_is_valid = True
#             for node in all_nodes:
#                 if not node.is_connected(new_G1):
#                     new_G1_is_valid = False
#                 if not node.is_connected(new_G2):
#                     new_G2_is_valid = False
#             if new_G1_is_valid:
#                 valid_graphs.append(new_G1)
#             if new_G2_is_valid:
#                 valid_graphs.append(new_G2)
#     return valid_graphs

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



# T = [[[4],[3],[2,1]],[[3,2],[3,1]],[[4,1]]]
# tree = build_fragmentation_tree(T)
# g = generate_first_graph(tree)
# sol = generate_all_graphs(tree)

## Real fragmentation trees
f = formula2degrees
x = formula2degrees('C9H16')
y = formula2degrees('C8H9N')
# T = [[[[x ,[4,2]], y], [4,2]], [2,1,1]]
T = [[[[x,[4,2]],[y]],[[[4,2]]]], [[[[2,1,1]]]]]
tree = build_fragmentation_tree(T)
sol = generate_all_graphs(tree)
