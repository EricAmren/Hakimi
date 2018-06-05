import cProfile
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
valence_dict['Cl'] = 7

def draw_mol(G):
    dic_colors={}
    dic_colors["black"]=[i for i in G if G.degree(i)==4]
    dic_colors["red"]=[i for i in G if G.degree(i)==2]
    dic_colors["white"]=[i for i in G if G.degree(i)==1]
    dic_colors["blue"]=[i for i in G if G.degree(i)==3]
    G.plot(vertex_colors=dic_colors).show()

def formula2degrees(formula):
    """
    From a chemical formula, return the corresponding degree sequence as a list.
    """
    splitted_formula = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
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

class Node:
    def __init__(self, atoms = {}):
        self.parent = None
        self.children = []
        self.atoms = list(atoms)
        self.depth = None

    def __gt__(self, other):
        try:
            return max(self.atoms) > max(other.atoms)
        except AttributeError :
            print("Those nodes are empty...")

    def __repr__(self):
        return str([atom.valence for atom in self.atoms])

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

    def rec_set_depth_of_nodes(self, depth):
        self.depth = depth
        depth += 1
        try:
            for child in self.children:
                child.rec_set_depth_of_nodes(depth)
        except TypeError:
            pass

class Atom:
    def __init__(self):
        self.vertex = None  # Index of this atom in general graph
        self.valence = None
        self.freeBonds = 0
        self.leaf = None

    def __gt__(self, other):
        return self.freeBonds > other.freeBonds

    def __repr__(self):
        return self.valence

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

def connect_node(node, G):
        free_children = node.children
        max_child = get_max_atom_of_nodes(node.children)[1]
        # max_child is a tuple containing the atom with the most of free bonds
        # and the node child containing this atom.
        connected_children = [max_child] ### Need to update this one
        free_children.remove(max_child)
        while free_children:
            max_connected_child = get_max_atom_of_nodes(connected_children)
            max_free_child = get_max_atom_of_nodes(free_children)
            if max_free_child[0].freeBonds == 0:
                raise ValueError("Not enough free atoms!")
            connect_2_atoms(G, max_connected_child[0], max_free_child[0])
            connected_children.append(max_free_child[1])
            free_children.remove(max_free_child[1])


def get_max_atom_of_nodes(node_list):
    max_atoms = []
    for node in node_list:
        max_atom = max(node.atoms)
        max_atoms.append((max_atom, node))
    return max(max_atoms)


# def connect_node2(node, G):
#     free_children = sorted(node.children, key=attrgetter('atoms'), reverse=True)
#     max_child = free_children.pop(0)
#     connected_children = [max_child]

#     while free_children:
#         max_connected_child = max(connected_children, key=attrgetter('atoms'))
#         if max(max_connected_child.atoms).freeBonds == 0:
#             raise ValueError("Not enough free atoms.")
#         max_free_child = free_children.pop(0)
#         try:
#             connect_2_atoms(G, max(max_connected_child.atoms), max(max_free_child.atoms))
#         except AssertionError:
#             print("Error in connect_node")
#             print(max(max_connected_child.atoms).valence, max(max_free_child.atoms).valence)
#             print(len(node.children[0].children[0].children[0].children[0].children))
#             exit(0)
#         connected_children.append(max_free_child)

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
            atom.leaf = leaf
            tree.atoms[index] = atom
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
        self.atoms = dict()



    def set_nodes(self):
        root = self.root
        branches = self.sequence
        self.rec_set_nodes(root, branches)
        self.upper_nodes = list(reversed(self.upper_nodes))

    def rec_set_nodes(self, node, branches):
        for branch in branches:
            if type(branch) == sage.rings.integer.Integer:
                new_atom = build_atom(branch)
                node.atoms.append(new_atom)
                new_atom.leaf = node
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

    def get_atom(self, index):
        return self.atoms[index]

    def set_depth_of_nodes(self):
        self.root.rec_set_depth_of_nodes(0)

def get_LCA(atom1, atom2):
    if atom1.leaf == atom2.leaf:
        return atom1.leaf
    else:
        node1 = atom1.leaf
        node2 = atom2.leaf
        return rec_get_LCA(node1, node2)

def rec_get_LCA(node1, node2):
    if node1 == node2:
        return node1
    else :
        if node1.parent == None:
            print(node1)
            print("is this root?")
            exit(0)
        if node1 == None:
            print("oijo!")
            exit(0)
        return rec_get_LCA(node1.parent, node2.parent)

def build_fragmentation_tree(sequence):
    tree = Fragmentation_Tree(sequence)
    tree.set_nodes()
    tree.set_depth_of_nodes()
    return tree


def relabel_graph(graph):
    keys = graph.vertices()
    values = graph.degree()
    dict1 = dict(zip(keys, values))
    for key, value in dict1.items():
        dict1[key] = str(key) + '-' + str(value)
    graph.relabel(dict1)

## Second step : Generate all graphs

def generate_all_graphs(tree):
    G = generate_first_graph(tree)
    valid_hashes = {graph2hash(G)}
    # valid_graphs = [G]
    stack = [G]
    while len(stack) > 0:
        new_G = stack.pop()
        rec_generate_all_graphs(tree, new_G, valid_hashes, stack)
        print(len(valid_hashes))
    return valid_hashes
    # return valid_graphs

def rec_generate_all_graphs(tree, G, valid_hashes, stack):
    all_pairs_of_edges = combinations(G.edges(), 2)
    with open("./graphs_creatine.txt", "a+") as f:
        f.write(graph2hash(G))
        f.write("\n")
        for pair in all_pairs_of_edges:
            edge1 = pair[0]
            edge2 = pair[1]
            if edge1 == edge2:
                pass
            atom11 = tree.get_atom(edge1[0])
            atom12 = tree.get_atom(edge1[1])
            atom21 = tree.get_atom(edge2[0])
            atom22 = tree.get_atom(edge2[1])
            # edges cannot share one vertex.
            if edge1[0] in edge2 or edge1[1] in edge2:
                pass
            # # If one edge has a hydrogene, there is no need to swap with an edge
            # # also containing a hydrogene.
            elif (atom11.valence == 1 or atom12.valence == 1) and (atom21.valence == 1 or atom22.valence == 1):
                pass
            else :
                last_common_ancestor = get_LCA_of_edges(tree, edge1, edge2)
                LCA1 = get_LCA(atom11, atom12)
                LCA2 = get_LCA(atom21, atom22)

                straight_switch(G, edge1, edge2)
                if LCA1.is_connected(G) and LCA2.is_connected(G):
                    if not graph2hash(G) in valid_hashes:
                        new_G = copy(G)
                        stack.append(new_G)
                        valid_hashes.add(graph2hash(new_G))
                        # valid_graphs.append(new_G)
                reverse_straight_switch(G, edge1, edge2)

                crossed_switch(G, edge1, edge2)
                if LCA1.is_connected(G) and LCA2.is_connected(G):
                    if not graph2hash(G) in valid_hashes:
                        new_G = copy(G)
                        stack.append(new_G)
                        valid_hashes.add(graph2hash(new_G))
                        f.write(graph2hash(new_G))
                        f.write("\n")
                        # valid_graphs.append(new_G)
                reverse_crossed_switch(G, edge1, edge2)

def get_LCA_of_edges(tree, edge1, edge2):
    atom11 = tree.get_atom(edge1[0])
    atom12 = tree.get_atom(edge1[1])
    atom21 = tree.get_atom(edge2[0])
    atom22 = tree.get_atom(edge2[1])
    LCA1 = get_LCA(atom11, atom12)
    LCA2 = get_LCA(atom21, atom22)
    while LCA1.depth != LCA2.depth:
        if LCA1.depth < LCA2.depth:
            LCA2 = LCA2.parent
        else:
            LCA1 = LCA1.parent
        assert LCA1 != None, "Last common ancestor is none!"
        assert LCA2 != None, "Last common ancestor is none!"
    while LCA1 != LCA2:
        LCA1 = LCA1.parent
        LCA2 = LCA2.parent
    return LCA1



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

def reverse_straight_switch(G, edge1, edge2):
    v1 = edge1[0]
    v2 = edge1[1]
    v3 = edge2[0]
    v4 = edge2[1]
    G.delete_edge((v1, v3, None))  # This may be dangerous is edge label are added...
    G.delete_edge((v2, v4, None))
    G.add_edge(v1,v2)
    G.add_edge(v3,v4)

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

def reverse_crossed_switch(G, edge1, edge2):
    v1 = edge1[0]
    v2 = edge1[1]
    v3 = edge2[0]
    v4 = edge2[1]
    G.delete_edge((v1, v4, None))
    G.delete_edge((v2, v3, None))
    G.add_edge(v1,v2)
    G.add_edge(v3,v4)

def graph2hash(graph):
    hash = graph.canonical_label().sparse6_string()
    return hash

def hash2graph(hash):
    assert isinstance(hash, str), "Wrong type : need a string"
    graph = Graph(hash, loops=False, multiedges=True, data_structure="sparse")
    return graph


def write_result_in_file(hashes, filename):
    with open (filename, w) as f:
        f.write(hashes)


# T = [[[4],[3],[2,1]],[[3,2],[3,1]],[[4,1]]]
# T = [[[4],[3],[2,1]],[[3,2],[3,1]],[4,1]]

# T = [[[4,1],[3,2]],[[4,2],[4,3,1]]]
# T = [[[4],[3,2]],[[4,2],[4,3]],[1],[1]]
## Real fragmentation trees
f = formula2degrees
# x = formula2degrees('C9H16')
# y = formula2degrees('C8H9N')
# T = [[[[x ,[4,2]], y], [4,2]], [2,1,1]]
# T = [[[[x,[4,2]],[y]],[[[4,2]]]], [[[[2,1]]]]]
# a = formula2degrees('C12H5ClN')
# T = [[[a, [1,7]],[[4,2]]], [[[4,1,1,1]]]]

# a = f('C9H17')
# b = f('C8H9N')
# T = [[[[a,[4,2]],[b]], [[[4,2]]]], [[[[2,1]]]]]

# a = f('C2H2O2')
# b = f('CH4N')
# c = f('CH2N2')
# T = [[a,b], [c]]

# Creatine real hash: ':Qa@jIaA_?JgHKkbCDNeFLN'
h = [[1]]
a = f('C2O2')
b = f('CN')
c = f('CN2')
T = [[a,b], [c], h,h,h,h,h,h,h,h,h]

# T = [[4,2],[1],[1],[1],[1]]

# tree1 = build_fragmentation_tree(T1)
# tree2 = build_fragmentation_tree(T2)
tree = build_fragmentation_tree(T)


# pr = cProfile.Profile()
# pr.enable()

sol = generate_all_graphs(tree)
# sol1 = generate_all_graphs(tree1)
# sol2 = generate_all_graphs(tree2)

# pr.disable()
# pr.print_stats(sort='time')
