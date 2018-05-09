from itertools import combinations


class Node:
    """
    A Node contains all conformations of a degree sequence.
    These conformations are defined by several graphs.
    """
    def __init__(self, degree = None, degreeSequence=None):
        self.degreeSequence = degreeSequence
        self.degree = degree
        self.parent = None
        self.children = []
        self.upper = None
        self.lower = None
        self.x = None
        self.freeBonds = None
        self.edges = []

    def add_children(self, list_of_children):
        self.children = list_of_children
        for child in list_of_children:
            child.parent = self
        self.set_upper()
        self.set_lower()

    def set_upper(self):
        if not self.children:
            self.upper = self.degree
        else:
            children_upper = []
            for child in self.children:
                children_upper.append(child.upper)
            upper = sum(children_upper) - 2 * (len(self.children) - 1)
            self.upper = upper

    def set_lower(self):
        if not self.children:
            self.lower = self.degree
        else:
            max_lower_child = max(self.children, key=operator.attrgetter('lower'))
            sum_of_upper = sum(child.upper for child in self.children)
            lower = max_lower_child.lower - sum_of_upper + max_lower_child.upper
            if lower >= 0 :
                self.lower = lower
            else:
                self.lower = lower % 2

    def is_valid(self):
        if self.lower == 0 & self.upper >=0 and self.upper % 2 == 0:
            return True
        else:
            return False

    # def set_x(self):
    #     if not self.parent:
    #         self.x = 0
    #     max_lower_child = max(self.children, key=operator.attrgetter('lower'))
    #     y = 2 * len(self.children) - 2 - sum(child.upper for child in self.children) - max_lower_child.upper
    #     x_max = max(max_lower_child.lower, y)
    #     if x_max <= 0:
    #         if x_max % 2 == 0: # If x_max is even : x_max = 2
    #             x_max = 2
    #         if x_max % 2 == 1: # If x_max is odd : x_max = 1
    #             x_max = 1
    #     for child in self.children:
    #         if child != max_lower_child:
    #             child.x = child.upper
    #         else:
    #             child.x = x_max

    def create_tree(self):
        #TODO
        if self.children:
            remaining_edges = len(self.children) - 1
            linked_atoms = []
            while remaining_edges > 0:
                sorted_children = sorted(self.children, key = lambda x: x.freeBonds, reverse=True)
                for i in range(len(sorted_children) - 1):
                    for j in range(i+1, len(sorted_children)):
                        if sorted_children[i] in linked_atoms and sorted_children[j] in linked_atoms:
                            pass
                        else:
                            self.create_virtual_edge(sorted_children[i], sorted_children[j])
                            remaining_edges -= 1
                            linked_atoms.append(sorted_children[i])
                            linked_atoms.append(sorted_children[j])

    def create_edge(self, atom1, atom2):
        if atom1 == atom2:
            raise ValueError("Loops aren't allowed.")
        self.edges.append(set([atom1,atom2]))
        # self.freeBonds -= 2 #TODO: need to do modifications on self.freeBonds
        atom1.freeBonds -= 1
        atom2.freeBonds -= 1

    def create_virtual_edge(self, node1, node2):
        if not node1.children:
            self.create_edge(node1, node2)
        # if not node1.children[0].children:
        #     max_atom1 = max(node1.children, key=operator.attrgetter('freeBonds'))
        #     max_atom2 = max(node2.children, key=operator.attrgetter('freeBonds'))
        #     self.create_edge(max_atom1, max_atom2)
        else:
            # sorted_children1 = sorted(node1.children, key = lambda x: x.freeBonds, reverse=True)
            max_child1 = max(node1.children, key=operator.attrgetter('freeBonds'))
            max_child2 = max(node2.children, key=operator.attrgetter('freeBonds'))
            self.create_virtual_edge(max_child1, max_child2)


    def remove_edge(self, atom1, atom2):
        try :
            self.edges.remove(set(atom1, atom2))
        except ValueError:
            raise ValueError("This edge doesn't exist.")
        self.freeBonds += 2
        atom1.freeBonds += 1
        atom2.freeBonds += 1

# def set_x(node):
#     if node.children:
#         max_lower_child = max(node.children, key=operator.attrgetter('lower'))
    def function_x(self):
        if self.parent == None:
            if self.is_valid():
                self.x = self.lower
                if self.children:
                    for child in self.children:
                        child.function_x()
            else :
                raise ValueError("Invalid fragmentation tree.")
        else:
            max_lower_node = max(self.parent.children, key=operator.attrgetter('lower'))
            if self == max_lower_node:
                y = 2 * len(self.parent.children) - 2 - sum(child.upper for child in self.parent.children) - self.upper
                x = max(self.lower, y)
                if x <= 0:
                    if x % 2 == 0: # If x is even : x = 2
                        x = 2
                    if x % 2 == 1: # If x is odd : x = 1
                        x = 1
                self.x = x
                for child in self.children:
                    child.function_x()

    def get_number_of_remaining_bonds(self):
        """
        This function compute the number of bonds needed for a node to be filled.
        """
        if not self.children:
            return self.x
        else:
            sum_of_x = 0
            for child in self.children:
                sum_of_x += child.x
            present_x = sum_of_x - 2 * (len(self.children) - 1)
            number_of_remaining_bonds = present_x - self.x
            return number_of_remaining_bonds

    # def add_remaining_bonds(self):
    #     remaining_bonds = self.get_number_of_remaining_bonds()
    #     if remaining_bonds == 0:
    #         self.create_tree()

    #     while remaining_bonds != 0:
    #         for child in self.children:
    #             if child.x == child.lower:
    #                 continue
    #             if child.x > child.lower:
    #                 if child.x - child.lower >= 2:

    def build_global_tree(self):
        self.create_tree()




def set_all_x_to_upper(root):
    root.x = root.upper
    if root.children:
        for child in root.children:
            set_all_x_to_upper(child)


# def set_all_x(node):
#     if node.children:
#         node.set_x()
#         for child in node.children:
#             set_all_x(child)
#     else:
#         pass


# def get_max_lower_child(node):
#     l = []
#     for child in node.children:
#         l.append((child.lower, child))
#     max_lower_child = max(l, key=operator.itemgetter(1))[1]
#     max_lower_value = max_lower_child.lower
#     # max_lower = max(l, key=operator.itemgetter(1))[0]
#     return max_lower_value


def set_atom(degree):
    atom  = Node(degree)
    atom.freeBonds = degree
    atom.set_upper()
    atom.set_lower()
    return atom



C1 = set_atom(4)
C2 = set_atom(4)
C3 = set_atom(4)
N1 = set_atom(3)
N2 = set_atom(3)
O1 = set_atom(2)
O2 = set_atom(2)
H1 = set_atom(1)
H2 = set_atom(1)

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()
n4 = Node()
n5 = Node()
n6 = Node()
n3.add_children([C1,H1])
n4.add_children([N1,O1])
n5.add_children([C2,O2])
n6.add_children([C3,N2,H2])
n1.add_children([n3,n4])
n2.add_children([n5,n6])
n0.add_children([n1,n2])

set_all_x_to_upper(n0)
n0.function_x()









