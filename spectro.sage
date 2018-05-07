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
        self.freeBonds = degree
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

    def set_x(self):
        max_lower_child = max(self.children, key=operator.attrgetter('lower'))
        y = 2 * len(self.children) - 2 - sum(child.upper for child in self.children) - max_lower_child.upper
        x_max = max(max_lower_child.lower, y)
        if x_max >= 0:   # TODO : what if x_max is lower than 0?
            if x_max %2 == 0:
                x_max += 2
            else:
                x_max += 1
        for child in self.children:
            if child != max_lower_child:
                child.x = child.upper
            else:
                child.x = x_max

    def create_tree(self):
        if self.children:
            remaining_edges = len(self.children) - 1
            while remaining_edges > 0:
                sorted_children = sorted(self.children, key = lambda x: x.freeBonds, reverse=True)
                print(sorted_children)
                self.create_virtual_edge(sorted_children[0], sorted_children[1])
                remaining_edges -= 1

    def create_edge(self, atom1, atom2):
        self.edges.append(set([atom1,atom2]))
        print("oij")
        self.freeBonds -= 2
        print("oij")
        atom1.freeBonds -= 1
        print("oij")
        atom2.freeBonds -= 1
        print("oij")

    def create_virtual_edge(self, node1, node2):
        # self must be parent of node1 and node2
        if not node1.parent == node2.parent:
            raise ValueError("Cannot create edge between node that aren't of same level.")
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




def set_all_x(node):
    if node.children:
        node.set_x()
        for child in node.children:
            set_all_x(child)
    else:
        pass


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
    atom.set_upper()
    atom.set_lower()
    return atom



C = set_atom(4)
N = set_atom(3)
O = set_atom(2)
H = set_atom(1)

n0 = Node()
n1 = Node()
n2 = Node()
n3 = Node()
n4 = Node()
n5 = Node()
n6 = Node()
n3.add_children([C,H])
n4.add_children([N,O])
n5.add_children([C,O])
n6.add_children([C,N,H])
n1.add_children([n3,n4])
n2.add_children([n5,n6])
n0.add_children([n1,n2])

set_all_x(n0)









