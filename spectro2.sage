def TODO():
    raise NotImplementedError("To be implemented")


class Fragmentation_Tree:
    def __init__(self):
        self.root = None
        self.leaves = []

    def find_all_conformations(self):
        self.verify_conditions()
        first_graph = self.generate_valid_graph()
        valid_conformations = {first_graph}
        stack = {first_graph}
        while stack != {}:


    def verify_conditions(self):
        TODO()

t = Fragmentation_Tree()





