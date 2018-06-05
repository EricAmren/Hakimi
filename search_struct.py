import mmap

def search_struct(mol_hash, filename):
    with open(filename, 'rb', 0) as file, \
    mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
        if s.find(mol_hash.encode()) != -1:
        # if s.find(b'cat') != -1:
            print('true')



creatine = ':Qa@jIaA_?JgHKkbCDNeFLN'
print(creatine)
example = ':Qi??C_?lHfGdEJaBKM`KLN'
python = 'Python anchietae'
search_struct(creatine, '../SageMath/graphs_creatine.txt')
# search_struct(example, '../SageMath/graphs_creatine.txt')
# search_struct(python, '../SageMath/exemple.txt')

