from sage.all import *

Liste_poids=[12.0107,1.00794,14.00674,15.9994]
Liste_atoms=["C","H","N","O"]
Liaisons=[4,1,3,2]

# def formula(N):

#     p = MixedIntegerLinearProgram(maximization=False, solver = "GLPK")
#     w = p.new_variable(integer=True, nonnegative=True)
#     if type(N)!=list:
#         p.add_constraint(p.sum([w[i]*Liste_poids[i] for i in range(len(Liste_poids))]) == N)
#     else:
#         n1=N[0]
#         n2=N[1]
#         p.add_constraint(p.sum([w[i]*Liste_poids[i] for i in range(len(Liste_poids))]) >= n1)
#         p.add_constraint(p.sum([w[i]*Liste_poids[i] for i in range(len(Liste_poids))]) <= n2)
#     p.add_constraint(p.sum([w[i]*Liaisons[i] for i in range(len(Liste_poids))])/2 >= p.sum([w[i] for i in range(len(Liste_poids))]) )
#     p.set_objective(None)
#     p.solve()
#     return p.get_values(w)





# def draw_formula(N):
#     dic=formula(N)
#     L=[i for i in dic if dic[i]!=0]
#     form=""
#     for i in L:
#         k=int(dic[i])
#         if k==1:
#             s=""
#         else:
#             s=str(k)
#         form+=Liste_atoms[i]+s
#     return form




def create_polyhedron_formula(n1,n2):
    ineqs=[]
    ineqs.append(tuple([-n1]+Liste_poids)) # weight larger than n1
    ineqs.append(tuple([n2]+[-i for i in Liste_poids])) # weight less than n2
    ineqs.append(tuple([1]+[(Liaisons[i]/2)-1 for i in range(len(Liste_poids))])) # The number of edges should be larger than the numebr of vertices -1
    for i in range(len(Liste_poids)):
        k=[0]*(len(Liste_poids)+1)
        k[i+1]=1
        ineqs.append(tuple(k))
    return Polyhedron(ieqs=ineqs)


def create_polyhedron_neutral_formula(n1,n2):
    ineqs=[]
    ineqs.append(tuple([-n1]+Liste_poids+[0])) # weight larger than n1
    ineqs.append(tuple([n2]+[-i for i in Liste_poids]+[0])) # weight less than n2
    ineqs.append(tuple([1]+[(Liaisons[i]/2)-1 for i in range(len(Liste_poids))]+[0])) # The number of edges should be larger than the numebr of vertices -1
    ineqs.append(tuple([0]+Liaisons+[-2])) # the degree sum should be even
    ineqs.append(tuple([0]+[-i for i in Liaisons]+[2])) # the degree sum should be even
    for i in range(len(Liste_poids)):
        k=[0]*(len(Liste_poids)+2)
        k[i+1]=1
        ineqs.append(tuple(k)) # all variables are positives
        #ineqs.append(tuple([0]+ Liaisons[0:i]+[-Liaisons[i]]+Liaisons[i+1:]+[0])) #to be modified
    return Polyhedron(ieqs=ineqs)


def create_polyhedron_neutral_formula2(n1,n2):
    ineqs=[]
    n=len(Liste_poids)
    ineqs.append(tuple([-n1]+Liste_poids+[0]*n+[0])) # weight larger than n1
    ineqs.append(tuple([n2]+[-i for i in Liste_poids]+[0]*n+[0])) # weight less than n2
    ineqs.append(tuple([1]+[(Liaisons[i]/2)-1 for i in range(n)]+[0]*n+[0])) # The number of edges should be larger than the numebr of vertices -1
    ineqs.append(tuple([0]+Liaisons+[0]*n+[-2])) # the degree sum should be even
    ineqs.append(tuple([0]+[-i for i in Liaisons]+[0]*n+[2])) # the degree sum should be even
    for i in range(n):
        k=[0]*(2*n+2)
        k[i+1]=1
        k2=[0]*(2*n+2)
        k2[n+i+1]=1
        ineqs.append(tuple(k)) # all variables are positives
        ineqs.append(tuple(k2)) # all variables are positives
        k3=[1]*(2*n+2)
        k3[n+i+1]=-1
        ineqs.append(tuple(k3))
        k4=[0]*(2*n+2)
        k4[n+i+1]=-1
        k4[i+1]=1
        ineqs.append(tuple(k4))
        k5=[0]*(2*n+2)
        k5[n+i+1]=1000
        k5[i+1]=-1
        ineqs.append(tuple(k5))
        k6=[0]+ Liaisons+ [0]*n +[0]
        k6[n+i+1]=-2*Liaisons[i]
        ineqs.append(tuple(k6)) #to be modified

    return Polyhedron(ieqs=ineqs)




def find_formulas(n1,n2,neutral=True):
    if  neutral:
        P=create_polyhedron_neutral_formula(n1,n2)
        S=P.integral_points()
        K=[]
        #KK=[i[:len(Liste_poids)] for i in S]
        for i in S:
            t=i[:len(Liste_poids)]
            f=[k for k in range(len(t)) if t[k]>0]
            m=max(Liaisons[k] for k in f)
            s=[k for k in f if Liaisons[k]==m][0]
            if t[s]>1 or sum(Liaisons[k]*t[k] for k in f[:s]+f[s+1:])>Liaisons[s]:
                K.append(t)

    else:
        P=create_polyhedron_formula(n1,n2)
        K=P.integral_points()


    return K


def formula_from_tuple(T):
    form=""
    for i in range(len(T)):
        if T[i]>0:
            if T[i]==1:
                s=""
            else:
                s=str(T[i])
            form+=Liste_atoms[i]+s
    return form

def html_from_tuple(T):
    form="$"
    for i in range(len(T)):
        if T[i]>0:
            if T[i]==1:
                form+=Liste_atoms[i]
            else:
                s=str(T[i])
                form+=Liste_atoms[i]+"_{"+s+"}"
    return html(form+"$")



def tuple_to_degree_seq(T):
    R=sorted(Liaisons,reverse=True)
    deg=[]
    for i in R:
        k=Liaisons.index(i)
        deg+=[Liaisons[k]]*T[k]
    return deg







# def build_neutral_graph(deg_seq):
#     n=len(deg_seq)
#     G=Graph(n,multiedges=True)

#     m=len([i for i in deg_seq if i>1])
#     part=range(m)
#     part.reverse()
#     G.add_edges([[i,i+1] for i in range(m-1)])
#     i=0
#     while i<n-m:
#         j=0
#         while j<len(part) and i<n-m:
#             G.add_edge(n-i-1,part[j])
#             i+=1
#             j+=1
#         part=[k for k in part if deg_seq[k] -G.degree(k)>0]
#     print G.degree_sequence()
#     for i in range(m):
#         t=deg_seq[i]-G.degree(i)
#         if t>0:
#             count=0
#             iter=0
#             while count<t:
#                 iter+=1
#                 if deg_seq[i+iter]-G.degree(i+iter)>0:
#                     G.add_edge(i,i+iter)
#                     count+=1
#     return G




def build_neutral_graph(deg_seq):
    # build first a not necessary connected graph
    n=len(deg_seq)
    G=Graph(n,multiedges=True)
    part=range(n)

    while len(part)>0:
        for i in range(min(deg_seq[part[0]]-G.degree(part[0]) , len(part)-1)):
            print("i")
            G.add_edge(part[0],part[1+ i])
        part=[i for i in part if deg_seq[i]>G.degree(i)]
        print("iii")
    # make the graph connected
    CC=G.connected_components()
    while len(CC)>1:
        #print len(CC), " connected components"
        found=False
        cou=0
        while (not found) and (cou<len(CC)):
            #print "Try find a non bridge"
            c=CC[cou]
            GG=G.subgraph(c)
            bridges=GG.bridges()
            if not (len(bridges)==GG.size()):
                e1=list(Set(GG.edges())-Set(bridges))[0]
                #print "non bridge found! in" + str(c)
                found=True
            cou+=1


        cou-=1
        for i in range(len(CC)):
            #print "Find another edge"
            if not i==cou:
                c2=CC[i]
                GG=G.subgraph(c2)
                if GG.size()>0:
                    e2=GG.edges()[0]
        #print e1,e2
        G.add_edge(e1[0],e2[0])
        G.add_edge(e1[1],e2[1])
        G.delete_edge(e1)
        G.delete_edge(e2)
        CC=G.connected_components()
    return G



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


def switch(G,v0,v1,v2,v3):
    """ v0v2 and v1v3 should be edges """
    G.add_edge(v0,v1)
    G.add_edge(v2,v3)
    G.delete_edge(v0,v2)
    G.delete_edge(v1,v3)




def fragmentation_tree(inter_list):
    max=0
    ind=0
    for i in range(len(inter_list)):
        if inter_list[i][0]>max:
            max =inter_list[i][0]
            ind=i
    Others=[]
    All=[]
    Max_formula=find_formulas(inter_list[ind][0],inter_list[ind][1],neutral=True)
    Max_formula=[tuple(i) for i in Max_formula]
    All.append(Max_formula)
    for i in range(len(inter_list)):
        if not i==ind:
            R=list(find_formulas(inter_list[i][0],inter_list[i][1],neutral=False))
            R=[tuple(i) for i in R]
            Others+=R
            All.append(R)
    print Max_formula+Others
    print All
    P=Poset(([tuple(i) for i in Max_formula+Others],lambda x,y:all(x[i]<=y[i] for i in range(len(x)) )))
    #P=Poset((Max_formula+Others,lambda x,y : all(x[i]<=y[i] for i in range(min(len(x),len(y))))))
    I=lambda z: [All.index(i) for i in All if z in i][0]

    Subelements=[]
    for z in Max_formula:
        if Set([I(x) for x in P.order_ideal([z])]).cardinality() ==len(inter_list):
            Subelements+=P.order_ideal([z])



    #[P.order_ideal([z]) for z in Max_formula if Set([I(x) for x in P.order_ideal([z])]).cardinality() ==len(inter_list) ]

    print Subelements

    return P.subposet(Subelements)#[tuple(i) for i in Max_formula+Others]
