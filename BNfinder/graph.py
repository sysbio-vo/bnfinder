# Copyright (c) 2004-2007 Bartek Wilczynski and Norbert Dojer; All Rights Reserved.
#
# This software is distributable under the terms of the GNU
# General Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl-2.0.txt. Installing, importing or otherwise
# using this module constitutes acceptance of the terms of this License.
#
# Disclaimer
# 
# This software is provided "as-is".  There are no expressed or implied
# warranties of any kind, including, but not limited to, the warranties
# of merchantability and fittness for a given application.  In no event
# shall the authors be liable for any direct, indirect, incidental,
# special, exemplary or consequential damages (including, but not limited
# to, loss of use, data or profits, or business interruption) however
# caused and on any theory of liability, whether in contract, strict
# liability or tort (including negligence or otherwise) arising in any way
# out of the use of this software, even if advised of the possibility of
# such damage.
#




class stack:
    def __init__(self):
        self.stack = []
    def put(self,item):
        self.stack.append(item)
        return self
    def get(self):
        res = self.stack[-1]
        self.stack=self.stack[:-1]
        return res
    def empty(self):
        return len(self.stack)==0

class queue(stack):
    def get(self):
        res = self.stack[0]
        self.stack=self.stack[1:]
        return res


class graph:
    def __init__(self):
        self.vertices = []
        self.edges = {}
        self.edge_labelling = {}
        self.vertice_labelling = {}

    def add_edge(self,st,end,label=" "):
        try:
            self.edges[st].append(end)
        except KeyError:
            if st not in self.vertices:
                raise "Wrong starting vertice"
            self.edges[st] =[end]
            
        self.edge_labelling[st,end]=label
        
    def add_vert(self,vert,label=" "):
        self.vertices.append(vert)
        self.vertice_labelling[vert]=label
        self.edges[vert]=[]
        
    def __repr__(self):
        rep = "Graph: \n"
        for v in self.vertices:
            rep += "\t%s"%v
            try:
                rep += "(%s) => "%str(self.vertice_labelling[v])
            except KeyError:
                rep += " => "
            for v2 in self.edges[v]:
                rep+= "%s"%v2
                try:
                     rep+="(%s), "%str(self.edge_labelling[(v,v2)])
                except KeyError:
                    rep+=", "
            rep +="\n"
        #print self.edge_labelling, self.vertice_labelling
        return rep
            
    def ford_bellman(self,s):
        """
        Bellman - Ford algorithm
        """
        D = {}
        inf = len(self.vertices)
        for x in self.vertices:
            if x in self.edges[s]:
                D[x]=1
            else:
                D[x]=inf
        D[s]=0
        for i in range(len(self.vertices)-2):
            for v in self.vertices:
                if v != s:
                    for u in self.vertices:
                        if v in self.edges[u]:
                            D[v]=min(D[v],D[u]+1)
        return D

    def avg_dist(self):
        s = 0
        for v in self.vertices:
            s+=sum(self.ford_bellman(v).values())
        return s*1.0/(len(self.vertices)**2-len(self.vertices))

    def make_undirected(self):
        for v in self.vertices:
            for w in self.vertices:
                if (w in self.edges[v]) and (not v in self.edges[w]):
                    self.edges[w].append(v)
                    
    def c_v(self,v, method="plain"):
        """clustering coefficient of a vertice"""
        nom = 0
        denom=0
        if method=="plain":
            for u in self.edges[v]: #for all neighbours of v
                for w in self.edges[u]: # for all neighbours of a given neighbour 
                    if w in self.edges[v]: # if it is a neighbour of v
                        nom+=1
            denom=(len(self.edges[v])*(len(self.edges[v])-1)) #divided by possible triangles
        elif method=="FFA": #feed-forward "master"
            for u in self.edges[v]: #for all neighbours of v
                for w in self.edges[u]: # for all neighbours of a given neighbour 
                    if w in self.edges[v]: # if it is a neighbour of v
                        nom+=1
            denom=len(self.edges[v])*(len(self.edges[v])-1) #divided by pairs of outgoing edges
        elif method=="FFB": #feed-forward "helper"
            for u in self.parents(v): #for all parents of v
                for w in self.edges[u]: # for all neighbours of a given parent
                    if w in self.edges[v]: # if it is a neighbour of v
                        nom+=1
            denom=len(self.parents(v))*len(self.edges[v]) #divided  outgoing * incoming edges
        elif method=="FFC": #feed-forward "target"
            for u in self.parents(v): #for all parents of v
                for w in self.edges[u]: # for all neighbours of a given parent
                    if v in self.edges[w]: # if it is a neighbour of v
                        nom+=1
            denom=len(self.parents(v))*(len(self.parents(v))-1) #  incoming edges squared
        elif method=="FB": #feed-back loop
            for u in self.edges[v]: #for all neighbors of v
                for w in self.edges[u]: # for all neighbours of a given neighbor
                    if v in self.edges[w]: # if v is a neighbour of it
                        nom+=1
            denom=len(self.parents(v))*len(self.edges[v]) #divided  outgoing * incoming edges
        
        try:
            return nom*1.0/denom
        except ZeroDivisionError:
            return 0.0
    def clust_sig(self,meths=["FFA","FFB","FFC","FB"]):
        lst=[]
        avg={}
        for m in meths:
            avg[m]=0.0
        for v in self.vertices:
            vec=map(lambda m: self.c_v(v,m),meths)
            s=sum(vec)
            if s==0.0:
                norm_vec=vec
            else:
                norm_vec=map(lambda x: x/s,vec)
            lst.append(norm_vec)
            for m,v in zip(meths,norm_vec):
                avg[m]+=v
        s_all=sum(avg.values())
        for m in meths:
            avg[m]/=s_all
        return lst,avg
    
    def random_regular(self,n,k):
        """
        create a random graph with n vertices, each of degree k
        """
        import random
        self.__init__()
        self.vertices=range(n)
        for v in self.vertices:
            self.edges[v]=random.sample(self.vertices,k)

    def random_of_size(self,n,m,e_labels=["+","-"],v_labels=["AND","OR"]):
        """
        create a random graph with n vertices and m edges
        """
        import random
        self.__init__()
        self.vertices=range(n)
        for v in self.vertices:
            self.vertice_labelling[v]=random.choice(v_labels)
            self.edges[v]=[]
        edge_pool = []
        for u in self.vertices:
            for v in self.vertices:
                edge_pool.append((u,v))

        sample = random.sample(edge_pool,m)

        for u,v in sample:
            self.edges[u].append(v)
            self.edge_labelling[u,v]=random.choice(e_labels)
        
        
        
    def search(self, start, fun, struct=stack,avoid=None,forward=True,backward=False):
        if not type(start)==type([]): # did we get a list of starting points ?
            start = [start]
            
        if avoid==None:
            avoid=[] # We can't set it to [] as default, because static parameter binding sucks !
            
        st = struct()
        for v in start: #put all the starting points on the stack/queue/whatever
            st.put(v)
            avoid.append(v)
            
        while not st.empty():
            vert = st.get()
            fun(vert)
            #print "visiting", vert,self.edges[vert]
            if forward: #do we search along the edge orientation 
                for v in self.edges[vert]:
                    if v not in avoid:
                        st.put(v)
                        avoid.append(v)
            if backward: #do we search "backwards" the edge orientation
                for v in self.parents(vert):
                    if v not in avoid:
                        st.put(v)
                        avoid.append(v)
                    
    def subgraph(self,source,forward=True,backward=False):
        import sys,copy
        res = graph()
        f = lambda v: res.vertices.append(v)#sys.stdout.write("%s\n"%v)
        self.search(source,f,queue,None,forward,backward)
        for v in res.vertices:
            res.vertice_labelling[v]=self.vertice_labelling[v]
            res.edges[v]=copy.copy(self.edges[v])
            for w in res.edges[v]: 
                if not forward: #we need to check if there are no spurious edges
                    if w not in res.vertices:
                        res.edges[v].remove(w)
                        continue
                res.edge_labelling[(v,w)]=self.edge_labelling[(v,w)]
        #pruning
        g=res
        for v in g.vertices:
            for v2 in g.edges[v]:
                if v2 not in g.vertices:
                    g.edges[v].remove(v2)
        return g
        
    def verticeByLabel(self,label):
        for k,v in self.vertice_labelling.items():
            if v==label:
                return k

    def edgeByLabel(self,label):
        for k,v in self.edge_labelling.items():
            if v==label:
                return k

    def toDot(self):
        import pydot
        edges = []
        edg_keys = []
        for g in self.vertices:
            g_lab = str(g) + " - " + (self.vertice_labelling[g])
            for g2 in self.edges[g]:
                g2_lab = str(g2) + " - " + (self.vertice_labelling[g2])
                edges.append((g_lab,g2_lab))
                edg_keys.append((g,g2))  
        #print edges,edg_keys
        g = pydot.graph_from_edges(edges,directed=True)
        for e in g.edge_list:
            i = edges.index((e.get_source(),e.get_destination()))
            if self.edge_labelling[edg_keys[i]]=="+":
                e.set_label("+")
            elif self.edge_labelling[edg_keys[i]]=="-":
                e.set_label("-")
                e.set_arrowhead("tee")
            else: # no info
                e.set_label(self.edge_labelling[edg_keys[i]])
                e.set_arrowhead("none")
        return g

    def show(self):
        import tempfile

        f,fn = tempfile.mkstemp(".jpg")
        #f.close()
        self.toDot().write_jpg(fn)
        import Image
        im = Image.open(fn)
        im.show()
        import os
        os.remove(fn)
    
    def parents(self,vert):
        par =[]
        for v in self.vertices:
            if vert in self.edges[v]:
                par.append(v)
        return par
        
    def fromParents(self,verts,pars):
        self.__init__()
        self.vertices=[]
        for v in verts:
            self.add_vert(v)
            self.edges[v]=[]
            
        for v,par in zip(verts,pars):
            for p in par:
                self.add_edge(p,v)
        return self
    
    def fromOIMfile(self,f):
        ln = f.readline()
        par = []
        while ln and ln.find("-")==-1:
            par.append(map(int,ln.strip().split()))
            ln = f.readline()
        self.fromParents(range(len(par)),par)

    def num_edges(self):
        return sum(map(len,self.edges.values()))

    def p_value(self,g):
        n2 = len(self.vertices)**2
        s1 = self.num_edges()
        s2 = g.num_edges()
        x = 0
        x2 = 0
        for v in self.vertices:
            for w in self.edges[v]:
                if w in g.edges[v]:
                    x+=1
                    try:
                        if self.edge_labelling[(v,w)]==g.edge_labelling[(v,w)]:
                            x2+=1
                    except KeyError:
                        pass
        print s1,s2,x,n2
        print s1,s2,x2,n2
        sens = x*1.0/s2
        spec = 1.0*(n2-s1-s2+x)/(n2-s2)
        return hyp_geom(float(s1),float(s2),float(x),float(n2)),binom_p(x,s1,s2,n2),binom_p(x2,s1,s2,n2,0.33),sens,spec
    
    def weighted_edges(self, weights,fraction=0.0):
        #print "XX",fraction
        res=graph()
        res.vertices=self.vertices
        res.vertice_labelling=self.vertice_labelling
        for u in self.vertices:
            res.edges[u]=[]
        for u in self.vertices:
            if weights: # we have weights for all edges
                s=0.0
                w_dict={}
                if weights[u]:
                    w0,l=weights[u][0]
                    for w,l in weights[u]: #calculate the weights for edges
                        #print w,l
                        w2=2**(w-w0)
                        s+=w2
                        for node in l:
                            try:
                                w_dict[node]+=w2
                            except KeyError:
                                w_dict[node]=w2
                    for v,w in w_dict.items():
                        #print w/s,fraction,v,u
                        if w/s>fraction:
                            res.add_edge(v,u,str(w/s))
        return res
                        
                            
    def to_SIF(self,weights=None):
        result=""
        if weights:
            g1=self.weighted_edges(weights)
        else:
            g1=self
            
        for u in g1.vertices:
            for v in g1.edges[u]:
                label=g1.edge_labelling[(u,v)]
                result+="%s\t%s\t%s\n"%(g1.vertice_labelling[u],label,g1.vertice_labelling[v])
        #print "XX",result,g1.edges
        return result
    
    def from_SIF(self,file):
        self.__init__()
        ln = file.readline()
        while ln:
            rec=ln.strip().split()
            v = self.verticeByLabel(rec[0])
            if v==None: #it's a first time we see vertice v
                v=len(self.vertices)
                self.add_vert(v,rec[0])
            label=rec[1]
            for target in rec[2:]:
                t = self.verticeByLabel(target)
                if t==None: #it's a first time we see vertice t
                    t=len(self.vertices)
                    self.add_vert(t,target)
                self.add_edge(v,t,label)
            ln= file.readline()
        return self

    def labels_2_verts(self):
        new_v=[]
        new_vl = {}
        new_e={}
        new_el={}
        for v in self.vertices:
            l= self.vertice_labelling[v]
            new_v.append(l)
            new_vl[l]=(v,l)
            new_e[l]=[]
            for v2 in self.edges[v]:
                l2 = self.vertice_labelling[v2]
                new_e[l].append(l2)
                new_el[(l,l2)]=self.edge_labelling[(v,v2)]
        g = graph()
        g.vertices=new_v
        g.vertice_labelling=new_vl
        g.edges=new_e
        g.edge_labelling=new_el
        return g

    def compare(self,other):
        #find the intersection of vertices
        g1 = self.labels_2_verts()
        g2 = other.labels_2_verts()
        l = filter(lambda x: x in g2.vertices, g1.vertices)
        g3 = g1.subgraph(l,forward=False)
        g4 = g2.subgraph(l,forward=False)
        #print g3,g4,l
        return  g3.p_value(g4)

    def to_GML(self,f,geom=None,ref=None):
        f.write("graph\t[\n")
        for iv,v in enumerate(self.vertices):
            f.write("\tnode\t[\n")
            f.write("\t\tid\t%d\n"%iv)
            if geom:
                for ln in geom[self.vertice_labelling[v]]:
                    f.write(ln)
            f.write('\t\tlabel\t"%s"\n'%self.vertice_labelling[v])
            f.write('\t]\n')
        for iv,v in enumerate(self.vertices):
            for iu,u in enumerate(self.edges[v]):
                f.write("\tedge\t[\n")
                f.write("\t\ttarget\t%d\n"%self.vertices.index(u))
                f.write("\t\tsource\t%d\n"%iv)
                if not ref:
                    label = ""
                else:
                    if u in ref.edges[v]:
                        label=":V"
                    else:
                        label=":X"
                f.write('\t\tlabel\t"%s"\n'%(self.edge_labelling[(v,u)]+label))
                f.write('\t]\n')
        f.write(']\n')
        
                


def hyp_geom(s1,s2,x,n):
    """
    calculating hypergeometric distribution value using the simplified formula
    """
    def s(s1,s2,i,n):
        res = 1
        for j in range(1,i+1):
            res*=(s1-i+j)*(s2-i+j)/((j)*(n-i+j))
    #        print res
        for j in range(1,s2-i+1):
            res*=(n-s1-s2+i+j)/(n-i-s2+i+j)
    #        print res
        return res
    if s1<s2:
        return hyp_geom(s2,s1,x,n)
    sum = 0
    for i in range(x,s1+1):
        sum +=s(s1,s2,i,n)
    return sum


def fact(n):
    f = 1
    for x in range(n):
        f*=x+1
    return f

def binom(x,n,p=0.5):
     res = 0
     for k in range(x,n+1):
         p_k = p**k*(1-p)**(n-k)
         for i in range(min(k,n-k)):
             p_k*=1.0*(n-i)/(i+1)
             #print n-i,i+1,p_k,k
         res+=p_k
     return res
 
def binom_p(x,s1,s2,N,p=0.5):
    return binom(2*x+N-s1-s2,N,p)

def read_gml_geometry(f):
    def process_node(f):
        ln = f.readline()
        lns=[]
        while ln:
            rec=ln.strip().split()
            if rec[0]=="label":
                label = eval(rec[1])
            if rec[0]=="]":
                return label,lns
            if rec[0]=="graphics":
                while rec[0]!="]":
                    lns.append(ln)
                    ln = f.readline()
                    rec=ln.strip().split()
                lns.append(ln)
            ln = f.readline()
    nodes={}
    ln = f.readline()
    while ln:
        rec=ln.strip().split()
        if rec[0]=="node":
            #next node
            lab,lns=process_node(f)
            nodes[lab]=lns
        ln = f.readline()
    return nodes
