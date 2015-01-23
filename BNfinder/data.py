#!/usr/bin/env python2.5
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

import graph,stats,math,fpconst
from itertools import chain
from BDE import BDE
from MDL import MDL
from MIT import MIT
import continuous
import util
from math import log
from random import shuffle


def eval_func(f_args):
    """Takes a tuple of a function and args, evaluates and returns result"""
    return f_args[0](*f_args[1:])

def learn_x(v):
    x,self,min_empty,min_optim,scr,verbose,n_min,limit,fpr_factor,max_tries,cores=v

    selected_data = self.select_1(x, sloops=scr.sloops)

    # dict comprehensions were introduced only in python 2.7
    import sys
    if sys.version_info >= (2,7,0):
        x_priors = {p : self.get_prior(x,p) for p in selected_data.parents}
    else:
        x_priors = {}
        for p in selected_data.parents:
            x_priors[p] = self.get_prior(x, p)

    selected_data.weight_parents(x_priors,fpr_factor,max_tries,scr)

    if min_empty:
        score_empty=scr.data_score(selected_data.subset([]))+scr.graph_score(len(selected_data.parents),x,[],len(selected_data))
        score_max=score_empty-math.log(min_empty,2)
    else:
        score_max=fpconst.PosInf
    if min_optim:
        score_delta=-math.log(min_optim,2)
    else:
        score_delta=fpconst.PosInf
    (par,sc),minlist = scr.learn_1(selected_data,verbose,n_min,limit,score_max,score_delta,cores)
    return x,par,minlist,sc



class experiment:
    """Class representing an experiment, possibly containing more than one time-point
    """
    def __init__(self,name,data,perturbed,points=[],static=False):
        """Constructor

        retains the argument data as self._data, while filling the self.data with double vectors
        (duplicates for static networks, consecutive pairs for dynamic)
        """
        self._data=data
        self.perturbed=perturbed
        self.data=[]
        self.name=name
        self.points=points
        if static:
            for i,d in enumerate(self._data):
                self.data.append(self._data[i]+d)
        else:
            # in format (data from time t - 1, data from time t)
            for i,d in enumerate(self._data[1:]):
                self.data.append(self._data[i]+d)

def unwrap_dataset_1_subset(arg, **kwarg):
    return dataset_1.subset(*arg, **kwarg)

class dataset_1:
    """Class representing data selected for learning parents of one vertex
    """
    def __init__(self,data,vertex,parents):
        self.data=data
        self.vertex=vertex
        self.parents=parents
        self.weights={}
        
    def weight_parents(self,priors,fpr_factor=None,max_tries=None,score=None):
        """Assigns weights to potential parents
        Required before learning
        """
        if self.parents:
            if fpr_factor: # calculate weights with fpr procedure
                n_parents = len(self.parents)
                n_data = len(self)
                no_tries = max(100, int(round( log(2+n_parents)*(max(priors.values())/fpr_factor) )))
                if max_tries:
                    no_tries = min(no_tries, max_tries)
                
                self_empty = self.subset()
                empty_score = score.data_score(self_empty) + score.graph_score(n_parents,self.vertex,[],n_data)
                sample_scores = [[score.lower_bound_for_data_score(self_empty)]*n_parents]
                
                # no_tries times shuffle vertex data and add scores to the distributions
                v_data = [d for d in self.data[-1]]  # create a copy before shuffling
                for i in range(no_tries):
                    shuffle(self.data[-1])
                    sample_scores.append([score.data_score(self.subset([par]))  for par in self.parents])
                distributions = [sorted(scores) for scores in zip(*sample_scores)]
                self.data[-1] = v_data
                
                for i,par in enumerate(self.parents):
                    # compute data_score thresholds
                    thr_ind = no_tries*fpr_factor/priors[par]
                    ind1 = min(no_tries, int(thr_ind))
                    ind2 = min(no_tries, ind1+1)
                    p = thr_ind-ind1
                    ds_thr = (1-p)*distributions[i][ind1] + p*distributions[i][ind2]

                    # compute required parents' graph_score
                    g_score = max(0, empty_score-ds_thr)
                    
                    # compute parent's weight
                    r = 2.001
                    r_score = score.graph_score(n_parents,self.vertex,[r],n_data)
                    while r_score < g_score:
                        r*=2
                        r_score = score.graph_score(n_parents,self.vertex,[r],n_data)
                    l = r/2
                    l_score = score.graph_score(n_parents,self.vertex,[l],n_data)
                    while r-l > 1e-9 and r_score-l_score > 1e-5:
                        w=(l+r)/2
                        w_score = score.graph_score(n_parents,self.vertex,[w],n_data)
                        if w_score < g_score:
                            l=w
                            l_score = w_score
                        else:
                            r=w
                            r_score = w_score
                    self.weights[par]=w
            
            else: # calculate weights based on priors only
                for par in self.parents:
                    self.weights[par]=par.base_weight()**priors[par]
            # sort parents and data according to weights
            order = [i for i,p in sorted(enumerate(self.parents), key = lambda (i,par) : self.weights[par])]
            self.parents = [self.parents[i] for i in order]
            self.data = [self.data[i] for i in order]+[self.data[-1]]
            
        
    def __len__(self):
        return len(self.data[-1])
        
    def signs(self):
        """Calculates signs of correlation coefficients between a child and its parents
        """
        res={}
        for i,p in enumerate(self.parents):
	    try:
		sign = stats.pearsonr(self.data[i],self.data[-1])[0]
	    except ValueError:
		print self.data[i],self.data[-1]
            if sign >=0:
                res[p]="+"
            else:
                res[p]="-"
        return res
        
    def subset(self,sub=[]):
        """Creates a new dataset_1 object with a subset of parents
        """
        ssub = set(sub)
        sub_ind=[i for i,p in enumerate(self.parents) if p in ssub]
        parents=[self.parents[i] for i in sub_ind]
        data = [self.data[i] for i in sub_ind]
        data.append(self.data[-1])
        d1 = dataset_1(data,self.vertex,parents)
        if self.weights:
            # dict comprehensions were introduced only in python 2.7
            import sys
            if sys.version_info >= (2,7,0):
                d1.weights = {par : self.weights[par] for par in parents}
            else:
                for par in parents:
                    d1.weights[par] = self.weights[par]
        return d1
        
#    def rm_sloops(self):
#        """Removes the child from parent set
#        Warning: modifies data
#        """
#        try:
#            i=self.parents.index(self.vertex)
#            del self.parents[i]
#            del self.data[i]
#        except ValueError:
#            pass
        
    def stats(self):
        """Counts frequencies of value vectors of
        parents and child (stats_all) and 
        parents only (stats_parents)
        """
        def stats_disc(data):
            stats_all = {}
            stats_par = {}
            for d in zip(*data):
                d_par = tuple(d[:-1])
                d_all = tuple(d)
                stats_par[d_par] = stats_par.get(d_par,0)+1
                stats_all[d_all] = stats_all.get(d_all,0)+1
            return stats_all,stats_par
  
        def strec(key,prob,d):
            first,first_disc=d[0]
            rest=d[1:]
            if rest==[]:
                stats_par[key] = stats_par.get(key,0)+prob
                if first_disc:
                    stats_all[key+(first,)] = stats_all.get(key+(first,),0)+prob
                else:
                    stats_all[key+(0,)] = stats_all.get(key+(0,),0)+prob*(1-first)
                    stats_all[key+(1,)] = stats_all.get(key+(1,),0)+prob*first
            else:
                if first_disc:
                    strec(key+(first,),prob,rest)
                else:
                    strec(key+(0,),prob*(1-first),rest)
                    strec(key+(1,),prob*first,rest)
                    
        d_disc=map(lambda p: p.n_disc,self.parents+[self.vertex])
        if 0 not in d_disc: #are all variables discrete?
            return stats_disc(self.data) # -yes
        # -no
        stats_all = {}
        stats_par = {}
        for d in zip(*self.data):
            strec((),1,zip(d,d_disc))
        return stats_all,stats_par
        
class gene:
    """Class corresponding to a single variable

    #TODO: rename to variable or smth... gene is inappropriate
    """
    def __init__(self,name,n_disc,index):
        self.name=name
        self.n_disc=n_disc
#        self.parents=[]
        self.values=[]
        self.index=index
        self.cpd=None
        
    def __str__(self):
        return self.name
    
    def disc_vals(self):
        """Returns discretized values of a variable

        works also for continuous variables
        """
        if self.values:
            try:
                return map(int,self.values)
            except ValueError:
                return self.values
        else:
            return [0,1] #This should take care of continuous variables... may need a change if we introduce more complex mixtures.

    def __index__(self):
        """returns the index of a variable.
        """
        return self.index

    def __add__(self,other):
        """Implements addition on indexes, so that we can write table[var+1] instad of table[var.index+1]
        
        """
        return other+self.index
    def __hash__(self):
        """returns the hash value for comparisons
        """
        return self.index

    def __eq__(self,other):
        if type(other)==int:
            return self.index==other
        else:
            return self.index==other.index    
            
    def base_weight(self):
        return max(1.5,self.n_disc)
        
class dataset:
    """
    A class representing a complete dataset.

    Implements input, output, data selection and also functions for learning networks given a scoring function
    """
    def __init__(self,name='no_name'):
        self.name=name
        self.n_series=0
        self.n_gen=0
        self.points=[]
        self.vertices= []
        self.vertice_names=[]
        self.static=False
        self.regulators=[]
        self.parents={}
        self.prior={}
        
    def fromFile(self,f,float_transform=continuous.BiNormMeans()):
        import warnings
        warnings.warn("The old data format is deprecated! Please use the new format.")
        #number of genes in the first line
        self.n_gen=int(f.readline())
        # experiments description
        ln = f.readline().strip().split()
        self.n_series = int(ln[0])
        series = []
        for txt in ln[1:]:
            if txt.find("[")==-1:
                series.append((int(txt),[]))
            else: #perturbed genes
                series.append((int(txt[:txt.find("[")]),
                               [int(txt[1+txt.find("["):txt.find("]")])]))
        lines = []
        self.vertices=[]
        for i in range(self.n_gen):
            ln = f.readline().strip().split()
            n_disc = max(0,int(ln[1]))
            self.vertices.append(gene(ln[0],n_disc,i))
            if n_disc==0:
                if float_transform: #use float_transformn if specified
                    lines.append(float_transform.estTrans(map(float,ln[2:])))
                else:
                    lines.append(map(float,ln[2:]))
            else:
                lines.append(map(int,ln[2:]))
        self.points=[]
        i=0
        for s,p in series:
            d = []
            for k in range(s): 
                d.append(map(lambda x:x[i],lines))
                i+=1
            name="serie%d"%len(self.points)
            self.points.append(experiment(name,d,p,range(len(d))))

#        for v in self.vertices:
#            v.parents=map(lambda v: (float(max(self.cont_weight,v.n_disc)),v),self.vertices)
        for e in self.points:
            e.perturbed=map(lambda i: self.vertices[i],e.perturbed)
        return self

    def toNewFile(self,file):
        #perturbations
        for s in self.points:
            for g in s.perturbed:
                file.write("#perturbed\t%s\t%s\n"%(s.name,self.vertices[g]))
        #condition names\
        file.write("conditions")
        for s in self.points:
            for i,d in enumerate(s._data):
                file.write("\t%s:%d"%(s.name,i))
        file.write("\n")
        #expr_data
        for i,v in enumerate(self.vertices):
            file.write(str(v))
            for s in self.points:
                for j,d in enumerate(s._data):
                    if v.n_disc>0: #discrete
                        format="\t%d"
                    else: #continuous
                        format="\t%f"
                    file.write(format%d[i])
            file.write("\n")

    def fromNewFile(self,file,float_transform=continuous.BiNormMeans()):
        pert=[]
        disc={}
        cpd={'and':[],'or':[]}
        contin=[]
        self.default=None
        rec=[]
        regul=[]
        while not rec:
            ln=file.readline()
            rec=ln.strip().split()
        #################
        #THE PREAMBLE
        ################
        #print 'preamble'
        while rec[0][0]=="#": #preamble
            if rec[0]=="#perturbed": #experiment serie name followed by list of its perturbed vertices
                pert.append((rec[1],rec[2:]))
            elif rec[0]=="#regulators": #list of potential regulators of all vertices (except specified in #parents command) not specified in previous or present #regulators command
                reg=filter(lambda x:True not in map(lambda r:x in r,regul),rec[1:])
                for i,r in enumerate(reg):
                    if i!=reg.index(r):
                        reg[i]=None
                reg=filter(None,reg)
                regul.append(reg)
            elif rec[0]=="#parents": #vertice name followed by list of its potential regulators
                reg=rec[2:]
                for i in range(len(reg)):
                    if i!=reg.index(reg[i]):
                        reg[i]=None
                reg=filter(None,reg)
                self.parents[rec[1]]=reg
            elif rec[0]=="#default": #list of default possible values or 'FLOAT' for default continuous values
                if len(rec)==2 and rec[1]=="FLOAT":
                    self.default=[]
                else:
                    self.default=rec[1:]
            elif rec[0]=="#discrete": #vertice name followed by list of all its possible values
                disc[rec[1]]=rec[2:]
            elif rec[0]=="#continuous": #list of vertices with continuous values
                contin.extend(rec[1:])
            elif rec[0]=="#and": #list of vertices with cpd noisy-and
                cpd['and'].extend(rec[1:])
            elif rec[0]=="#or": #list of vertices with cpd noisy-or
                cpd['or'].extend(rec[1:])
            elif rec[0]=="#priorvert": #prior weight of all outgoing edges (except specified in #prioredge command) followed by list of related vertices
                val=float(rec[1])
                if val <= 0:
                    raise Exception("Non-positive prior weight in dataset")
                for name in rec[2:]:
                    self.prior[name]=val
            elif rec[0]=="#prioredge": #prior weight of edges preceded by edge target (common for all specified edges) and followed by list of edge sources
                target_name=rec[1]
                val=float(rec[2])
                if val <= 0:
                    raise Exception("Non-positive prior weight in dataset")
                for name in rec[3:]:
                    self.prior[(target_name,name)]=val
            else:
                pass #ignoring other stuff
            rec=[]
            while not rec:
                ln=file.readline()
                rec=ln.strip().split()
        self.regulators=regul
        
        ################
        #condition names
        ################
        #print 'conditions'
        if rec[0]!='conditions':
            self.name=rec[0]
        #dynamic or static data?
        if len(rec[1].split(":"))==1:
            self.static=True
        #else:
        #    self.static=False
        if self.static:
            conds=map(lambda x : [x],rec[1:])
            names=conds
        else:
            # todo:pbd serie not used any more?
            serie=rec[1].split(":")[0]
            # pbd: explain why conds is a singleton and not just a single standalone element
            cond=[rec[1].split(":")[1]]
            conds=[cond]
            names=[serie]
            for c in rec[2:]:
                r=c.split(":")
                if r[0]==serie:
                    cond.append(r[1])
                else:
                    serie=r[0]
                    names.append(serie)
                    cond=[r[1]]
                    conds.append(cond)
            # pbd: after that conds looks like [['1', '2', '3', ...]]
        
        ln=file.readline()
        ###################
        # EXPRESSION DATA #
        ###################
        #print 'expression'
        self.vertices=[]
        self.vertice_names=[]
        #prepare the lists for expression data
        n_cond=sum(map(len,conds))
        exp_data=map(lambda l: [ [] for el in l],conds)
        geneindex=0
        while ln: #expression data
            rec=ln.strip().split()
            if not rec:
                ln=file.readline()
                continue
            name=rec[0]
            vals=[]
            #are the values discrete or not?
            if name in contin:#did we say it's continuous?
                discrete=False
            elif name in disc.keys(): #did we say is's discrete?
                discrete=True
                vals=disc[name]
            elif self.default!=None:# did we specify a default
                if self.default==[]:
                    discrete=False
                else:
                    discrete=True
                    vals=self.default
            else: #We did not specify the type of values for this gene, try to guess it
                try:
                    i = int(rec[1])
                except ValueError: #it's not an integer
                    try:
                        i=float(rec[1])
                    except ValueError: #it's not a float -> it's an alphanumerical value
                        discrete=True
                    else:
                        discrete=False # it's a float
                else:
                    discrete=True #it's an int, let's map the strings to ints
            if discrete:
                line=rec[1:]
                if vals==[]:#we don't know what are the possible values
                    for v in line:
                        if v not in vals:
                            vals.append(v)
                    #let's sort the values
                    vals.sort()
            else: #not discrete
                line=map(float,rec[1:])
                float_params=[]
                if float_transform: #use float_transform if specified
                    float_params=float_transform.estimate(line)
                    line=float_transform.transform(line)
                
            self.vertices.append(gene(name,len(vals),geneindex))
            self.vertices[-1].values=vals
            if discrete==False:
                self.vertices[-1].floatParams=float_params
            else:
                self.vertices[-1].floatParams=None
            self.vertice_names.append(name)
            #append this line to the expression data
            for el,l in zip(line,chain(*exp_data)):
                if discrete:
                    l.append(vals.index(el))
                else:
                    l.append(el)
            ln=file.readline()
            geneindex+=1
        ################
        #POSTPROCESSING#
        ################
        #print 'postprocessing'
        #process regulatory constraints
        for v in self.vertices:
            #specify cpd
            if v.name in cpd['and']:
                v.cpd='and'
            elif v.name in cpd['or']:
                v.cpd='or'
        #process the lists of perturbed genes and make the list of experiments
        for i,(c,n) in enumerate(zip(conds,names)):
            pg = []
            for x,p in pert:
                if x==n:
                    pg.extend(map(lambda n: self.vertices[self.vertice_names.index(n)],p))# map gene names to indexes
            
            self.points.append(experiment(n,exp_data[i],pg,c,self.static))

        #set the number of genes
        self.n_gen=geneindex
        return self


    def select_1(self,vertex,parents=None, sloops=True):
        """
        given a vertice and a list of his parents,
        select parents values at time t
        and child value at time t+1,
        omitting the experiments where the child was perturbed.
        """
        if parents is None:
            parents = self.get_potential_parents(vertex, sloops)
        res = []
        for expment in self.points:
            if vertex not in expment.perturbed:
                for sample in expment.data:
                    res.append([sample[par] for par in parents]+[sample[vertex+self.n_gen]])
        return dataset_1(map(list,zip(*res)),vertex,parents)


    def get_potential_parents(self, node, sloops=True):
        """selects potential parents for node"""
        if self.parents.has_key(node.name): # potential parents of a node are explicitly specified
            parents=map(lambda n: self.vertices[self.vertice_names.index(n)],self.parents[node.name])
        elif self.regulators: # regulators are specified
            parents=[]
            for r in self.regulators:
                """ !!! """
                if node.name in r:
                    break
                else:
                    parents.extend(map(lambda n: self.vertices[self.vertice_names.index(n)],r))
        else: # all vertices are potential parents
            parents=[v for v in self.vertices]
        if not sloops:
            try:
                i=parents.index(node)
                del parents[i]
            except ValueError:
                pass
        return parents

 
    def get_prior(self,vertex,parent):
        return self.prior.get((vertex.name,parent.name),self.prior.get(parent.name,1))
 
 
    def learn(self,score,data_factor, prior=None,cores=False, subset=False, \
            verbose=None,n_min=1,limit=None,min_empty=None,min_optim=None,fpr=None,max_tries=None):
        
        scr=score
        distribs = []

        # calculate fpr_factor
        if fpr:
            normalizer = 0.
            n_edge = 0
            for vertex in self.vertices:
                for parent in self.get_potential_parents(vertex,scr.sloops):
                    normalizer += 1./self.get_prior(vertex,parent)
                    n_edge +=1
            try:
                fpr_factor = fpr*n_edge/normalizer
#                fpr_factor = fpr*self.n_gen/normalizer
#                fpr_factor = fpr/normalizer
            except ZeroDivisionError:
                raise Exception("Potential parents set is empty for each vertex")
        else:
            fpr_factor = None

        # Handle case of computing subsets of genes
        vertices = self.vertices
        if (subset) and (subset != "concat"):
            filename = subset
            with open(subset, "r") as subset_file:
                subset = subset_file.readlines()
            subset = ''.join(subset).split()
            vertices = []
            for x in self.vertices:
                if x.name in subset:
                    vertices.append(x)

        if cores:
            # Some preparatory work
            try:
                import multiprocessing
                import multiprocessing.pool
            except ImportError:
                print "Problem invoking multiprocessing module. Running on a single CPU"
                cores=False
            else:
                # Daemonic processes are not allowed to have children
                class NoDaemonProcess(multiprocessing.Process):
                # make 'daemon' attribute always return False
                    def _get_daemon(self):
                        return False
                    def _set_daemon(self, value):
                       pass
                    daemon = property(_get_daemon, _set_daemon)

                class MyPool(multiprocessing.pool.Pool):
                    Process = NoDaemonProcess                
                
                # We need to define how to distribute cores among the tasks
                # In case of limit<=3 it is efficient to use hybrid alg
                # There are two situations: number of cores >= number of vertices and opposite one.
                # So, we need to handle it separately to distritube cores equaly
                pool=MyPool(1)
                if (limit is not None) and (int(limit)<=3) and (cores>int(limit)):
                    if (cores>=len(vertices)):
                        pool=MyPool(len(vertices))
                        for counter in range(1, len(vertices)+1):
                            if (counter<=(cores%len(vertices))):
                                distribs.append(cores/len(vertices) + 1)
                            else:
                                distribs.append(cores/len(vertices))
                    else:
                        count = cores/int(limit)
                        if (cores%int(limit)>0):
                            count = count + 1
                            
                        for counter in range(1, count+1):
                            if (counter<=(cores%count)):
                                distribs.append(cores/count + 1)
                            else:
                                distribs.append(cores/count)

                        temp = []
                        for i in range (1, len(vertices)+1):
                            temp.append(distribs[i%count])
                        distribs = temp 
                        pool=MyPool(count)
                
                # Computational part
                if verbose:
                    print "Computing in parallel"

                from itertools import izip,repeat
                import timeit
                start = timeit.default_timer()

                # Need to check if we just want to concatenate previously obtained results
                if (subset==False) or (subset != "concat"):
                    # Hybrid alg
                    if (limit is not None) and (int(limit)<=3) and (cores>int(limit)):
                        result=pool.map(learn_x,[(x,self,min_empty,min_optim,scr,verbose,n_min,limit,fpr_factor,max_tries, y) for x, y in zip(vertices, distribs)])
                    # Simple alg    
                    else:
                        result = map(learn_x, [(x,self,min_empty,min_optim,scr,verbose,n_min,limit,fpr_factor,max_tries, cores) for x in vertices])

                pool.close()
                pool.join()

                stop = timeit.default_timer()
                s = str(stop - start)
                if verbose:
                    print "----------------------------------------"
                    print "Time, secs: ", s
                    # Save computational time to file
                    with open("time"+str(cores)+".txt", 'w') as file:
                        file.write(s)

                # Save computed subset of genes to file
                if (subset) and (subset != "concat"):
                    file_pi = open("genes_"+filename+".obj", 'w') 
                    import pickle
                    pickle.dump(result, file_pi)
                    exit()
                 
        # Computing without any parallelizing                   
        else:
            result=map(learn_x,[(x,self,min_empty,min_optim,scr,verbose,n_min,limit,fpr_factor,max_tries, cores) for x in self.vertices])

        # Concatenating previously computed subsets of genes
        if subset == "concat":
            import pickle 
            result = []
            names = []
            import os
            for root, dirs, files in os.walk(os.getcwd()):
                for name in files:
                    if name.startswith("genes"):
                        names.append(name)
            
            for name in names:
                file_pi = open(name)
                result.extend(pickle.load(file_pi))

        #print result
        self.vertices=[r[0] for r in result]
        #TODO do we need this line?^

        total_score = 0.0
        pars ={}
        subpars ={}

        for x,par,minlist,sc in result:
            pars[x]=par
            subpars[x]=minlist
            total_score+=sc
        #print pars
        #move parents from a dictionary to a list
        par_list = []
        for v in self.vertices:
            try:
                par_list.append(pars[v.__index__()])
            except KeyError:
                print "ERRR",par_list,v.name,id(v),v.__index__(),[map(lambda x:(x.name,id(x)),l) for l in pars.values()]
        pars=par_list

        g = graph.graph()
        g.fromParents(self.vertices ,pars)
        for v,par in zip(self.vertices,pars):
            g.vertice_labelling[v]=v.name
            v_signs=self.select_1(v,par).signs()
            for p in par:
                g.edge_labelling[p,v]=v_signs[p]
        return total_score,g,subpars
   
    def get_stats(self,g):
        stats={}
        for v in self.vertices:
            stats[v] = self.select_1(v,g.parents(v)).stats()
        return stats
        
    def write_txt(self,subpars,file_name):
        """Outputs suboptimal parents sets into a text file
        """

        f=open(file_name,"w")
        for v,minlist in subpars.items():
            f.write('\n%s'% v.name)
            for prob, pars in minlist:
                prob_s=util.safe_exponent(prob)
                f.write('\n %s '% prob_s)
                for p in pars:
                    f.write(' %s'% p.name)
        f.close()

    def cpd_andor(self,n_par,stats_all,stats_par,prod_in,eps=10**(-7)):
        n_all=n_par+1
        epsilon=eps*n_all*10
        v_all=0.0
        v_in=0.0
        for all,val in stats_all.items():
            v_all+=val
            if all[-1]==prod_in:
                v_in+=val
        if v_all-v_in<=eps:
            return [0.0]*n_all
        elif v_in<=eps:
            return [1.0]*n_all
        p=[1.0-eps]*n_all
        delta=epsilon+1.0
        while delta>epsilon:
            delta=0.0
            for i in range(n_all):
                p_old=p[i]
                delt=eps+1.0
                while abs(delt)>eps:
                    sumf=0.0
                    sumfprim=0.0
                    for par in stats_par.keys():
                        all=par+(prod_in,)
                        if all[i]==prod_in:
                            if all in stats_all.keys():
                                prod=1.0
                                for j in range(n_all):
                                    if all[j]==prod_in:
                                        prod*=p[j]
                                q=stats_all[all]/(1.0-prod)
                                sumf+=q-stats_par[par]
                                sumfprim+=(q/(1.0-prod))*(prod/p[i])
                            else:
                                sumf-=stats_par[par]
                    delt=sumf/sumfprim
                    p[i]-=delt
                    if p[i]>=1.0:
                        delt=1.0-eps-p[i]-delt
                        p[i]=1.0-eps
                delta+=abs(p[i]-p_old)
        return p

    def to_cpd(self,g,dirichlet=None):
        if dirichlet==None:
            dirichlet=1.0
        else:
            dirichlet=float(dirichlet)
        result={}
        for v in self.vertices:
            if True:#if v.name not in self.regulators[0]:
                n_disc=max(2,v.n_disc)
                stats_all,stats_par = self.select_1(v,g.parents(v)).stats()
                var_d={}
                result[str(v)]=var_d
                var_d["vals"]=v.values
                var_d["pars"]=map(str,g.parents(v))
                var_d["floatParams"]=str(v.floatParams)
                if v.cpd=='and':
                    var_d["type"]="and"
                elif v.cpd=='or':
                    var_d["type"]="or"
                var_cpds={}
                var_d["cpds"]=var_cpds
                if v.cpd:
                    p=self.cpd_andor(len(g.parents(v)),stats_all,stats_par,v.cpd=='or')
                    for (i,pv) in enumerate(p[:-1]):
                        var_cpds[str(g.parents(v)[i])]=pv
                    var_cpds[None]=p[-1]
                else:
                    stats_values={}
                    for par_val in stats_par.keys():
                        stats_values[par_val]=[]
                    for all_val in stats_all.keys():
                        stats_values[all_val[:-1]].append(all_val[-1])
                    for par_val in stats_par.keys():
                        d={}
                        var_cpds[par_val]=d
                        for val in stats_values[par_val]:
                            d[val]=(stats_all[par_val+(val,)]+dirichlet)/(stats_par[par_val]+dirichlet*n_disc)
                        d[None]=dirichlet/(stats_par[par_val]+dirichlet*n_disc)
                    var_cpds[None]= 1.0/n_disc
        return result

        
        
    def write_cpd(self,g,f,dirichlet=None):
        if dirichlet==None:
            dirichlet=1.0
        else:
            dirichlet=float(dirichlet)
        f.write('{\n')
        for v in self.vertices:
            n_disc=max(2,v.n_disc)
            stats_all,stats_par = self.select_1(v,g.parents(v)).stats()
            f.write('\''+str(v)+'\' : {\n')
            f.write('  \'vals\' : '+str(v.values)+' ,\n')
            f.write('  \'pars\' : '+str(map(str,g.parents(v)))+' ,\n')
            if v.cpd=='and':
                f.write('  \'type\' : \'and\' ,\n')
            elif v.cpd=='or':
                f.write('  \'type\' : \'or\' ,\n')
            f.write('  \'cpds\' : {\n')
            if v.cpd:
                p=self.cpd_andor(len(g.parents(v)),stats_all,stats_par,v.cpd=='or')
                for (i,pv) in enumerate(p[:-1]):
                    f.write('     '+str(g.parents(v)[i])+' : '+str(pv)+'\n')
                f.write('     '+str(v)+' : '+str(p[-1])+' } } ,\n')
            else:
                stats_values={}
                for par_val in stats_par.keys():
                    stats_values[par_val]=[]
                for all_val in stats_all.keys():
                    stats_values[all_val[:-1]].append(all_val[-1])
                for par_val in stats_par.keys():
                    f.write('     '+str(par_val)+' : { ')
                    for val in stats_values[par_val]:
                        f.write(str(val)+' : '+str((stats_all[par_val+(val,)]+dirichlet)/(stats_par[par_val]+dirichlet*n_disc))+' , ')
                    f.write('None : '+str(dirichlet/(stats_par[par_val]+dirichlet*n_disc))+' } ,\n')
                f.write('     None : '+str(1.0/n_disc)+' } } ,\n')
        f.write('}\n')


    def write_bif(self,g,file_name,dirichlet=None,comments=[]):
        if dirichlet==None:
            dirichlet=1.0
        else:
            dirichlet=float(dirichlet)
        f=open(file_name,"w")
        f.write('\\\\ File generated by BNfinder\n')
        import time
        f.write('\\\\ '+time.strftime('%x %X')+'\n')
        for com in comments:
            f.write('\\\\ '+com+'\n')
        f.write('\\\\ Conditional probability distributions generated with total pseudocounts number %f\n' % dirichlet)
        f.write('\n')
        f.write('network \"%s\" {}\n' % self.name)
        f.write('\n')
        for v in self.vertices:
            n_disc=max(2,v.n_disc)
            stats_all,stats_par = self.select_1(v,g.parents(v)).stats()
            f.write('variable \"%s\" {\n' % str(v))
            f.write('   type discrete[%d] { \"%s\" }\n'% (n_disc,'\" \"'.join(map(str,v.disc_vals()))))
            f.write('    }\n')
            if g.parents(v):
                f.write('probability ( \"%s\" | \"%s\" ) {\n'% (str(v),'\" \"'.join(map(str,g.parents(v)))))
                f.write('     default %s ;\n' % ' '.join([str(1.0/n_disc)]*n_disc))
            else:
                f.write('probability ( \"%s\" ) {\n'% str(v))
            stats_values={}
            for par_val in stats_par.keys():
                stats_values[par_val]=[]
            for all_val in stats_all.keys():
                stats_values[all_val[:-1]].append(all_val[-1])
            for par_val in stats_par.keys():
                if par_val:
                    f.write('     ( \"%s\" ) ' % '\" \"'.join(map(str,par_val)))
                else:
                    f.write('     table ')
                for val in v.disc_vals():
                    if val in stats_values[par_val]:
                        f.write(str((stats_all[par_val+(val,)]+dirichlet)/(stats_par[par_val]+dirichlet*n_disc))+' ')
                    else:
                        f.write(str(dirichlet/(stats_par[par_val]+dirichlet*n_disc))+' ')
                f.write(';\n')
            f.write('    }\n\n')
        f.close()

    def point2dict(self,i):
        """returns a dictionary representing one (i-th) data-point
        """
	#
	return dict([(x,self.points[i].data[0][self.get_vertice_by_name(x)]) for x in self.vertice_names])
        #return dict(zip(self.vertice_names+map(lambda x: x+"'",self.vertice_names),self.points[i].data[0]))
	
    def get_vertice_by_name(self,name):
        """Returns vertex of self corresponding to the given name
        """
	for v in self.vertices:
	    if v.name==name:
		return v
        #return self.vertices[self.vertice_names.index(name)]
    
    def get_probs(self,d,cpd,pars,transformed=False):
        """returns list of posterior probabilities of all parents being in all possible states.

        for discrete variables gives singular distributions, important for continuous vars.
	if transformef==True then we assume that the dataset is already transformed into probabilities, otherwise apply continuous transformation
        """
        res=[]
        import collections
        for p in pars:
            v=self.get_vertice_by_name(p)
            if str(d[p]) in v.values: #discrete?
                dd=collections.defaultdict(float)
                dd[d[p]]=1.0
                res.append(dd) #singular distribution for discrete vars
            else: #continuous vars need care
		if transformed:
		    x=d[p]
		else:
		    model=continuous.BiNormMeans()
		    model.fromParams(eval(cpd[p]["floatParams"]))
		    x=model.transform([d[p]])[0] # x -> probability of parent being 1
                res.append({1:x,0:1-x})
#        print res
#        a=b
        return res
                
            
    def classify(self, cpd,f,ML=False,prob=None,transformed=False):
        """
        Output the classification based on the regulators' values from dataset and probabilities from a cpd.

        if ML==True:
           output only the most likely value
        else: (Default)
           output posterior probabilities
        if prob:
           output only the probability of a single class
        """
        f.write("classes")
        for p in self.points:
            f.write("\t%s"%p.name[0])
        f.write("\n")
        dicts=map(self.point2dict,range(len(self.points)))
        for v in cpd.keys():
            if v not in self.vertice_names:
                f.write("%s\t"%v)
                c=cpd[v]
                for i,p in enumerate(self.points):
                    total_prob={}
                    for val in c["vals"]:
                        total_prob[int(val)]=0.0
                    probs=self.get_probs(dicts[i],cpd,c["pars"],transformed=transformed)
                    for possib,poster in c["cpds"].items(): #possible parent values
                        if possib!=None:
                            r=reduce(lambda x,y:x*y,[d[x] for (x,d) in zip(possib,probs)],1)
                            for val in total_prob.keys():
                                try:
                                    total_prob[val]+=r*poster[val]
                                except KeyError:
                                    pass
                                
                    ret_d=total_prob
		    if i==5:
			#print i,probs,total_prob,c["vals"],dicts[i][c["pars"][0]],c["pars"],self.points[i].data[0][self.get_vertice_by_name(c["pars"][0])]
			#print self.get_vertice_by_name(c["pars"][0]).index
			#print dicts[i].items()
			#print [(self.get_vertice_by_name(x).index,self.points[i].data[0][self.get_vertice_by_name(x)]) for x in self.vertice_names]
			for v in self.vertices:
			    #print v.name[:10],v.index,
			    pass
                    #print ret_d
                    if ML: # Max likelihod
                        f.write("%s\t"%(str(sorted(zip(*reversed(zip(*(ret_d.items())))))[-1][-1]))) # pick the arg max from ret_d
                    elif prob!=None: # the probability of a given value
                        try:
                            f.write("%s\t"%(str(ret_d[prob])))
                        except KeyError:
                            f.write("%s\t"%"-0.0")
                    else: # output whole posterior
                        f.write("%s\t"%(str(ret_d)))
                    #par_vals=tuple(map(lambda x: dicts[i][x],c["pars"]))
                    #try:
                    #    ret_d=c["cpds"][par_vals]
                    #except KeyError:
                        #testing, to be removed
                        #print"*",
                        #model=continuous.BiNormMeans()
                        #print cpd[c["pars"][0]]["floatParams"]
                        #model.fromParams(eval(cpd[c["pars"][0]]["floatParams"]))
                        #x=model.transform(par_vals)[0] # x -> probability of parent being 1
                        #p=x*c["cpds"][1,][1]+(1-x)*c["cpds"][0,][1]
                        #f.write("%f\t"%(p))
                    #    f.write("%s\t"%"0.0")
                    #else:
#                         if ML: # Max likelihod
#                             f.write("%s\t"%(str(sorted(zip(*reversed(zip(*(ret_d.items())))))[-1][-1]))) # pick the arg max from ret_d
#                         elif prob!=None: # the probability of a given value
#                             try:
#                                 f.write("%s\t"%(str(ret_d[prob])))
#                             except KeyError:
#                                 f.write("%s\t"%"0.0")
#                         else: # output whole posterior
#                             f.write("%s\t"%(str(ret_d)))
                f.write("\n")
    def classify_list(self, cpd,ML=False,prob=None,transformed=False):
        """
        Output the classification based on the regulators' values from dataset and probabilities from a cpd.

        if ML==True:
           output only the most likely value
        else: (Default)
           output posterior probabilities
        if prob:
           output only the probability of a sigle class
        """
        result={}
        dicts=map(self.point2dict,range(len(self.points)))
        for v in cpd.keys():
            if v not in self.vertice_names:
                c=cpd[v]
                for i,p in enumerate(self.points):
                    total_prob={}
                    for val in c["vals"]:
                        total_prob[int(val)]=0.0
                    probs=self.get_probs(dicts[i],cpd,c["pars"],transformed=transformed)
                    for possib,poster in c["cpds"].items(): #possible parent values
                        if possib!=None:
                            r=reduce(lambda x,y:x*y,[d[x] for (x,d) in zip(possib,probs)],1)
                            for val in total_prob.keys():
                                try:
                                    total_prob[val]+=r*poster[val]
                                except KeyError:
                                    pass
                                
                    ret_d=total_prob
                    #print p.name[0],v
                    if ML: # Max likelihod
                        result[p.name[0],v]=sorted(zip(*reversed(zip(*(ret_d.items())))))[-1][-1] # pick the arg max from ret_d
                    elif prob!=None: # the probability of a given value
                        try:
                            result[p.name[0],v]=(ret_d[prob])
                        except KeyError:
                            result[p.name[0],v]=0.0
                    else: # output whole posterior
                        result[p.name[0],v]=ret_d
        return result
        
    
if __name__=="__main__":
    import sys
    try:
        new=open(sys.argv[2],"w")
        old=open(sys.argv[1])
        dataset().fromFile(old,None).toNewFile(new)
        new.close()
        old.close()
    except:
        print "As a quick and dirty solution for transforming old-style files into new style input files you may just run:\n python data.py old_file.txt new_file.txt"
