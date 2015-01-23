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

import math,stats
from score import score

def log(x):
    return math.log(x,2)

class BDE(score):
    def __init__(self,prior=None,*args,**kwds):
        score.__init__(self,*args,**kwds)
        self.prior=prior
        
    def HP(self,v,par):
        if self.prior==None:
            return max(2,v.n_disc)
        else:
            hp=float(self.prior)
            for p in par:
                hp/=max(2,p.n_disc)
            return hp
    
    def H(self,v,par):
        if self.prior==None:
            return 1.0
        else:
            return self.HP(v,par)/max(2,v.n_disc)
    
    def graph_score(self,number_of_potential_parents,gene_vertex,weights_of_parents,number_of_data_points):
        return sum(map(log,weights_of_parents))*log(number_of_data_points+1)*min(1,gene_vertex.n_disc+0.1)
    
    def lower_bound_for_data_score(self,selected_data_empty):
        if self.prior!=None:
            return 0.0
        HP=self.HP(selected_data_empty.vertex,[])
        H=self.H(selected_data_empty.vertex,[])
        stats_all,stats_par = selected_data_empty.stats()
        if selected_data_empty.vertex.n_disc:
            s = 0
            for a,cv in stats_all.items():
                for i in range(0,cv):
                    s-=log(H+i)
                    s+=log(HP+i)
        else:
            return 0.0
#            gh = stats.gammln(H)
#            s = stats.gammln(HP+stats_par[()])-stats.gammln(HP)
#            for a,cv in stats_all.items():
#                s+=gh-stats.gammln(H+cv)
        return s*self.data_factor
        
    def data_score(self,selected_data):
        HP=self.HP(selected_data.vertex,selected_data.parents)
        H=self.H(selected_data.vertex,selected_data.parents)
        stats_all,stats_par = selected_data.stats()
        s = 0
        if 0 not in map(lambda p: p.n_disc,selected_data.parents+[selected_data.vertex]):
            par_counted = set()
            for a,cv in stats_all.items():
                numbers_of_parents = a[:-1]
                cp = stats_par[numbers_of_parents]
                if numbers_of_parents not in par_counted:
                    tmp =0
                    for i in range(0,cp):
                        tmp+=log(HP+i)
                    s+=tmp
                    par_counted.add(numbers_of_parents)
                for i in range(0,cv):
                    s-=log(H+i)
        else:
            gh = stats.gammln(H)
            ghp = stats.gammln(HP)
            for a,cv in stats_all.items():
                s+=gh-stats.gammln(H+cv)
            for a,cp in stats_par.items():
                s+=stats.gammln(HP+cp)-ghp

        return s*self.data_factor

