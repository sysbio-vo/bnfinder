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

import math,fpconst
from score import score
from scipy.stats import chi2

class MIT(score):
    def __init__(self, *args,**kwds):
        score.__init__(self,*args,**kwds)
        self.alpha=kwds["chi_alpha"]
    
    def graph_score(self,number_of_potential_parents,gene_vertex,weights_of_parents,number_of_data_points): # g(Pa)
        
        # chisquare
        #print list_of_parents
            
#        sorted_parents = sorted(list_of_parents, key=lambda gene: gene.n_disc,reverse=True)
        #print "parent.size: ", len(list_of_parents)
        sum_of_l_values = 0.0
        
        for i,wi in enumerate(weights_of_parents):            
            l_i_sigma = (gene_vertex.base_weight() - 1)*(wi - 1)
            for w in weights_of_parents[:i]:
                l_i_sigma *= w
            chisquare_value = chi2.isf(1-self.alpha, l_i_sigma)
            #print "chisquare:", chisquare_value
            sum_of_l_values += chisquare_value

        #print [x.name for x in list_of_parents],gene_vertex.name
        #print "g_score: ", sum_of_l_values
        return sum_of_l_values    
            
    def lower_bound_for_data_score(self,selected_data):
        return 0.


    def data_score(self,selected_data):
        stats_all,stats_parents = selected_data.stats()
        stats_for_empty_parents, stats_for_empty_parents_fool = selected_data.subset().stats()
        #print "all",stats_all
        #print "parent",stats_parents
        #print "empty",stats_for_empty_parents

    
   
        score = 0.0
        #print "vertices",[v.name for v in parents_sequence]
        
        if selected_data.vertex.cpd in ['or','and']: # continuous
            raise Exception
        else:
            number_of_effective_observations = len(selected_data)
            #compute myEntropy
            H=0.0
            for k,v in stats_for_empty_parents.items():
                p=float(v)/number_of_effective_observations
                H-=p*math.log(p,2)
            #print H
            if len(selected_data.parents)==0: # CMI(X,[])==0
                score=H
            else:
                #compute CMI
                CMI=0.
                for obs_values,count in stats_all.items():
                    obs_parents=obs_values[:-1]
                    obs_vert=obs_values[-1]
                    p_xy=1.*count/number_of_effective_observations
                    p_x=1.*stats_parents[obs_parents]/number_of_effective_observations
                    p_y=1.*stats_for_empty_parents[(obs_vert,)]/number_of_effective_observations
                    CMI+=p_xy*math.log(p_xy/(p_x*p_y),2)
                    #print obs_parents,obs_vert,p_xy,p_x,p_y,p_xy*(math.log(p_xy/(p_x*p_y))),CMI
                score=H-CMI
            
        return 2*number_of_effective_observations*score*self.data_factor
        #print "d_score: ", score
        #return score


