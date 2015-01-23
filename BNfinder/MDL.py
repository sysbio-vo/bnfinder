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

import math
from score import score

class MDL(score):
    
    def graph_score(self,number_of_potential_parents,gene_vertex,weights_of_parents,number_of_data_points): # g(Pa)
        graph_structure=(1+len(weights_of_parents))*math.log(number_of_potential_parents+1,2) # = |Pa| * log n
        graph_cpd=0.5*(max(1.1,gene_vertex.n_disc)-1)*math.log(number_of_data_points+1,2) # = 1/2 * log N * (k_x - 1)
        if gene_vertex.cpd in ['or','and','preproc']:
            graph_cpd*=1+sum(weights_of_parents)/2.0
        else:
            for weight in weights_of_parents:
                graph_cpd*=weight
        return graph_structure+graph_cpd
    
    
    def lower_bound_for_data_score(self,selected_data):
        if selected_data.vertex.n_disc:
            return 0.0
        score=0.0
        for d in selected_data.data:
            try:
                score -= d[0]*math.log(d[0],2) + (1-d[0])*math.log(1-d[0],2)
            except:
                pass
        return score*self.data_factor

    def data_score(self,selected_data):
        stats_all,stats_parents = selected_data.stats()
        score = 0.0
        if selected_data.vertex.cpd in ['or','and']: # continuous
            prod_in=(selected_data.vertex.cpd=='or')
            prod_out=1-prod_in
            probabilities=self.cpd_andor(len(selected_data.parents),stats_all,stats_parents,prod_in)
            for a in stats_parents.keys():
                prob_a=probabilities[-1] # last element
                for i,av in enumerate(a):
                    if av==prod_in:
                        prob_a*=probabilities[i]
                try:
                    score-=stats_all[a+(prod_out,)]*math.log(prob_a,2)
                except:
                    pass
                try:
                    score-=stats_all[a+(prod_in,)]*math.log(1-prob_a,2)
                except:
                    pass
        else:
            for a,cv in stats_all.items():
                numbers_of_parents = a[:-1] # except the last one
                cp = stats_parents[numbers_of_parents]
                try:
                    sc = math.log(1.0*cp/cv,2) # - log(cv/cp)
                    score+=cv*sc
                except ZeroDivisionError:
                    pass
        return score*self.data_factor


