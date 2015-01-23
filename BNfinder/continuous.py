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

import stats,math,util

def linspace(minv,maxv,num):
    """
    """
    res=[]
    step=1.0*(maxv-minv)/(num-1)
    v=minv
    while v<=maxv:
        res.append(v)
        v+=step
    return res

class Transform:
    """general (vritual) class for continuous data transformations
    """
    def __init__(self,*args,**kwds):
        """Initialize general parameter of a transformation
        """
        pass
    def estimate(self,vector):
        """Estimate parameters of transformation dependent on input data

        estimates parameters and seves them internally.
        Returns serialized parameters for use with .fromParams method
        """
        pass
    
    def transform(self,vector):
        """Transform a vector of data using internal parameters
        """
        pass
    
    def estTrans(self,vector):
        """Estimate transformation parameters and transform data in one go.

        Only the transformed data are returned. 
        """
        self.estimate(vector)
        return self.transform(vector)
    
    def fromParams(self,params):
        """Fill in internal parameters using serialized data.
        """
        pass
    def getParams(self):
        """Return serialized parameters.
        """
        pass

class BiNormTransform(Transform):
    def __init__(self,t=1.0):
        self.t=t
        
    def transform(self,l):
        trans=lambda x: 1.0/(1.0+(self.p0/self.p1)*math.exp((self.mu1**2-self.mu0**2+2*(self.mu0-self.mu1)*x)/(2*self.t*self.sigma**2)))
        ret=[]
        for x in l:
            try:
                ret.append(trans(x))
            except OverflowError:
                #print "ERROR",x
                ret.append(0.0)
        return ret

    def fromParams(self,p):
        self.mu0,self.mu1,self.sigma,self.p0,self.p1=p
        
    def getParams(self):
        return self.mu0, self.mu1, self.sigma, self.p0, self.p1

    
class BiNormSimple(BiNormTransform):
    def estimate(self,l):
        mu=stats.mean(l)
        self.sigma=stats.stdev(l)/2
        self.mu0=mu-self.sigma
        self.mu1=mu+self.sigma
        self.p0=self.p1=0.5
    
        return self.mu0, self.mu1, self.sigma, self.p0, self.p1

class BiNormMeans(BiNormTransform):
    """Simple estimation of mixture of two gaussians

    Adapted from previous code in 
    """
    def estimate(self,l):
        def sumvar(l,m):
            return stats.sum(map(lambda x:(x-m)**2,l))
        med=stats.lmedianscore(l)
        oldmed=med+1
        while med!=oldmed:
            l0=filter(lambda x: x<med,l)
            l1=filter(lambda x: x>med,l)
            lm=filter(lambda x: x==med,l)
            mu0=stats.mean(l0*2+lm)
            mu1=stats.mean(l1*2+lm)
            (oldmed,med)=(med,(mu0+mu1)/2)
            #print mu0,mu1, oldmed,med#,l0,l1
        self.sigma=max(0.1**10, math.sqrt((sumvar(l0+lm,mu0)+sumvar(l1,mu1)+(mu1-med)**2)/(len(l))))
        self.p0=(len(l0)+len(lm)*0.5)/len(l)
        self.p1=(len(l1)+len(lm)*0.5)/len(l)
        self.mu0=mu0
        self.mu1=mu1
        return self.mu0, self.mu1, self.sigma, self.p0, self.p1


class GaussianMixture(Transform):
    """ Implements a mixture of N Gaussian variables. Estimation done with EM algorithm.

    default N==2, we assume a common variance"""
    def __init__(self,n=2):
        self.n=n
        

    def _Estep(self,mus,vector):
        #assign the nearest mu to each observation
        import collections
        obs=collections.defaultdict(list)
        for v in vector:
            m,val=util.bin_search(v,mus)
            obs[m].append(v)
        return obs
        
    def _Mstep(self,obs):
        #find best self.mus given self.assign
        mus=[]
        for i in range(self.n):
            mus.append(stats.mean(obs[i]))
        return mus
    
    def estimate(self,vector):
        sv=sorted(vector)
        minv=sv[0]
        maxv=sv[-1]
        self.mus=linspace(minv,maxv,self.n+1)
        self.mus=[(x+self.mus[i+1])/2 for i,x in list(enumerate(self.mus))[:-1]]
        self.obs=self._Estep(self.mus,sv)
        mus=self._Mstep(self.obs)
        while mus!=self.mus:
            self.mus=mus #update means
            self.obs=self._Estep(mus,sv) #update observation assignment
            mus=self._Mstep(self.obs) # calculate new mus
        #estimate variance
        self.var=0.0
        for i,mu in enumerate(self.mus):
            for x in self.obs[i]:
                self.var+=(x-mu)**2
        self.var/=len(vector)
        #estimate priors
        self.priors={}
        for x,l in self.obs.items():
            self.priors[x]=len(l)*1.0/len(vector)
        return self.mus,self.var,self.priors
            
    def transform(self,vector):
        """Transformation using multiple gaussian mixture. as of yet unfinished, requires code changes in data.py and score.py
        #TODO: adapt data.py and score.py to allow more general gaussian mixtures
        """

        result=[]
        for v in vector:
            pass
        #trans=lambda x: 1.0/(1.0+(self.p0/self.p1)*math.exp((self.mu1**2-self.mu0**2+2*(self.mu0-self.mu1)*x)/(2*self.t*self.sigma**2)))
        #return map(trans,l)
