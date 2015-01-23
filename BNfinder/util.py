# Copyright (c) 2004-2009 Bartek Wilczynski and Norbert Dojer; All Rights Reserved.
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

def safe_exponent(sc):
    """Gives the safe (against overflows) exponent of sc (base 2).

    In case of overflow, give the approximate value.
    """
    try:
        prob=str(2**sc)
    except OverflowError:
        sc10=math.log(2,10)*sc
        prob=str(10**(sc10-int(sc10)))+"e+"+str(int(sc10))
    return prob

def write_indent(file,txt,ind,ln=False):
    """write a string to a file with indent
    """
    file.write(ind)
    file.write(txt)
    if ln:
        file.write("\n")

def writeln_indent(file,txt,ind):
    """write an indented string followed by EOL
    """
    write_indent(file,txt,ind,ln=True)

def pprint(file,obj,level=0,indent="  "):
    """Pretty print a dictionary
    """
    try:
        items=obj.items()
    except AttributeError:
        file.write(repr(obj))
    else:
        writeln_indent(file,"{",level*indent)
        for k,val in sorted(items):
            write_indent(file,repr(k)+" : ",(level+1)*indent)
            pprint(file,val,level+1)
            writeln_indent(file,",","")
        write_indent(file,"}",level*indent)

def bin_search(x,l,b=0,e=None):
    """Finds in a _sorted_ list l the value nearest to x""" 
    if e==None:
        e=len(l)-1
    while b<e-1:
        i=(b+e)/2
        if l[i]<x:
            b=i
        elif l[i]>x:
            e=i
        else: #l[i]==x
            return i,l[i]

    #b==e find the nearer
    if abs(l[b]-x)>abs(l[b+1]-x):
        return b+1,l[b+1]
    else:
        return b,l[b]
        
                               
def rand_exp(n_genes,n_exp,n_disc):
    """Gives random expresssion data.

    n_genes : number of genes
    n_exp : number of experiments
    n_disc : nuimber of discrete levels
    """
    import random
    res = []
    for i in range(n_exp):
        exp = []
        for j in range(n_genes*2):
            exp.append(random.randint(0,n_disc-1))
        res.append(exp)
    return res

def test(n,m,d):
    """Testing function. Old code...
    """
    e = rand_exp(n,m,d)
    return d.learn_1(MDL(),0,range(n),e)
