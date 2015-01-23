import random,collections

def combine_fold(folds,i):
    testing=folds[i]
    del folds[i]
    training=reduce(lambda x,y:x+y,folds)
    folds.insert(i,testing)
    return training,testing

def cv_folds(elems,k=10,classes=None,state=None):
    if state!=None: #save (or set) the state of random generator
        random.setstate(state)
    else:
       state=random.getstate()
        
    if classes==None: # if we have a class variable, split list into lists of elements with the same class variable
        lists=[elems]
    else:
        d=collections.defaultdict(list)
        for el,cl in zip(elems,classes):
            d[cl].append(el)
        lists=d.values()

    #shuffle the elements
    shuffled=[]
    for l in lists:
        random.shuffle(l)
        shuffled.extend(l)

    #create folds
    folds=[[] for i in range(k)]
    for i,e in enumerate(shuffled):
        folds[i%k].append(e)    
        
    return state,(combine_fold(folds,i) for i in range(k))
        
