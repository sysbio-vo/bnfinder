import rpy2.robjects as robjects

robjects.r("library(ROCR)")
robjects.r('auc<-function(p,s){performance(p,"auc",fpr.stop=s)}')

def device(filename):
    robjects.r.pdf(file=filename,width=5,height=5)

def close():
    robjects.r("dev.off()")

def aucROC(predictions,labels,part=1.):
    try:
        l= robjects.IntVector(labels)
    except:
        print "Problem with data",labels
    else:
        try:
            p= robjects.FloatVector(predictions)
        except:
            print "Problem with data",predictions
        else:
            pred = robjects.r.prediction(p,l)
            
            return robjects.r.auc(pred,part).do_slot("y.values")[0][0]

def simpleROC(predictions,labels,title,y="tpr",x="fpr",diag=True):
    try:
        l= robjects.IntVector(labels)
    except:
        print "Problem with data",labels
    else:
        try:
            p= robjects.FloatVector(predictions)
        except:
            print "Problem with data",predictions
        else:
            pred = robjects.r.prediction(p,l)
            perf = robjects.r.performance(pred,y,x)
            auc  = robjects.r.performance(pred,"auc").do_slot("y.values")[0][0]
            print "AUC",auc
            title="%s (AUC=%f)"%(title,auc)
            try:
                robjects.r.plot(perf,main=title,colorize=True)
            except:
                print "problem with plotting"
            if diag:
                robjects.r.abline(a=0,b=1)
    
    
def cvROC(predictions,labels,title,y="tpr",x="fpr",diag=True):
    #remove folds with only positive or only negative labels
    bad_folds=[]
    for i in range(len(labels)):
        if (True in labels[i]) and (False in labels[i]):
            pass
        else:
            bad_folds.append(i)
    print bad_folds
    for i in reversed(bad_folds): #remove in reversed order not to break the indexing
        del predictions[i]
        del labels[i]
        
    #Now do the ROC curve
    try:
        rlabs=map(robjects.IntVector,labels) #make IntVector out of each CV fold labels
        ll=robjects.r.list(*rlabs)
    except:
        print "Problem with data",labels
    else:
        try:
            rpreds=map(robjects.FloatVector,predictions)
            pp=robjects.r.list(*rpreds)
        except:
            print "Problem with data",predictions
        else:
            
            try:
                pred = robjects.r.prediction(pp,ll)
                perf = robjects.r.performance(pred,y,x)
            except:
                print "problem with ROCR",predictions,labels
            else:
                try:
                    robjects.r.plot(perf,main=title,col="grey82")
                    robjects.r.plot(perf,lwd=3,avg="vertical",add=True)
                    if diag:
                        robjects.r.abline(a=0,b=1)
                except:
                    print "problem with plotting"


def top_n(predictions,labels,title):
    try:
        l= robjects.IntVector(labels)
    except:
        print "Problem with data",labels
    else:
        try:
            p= robjects.FloatVector(predictions)
        except:
            print "Problem with data",predictions
        else:
            pred = robjects.r.prediction(p,l)
            perf = robjects.r.performance(pred,"prec","rec")
            xvs=perf.do_slot("x.values")[0]
            yvs=perf.do_slot("y.values")[0]
            print "AUC",auc
            title="%s (AUC=%f)"%(title,auc)
            try:
                robjects.r.plot(perf,main=title,colorize=True)
            except:
                print "problem with plotting"
            robjects.r.abline(a=0,b=1)
