from operator import itemgetter

__author__ = 'juliewe'

import numpy as np, scipy.sparse as sparse,sys
from composition import getorder,getpathtype

def isAny(token):
    return True

def profile(featdict,minorder=0,maxorder=10):

        paths={}
        totalweight=0
        thisorderweight=0
        for feat in featdict.keys():
            path = getpathtype(feat)
            order = getorder(feat)
            weight= featdict[feat]
            sofar=paths.get(path,0)
            if order>=minorder and order<=maxorder:
                paths[path]=sofar+weight
                thisorderweight+=weight
            totalweight+=weight

        #print "total weight of features",totalweight
        print "total weight of required order features",thisorderweight
        profile=sorted(paths.items(),key=itemgetter(1),reverse=True)

        print profile

class WordVector:

    def __init__(self,token):
        self.name=token
        self.features={}  #dictionary representation
        self.array=None #array representation
        self.lgth=-1
        self.wdth=-1
        self.total=0

    def addfeature(self,f,sc):
        self.features[f]=sc
        self.total+=sc

    def makearray(self,fk_idx):
        temparray=np.zeros(len(fk_idx))
        for feature in self.features.keys():
            temparray[fk_idx[feature]]=self.features[feature]

        self.array=sparse.csr_matrix(temparray)

    def width(self):
        if self.wdth==-1:
            self.wdth=len(self.features.keys())
        return self.wdth

    def length(self):
        if self.lgth==-1:
            self.lgth=pow(self.dotprod(self.array),0.5)
        return self.lgth

    def dotprod(self,anarray):
        return self.array.multiply(anarray).sum()



    def profile(self,minorder=0,maxorder=10):
        profile(self.features,minorder,maxorder)

    def cosine(self,avector):

        if (self.length()*avector.length())==0:
            sim =0
            print "Warning: zero length vectors"
        else:
            dotprod=self.dotprod(avector.array)
            sim = dotprod/(self.length()*avector.length())

        return sim

    def intersection(self,avector):
        intersect=0
        rtot=0
        ptot=0
        pset={}
        rset={}
        diff={}
        for feat in self.features.keys():
            if feat in avector.features.keys():
                intersect+=1
                pweight=avector.features[feat]
                rweight=self.features[feat]
                ptot+=pweight
                rtot+=rweight
                pset[feat]=pweight
                rset[feat]=rweight
            else:
                diff[feat]=self.features[feat]
        #profile(pset)
        print "Profiling intersection for ",self.name,avector.name
        profile(rset)
        profile(diff)
        return intersect,ptot,rtot

    def jaccard(self,avector):
        intersect= self.intersection(avector)[0]
        if intersect>0:
            den =self.width() + avector.width() - intersect
            return float(intersect)/float(den)
        else:
            return 0

    def precision(self,avector): #precision in predicting self by avector
        intersect=self.intersection(avector)[1]
        if intersect>0:
            #profile(avector.features)
            return float(intersect)/float(avector.total)
        else:
            return 0

    def recall(self,avector): #recall of self by avector

        intersect=self.intersection(avector)[2]
        if intersect>0:
            #profile(self.features)
            return float(intersect)/float(self.total)
        else:
            return 0

    def harmonicmean(self,avector):
        p=self.precision(avector)
        r=self.recall(avector)
        if p+r>0:
            return 2*p*r/(p+r)
        else:
            return 0

    def arithmeticmean(self,avector,beta=0.5):
        p=self.precision(avector)
        r=self.recall(avector)
        return beta*p+(1-beta)*r

    def sim(self,avector,measure):
        if measure =="cosine":
            return self.cosine(avector)
        elif measure =="jaccard":
            return self.jaccard(avector)
        elif measure =="recall":
            return self.recall(avector)
        elif measure=="precision":
            return self.precision(avector)
        elif measure=="harmonicmean":
            return self.harmonicmean(avector)
        elif measure=="arithmeticmean":
            return self.arithmeticmean(avector)
        else:
            print "Unknown similarity measure ",measure
            exit(-1)

    def reducesaliency(self,saliency,saliencyperpath=False):
        if saliency==0:
            return
        else:

            feats=sorted(self.features.items(),key=itemgetter(1),reverse=True)
            self.features={}
            donetypes={}
            all=0
            for tuple in feats:
                feature=tuple[0]
                pathtype=getpathtype(feature)
                done=donetypes.get(pathtype,0)
                if (saliencyperpath and done<saliency)or(not saliencyperpath and all<saliency):
                    self.features[feature]=tuple[1]
                    donetypes[pathtype]=done+1
                    all+=1


class SimEngine():
    #holds dictionary of vectors - manages conversion to sparse arrays and similarity calculations
    #requires a dictionary of filenames (anyRefKey->filename) for initialisation
    #optionally also an include_function - otherwise all vectors will be loaded from all files
    #allpairs similarities is cosine similarity within each file
    #pointwise similarities is for same key from different files e.g., to compare composed with observed vectors

    minorder=1 #minimum order to be used in similarity calculations
    maxorder=1 #maximum order to be used in simiarity calculations
    #paths_to_include=['_dobj','amod','nn','_nn','_nsubj']  #empty to include all
    paths_to_include=[]
    blacklist=[]
    matrix_sims=["cosine"]

    def __init__(self,filename_dict,include_function=isAny,pathdelim="\xc2\xbb",saliency=0,saliencyperpath=False):
        self.filenames=filename_dict
        self.vectors={} #dictionaries of vectors
        self.allfeatures={} #dictionary of all features observed for matrix generation
        self.fk_idx={} #feature --> dimension
        self.include_fn=include_function
        self.pathdelim=pathdelim
        self.saliency=saliency
        self.saliencyperpath=saliencyperpath
        for type in self.filenames.keys():
            self.vectors[type]={}
        for type in self.filenames.keys():
            self.load(type)

        self.madematrix=False



    def load(self,type):

        vectorfile=self.filenames[type]

        with open(vectorfile) as instream:
            for line in instream:
                line=line.rstrip()
                fields=line.split('\t')
                token=fields[0]
                featurelist=fields[1:]
                if self.include(token):
                    print "Loading ",token
                    try:
                        self.createvector(token,featurelist,type)
                    except:
                        print "Failed to load", featurelist
                        exit(-1)
                else:
                    #print "Ignoring : "+token
                    pass
        print "Loaded %s vectors from file %s with key: %s" %(str(len(self.vectors[type].keys())),vectorfile, type)


    def addfile(self,key, filename):
        self.filenames[key]=filename
        self.vectors[key]={}
        self.load(key)

        self.madematrix=False  #matrix must be remade

    def include(self,token):
        return self.include_fn(token)

    def includepath(self,feat):
        if len(SimEngine.paths_to_include)>0:
            return getpathtype(feat) in SimEngine.paths_to_include
        elif len(SimEngine.blacklist)>0:
            return getpathtype(feat) not in SimEngine.blacklist
        else:
            return True

    def createvector(self,token,featurelist,type):
        self.vectors[type][token]=WordVector(token)
        featurelist.reverse() #reverse list so can pop features and scores off

        while (len(featurelist)>0):
            f=featurelist.pop()
            sc=float(featurelist.pop())
            forder=getorder(f,delim=self.pathdelim)

            if forder>=SimEngine.minorder and forder<=SimEngine.maxorder: #filter features by path length
                if self.includepath(f):  #filter by path type
                    self.vectors[type][token].addfeature(f,sc)
                    self.allfeatures[f]=1

        self.vectors[type][token].reducesaliency(self.saliency,saliencyperpath=self.saliencyperpath)

    def makematrix(self):
        self.setup_matrix()

        for type in self.vectors.keys():
            for wordvector in self.vectors[type].values():
                wordvector.makearray(self.fk_idx)
        print "Completed matrix generation"


    def setup_matrix(self):
        print "Converting to matrix form"
        fkeys=self.allfeatures.keys()
        fkeys.sort()  #don't actually need to sort - but makes the indexing more predictable
        for i in range(len(fkeys)):
            self.fk_idx[fkeys[i]]=i

        del self.allfeatures #can now delete this dictionary if memory is an issue

        dim=len(self.fk_idx)
        self.madematrix=True
        print "Dimensionality is " + str(dim)


    def allpairs(self,outstream=None,simmetric="cosine"):

        if not self.madematrix and simmetric in SimEngine.matrix_sims:
            self.makematrix()

        todo=0
        for typeA in self.vectors.keys():
            todo+=pow(len(self.vectors[typeA].keys()),2)

        done=0

        for type in self.vectors.keys():
            for wordA in self.vectors[type].keys():
                for wordB in self.vectors[type].keys():
                    sim=self.vectors[type][wordA].sim(self.vectors[type][wordB],measure=simmetric)
                    if outstream==None:
                        print "%s(%s,%s) = %s [%s]"%(simmetric,wordA,wordB,str(sim),type)
                    else:
                        outstream.write("%s\t%s\t%s\t%s\n"%(type,wordA,wordB,str(sim)))
                    done+=1
                    if (done%10000)==0:
                        percentage=done*100.0/todo
                        print "Completed %s calculations = %s percent"%(str(done),str(percentage))


    def pointwise(self,outstream=None,simmetric="cosine"):

        if not self.madematrix and simmetric in SimEngine.matrix_sims:
            self.setup_matrix()

        todo=0
        for typeA in self.vectors.keys():
            todo+=len(self.vectors[typeA].keys())

        done=0
        for typeA in self.vectors.keys():
            for typeB in self.vectors.keys():
                if typeA !=typeB:
                    for wordA in self.vectors[typeA].keys():
                        vectorB=self.vectors[typeB].get(wordA,None)
                        if vectorB==None:
                            sim=0.0
                        else:
                            vectorA=self.vectors[typeA][wordA]
                            if simmetric in SimEngine.matrix_sims:
                                vectorA.makearray(self.fk_idx)
                                vectorB.makearray(self.fk_idx)
                            sim = vectorA.sim(vectorB,measure=simmetric)
                            vectorA.array = None
                            vectorB.array = None
                        if outstream==None:
                            print "%s(%s_[%s],%s_[%s]) = %s"%(simmetric,wordA,typeA,wordA,typeB,str(sim))
                        else:
                            outstream.write("%s\t%s\t%s\t%s\n"%(wordA,typeA,typeB,str(sim)))

                        done+=1
                        if (done%10000)==0:
                            percentage=done*100.0/todo
                            print "Completed %s calculations = %s percent"%(str(done),str(percentage))

    def selectedSims(self,pairlist,outstream=None,simmetric="cosine"):
        if not self.madematrix and simmetric in SimEngine.matrix_sims:
            self.setup_matrix()
        todo = len(pairlist)


        for type in self.vectors.keys():
            for pair in pairlist:
                wordA=pair(0)+"/"+type  #this may need to be generalised as type is not always POS but it is for MEN exps
                wordB=pair(1)+"/"+type
                vectorA=self.vectors[type].get(wordA,None)
                vectorB=self.vectors[type].get(wordB,None)

                if vectorA==None or vectorB==None:
                    sim=0
                else:
                    if simmetric in SimEngine.matrix_sims:
                        vectorA.makearray(self.fk_idx)  #TODO: finish

if __name__=="__main__":

    filename_dict={}
    key=1
    for filename in sys.argv[1:]:
        filename_dict[str(key)]=filename
        key+=1
    mySimEngine=SimEngine(filename_dict)

    outfilename="testout"
    outfilename=""
    measure="recall"

    if outfilename!="":
        outstream=open(outfilename,"w")
    else:
        outstream=None

    if len(filename_dict.keys())==1:
        mySimEngine.allpairs(outstream,simmetric=measure)
    else:
        mySimEngine.pointwise(outstream,simmetric=measure)


    if outstream!=None:
        outstream.close()