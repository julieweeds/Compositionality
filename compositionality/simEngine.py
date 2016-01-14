from operator import itemgetter

__author__ = 'juliewe'

import numpy as np, scipy.sparse as sparse,sys,math
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

#----
#get the path prefix / dependency path of a given feature
#self.getpathtype("amod:red") = amod
#self.getpathtype("_dobj>>amod:red") = _dobj>>amod
#----
def getpathtype(feature):
    #get the path of a given feature
    fields=feature.split(":")
    return fields[0]

class WordVector:

    def __init__(self,token):
        self.name=token
        self.features={}  #dictionary representation
        self.array=None #array representation
        self.lgth=-1
        self.wdth=-1
        self.total=0
        self.pathtotals={}

    def addfeature(self,f,sc):
        self.features[f]=sc
        self.total+=sc

    def updateweights(self,featdict):
        self.features=dict(featdict)
        self.total=0
        for f in self.features:
            self.total+=self.features[f]

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

    def computepathtotals(self):

        self.pathtotals={}
        for feature in self.features.keys():
            pathtype=getpathtype(feature)
            sofar=self.pathtotals.get(pathtype,0.0)
            self.pathtotals[pathtype]=sofar+float(self.features[feature])

    def reweight(self,weighting,feattots,typetots,grandtot=0,ppmithreshold=0):
         for feature in self.features.keys():
            freq=float(self.features[feature])  # C<w1,p,w2>
            try:
                total=float(self.pathtotals[getpathtype(feature)]) # C<w1,p,*>
            except:
                total=0.0001
                print "Warning: no path total for %s: %s"%(feature,getpathtype(feature))
            feattot=float(feattots[feature]) #C<*,p,w2>
            typetot=float(typetots[getpathtype(feature)]) #C<*,p,*>
            entrytotal=float(self.total) # C<w1,*,*>

            if "ttest" in weighting:
                expected = (total*feattot)/(typetot*typetot)
                obs=freq/typetot
                score= (obs-expected)/math.pow(expected,0.5)
                if score>ppmithreshold:
                    self.features[feature]=score
            else:

                try:
                    if "gof_ppmi" in weighting:

                        pmi=math.log10((freq*grandtot)/(feattot*entrytotal))
                    else:
                        pmi=math.log10((freq*typetot)/(feattot*total))
                except:
                    pmi=0
                shifted_pmi=pmi-ppmithreshold
                if shifted_pmi>0:
                    if "pnppmi" in weighting:

                        shifted_pmi=shifted_pmi * total/entrytotal

                    if "plmi" in weighting:
                        shifted_pmi=shifted_pmi * freq/typetot
                    self.features[feature]=shifted_pmi


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
        self.coltots_loaded={}
        self.rowtots_loaded={}
        for type in self.filenames.keys():
            self.vectors[type]={}
            self.coltots_loaded[type]=False
            self.rowtots_loaded[type]=False
        for type in self.filenames.keys():
            self.load(type)


        self.entrytotals={}
        self.feattots={}
        self.typetots={}
        #self.pathtots={}
        self.madematrix=False
        self.ppmithreshold=0


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

    def load_rowtotals(self,type):

        if not self.rowtots_loaded[type]:
            rowtotals=self.filenames[type]+".rtot"
            self.entrytotals={}
            print "Loading entry totals from: "+rowtotals
            with open(rowtotals) as instream:
                for line in instream:
                    line=line.rstrip()
                    fields=line.split("\t")
                    if self.include(fields[0]):
                        self.vectors[type][fields[0]].total=float(fields[1])

                        #sofar=self.entrytotals.get(fields[0],0)
                        #self.entrytotals[fields[0]]=sofar+float(fields[1])
            print "Loaded "+str(len(self.entrytotals.keys()))

            self.setrowsloaded(type)

    def setrowsloaded(self,type):
        for t in self.rowtots_loaded.keys():
            self.rowtots_loaded[t]=False
        self.rowtots_loaded[type]=True

    def setcolsloaded(self,type):
        for t in self.coltots_loaded.keys():
            self.coltots_loaded[t]=False
        self.coltots_loaded[type]=True


    def load_coltotals(self,type,cds=False):

        if not self.coltots_loaded[type]:
            coltotals=self.filenames[type]+".ctot"
            self.feattots={}
            print "Loading feature totals from: "+coltotals
            with open(coltotals) as instream:
                for line in instream:
                    line=line.rstrip()
                    fields=line.split("\t")
                    if fields[0] != "___FILTERED___":
                        feat=fields[0]
                        sofar=self.feattots.get(feat,0)
                        if cds:
                            self.feattots[feat]=pow(float(fields[1]),0.75)+sofar
                        else:
                            self.feattots[feat]=float(fields[1])+sofar
            print "Loaded "+str(len(self.feattots.keys()))
            self.setcolsloaded(type)

    #---
    #use the totals for each feature to compute grand totals for each feature type (e.g., "amod")
    #this is C<*,t,*> and is computed by S_f C<*,t,f>
    #----
    def compute_typetotals(self,type,cds=False):
        #compute totals for different paths over all entries (using column totals given in feattots)

        if not self.coltots_loaded[type]:
            self.load_coltotals(type,cds)
        print "Computing path totals C<*,t,*>"
        self.typetots={}
        for feature in self.feattots.keys():
            pathtype=getpathtype(feature)
            sofar=self.typetots.get(pathtype,0.0)
            self.typetots[pathtype]=sofar+float(self.feattots[feature])



    #--
    #compute path totals for each entry
    #i.e., C<w1,t,*>
    #this is needed in the standard PPMI calculation we use
    #----
    def compute_entrypathtotals(self,type):
        #compute totals for the different paths for each entry
        print "Computing path totals for each entry C<w1,t,*>"
        #self.pathtots={}
        for entry in self.vectors[type].keys():
            #totalvector={}
            self.vectors[type][entry].computepathtotals()


    def reweight(self,type,weighting=["ppmi"],ppmithreshold=0):

        #self.load_rowtotals(type)
        self.compute_typetotals(type,cds=("smooth_ppmi" in weighting))
        self.compute_entrypathtotals(type)
        grandtot=0.0
        if "pnppmi" in weighting:
            print "Computing pnppmi"
        elif "gof_ppmi" in weighting:
            print "Computing gof_ppmi"
            for type in self.typetots.keys():
                grandtot+=float(self.typetots[type])
        elif "plmi" in weighting:
            print "Computing localised PPMI"

        elif "ttest" in weighting:
            print "Computing ttest values"
        else:
            print "Computing ppmi"
        done =0
        todo=len(self.vectors[type].keys())

        for entry in self.vectors[type].keys():

            self.vectors[type][entry].reweight(weighting,self.feattots,self.typetots,grandtot=grandtot,ppmithreshold=ppmithreshold)
            done+=1
            if done%1000==0:
                percent=done*100.0/todo
                print "Completed "+str(done)+" vectors ("+str(percent)+"%)"



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

        done=0
        print self.vectors.keys()
        print pairlist
        results=[]
        for type in self.vectors.keys():
            for pair in pairlist:
                wordA=pair[0]+"/"+type  #this may need to be generalised as type is not always POS but it is for MEN exps
                wordB=pair[1]+"/"+type
                vectorA=self.vectors[type].get(wordA,None)
                vectorB=self.vectors[type].get(wordB,None)

                if vectorA==None or vectorB==None:
                    sim=-1
                else:
                    if simmetric in SimEngine.matrix_sims:
                        vectorA.makearray(self.fk_idx)
                        vectorB.makearray(self.fk_idx)
                    sim=vectorA.sim(vectorB,measure=simmetric)
                    vectorA.array=None
                    vectorB.array=None
                results.append(sim)
                if outstream==None:
                    print "%s(%s,%s) = %s"%(simmetric,wordA,wordB,str(sim))
                else:
                    outstream.write("%s\t%s\t%s\n"%(wordA,wordB,str(sim)))

                done+=1
                if (done%10000)==0:
                    percentage=done*100.0/todo
                    print "Completed %s calculations = %s percent"%(str(done),str(percentage))
        return results

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