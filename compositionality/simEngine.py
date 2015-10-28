__author__ = 'juliewe'

import numpy as np, scipy.sparse as sparse,sys

def isAny(token):
    return True

class WordVector:

    def __init__(self,token):
        self.name=token
        self.features={}  #dictionary representation
        self.array=None #array representation
        self.lgth=-1
        self.wdth=-1

    def addfeature(self,f,sc):
        self.features[f]=sc

    def makearray(self,fk_idx):
        temparray=np.zeros(len(fk_idx))
        for feature in self.features.keys():
            temparray[fk_idx[feature]]=self.features[feature]

        self.array=sparse.csr_matrix(temparray)

    def width(self):
        if self.wdth==-1:
            self.width=len(self.features.keys())
        return self.wdth

    def length(self):
        if self.lgth==-1:
            self.lgth=pow(self.dotprod(self.array),0.5)
        return self.lgth

    def dotprod(self,anarray):
        return self.array.multiply(anarray).sum()

    def cosine(self,avector):

        if (self.length()*avector.length())==0:
            sim =0
            print "Warning: zero length vectors"
        else:
            dotprod=self.dotprod(avector.array)
            sim = dotprod/(self.length()*avector.length())

        return sim


class SimEngine():
    #holds dictionary of vectors - manages conversion to sparse arrays and similarity calculations
    #requires a dictionary of filenames (anyRefKey->filename) for initialisation
    #optionally also an include_function - otherwise all vectors will be loaded from all files
    #allpairs similarities is cosine similarity within each file
    #todo: pointwise similarities is for same key from different files e.g., to compare composed with observed vectors

    def __init__(self,filename_dict,include_function=isAny):
        self.filenames=filename_dict
        self.vectors={} #dictionaries of vectors
        self.allfeatures={} #dictionary of all features observed for matrix generation
        self.fk_idx={} #feature --> dimension
        self.include_fn=include_function
        for type in self.filenames.keys():
            self.vectors[type]={}
        for type in self.filenames.keys():
            self.load(type)
            print "Loaded %s vectors from file %s" %(str(len(self.vectors[type].keys())),type)
        print "Converting to matrix form"
        self.makematrix() #sparse array generation

    def load(self,type):

        vectorfile=self.filenames[type]

        with open(vectorfile) as instream:
            for line in instream:
                line=line.rstrip()
                fields=line.split('\t')
                token=fields[0]
                featurelist=fields[1:]
                if self.include(token):
                    self.createvector(token,featurelist,type)
                else:
                    #print "Ignoring : "+token
                    pass

    def include(self,token):
        return self.include_fn(token)

    def createvector(self,token,featurelist,type):
        self.vectors[type][token]=WordVector(token)
        featurelist.reverse() #reverse list so can pop features and scores off

        while (len(featurelist)>0):
            f=featurelist.pop()
            sc=featurelist.pop()
            self.vectors[type][token].addfeature(f,sc)
            self.allfeatures[f]=1

    def makematrix(self):
        fkeys=self.allfeatures.keys()
        fkeys.sort()  #don't actually need to sort - but makes the indexing more predictable
        for i in range(len(fkeys)):
            self.fk_idx[fkeys[i]]=i

        #del self.allfeatures #can now delete this dictionary if memory is an issue

        dim=len(self.fk_idx)
        print "Dimensionality is " + str(dim)

        for type in self.vectors.keys():
            for wordvector in self.vectors[type].values():
                wordvector.makearray(self.fk_idx)

    def allpairs(self):
        for type in self.vectors.keys():
            for wordA in self.vectors[type].keys():
                for wordB in self.vectors[type].keys():
                    sim=self.vectors[type][wordA].cosine(self.vectors[type][wordB])
                    print "cosine(%s,%s) = %s [%s]"%(wordA,wordB,str(sim),type) #probably want to pipe this to a file


if __name__=="__main__":

    filename_dict={}
    key=1
    for filename in sys.argv[1:]:
        filename_dict[str(key)]=filename
        key+=1
    mySimEngine=SimEngine(filename_dict)
    mySimEngine.allpairs()