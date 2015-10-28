__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds
import numpy as np, scipy.sparse as sparse

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

    def __init__(self,filename_dict,include_function=isAny):
        self.filenames=filename_dict
        self.vectors={} #dictionaries of vectors
        self.allfeatures={} #dictionary of all features observed for matrix generation
        self.fk_idx={} #feature --> dimension
        self.include_fn=include_function
        for type in self.filenames.keys():
            self.vectors[type]={}


    def load(self,type,exp_type="simple"):

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
                    print "cosine(%s,%s) = %s [%s]"%(wordA,wordB,str(sim),type)


    def run(self):
        for type in self.filenames.keys():
            self.load(type)
            print "Loaded %s vectors" %(str(len(self.vectors['observed'].keys())))
        print "Converting to matrix form"
        self.makematrix() #sparse array generation
        self.allpairs()

class Comparator():

    def __init__(self,configfile):
        self.config=ConfigParser.RawConfigParser()
        self.config.read(configfile)

        self.exp_type=self.config.get('default','exp_type') # set this to anything other than 'compounds' if not wanting to load the phrasal compounds and use other config options
        self.filenames={}
        self.filenames['observed']=self.config.get('default','observedfile')
        if self.exp_type=="compounds":
            self.setup_compounds_exp(configfile)
        self.mySimEngine=self.generate_SimEngine()

    def setup_compounds_exp(self,configfile):
        self.compounder=compounds.Compounder(configfile)
        self.composer = nouncompounds.NounCompounder(["config",configfile])
        self.rels=ast.literal_eval(self.config.get('default','rels'))
        self.testcompoundfile=self.config.get('compounder','compound_file')

        self.reducestring={}
        self.reducestring['observed']=".nouns.reduce_1_1"
        self.normstring=".filtered"
        if self.composer.normalised:
            self.normstring+=".norm"
        if self.composer.weighting in ['smooth_ppmi','ppmi','pnppmi','gof_ppmi']:
            self.weightingstring="."+self.composer.weighting
        else: self.weightstring=""

        for type in self.filenames.keys():
            self.filenames[type]=self.filenames[type]+self.reducestring.get(type,"")+self.normstring+self.weightingstring

    def generate_SimEngine(self):

        simEngine=SimEngine(self.filenames)

        if self.exp_type==('compounds'):
            simEngine=SimEngine(self.filenames,self.isListedCompound)
        elif self.exp_type==('simple_compounds'):
            simEngine=SimEngine(self.filenames,self.isCompound)

        return simEngine

    def isCompound(self,token):
        return len(token.split('|'))==3

    def isListedCompound(self,token):
        return len(token.split('|'))==3 and token in self.compounder.generated_compounds


    def run(self):
        if self.exp_type=='compounds':
            self.compounder.generate(self.rels,outfile=self.testcompoundfile)
            print self.compounder.generated_compounds
        self.mySimEngine.run()



if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()