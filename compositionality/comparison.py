__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds
import numpy as np, scipy.sparse as sparse

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



class Comparator():

    def __init__(self,configfile):
        self.config=ConfigParser.RawConfigParser()
        self.config.read(configfile)
        self.vectors={} #dictionaries of vectors
        self.vectors['observed']={}
        self.allfeatures={} #dictioray of all features observed - needed for matrix generation
        self.fk_idx={} #feature--> dimension

        self.exp_type=self.config.get('default','exp_type') # set this to anything other than 'compounds' if not wanting to load the phrasal compounds and use other config options
        self.filenames={}
        self.filenames['observed']=self.config.get('default','observedfile')

        if self.exp_type=="compounds":
            self.setup_compounds_exp(configfile)

    def setup_compounds_exp(self,configfile):
        self.compounder=compounds.Compounder(configfile)
        self.composer = nouncompounds.NounCompounder(["config",configfile])
        self.rels=ast.literal_eval(self.config.get('default','rels'))
        self.testcompoundfile=self.config.get('compounder','compound_file')
        self.vectors['composed']={}

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

    def load(self,type):

        vectorfile=self.filenames[type]

        with open(vectorfile) as instream:
            for line in instream:
                line=line.rstrip()
                fields=line.split('\t')
                token=fields[0]
                featurelist=fields[1:]
                if self.include(token,type):
                    self.createvector(token,featurelist,type)
                else:
                    #print "Ignoring : "+token
                    pass

    def include(self,token,type):

        if self.exp_type==('compounds'):

            if len(token.split('|'))==3:
                if token in self.compounder.generated_compounds:
                    return True
            return False
        elif self.exp_type==('simple_compounds'):
            if len(token.split('|'))==3:
                return True
            else:
                return False

        return True



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
        fkeys.sort()
        for i in range(len(fkeys)):
            self.fk_idx[fkeys[i]]=i

        #del self.allfeatures #can now delete this dictionary if memory is an issue

        dim=len(self.fk_idx)
        print "Dimensionality is " + str(dim)

        for type in self.vectors.keys():
            for wordvector in self.vectors[type].values():
                wordvector.makearray(self.fk_idx)


    def test1(self):

        for wordA in self.vectors["observed"].keys():
            for wordB in self.vectors["observed"].keys():
                sim=self.vectors["observed"][wordA].cosine(self.vectors["observed"][wordB])
                print "Cosine("+wordA+","+wordB+") = "+str(sim)

    def run(self):
        if self.exp_type=='compounds':
            self.compounder.generate(self.rels,outfile=self.testcompoundfile)
            print self.compounder.generated_compounds

        self.load('observed')
        print "Loaded %s vectors, converting to matrix form" %(str(len(self.vectors['observed'].keys())))
        self.makematrix() #sparse array generation
        self.test1()


if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()