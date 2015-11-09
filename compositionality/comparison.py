__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds
from simEngine import SimEngine

class Comparator():
    key1="observed"
    key2="composed"

    def __init__(self,configfile):
        self.config=ConfigParser.RawConfigParser()
        self.config.read(configfile)

        self.exp_type=self.config.get('default','exp_type') # set this to anything other than 'compounds' if not wanting to load the phrasal compounds and use other config options
        self.filenames={}
        self.filenames[Comparator.key1]=self.config.get('default','observedfile')
        if self.exp_type=="compounds":
            self.setup_compounds_exp(configfile)
            self.compounder.generate(self.rels,outfile=self.testcompoundfile) #generate list of compounds from observed file
            print self.compounder.generated_compounds

        self.mySimEngine=self.generate_SimEngine()  #will load observed vectors

    def setup_compounds_exp(self,configfile):
        self.compounder=compounds.Compounder(configfile)
        self.composer = nouncompounds.NounCompounder(["config",configfile])
        self.rels=ast.literal_eval(self.config.get('default','rels'))
        self.testcompoundfile=self.config.get('compounder','compound_file')

        self.reducestring={}
        self.reducestring[Comparator.key1]=".nouns.reduce_1_1"
        self.normstring=".filtered"
        if self.composer.normalised:
            self.normstring+=".norm"
        if self.composer.weighting in ['smooth_ppmi','ppmi','pnppmi','gof_ppmi']:
            self.weightingstring="."+self.composer.weighting
        if self.composer.ppmithreshold>0:
            self.weightingstring+="_"+str(self.composer.ppmithreshold)
        else: self.weightstring=""

        for type in self.filenames.keys():
            self.filenames[type]=self.filenames[type]+self.reducestring.get(type,"")+self.normstring+self.weightingstring

    def generate_SimEngine(self):

        if self.exp_type==('compounds'):
            simEngine=SimEngine(self.filenames,self.isListedCompound,pathdelim=self.composer.pathdelim)
        elif self.exp_type==('simple_compounds'):
            simEngine=SimEngine(self.filenames,self.isCompound,pathdelim=self.composer.pathdelim)

        return simEngine

    def isCompound(self,token):
        return len(token.split('|'))==3

    def isListedCompound(self,token):
        return len(token.split('|'))==3 and token in self.compounder.generated_compounds


    def correlate(self,instream):

        for line in instream:
            line=line.rstrip()
            fields=line.split('\t')
            if fields[1]==Comparator.key1 and fields[2]== Comparator.key2:
                self.compounder.addAutoSim(fields[0],fields[3])

        self.compounder.correlate()


    def run(self):
        if self.exp_type=='compounds':
            self.composer.run()  #run composer to create composed vectors
            self.mySimEngine.addfile(Comparator.key2,self.composer.outfile)  #add composed vector file to SimEngine
            with open("testout","w") as outstream:
                self.mySimEngine.pointwise(outstream)

            with open("testout",'r') as instream:
                self.correlate(instream)
        else:
            self.mySimEngine.allpairs()



if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()