__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds
from simEngine import SimEngine

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
        self.mySimEngine.allpairs()



if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()