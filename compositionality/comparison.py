__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds
from simEngine import SimEngine

class Comparator():
    key1="observed"
    key2="composed"
    offsetting=1.0 #offset the dependency vector before composition (False/0 for baseline)

    def __init__(self,configfile):
        self.config=ConfigParser.RawConfigParser()
        self.config.read(configfile)

        self.exp_type=self.config.get('default','exp_type') # set this to anything other than 'compounds' if not wanting to load the phrasal compounds and use other config options
        try:
            self.parentdir=self.config.get('default','parentdir')
        except:
            self.parentdir=""
        self.filenames={}
        self.filenames[Comparator.key1]=self.parentdir+self.config.get('default','observedfile')
        if self.exp_type=="compounds":
            self.setup_compounds_exp(configfile)
            self.compounder.generate(self.rels,outfile=self.testcompoundfile) #generate list of compounds from observed file
            print self.compounder.generated_compounds
            print self.compounder.firstindex.keys()
            print self.compounder.secondindex.keys()
            if self.crossvalidate:
                self.compounder.setup_folds(self.nfolds)

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


        self.freqfile=self.filenames[Comparator.key1]+self.reducestring[Comparator.key1]+".rtot"
        for type in self.filenames.keys():
            self.filenames[type]=self.filenames[type]+self.reducestring.get(type,"")+self.normstring+self.weightingstring

        try:
            self.nfolds=int(self.config.get('default','nfolds'))
            trialp=ast.literal_eval(self.config.get('default','trialp'))
            self.crossvalidate=True
            self.paramdict={}
            self.paramdict["offsetting"]=trialp
        except:
            self.nfolds=0
            self.paramdict={}
            self.crossvalidate=False

            try:
                offsetting=float(self.config.get('default','offsetting'))
            except:
                offsetting=Comparator.offsetting
            self.paramdict["offsetting"]=offsetting


        if self.crossvalidate:
            print "Cross-valiation: number of folds = "+str(self.nfolds)
            print self.paramdict
        else:
            print "No cross-validation"
            print self.paramdict

    def generate_SimEngine(self):

        if self.exp_type==('compounds'):
            simEngine=SimEngine(self.filenames,self.isListedCompound,pathdelim=self.composer.pathdelims[0],saliency=self.composer.saliency,saliencyperpath=self.composer.saliencyperpath)
        elif self.exp_type==('simple_compounds'):
            simEngine=SimEngine(self.filenames,self.isCompound,pathdelim=self.composer.pathdelims[0],saliency=self.composer.saliency,saliencyperpath=self.composer.saliencyperpath)

        return simEngine

    def isCompound(self,token):
        return len(token.split('|'))==3

    def isListedCompound(self,token):
        return len(token.split('|'))==3 and token in self.compounder.generated_compounds


    def isConstituent(self,token):
        lex =token.split('/')[0]
        return lex in self.compounder.firstindex.keys() or lex in self.compounder.secondindex.keys()

    def loadFreqs(self):
        print("Loading "+self.freqfile+" for frequency analysis")
        with open(self.freqfile) as instream:
            for line in instream:
                line=line.rstrip()
                fields=line.split('\t')
                if len(fields[0].split('|'))==3:
                    self.compounder.addFreq(fields[0],float(fields[1]))




    def calcInternalSims(self):
        filenames={Comparator.key1:self.filenames[Comparator.key1]}
        print "Starting calculation of constituent similarities"
        aSimEngine=SimEngine(filenames,include_function=self.isConstituent)
        with open("intsims","w") as outstream:
            aSimEngine.allpairs(outstream=outstream)
        with open("intsims","r") as instream:
            for line in instream:
                line=line.rstrip()
                fields=line.split('\t')
                self.compounder.addIntSim(fields[1],fields[2],float(fields[3]))

    def correlate(self,instream,parampair=('','')):

        for line in instream:
            line=line.rstrip()
            fields=line.split('\t')
            if fields[1]==Comparator.key1 and fields[2]== Comparator.key2:
                self.compounder.addAutoSim(fields[0],fields[3])

        self.compounder.correlate()
        if self.crossvalidate:
             self.compounder.crossvalidate(self.nfolds,p=str(parampair[0])+":"+str(parampair[1]))


    def run(self):
        if self.exp_type=='compounds':

            for key in self.paramdict.keys():
                for value in self.paramdict[key]:
                    self.composer.run(parampair=(key,value))  #run composer to create composed vectors
                    self.mySimEngine.addfile(Comparator.key2,self.composer.outfile)  #add composed vector file to SimEngine
                    with open("testout","w") as outstream:
                        self.mySimEngine.pointwise(outstream)

                    self.loadFreqs()
                    self.calcInternalSims()
                    with open("testout",'r') as instream:
                        self.correlate(instream,parampair=(key,value))


        else:
            self.mySimEngine.allpairs()



if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()