__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds, numpy as np, composition,math
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
            #self.compounder.generate(self.rels,outfile=self.testcompoundfile) #generate list of compounds from observed file
            self.compounder.readcompounds()
            self.loadFreqs(self.rels,outfile=self.testcompoundfile)
            print self.compounder.generated_compounds
            if self.crossvalidate:
                self.compounder.setup_folds(self.nfolds)
        print "Revectorising observed phrasal vectors"
        self.revectorise_observed(configfile,self.compounder.generated_compounds)
        print "Reloading observed phrasal vectors"
        self.mySimEngine=self.generate_SimEngine()  #will load observed vectors

    def setup_compounds_exp(self,configfile):
        self.compounder=compounds.Compounder(configfile)
        self.composer = nouncompounds.NounCompounder(["config",configfile])
        self.rels=ast.literal_eval(self.config.get('default','rels'))
        self.testcompoundfile=self.config.get('compounder','compound_file')

        self.reducestring={}
        self.reducestring[Comparator.key1]=".nouns.reduce_0_2"
        self.normstring=".filtered"
        if self.composer.normalised:
            self.normstring+=".norm"


        if self.composer.weighting in ['smooth_ppmi','ppmi','pnppmi','gof_ppmi']:
            self.weightingstring="."+self.composer.weighting
        if self.composer.ppmithreshold>0:
            self.weightingstring+="_"+str(self.composer.ppmithreshold)
        else: self.weightstring=""
        #self.weightingstring=""

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
            self.paramdict["offsetting"]=[offsetting]


        if self.crossvalidate:
            print "Cross-valiation: number of folds = "+str(self.nfolds)
            print self.paramdict
        else:
            print "No cross-validation"
            print self.paramdict

    def revectorise_observed(self,configfile,phraselist):
        vectoriser=composition.Composition(["config",configfile])
        vectoriser.options=['revectorise']
        vectoriser.run(phraselist)

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
        return lex in self.composer.getLeftIndex() or lex in self.composer.getRightIndex()

    def loadFreqs(self,rel_list,outfile): #should be part of compounder
        print("Loading "+self.freqfile+" for frequency analysis")
        with open(outfile,"w") as outstream:
            self.compounder.generated_compounds=[]
            with open(self.freqfile) as instream:
                for line in instream:
                    line=line.rstrip()
                    fields=line.split('\t')
                    parts=fields[0].split('|')
                    if len(parts)==3 and parts[1] in rel_list:
                        posparts=parts[2].split('/')
                        if len(posparts)==2:

                            self.compounder.addFreq(fields[0],float(fields[1]))
                            self.compounder.generated_compounds.append(fields[0])
                            outstream.write(fields[0]+"\n")




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

        self.compounder.correlate(show_graph=(not self.crossvalidate))
        if self.crossvalidate:
            return self.compounder.crossvalidate(self.nfolds,p=str(parampair[0])+":"+str(parampair[1]))

        else:
            return []

    def analyse(self,cv_matrix):
        print cv_matrix

        testrs=[]
        testps=[]
        #analyse training performance
        #for each fold find best parameter
        #for that fold and parameter collect test performance
        for i in range(0,self.nfolds):
            maxtraining=0
            maxindex=-1
            for index,line in enumerate(cv_matrix):
                if line[1]==i:
                    if line[2]>maxtraining:
                        maxtraining=line[2]
                        maxindex=index
            testrs.append(cv_matrix[maxindex][3])
            testps.append(cv_matrix[maxindex][0])

        perf=np.mean(testrs)
        error=np.std(testrs)/math.sqrt(self.nfolds)
        print "Cross-validated performance",str(perf),str(error)
        print "Chosen parameter settings: ",testps

    def run(self):
        if self.exp_type=='compounds':
            cv_matrix=[]
            for key in self.paramdict.keys():
                for value in self.paramdict[key]:
                    self.composer.run(parampair=(key,value))  #run composer to create composed vectors
                    self.mySimEngine.addfile(Comparator.key2,self.composer.outfile)  #add composed vector file to SimEngine
                    with open("testout","w") as outstream:
                        self.mySimEngine.pointwise(outstream)

                    #self.calcInternalSims()
                    with open("testout",'r') as instream:
                        m=self.correlate(instream,parampair=(key,value))
                    if len(m)>0:
                        for line in m:
                            cv_matrix.append(line)


            if len(cv_matrix)>0:
                self.analyse(cv_matrix)

        else:
            self.mySimEngine.allpairs()



if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()