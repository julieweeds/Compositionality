__author__ = 'juliewe'
#compare observed and composed vectors, correlate with compositionality judgements

import compounds,sys, ConfigParser,ast, nouncompounds, numpy as np, composition,math
from simEngine import SimEngine

def getValue(text):
    # extract 0.25 from offsetting:0.25
    fields=text.split(":")
    return float(fields[1])


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
        try:
            self.simmetric=self.config.get('default','simmetric')
        except:
            self.simmetric="cosine"
        if self.exp_type=="compounds":
            self.setup_compounds_exp(configfile)

            if 'observed' not in self.skip:
                self.compounder.readcompounds()
                self.loadFreqs(self.rels,outfile=self.testcompoundfile)
            else:
                self.compounder.generate(self.rels,outfile=self.testcompoundfile) #generate list of compounds from observed file
            print len(self.compounder.generated_compounds),self.compounder.generated_compounds
            #if self.crossvalidate:
            #    self.compounder.setup_folds(self.nfolds)
            # do this later

        if 'revectorise' not in self.skip:
            print "Revectorising observed phrasal vectors"
            self.revectorise_observed(configfile,self.compounder.generated_compounds)

    def setup_compounds_exp(self,configfile):

        try:
            self.skip=ast.literal_eval(self.config.get('default','skip'))
        except:
            self.skip=[]
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
            self.offsetting=float(self.config.get('default','offsetting'))
        except:
            self.offsetting=Comparator.offsetting

        try:
            self.nfolds=int(self.config.get('default','nfolds'))
            trialp=ast.literal_eval(self.config.get('default','trialp'))
            self.crossvalidate=True
            self.paramdict={}

            try:
                self.cv_param=self.config.get('default','cv_param')
            except:
                self.cv_param="offsetting"
            self.paramdict[self.cv_param]=trialp
            try:
                self.repetitions=int(self.config.get('default','repetitions'))
            except:
                self.repetitions=1
        except:
            self.nfolds=0
            self.paramdict={}
            self.crossvalidate=False
            self.paramdict["offsetting"]=[self.offsetting]


        if self.crossvalidate:
            print "Cross-validation: number of folds = "+str(self.nfolds)
            print "Number of repetitions = "+str(self.repetitions)
            print self.paramdict
            print "Default off-setting: ",self.offsetting

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

                            if self.compounder.addFreq(fields[0],float(fields[1])):
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
            reps=self.repetitions
            m=[]
            while reps>0:
                reps=reps-1
                m+=self.compounder.crossvalidate(self.nfolds,p=str(parampair[0])+":"+str(parampair[1]),rep=reps)
            return m
        else:
            return []

    def analyse(self,cv_matrix):
        #print cv_matrix

        testrs=[]
        testps=[]
        #analyse training performance
        #for each fold find best parameter
        #for that fold and parameter collect test performance
        folds = self.nfolds*self.repetitions
        for i in range(0,folds):
            besttraining=0  #make 1 for worst, 0 for best
            bestindex=-1
            for index,line in enumerate(cv_matrix):
                if line[1]==i:
                    if line[2]>besttraining:  #make < for worst, > for best
                        besttraining=line[2]
                        bestindex=index
            testrs.append(cv_matrix[bestindex][3])
            testps.append(getValue(cv_matrix[bestindex][0]))

        perf=np.mean(testrs)
        error=np.std(testrs)/math.sqrt(folds)
        print "Cross-validated performance over %s repetitions is %s with error %s"%(str(len(testrs)),str(perf),str(error))
        mp=np.mean(testps)
        msd=np.std(testps)
        print "Mean Chosen parameter settings: ",str(mp),str(msd)


    def run(self):
        if self.exp_type=='compounds':
            cv_matrix=[]
            for key in self.paramdict.keys():
                for value in self.paramdict[key]:
                    if 'compose' not in self.skip:
                        print "Running composer"
                        self.composer.run(parampair=(key,value))  #run composer to create composed vectors
                        self.composer.close()
                    else:
                        self.composer.outfile=self.composer.getComposedFilename(parampair=(key,value))

                    simfile=self.composer.outfile+".sims"

                    if 'sim' not in self.skip:
                        print "Running sim engine"
                        print "Reloading observed phrasal vectors"
                        self.mySimEngine=self.generate_SimEngine()  #will load observed vectors

                        self.mySimEngine.addfile(Comparator.key2,self.composer.outfile)  #add composed vector file to SimEngine
                        with open(simfile,"w") as outstream:
                            self.mySimEngine.pointwise(outstream,simmetric=self.simmetric)

                    #self.calcInternalSims()
                    if 'correlate' not in self.skip:
                        print "Running correlation"
                        with open(simfile,'r') as instream:
                            m=self.correlate(instream,parampair=(key,value))
                        if len(m)>0:
                            for line in m:
                                cv_matrix.append(line)


            if len(cv_matrix)>0:
                self.analyse(cv_matrix)

        else:
            print "Reloading observed phrasal vectors"
            self.mySimEngine=self.generate_SimEngine()  #will load observed vectors
            self.mySimEngine.allpairs()



if __name__=="__main__":
    myComparator=Comparator(sys.argv[1])
    myComparator.run()