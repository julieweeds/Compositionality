__author__ = 'juliewe'
#contains a number of tools useful for working with APT vector files - started March 2015
#split : by POS
#maketotals : for PPMI calculation/Byblo
#filter : only words in lists (nouns and adjectives) and/or with particular frequency
#compose : words in lists using simple add composition
#reduceorder : only retain features with given orders
#revectorise : output PPMI vectors or PNPPMI vectors
#normalise: convert counts to probabilities

#vectors are displayed via their most salient features



import sys,math,gzip,ast
import ConfigParser

from operator import itemgetter

try:
    import graphing
except IOError,ImportError:
    print "Warning: Unable to import graphing module"

try:
    import yaml
except ImportError,IOError:
    print "Warning: Unable to import yaml for reading composition pair file"


#functions on features
#----
#get the path prefix / dependency path of a given feature
#self.getpathtype("amod:red") = amod
#self.getpathtype("_dobj>>amod:red") = _dobj>>amod
#----
def getpathtype(feature):
    #get the path of a given feature
    fields=feature.split(":")
    return fields[0]

#----
#get the path value of a given feature
#self.getpathvalue("amod:red")=red
#self.getpathvalue("_dobj>>amod:red") = red
#----
def getpathvalue(feature):
    fields=feature.split(":")
    if len(fields)>1:
        return ":"+fields[1]
    else:
        #print "No feature value for "+feature
        return ""

#----
#get the order of a given feature
#----
def getorder(feature,delim="\xc2\xbb"):

    path=getpathtype(feature)

    if path=="":
        order=0
    else:

        fields=path.split(delim)
        order=len(fields)
    return order

def convert(feature,delims=[]):
    if len(delims)<2:
        return feature
    else:
        for delim in delims[1:]:
            #print "Attempting to convert ",feature,delim
            feature = feature.replace(delim,delims[0])
            #print feature
        return feature

#---
#split a higher order feature / find path prefix
#0th order e.g., :red => return ("","")
#1st order e.g., amod:red =>  return ("amod:red","")
#2nd order e.g., _dobj>>amod:red => return ("_dobj","amod:red:)
#3rd order e.g., nsubj>>_dobj>>amod:red => return ("nsubj","dobj>>amod:red")
#----
def splitfeature(feature,delim="\xc2\xbb"):


    path=getpathtype(feature)

    if path=="":
        return "",""
    else:
        path=getpathtype(feature)
        fields=path.split(delim)

        if len(fields)>1:
            text=fields[1]
            if len(fields)>2:
                for field in fields[2:]:
                    text+=delim+field
            return fields[0],text
        else:
            return fields[0],""

#---
#turn alist into a string concatenated using the achar
#---
def join (alist,achar):
    if len(alist)>1:
        astring=alist[0]
        for element in alist[1:]:
            astring+=achar+element
        return astring

    elif len(alist)==1:
        return alist[0]
    else:
        return ""


def checkphraseformat(line):

    entry=line.split("\t")[0]
    phrasefields=entry.split('|')
    if len(phrasefields)==3:
        depparts=phrasefields[0].split("/")
        if len(depparts)==2:
            restline=line[len(entry):]
            entry=depparts[0]+"|"+phrasefields[1]+"-"+depparts[1]+"|"+phrasefields[2]
            line=entry+restline
            return line
        else:
            return line

    else:
        return line

class Composition:

    nouns=[]
    adjectives=[]
    verbs=[]
    adverbs=[]
    others=[]

    includedtypes=[]  #all types of relation will be considered if this is empty
    featmax=3  #how many features of each path type to display when showing most salient features
    display=True
    allphrases=True  #include allphrases regardless of frequency

    headp=0.5
    reversed=False #reverse the order of the constituents before composition (True for order baseline)
    compop="add"
    offsetting=1.0 #offset the dependency vector before composition (False/0 for baseline)

    ppmithreshold=0
    filterfreq=1000
    saliency=0
    saliencyperpath=False

    headPoS={"nn":"N","amod":"N","mod":"N","nn-n":"N","amod-j":"N","amod-v":"N"}
    depPoS={"nn":"N","amod":"J","mod":"J","nn-n":"N","amod-j":"J","amod-v":"V"}
    basicRel={"nn":"nn","amod":"amod","mod":"mod","nn-n":"nn","amod-j":"amod","amod-v":"amod"}


    def __init__(self,options):

        print options
        if options[0] == "config":
            self.configure(options[1]) # configure via configuration file
        else:
            #configure via command line options
            # parameter0 = function
            self.options=[options[0]]

            # parameter1 or default = original input file name (in current working directory)
            if len(options)>1:
                self.inpath=options[1]
            else:
                print "Requires the base filename as second input"

            #parameter2 = pos to be considered (necessary for all functions other than split).  If not given = N
            if len(options)>2:
                self.pos=options[2]
            else:
                self.pos="N"

            #parameter3 = minimum order of dependency relation to include
            #parameter4 = maximum order of dependency relation to include
            #parameter3 can = X if you want to set further options without setting this (i.e., work with all orders given)

            if len(options)>3 and not options[3]=="X":

                self.minorder=int(options[3])
                self.maxorder=int(options[4])
                self.reducedstring=".reduce_"+str(self.minorder)+"_"+str(self.maxorder)
            else:
                self.minorder=0
                self.maxorder=2
                self.reducedstring=""

            #optional parameters.  Can be anywhere from parameter4 onwards
            #type of ppmi calculation: default = ppmi, alternatives are gof_ppmi (where probabilities are calculated over all types rather than on a path type basis)
            #   and pnppmi (where standard ppmi calculation is multiplied by path probability)
            #normalised: this is required so that the result of the normalisation function can be used as input to a future function e.g. revectorisation or composition

            self.pp_normal= "pp_normalise" in options or "pnppmi" in options #include one of these flags in order for PPMI values to be multiplied by path probability in final vectors
            self.gof_ppmi="gof_ppmi" in options
            self.normalised = "normalise" in options or "normalised" in options #this may be the main option (to carry out normalisation) or be included as one of the optional options so that normalised counts are used

            self.ppmithreshold=Composition.ppmithreshold
            self.rfilterfreq=Composition.filterfreq
            self.cfilterfreq=Composition.filterfreq
            self.saliency=Composition.saliency
            self.saliencyperpath=Composition.saliencyperpath
            self.display=Composition.display
            self.parentdir=""
            self.filenames=[]
            self.pathdelims=['\xc2\xbb']

          #suffixes for pos
        self.filesbypos={"N":".nouns","V":".verbs","J":".adjs","R":".advs","F":".others","ANS":".ans"}

        #these are dictionaries which will hold vectors and totals
        self.vecsbypos={}
        self.totsbypos={}
        self.feattotsbypos={}
        self.pathtotsbypos={}
        self.typetotsbypos={}

        for pos in self.filesbypos.keys():
            self.vecsbypos[pos]={}
            self.totsbypos[pos]={}
            self.feattotsbypos[pos]={}
            self.pathtotsbypos[pos]={}
            self.typetotsbypos[pos]={}



        self.includedtypes=Composition.includedtypes


    def configure(self,filename):
        #load and configure
        print "Reading configuration from "+filename
        self.config=ConfigParser.RawConfigParser()
        self.config.read(filename)

        self.options=ast.literal_eval(self.config.get('default','options'))
        try:
            self.parentdir=self.config.get('default','parentdir')
        except:
            self.parentdir=""
        try:
            filenames=ast.literal_eval(self.config.get('default','filenames'))
            #print filenames
            self.inpaths=[self.parentdir+p for p in filenames]
            self.inpath=self.parentdir+"combo"

        except:
            self.inpaths=[self.parentdir+self.config.get('default','filename')]
            self.inpath=self.inpaths[0]
        print "Paths to process: ",self.inpaths
        self.pos=self.config.get('default','pos')
        mini = self.config.get('default','minorder')
        maxi = self.config.get('default','maxorder')
        if mini == "X":
            self.minorder=0
            self.maxorder=2
            self.reducedstring=""
        else:
            self.minorder=int(mini)
            self.maxorder=int(maxi)
            self.reducedstring=".reduce_"+str(mini)+"_"+str(maxi)

        self.weighting=self.config.get('default','weighting')
        self.pp_normal=(self.weighting=="pnppmi" or self.weighting=="pp_normalise")
        self.gof_ppmi=(self.weighting=="gof_ppmi")
        self.smooth_ppmi=(self.weighting=="smooth_ppmi" or self.weighting=="smoothed_ppmi")
        self.normalised=(self.config.get('default','normalised')=="True") or self.options[0]=="normalise"
        self.ppmithreshold=float(self.config.get('default','wthreshold'))
        self.saliency=int(self.config.get('default','saliency'))
        self.saliencyperpath=self.config.get('default','saliencyperpath')
        try:
            self.rfilterfreq=int(self.config.get('default','rfthreshold'))
            self.cfilterfreq=int(self.config.get('default','cfthreshold'))
        except:

            self.rfilterfreq=int(self.config.get('default','fthreshold'))
            self.cfilterfreq=int(self.config.get('default','fthreshold'))
        self.comppairfile=self.config.get('default','comppairfile')
        self.filterfile=self.config.get('default','filterfile')
        try:
            self.display=(self.config.get('default','display')=='True')
        except:
            self.display=Composition.display


        try:
            self.pathdelims=ast.literal_eval(self.config.get('default','path_delims'))
        except:
            self.pathdelims=[]
            try:
                self.pathdelims=[self.config.get('default','path_delim')]
            except:
                self.pathdelims=['\xc2\xbb']

        try:
            self.allphrases=(self.config.get('default','allphrases')=='True')
        except:
            self.allphraser=Composition.allphrases



        try:
            self.headp=float(self.config.get('default','headp'))
        except:
            self.headp=Composition.headp

        try:
            self.reversed=(self.config.get('default','reversed')=='True')
        except:
            self.reversed=Composition.reversed

        try:
            self.compop=self.config.get('default','compop')
        except:
            self.compop=Composition.compop
        try:
            self.offsetting=float(self.config.get('default','offsetting'))
        except:
            self.offsetting=Composition.offsetting

        try:
            self.distinguish=(self.config.get('default','distinguish')=='True')
        except:
            self.distinguish=False

        print "Composition offsetting: ",self.offsetting

    #----HELPER FUNCTIONS

    #-----
    # set the words of interest
    #1) if self.comppairfile set, add the appropriate member of pair form self.comppairlist to self.words
    #2) elseif filterfile not present, add any nouns/adjectives in defaults
    # 3) elseif filterfile is present, add these words to self.words
    #-----
    def set_words(self,words=[]):


        if words==[]:
            if self.comppairfile!="":
                if self.pos=="J":
                    index=2
                else:
                    index=0
                self.words=[]
                for pair in self.comppairlist:
                    if not pair[index] in self.words:
                        self.words.append(pair[index])


            elif self.filterfile=="":
                if self.pos=="N":
                    self.words=Composition.nouns
                elif self.pos=="J":
                    self.words=Composition.adjectives
                elif self.pos=="V":
                    self.words=Composition.verbs
                elif self.pos=="R":
                    self.words=Composition.adverbs
                else:
                    self.words=Composition.others
            else:
                with open(self.filterfile) as fp:
                    self.wordlistlist=yaml.safe_load(fp)
                self.words=[]
                for wordlist in self.wordlistlist:
                    self.words+=wordlist
        else:
            self.words=words
        print "Setting words of interest: ",self.words


    #----
    #boolean function as to whether a word is in self.words or self.words=[]
    #-----
    def include(self,word):
        if len(self.words)==0:
            return True
        elif word in self.words:
            return True
        else:
            return False


    #----
    #boolean function to decide whether an entry is a phrase
    #----
    def phraseinclude(self,entry):
        #do not want to filter out phrases even if low frequency
        if self.allphrases:
            if len(entry.split('|'))==3:
                return True
            else:
                return False
        else:
            return False #don't use this functionality

    #---
    #boolean function as to whether a pathtype is in self.includedtypes or self.includedtypes=[]
    #----
    def typeinclude(self,pathtype):
        if self.includedtypes==[]:
            return True
        else:
            return pathtype in self.includedtypes

    #---
    #generate the input file string according to POS
    #---
    def selectpos(self,path=""):
        if path =="":
            path=self.inpath
        return path+self.filesbypos.get(self.pos,self.filesbypos["N"])


    #----
    #get the path prefix / dependency path of a given feature
    #self.getpathtype("amod:red") = amod
    #self.getpathtype("_dobj>>amod:red") = _dobj>>amod
    #----
    def getpathtype(self,feature):
        return getpathtype(feature)

    #----
    #get the path value of a given feature
    #self.getpathvalue("amod:red")=red
    #self.getpathvalue("_dobj>>amod:red") = red
    #----
    def getpathvalue(self,feature):
       return getpathvalue(feature)

    #----
    #get the order of a given feature
    #----
    def getorder(self,feature):
        return getorder(feature,delim=self.pathdelims[0])

    #---
    #split a higher order feature / find path prefix
    #0th order e.g., :red => return ("","")
    #1st order e.g., amod:red =>  return ("amod:red","")
    #2nd order e.g., _dobj>>amod:red => return ("_dobj","amod:red:)
    #3rd order e.g., nsubj>>_dobj>>amod:red => return ("nsubj","dobj>>amod:red")
    #----
    def splitfeature(self,feature):
        return splitfeature(feature,delim=self.pathdelims[0])

    #---
    #turn alist into a string concatenated using the achar
    #---
    def join (self,alist,achar):
        return join(alist,achar)

    #----MAIN FUNCTIONS

    #----
    #SPLIT
    #take the original gzipped file and split it by POS
    #----
    def splitpos(self):

        nouns=open(self.inpath+self.filesbypos["N"],"w")
        verbs=open(self.inpath+self.filesbypos["V"],"w")
        adjs=open(self.inpath+self.filesbypos["J"],"w")
        advs=open(self.inpath+self.filesbypos["R"],"w")
        others=open(self.inpath+self.filesbypos["F"],"w")


        infile=self.inpath+".gz"

        try:
            instream=gzip.open(infile)
        except:
            instream=open(self.inpath)

        for lines,line in enumerate(instream):

            line=line.rstrip()
            entry=line.split("\t")[0]

            try:
                pos=entry.split("/")[-1].lower()
            except:
                print "Cannot split "+entry+" on line "+str(lines)
                pos=""

            #pos=entry.split("/")[-1].lower()
            if pos.startswith("n"):
                nouns.write(checkphraseformat(line)+"\n")
            elif pos.startswith("v"):
                verbs.write(line+"\n")
            elif pos.startswith("j"):
                adjs.write(line+"\n")
            elif pos.startswith("r"):
                advs.write(line+"\n")
            else:
                others.write(line+"\n")
            if lines % 1000==0:print "Processed "+str(lines)+" lines"

        nouns.close()
        verbs.close()
        adjs.close()
        advs.close()
        others.close()
        instream.close()
        return

    #----
    #REDUCEORDER
    #generally used after SPLIT
    #remove dependencies which do not have an order within the thresholds set by minorder and maxorder
    #note that
    # :gunman/N etc is considered to be a 0th order dependency
    #_dobj:shoot/V is a 1st order dependency
    # _dobj:nsubj:man/N is a 2nd order dependency
    #-------
    def reduceorder(self):

        infile=self.selectpos()
        outfile=infile+self.reducedstring
        with open(outfile,"w") as outstream:
            with open(infile) as instream:
                for lines, line in enumerate(instream):
                    line=line.rstrip()
                    if lines % 1000==0:print "Processing line "+str(lines)
                    fields=line.split("\t")
                    entry=fields[0]
                    features=fields[1:]
                    outline=entry
                    nofeats=0
                    while len(features)>0:
                        freq=features.pop()
                        feat=convert(features.pop(),delims=self.pathdelims)
                        forder=self.getorder(feat)

                        if forder>=self.minorder and forder<=self.maxorder:
                            outline+="\t"+feat+"\t"+freq
                            nofeats+=1
                    #print entry, nofeats
                    if nofeats>0:
                        outstream.write(outline+"\n")
                        #print "written out"


    #----
    #MAKETOTALS
    #calculate row and column totals
    #this is usually done before filtering and also after filtering and normalisation
    #------

    def maketotals(self):

        if self.normalised:
            infile=self.selectpos()+self.reducedstring+".filtered"+".norm"
        else:
            infile= self.selectpos()+self.reducedstring
        rowtotals=infile+".rtot"
        coltotals=infile+".ctot"

        rows=open(rowtotals,"w")
        cols=open(coltotals,"w")

        featuretotals={}
        with open(infile) as instream:
            for lines,line in enumerate(instream):

                if lines%1000==0:print "Processing line "+str(lines)
                rowtotal=0.0
                line=line.rstrip()
                fields=line.split("\t")
                entry=fields[0]
                features=fields[1:]

                index=0
                while len(features)>0:
                    index+=1

                    freq=features.pop()
                    feat=convert(features.pop(),delims=self.pathdelims)

                    #print str(index)+"\t"+feat+"\t"+str(freq)
                    try:
                        freq=float(freq)
                        rowtotal+=freq
                        current=featuretotals.get(feat,0.0)
                        featuretotals[feat]=current+freq
                    except ValueError:
                        print "Error: "+str(index)+"\t"+feat+"\t"+str(freq)+"\n"
                        features=features+list(feat)



                rows.write(entry+"\t"+str(rowtotal)+"\n")

        for feat in featuretotals.keys():
            cols.write(feat+"\t"+str(featuretotals[feat])+"\n")


        rows.close()
        cols.close()

    #---
    #subsequenct functions in pipeline can load pre-calcualated row totals using this function
    #---
    def load_rowtotals(self):
        for path in self.inpaths:
            infile= self.selectpos(path=path)+self.reducedstring
            if self.normalised and not self.option=="normalise":
                infile+=".filtered.norm"
            rowtotals=infile+".rtot"
            totals={}
            print "Loading entry totals from: "+rowtotals
            with open(rowtotals) as instream:
                for line in instream:
                    line=line.rstrip()
                    fields=line.split("\t")
                    if self.normalised or float(fields[1])>self.rfilterfreq or self.phraseinclude(fields[0]):
                        sofar=totals.get(fields[0],0)
                        totals[fields[0]]=sofar+float(fields[1])
            print "Loaded "+str(len(totals.keys()))

        return totals

    #---
    #subsequent functions in pipeline can load pre-calculated column totals using this function
    #----
    def load_coltotals(self,cds=False):
        #set cds =True to perform context distribution smoothing

        for path in self.inpaths:
            infile=self.selectpos(path=path)+self.reducedstring
            if self.normalised and not self.option=="normalise":
                infile+=".filtered.norm"
            coltotals=infile+".ctot"
            totals={}
            print "Loading feature totals from: "+coltotals
            with open(coltotals) as instream:
                for line in instream:
                    line=line.rstrip()
                    fields=line.split("\t")
                    if self.normalised or float(fields[1])>self.cfilterfreq:
                        feat=convert(fields[0],delims=self.pathdelims)
                        sofar=totals.get(feat,0)
                        if cds:
                            totals[feat]=pow(float(fields[1]),0.75)+sofar
                        else:
                            totals[feat]=float(fields[1])+sofar
            print "Loaded "+str(len(totals.keys()))
        return totals


    #---
    #FILTER
    #filter by frequency and by words of interest
    #make sure filterfile and comppairfile etc are undefined if you want vectors for all words above threshold frequency
    #no threshold on individual event frequencies, only on entry/row and feature/column totals.
    #----
    def filter(self):

        infile = self.selectpos()+self.reducedstring
        outfile=infile+".filtered"

        coltotals= self.load_coltotals()
        savereducedstring=self.reducedstring
        self.reducedstring=".reduce_1_1" #always use same rowtotals for filtering whatever the reduction
        rowtotals = self.load_rowtotals()
        self.reducedstring=savereducedstring
        outstream=open(outfile,"w")
        print "Filtering for words ",self.words
        print "Filtering for frequency ",self.rfilterfreq, self.cfilterfreq
        todo=len(rowtotals)
        with open(infile) as instream:
            done=0
            for lines,line in enumerate(instream):
                line = line.rstrip()
                if lines%1000==0:
                    percent=done*100.0/todo
                    print "Processing line "+str(lines)+"("+str(percent)+"%)"
                fields=line.split("\t")
                #entry=fields[0].lower()
                entry=fields[0]
                features=fields[1:]
                entrytot=rowtotals.get(entry,0)
                nofeats=0
                if self.phraseinclude(entry) or( entrytot>self.rfilterfreq and self.include(entry)):
                    done+=1
                    outline=entry
                    #print "Filtering entry for "+entry
                    while len(features)>0:
                        freq=features.pop()
                        #feat=features.pop().lower()
                        feat=convert(features.pop(),delims=self.pathdelims)
                        feattot=float(coltotals.get(feat,0))
                        #print feat+"\t"+str(feattot-self.filterfreq)

                        if feattot>self.cfilterfreq:
                            outline+="\t"+feat+"\t"+freq
                            nofeats+=1

                    if nofeats>0:
                        outstream.write(outline+"\n")
                else:
                    pass
                    #print "Ignoring "+entry+" with frequency "+str(entrytot)

        outstream.close()

    #----
    #NORMALISE
    #this option normalises vectors so that they "sum to 1"
    #however they won't actually sum to 1 as row totals are computed before filtering
    #so necessary to call maketotals again after normalise
    #-----
    def normalise(self):
        rowtotals=self.load_rowtotals()
        infile= self.selectpos()+self.reducedstring+".filtered"
        outfile=infile+".norm"

        print "Normalising counts => sum to 1"
        outstream=open(outfile,"w")

        todo=len(rowtotals.keys())
        print "Estimated total vectors to do = "+str(todo)
        with open(infile) as instream:
            for lines,line in enumerate(instream):

                line = line.rstrip()
                fields=line.split("\t")
                entry=fields[0]
                features=fields[1:]
                entrytot=rowtotals[entry]
                outline=entry
                while len(features)>0:
                    weight=float(features.pop())
                    feat=convert(features.pop(),delims=self.pathdelims)
                    weight = weight/entrytot
                    outline+="\t"+feat+"\t"+str(weight)
                outline+="\n"
                outstream.write(outline)
                if lines%1000==0:
                    percent=lines*100.0/todo
                    print "Completed "+str(lines)+" vectors ("+str(percent)+"%)"
        outstream.close()
        self.normalised=True

    #---
    #load in pre-filtered and (optionally) normalised vectors
    #----
    def load_vectors(self,inpath=""):
        if inpath=="":
            infiles=[self.selectpos(p)+self.reducedstring+".filtered" for p in self.inpaths]
            if self.normalised and not self.option=="normalise":
                infiles=[p+".norm" for p in infiles]
        else:
            infiles=[inpath]


        vecs={}
        for infile in infiles:
            print "Loading vectors from: "+infile
            print "Words of interest: ",self.words
            with open(infile) as instream:

                for lines,line in enumerate(instream):
                    if lines%1000==0: print "Reading line "+str(lines)

                    line=line.rstrip()
                    fields=line.split("\t")
                    entry=fields[0]
                    #print entry
                    if self.include(entry):
                        vector={}
                        features=fields[1:]

                        index=0
                        while len(features)>0:
                            index+=1

                            freq=features.pop()
                            feat=convert(features.pop(),delims=self.pathdelims)
                            #print feat

                            #print str(index)+"\t"+feat+"\t"+str(freq)
                            try:
                                freq=float(freq)
                                vector[feat]=freq
                            except ValueError:
                                print "Error: "+str(index)+"\t"+feat+"\t"+str(freq)+"\n"
                                features=features+list(feat)
                        if entry in vecs.keys():
                            vecs[entry]=self.add(vecs[entry],vector)
                        else:
                            vecs[entry]=vector

            print "Loaded "+str(len(vecs.keys()))+" vectors"
        return vecs

    #----
    #write a set of vectors to file
    #these could be raw or PPMI vectors
    #----
    def output(self,vectors,outfile,append=False):
        #write a set of vectors to file
        print "Writing vectors to output file: "+outfile
        if append:
            outstream=open(outfile,"a")
        else:
            outstream=open(outfile,"w")
        for entry in vectors.keys():
            vector=vectors[entry]
            #print entry
            #print vector
            if len(vector.keys())>0:
                outstring=entry
                ignored=0
                nofeats=0
                for feat in vector.keys():
                    #print feat
                    forder=self.getorder(feat)

                    if forder>=self.minorder and forder<=self.maxorder:

                        try:
                            outstring+="\t"+feat+"\t"+str(vector[feat])
                            nofeats+=1
                        except:
                            ignored+=1
                #print "Ignored "+str(ignored)+" features"
                if nofeats>0:
                    outstream.write(outstring+"\n")

        outstream.close()



    #----SALIENCY FUNCTIONS

    #---
    #use the totals for each feature to compute grand totals for each feature type (e.g., "amod")
    #this is C<*,t,*> and is computed by S_f C<*,t,f>
    #----
    def compute_typetotals(self,feattots):
        #compute totals for different paths over all entries (using column totals given in feattots)
        print "Computing path totals C<*,t,*>"
        typetots={}
        for feature in feattots.keys():
            pathtype=self.getpathtype(feature)
            sofar=typetots.get(pathtype,0.0)
            typetots[pathtype]=sofar+float(feattots[feature])

        return typetots

    #--
    #compute path totals for each noun
    #i.e., C<w1,t,*>
    #this is needed in the standard PPMI calculation we use
    #----
    def compute_nounpathtotals(self,vectors):
        #compute totals for the different paths for each entry
        print "Computing path totals for each entry C<w1,t,*>"
        pathtotals={}
        for entry in vectors.keys():
            totalvector={}
            vector=vectors[entry]
            for feature in vector.keys():
                pathtype=self.getpathtype(feature)
                sofar=totalvector.get(pathtype,0.0)
                totalvector[pathtype]=sofar+float(vector[feature])

            pathtotals[entry]=totalvector
            #print entry
            #print pathtotals[entry]
            #print "nn:" + str(pathtotals[entry].get('nn',"not present"))
        return pathtotals

    #----
    #compute ppmi (or similar) for a set of vectors and return new set of vectors
    #@vecs: dict of dicts representing a set of vectors for which PPMI calculations to be carried out on i.e., vectors[w1][p:w2]=f  => C<w1,p,w2>=f
    #@pathtotals: dict of dicts representing a set of totals indexed by entry and path i.e., pathtots[w1][p] = f => C<w1,p,*> = f
    #@feattots: dict where feattots[p:w2]=f => C<*,p,w2>=f
    #@typetots: dict where typetots[p]=f => C<*,p,*>=f
    #@entrytots: dict where entrytots[w1]=f => C<w1,*,*>=f
    #
    #TODO: play with PPMI threshold and/or number of features
    #-----

    def computeppmi(self,vecs,pathtots,feattots,typetots,entrytots):

        ppmivecs={}
        grandtot=0.0
        if self.pp_normal:
            print "Computing pnppmi"
        elif self.gof_ppmi:
            print "Computing gof_ppmi"
            for type in typetots.keys():
                grandtot+=float(typetots[type])
            if self.smooth_ppmi:
                grandtot=math.pow(grandtot,0.75)

                #print type, grandtot
        else:
            print "Computing ppmi"
        done =0
        todo=len(vecs.keys())

        for entry in vecs.keys():
            #print entry
            ppmivector={}

            vector=vecs[entry]
            for feature in vector.keys():
                freq=float(vector[feature])  # C<w1,p,w2>
                try:
                    total=float(pathtots[entry][self.getpathtype(feature)]) # C<w1,p,*>
                except:
                    total=0.0001
                    print "Warning: no path total for %s: %s"%(feature,self.getpathtype(feature))
                feattot=float(feattots[feature]) #C<*,p,w2>
                typetot=float(typetots[self.getpathtype(feature)]) #C<*,p,*>
                entrytotal=float(entrytots[entry]) # C<w1,*,*>

              #  if self.smooth_ppmi:
              #      feattot=math.pow(feattot,0.75)
              #      typetot=math.pow(typetot,0.75)  # NO - not the same - need totals computed from smoothed values - need to smooth before computing type totals


                try:
                    if self.gof_ppmi:

                        pmi=math.log10((freq*grandtot)/(feattot*entrytotal))
                    else:
                        pmi=math.log10((freq*typetot)/(feattot*total))
                except:
                    pmi=0
                shifted_pmi=pmi-self.ppmithreshold
                if shifted_pmi>0:
                    if self.pp_normal:

                        shifted_pmi=shifted_pmi * total/entrytotal
                    ppmivector[feature]=shifted_pmi

            done+=1
            if done%1000==0:
                percent=done*100.0/todo
                print "Completed "+str(done)+" vectors ("+str(percent)+"%)"



            ppmivecs[entry]=self.mostsalient_vector(ppmivector)
            #print ppmivector
        return ppmivecs

    #-----
    #REVECTORISE
    #load the appropriate vectors and totals files, compute more totals, compute PPMI and output the returned vectors
    #------
    def revectorise(self):

        if self.normalised:
            suffix=".norm"
        else:
            suffix=""
        if self.pp_normal:
            suffix += ".pnppmi"
        elif self.gof_ppmi:
            suffix += ".gof_ppmi"
        elif self.smooth_ppmi:
            suffix += ".smooth_ppmi"
        else:
            suffix += ".ppmi"
        if self.ppmithreshold>0:
            suffix+="_"+str(self.ppmithreshold)
        if self.saliency>0:
            if self.saliencyperpath:
                suffix+=".spp_"+str(self.saliency)
            else:
                suffix+=".sal_"+str(self.saliency)
        outfile=self.selectpos()+self.reducedstring+".filtered"+suffix
        self.vecsbypos[self.pos]=self.load_vectors()
        self.feattotsbypos[self.pos]=self.load_coltotals(cds=self.smooth_ppmi)
        self.totsbypos[self.pos]=self.load_rowtotals()
        self.pathtotsbypos[self.pos]=self.compute_nounpathtotals(self.vecsbypos[self.pos])
        self.typetotsbypos[self.pos]=self.compute_typetotals(self.feattotsbypos[self.pos])

        ppmivecs=self.computeppmi(self.vecsbypos[self.pos],self.pathtotsbypos[self.pos],self.feattotsbypos[self.pos],self.typetotsbypos[self.pos],self.totsbypos[self.pos])
        self.output(ppmivecs,outfile)


    #---
    #use POS to determine which vectors/totals to supply to self.mostsalientvecs
    #----
    def mostsalient(self):
        return self.mostsalientvecs(self.vecsbypos[self.pos],self.pathtotsbypos[self.pos],self.feattotsbypos[self.pos],self.typetotsbypos[self.pos],self.totsbypos[self.pos])

    #---
    #compute PPMI and then only retain the most salient features (up to featmax for each includedtype)
    #does not modify ppmivectors
    #primary purpose has been to compute complete vectors to output to file but display the most salient ones for inspection
    #-----
    def mostsalientvecs(self,vecs,pathtots,feattots,typetots,entrytots):

        ppmivecs=self.computeppmi(vecs,pathtots,feattots,typetots,entrytots)
        if self.display:
            for entry in ppmivecs.keys():
                print "Most salient features for "+entry+" , width: "+str(len(vecs[entry].keys()))+", "+str(len(ppmivecs[entry].keys()))
                vector=ppmivecs[entry]
                #print vector
                feats=sorted(vector.items(),key=itemgetter(1),reverse=True)

                donetypes={}

                for tuple in feats:
                    feature=tuple[0]
                    pathtype=self.getpathtype(feature)
                    done=donetypes.get(pathtype,0)
                    if done<Composition.featmax and self.typeinclude(pathtype):
                        print feature+" : "+str(tuple[1])+" ("+str(vecs[entry][feature])+")"
                    donetypes[pathtype]=done+1

                print donetypes
                print "-----"
        return ppmivecs

    #-----
    #take a vector and retain only the most highly weighted features
    #-----
    def mostsalient_vector(self,ppmivector):

        if self.saliency>0:
            newvector={}
            feats=sorted(ppmivector.items(),key=itemgetter(1),reverse=True)
            donetypes={}
            all=0
            for tuple in feats:
                feature=tuple[0]
                pathtype=self.getpathtype(feature)
                done=donetypes.get(pathtype,0)
                if self.typeinclude(pathtype) and ((self.saliencyperpath and done<self.saliency)or(not self.saliencyperpath and all<self.saliency)):
                    newvector[feature]=tuple[1]
                    donetypes[pathtype]=done+1
                    all+=1
            return newvector
        else:
            return ppmivector
    #----
    #INSPECT
    #display the path distribution graph for a set of noun vectors and the most salient feature for those vectors
    #-----
    def inspect(self):
        self.pos="N"
        self.set_words()
        self.feattotsbypos[self.pos]=self.load_coltotals(cds=self.smooth_ppmi)
        self.totsbypos[self.pos]=self.load_rowtotals()
        self.vecsbypos[self.pos]= self.load_vectors()
        self.pathtotsbypos[self.pos]=self.compute_nounpathtotals(self.vecsbypos[self.pos])
        self.typetotsbypos[self.pos]=self.compute_typetotals(self.feattotsbypos[self.pos])
        print self.typetotsbypos[self.pos]
        graphing.display_bargraph(self.typetotsbypos[self.pos],title="Path Distribution over all Nouns")
        for entry in self.vecsbypos[self.pos].keys():
            title="Path Distribution for "+entry
            graphing.display_bargraph(self.pathtotsbypos[self.pos][entry],title)

        self.mostsalient()

    #----COMPOSITION FUNCTIONS

    #----
    #COMPOSE
    #load appropriate vectors, display most salient features for each vector, then runANcomposition and output to file
    #----
    def getComposedFilename(self,parampair):

        if self.normalised:
            suffix=".norm"
        else:
            suffix=""
        if self.pp_normal:
            suffix += ".pnppmi"
        elif self.gof_ppmi:
            suffix += ".gof_ppmi"
        elif self.smooth_ppmi:
            suffix += ".smooth_ppmi"
        else:
            suffix += ".ppmi"
        if self.ppmithreshold>0:
            suffix+="_"+str(self.ppmithreshold)
        if self.saliency>0:
            if self.saliencyperpath:
                suffix+=".spp_"+str(self.saliency)
            else:
                suffix+=".sal_"+str(self.saliency)

        key=str(parampair[0])+"-"+str(parampair[1])
        if parampair[0]=="hp":
            key1="_offsetting-"+str(self.offsetting)
        else:
            key1="_hp-"+str(self.headp)
        key2="_"+str(self.compop)
        outfile=self.selectpos()+self.reducedstring+".composed_"+key+key1+key2+suffix

        return outfile

    def compose(self,parampair=('','')):

        self.outfile=self.getComposedFilename(parampair)

        for pos in ["N","J","V"]:
            self.pos=pos
            self.set_words()
            self.feattotsbypos[pos]=self.load_coltotals(cds=self.smooth_ppmi)
            self.totsbypos[pos]=self.load_rowtotals()
            self.vecsbypos[pos]= self.load_vectors()
            self.pathtotsbypos[pos]=self.compute_nounpathtotals(self.vecsbypos[pos])
            self.typetotsbypos[pos]=self.compute_typetotals(self.feattotsbypos[pos])

            #self.mostsalient()


        self.output(self.runANcomposition(parampair=parampair),self.outfile)



    #----
    #run ANcompose for each adjective, noun pair of interest
    #then run mostsalientvecs on ANvecs which cause PPMI to be computed, most salient features displayed and PPMI vectors returned for output
    #----
    def runANcomposition(self):
        #this overridden in nouncompounds.py

        self.ANfeattots=self.addAN(self.feattotsbypos["J"],self.feattotsbypos["N"])  #C<*,t,f>
        self.ANtypetots=self.addAN(self.typetotsbypos["J"],self.typetotsbypos["N"])  #C<*,t,*>
       #print ANtypetots

        self.ANvecs={}
        self.ANtots={}
        self.ANpathtots={}

        if self.comppairfile:
            for comppair in self.comppairlist:
                 self.ANcompose(comppair[2],comppair[0])



        else:

            for adj in self.adjectives:
                for noun in self.nouns:
                    self.ANcompose(adj,noun)  #C<an,t,f>

                    #print ANvecs,ANtots


        return self.mostsalientvecs(self.ANvecs,self.ANpathtots,self.ANfeattots,self.ANtypetots,self.ANtots)

    #----
    #for a given adjective and noun, compute the compsed vector using addAN and the appropriate composed totals
    #add these to the dicts for ANs
    #----
    def ANcompose(self,adj,noun):
        self.CompoundCompose(adj,noun,"mod")

    def CompoundCompose(self,dep,head,rel,hp=0.5,compop="add",offsetting=1):
        hdpos=Composition.headPoS.get(rel,"N")
        dppos=Composition.depPoS.get(rel,"J")


        headvector=self.vecsbypos[hdpos][head]
        headpathtots=self.pathtotsbypos[hdpos][head]
        headtot=self.totsbypos[hdpos][head]

        depvector=self.vecsbypos[dppos][dep]
        deppathtots=self.pathtotsbypos[dppos][dep]
        deptot=self.totsbypos[dppos][dep]

        if self.distinguish:
            tag="composed:"
        else:
            tag=""

        entry=dep.split("/")[0]+"|"+tag+rel+"|"+head
        print "Composing vectors for "+entry
        self.ANvecs[entry]=self.doCompound(depvector,headvector,rel,hp=hp,op=compop,offsetting=offsetting)
        print "Composing path totals"
        self.ANpathtots[entry]=self.doCompound(deppathtots,headpathtots,rel,hp=hp,op=compop,offsetting=offsetting)

        # print self.ANpathtots[entry]
        # print "nn: "+str(self.ANpathtots[entry].get('nn',"not present"))
        self.ANtots[entry]=float(deptot)+float(headtot)


    #----
    #handle composition operations such as offsetting (required or not) and merge operation (add, min etc)
    #-----
    def doCompound(self,depvector,headvector,rel,nntest=False,hp=0.5,op="add",offsetting=1):
        if self.reversed:
            temp =depvector
            depvector=headvector
            headvector=temp

        if offsetting<=0:
            #print "using depvector"
            offsetvector=dict(depvector)
        elif offsetting>=1:
            #print "using offsetvector"
            offsetvector=self.offsetVector(depvector,Composition.basicRel[rel])
        else:
            #print "using mixture",offsetting
            offsetvector=self.add(self.offsetVector(depvector,Composition.basicRel[rel]),depvector,weight=offsetting)


        if nntest:
            print offsetvector
            print "Offset vector nn: "+str(offsetvector.get('nn',"not present"))
            print "Offset vector '': "+str(offsetvector.get('',"not present"))
            print headvector
            print "Head vector '': "+str(headvector.get('',"not present"))

        if op=="add":
            return self.addCompound(offsetvector,headvector,hp)
        elif op=="min":
            return self.minCompound(offsetvector,headvector)
        elif op=="smooth_mult":
            return self.smooth_multCompound(offsetvector,headvector)
        else:
            print "Error unknown composition operation: "+op
            exit(-1)

    #----
    #add an adjective vector to a noun vector (may be feature vectors or path vectors)
    #do this by offsetting the adjective vector so that it is aligned with the noun vector
    #then add noun vector to adjective vector
    #----
    def addCompound(self,offsetvector,headvector,hp=0.5):

        COMPOUNDvector={}
        #print "Processing noun features "+str(len(nounvector.keys()))
        count=0
        intersect=[]
        #print nounvector
        #print adjvector
        for feature in headvector.keys():
            count+=1
            if feature in offsetvector:
                COMPOUNDvector[feature]=hp*float(headvector[feature])+(1-hp)*float(offsetvector[feature])
                intersect.append(feature)
                offsetvector.__delitem__(feature)
            else:
                COMPOUNDvector[feature]=hp*headvector[feature]
            #if count%10000==0:print"Processed "+str(count)

        print "Intersecting features: "+str(len(intersect))
        #print "Processing remaining adj features "+str(len(adjvector.keys()))+" : reduced to : "+str(len(offsetvector.keys()))
        for feature in offsetvector.keys():
            COMPOUNDvector[feature]=(1-hp)*offsetvector[feature]
        #print "Complete"
        return COMPOUNDvector

    def minCompound(self,offsetvector,headvector):
        COMPOUNDvector={}
        intersect=[]
        for count,feature in enumerate(headvector.keys()):
            if feature in offsetvector:
                COMPOUNDvector[feature]=min(float(headvector[feature]),float(offsetvector[feature]))
                intersect.append(feature)

        print "Intersecting features: "+str(len(intersect))
        return COMPOUNDvector

    def smooth_multCompound(self,offsetvector,headvector):

        COMPOUNDvector={}
        intersect=[]
        for count,feature in enumerate(headvector.keys()):
            if feature in offsetvector:
                COMPOUNDvector[feature]=(float(headvector[feature])+1)*(float(offsetvector[feature])+1)
                intersect.append(feature)
                offsetvector.__delitem__(feature)
            else:
                COMPOUNDvector[feature]=float(headvector[feature])+1

        for feature in offsetvector.keys():
            COMPOUNDvector[feature]=float(offsetvector[feature])+1

        print "Intersecting features: "+str(len(intersect))

        return COMPOUNDvector


    def addAN(self,adjvector,nounvector):
        return self.doCompound(adjvector,nounvector,"mod")



    #----
    #offset an adjective vector so that it aligns with the noun vector it is modifying
    #----
    def offsetAN(self,adjvector):
        return self.offsetVector(adjvector,"mod")

    def offsetVector(self,depvector,rel):
        depPREFIX="_"+rel
        headPREFIX=rel

        offsetvector={}
        incomp=0
        for feature in depvector.keys():
            (prefix,suffix)= self.splitfeature(feature)
            #print depPREFIX,prefix,suffix
            if prefix==depPREFIX:
                newfeature=suffix+self.getpathvalue(feature)

            elif prefix.startswith("_"):
                #incompatible feature for composition
                #print "Incompatible feature for composition: "+feature
                incomp+=1
                newfeature="ignore"
            #elif feature.startswith(":"):
            elif prefix=="":
                newfeature=headPREFIX+feature
            else:
                newfeature=headPREFIX+self.pathdelims[0]+feature
            #print newfeature
            if not newfeature == "ignore":
                offsetvector[newfeature]=depvector[feature]
        #print "Features in original adj vector: "+str(len(adjvector.keys()))
        #print "Incompatible features in adjective vector: "+str(incomp)
        #print "Features in offset adj vector: "+str(len(offsetvector.keys()))
        return offsetvector

    #---
    #add two vectors
    #not used I think
    #---
    def add(self,avector,bvector,weight=0.5):
        rvector={}
        for feature in bvector.keys():
            if feature in avector:
                rvector[feature]=weight*float(avector[feature])+(1-weight)*float(bvector[feature])
                avector.__delitem__(feature)
            else:
                rvector[feature]=(1-weight)*float(bvector[feature])

        for feature in avector.keys():
            rvector[feature]=weight*float(avector[feature])

        return rvector
        #rvector=dict(avector)
        #for feat in bvector.keys():
        #    rvector[feat]=rvector.get(feat,0)+bvector[feat]
        #return rvector


    #----OTHER FUNCTIONS

    def intersect(self):

        self.nounfeattots=self.load_coltotals()
        self.nountots=self.load_rowtotals()
        self.nountypetots=self.compute_typetotals(self.nounfeattots)
        self.nounvecs= self.load_vectors()
       # self.nounpathtots=self.compute_nounpathtotals(self.nounvecs)

        intersectedvecs=self.intersectall()
        self.nounpathtots=self.compute_nounpathtotals(intersectedvecs)
        self.mostsalientvecs(intersectedvecs,self.nounpathtots,self.nounfeattots,self.nountypetots,self.nountots)

    def intersectall(self):

        intersected={}
        for wordlist in self.wordlistlist:
            name=self.join(wordlist,'_')
            vector=self.nounvecs[wordlist[0]]
            for aword in wordlist[1:]:
                vector = self.intersecteach(vector,self.nounvecs[aword])
            intersected[name]=vector
            total=0
            for value in vector.values():
                total+=value
            self.nountots[name]=total
        return intersected

    def intersecteach(self,avector,bvector):
        newvector={}
        for feat in avector.keys():
            value=min(avector[feat],bvector.get(feat,0))
            if value>0:
                newvector[feat]=value
        return newvector




    def rewrite(self):
        self.output(self.load_vectors(self.inpath),self.inpath+".new")

    #----main run function
    def run(self,words=[]):

         #if present load phrases for composition
        # and set words/paths of interest

        if self.comppairfile!="":
            with open(self.comppairfile) as fp:
                self.comppairlist = yaml.safe_load(fp)
        else:
            self.comppairlist=[]
        self.set_words(words=words)

        while len(self.options)>0:
            self.option=self.options[0]
            self.options=self.options[1:]

            print "Stage: "+self.option
            if self.option=="split":
                self.splitpos()
            elif self.option=="reduceorder":
                self.reduceorder()
            elif self.option=="maketotals":
                self.maketotals()
            elif self.option =="filter":
                self.filter()
            elif self.option == "normalise":
                self.normalise()
            elif self.option=="compose":
                self.compose()
            elif self.option=="inspect":
                self.inspect()
            elif self.option=="revectorise":
                self.revectorise()
            elif self.option=="intersect":
                self.intersect()
            elif self.option=="rewrite":
                self.rewrite()


            else:
                print "Unknown option: "+self.option

    def close(self):
        #release memory by deleting stored vectors
        self.vecsbypos={}
        self.totsbypos={}
        self.feattotsbypos={}
        self.pathtotsbypos={}
        self.typetotsbypos={}
        self.ANfeattots={}
        self.ANpathtots={}
        self.ANtots={}
        self.ANvecs={}
        self.ANtypetots={}


if __name__=="__main__":

    #----
    #example runs:
    #python composition.py split filename
    #python composition.py reduceorder filename N 1 2
    #python composition.py maketotals filename N 1 2
    #python composition.py filter filename N 1 2
    #python composition.py normalise filename N 1 2
    #python composition.py revectorise filename N 1 2 ppmi
    #python composition.py compose filename AN 0 2 normalised pnppmi
    #-----

    myComposer = Composition(sys.argv[1:])
    myComposer.run()