__author__ = 'juliewe'
#take in the list of noun compounds and generate other formats for composition etc

import sys, ConfigParser, numpy as np, scipy.stats as stats
try:
    import graphing
    graphing_loaded=True
except:
    print "Warning: unable to load graphing module"
    graphing_loaded=False


def getPos(word):
    try:
        return word.split('/')[-1][0]
    except:
        return word.split('/')[1]

def getLex(word,delim='/',postagged=True,tdelim=' '):
    if postagged:
        tokens=word.split(tdelim)
        lex=tokens[0].split(delim)[0].lower()
        if len(tokens)>1:
            for token in tokens[1:]:
                lex += tdelim+token.split(delim)[0].lower()
        return lex
    else:
        return word

def matchmake(phrase):
    tokens=phrase.split('|')
    return tokens[0]+" "+tokens[2],tokens[1]



class Compound:

    #root=("root",-1,0)
    #pos=("?",-1)

    def __init__(self,name,ctype="farahmand"):
        self.name=name
        self.judgements=[]
        self.rels=[]
        #self.PoSes=[]
        self.counts=[]
        self.ctype=ctype
        self.autosims={}
        self.frequency=-1

    def addJudgement(self,NC):
        self.judgements.append(NC)

    def setJudgements(self,NCs):
        self.judgements=[float(item) for item in NCs]

    def setFreq(self,freq):
        self.freq=freq

    def getFreq(self):
        return np.log(self.freq)

    def addAutoSim(self,key,sim):
        self.autosims[key]=float(sim)

    def addRel(self,rellist):
        found =False
        for index,rel in enumerate(self.rels):
            if rellist==rel:
                self.counts[index]+=1
                found=True
                break
        if not found:
            self.rels.append(rellist)
            self.counts.append(1)

    def getFirst(self):
        return self.name.split(" ")[0]

    def getSecond(self):
        return self.name.split(" ")[1]

    def match(self,dep,head, rel):

        if (head==self.getSecond() and dep ==self.getFirst()):
            self.addRel([rel[0],1,2,rel[1],rel[2]])
            return True
        # elif (head==self.getFirst() and dep ==self.getSecond()):
        #     self.addRel([rel[0],2,1,rel[2],rel[1]])
        #     return True
        else:
            return False

    def filter(self,freqthresh,blacklist):

        newrels=[]
        newcounts=[]
        for index,count in enumerate(self.counts):
            if count>=freqthresh and self.rels[index] not in blacklist:
                newrels.append(self.rels[index])
                newcounts.append(count)

        if len(newcounts)==0:
            self.filter(0,blacklist)
            if len(self.counts)>0:
                self.pickBest()
        else:
            self.counts=newcounts
            self.rels=newrels
        return len(newcounts)>1

    def pickBest(self):
        max=0
        mindex=-1
        for index,count in enumerate(self.counts):
            if count>max:
                max=count
                mindex=index
        self.counts=[self.counts[mindex]]
        self.rels=[self.rels[mindex]]

    def getScore(self):
        if self.ctype=="reddy":
            return self.judgements[2]
        else:
            return sum(self.judgements)

    def getAutoSim(self):
        res=""
        for key in self.autosims.keys():
            res+=key+":"+str(self.autosims[key])+"\t"

        return res

    def getSimScore(self):
        if len(self.autosims.values())>0:
            return np.max(self.autosims.values())
        else:
            return -1

    def getTotal(self):
        return sum(self.counts)

    def getConll(self):
        tokens=self.name.split(" ")
        lines=[]

        rels_to_use=self.rels
        counts_to_use=self.counts
        if len(rels_to_use)==0:
            rels_to_use=[("?",1,2,"N","N")]
            counts_to_use=[0]

        for index,rel in enumerate(rels_to_use):
            for(id,token) in enumerate(tokens):
                sid=id+1

                if sid==rel[1]: #dependent
                    line=str(sid)+"\t"+token+"/"+rel[3]+"\t"+str(rel[2])+"\t"+rel[0]
                elif sid==rel[2]: # head
                    line=str(sid)+"\t"+token+"/"+rel[4]+"\t0\troot"
                else:
                    print "Too many tokens"
                lines.append(line)
            lines.append("counts: "+str(counts_to_use[index]))

        return lines

    def generate(self,rel_list):
        self.compounds={}
        for rel in rel_list:
            acompound=self.getFirst()+"|"+rel+"|"+self.getSecond()+"/N"
            self.compounds[rel]=acompound

    def display(self):
        for line in self.getConll():
            print line
        print "Non-Compositionality",self.getScore()
        print "Automatic Similarity",self.getAutoSim()
        print "-----"

    def display_compounds(self):
        for value in self.compounds.values():
            print value

    def write_compounds_to_file(self,outpath):
        for value in self.compounds.values():
            outpath.write(value+"\n")

class Compounder:

    def __init__(self,filename):
        self.configfile=filename
        #with open(self.configfile) as fp:
        #    self.configured=yaml.safe_load(fp)
        self.config=ConfigParser.RawConfigParser()
        self.config.read(self.configfile)


#        print self.configured
        self.compounds={} #key is name, value is list of compounds
        self.firstindex={} #key is first word, value is name
        self.secondindex={} #key is second word, value is name
        self.ctype=self.config.get('default','ctype')
        try:
            self.compound_file=self.config.get('default','compound_file')

        except:
            print "Error: problem with configuration"

    def readcompounds(self):
        if self.ctype=="reddy":
            self.readreddy()
        else:
            with open(self.compound_file) as fp:
                for line in fp:
                    line=line.rstrip()
                    fields=line.split(",")
                    aCompound=Compound(fields[0])
                    aCompound.setJudgements(fields[1:])

                    if self.compounds.get(aCompound.name,None) == None:
                        self.compounds[aCompound.name]=aCompound
                        sofar= self.firstindex.get(aCompound.getFirst(),[])
                        sofar.append(aCompound.name)
                        self.firstindex[aCompound.getFirst()]=sofar
                        sofar=self.secondindex.get(aCompound.getSecond(),[])
                        sofar.append(aCompound.name)
                        self.secondindex[aCompound.getSecond()]=sofar
                    else:
                        print "Error: Duplicate noun compound"


    def readreddy(self):

        with open(self.compound_file) as fp:
            for line in fp:
                line=line.rstrip()
                if not line.startswith("#"):
                    fields=line.split('\t')
                    text=fields[0]

                    #print text, getLex(text,'-')
                    aCompound=Compound(getLex(text,'-'),ctype=self.ctype)
                    scores=fields[1].split(' ')
                    aCompound.setJudgements([scores[0],scores[2],scores[4]])

                    if self.compounds.get(aCompound.name,None) == None:
                        self.compounds[aCompound.name]=aCompound
                        sofar= self.firstindex.get(aCompound.getFirst(),[])
                        sofar.append(aCompound.name)
                        self.firstindex[aCompound.getFirst()]=sofar
                        sofar=self.secondindex.get(aCompound.getSecond(),[])
                        sofar.append(aCompound.name)
                        self.secondindex[aCompound.getSecond()]=sofar
                    else:
                        print "Error: Duplicate noun compound"


    def generate(self,rel_list,outfile=""):
        self.readcompounds()
        self.generated_compounds=[]
        if outfile != "":
            outpath=open(outfile,"w")

        for compound in self.compounds.values():
            compound.generate(rel_list)
            self.generated_compounds+=compound.compounds.values()
            if outfile=="":
                compound.display_compounds()
            else:
                compound.write_compounds_to_file(outpath)

        if outfile != "":
            outpath.close()

    def addAutoSim(self,compound,sim):
        phrase,type = matchmake(compound.split('/')[0])
        #print phrase
        if phrase in self.compounds.keys():
            self.compounds[phrase].addAutoSim(type,sim)

    def addFreq(self,compound,freq):
        phrase,type = matchmake(compound.split('/')[0])
        if phrase in self.compounds.keys():
            self.compounds[phrase].setFreq(freq)

    def correlate(self):
        listX=[]
        listY=[]
        listZ=[]
        for compound in self.compounds.values():
            compound.display()
            if compound.getSimScore()>-1:
                listX.append(compound.getScore())
                listY.append(compound.getSimScore())
                listZ.append(compound.getFreq())

        #print listX
        #print listY

        print "Mean Human Judgement score: ",np.mean(listX)
        print "Mean AutoSim score: ",np.mean(listY)
        print "Mean Frequency: ",np.mean(listZ)
        print "Spearman's Correlation Coefficient and p'value for Human Judgements vs Automatic Similarity over %s values: "%(str(len(listX))),stats.spearmanr(np.array(listX),np.array(listY))
        if graphing_loaded:
            graphing.makescatter(listX,listY)
        print "Spearman's Correlation Coefficient and p'value for Human Judgments vs Frequency: ",stats.spearmanr(np.array(listX),np.array(listZ))
        #if graphing_loaded:
        #    graphing.makescatter(listX,listZ)
        print "Spearman's Correlation Coefficient and p'value for Automatic Similarity vs Frequency: ", stats.spearmanr(np.array(listY),np.array(listZ))
        #if graphing_loaded:
        #    graphing.makescatter(listY,listZ)

    def run(self):

        self.readcompounds()
        print self.compounds.keys()
        for compound in self.compounds.values():
                compound.display()
        return

if __name__=="__main__":
    #first argument should be name of config file
    myCompounder=Compounder(sys.argv[1])
    myCompounder.run()
