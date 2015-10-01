__author__ = 'juliewe'
#take in the list of noun compounds and generate other formats for composition etc

import sys, yaml

def getPos(word):
    try:
        return word.split('/')[-1][0]
    except:
        return word.split('/')[1]

def getLex(word,delim='/',postagged=True):
    if postagged:
        tokens=word.split(' ')
        lex=tokens[0].split(delim)[0].lower()
        if len(tokens)>1:
            for token in tokens[1:]:
                lex += " "+token.split(delim)[0].lower()
        return lex
    else:
        return word

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

    def addJudgement(self,NC):
        self.judgements.append(NC)

    def setJudgements(self,NCs):
        self.judgements=[float(item) for item in NCs]

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

    def display(self):
        for line in self.getConll():
            print line
        print "Non-Compositionality",self.getScore()
        print "-----"


class Compounder:

    def __init__(self,filename):
        self.configfile=filename
        with open(self.configfile) as fp:
            self.configured=yaml.safe_load(fp)
        print self.configured
        self.compounds={} #key is name, value is list of compounds
        self.firstindex={} #key is first word, value is name
        self.secondindex={} #key is second word, value is name
        self.ctype=self.configured.get("ctype","farahmand")
        try:
            self.compound_file=self.configured["compound_file"]

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
