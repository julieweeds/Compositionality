__author__ = 'juliewe'

import sys
from compounds import Compounder, Compound, getLex, getPos

class CompoundFinder(Compounder):

    freqthresh=10
    blacklist=["conj","dobj","nsubj","poss","appos"]
    def __init__(self,configfile):
        Compounder.__init__(self,configfile)
        self.corpusfile=self.configured.get("corpus")
        self.contiguous=self.configured.get("contiguous",False)
        self.ptype=self.configured.get("ptype","wiki")
        self.readcompounds()
        self.counts={}
        self.countpos={}
        self.rels={}

        if self.ptype=="nyt":
            self.lex=0
            self.headpos=3
            self.relname=4
            self.pos=2
            self.postagged=False
            self.arclength=5
        else:
            self.headpos=1
            self.relname=2
            self.lex=0
            self.postagged=True
            self.arclength=3

        for comp in self.compounds.keys():
            self.counts[comp]=0
            self.countpos[comp]=0

    def process_corpus(self):
        print self.compounds.keys()
        with open(self.corpusfile) as fp:
            sentencebuffer={}
            lines=0
            for line in fp:
                #need to load in sentence at a time and process each sentence
                #need to be careful in case there are any dependencies on dependent
                line=line.rstrip()
                fields=line.split('\t')
                if len(fields)==4 or len(fields)==6:
                    sentencebuffer[fields[0]]=fields[1:]
                else:
                    if self.contiguous:
                        self.process_sentence_contiguous(sentencebuffer)
                    else:
                        self.process_sentence(sentencebuffer)
                    sentencebuffer={}
                lines+=1
                if lines%1000000==0:
                    print "Processed "+str(lines)+" lines"
                    print "Found "+str(self.non_zero())+" out of "+str(len(self.counts.keys()))+" compounds"
        #analyse counts
        if self.contiguous:
            self.analyse_contiguous()
        else:
            self.analyse()

    def non_zero(self):
        sum=0
        for v in self.counts.values():
            if v>0:sum+=1
        return sum

    def analyse_contiguous(self):
        sum=0
        sumpos=0
        for candkey in self.counts.keys():
            if self.counts[candkey]>0:
                sum+=1
                if self.countpos[candkey]>0:sumpos+=1
                else: print "Not found with correct PoS: "+candkey
            else: print "Not found: "+candkey


        print "Found "+str(sum)+" contiguous compounds"
        print "Found "+str(sumpos)+" with correct PoS tags"
        print "Min count is "+str(min(self.counts.values()))
        print "Max count is "+str(max(self.counts.values()))
        print "Relation set observed ",self.rels

    def analyse(self):
        sum=0
        multiple=0
        counts=[]
        rels=[]
        freqmultiple=0
        for comp in self.compounds.values():
            if len(comp.counts)>0:
                sum+=1

                if len(comp.counts)>1:
                    multiple+=1
                    if comp.filter(CompoundFinder.freqthresh,CompoundFinder.blacklist):
                        freqmultiple+=1
                        comp.display()
                        #comp.pickBest()
                #comp.display()
                counts.append(comp.getTotal())
                rel =comp.rels[0][0]
                if rel not in rels: rels.append(rel)
        print "Found "+str(sum)+" compounds"
        print "Multiple analyses: "+str(multiple)
        print "Multiple analyses passing filters: "+str(freqmultiple)
        print "Relations included: ",rels
        mini = min(counts)
        maxi = max(counts)
        print "Min count is "+str(mini)
        print "Max count is "+str(maxi)

    def process_sentence_contiguous(self,sentence):
        #just interested in contiguous nouns don't care about relations

        #print sentence
        for i in sentence.keys():
            try:
                sid =int (i)
                if sid<len(sentence.values())-1:
                    arc=sentence[i]
                    canddep=getLex(arc[self.lex]).lower()
                    candhead=getLex(sentence[str(sid+1)][self.lex]).lower()

                    candkey=canddep+" "+candhead
                    #print candkey
                    if candkey in self.compounds.keys():
                        self.counts[candkey]+=1
                        if self.postagged:
                            deppos=getPos(arc[self.lex])
                            dephead=getPos(sentence[str(sid+1)][self.lex])
                        else:
                            deppos=arc[self.pos]
                            dephead=sentence[str(sid+1)][self.pos]

                        if deppos=="N" and dephead=="N":
                            self.countpos[candkey]+=1
                        if arc[self.headpos]==str(sid+1): #dependency relationship
                            sofar=self.rels.get(arc[self.relname],0)
                            self.rels[arc[self.relname]]=sofar+1
                        else:
                            print candkey, " :contiguous but no dependency: ",i,arc,sentence[str(sid+1)],sentence[arc[self.headpos]]
            except:
                pass

    def process_sentence(self,sentence):
        #do not know relation name or Pos tags
        #looking for any combination of the compound words
        #print sentence
        for i,arc in enumerate(sentence.values()):
            if len(arc)==self.arclength:
                canddep=getLex(arc[self.lex],postagged=self.postagged).lower()
                #print self.relname,arc
                candrel=arc[self.relname]
                #print arc
                if int(arc[self.headpos])==i+2 and candrel not in CompoundFinder.blacklist:
                    candhead= getLex(sentence[arc[self.headpos]][0]).lower()

                    candkey = canddep+" "+candhead
                    #print candkey
                    if candkey in self.compounds.keys():
                        if self.ptype=="nyt":
                            self.compounds[candkey].match(canddep,candhead,(candrel,arc[self.pos][0],sentence[arc[self.headpos]][self.pos][0]))
                        else:
                            self.compounds[candkey].match(canddep,candhead,(candrel,getPos(arc[0]),getPos(sentence[arc[1]][0])))
                        self.counts[candkey]+=1
                            #print "Match found for: ", arc
                            #self.compounds[candkey].display()




    def run(self):

        self.process_corpus()


if __name__=="__main__":
    myCompoundFinder=CompoundFinder(sys.argv[1])
    myCompoundFinder.run()