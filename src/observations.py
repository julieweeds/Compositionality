__author__ = 'juliewe'

import sys
from compounds import Compounder, Compound


def getPos(word):
    try:
        return word.split('/')[-1][0]
    except:
        return word.split('/')[1]

def getLex(word):
    return word.split('/')[0].lower()

class CompoundFinder(Compounder):

    freqthresh=10
    blacklist=["conj","dobj","nsubj","poss","appos"]
    def __init__(self,configfile):
        Compounder.__init__(self,configfile)
        self.corpusfile=self.configured.get("corpus")
        self.readcompounds()
        self.counts={}

        # for comp in self.compounds.keys():
        #     self.counts[comp]=0


    def process_corpus(self):
        with open(self.corpusfile) as fp:
            sentencebuffer={}
            lines=0
            for line in fp:
                #need to load in sentence at a time and process each sentence
                #need to be careful in case there are any dependencies on dependent
                line=line.rstrip()
                fields=line.split('\t')
                if len(fields)==4:
                    sentencebuffer[fields[0]]=fields[1:]
                else:
                    self.process_sentence(sentencebuffer)
                    sentencebuffer={}
                lines+=1
                if lines%1000000==0: print "Processed "+str(lines)+" lines"
        #analyse counts

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


    def process_sentence(self,sentence):
        #do not know relation name or Pos tags
        #looking for any combination of the compound words
        for i,arc in enumerate(sentence.values()):
            canddep=getLex(arc[0]).lower()
            candrel=arc[2]
            #print arc
            if int(arc[1])==i+2 and candrel not in CompoundFinder.blacklist:
                candhead= getLex(sentence[arc[1]][0]).lower()

                candkey = canddep+" "+candhead

                if candkey in self.compounds.keys():
                    self.compounds[candkey].match(canddep,candhead,(candrel,getPos(arc[0]),getPos(sentence[arc[1]][0])))
                        #print "Match found for: ", arc
                        #self.compounds[candkey].display()




    def run(self):

        self.process_corpus()


if __name__=="__main__":
    myCompoundFinder=CompoundFinder(sys.argv[1])
    myCompoundFinder.run()