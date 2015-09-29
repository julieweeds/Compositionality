__author__ = 'juliewe'

import sys
from compounds import Compounder, Compound


def getPos(word):
    return word.split('/')[1]

def getLex(word):
    return word.split('/')[0].lower()

class CompoundFinder(Compounder):

    def __init__(self,configfile):
        Compounder.__init__(self,configfile)
        self.corpusfile=self.configured.get("corpus")
        self.readcompounds()
        self.counts={}
        for head in self.compounds.keys():
            for comp in self.compounds[head]:
                self.counts[comp]=0


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

        mini = min(self.counts.values())
        maxi = max(self.counts.values())
        print "Min count is "+str(mini)
        print "Max count is "+str(maxi)
        print self.counts.keys()

        minicounts=0
        for comp in self.counts.keys():
            if self.counts[comp]<mini+5:
                minicounts+=1

                for line in comp.getConll():
                    print line
                print self.counts[comp]
            elif self.counts[comp]>maxi-5:

                for line in comp.getConll():
                    print line
                print self.counts[comp]

        print "Low counts: "+str(minicounts)
        #for key in self.counts.keys():
         #   for line in key.getConll():
          #      print line
           # print self.counts[key]

    def process_sentence(self,sentence):

        rel_to_match=self.rel[0]  #assuming there is only one rel in list at the moment
        dp=str(self.rel[1])
        hp=str(self.rel[2])

        for arc in sentence.values():
            if arc[2]==rel_to_match:

                dep = arc[0]
                head=sentence[arc[1]][0]
                if getPos(dep).startswith(self.posmap[dp]) and getPos(head).startswith(self.posmap[hp]):

                    dep=getLex(dep)
                    head=getLex(head)

                    for comp in self.compounds.get(head,[]):

                        if comp.match(head,dep):
                            sofar=self.counts[comp]

                            self.counts[comp]=sofar+1
        #print self.compounds.keys()


    def run(self):

        self.process_corpus()


if __name__=="__main__":
    myCompoundFinder=CompoundFinder(sys.argv[1])
    myCompoundFinder.run()