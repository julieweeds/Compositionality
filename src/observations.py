__author__ = 'juliewe'

import sys
from compounds import Compounder, Compound

class CompoundFinder(Compounder):

    def __init__(self,configfile):
        Compounder.__init__(self,configfile)
        self.corpusfile=self.configured.get("corpus_file")

        frequencies={}
        for comp in self.compounds:
            frequencies[comp]=0

    def process_corpus(self):
        with open(self.corpusfile) as fp:
            for line in fp:
                #need to load in sentence at a time and process each sentence
                #need to be careful in case there are any dependencies on dependent
                line=line.rstrip()
        return

    def run(self):
        self.process_corpus()


if __name__=="__main__":
    myCompoundFinder=CompoundFinder(sys.argv[1])
    myCompoundFinder.run()