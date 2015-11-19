__author__ = 'juliewe'

import sys
import os
import ast

import numpy as np, gzip

from compounds import Compounder, getLex, getPos


class CompoundFinder(Compounder):

    freqthresh=10
    blacklist=["conj","dobj","nsubj","poss","appos"]

    postags=["N","V","J","R","F"]

    def __init__(self,configfile):
        Compounder.__init__(self,configfile)
        self.corpusfile=self.config.get('default','corpus')
        self.contiguous=(self.config.get('default','contiguous')=='True')
        self.ptype=self.config.get('default','ptype')
        self.convert=(self.config.get('default','convert')=='True')
        self.clean=(self.config.get('default','clean')=='True')
        self.dir=(self.config.get('default','dir')=='True')
        self.readcompounds()
        self.counts={}
        self.countpos={}
        self.countrel={}
        self.rels={}


        self.lex=int(self.config.get(self.ptype,'lex'))
        self.lemma=int(self.config.get(self.ptype,'lemma'))
        self.headpos=int(self.config.get(self.ptype,'headpos'))
        self.relname=int(self.config.get(self.ptype,'relname'))
        self.pos=int(self.config.get(self.ptype,'pos'))
        self.postagged=(self.config.get(self.ptype,'postagged')=='True')
        self.arclength=int(self.config.get(self.ptype,'arclength'))
        self.firstindex=int(self.config.get(self.ptype,'firstindex'))
        self.minlength=int(self.config.get(self.ptype,'minlength'))
        self.extra=ast.literal_eval(self.config.get(self.ptype,'extra'))
        self.erased=ast.literal_eval(self.config.get(self.ptype,'erased'))
        try:
            self.uselemma=(self.config.get(self.ptype,'uselemma')=='True')
        except:
            self.uselemma=False

        #if self.uselemma and self.lemma>-1:
        #    self.lex=self.lemma  #replace lexeme with lemma
        #    self.lemma=-1 #ignore lemma

        for comp in self.compounds.keys():
            self.counts[comp]=0
            self.countpos[comp]=0
            self.countrel[comp]=0

        self.cont_nodep=0
        self.cont=0
        self.cont_revdep=0
        self.skippedfiles=[]

    def process_file(self,file=""):
        if file=="":file=self.corpusfile

        if file.endswith(".gz"):
            fp = gzip.open(file)
            if self.convert:
                self.outstream=gzip.open(file[0:-3]+".compounds"+".gz","w")
            if self.clean:
                self.cleanstream=gzip.open(file[0:-3]+".clean"+".gz","w")

        else:
            try:
                fp = open(file)
            except:
                try:
                    fp = gzip.open(file+".gz")
                except:
                    self.skippedfiles.append(file)
            if self.convert:
                self.outstream = open(file+".compounds.gz","w")
            if self.clean:
                self.cleanstream = open(file+".clean.gz","w")


        sentencebuffer={}

        for lines,line in enumerate(fp):
            #need to load in sentence at a time and process each sentence
            #need to be careful in case there are any dependencies on dependent
            line=line.rstrip()
            fields=line.split('\t')
            if len(fields)==self.arclength+1:
                sentencebuffer[fields[0]]=fields[1:]
            else:
                if len(fields)==self.minlength+1:
                    arc=fields[1:]+self.extra
                    sentencebuffer[fields[0]]=arc
                else:
                    if self.contiguous:
                        self.process_sentence_contiguous(sentencebuffer)
                    else:
                        self.process_sentence(sentencebuffer)
                    sentencebuffer={}
            lines+=1
            if lines%1000000==0:
                print "Processed "+str(self.lines)+" lines"
                print "Found "+str(self.non_zero())+" out of "+str(len(self.counts.keys()))+" compounds"
                print "Found "+str(self.non_zero(list=self.countrel.values()))+" with dependency relationship"

        fp.close()

        if self.convert:
            self.outstream.close()
        if self.clean:
            self.cleanstream.close()

    def process_corpus(self):

        print self.compounds.keys()
        self.lines=0
        if self.dir:
            indir=self.corpusfile
            print "Processing ",indir
            os.chdir(indir)
            for datafile in [df for df in os.listdir(indir)]:
                print "Processing ", datafile
                if datafile.endswith(".clean") or datafile.endswith(".compounds") or datafile.endswith(".clean.gz") or datafile.endswith(".compounds.gz"):
                    pass
                else:
                    self.process_file(file=datafile)



        else: self.process_file()

        if len(self.skippedfiles)>0:
            todo=self.skippedfiles
            self.skippedfiles=[]
            print "Retrying skipped files: ",todo
            for datafile in self.todo:
                print "Processing ", datafile
                self.process_file(file=datafile)
            if len(self.skippedfiles)>0:
                print "Really skipped "+str(len(self.skippedfiles))+" files: ",self.skippedfiles

        #analyse counts
        if self.contiguous:
            self.analyse_contiguous()
        else:
            self.analyse()


    def non_zero(self,list=[]):
        if len(list)==0:list=self.counts.values()
        sum=0
        for v in list:
            if v>0:sum+=1
        return sum

    def analyse_contiguous(self):
        sum=0
        sumpos=0
        sumrel=0
        print "In total found "+str(self.cont)+" contiguous compounds"
        print "In total "+str(self.cont_revdep)+" with reverse dependency"
        print "In total "+str(self.cont_nodep)+" without dependency"
        for candkey in self.counts.keys():
            if self.counts[candkey]>0:
                sum+=1
                if self.countpos[candkey]>0:sumpos+=1
                else: print "Not found with correct PoS: "+candkey
                if self.countrel[candkey]>0:sumrel+=1
                else: print "Not found with dependency relation: "+candkey
            else: print "Not found: "+candkey


        print "Found "+str(sum)+" contiguous compounds"
        print "Found "+str(sumpos)+" with correct PoS tags"
        print "Found "+str(sumrel)+" with dependency relationship"
        print "Min count is "+str(min(self.counts.values()))
        print "Min count with dependency is "+str(min(self.countrel.values()))
        print "Max count is "+str(max(self.counts.values()))
        print "Max count with dependency is "+str(max(self.countrel.values()))
        print "Mean count is "+str(np.mean(self.counts.values()))
        print "Mean count with dependency is "+str(np.mean(self.countrel.values()))

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
        if self.clean:
            self.output_sentence(sentence,self.cleanstream)
        mxindex=len(sentence.values())-1+self.firstindex
        lemmamatch=False
        for i in sentence.keys():
            try:
                sid =int (i)
                if sid<mxindex:
                    arc=sentence[i]
                    canddep=getLex(arc[self.lex]).lower()
                    candhead=getLex(sentence[str(sid+1)][self.lex]).lower()

                    candkey=canddep+" "+candhead
                    if self.uselemma:
                        candkey2=getLex(arc[self.lemma]).lower()+" "+getLex(sentence[str(sid+1)][self.lemma]).lower()
                        candkey3=getLex(arc[self.lex].lower())+" "+getLex(sentence[str(sid+1)][self.lemma]).lower()
                        if candkey2 in self.compounds.keys():
                            lemmamatch=True
                            candkey=candkey2

                        elif candkey3 in self.compounds.keys():
                            lemmamatch=True
                            candkey=candkey3
                        else:
                            lemmamatch=False


                    #print candkey
                    if lemmamatch or (candkey in self.compounds.keys()):
                        #print "Found",candkey
                        self.counts[candkey]+=1
                        self.cont+=1
                        if self.counts[candkey]==1:
                            print "First found compound: line "+str(self.lines)
                            print i, arc, sentence[str(sid+1)],sentence[arc[self.headpos]]

                        if self.postagged:
                            deppos=getPos(arc[self.lex])
                            dephead=getPos(sentence[str(sid+1)][self.lex])
                        else:
                            deppos=arc[self.pos][0]
                            dephead=sentence[str(sid+1)][self.pos][0]

                        if deppos=="N" and dephead=="N":
                            self.countpos[candkey]+=1
                        if arc[self.headpos]==str(sid+1): #dependency relationship
                            sofar=self.rels.get(arc[self.relname],0)
                            self.rels[arc[self.relname]]=sofar+1

                            self.countrel[candkey]+=1
                            if self.convert: sentence=self.convert_compounds(sentence,sid)

                        else:
                            if sentence[str(sid+1)][self.headpos]==str(sid): #dependency relationship in other direction - ignore this so as to not create cyclic dependencies
                            #if False:
                                print "Not converting reverse dependency into compound"
                                self.cont_revdep+=1
                            else:
                                print candkey, " :contiguous but no dependency: ",i,arc,sentence[str(sid+1)],sentence[arc[self.headpos]]
                                self.cont_nodep+=1
                                if self.convert: sentence=self.convert_compounds(sentence,sid)
            except:
                print "Warning: error ignored"
                pass

        if self.convert: self.output_sentence(sentence, self.outstream)

    def convert_compounds(self,sentence,sid):
        if self.ptype =="conll7" or self.ptype=="nyt":
            deppos = self.getPosTag(sentence[str(sid)][self.pos])
        else:
            deppos = getPos(sentence[str(sid)][self.lex])

        if self.uselemma:
            token_to_use=self.lemma
        else:
            token_to_use=self.lex
        cmpd=sentence[str(sid)][token_to_use]+"|"+sentence[str(sid)][self.relname]+"-"+deppos+"|"+sentence[str(sid+1)][token_to_use]
        sentence[str(sid+1)][token_to_use]=cmpd
        if self.lemma >-1 and not self.uselemma:
            sentence[str(sid+1)][self.lemma]=sentence[str(sid)][self.lemma]+"|"+sentence[str(sid)][self.relname]+"|"+sentence[str(sid+1)][self.lemma]

        sentence[str(sid)]=self.erased
        #print "Converting compounds: ",cmpd,sid
        #print sentence.keys()
        for index in sentence.keys():
            if len(sentence[index])==self.arclength:

                #print index,sentence[index]
                if sentence[index][self.headpos]==str(sid):
            #      print "Found dependency on dependency"
                   sentence[index][self.headpos]=str(sid+1)
        #print "Complete"

        return sentence

    def output_sentence(self,sentence,out):
        slength=len(sentence.keys())
        index_adj=1-self.firstindex
        for i in range(self.firstindex,slength+self.firstindex):
            index=str(i)

            arc=sentence.get(index,None)
            if arc!=None:



                out.write(str(i+index_adj))
                for value in self.aptTransform(arc):
                    out.write("\t"+str(value))
                out.write("\n")

        out.write("\n")

    def aptTransform(self,arc):
    #ptype=conll7: arc = [form,lemma,POS,NER,head,rel]

    #aptInput requires: arc = [form/POS,head,rel]
        if self.uselemma:
            token_to_use = self.lemma
        else:
            token_to_use = self.lex

        index_adj=1-self.firstindex
        if self.ptype=="conll7" or self.ptype=="nyt":
            return [arc[token_to_use].lower()+"/"+self.getPosTag(arc[self.pos]),int(arc[self.headpos])+index_adj,arc[self.relname]]
        else:
            return [getLex(arc[token_to_use],tdelim='|')+"/"+getPos(arc[token_to_use]),int(arc[self.headpos])+index_adj,arc[self.relname]]

    def getPosTag(self,tag):
        newtag=tag[0]
        if newtag in CompoundFinder.postags: return newtag
        else: return "F"

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
                        if self.ptype=="nyt" or self.ptype=="conll7":
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