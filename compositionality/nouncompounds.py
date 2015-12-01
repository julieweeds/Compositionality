__author__ = 'juliewe'
#support composition of compounds and comparison with observed vectors

from composition import Composition
import sys

def union(list1,list2):

    for value in list2:
        if value not in list1:
            list1.append(value)

    return list1

class Compound:
    leftRels={"J":["amod","mod","amod-j"],"N":["nn","nn-n"],"V":["amod-v"]}
    rightRels={"N":["amod","nn","mod","amod-j","amod-v","nn-n"],"J":[],"V":[]}

    def __init__(self,text):
        self.text=text
        if not self.verify():
            print "Error: incorrect format for compound "+self.text
            exit(-1)

    def verify(self):
        if len(self.text.split('|'))==3:
            if len(self.text.split('|')[2].split('/'))==2:
                return True
        return False


    def getLeftLex(self):
        return self.text.split('|')[0]

    def getRel(self):
        return self.text.split('|')[1]

    def getRightLex(self):
        return self.text.split('|')[2].split('/')[0]

    def getPos(self):
        return self.text.split('|')[2].split('/')[1]

    def getWordsByPos(self,pos):
        #check not lower case for pos
        words=[]
        if self.getRel() in Compound.leftRels.get(pos,[]):
            words.append(self.getLeftLex()+"/"+pos)
        if self.getRel() in Compound.rightRels.get(pos,[]):
            words.append(self.getRightLex()+"/"+pos)
        return words

    def toString(self):
        return self.getRel()+":"+self.getLeftLex()+"\t"+self.getRightLex()+"\t"+self.getPos()

class DepCompounder:

    def __init__(self,config):

        try:
            self.datadir=config.get('compounder','datadir')
        except:
            self.datadir=""

        self.compoundfile=self.datadir+config.get('compounder','compound_file')
        self.leftindex={}
        self.relindex={}
        self.rightindex={}
        self.wordsByPos={"J":[],"N":[],"V":[]}

    def readcompounds(self):
        print "Reading "+self.compoundfile
        with open(self.compoundfile) as fp:
            for line in fp:
                line=line.rstrip()
                acompound=Compound(line)

                self.leftindex=self.addtoindex(acompound.getLeftLex(),self.leftindex,acompound)
                self.rightindex=self.addtoindex(acompound.getRightLex(),self.rightindex,acompound)
                self.relindex=self.addtoindex(acompound.getRel(),self.relindex,acompound)

                for pos in self.wordsByPos.keys():
                    self.wordsByPos[pos]=union(self.wordsByPos[pos],acompound.getWordsByPos(pos))


    def addtoindex(self,key,index,acompound):
        """

        :type index: dict
        """
        sofar=index.get(key,[])
        sofar.append(acompound)
        index[key]=sofar
        return index

    def run(self):
        self.readcompounds()
        print "Compounder stats... "
        print "Left index: "+str(len(self.leftindex.keys()))
        print "Right index: "+str(len(self.rightindex.keys()))
        print "Rel index: "+str(len(self.relindex.keys()))



class NounCompounder(Composition):
    left=Composition.depPoS
    right=Composition.headPoS
    basicRel=Composition.basicRel
    rels_to_include=[]

    def set_words(self):
        self.words=self.myCompounder.wordsByPos[self.pos]
        if self.words==[]:self.words=['-']


    def includeRel(self,rel):
        return NounCompounder.rels_to_include==[] or rel in NounCompounder.rels_to_include


    def compose(self,parampair=('','')):

        self.outfile=self.getComposedFilename()

        for pos in ["N","J","V"]:
            self.pos=pos
            self.set_words()
            self.feattotsbypos[pos]=self.load_coltotals(cds=self.smooth_ppmi)
            self.totsbypos[pos]=self.load_rowtotals()
            self.vecsbypos[pos]= self.load_vectors()
            self.pathtotsbypos[pos]=self.compute_nounpathtotals(self.vecsbypos[pos])
            self.typetotsbypos[pos]=self.compute_typetotals(self.feattotsbypos[pos])


        append=False
        for rel in self.myCompounder.relindex.keys():
            if self.includeRel(rel):
                self.runANcompositionByRel(rel,parampair=parampair)

            self.output(self.runANcomposition(parampair=parampair),self.outfile,append=append)
            append=True


    def runANcompositionByRel(self,rel,parampair=('','')):
        myvectors={}
        if parampair[0]=="offsetting":
            offsetting=float(parampair[1])
        else:
            offsetting=self.offsetting

        print "Composing type totals for "+rel
        self.ANtypetots=self.doCompound(self.typetotsbypos[NounCompounder.left[rel]],self.typetotsbypos[NounCompounder.right[rel]],rel,hp=self.headp,op=self.compop,offsetting=offsetting)  #C<*,t,*>
        print "Composing feature totals for "+rel
        self.ANfeattots=self.doCompound(self.feattotsbypos[NounCompounder.left[rel]],self.feattotsbypos[NounCompounder.right[rel]],rel,hp=self.headp,op=self.compop,offsetting=offsetting)  #C<*,t,f>

        self.ANvecs={}
        self.ANtots={}
        self.ANpathtots={}

        for compound in self.myCompounder.relindex[rel]:
            #should check not lower case for pos
            try:
                #print "Composing: "+compound.text
                self.CompoundCompose(compound.getLeftLex()+"/"+NounCompounder.left[rel],compound.getRightLex()+"/"+NounCompounder.right[rel],rel,hp=self.headp,compop=self.compop,offsetting=offsetting)

            except KeyError:
                pass
                #print "Warning: 1 or more vectors not present for "+compound.text

        myvectors.update(self.mostsalientvecs(self.ANvecs,self.ANpathtots,self.ANfeattots,self.ANtypetots,self.ANtots)) #compute ppmi vectors and store in myvectors
        return myvectors

    def getLeftIndex(self):
        return self.myCompounder.leftindex.keys()

    def getRightIndex(self):
        return self.myCompounder.rightindex.keys()

    def run(self,parampair=('','')):
        self.option=self.options[0]
        self.myCompounder=DepCompounder(self.config)
        self.myCompounder.run()
        self.compose(parampair=parampair)


if __name__=="__main__":
    myCompounder=NounCompounder(["config",sys.argv[1]])
    myCompounder.run()