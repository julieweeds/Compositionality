__author__ = 'juliewe'
#take in the list of noun compounds and generate other formats for composition etc

import sys, yaml

class Compound:

    root=("root",-1,0)
    pos=("?",-1)

    def __init__(self,name):
        self.name=name
        self.judgements=[]
        self.relmap={}
        self.deps=[]
        self.heads=[]
        self.posmap={}

    def addJudgement(self,NC):
        self.judgements.append(NC)

    def setJudgements(self,NCs):
        self.judgements=[int(item) for item in NCs]

    def setRel(self,rel):

        if len(rel)>0:
            dep=rel[1]
            self.relmap[dep]=rel
            for rel in self.relmap.values():
                self.deps.append(rel[1])
                if rel[2] not in self.deps:self.heads.append(rel[2])


    def setPos(self,poses):
        self.posmap=poses

    def getFirst(self):
        return self.name.split(" ")[0]

    def getSecond(self):
        return self.name.split(" ")[1]

    def getHead(self):

        if len(self.heads)==1:
            return self.name.split(" ")[self.heads[0]-1]
        else:
            print "Error: multiple or no heads"
            return self.name.split(" ")[self.heads[0]-1]

    def getDeps(self):
        return [self.name.split(" ")[dep-1] for dep in self.deps]

    def match(self,head,dep):

        return head == self.getHead() and dep in self.getDeps()

    def getScore(self):
        return sum(self.judgements)

    def getConll(self):
        tokens=self.name.split(" ")
        lines=[]
        for(id,token) in enumerate(tokens):
            sid=id+1
            rel=self.relmap.get(sid,Compound.root)
            pos=self.posmap.get(sid,Compound.pos)
            line=str(sid)+"\t"+token+"/"+pos[0]+"\t"+str(rel[2])+"\t"+rel[0]
            lines.append(line)

        return lines

    def display(self):
        for line in self.getConll():
            print line
        print self.getScore()


class Compounder:

    def __init__(self,filename):
        self.configfile=filename
        with open(self.configfile) as fp:
            self.configured=yaml.safe_load(fp)
        print self.configured
        self.compounds={} #key is name, value is list of compounds
        self.firstindex={} #key is first word, value is name
        self.secondindex={} #key is second word, value is name
        try:
            self.compound_file=self.configured["compound_file"]

        except:
            print "Error: problem with configuration"

    def readcompounds(self):

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





    def run(self):
        self.readcompounds()
        for compound in self.compounds.values():
                compound.display()
        return

if __name__=="__main__":
    #first argument should be name of config file
    myCompounder=Compounder(sys.argv[1])
    myCompounder.run()
