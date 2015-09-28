__author__ = 'juliewe'
#take in the list of noun compounds and generate other formats for composition etc

import sys, yaml

class Compound:

    root=("root",0,-1)
    pos=("N",-1)

    def __init__(self,name):
        self.name=name
        self.judgements=[]

    def addJudgement(self,NC):
        self.judgements.append(NC)

    def setJudgements(self,NCs):
        self.judgements=[int(item) for item in NCs]

    def setRel(self,rels):
        self.relmap={}
        for rel in rels:
            dep=rel[2]
            self.relmap[dep]=rel

    def setPos(self,poses):
        self.posmap={}
        for pos in poses:
            self.posmap[pos[1]]=pos


    def getScore(self):
        return sum(self.judgements)

    def getConll(self):
        tokens=self.name.split(" ")
        lines=[]
        for(id,token) in enumerate(tokens):
            sid=id+1
            rel=self.relmap.get(sid,Compound.root)
            pos=self.posmap.get(sid,Compound.pos)
            line=str(sid)+"\t"+token+"/"+pos[0]+"\t"+str(rel[1])+"\t"+rel[0]
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
        self.compounds=[]

        try:
            self.rel=self.configured.get("rels",[("nn",2,1)])
            self.pos=self.configured.get("pos",[])
            self.compound_file=self.configured["compound_file"]
        except:
            print "Error: problem with configuration"

    def readcompounds(self):

        try:
            with open(self.compound_file) as fp:
                for line in fp:
                    line=line.rstrip()
                    fields=line.split(",")
                    aCompound=Compound(fields[0])
                    aCompound.setRel(self.rel)
                    aCompound.setPos(self.pos)
                    aCompound.setJudgements(fields[1:])
                    self.compounds.append(aCompound)

        except:
            print "Error: problem reading compound file"

    def run(self):
        self.readcompounds()
        for compound in self.compounds:
            compound.display()
        return

if __name__=="__main__":
    #first argument should be name of config file
    myCompounder=Compounder(sys.argv[1])
    myCompounder.run()
