
from operator import itemgetter
from composition import getorder, getpathtype, getpathvalue
import sys

class Vector:

    def __init__(self,name,features):
        self.name=name
        self.features={}
        self.addall(features)

    def addall(self,features):

        while(len(features)>0):
            score=features.pop()
            feat=features.pop()
            self.features[feat]=float(score)

    def lookup(self,feature):
        return self.features.get(feature,-1)

    def profile(self,minorder=0,maxorder=10):

        paths={}
        totalweight=0
        thisorderweight=0
        for feat in self.features.keys():
            path = getpathtype(feat)
            order = getorder(feat)
            weight= self.features[feat]
            sofar=paths.get(path,0)
            if order>=minorder and order<=maxorder:
                paths[path]=sofar+weight
                thisorderweight+=weight
            totalweight+=weight

        print "total weight of features",totalweight
        print "total weight of required order features",thisorderweight
        profile=sorted(paths.items(),key=itemgetter(1),reverse=True)

        print profile

    def view(self,minorder=0,maxorder=10,cutoff=100):
        features={}
        for feat in self.features.keys():
            order = getorder(feat)
            if order>=minorder and order<=maxorder:
                features[feat]=self.features[feat]

        aview=sorted(features.items(),key=itemgetter(1),reverse=True)
        line=""
        for i,pair in enumerate(aview[:cutoff-1]):
            if i%5==4:
                line += pair[0]+" "+str(pair[1])
                print line
                line=""
            else:
                line += pair[0]+" "+str(pair[1])+"\t"


    def showsuffix(self,path,minorder=1,maxorder=1):

        feats={}
        for feat in self.features.keys():
            thispath = getpathtype(feat)
            order = getorder(feat)
            if thispath.endswith(path) and order >=minorder and order <=maxorder:
                feats[feat]=self.features[feat]

        values = sorted(feats.items(),key=itemgetter(1),reverse=True)
        print values



class Viewer:

    def __init__(self,fileoption="filtered"):
        self.configure(fileoption)
        self.vectors={}

    def configure(self,fileoption="filtered"):
        self.parentdir="/home/j/ju/juliewe/Documents/workspace/SentenceCompletionChallenge/data/apt/"
        self.filename="wiki_rbt/wiki_rbt_lemma_reddy_2.tsv"
        self.reducestring="reduce_0_2."
        self.fileoption=fileoption
        self.suffix=".norm.smooth_ppmi"
        self.filesByPos={"N":"nouns","J":"adjs","V":"verbs","R":"advs","F":"other"}

    def loadvectors(self,pos,alist):

        infile=self.parentdir+self.filename+"."+self.filesByPos[pos]+"."+self.reducestring+self.fileoption+self.suffix

        found=0
        with open(infile) as instream:
            print "Searching", infile

            for line in instream:
                line=line.rstrip()
                fields=line.split('\t')
                if fields[0] in alist:
                    self.vectors[fields[0]]=Vector(fields[0],fields[1:])
                    found+=1

        if found == len(alist):
            return True

        else:
            print "Found "+str(found)+" out of "+str(len(alist))
            return False

    def load(self,entry):
        parts=entry.split('/')
        if entry not in self.vectors.keys():
            found = self.loadvectors(parts[1],[entry])
        else:
            found = True

        return found

    def lookup(self,entry,feature):

        if self.load(entry):
            return self.vectors[entry].lookup(feature)
        else:
            return "Error: not found vector"

    def profile(self,entry,minorder=0,maxorder=10):
        if self.load(entry):
            self.vectors[entry].profile(minorder=minorder,maxorder=maxorder)
        else:
            return "Error: not found vector"

    def view(self,entry,minorder=0,maxorder=10):
        if self.load(entry):
            self.vectors[entry].view(minorder=minorder,maxorder=maxorder)
        else:
            return "Error: not found vector"

    def showsuffix(self,entry,path,minorder=1,maxorder=1):
        if self.load(entry):
            self.vectors[entry].showsuffix(path,minorder=minorder,maxorder=maxorder)
        else:
            return "Error: not found vector"

if __name__=="__main__":
    cont=True
    if len(sys.argv)>1:
        myviewer=Viewer(sys.argv[1])
    else:
        myviewer=Viewer()
    while cont==True:
        query=raw_input("Enter query: ")
        if query=="quit":
            cont=False
        else:
            parts=query.split(" ")
            if parts[0]=="lookup":
                result=myviewer.lookup(parts[1],parts[2])
                print "Result: "+str(result)
            elif parts[0]=="profile":
                if len(parts)>=4:
                    min=int(parts[2])
                    max=int(parts[3])

                    myviewer.profile(parts[1],min,max)
                else:
                    myviewer.profile(parts[1])
            elif parts[0]=="suffix":
                if len(parts)>=4:
                    min = int(parts[3])
                    max=int(parts[4])
                    myviewer.showsuffix(parts[1],parts[2],min,max)
                else:
                    myviewer.showsuffix(parts[1],parts[2])

            elif parts[0]=="view":
                if len(parts)>=4:
                    min=int(parts[2])
                    max=int(parts[3])

                    if len(parts)>=5:
                        cutoff=int(parts[4])
                        myviewer.view(parts[1],min,max,cutoff)
                    else:
                        myviewer.view(parts[1],min,max)
                else:
                    myviewer.view(parts[1])

