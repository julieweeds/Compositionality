
from operator import itemgetter
from composition import getorder, getpathtype, getpathvalue

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

    def suffix(self,path):

        feats={}
        for feat in self.features.keys():
            thispath = getpathtype(feat)
            if thispath.endswith(path):
                feats[feat]=self.features[feat]

        values = sorted(feats.items(),key=itemgetter(1),reverse=True)
        print values



class Viewer:

    def __init__(self):
        self.configure()
        self.vectors={}

    def configure(self):
        self.parentdir="/home/j/ju/juliewe/Documents/workspace/SentenceCompletionChallenge/data/apt/"
        self.filename="wiki_rbt/wiki_rbt_reddy_lemma_2.tsv"
        self.reducestring="reduce_0_2.filtered.norm"
        self.filesByPos={"N":"nouns","J":"adjs","V":"verbs","R":"advs","F":"other"}

    def loadvectors(self,pos,alist):

        infile=self.parentdir+self.filename+"."+self.filesByPos[pos]+"."+self.reducestring

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

    def suffix(self,entry,path):
        if self.load(entry):
            self.vectors[entry].suffix(path)
        else:
            return "Error: not found vector"

if __name__=="__main__":
    cont=True
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
                myviewer.suffix(parts[1],parts[2])

