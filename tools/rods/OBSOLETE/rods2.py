#!/usr/bin/python

import sys
import random as rand
import time
import math
import cProfile, pstats, StringIO


class Rod:
    def __init__(self, Np, Lrod, pos, director):
        self.rod = []
        self.pos = pos
        
        pp=[1,1, 0]
        pp[2]=(-director[0]  - director[1]) / director[2]
        ppn = math.sqrt(sum(abs(d*d) for d in pp))
        pp = [ d/ppn for d in pp ]

        for i in range(0,Np):
            for j in range(0,3):
                self.rod.append(0.0) 
                
        for i in range(0,Np):
            for j in range(0,3):
                self.rod[3*i+j] = pos[j]+(2.0*i/Np-1)*0.5*Lrod*director[j]
        for j in range(0,3):
            self.rod[3*(Np-1)+j] = pos[j]-0.5*Lrod*director[j] +pp[j]
        for j in range(0,3):
            self.rod[3*(Np-2)+j] = pos[j]-0.5*Lrod*director[j] -pp[j]

    def collides(self, other_rod):
        Np = len(self.rod)/3;
        for i in range(0,Np):
            for j in range(0,Np):
                dist2 = sum((a-b)**2 for a,b in zip(self.rod[3*i:3*i+3], other_rod.rod[3*j:3*j+3]))
                if(dist2<4.0):
                    return True
        return False

    def plist(self):
        buff=""
        for i in range(0,Np):
            buff+='%.3f %.3f %.3f %d\n' % (self.rod[3*i], self.rod[3*i+1], self.rod[3*i+2], (i/(Np-2)))
        return buff

    def bondlist(self, irod, k, r0):
        buff=""
        buff+='%d %d %d %.3f %.3f \n' % (Np*irod+Np-2, Np*irod, Np*irod+Np-1, k, r0)
        for i in range(1,Np-3):
            buff+='%d %d %d %.3f %.3f \n' % (Np*irod+i-1, Np*irod+i, Np*irod+i+1, k, r0)
        return buff
        


def random_rod(Np, Lrod, Lbox):
    dire = [rand.random()-0.5 for i in range(0,3) ]
    dirnorm = math.sqrt(sum(abs(d*d) for d in dire))
    dire = [ d/dirnorm for d in dire ]
    
#    Lbox-=Lrod*0.5;
    pos = [(rand.random()-0.5)*Lbox for i in range(0,3)]
    return Rod(Np, Lrod, pos, dire);
    
    
f = open('read.in','r');

line = f.readline().split();

Nrods = int(line[0]);
Np = int(line[1]);
Lrod = float(line[2]);

line = f.readline().split();
L = float(line[0])

line = f.readline().split();
k = float(line[0])
r0 = Lrod/Np;

f.close();


rand.seed(time.time());
rods = []
rods.append(random_rod(Np, Lrod, L))
trycount = 0


#pr = cProfile.Profile()
#pr.enable()
count = 0;
while(len(rods)<Nrods and trycount<1000 and count<10000):
    count+=1 
    accepted=True
    if(trycount>500):
        print len(rods),trycount
    a_rod = random_rod(Np, Lrod, L);
    for i in range(0,len(rods)):
        if(rods[i].collides(a_rod)):
            trycount+=1
            accepted=False
            break;
    if(accepted):
        rods.append(a_rod)
        trycount = 0

#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()

if(trycount>0):
    print "ERROR!: Could not achieve desired density, only "+str(len(rods))+" were written"


#Np = len(rods);
fout = open('rods.pos','w');
print >> fout, Np*Nrods
for i in range(0,Nrods):
   print >> fout, rods[i].plist(),
fout.close();

fbond = open('rods.3bonds','w');
print >> fbond, Nrods*(Np-3)
for i in range(0,Nrods):
    print >> fbond, rods[i].bondlist(i,k, r0),
fbond.close();    





