#! /usr/bin/env python
#A. Bayegan

import sys,os,shlex,time
from subprocess import PIPE,Popen
from glob import glob

NUC = ['A','C','G','U','T','a','c','g','u','t']

def pipeCommand(cmd1,cmd2="",cmd3=""):
	c1 = Popen(shlex.split(cmd1),stdout=PIPE,stderr=PIPE)
	if cmd2=="":
		out,err = c1.communicate()
	else:
		c2 = Popen(shlex.split(cmd2),stdin=c1.stdout,stdout=PIPE,stderr=PIPE)
		if cmd3=="":
			out,err = c2.communicate()
		else:
			c3==Popen(shlex.split(cmd3),stdin=c2.stdout,stdout=PIPE,stderr=PIPE)
			out,err=c3.communicate()
	if(err):
		print("Error in %s|%s|%s\n" % (cmd1,cmd2,cmd3))
		print err
		sys.exit(1)
	return out


def parseInput(infile):
	infile = open(infile,'r')
	lines = infile.readlines()
	s = lines[0].strip()
	t = lines[1].strip()
	seq=""
	if len(lines)>2 and (lines[2][0] in NUC):
		seq=lines[2].strip()
	return s,t,seq

def getBpList(s):
          BpList = []; q = []; BaseList = []
          for i in range(len(s)):
            if s[i]=="(":
              q.append(i)
            elif s[i]==")":
              p = q.pop()
              BpList.append((p,i))
              BaseList.append(p); BaseList.append(i);
          #print BpList[::-1], pBaseList[::-1]
          return BpList[::-1]

def getBPdistance(s,t):
	bps = getBpList(s)
	bpt = getBpList(t)
	return len(set(bps)-set(bpt)) + len(set(bpt)-set(bps))	

def getBpArray(s):
	BpList = {}; q = [];
	for i in range(len(s)):
		BpList[i]=i
	for i in range(len(s)):
		if s[i]=="(":
			q.append(i)
		elif s[i]==")":
			p = q.pop()
			BpList[p]=i
			BpList[i]=p
	return BpList

def getHammingDistance(s,t):
	sbp = getBpArray(s)
	tbp = getBpArray(t)
	hdist=0
	for i in range(len(s)):
		if sbp[i]!=tbp[i]:
			hdist+=1
	return hdist
	
def main(s,t,seq,fname):
	bpdist = getBPdistance(s,t)
	hdist = getHammingDistance(s,t)
	bps = getBpList(s)
	bpt = getBpList(t)
	print '%s\t%d\t%d\t%d\t%d\t%d'%(fname,len(s),len(bps),len(bpt),bpdist,hdist)

if __name__ == '__main__':
	if len(sys.argv) < 2 :
		printUsage()
	args = sys.argv
	oname=""
	i=1
	while i<len(args):
		arg = args[i]
		if arg=='-i' and args[i+1][0]!='-':
			fname=args[i+1]
			i=i+1
		else:
			printerr('Error in the input arguments! Try -h for help.')
		i=i+1
	s,t,seq = parseInput(fname)
	main(s,t,seq,os.path.basename(fname))

