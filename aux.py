#A. Bayegan and P. Clote

import sys,itertools,os,time,ast,bisect,shlex
import networkx as nx
from collections import namedtuple 
from subprocess import PIPE,Popen
from edgedfs import helper_funcs, edge_dfs
import cycles
import ast,bisect

outprefix=""
NUC = ['A','C','G','U','T','a','c','g','u','t']
DEBUG =0
RCHIE_PATH = "/home/clotelab/Amir/ms2path/src/rchie.R"

def getBpList(s):
	n = len(s)
	L = [-1]*n; q=[]
	for i in range(1,n):
		if s[i]=="(":
			q.append(i)
		elif s[i]==")":
			p = q.pop()
			L[i] =p
			L[p]=i
	return L

def plotArcsR4RNA(s,t,outname):
	f1 = open('f1.fa','w')
	f1.write(">1\n%s" % t)
	f1.close()
	f2 = open('f2.fa','w')
	f2.write(">2\n%s" % s)
	f2.close()
	cmd = "%s --format1=vienna --format2=vienna --colour1=\"blue\" --colour2=\"red\" --quiet f1.fa f2.fa --pdf --output=%s.pdf" %(RCHIE_PATH,outname)
	os.popen(cmd)
	os.remove('f1.fa')
	os.remove('f2.fa')


#partition positions in s->t path to 4 sets and simultaneously construct node set V
#A: base-paired in both s and t, but the base pairs in s and t are not identical
#B1: base-paired in s not in t (base pairs are added to rmList)
#B2: base-paired in t not in s (base pairs are added to addList)
#C: base-paired in neither s nor t
#D: base-paired to the same partner in both s and t
#V: nodes of the conflict graph. (x,y,z) s.t {x,y} in t and {y,z} in s
def partitionPositions(sbp,tbp,loc):
	if(len(sbp) != len(tbp)):
		printerr("Error: source and target structures must have the same legnth")
	A=[];B1=[];B2=[];C=[];D=[];V=[];v_loc=[]
	addList=[];rmList=[]
	n = len(sbp)
	for i in range(1,n):
		if sbp[i]!=-1 and tbp[i]!=-1:
			if sbp[i]!=tbp[i]: 
				A.append(i)
				#print i,sbp[i],tbp[i]
				if(abs(sbp[i]-tbp[i])<=loc):
					V.append((tbp[i],i,sbp[i]))
				else:
					rmList.append((min(i,sbp[i]),max(i,sbp[i])))#WARNING:might remove a bp twice here
					addList.append((min(i,tbp[i]),max(i,tbp[i])))#WARNING:might add a bp twice here
			elif sbp[i]==tbp[i]:
				D.append(i)
		elif sbp[i]!=-1 and tbp[i]==-1:
			B1.append(i)
		elif sbp[i]==-1 and tbp[i]!=-1:
			B2.append(i)
		elif sbp[i]==-1 and tbp[i]==-1:
			C.append(i)
		if sbp[i]>i and tbp[i]==-1 and  tbp[sbp[i]]==-1:
			rmList.append((i,sbp[i]))
			#print "rm*",i,sbp[i]
		elif sbp[i]==-1 and tbp[i]>i and  sbp[tbp[i]]==-1:
			#print "add*",i,tbp[i]
			addList.append((i,tbp[i]))
	return A,B1,B2,C,D,V,set(addList),set(rmList)


def getBPs(s):
	bplist=[];stack=[]
	for i in range(len(s)):
		if s[i]=='(':	
			stack.append(i)	
		elif s[i]==')':
			j =stack.pop()
			bplist.append([j,i])
	return bplist

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


def checkSeq(seq):
	for nuc in seq:
		if nuc not in NUC:
			printerr("Input sequence should include only {A,C,G,T,U}")

def isCrossing(bp1,bp2):
	a,b = sorted(bp1)
	c,d = sorted(bp2)
	assert(a>-1 and b>-1 and c>-1 and d>-1)
	return (a<c<b<d or c<a<d<b)

#given a structure and a base pair x checks if adding x creates psuedoknots or triples
def isCompatibleBP(s,bp):
	s=' '+s
	bps = getBPs(s)
	for bp1 in bps:
		if isCrossing(bp1,bp):
			print 'Crossing between:',bp1,bp
			return False
	return True

def drawGraphWithDot(V,E,oname):
	txt =''
	V1=[]
	for x,y in E:
		V1.extend([x,y])
	isolated = set(V) - set(V1)
	for e in E:
		txt += '\"%s\"->\"%s\" \n' %(e[0],e[1])
	for v in isolated:
		txt+='\"(%d,%d,%d)\"'%(v[0],v[1],v[2])
	dot = """digraph G {
layout = "circo"
%s
node [shape=plaintext, fontsize=16];
}""" %txt
	outf = open(oname+'.dot','w')
	outf.write(dot)
	outf.close()
	#cmd = "dot -Teps %s.dot  -o %s.eps" %(oname,oname)
	#os.popen(cmd)


def printerr(msg):
	print msg
	sys.exit(1)
	
def parseInput(infile):
	infile = open(infile,'r')
	lines = infile.readlines()
	s = lines[0].strip()
	t = lines[1].strip()
	seq=""
	if len(lines)>2 and (lines[2][0] in NUC):
		seq=lines[2].strip()
	return s,t,seq

def getBPs(s):
	bplist=[];stack=[]
	for i in range(len(s)):
		if s[i]=='(':	
			stack.append(i)	
		elif s[i]==')':
			j =stack.pop()
			bplist.append([j,i])
	return bplist

def checkPath(s,t,rmlist,addlist,shlist):
	bps = getBPs(s)
	slist=list(s)
	for x,y in rmlist:
		slist[x]='.'
		slist[y]='.'
		sd = ''.join(slist).strip()
	for x,y,z in shlist:
		slist[z]='.'
		slist[y]='.'
		if(not isCompatibleBP(''.join(slist).strip(),(x,y))):
			print "ERROR IN ",s,t
			sys.exit(1)
		if(x<y):
			slist[x]='('
			slist[y]=')'
		else:
			slist[x]=')'
			slist[y]='('
	for x,y in addlist:
		if(not isCompatibleBP(''.join(slist).strip(),(x,y))):
			print "ERROR IN",s,t
			
		if(x<y):
			slist[x]='('
			slist[y]=')'
		else:
			slist[x]=')'
			slist[y]='('

def printPath(s,t,seq,rmlist,addlist,shlist):
	if len(seq)>0:
		print seq
		for i in range(1,len(seq)+1):
			sys.stdout.write("%d"%(i%10))
		sys.stdout.write("\n\n")
	bps = getBPs(s)
	slist=list(s)
	print s.strip(),"initial structure"
	for x,y in rmlist:
		slist[x]='.'
		slist[y]='.'
		ss = ''.join(slist).strip()
		print "%s remove (%d,%d)"%(ss,x,y)
	for x,y,z in shlist:
		slist[z]='.'
		slist[y]='.'
		assert(isCompatibleBP(''.join(slist).strip(),(x,y)))
		if(x<y):
			slist[x]='('
			slist[y]=')'
		else:
			slist[x]=')'
			slist[y]='('
		ss =  ''.join(slist).strip()
		print "%s shift (%d,%d) -> (%d,%d)"%(ss,min(y,z),max(y,z),min(x,y),max(x,y))
	for x,y in addlist: 
		assert(isCompatibleBP(''.join(slist).strip(),(x,y)))
		if(x<y):
			slist[x]='('
			slist[y]=')'
		else:
			slist[x]=')'
			slist[y]='('
		ss =  ''.join(slist).strip()
		print "%s add (%d,%d)"%(ss,x,y)
	assert(t[1:]==ss)
	numadd = len(addlist)
	numrm = len(rmlist)
	numshift = len(shlist)
	print "Number of base pair removals: %d\nNumber of base pair additions: %d\nNumber of base pair shifts: %d\nMS2 Distance: %d\n "% \
	(numrm,numadd,numshift,numadd+numrm+numshift)


#Given a component, return a list of base pair belonging to t
#and a list of bases (NOT pairs) belonging to s
def partitionPaths(path,ptype):
	if ptype==1:
		s_portion = sorted(path)
		t_portion = sorted(path[1:-1])
	elif ptype==2:
		s_portion = sorted(path[:-1])
		t_portion = sorted(path[1:])
	elif ptype==3:
		s_portion = sorted(path[1:])
		t_portion = sorted(path[:-1])
	elif ptype==4:
		s_portion = sorted(path[1:-1])
		t_portion = sorted(path)
	elif ptype==5:#The first edge in cycles(type 5) belongs t
		s_portion = sorted(path)
		t_portion = sorted(path)
	return s_portion,t_portion

#A is the list of pivot points
#output is a list of namedtuple Components ordered by the pivot points 
def findComponents(A,sbp,tbp):
	Component = namedtuple("Component","path type t_bases s_bases")
	complist=[]
	Q = set(A)
	while (len(Q)!=0):
		x0=min(Q);pathlist=[x0];nf=0;nb=0
		Q.discard(x0)
		x=tbp[x0]
		while (x!=-1 and x!=x0):
			nf+=1
			pathlist.append(x)
			Q.discard(x)
			if nf%2==1:
				x=sbp[x]
			else:
				x=tbp[x]
		nf+=1
		# assert(len(pathlist)==nf)
		if x==x0:#cycle instead of a path
			pathlist.append(x)
			nf+=1
			s_b,t_b = partitionPaths(pathlist,5)
			complist.append(Component(path=pathlist,type=5,t_bases=t_b,s_bases=s_b))
			continue
		x=sbp[x0]
		while(x!=-1):
			nb+=1
			pathlist.insert(0, x)
			Q.discard(x)
			if nb%2==1:
				x = tbp[x]
			else:
				x = sbp[x]
		# assert(len(pathlist)==nb+nf)
		ptype = findPathType(pathlist,sbp,tbp)
		s_b,t_b = partitionPaths(pathlist,ptype)
		complist.append(Component(path=pathlist,type=ptype,t_bases=t_b,s_bases=s_b))
	return complist 

#Determine 'path' is which of the 4 types:
#type 1: starts in s ends in s
#type 2: starts in s ends in t
#type 3: starts in t ends in s
#type 4: starts in t ends in t
def findPathType(path,sbp,tbp):
	a,b,x,y = path[0],path[1],path[-2],path[-1]
	if(sbp[a]==b and sbp[x]==y):
		ptype=1
	elif (sbp[a]==b and tbp[x]==y):
		ptype=2
	elif (tbp[a]==b and sbp[x]==y):
			ptype=3
	elif (tbp[a]==b and tbp[x]==y):
		ptype=4
	else:
		print "Error in findPathType(...)!"
		sys.exit(1)
	return ptype

def buildGraph(V):
	E = getEdges(V)
	G = nx.DiGraph()
	G.add_nodes_from(V)
	G.add_edges_from(E)
	return G
	
def getEdges(V):
	E=[]
	for v1,v2 in itertools.combinations(V,2):
		intersect = set(v1).intersection(v2)
		t1,p1,s1=v1
		t2,p2,s2=v2
		if len(intersect)<=1:
			if s1==t2 or isCrossing(sorted((s1,p1)),sorted((p2,t2))):
				E.append((v1,v2))
			if s2==t1 or isCrossing(sorted((s2,p2)),sorted((p1,t1))):
				E.append((v2,v1))
	return E

def isValidStructure(s):
	bplist=[];stack=[]
	ok = True
	if s.count('(')!=s.count(')'):
		print 'unbalanced!'
		return False					
	for i in range(len(s)):
		if s[i]=='(':	
			stack.append(i)	
		elif s[i]==')':
			j =stack.pop()
			bplist.append((j,i))
	n = len(bplist)
	for i in range(n):
		for j in range(i+1,n):
			if isCrossing(bplist[i],bplist[j]):
				print 'crossing'
				ok = False
	return ok
