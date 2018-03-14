#! /usr/bin/env python
#A. Bayegan and P. Clote

import sys,itertools,os
import networkx as nx
from pk_distance import *
from edgedfs import helper_funcs, edge_dfs
from gurobipy import *
from near_optimal import near_optimal
from ms1_distance import getBPdistance,getHammingDistance
from aux import *


outprefix=""
NUC = ['A','C','G','U','T','a','c','g','u','t']
DEBUG =0

RCHIE_PATH = "./rchie.R"

# c is the output of simple_cycles() which does not include the last edge
#convert list of vertices in a c to list of edges
def getEdgesFromCycle(c):
	e = []
	for i in range(len(c)):
		if (i==len(c)-1):
			e.append((c[i],c[0]))
		else:
			e.append((c[i],c[i+1]))
	return e


#this is used to get proper names in the lp problem
def vertexList2strList(e):
	return map(lambda x: str(x).replace(' ' ,'') ,e)

#convert back the lp variable names to tuples
def strList2vertexList(v):
	l = v.strip('vertex_(').strip(')').split(',')
#	print (int(l[0]),int(l[1]),int(l[2]))
	return ((int(l[0]),int(l[1]),int(l[2])))

#give graph G, returns the minimum list of vertecies needs to be removed
#from G to make it acyclic.
def solve_MFVS_ILP(G,cycList,s,t):
	sbplist = getBPs(s)
	tbplist= getBPs(t)
	#define the LP problem
	model = Model("Minimum feedback vertex set problem model for the input file "+outprefix)
	#define variables
	V = G.nodes()
	n = len(V)
	vname2idx ={}
	for i in range(n):
		vname2idx[V[i]] = i
	
	nodevars = model.addVars(n, vtype=GRB.BINARY)
	for i in range(n):
		nodevars[i].setAttr(GRB.Attr.VarName,"vertex_(%d,%d,%d)"%(V[i][0],V[i][1],V[i][2]))
	#define objective
	model.update()

	objexp = LinExpr()
	for i in range(n):
		objexp += nodevars[i]
	model.setObjective(objexp,GRB.MAXIMIZE)
		
	#add constraints
	for i in range(len(cycList)):
		objexp = LinExpr()
		cyc=cycList[i]
		m = len(cyc)
		for v in cyc:
			objexp += nodevars[vname2idx[v]]
		model.addConstr(objexp <= m-1, "cycle_%d"%i)
	cnt=1
	#constraints on s base pairs for which 2 shift moves are possible
	for v1,v2 in itertools.combinations(G.nodes(),2):
		objexp = LinExpr()
		if len(set(v1).intersection(v2))==2:
			objexp = nodevars[vname2idx[v1]] + nodevars[vname2idx[v2]]
			model.addConstr(objexp<=1,'triplet_%d'%cnt)
			
	model.params.Outputflag=0
	model.optimize()
	if DEBUG:
		# Each of the variables is printed with it's resolved optimum value
		for v in model.getVars():
			 print('%s %g' % (v.varName, v.x))

	
	#----------------processing ilp output
	L1=[] #list of edges needs to be removed from G to make it acyclic
	for v in model.getVars():
		if v.x ==0:
			L1.append(strList2vertexList(v.varName))
	remove=[];add=[] #s base pairs need to be added or removed from the structures
	scovered=[];tcovered=[]
	for v in model.getVars():
		x = strList2vertexList(v.varName)
		if v.x==1:
			tcovered.append(sorted((x[0],x[1])))
			scovered.append(sorted((x[1],x[2])))
	uniq_sbplist =[x for x in sbplist if x not in tbplist ]
	uniq_tbplist =[x for x in tbplist if x not in sbplist ]
	for x in uniq_sbplist:
		if x not in scovered:
			if x not in remove:
				remove.append(x)
	for x in uniq_tbplist:
		if x not in tcovered:
			if x not in add:
				add.append(x)
	return L1,add,remove


def runms2(s,t,seq,outprefix,loc):
	rmList=[];addList=[];shiftList=[]
	#plotArcsR4RNA(s,t,outprefix)

	#use 1-indexed
	s = ' ' +s
	t = ' '+ t
	sbp = getBpList(s)
	tbp = getBpList(t)
	n = len(s)
	A,B1,B2,C,D,V,addList,rmList= partitionPositions(sbp,tbp,loc)
	if len(V)==0:
		printPath(s,t,seq,rmList,addList,shiftList)
		return

	G = buildGraph(V)
	num_edges = nx.number_of_edges(G)
	num_nodes = nx.number_of_nodes(G)
	v1=[]
	for e1,e2 in G.edges():
		v1.extend([e1,e2])
	num_isolated = num_nodes - len(set(v1))
	drawGraphWithDot(V,G.edges(),outprefix)
	
	#start_time = time.time()
	cycList = list(cycles.simple_cycles(G))#Call to a modified version of nx.simple_cycles(G). Modification allows a threshold on the number of cycles
	#cyctime = time.time() - start_time
	if len(cycList)!=0 and cycList[-1]==-1:#if threshold is reached in simple_cycles then cycList[-1]==-1
		printerr("Number of cycles in the conflict graph is greater than %d! Try near-optimal algorithm!"%(CYCTHRESH))
	else:
		numcyc = len(cycList)
		if (numcyc>IPTHRESH):
			printerr("Number of cycles in the conflict graph is %d! Try near-optimal algorithm!"%(numcyc))	
		else:
			l1,addilp,rmilp = solve_MFVS_ILP(G,cycList,s,t)
			G.remove_nodes_from(l1)
			rmList = rmilp
			addList = addilp
			sortedNodes= nx.topological_sort(G)
			for node in sortedNodes:
				sp,tp= sorted((node[2],node[1])),sorted((node[0],node[1]))
				shiftList.append(node)

			print 'Total number of nodes:',num_nodes
			print 'Number of isolated nodes:',num_isolated
			print 'Number of edges:',num_edges
			print 'Number of cycles:',numcyc
			print
	printPath(s,t,seq,rmList,addList,shiftList)
	return

def printUsage():
	print """Usage:%s -i inputFile [Options]
	
	-i <string>: Input file containing the initial structure in the first and the target structure in the second line. RNA sequence in line 3 is optional.
				 NOTE: The user is responsible for checking the Watson-Crick compatibility of base paring for the given sequence(if any).  If the given initial and target structures are watson-crick then intermediate structures would satisfy Watson-Crick, as well. However, the program is not limitted to Watson-Crick.
	Options:
	
	-pk  			Allow formation of psuedoknots
	-nopt			Near-optimal algorithm
	-o <string>	 	Prefix for the output names. Defualt is the input prefix.
	-k <int> 		Threshold for locality of shift moves
	-h  			Print help
	
	output:	
			1) MS2 path is printed in the standard output.
			2) dot file of the conflict digraph""" % (sys.argv[0])

	sys.exit(1)

if __name__ == '__main__':
	if len(sys.argv) < 2 :
		printUsage()
	args = sys.argv
	oname=""
	PK=0
	NOPT=0
	CYCTHRESH=10000000
	IPTHRESH = 10000000
	i=1
	loc=sys.maxint
	while i<len(args):
		arg = args[i]
		if arg=='-i' and args[i+1][0]!='-':
			fname=args[i+1]
			i=i+1
		elif arg=='-o' and args[i+1][0]!='-':
			oname = args[i+1]
			i=i+1
		elif arg=='-pk':
			PK=1
		elif arg=='-nopt':
			NOPT=1
		elif arg=='-k':
			loc=int(args[i+1])
			i=i+1
		elif arg=='-h':
			printUsage()
		else:
			printerr('Error in the input arguments! Try -h for help.')
		i=i+1
	s,t,seq = parseInput(fname)
	if oname=="":
		o = fname.split('/')
		oname = ''.join(o[:-1]) + o[-1].split('.')[0]
	if PK==0 and NOPT==0:
		runms2(s,t,seq,oname,loc)
	elif PK==1:
		pk_distance(s,t,seq,oname,loc)
	elif NOPT==1:
		near_optimal(s,t,seq,oname,loc)
	else:
		printerr('Either -pk of -nopt must be set.')
	print "MS1 distance:", getBPdistance(s,t)
	print "Hamming distance:", getHammingDistance(s,t)	
	print
