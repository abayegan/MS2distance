#A. Bayegan and P. Clote

import sys,itertools,os
import networkx as nx
from collections import namedtuple,Counter 
from gurobipy import *
from aux import *

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
	start_time = time.time()
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
		
	#constraints on s base pairs for which 2 shift moves are possible
	cnt=1
	for v1,v2 in itertools.combinations(G.nodes(),2):
		objexp = LinExpr()
		if len(set(v1).intersection(v2))==2:
			objexp = nodevars[vname2idx[v1]] + nodevars[vname2idx[v2]]
			model.addConstr(objexp<=1,'triplet_%d'%cnt)
	
	#model.write(outprefix+'.lp')
	model.params.OutputFlag=0
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
	#print add,remove
	return L1,add,remove


#A is the list of pivot points
#output is a list of namedtuple Components ordered by the pivot points 
def find_equivalence_classes(A,sbp,tbp):
	Component = namedtuple("Component","path type sportion tportion")
	complist=[]
	eqt=[];eqs=[]
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
				eqt.append(sorted((x,tbp[x])))
				x=sbp[x]
			else:
				eqs.append(sorted((x,sbp[x])))
				x=tbp[x]
		nf+=1
		# assert(len(pathlist)==nf)
		if x==x0:#cycle instead of a path
			pathlist.append(x)
			nf+=1
			eqs.append(sorted((x,sbp[x])))
			complist.append(Component(path=pathlist,type=5,tportion=eqt,sportion=eqs))
			eqt=[];eqs=[]
			continue
		x=sbp[x0]
		while(x!=-1):
			nb+=1
			pathlist.insert(0,x)
			Q.discard(x)
			if nb%2==1:
				eqs.append(sorted((x,sbp[x])))
				x = tbp[x]
			else:
				eqt.append(sorted((x,tbp[x])))
				x = sbp[x]
		# assert(len(pathlist)==nb+nf)
		ptype = findPathType(pathlist,sbp,tbp)
		complist.append(Component(path=pathlist,type=ptype,tportion=eqt,sportion=eqs))
		eqt=[];eqs=[]
	return complist 

def getEdges(V):
	E=[]
	for v1,v2 in itertools.combinations(V,2):
		intersect = set(v1).intersection(v2)
		t1,p1,s1=v1
		t2,p2,s2=v2
		# print v1,v2,intersect
		if len(intersect)<=1:
			if s1==t2 or isCrossing(sorted((s1,p1)),sorted((p2,t2))):
				E.append((v1,v2))
			if s2==t1 or isCrossing(sorted((s2,p2)),sorted((p1,t1))):
				E.append((v2,v1))
	return E


def get_edges_bw_eq_classes(s2eq,t2eq):
	E=[];Ew=[]
	crossDict={}
	for sbp,eqs in s2eq.items():
		for tbp,eqt in t2eq.items():
			if isCrossing(tbp,sbp) and eqs!=eqt:	
				E.append((eqs,eqt))
				if (eqs,eqt) not in crossDict.keys():
					crossDict[(eqs,eqt)]=[(sbp,tbp)]
				else:
					crossDict[(eqs,eqt)].append((sbp,tbp))
	#compute the weights
	d = Counter(E)
	for e in d.keys():
		Ew.append((e[0],e[1],d[e]))
	return Ew,crossDict	

def solve_feedback_arc(G,cyclist):
	model = Model("Minimum feedback arc problem"+outprefix)
	E=G.edges()
	Edata = G.edges(data=True)
	edge2idx = {}
	n = len(E)
	for i in range(n):
		edge2idx[E[i]]=i
	edgeVars = model.addVars(n,vtype = GRB.BINARY)
	for i in range(n):
		edgeVars[i].setAttr(GRB.Attr.VarName,"(%d,%d)"%(E[i][0],E[i][1]))
	objexp = LinExpr()
	for i in range(n):
		#print Edata[i][2]['weight']
		objexp += Edata[i][2]['weight'] * edgeVars[i]
	model.setObjective(objexp,GRB.MAXIMIZE)
	for i in range(len(cyclist)):
		cyc = cyclist[i]
		elist = getEdgesFromCycle(cyc)
		m = len(elist)
		objexp = LinExpr()
		for j in range(m):
			#print elist[j],edge2idx[elist[j]],'*'
			objexp += edgeVars[edge2idx[elist[j]]]
		model.addConstr(objexp<=m-1,'cycle_%d'%i)
	#model.write('a.lp')
	model.params.OutputFlag=0
	model.optimize()	
	L1=[] #list of edges needs to be removed from G to make it acyclic
	for v in model.getVars():
		if v.x==0:
			L1.append(ast.literal_eval(v.varName))
	#print L1
	return L1

def drawGraphWithDotw(V,E,oname,weighted):
	txt =''
	ee = set()
	for e0,e1,e2 in E:
		ee.update([e0,e1])
	isolated = set(V) - ee
	for e in E:
		txt += '\"%s\"->\"%s\"' %(e[0],e[1])
		if weighted:
			txt += '[label=%d,weight=\"%d\"];\n'%(e[2],e[2])
		else:
			txt+='\n'
	for v in isolated:
		txt+='\"%d\"'%(v)
	dot = """digraph G {
layout = "circo"
%s
node [shape=plaintext, fontsize=16];
}""" %txt
	outf = open(oname+'.dot','w')
	outf.write(dot)
	outf.close()
	cmd = "dot -Teps %s.dot  -o %s.eps" %(oname,oname)
	os.popen(cmd)

def near_optimal(s,t,seq,outprefix,loc):
	rmList=[];addList=[];shiftList=[]
	s2eq={};t2eq={}
	rmlist2=[]
	rmIp=[]
	plotArcsR4RNA(s,t,outprefix)

	#use 1-indexed
	s = ' ' +s
	t = ' '+ t
	sbp = getBpList(s)
	tbp = getBpList(t)
	n = len(s)
	A,B1,B2,C,D,V,addList,rmList= partitionPositions(sbp,tbp,loc)
	
	#phase 1: construct the conflict graph between equivalence classes
	eqlist = find_equivalence_classes(A,sbp,tbp)

	#for each bp determine the equivalance class
	for i in range(len(eqlist)):
		eq = eqlist[i]
		for s0 in eq.sportion:
			s2eq[tuple(s0)]=i
		for t0 in eq.tportion:
			t2eq[tuple(t0)]=i

	E,crossDict = get_edges_bw_eq_classes(s2eq,t2eq)

	G_eq = nx.DiGraph()
	V = range(len(eqlist))
	G_eq.add_nodes_from(V)
	G_eq.add_weighted_edges_from(E)
	drawGraphWithDotw(V,E,'equivalance_'+outprefix,True)
	#print '#number of equivalance classes:',len(V)
	#print '#number of crossings between eq classes:',len(E)
	cyclist = list(nx.simple_cycles(G_eq))
	vnum1 ,enum1,cnum1 = len(V),len(E),len(cyclist)
	if cnum1>0:
		rmIp =solve_feedback_arc(G_eq,cyclist)
	
	print '#number of nodes in equivalance graph:',vnum1
	print '#number of edges in equivalance graph:',enum1
	print '#number of cycles in equivalance graph:',cnum1
	print
	#remove all the s base pairs for edges found by IP
	rm1=[]
	ll = list(s)
	for e in rmIp:
		for s0,t0 in crossDict[e]:
			if s0 not in rm1:
				rm1.append(s0)
				sbp[s0[1]]=-1
				sbp[s0[0]]=-1
				ll[s0[0]]='.'
				ll[s0[1]]='.'
	new_s = ''.join(ll)
	#new equivalance classes can be computed from the previous ones
	#THIS SECTION SHOUlD BE IMPROVED
	A,B1,B2,C,D,V,addListprime,rmListprime= partitionPositions(sbp,tbp,loc)
	rmList=rm1
	#Find the order of solving the components using the conflict graph
	G = buildGraph(V)
	drawGraphWithDot(V,G.edges(),'conflict_'+outprefix)
	cyclist2 = list(nx.simple_cycles(G))
	enum2,vnum2,cnum2 = nx.number_of_edges(G),nx.number_of_nodes(G),len(cyclist2)
	l1,addilp,rmilp = solve_MFVS_ILP(G,cyclist2,new_s,t)
	G.remove_nodes_from(l1)
	rmList.extend(rmilp)
	addList = addilp
	sortedNodes= nx.topological_sort(G)
	for node in sortedNodes:
		sp,tp= sorted((node[2],node[1])),sorted((node[0],node[1]))
		shiftList.append(node)
	
	print '#number of nodes in conflict graph:',vnum2
	print '#number of edges in conflict graph:',enum2
	print '#number of cycles in conflict graph:',cnum2
	print
	#print 'Total number of nodes:',vnum2
	#print 'Number of edges:',enum2
	#print 'Number of cycles:',cnum2
	#print
	printPath(s,t,seq,rmList,addList,shiftList)
