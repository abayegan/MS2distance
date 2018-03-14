#A. Bayegan and P. Clote

from aux import *

def printPath_pk(s,t,seq,rmlist,addlist,shlist):
	if len(seq)>0:
		print seq
		for i in range(1,len(seq)+1):
			sys.stdout.write("%d"%(i%10))
		sys.stdout.write("\n\n")
	bps = getBPs(s)
	slist=list(s)
	print s.strip(),'initial strucutre'
	for x,y in rmlist:
		slist[x]='.'
		slist[y]='.'
		ss = ''.join(slist).strip()
		print "%s remove (%d,%d)"%(ss,x,y)
	for x,y,z in shlist:
		slist[z]='.'
		slist[y]='.'
		if(x<y):
			slist[x]='('
			slist[y]=')'
		else:
			slist[x]=')'
			slist[y]='('
		ss =  ''.join(slist).strip()
		print "%s shift (%d,%d) -> (%d,%d)"%(ss,min(y,z),max(y,z),min(x,y),max(x,y))
	for x,y in addlist: 
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

def getComponentMoves(c,loc):
	path=c.path
	m = len(path)
	shList=[];rmList=[];adList=[] # shList is stored as (t-portion,pivot,s-portion)
	if c.type==1:
		rmList.append(sorted((path[0],path[1])))
		for i in range(1,m-2,2):
			if(abs(path[i]-path[i+2])<=loc):
				shList.append((path[i],path[i+1],path[i+2]))
	elif c.type==2:
		for i in range(m-1,1,-2): #[m-1,2]
			if(abs(path[i]-path[i-2])<=loc):
				shList.append((path[i],path[i-1],path[i-2]))
	elif c.type==3:
		for i in range(0,m-2,2):
			if(abs(path[i]-path[i+2])<=loc):
				shList.append((path[i],path[i+1],path[i+2]))
	elif c.type==4:
		adList.append(sorted((path[m-2],path[m-1])))
		for i in range(0,m-2,2):
			if(abs(path[i]-path[i+2])<=loc):
				shList.append((path[i],path[i+1],path[i+2]))

	#Note: these cycles have nothing to do with the cycles of the conflict graph
	elif c.type==5:
		rmList.append(sorted((path[m-2],path[m-1])))
		adList.append(sorted((path[m-3],path[m-2])))
		for i in range(0,m-4,2):
			if(abs(path[i]-path[i+2])<=loc):
				shList.append((path[i],path[i+1],path[i+2])) 
	return adList,rmList,shList

def pk_distance(s,t,seq,outprefix,loc):
	rmList=[];addList=[];shiftList=[]
	#plotArcsR4RNA(s,t,outprefix)

	#use 1-indexed
	s = ' ' +s
	t = ' '+ t
	sbp = getBpList(s)
	tbp = getBpList(t)
	n = len(s)
	A,B1,B2,C,D,V,addList1,rmList1= partitionPositions(sbp,tbp,loc)
	complist = findComponents(A,sbp,tbp)
	for c in complist:
		ad,rm,sh = getComponentMoves(c,loc)
		addList.extend(ad)
		rmList.extend(rm)
		shiftList.extend(sh)
	scovered=[];tcovered=[]
	for v in shiftList:
		tcovered.append(sorted((v[0],v[1])))
		scovered.append(sorted((v[1],v[2])))

	for x in getBPs(s):
		if x not in scovered:
			if x not in rmList:
				rmList.append(x)
	for x in getBPs(t):
		if x not in tcovered:
			if x not in addList:
				addList.append(x)
	printPath_pk(s,t,seq,rmList,addList,shiftList)
	return
