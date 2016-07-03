import sys


def computekmers(k, reads):
	kmer = dict()
		
	for r in reads:
		for i in range(len(r)-k+1):
			if r[i:i+k] in kmer:
				kmer[r[i:i+k]] += 1
			else:
				kmer[r[i:i+k]] = 1
			
	return kmer


#def addreverse(kmerlist):
	
def errorcorrect(kmerlist):
	
	error = []
	
	for k in kmerlist:
		if kmerlist[k] == 1:
			if k not in error:
				error.append(k)
	
	for e in error:
		del kmerlist[e]
		
	return kmerlist
			
	
		
	
def getneighbor(node, kmers, k, currlen):
	ns = node
	s1 = node[:k-1]
	s2 = node[-(k-1):]
	found1 = False
	found2 = False
	#print s1, s2
	for c in ["A", "C", "G", "T"]:
		#print s2+c, c+s1
		if s2+c in kmers:
			if not found1:
				ns += c
				found2 = True
		if c+s1 in kmers:
			if not found2:
				ns = c + ns
				found2 = True
	
	#print ns

	if len(node) == len(ns):
		return ns
	else:
		ns = getneighbor(ns, kmers, k, len(ns))
		
	return ns
	

	
#main

#first argument is the input file
#second argument is the name of the output file
reads = []

with open(str(sys.argv[1])) as f:
	for line in f:
		if line == "\n":
			break
			
		if line[0] == ">":
			continue
		else:
			reads.append(line)

k = 6

assemb = []

kmerlist = computekmers(k, reads)
print kmerlist

#kmerlist = addreverse(kmerlist)

kmerlist = errorcorrect(kmerlist)
print kmerlist

if len(kmerlist) == 0:
	exit()

start = kmerlist.items()[0][0]
print start

finals = getneighbor(start, kmerlist, k, len(start))
print finals
#print assembl(start, kmers)

