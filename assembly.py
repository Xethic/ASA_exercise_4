import sys
from string import maketrans

alphabet = ['A', 'C', 'G', 'T']
original = 'ACGT'
reversecomplement = 'TGCA'
revtable = maketrans(original, reversecomplement)


def reverseComplement(sequence):
	return sequence.translate(revtable)


def addreverse(reads):
	revreads = []
	for read in reads:
		revreads.append(reverseComplement(read)[::-1])
	reads.extend(revreads)


def computekmers(k, reads):
	kmer = dict()
	count = 0
	for r in reads:
		count+=1
		for i in range(len(r)-k+1):
			if r[i:i+k] in kmer:
				kmer[r[i:i+k]] += 1
			else:
				kmer[r[i:i+k]] = 1
	print 'Parsed ' + str(count) + ' reads'
	return kmer
	
	
def bestHammingNeighbor(v, cutoff, kmers):
	bestSeq = v
	bestCount = cutoff-1
	for i in range(0, len(v)):
		for a in alphabet:
			new = list(v)
			new[i] = a
			new = ''.join(new)
			if(new in kmers):
				if(kmers[new] > bestCount):
					bestSeq = new
					bestCount = kmers[new]
	return bestSeq
	
	
def errorcorrect(reads, cutoff, k):
	
	kmers = computekmers(k, reads)
	print 'Correcting errors....'
	for read in reads:
		for i in range(len(read)-k+1):
			v = read[i:i+k]
			if(not v in kmers):
				kmers[v] = 1
			if(kmers[v] < cutoff):
				new = bestHammingNeighbor(v, cutoff, kmers)
				if new != v:
					kmers[new] += 1
					read = read[0:i] + new + read[i+k: len(read)]
	print 'Finished!'
			
	
def getneighbor(node, kmers, k):
	temp = dict(kmers)
	
	converged = False
	olds1 = ""
	olds2 = ""
	while not converged:
		s1 = node[:k-1]
		s2 = node[-(k-1):]
		#print s2
		#print s1
		found1 = False
		found2 = False
		oldlength = len(node)
		#print s1, s2
		for c in ['A', 'C', 'G', 'T']:
			#print s2+c, c+s1
			if s2+c in temp:
				if not temp[s2+c] == 0:
					if not found1:					
						node += c
						found1 = True
						temp[s2+c] = 0
	
			if c+s1 in temp:
				if not temp[c+s1] == 0:
					if not found2:
						node = c + node
						found2 = True
						temp[c+s1] = 0
		
		if oldlength == len(node):
			converged = True
			
				
	return node


def countKmersOfCountOne(kmers, k):
	count = 0
	for kmer in kmers:
		if kmers[kmer] == 1:
			count +=1
	
	print 'Found ' + str(count) + ' ' +str(k) + '-mers that were present once' 
	
	
def assemble(kmers, k):
	finals = []
	for kmer in kmers:
		f = getneighbor(kmer, kmerlist, k)
		if f not in finals:
			finals.append(f)
		#exit()
	return finals


#######
#	MAIN
#
#	first argument is the input file*
#	second argument is the name of the output file
#######

reads = []
readset = set()

with open(str(sys.argv[1])) as f:
	for line in f:
		if line == "\n":
			break
			
		if line[0] == ">":
			continue
		else:
			read = line.replace('\n', '')
			reads.append(read)
			readset.add(read)

k = 15
errorcutoff = 2
assemb = []

addreverse(reads)

errorcorrect(reads, errorcutoff, k)

kmerlist = computekmers(k, reads)
#print kmerlist

print 'Generated ' + str(len(kmerlist)) + ' many ' + str(k) + '-mers'
countKmersOfCountOne(kmerlist, k)

if len(kmerlist) == 0:
	exit()

finals = assemble(kmerlist, k)
#print finals

with open(str(sys.argv[2]), "w") as out:
	for f in finals:
		out.write(str(f)+"\n")
		
#print len(finals)
#print max(finals)
#finals = list(finals)
#print str(max(finals))
