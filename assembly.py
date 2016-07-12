import sys
from string import maketrans

alphabet = ['A', 'C', 'G', 'T']
original = 'ACGT'
complement = 'TGCA'
revtable = maketrans(original, complement)


def getComplement(sequence):
	return sequence.translate(revtable)


def addreverse(reads):
	revreads = []
	for read in reads:
		revreads.append(getComplement(read)[::-1])
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
			
	
def getneighbor(node, kmers, k, dispatched, maxrepeats):
	converged = False
	dispatched[node] = 1
	breakright = False
	breakleft = False
	while not converged:
		s1 = node[:k-1]
		s2 = node[-(k-1):]
		oldlength = len(node)
		rightextensions = []
		leftextensions = []
		for c in ['A', 'C', 'G', 'T']:
			if (not breakright) and s2+c in kmers:
				if s2+c not in dispatched :
					rightextensions.append(c)
			if (not breakleft) and c+s1 in kmers:
				if c+s1 not in dispatched :
					leftextensions.append(c)
		
		if(len(rightextensions) == 1):
			c = rightextensions[0]
			dispatched[s2+c] = 1
			node += c
		elif len(rightextensions) > 1:
			breakright = True
		
		if(len(leftextensions) == 1):
			c = leftextensions[0]
			dispatched[c+s1] = 1
			node = c + node
		elif len(leftextensions) > 1:
			breakleft = True
						
		if oldlength == len(node):
			converged = True
				
	return node


def countKmersOfCountOne(kmers, k):
	count = 0
	for kmer in kmers:
		if kmers[kmer] == 1:
			count +=1
	
	print 'Found ' + str(count) + ' ' +str(k) + '-mers that were present once' 
	
	
def assemble(kmers, k, maxrepeats):
	finals = []
	total = len(kmers)
	count = 0
	lastp = -10
	percentage = 0
	dispatched = dict()
	for kmer in kmers:
		if(kmer not in dispatched):
			f = getneighbor(kmer, kmerlist, k, dispatched, maxrepeats)
			if f not in finals:
				finals.append(f)
		count += 1
		percentage = float(count)/float(total) * 100
		if abs(percentage - lastp) >= 10:
			lastp = percentage
			print "{0:.0f}%".format(lastp)
			
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

k = 23
errorcutoff = 2
maxrepeats = 1
assemb = []

addreverse(reads)

errorcorrect(reads, errorcutoff, k)

kmerlist = computekmers(k, reads)

print 'Generated ' + str(len(kmerlist)) + ' many ' + str(k) + '-mers'
countKmersOfCountOne(kmerlist, k)

if len(kmerlist) == 0:
	exit()

finals = assemble(kmerlist, k, maxrepeats)

with open(str(sys.argv[2]), "w") as out:
	count = 0
	for f in finals:
		count += 1
		out.write('> contig '+ str(count)+'\n')
		out.write(str(f)+"\n")
