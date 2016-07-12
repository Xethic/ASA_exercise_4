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
	kmers = dict()
	for r in reads:
		for i in range(len(r)-k+1):
			if r[i:i+k] in kmers:
				kmers[r[i:i+k]] += 1
			else:
				kmers[r[i:i+k]] = 1
	return kmers
	
	
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
			
	
def getneighbor(node, kmers, k, startnodes):
	dispatched = dict()
	converged = False
	dispatched[node] = 1
	breakright = False
	
	while not converged:
		s2 = node[-(k-1):]
		oldlength = len(node)
		rightextensions = []
		for c in ['A', 'C', 'G', 'T']:
			if (not breakright) and (s2+c in kmers):
				if s2+c not in dispatched :
					rightextensions.append(c)
		
		if(len(rightextensions) == 1):
			c = rightextensions[0]
			dispatched[s2+c] = 1
			node += c
		elif len(rightextensions) > 1:
			for c in rightextensions:
				if s2+c not in startnodes:
					startnodes.append(s2+c)
			breakright = True
						
		if oldlength == len(node):
			converged = True
				
	return node


def countKmersOfCountOne(kmers, k):
	count = 0
	for kmer in kmers:
		if kmers[kmer] == 1:
			count +=1
	
	print 'Found ' + str(count) + ' ' +str(k) + '-mers that were present once'


def searchforstartnodes(kmers,k):
	startnodes = []
	for kmer in kmers:
		s1 = kmer[:k-1]
		found = False
		for c in ['A', 'C', 'G', 'T']:
			if c+s1 in kmers:
				found = True
				break
		
		if not found:
			startnodes.append(kmer)
	return startnodes

	
def assemble(kmers, k, maxrepeats):
	finals = []
	count = 0
	lastp = -10
	percentage = 0
	
	startnodes = searchforstartnodes(kmers, k)
	total = len(startnodes)
	print 'Found ' + str(len(startnodes)) + ' start-nodes'
	for kmer in startnodes:
		total = len(startnodes)
		f = getneighbor(kmer, kmers, k, startnodes)
		if f not in finals:
			finals.append(f)
		count += 1
		percentage = float(count)/float(total) * 100
		if abs(percentage - lastp) >= 1:
			lastp = percentage
			print "{0:.0f}%".format(lastp)
			
	return finals


#######
#	MAIN
#
#	first argument is the input file*
#	second argument is the name of the output file
#######
def main():

	reads = []
	readset = set()

	with open(str(sys.argv[1])) as f:
		duplicated = 0
		count = 0
		for line in f:
			if line == "\n":
				break
				
			if line[0] == ">":
				continue
			else:
				read = line.replace('\n', '')
				reads.append(read)
				count += 1
				if read in readset:
					duplicated += 1
				else:
					readset.add(read)
		#reads = list(readset)	
		print 'Parsed ' + str(count) + ' reads'
		print 'Found ' + str(duplicated) +' duplicated reads'

	k = 21
	errorcutoff = 2
	maxrepeats = 1
	assemb = []

	addreverse(reads)

	errorcorrect(reads, errorcutoff, k)

	kmerlist = computekmers(k, reads)

	print 'Generated ' + str(len(kmerlist.keys())) + ' many ' + str(k) + '-mers'
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
		
main()
