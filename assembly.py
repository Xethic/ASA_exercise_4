import sys
from string import maketrans

#function for compution of the reverse complement of a read
alphabet = ['A', 'C', 'G', 'T']
original = 'ACGT'
complement = 'TGCA'
revtable = maketrans(original, complement)


def getComplement(sequence):
	return sequence.translate(revtable)


#addition of the reverse complement of a read to the final list of reads
def addreverse(reads):
	revreads = []
	for read in reads:
		revreads.append(getComplement(read)[::-1])
	reads.extend(revreads)


#computation of kmers of a read
#stored in dictionary, kmer as key, coverage as value
def computekmers(k, reads):
	kmers = dict()
	for r in reads:
		for i in range(len(r)-k+1):
			if r[i:i+k] in kmers:
				kmers[r[i:i+k]] += 1
			else:
				kmers[r[i:i+k]] = 1
	return kmers
	
	
#computation of the most similar sequence using hamming distance
#the cutoff gives the maximum distance of the two sequences that is allowed
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
	
	
#correction of the putative read errors
#for each kmer with coverage 1, search for a sequence with hamming distance < cutoff
#if there exist one with a higher coverage then the kmer is correcet by the better one 
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
			
	
#searching for extensions of kmers
#if there are more than one possible extension, the contig breaks and starts with each branch again as a start node
#the dictionary dispatched stores the kmers that are used at most 5 times, otherwise the kmers are not used again for extension
#returns the extended kmer
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
				if s2+c not in dispatched or dispatched[s2+c] < 6:
					rightextensions.append(c)
		
		if(len(rightextensions) == 1):
			c = rightextensions[0]
			if s2+c in dispatched:
				dispatched[s2+c] += 1
			else:
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

#counting of the kmers that occur only once
def countKmersOfCountOne(kmers, k):
	count = 0
	for kmer in kmers:
		if kmers[kmer] == 1:
			count +=1
	
	print 'Found ' + str(count) + ' ' +str(k) + '-mers that were present once'


#searching for nodes that can only be extended on the right site
#stored as start nodes of finding neigbors
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


#for all startnode extensions are searched
#if there is no more extension, the contigs are stored in a final list
def assemble(kmers, k):
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
	if len(sys.argv) > 3:
		print "Too many arguments, only two allowed!"
		exit()
		
	reads = []
	readset = set()

	#reading of the input file of reads
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
		
		print 'Parsed ' + str(count) + ' reads'
		print 'Found ' + str(duplicated) +' duplicated reads'

	#used k for the kmers
	k = 25
	#cutoff for hamming distance
	errorcutoff = 2
	
	assemb = []

	addreverse(reads)

	errorcorrect(reads, errorcutoff, k)

	kmerlist = computekmers(k, reads)

	print 'Generated ' + str(len(kmerlist.keys())) + ' many ' + str(k) + '-mers'
	countKmersOfCountOne(kmerlist, k)

	if len(kmerlist) == 0:
		exit()

	finals = assemble(kmerlist, k)

	#writing of the output file in a fasta format 
	with open(str(sys.argv[2]), "w") as out:
		count = 0
		for f in finals:
			count += 1
			out.write('> contig '+ str(count)+'\n')
			out.write(str(f)+"\n")
		
main()
