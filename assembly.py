import sys


def kmers(k, reads):
	kmer = dict()
	
	for r in reads:
		for i in range(len(r)-k):
			if r[i:i+k+1] in kmer:
				kmer[r[i:i+k]] += 1
			else:
				kmer[r[i:i+k]] = 1
			
	return kmer


def graph(kmers):
	g = []

	for k in kmers:
		for j in kmers:
			if k[1:] == j[:-1]:
				g.append((k,j))		
	return g
		

def route(g):
	
	seq = g[0][0]
	
	first = g[0][0]
	
	for i in range(len(g)):		
		for node in g:
			if node[1] == first:
				seq += node[1][:-1]
				first = node[0]
				g.remove(node)
	
				
	return seq
		
		
	
	
	
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
			reads.append(line[:-1])

print reads

kmerlist = kmers(8, reads)
	
g = graph(kmerlist)

print g

print route(g)
