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

	
def find(kmerlist):	
	assemb = []
	
	for key in kmerlist:
		s = ""
		found = True		
		query = []		
		start = key
			 
		s += start
		
		query.append(start)
		
		f = 0
		
		while len(query) != 0:
			if not found:
				if s not in assemb:
					assemb.append(s)
				s = query[-1]		
			
			first = query[-1][1:]
			print first
					
			f = 0
			s1 = first + "A"
			s2 = first + "C"
			s3 = first + "G"
			s4 = first + "T"
			
			if s1 in kmerlist:
				s += "A"
				f += 1
				found = True
				query += [s1]
			else:
				found = False
			if s2 in kmerlist:
				if f == 0: 
					s += "C"
					f += 1
					found = True
					query += [s2]
				else:
					query += [s2]
					found = True
			else:
				if f == 0:
					found = False
			if s3 in kmerlist:
				if f == 0: 
					s += "G"
					f += 1
					found = True
					query += [s3]
				else:
					query += [s3]
					found = True
			else:
				if f == 0:
					found = False
			if s4 in kmerlist:
				if f == 0: 
					s += "T"
					f += 1
					found = True
					query += [s4]
				else:
					query += [s4]
					found = True
			else:
				if f == 0:
					found = False
						
			query.pop(0)
			
	
		if s not in assemb:		
			assemb.append(s)	
	

	return assemb
	
	
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

k = 3

kmerlist = kmers(k, reads)
print kmerlist

seq = find(kmerlist)
print seq
