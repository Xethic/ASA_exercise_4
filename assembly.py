import sys


def kmers(k, reads):
	kmer = dict()
		
	for r in reads:
		for i in range(len(r)-k+1):
			if r[i:i+k] in kmer:
				kmer[r[i:i+k]] += 1
			else:
				kmer[r[i:i+k]] = 1
			
	return kmer


def find(k, kmerlist, longest):
	print longest
	s = k
	
	if s[1:] == longest.items()[0][0][:3]:
		longest[s] = s[1:]+longest.items()[0][0][:3]
	
	first = k
	
	for k2 in kmerlist:
			
		s1 = first[1:] + "A"
		s2 = first[1:] + "C"
		s3 = first[1:] + "G"
		s4 = first[1:] + "T"
		
		if k2 in longest:
			s += longest[k2]
		if k2 == s1:
			s += "A"
			first = s1
			found = True
		elif k2 == s2:
			s += "C"
			first = s2
			found = True
		elif k2 == s3:
			s += "G"
			first = s3
			found = True
		elif k2 == s4:
			s += "T"
			first = s4
			found = True
		else:
			found = False
		
		if not found:
			return s
				
	
	'''
def find(kmerlist):	
	assemb = []
	startnode = ""
	

	
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
	
	'''
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

k = 4

longest = dict()
assemb = []

kmerlist = kmers(k, reads)
print kmerlist

for km in kmerlist:
	assemb.append(find(km, kmerlist, longest))
	
print assemb
