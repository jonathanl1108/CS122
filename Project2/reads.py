def reads(filename = 'reads_hw1_W_2_chr_1.txt'):
	"""Read in reads as pair-ended read"""
	file = open(filename,'r')
	first_line = True
	reads = []
	for line in file:
		if(first_line):
			first_line = False
			continue
		reads.append(line.strip())
	file.close()
	return reads

def split_pair(pair_read):
	"""Split pair-ended read """
	return pair_read.split(',')

def kmer_read(read, k = 10):
	"""Split read into k-mers"""
	kmers = []
	for i in range(0,len(read),k):
		kmers.append(read[i:i+k])
	return kmers

def generate_single_reads(filename = 'reads_hw1_W_2_chr_1.txt'):
	"""Generate single-end reads from file"""
	paired = reads(filename)
	single_reads = []
	for read in paired:
		single = split_pair(read)
		single_reads.extend(single)
	return single_reads

def generate_pair_reads(filename):
	paired = reads(filename)
	return paired

def HammingDistance(read1,read2):
	"""Calculate Hamming Distance between two n-length strings (mismatches)"""
	assert len(read1) == len(read2), "HammingDistance Error: Not Same Length: {} vs. {}".format(len(read1),len(read2))
	count = 0
	for i in range(len(read1)):
		if (read1[i] != read2[i]):
			count += 1
	return count

"""
r = generate_single_reads()
for item in r:
	print item
	kmers = kmer_read(item)
	for fragment in kmers:
		print "\t{}".format(fragment)
"""

#print HammingDistance('ATT','ATT')