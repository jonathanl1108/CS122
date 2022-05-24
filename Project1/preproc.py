from time import clock
import cPickle

def read_genome(filename):
	"""Read in genome file to string """
	file = open(filename,'r')
	first_line = True
	genome = ''
	for line in file:
		if(first_line):
			first_line = False
			continue
		genome = genome + line.strip()
	file.close()
	return genome

def index_genome(genome, read_length = 50, mismatches = 4):
	"""Build genome index/hash table using dictionary implementation"""
	fragments = mismatches + 1
	k = int(round(read_length/fragments))
	dict = {}
	for i in range(len(genome)-k):
		kmer = genome[i:i+k]
		if kmer in dict:
			dict[kmer].append(i)
		else:
			dict[kmer] = [i]
	return dict

def serialize(genome, filename = 'genome_index.txt'):
	"""Serialize genome index for future use in .txt file"""
	ser = cPickle.dumps(genome)
	f = open(filename,'w')
	f.write(ser)
	f.close()

def deserialize(file = 'genome_index.txt'):
	"""Deserialize genome index from .txt file"""
	f = open(file,'r')
	ser = f.read().strip()
	return cPickle.loads(ser)

"""
start1 = clock()
genome = read_genome('ref_hw1_W_2_chr_1.txt')
end1 = clock()
print "Genome reading: {} ms".format(end1-start1)

start2 = clock()
ind = index_genome(genome)
end2 = clock()
print "Indexing: {} ms".format(end2-start1)

start3 = clock()
serialize(ind)
end3 = clock()
print "Serializing: {} ms".format(end3-start3)

start4 = clock()
genome2 = deserialize()
end4 = clock()
print "Deserialzing: {} ms".format(end4-start4)

f = open('index.txt','w')
for key in genome2.keys():
	f.write("{} : {} \n".format(key,genome2[key]))
f.close()
"""
