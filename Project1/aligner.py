from output import *
from preproc import *
from reads import *

reference = read_genome('ref_hw1_W_2_chr_1.txt')

ref_index = index_genome(reference,50,4)
serialize(ref_index)

single_reads = generate_single_reads('reads_hw1_W_2_chr_1.txt')

#account for reverse order reads from pair ending
for each in single_reads[:]:
	single_reads.append(each[::-1])

mismatches = 4

SNP_candidates = {} #stores tuple (Reference,Variant,Position)

for read in single_reads:
	kmers = kmer_read(read, 10) #split read into k-mer
	for each in kmers:
		if each in ref_index.keys(): #if k-mer found in reference index
			positions = ref_index[each] #positions equal to reference index
			for pos in positions:
				index = read.find(each) #find position of k-mer in read to align middle/end k-mers
				if(pos-index < 0 or pos-index+len(read) > len(reference)):
					continue #skip if read is out of bounds in reference
				differences = HammingDistance(read, reference[pos-index:pos-index+len(read)])
				if (differences > mismatches):
					continue #discard if differences are greater than our allowed mismatches
				else:
					for i in range(len(read)):
						if (read[i] != reference[pos-index+i]): #if mismatch record tuple
							s = (reference[pos-index+i],read[i],pos-index+i)
							if s in SNP_candidates: #increase count
								SNP_candidates[s] += 1
							else:
								SNP_candidates[s] = 1

for key in SNP_candidates.keys():
	if (SNP_candidates[key] <= 27):
		del SNP_candidates[key] #remove SNP if not >= 90%

print "SNP's Found:", len(SNP_candidates)

list_of_SNPS = []

for key in SNP_candidates.keys():
	s = "{},{},{}".format(key[0],key[1],key[2])
	list_of_SNPS.append(s)

generate_file(SNP=list_of_SNPS)