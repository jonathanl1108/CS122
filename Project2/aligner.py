from output import *
from preproc import *
from reads import *
from time import clock
import concurrent.futures
from multiprocessing import Manager
import logging as logger

#speed ups = larger key length
#don't check both of the read pair in reverse (do 3 checks instead of 4 per read pair)
	#use boolean to check which is in reverse
#double loop when checking for HammingDistance can be reduced to single loop
	#keep track of SNP's while doing first loop for mistmatch count
#read in less reads at a time (MICHAEL SAYS USE ALL READS)

logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

start = clock()

reference = read_genome('ref_hw2undergrad_E_2_chr_1.txt')

ref_index = index_genome(reference,50,4)
#serialize(ref_index)
#ref_index = deserialize('genome_index.txt')

single_reads = generate_single_reads('reads_hw2undergrad_E_2_chr_1.txt')

for each in single_reads[:]:
	single_reads.append(each[::-1])

#account for reverse order reads from pair ending

# single_reads = []

# pairs = generate_pair_reads('reads_hw1_W_2_chr_1.txt')

# logger.warn("Read pairs")
# for each in pairs:
# 	split = each.split(',')
# 	for i in range(len(split)):
# 		if (split[0] in ref_index):
# 			single_reads.append(split[0])
# 			single_reads.append(split[1][::-1])
# 		else:
# 			single_reads.append(split[1])
# 			single_reads.append(split[0][::-1])

mismatches = 4

manager = Manager()

SNP_candidates = manager.dict() #stores tuple (Reference,Variant,Position)

def align(read):
	kmers = [read[:10]]
	#kmers = kmer_read(read, 10) #split read into k-mer
	for each in kmers:
		#if(each in ref_index.keys()):
		if (ref_index.has_key(each)): #if k-mer found in reference index
			#positions = ref_index[each] #positions equal to reference index
			for pos in ref_index[each]:
				index = read.find(each) #find position of k-mer in read to align middle/end k-mers
				if(pos-index < 0 or pos-index+len(read) > len(reference)):
					continue #skip if read is out of bounds in reference
				#differences = HammingDistance(read, reference[pos-index:pos-index+len(read)])
				differences = 0
				ls = []
				for i in range(len(read)):
					if (read[i] != reference[pos-index+i]): #if mismatch record tuple
						s = (reference[pos-index+i],read[i],pos-index+i)
						differences += 1
						ls.append(s)
				if (differences > mismatches):
					continue
				else:
					for s in ls:
						if SNP_candidates.has_key(s):
						#if s in SNP_candidates: #increase count
							SNP_candidates[s] += 1
						else:
							SNP_candidates[s] = 1
				# if (differences > mismatches):
				# 	continue #discard if differences are greater than our allowed mismatches
				# else:
				# 	for i in range(len(read)):
				# 		if (read[i] != reference[pos-index+i]): #if mismatch record tuple
				# 			s = (reference[pos-index+i],read[i],pos-index+i)
				# 			if s in SNP_candidates: #increase count
				# 				SNP_candidates[s] += 1
				# 			else:
				# 				SNP_candidates[s] = 1

logger.warn("multiprocessing:")
executor = concurrent.futures.ProcessPoolExecutor(10)
futures = [executor.submit(align,read) for read in single_reads]
concurrent.futures.wait(futures)

list_of_SNPS = []

logger.warn("Filtering:")
for key in SNP_candidates.keys():
	if (SNP_candidates[key] <= 10):
		del SNP_candidates[key] #remove SNP if not >= 90%
	else:
		s = "{},{},{}".format(key[0],key[1],key[2])
		list_of_SNPS.append(s)

print "SNP's Found:", len(SNP_candidates)

#for key in SNP_candidates.keys():
	# s = "{},{},{}".format(key[0],key[1],key[2])
	# list_of_SNPS.append(s)

generate_file(header='hw2undergrad_E_2_chr_1',SNP=list_of_SNPS)
end = clock()
logger.warn("End")
print "Time:",end-start