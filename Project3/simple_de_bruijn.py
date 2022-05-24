from os.path import join
import sys
import time
from collections import defaultdict, Counter
import zipfile
import logging as logger
from random import choice

logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

file_sets = ['practice_W_3_chr_1','practice_E_1_chr_1','hw2undergrad_E_2_chr_1','hw3all_A_3_chr_1']

def read_reads(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        read = line.split(',')  # Only take the first read.
        #possible = line.split(',')
        #possible1 = [possible[0],possible[0][::-1]]
        #possible2 = [possible[1],possible[1][::-1]]
        #read = choice(possible1)
        #read2 = choice(possible2)
        # Clearly, there is room for improvement here.
        all_reads.append(read[0])
        all_reads.append(read[1])
    return all_reads


def simple_de_bruijn(sequence_reads, k):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the set
                of k-mers that
    """
    de_bruijn_counter = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read in sequence_reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 2}
                       for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}

    out_degree = defaultdict(int)
    for key in de_bruijn_counter:
        for k in de_bruijn_counter[key]:
            out_degree[key] += de_bruijn_counter[key][k]

    in_degree = defaultdict(int)
    for value in de_bruijn_counter.values():
        for kmer in value:
            in_degree[kmer] += value[kmer]

    return de_bruijn_graph, in_degree, out_degree


def de_bruijn_reassemble(de_bruijn_graph, in_degree, out_degree):
    """
    Traverses the DeBruijn Graph created by simple_de_bruijn and
    returns contigs that come from it.
    :param de_bruijn_graph: A De Bruijn Graph
    :return: a list of the
    """

    assembled_strings = []
    #counter = 0
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        #print n_values
        if n_values == 0:
            break
        #good_starts = [k for k in de_bruijn_graph if in_degree[k] == 0]
        good_starts = [k for k in de_bruijn_graph if out_degree[k] > in_degree[k] and de_bruijn_graph[k]]
        # You may want to find a better start
        # position by looking at in and out-degrees,
        # but this will work okay.
        if len(good_starts) == 0:
            good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]
        current_point = good_starts[0]
        assembled_string = current_point
        while True:
            try:
                next_values = de_bruijn_graph[current_point]
                next_edge = next_values.pop()
                #counter += 1
                #print "Reads poppped:", counter
                assembled_string += next_edge[-1]
                de_bruijn_graph[current_point] = next_values
                current_point = next_edge
            except KeyError:
                assembled_strings.append(assembled_string)
                break
    return assembled_strings


if __name__ == "__main__":
    logger.warn('START:')
    chr_name = file_sets[1]
    input_folder = './{}'.format(chr_name)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(chr_name))
    reads = read_reads(reads_fn)
    logger.warn('DEBRUIJN:')
    db_graph, in_degs, out_degs = simple_de_bruijn(reads, 35)
    # for k in db_graph.keys()[:40]:
    #     print k, db_graph[k]
    logger.warn('ASSEMBLE:')
    output = de_bruijn_reassemble(db_graph,in_degs,out_degs)
    output_fn_end = 'assembled_{}.txt'.format(chr_name)
    output_fn = join(input_folder, output_fn_end)
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))
    zip_fn = join(input_folder, 'assembled_{}.zip'.format(chr_name))
    with zipfile.ZipFile(zip_fn,'w') as myzip:
        myzip.write(output_fn)
    logger.warn('END:')
