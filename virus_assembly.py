### virus_assembly.py
# Given several different viral DNA strains and a set of reads, calculates the relative frequencies of each strain using least squares


from os.path import join
import numpy as np
from collections import Counter, defaultdict


# Skips first line of file and reads in remainder of lines into a list
def read_data(read_fn):
    '''
    :param read_fn: Name of file
    :return reads: Returns each line of the file in a list
    '''
    reads = []
    with open(read_fn, 'r') as f:
        next(f)
        for line in f:
            line = line.strip()
            reads.append(line)
    return reads


# Creates the SNP matrix of 1s and 0s from a list of strains. The reference sequence is taken to be the consensus sequence of all the strains.
def generate_snp_matrix_from_strains(strains):
    '''
    :param strains: List of strain sequences
    :return reference: The consensus sequence
    :return snp_matrix: The SNP matrix
    :return snp_inds: The indices of the SNPs in the reference
    '''
    snp_inds = []
    snp_matrix = []
    reference = ''
    for ind, base in enumerate(zip(*strains)):
        base = Counter(base)
        cons = base.most_common(1)
        if cons[0][1] < len(strains):
            snp_inds.append(ind)
        reference += cons[0][0]
    for strain in strains:
        snp_row = [0 if reference[j] == strain[j] else 1 for j in snp_inds]
        snp_matrix.append(snp_row)
    return reference, snp_matrix, snp_inds


# Builds a hash table of indices from the reference
def build_index(ref, length):
    '''
    :param ref: The reference sequence
    :param length: The length of the keys
    :return index: The hash table of indices
    '''
    index = defaultdict(list)
    for i in range(0, len(ref)-length):
        subset = ref[i:i+length]
        index[subset].append(i)
    return index


# Maps a read to the reference
def compare_reads(ref, read, index, length_index, max_mismatches, num_pieces):
    '''
    :param ref: The reference sequence
    :param read: The read to be mapped
    :param index: The hashed indices generated from build_index
    :param length_index: The length of the keys in the hashed indices
    :param max_mismatches: The maximum number of mismatches to be allowed
    :param num_pieces: The number of pieces the read is divided into so that the length of each piece equals the length of the index keys
    :return start: Returns the index of the match if there is a match, otherwise returns None
    '''
    for sub in range(num_pieces):
        subsection = read[sub * length_index:(sub + 1) * length_index]
        if subsection in index:
            for ind in index[subsection]:
                start = ind - sub * length_index
                if start + len(read) <= len(ref):
                    mismatches = [1 if read[j] != ref[start + j] else 0 for j in range(len(read))]
                    n_mismatches = sum(mismatches)
                    if n_mismatches <= max_mismatches:
                        return start
    return None


# Maps all reads to a reference and generates pileup
def alignment_algorithm(reads, ref, max_mismatches, num_pieces):
    '''
    :param reads: List of reads
    :param ref: Reference sequence
    :param max_mismatches: Maximum amount of mismatches betwen read and reference allowed
    :param num_pieces: Number of subsections read is divided to when looking up in the hashed indices
    :return pileup: Returns a dictionary containing the total number of each base at each index in the reference
    '''
    length_read = len(reads[0])
    length_index = int(length_read/num_pieces)
    index = build_index(ref, length_index)
    pileup = defaultdict(Counter)

    for read in reads:
        ind = compare_reads(ref, read, index, length_index, max_mismatches, num_pieces)
        if ind is None:
            rev_read = read[::-1]
            rev_ind = compare_reads(ref, rev_read, index, length_index, max_mismatches, num_pieces)
            if rev_ind is not None:
                for i, base in enumerate(rev_read):
                    pileup[rev_ind + i].update(rev_read[i])
        else:
            for i, base in enumerate(read):
                pileup[ind + i].update(read[i])

    return pileup


# Calculates strain frequencies
def calculate_frequencies(pileup, snp_inds, snp_matrix, reference):
    '''
    :param pileup: The pileup generated from alignment_algorithm
    :param snp_inds: The indices of the SNPs in the reference
    :param snp_matrix: The SNP matrix
    :param reference: The reference sequence
    :return strains_freqs: Returns a list of strain frequencies
    '''
    snp_freqs = []
    for i in snp_inds:
        base = reference[i]
        variant = next(x for x in list(pileup[i].elements()) if x != base)
        variant_freq = float(pileup[i][variant])/(pileup[i][variant] + pileup[i][base])
        snp_freqs.append(variant_freq)
    snp_freqs = np.array(snp_freqs)
    snp_matrix = np.array(snp_matrix)
    strain_freqs = np.linalg.lstsq(snp_matrix.T, snp_freqs)
    return strain_freqs[0]


if __name__ == '__main__':
    data_folder = 'hw4_W_1'
    input_folder = join('/Users/douglasyao/Downloads/', data_folder)
    reads_fn = join(input_folder, '{}_reads.txt'.format(data_folder))
    strains_fn = join(input_folder, '{}_strains.txt'.format(data_folder))
    output_fn = join(input_folder, '{}_ans.txt'.format(data_folder))

    reads = read_data(reads_fn)
    strains = read_data(strains_fn)
    reference, snp_matrix, snp_inds = generate_snp_matrix_from_strains(strains)
    pileup = alignment_algorithm(reads, reference, 2, 2)
    frequencies = calculate_frequencies(pileup, snp_inds, snp_matrix, reference)
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + data_folder + '\n')
        for i in range(len(frequencies)):
            output_file.write('{},{}\n'.format(frequencies[i], strains[i]))


