### find_copy_number.py
# Finds copy number of donor sequence from reference sequence and reads


import numpy as np
from os.path import join
from collections import defaultdict, OrderedDict, Counter
from itertools import groupby, combinations
import zipfile
import re
from multiple_sequence_aligner_undergrad_assignment import read_reads, read_reference, build_index
from math import ceil
from find_strs import read_ans


# Given a dict of coverage of each section of reference aligned to itself, return start of sections where coverage exceeds 1x
def collapse_dict(all_starts):
    '''
    :param all_starts: Dict of reference coverage
    :return: locations of sections
    '''
    dups = []
    iter_starts = all_starts.iteritems()
    prev, b = next(iter_starts)
    for k, v in iter_starts:
        if v != 1:
            if k - prev > 10:
                dups.append(k)
            prev = k
    return dups


# Finds copy number in reference
def find_copy_number_reference(reference, index):
    '''
    :param reference: Reference sequence
    :param index: Indexed reference sequence as a dict
    :return: List of copy number sequences and their locations in the reference
    '''
    all_starts = OrderedDict([(10*i, 0) for i in range(1000)])
    seqs = []
    for i in range(1000):
        subsection = reference[i * 10:(i + 1) * 10]
        if subsection in index:
            start = round_down(index[subsection][0], 10)
            all_starts[start] += 1
    dups = collapse_dict(all_starts)
    for i in dups:
        sec = reference[i:i+10]
        forward = 0
        backward = 0
        matches = re.findall(sec, reference)
        if len(matches) > 1:
            forward += 1
            sec = reference[i:i+10+forward]
            new_match = re.findall(sec, reference)
            while len(new_match) == len(matches):
                forward += 1
                sec = reference[i:i+10+forward]
                new_match = re.findall(sec, reference)
            backward += 1
            sec = reference[i-backward:i+9+forward]
            new_match = re.findall(sec, reference)
            while len(new_match) == len(matches):
                backward += 1
                sec = reference[i-backward:i+9+forward]
                new_match = re.findall(sec, reference)
            sec = sec[1:]
        indices = [x.start() for x in re.finditer(sec, reference)]
        seqs.append([sec, indices])
    seqs = [[sec, indices] for sec, indices in seqs if len(sec) > 15 and len(indices) > 1]
    seqs.sort()
    seqs = [seqs for seqs, x in groupby(seqs)]
    if len(seqs) > 1:
        for combo in combinations([x[0] for x in seqs], len(seqs)):
            if combo[0] in combo[1] or combo[1] in combo[0]:
                to_compare = [x for x in seqs if x[0] in combo]
                to_remove = min(to_compare, key = lambda x: len(x[1]))
                seqs.remove(to_remove)
    return seqs


# Round a number down to the nearest multiple of a given number
def round_down(num, divisor):
    return num - (num%divisor)


# Aligns reads to reference and identifies regions where coverage is some multiple of normal coverage.
# Perform assembly at copy number altered regions to identify copy number in donor sequences
def calculate_coverage(paired_reads, index, length_index, num_pieces, reference_copy_number, reference):
    '''
    :param paired_reads: List of paired end reads
    :param index: Indexed reference sequence as a dict
    :param length_index: Length of each key in the reference index
    :param num_pieces: Number of pieces to divide each read into
    :param reference_copy_number: List of reference copy number sequences and their locations in the reference
    :param reference: The reference sequence
    :return new_insertions: Locations in reference where a copy is inserted
    :return deletions: Locations in reference where a copy is deleted
    '''
    nreads = [item for sublist in paired_reads for item in sublist]
    insertions = defaultdict(set)
    deletions = []
    all_starts = OrderedDict([(10*i, 0) for i in range(int(ceil(len(reference)/10.0)))])
    for read_pair in paired_reads:
        for read in read_pair:
            ind = compare_reads_coverage(read, index, length_index, num_pieces)
            if ind.count(None) > 2:
                rev_read = read[::-1]
                rev_ind = compare_reads_coverage(rev_read, index, length_index, num_pieces)
                if rev_ind.count(None) < 3:
                    rev_ind = [round_down(x, 10) for x in rev_ind if x is not None]
                    for i in rev_ind:
                        all_starts[i] += 1
            else:
                ind = [round_down(x, 10) for x in ind if x is not None]
                for i in ind:
                    all_starts[i] += 1
    ave_cov = np.mean([all_starts.values()])
    locs = OrderedDict([(k, v) for k, v in all_starts.iteritems() if v > 1.5*ave_cov or v < ave_cov/4])
    locs = collapse_dict(locs)

    # Perform assembly
    for i in locs:
        print i
        begin = -20
        matching_reads = assemble_reads(i, reference, begin, nreads)
        if matching_reads:
            assembled_read = debruijn_assemble(matching_reads)
            if assembled_read is not None:
                match_ref = compare_cnv(reference_copy_number, assembled_read)
                while assembled_read[:-20] in reference and match_ref is None:
                    begin += 10
                    matching_reads = assemble_reads(i, reference, begin, nreads)
                    if not matching_reads:
                        break
                    else:
                        assembled_read = debruijn_assemble(matching_reads)
                        match_ref = compare_cnv(reference_copy_number, assembled_read)
                if match_ref is not None and assembled_read[:15] in reference:
                    read_start = re.search(assembled_read[:15], reference).start()
                    cnv_start = re.search(match_ref[:10], assembled_read).start()
                    cnv_insertion = read_start + cnv_start
                    insertions[match_ref].add(cnv_insertion)
                else:
                    add = 0
                    seq = assembled_read[:10]
                    while seq in reference:
                        add += 1
                        seq = assembled_read[:10+add]
                        if 10+add > len(seq):
                            break
                    seq = seq[:-1]
                    if len(seq) > 15:
                        pos = add + re.search(seq, reference).start() + 9
                        remainder = assembled_read[9+add:]
                        if remainder in reference:
                            to_add = compare_cnv_indices(reference_copy_number, pos)
                            if to_add is not None:
                                deletions.append(to_add)
    new_insertions = defaultdict(list)
    ref_dict = OrderedDict(reference_copy_number)
    for k, v in insertions.iteritems():
        for ind in v:
            if ind not in ref_dict[k]:
                new_insertions[k].append(ind)
    return new_insertions, deletions


# At a given location in reference, find reads that match to region around location
def assemble_reads(loc, reference, start, reads):
    '''
    :param loc: Location in reference
    :param reference: Reference sequence
    :param start: Where reads should map relative to location
    :param reads: List of paired end reads
    :return: List of reads that match to location
    '''
    match = reference[loc + start:loc + start + 15]
    matching_reads = [x for x in reads if match in x]
    matching_reads.extend([x[::-1] for x in reads if match in x[::-1]])
    return matching_reads


# Assemble reads using a De Bruijn graph
def debruijn_assemble(reads):
    '''
    :param reads: List of reads to be assembled
    :return: Assembed sequence
    '''
    debruijn = simple_de_bruijn(reads, 25)
    debruijn = {k: list(v)[0] for k, v in debruijn.iteritems()}
    if len(debruijn.keys()) == 0:
        return None
    else:
        start = [x for x in debruijn.keys() if x not in debruijn.values()][0]
        current_point = start
        while True:
            try:
                next_values = debruijn[current_point]
                next_edge = next_values
                start += next_edge[-1]
                current_point = next_edge
            except KeyError:
                break
        return start


# Detects whether any reference copy number sequences overlap with given sequence
def compare_cnv(reference_copy_number, cmp):
    '''
    :param reference_copy_number: List of reference copy number sequences and locations
    :param cmp: Sequence to compare with
    :return: Reference sequence if it overlaps with the given sequence, otherwise None
    '''
    ref = [x[0] for x in reference_copy_number]
    for loc in ref:
        if loc[:10] in cmp:
            return loc
    return None


# Detects whether a given location overlaps with reference copy number locations
def compare_cnv_indices(reference_copy_number, index):
    '''
    :param reference_copy_number: List of reference copy number sequences and locations
    :param index: Location to be compared
    :return: Reference sequence index if it is within 3 bp of the given index, otherwise None
    '''
    indices = [x[1] for x in reference_copy_number]
    indices = [x for y in indices for x in y]
    for ind in indices:
        if abs(index - ind) < 3:
            return ind
    return None


# Creates a De Bruijn graph from a set of reads
def simple_de_bruijn(sequence_reads, k):
    '''
    :param sequence_reads: List of reads from which to generate De Bruijn graph
    :param k: length of k-mer to generate read spectrum
    :return: the De Bruijn graph
    '''
    de_bruijn_counter = defaultdict(Counter)
    for read in sequence_reads:
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 5} for key in de_bruijn_counter}
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    return de_bruijn_graph


# Find locations where read matches to indexed reference
def compare_reads_coverage(read, index, length_index, num_pieces):
    '''
    :param read: Read sequence
    :param index: Indexed reference as a dict
    :param length_index: Length of keys in indexed reference
    :param num_pieces: Number of pieces to divide read into
    :return: All positions where read matches
    '''
    positions = []
    for sub in range(num_pieces):
        subsection = read[sub * length_index:(sub + 1) * length_index]
        if subsection in index:
            positions.append(index[subsection][0])
        else:
            positions.append(None)
    return positions


# Reads in existing answer file and modifies accordingly to add CNVs. CNVs are also added to insertions and deletions where appropriate.
def prepare_output(ans_fn, cnv, new_ins, new_dels):
    '''
    :param ans_fn: File name of answer file
    :param cnv: List of copy number sequences and locations
    :param new_ins: List of copy number insertions
    :param new_dels: List of copy number deletions
    :return: Updated SNPs, insertions, deletions, and CNVs
    '''
    new_cnv = defaultdict(list)
    final_cnv = []
    all_dels = defaultdict(list)
    for k, v in cnv:
        for ind in v:
            if ind in new_dels:
                all_dels[k].append(ind)
                continue
            new_cnv[k].append(ind)
    snps, ins, dels, cnvs, strs = read_ans(ans_fn)
    for k, v in new_ins.iteritems():
        for ind in v:
            ins.append([k, str(ind)])
        val = new_cnv[k]
        val.extend(v)
        val = sorted(val)
        new_cnv[k] = val
    if len(all_dels.keys()) > 0:
        for k, v in all_dels.iteritems():
            for ind in v:
                dels.append([k, str(ind)])
    ins.sort()
    ins = [ins for ins, x in groupby(ins)]
    ins.sort(key=lambda x: int(x[1]))
    dels.sort()
    dels = [dels for dels, x in groupby(dels)]
    dels.sort(key=lambda x: int(x[1]))
    for k, v in new_cnv.iteritems():
        to_add = [k]
        to_add.extend(v)
        final_cnv.append(to_add)
    final_cnv.sort(key=lambda x:x[1])
    return snps, ins, dels, final_cnv


if __name__ == "__main__":

    data_folder = 'hw1_W_2'
    input_folder = '/Users/douglasyao/Downloads/' + data_folder
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    input_reads = read_reads(reads_fn)
    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    ans_fn = join(input_folder, 'snps_{}.txt'.format(f_base))
    index = build_index(reference, 10)
    cn = find_copy_number_reference(reference, index)
    ins, dels = calculate_coverage(input_reads, index, 10, 5, cn, reference)
    snps, ins, dels, cnv = prepare_output(ans_fn, cn, ins, dels)
    zip_fn = join(input_folder, 'snps_{}.zip'.format(f_base))

    with open(ans_fn, 'w') as output_file:
        output_file.write('>' + f_base + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>CNV\n')
        for x in cnv:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        tails = ('>' + x for x in ('STR', 'ALU', 'INV\n'))
        output_file.write('\n'.join(tails))
        output_file.write('>INS\n')
        for x in ins:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in dels:
            output_file.write(','.join([str(u) for u in x]) + '\n')

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(ans_fn)



