### multiple_sequence_aligner_graduate_assignment.py
# Detects SNPs and indels given a reference sequence and reads from a donor sequence.
# First maps reads to reference using a hashed index of the reference. Reads that map poorly are re-aligned using
# the Smith-Waterman algorithm.
#
# Same as multiple_sequence_aligner_undergrad_assignment.py except designed to handle much larger datasets
# (ref sequence ~100 million base pairs, ~60 million reads from donor sequence)




from os.path import join
import time
from collections import defaultdict, OrderedDict
import zipfile
from itertools import groupby
import numpy as np

# Reads reference from input file
def read_reference(ref_fn):
    '''
    :param ref_fn: name of the file containing the reference sequence
    :return: Returns the reference sequence
    '''
    with open(ref_fn, 'r') as f:
        next(f)
        output_reference = ''
        for line in f:
            line = line.strip()
            output_reference += line  # We append each line to the output reference string.
    return output_reference


# Build an index from a reference sequence
def build_index(ref, length):
    '''
    :param ref: the reference sequence
    :param length: length of each index
    :return: the index as a dict
    '''
    index = defaultdict(list)
    for i in range(0, len(ref)-length):
        subset = ref[i:i+length]
        index[subset].append(i)
    return index


# Aligns a read to a reference sequence based on an index. Returns the location and total number of mismatches for the best alignment.
# If a potential indel is detected, also returns the location of the indel.
def compare_reads(ref, read, index, length_index, max_mismatches, num_pieces):
    '''
    :param ref: the reference sequence
    :param read: the read
    :param index: the indexed reference sequence as a dict
    :param length_index: the length of each key in the index
    :param max_mismatches: maximum number of mismatches to allow
    :param num_pieces: number of subsections to split the read into
    :return ind: the location on the reference where the read matches with the least amount of mismatches
    :return indel_pos: location of a possible indel
    :return least_mismatches: the number of mismatches between the read and the reference at the best location
    '''
    best_ind = None
    curr_sub = None
    indel_pos = None
    least_mismatches = len(read)

    # Divides the read into subsections and looks up each section in the index
    for sub in range(num_pieces):
        subsection = read[sub * length_index:(sub + 1) * length_index] # Subset read
        if subsection in index:
            for ind in index[subsection]:
                start = ind - sub * length_index # Find start of read
                if start < 0: # If read starts before the beginning of the reference, then truncate read
                    read = read[abs(start):]
                    mismatches = [1 if read[j] != ref[j] else 0 for j in range(len(read))]
                    n_mismatches = sum(mismatches)
                elif start + len(read) >= len(ref): # If read ends after the end of the reference, then truncate read
                    read = read[:(len(ref) - start)]
                    mismatches = [1 if read[j] != ref[start + j] else 0 for j in range(len(read))]
                    n_mismatches = sum(mismatches)
                else:
                    mismatches = [1 if read[j] != ref[start + j] else 0 for j in range(len(read))]
                    n_mismatches = sum(mismatches)
                if n_mismatches <= max_mismatches: # If few mismatches, then return location. No indels.
                    best_ind = start
                    return (best_ind, None, n_mismatches)
                elif n_mismatches < least_mismatches:
                    best_ind = start
                    curr_sub = sub
                    least_mismatches = n_mismatches
    if best_ind == None: # If read subset not in index
        return (None, None, least_mismatches)

    # The following code looks for potential indels. The idea is that a read that covers an indel will match (almost) perfectly
    # on the reference up until the location of the indel, after which it will match very poorly. We can therefore start from
    # the end of the subsection that matches the index and work our way forwards/backward until we hit a mismatch, then record
    # this position as a potential indel.

    # Assuming the read is divided into two subsets
    if curr_sub == 0: # If the first subset matches the index
        read_indel_pos = length_index
        indel_pos = best_ind + length_index
        while read_indel_pos < len(read) and read[read_indel_pos] == ref[indel_pos]: # Work our way forward until we hit a mismatch
            read_indel_pos += 1
            indel_pos += 1

    elif curr_sub == 1: # If the second subset matches the index
        read_indel_pos = length_index
        indel_pos = best_ind + length_index
        while read_indel_pos >= 0 and read[read_indel_pos] == ref[indel_pos]: # Work our way backward until we hit a mismatch
            read_indel_pos -= 1
            indel_pos -= 1

    return (best_ind, indel_pos, least_mismatches)


# Compare consensus to reference and identify SNPs
def snp_calls(ref_string, consensus_string, start_index):
    """
    :param ref_string: A piece of the reference string
    :param consensus_string: A piece of the consensus string
    :param start_index: The start
    :return: Correctly formatted SNPs for output to the herokuapp server.
    """
    snps = []
    for i in range(len(ref_string)):
        if ref_string[i] != consensus_string[i]:
            snps.append([ref_string[i], consensus_string[i], start_index + i])
    return snps


# Get the number of lines in a file
def file_len(fname):
    '''
    :param fname: name of the file
    :return: number of lines in the file
    '''
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# Aligns reads to a reference sequence and returns potential snps/indels that will be processed further. Also has the option to split
# the reference into multiple sections and perform alignment of each section individually in case the overall reference is too long to index.
def alignment_algorithm(read_fn, ref_fn, length_ref, divisions, max_mismatches, num_pieces):
    '''
    :param read_fn: Name of the file containing the reads
    :param ref_fn: Name of the file containing the reference
    :param length_ref: Length of the reference
    :param divisions: Number of times we want to split up the reference. Each piece can be around 5 million bp long.
    :param max_mismatches: Maximum number of mismatches per read
    :param num_pieces: Number of times we want to split up each read.
    :return all_snps: Location and bps of each snp. Will contain errors because indels not accounted for.
    :return all_indels: Location of each potential indel as well as the reads that map to that location as a dict.
    '''
    lines_ref = file_len(ref_fn) - 1
    all_snps = []
    all_indels = defaultdict(list)

    # Divides the reference into i pieces and performs alignment of all reads against each piece separately, then combining the results.
    for i in range(divisions):
        print 'Iteration {} out of {}'.format(i+1, divisions)
        output_reference = ''
        consensus = ''
        with open(ref_fn, 'r') as f: # Read in the ith piece of the reference
            next(f)
            for j in range(0, lines_ref/divisions*i):
                next(f)
            count = 0
            for line in f:
                line = line.strip()
                output_reference += line  # We append each line to the output reference string.
                count += 1
                if count >= lines_ref/divisions:
                    break
        length_read = 50
        length_index = 25
        index = build_index(output_reference, length_index) # Build index from reference subset
        count = 0
        start = time.clock()
        each_location = OrderedDict()

        # Reads are read in line by line. In order to generate the consensus sequence, the reads are first aligned, then each base
        # pair is added to a dict containing the number of each base at each location of the reference. Afterward, we can look up
        # the most common base at each position
        with open(read_fn, 'r') as f:
            next(f)
            for k in range(length_ref/divisions):
                bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                each_location[k] = bases
            for line in f:
                line = line.strip()
                paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
                read_alignment_locations = []
                output_read_pair = []
                indels = []

                for read in paired_end_reads:
                    count += 1
                    if count % 1000000 == 0:
                        time_passed = (time.clock() - start) / 60
                        print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)

                    ind, indel, mm = compare_reads(output_reference, read, index, length_index, max_mismatches, num_pieces) # Compare paired end read to reference
                    if mm > max_mismatches:
                        rev_read = read[::-1] # If forward read doesn't match well, try reversing the read and compare again
                        rev_ind, rev_indel, rev_mm = compare_reads(output_reference, rev_read, index, length_index, max_mismatches, num_pieces)
                        if rev_mm < mm and rev_ind is not None: # if reverse matches better than forward
                            read_alignment_locations.append(rev_ind)
                            output_read_pair.append(rev_read)
                            indels.append(rev_indel)
                        elif ind is not None: # if forward read still matches better
                            read_alignment_locations.append(ind)
                            output_read_pair.append(read)
                            indels.append(indel)
                    else: # if forward read matches well
                        read_alignment_locations.append(ind)
                        output_read_pair.append(read)
                        indels.append(indel)
                if len(read_alignment_locations) == 2:
                    for l in range(2):
                        if read_alignment_locations[l] < 0: # if read starts before the beginning of the reference, truncate the read
                            output_read_pair[l] = output_read_pair[l][abs(read_alignment_locations[l]):]
                            for m in range(50 + read_alignment_locations[l]): # store bp counts at each location
                                to_add = output_read_pair[l][m]
                                each_location[m][to_add] += 1
                        else:
                            for m in range(50): # store bp counts at each location
                                if read_alignment_locations[l] + m >= len(output_reference): # if read ends after the end of the reference, truncate the read
                                    break
                                to_add = output_read_pair[l][m]
                                each_location[read_alignment_locations[l] + m][to_add] += 1
                if len(indels) > 0:
                    for n in range(len(indels)):
                        if indels[n] is not None:
                            all_indels[(indels[n] + length_ref/divisions*i)].append(output_read_pair[n]) # store potential indel location and the read that maps there

        # Generate consensus sequence
        for ind, bases in each_location.iteritems():
            cons = max(bases, key=bases.get)
            dups = [k for k,v in bases.iteritems() if v == max(bases.values())]
            if cons in dups and len(dups) > 1: # If there is a tie between bases, favor the reference
                cons = output_reference[ind]
            consensus += cons
        snps = snp_calls(output_reference, consensus, length_ref/divisions*i) # Generate snps from consensus
        all_snps.extend(snps)
    return all_snps, all_indels


# Performs local alignment between a donor and reference sequence according to the Smith-Waterman algorithm.
def edit_distance_matrix(ref, donor):
    '''
    :param ref: The reference sequence
    :param donor: The donor sequence
    :return output_matrix: The distance matrix containing the scores at each location
    :return best_moves_matrix: A matrix containing the best move at each location. Each entry in the matrix is a set of coordinates
    that indicates whether to go down (for insertion), left (for deletion), or diagonally (for identity/substitution)
    '''
    output_matrix = np.zeros((len(ref), len(donor)))
    best_moves_matrix = np.zeros((len(ref), len(donor)), dtype = (int, 2))

    for j in range(1, len(donor)):
        for i in range(1, len(ref)): # scores can be tweaked to give different results
            deletion = output_matrix[i - 1, j] - 3
            insertion = output_matrix[i, j - 1] - 3
            identity = output_matrix[i - 1, j - 1] + 5 if ref[i] == donor[j] else -np.inf
            substitution = output_matrix[i - 1, j - 1] - 2 if ref[i] != donor[j] else -np.inf
            start = 0
            best_score = max(insertion, deletion, identity, substitution, start)
            output_matrix[i, j] = best_score
            if best_score == insertion:
                best_moves_matrix[i, j] = (0, -1)
            elif best_score == deletion:
                best_moves_matrix[i, j] = (-1, 0)
            elif best_score == identity or best_score == substitution:
                best_moves_matrix[i, j] = (-1, -1)
    return output_matrix, best_moves_matrix


# Compares a donor to a reference sequence according the Smith-Waterman algorithm. Performs a backtrace-based re-alignment
# of the donor to the reference and identifies SNPs, insertions, and deletions.
def identify_changes(ref, donor, offset):
    '''
    :param ref: The reference sequence
    :param donor: The donor sequence
    :param offset: Where the reference starts, if the inputted reference sequence is a subset of the entire reference
    :return ref_start: Where alignment starts on the reference
    :return ref_end: Where alignment ends on the reference
    :return read_start: Where alignment starts on the donor
    :return read_end: Where alignment ends on the donor
    :return changes: A list of SNPs, insertions, and deletions, and their locations
    '''
    ref = '${}'.format(ref)
    donor = '${}'.format(donor)
    edit_matrix, best_moves_matrix = edit_distance_matrix(ref=ref, donor=donor) # Create distance matrix
    current_row, current_column = np.unravel_index(edit_matrix.argmax(), edit_matrix.shape) # Find starting location which is the highest score
    ref_end = offset + current_row + 1
    read_end = current_column + 1
    current_move = best_moves_matrix[current_row, current_column]
    changes = []
    while not (current_move[0] == 0 and current_move[1] == 0):
        ref_index = current_row - 1
        if current_move[0] == -1 and current_move[1] == 0:
            if len(changes) > 0 and changes[-1][0] == 'DEL' and changes[-1][-1] == offset + ref_index + 1:
                changes[-1] = ['DEL', ref[current_row] + changes[-1][1], offset + ref_index]
            else:
                changes.append(['DEL', ref[current_row], offset + ref_index])
        elif current_move[0] == 0 and current_move[1] == -1:
            if len(changes) > 0 and changes[-1][0] == 'INS' and changes[-1][-1] == offset + ref_index + 1:
                changes[-1][1] = donor[current_column] + changes[-1][1]
            else:
                changes.append(['INS', donor[current_column], offset + ref_index + 1])
        elif current_move[0] == -1 and current_move[1] == -1 and ref[current_row] != donor[current_column]:
            changes.append(['SNP', ref[current_row], donor[current_column], offset + ref_index])
        current_row += current_move[0]
        current_column += current_move[1]
        current_move = best_moves_matrix[current_row, current_column]
    changes = sorted(changes, key=lambda change: change[-1])
    ref_start = offset + current_row
    read_start = current_column
    return ref_start, ref_end, read_start, read_end, changes


# Given a dictionary containing all possible indel locations and reads that map there, collapse the entries if they are represented
# by at least 2 reads or if they fall within 10 bp of the previous indel location. In the latter case, the locations are collapsed to
# the average of all locations.
def collapse_indels(indels):
    '''
    :param indels: A dictionary containing all possible indel locations and reads that map there
    :return: A collapsed dictionary of indel locations and reads
    '''
    total_indels = OrderedDict(sorted(indels.items(), key=lambda t: t[0]))
    current_indel = total_indels.keys()[0]
    grouped_indel_locations = []
    grouped_indels = []
    total_grouped_indels = {}
    for k, v in total_indels.iteritems():
        if k - current_indel <= 10:
            grouped_indel_locations.append(k)
            grouped_indels.extend(v)
        elif len(grouped_indels) >= 2:
            ave_location = int(np.mean(grouped_indel_locations))
            total_grouped_indels[ave_location] = grouped_indels
            grouped_indel_locations = [k]
            grouped_indels = v
        else:
            grouped_indel_locations = [k]
            grouped_indels = v
        current_indel = k
    return total_grouped_indels


# Given a dictionary containing collapsed indel locations and read that map there, compare each read against the reference sequence
# near the location of the indel using Smith-Waterman. If a SNP/insertion/deletion is found in most reads, then we can be confident
# that it is real.
def check_indels(total_grouped_indels, ref):
    '''
    :param total_grouped_indels: A collapsed dictionary of indels and reads
    :param ref: The reference sequence
    :return: SNPs, insertions, and deletions, and their locations
    '''
    all_changes = []
    count = 0
    new_start = time.clock()
    number_indels = len(total_grouped_indels)
    for location, reads in total_grouped_indels.iteritems():
        count += 1
        if count % 100 == 0:
            time_passed = (time.clock() - new_start) / 60
            print '{} out of {} locations re-aligned in {:.3} minutes'.format(count, number_indels, time_passed)
            remaining_time = time_passed / count * (number_indels - count)
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)
        ref_start = location - 30 # perform Smith-Waterman re-alignment on the section surrounding the potential indel
        ref_end = location + 30
        if ref_start < 0:
            ref_start = 0
        if ref_end > len(ref) - 1:
            ref_end = len(ref) - 1
        ref_subset = ref[ref_start:ref_end]
        for read in reads:
            start, end, read_start, read_end, changes = identify_changes(ref_subset, read, ref_start)

            # One issue of Smith-Waterman is that is assumes the independence of insertion/deletion locations, which is not true.
            # Indels occur in contiguous pieces. It is much more likely to have an indel of length of 5 than two indels of length 2 and 3
            # separated by one bp. The following detects indels that are within 3 bp of each other and connects then together.
            if len(changes) > 1:
                if sum(x.count('INS') for x in changes) + sum(x.count('DEL') for x in changes) > 2: # Ignore reads w/ too many indels/deletions
                    continue
                elif sum(x.count('INS') for x in changes) == 2:
                    indels = [x for x in changes if x[0] == 'INS']
                    changes = [x for x in changes if x[0] != 'INS']
                    if indels[1][-1] - indels[0][-1] <= 3:
                        combine_len = len(indels[0][1]) + len(indels[1][1])
                        indels[0][1] = read[indels[0][-1] - start + read_start:indels[0][-1] - start + read_start + combine_len]
                        del indels[1]
                        all_changes.extend(indels)
                elif sum(x.count('DEL') for x in changes) == 2:
                    indels = [x for x in changes if x[0] == 'DEL']
                    changes = [x for x in changes if x[0] != 'DEL']
                    if indels[1][-1] - indels[0][-1] <= 3:
                        combine_len = len(indels[0][1]) + len(indels[1][1])
                        indels[0][1] = ref[indels[0][-1]:indels[0][-1] + combine_len]
                        del indels[1]
                        all_changes.extend(indels)
            all_changes.extend(changes)

    # Changes that occur in many reads are more likely to be true, so only keep these
    all_changes.sort(key = lambda x: (x[-1], x[0], x[1], x[2]))
    snps = [x for x in all_changes if x[0] == 'SNP']
    snps = [key for key, group in groupby(snps) if len(list(group)) > 5]

    ins = [x for x in all_changes if x[0] == 'INS']
    ins = [key for key, group in groupby(ins) if len(list(group)) > 2]

    dels = [x for x in all_changes if x[0] == 'DEL']
    dels = [key for key, group in groupby(dels) if len(list(group)) > 2]

    return snps, ins, dels

# The SNPs generated from the alignment_algorithm function will contain errors, since indels will be incorrectly reported as SNPs.
# This function removes SNPs that are either within 3 bp of each other or fall within a certain distance to the potential indel locations.
def clean_snps(new_snps, grouped_indels, distance_to_indel):
    '''
    :param new_snps: List of SNPs
    :param grouped_indels: List of collapsed potential indels locations
    :param distance_to_indel: SNPs within this distance to the indels will be removed
    :return:
    '''
    previous_snp_location = new_snps[0][2]
    itersnp = iter(new_snps)
    next(itersnp)
    for snp in itersnp:
        if snp[2] - previous_snp_location <= 3:
            snp[0] = 'TODELETE'
            previous_snp_location = snp[2]
            continue
        for index in grouped_indels.keys():
            if abs(snp[2] - index) < distance_to_indel:
                snp[0] = 'TODELETE'
                previous_snp_location = snp[2]
                continue
        previous_snp_location = snp[2]

    new_snps = [x for x in new_snps if x[0] != 'TODELETE']
    return new_snps

if __name__ == "__main__":
  #  data_folder = 'hw2grad_M_1_copy2'
  #  input_folder = join('/u/scratch/d/douglasy/', data_folder)

    data_folder = 'hw2undergrad_E_2'
    input_folder = join('/Users/douglasyao/Downloads/', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    new_snps, indels = alignment_algorithm(reads_fn, reference_fn, 1000000, 5, 3, 2)
  #  new_snps, indels = alignment_algorithm(reads_fn, reference_fn, 100000000, 20, 3, 2)

    grouped_indels = collapse_indels(indels)

    reference = read_reference(reference_fn)
    snps, ins, dels = check_indels(grouped_indels, reference)

    snps = [x[1:] for x in snps]

    snps = clean_snps(snps, grouped_indels, 20)
    new_snps = clean_snps(new_snps, grouped_indels, 10)

    new_snps.extend(snps)
    new_snps.sort(key=lambda x: x[2])

    final_output_fn = join(input_folder, 'snps_{}.txt'.format(f_base))
    zip_fn = join(input_folder, 'snps_{}.zip'.format(f_base))

    with open(final_output_fn, 'w') as output_file:
        output_file.write('>' + f_base + '\n>SNP\n')
        for x in new_snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV\n'))
        output_file.write('\n'.join(tails))
        output_file.write('>INS\n')
        for x in ins:
            output_file.write(','.join([str(u) for u in x[1:]]) + '\n')
        output_file.write('>DEL\n')
        for x in dels:
            output_file.write(','.join([str(u) for u in x[1:]]) + '\n')

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(final_output_fn)
