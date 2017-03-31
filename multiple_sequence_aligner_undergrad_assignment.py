### multiple_sequence_aligner_undergrad_assignment.py
# Detects SNPs and indels given a reference sequence and reads from a donor sequence.
# First maps reads to reference using a hashed index of the reference. Reads that map poorly are re-aligned using
# the Smith-Waterman algorithm.


import numpy as np
from os.path import join
import time
from collections import defaultdict, OrderedDict
from itertools import groupby
import zipfile

# Reads paired end reads from input file
def read_reads(read_fn):
    '''
    :param read_fn: name of the file containing the paired end reads
    :return: Returns the reads as a list with each pair of paired end reads in a sublist
    '''
    all_reads = []
    with open(read_fn, 'r') as f:
        next(f)
        for line in f:
            line = line.strip()
            paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
            all_reads.append(paired_end_reads)
    return all_reads


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
    :param max_mismatches: maximum amount of mismatches to allow
    :param num_pieces: number of subsections to split read into
    :return:
    '''
    best_ind = None
    curr_sub = None
    indel_pos = None
    least_mismatches = len(read)

    # Divides the read into subsections and looks up each section in the index
    for sub in range(num_pieces):
        subsection = read[sub * length_index:(sub + 1) * length_index]
        if subsection in index:
            for ind in index[subsection]:
                start = ind - sub * length_index
                if start + len(read) <= len(ref):
                    mismatches = [1 if read[j] != ref[start + j] else 0 for j in range(len(read))]
                    n_mismatches = sum(mismatches)
                    if n_mismatches <= max_mismatches:
                        best_ind = start
                        return (best_ind, None, n_mismatches)
                    elif n_mismatches < least_mismatches:
                        best_ind = start
                        curr_sub = sub
                        least_mismatches = n_mismatches
    if best_ind == None: # if read subset not in index
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


# Prints aligned reads into a txt file with each read displayed in 100 bp windows
# Subsequent SNP and indel calling will read from this file
def pretty_print_aligned_reads_with_ref(genome_oriented_reads, read_alignments, ref, read_length=50,
                                        line_length=100, read_sep=100, buffer=50):
    """
    :param genome_oriented_reads: oriented reads generated by an alignment algorithm
    :param read_alignments: alignments generated from an alignment algorithm
    :param ref: reference generated by read_ref
    :return: Returns nothing, but prints the reads aligned to the genome to
     show you what pileup actually *LOOKS* like. You should be able to call SNPs
     by eyeballing the output. However, there are some reads that will not align.
     In the future you'll want to re-check why these reads aren't aligning--the cause
     is usually a structural variation, like an insertion or deletion.
    """
    output_str = ''
    good_alignments = [read_sep + read_length - buffer < x[1] - x[0] <
                       read_sep + read_length + buffer for x in read_alignments]
    # There should be read_length + x (90 < x < 110) p between the reads, and we give a little
    # extra space in case there's been a deletion or insertion.  Depending on the type of
    # deletions/insertions

    best_reads = [genome_oriented_reads[i] for i in range(len(good_alignments))
                  if good_alignments[i]]
    # Remove the reads that do not have a good alignment, or a good reverse alignment.
    best_alignments = [read_alignments[i] for i in range(len(read_alignments))
                       if good_alignments[i]]
    # Take their corresponding alignments
    aligned_reads = [best_reads[i][0] + '.' * (best_alignments[i][1] - best_alignments[i][0] - read_length)
                     + best_reads[i][1] for i in range(len(best_reads))]
    # This turns the reads into strings oriented towards the genome.
    # We get the first read, followed by the correct number of dots to join the first and second reads,
    # and then the second read.

    first_alignment = [x[0] for x in best_alignments]
    alignment_indices = np.argsort(first_alignment)
    sorted_reads = np.array([aligned_reads[i] for i in alignment_indices])
    sorted_alignments = np.array([best_alignments[i] for i in alignment_indices])

    # You don't need to worry too much about how the code block below works--its job is to make it so
    # that a read that starts printing in the third row will continue printing in the third row of the
    # next set of lines.
    active_reads = []
    output_str += '\n\n' + '-' * (line_length + 6) + '\n\n'
    read_indices = np.array([sorted_alignments[j][0]/line_length for j in range(len(sorted_alignments))])

    for i in range(len(ref) / line_length):
        next_ref = ref[i * line_length: (i + 1) * line_length]
        read_mask = (read_indices == i)
        new_alignments = sorted_alignments[read_mask]
        new_reads = sorted_reads[read_mask]
        space_amounts = [_[0] % line_length for _ in new_alignments]
        new_reads_with_spaces = [' ' * space_amounts[j] + new_reads[j] for j in range(len(new_reads))]
        empty_active_read_indices = [index for index in range(len(active_reads)) if active_reads[index] == '']
        for j in range(min(len(new_reads_with_spaces), len(empty_active_read_indices))):
            active_reads[empty_active_read_indices[j]] = new_reads_with_spaces[j]

        if len(new_reads_with_spaces) > len(empty_active_read_indices):
            active_reads += new_reads_with_spaces[len(empty_active_read_indices):]
        printed_reads = ['Read: ' + read[:line_length] for read in active_reads]
        active_reads = [read[line_length:] for read in active_reads]
        while len(active_reads) > 0:
            last_thing = active_reads.pop()
            if last_thing != '':
                active_reads.append(last_thing)
                break
        output_lines = ['Ref:  ' + next_ref] + printed_reads
        output_str += 'Reference index: ' + str(i * line_length) + \
                      '\n' + '\n'.join(output_lines) + '\n\n' + '-' * (line_length + 6) + '\n\n'
    # print output_str
    return output_str


# Aligns reads to a reference sequence and returns potential snps/indels that will be processed further
def alignment_algorithm(paired_end_reads, ref, max_mismatches, num_pieces):
    '''
    :param paired_end_reads: List of paired of reads
    :param ref: Reference sequence
    :param max_mismatches: Maximum amount of mismatches to allow
    :param num_pieces: Number of pieces to divide each read into
    :return all_read_alignment_locations: Each location of the reads
    :return output_read_pairs: Each read pair in the properly oriented direction
    :return all_indels: Each location of potential indels
    '''
    length_read = len(paired_end_reads[0][0])
    length_index = int(length_read/num_pieces)
    index = build_index(ref, length_index)
    all_read_alignment_locations = []
    output_read_pairs = []
    all_indels = []
    count = 0
    start = time.clock()


    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        indels = []
        if count % 1000 == 0:
            time_passed = (time.clock() - start) / 60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed / count * (len(paired_end_reads) - count)
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)
        for read in read_pair:
            ind, indel, mm = compare_reads(ref, read, index, length_index, max_mismatches, num_pieces)
            if mm > max_mismatches:
                rev_read = read[::-1]
                rev_ind, rev_indel, rev_mm = compare_reads(ref, rev_read, index, length_index, max_mismatches, num_pieces)
                if rev_mm < mm and rev_ind is not None:
                    read_alignment_locations.append(rev_ind)
                    output_read_pair.append(rev_read)
                    indels.append(rev_indel)
                elif ind is not None:
                    read_alignment_locations.append(ind)
                    output_read_pair.append(read)
                    indels.append(indel)
            else:
                read_alignment_locations.append(ind)
                output_read_pair.append(read)
                indels.append(indel)
        if len(read_alignment_locations) == 2:
            all_read_alignment_locations.append(read_alignment_locations)
            output_read_pairs.append(output_read_pair)
            all_indels.append(indels)
    return all_read_alignment_locations, output_read_pairs, all_indels


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
        for i in range(1, len(ref)):
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
    edit_matrix, best_moves_matrix = edit_distance_matrix(ref=ref, donor=donor)
    current_row, current_column = np.unravel_index(edit_matrix.argmax(), edit_matrix.shape)
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
    total_indels = defaultdict(list)
    for i in range(len(indels)):
        if indels[i][0] is not None:
            total_indels[indels[i][0]].append(reads[i][0])
        if indels[i][1] is not None:
            total_indels[indels[i][1]].append(reads[i][1])

    total_indels = OrderedDict(sorted(total_indels.items(), key=lambda t: t[0]))
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
    for location, reads in total_grouped_indels.iteritems():
        ref_start = location - 30
        ref_end = location + 30
        if ref_start < 0:
            ref_start = 0
        if ref_end > len(ref) - 1:
            ref_end = len(ref) - 1
        ref_subset = ref[ref_start:ref_end]
        for read in reads:
            start, end, read_start, read_end, changes = identify_changes(ref_subset, read, ref_start)
            if len(changes) > 1:
                if sum(x.count('INS') for x in changes) + sum(x.count('DEL') for x in changes) > 2:
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
    all_changes.sort(key = lambda x: (x[-1], x[0], x[1], x[2]))
    snps = [x for x in all_changes if x[0] == 'SNP']
    snps = [key for key, group in groupby(snps) if len(list(group)) > 5]

    ins = [x for x in all_changes if x[0] == 'INS']
    ins = [key for key, group in groupby(ins) if len(list(group)) > 2]

    dels = [x for x in all_changes if x[0] == 'DEL']
    dels = [key for key, group in groupby(dels) if len(list(group)) > 2]

    return snps, ins, dels

# Generate consensus from file generated from pretty_print_aligned_reads_with_ref function
def generate_consensus(aligned_fn):
    """
    :param aligned_fn: The filename of the saved output of the basic aligner
    :return: SNPs (the called SNPs for uploading to the herokuapp server)
             output_lines (the reference, reads, consensus string, and diff string to be printed)
    """
    with open(aligned_fn, 'r') as input_file:
        line_count = 0
        lines_to_process = []
        SNPs = []
        output_lines = []
        for line in input_file:
            line_count += 1
            line = line.strip()
            if line_count <= 4 or line == '':  # The first 4 lines need to be skipped
                output_lines.append(line)
                continue
            if len(line) > 0 and all(x == '-' for x in line):  # The different pieces of the genome are set off
                                                               # with lines of all dashes '--------'
                new_snps, new_output_lines = process_lines(lines_to_process)
                lines_to_process = []
                SNPs += new_snps
                output_lines += new_output_lines
                output_lines.append(line)
            else:
                lines_to_process.append(line)
        return SNPs, output_lines


def process_lines(genome_lines):
    """

    :param genome_lines: Lines in between dashes from the saved output of the basic_aligner
    :return: snps (the snps from this set of lines)
             output_lines (the lines to print, given this set of lines)
    """
    line_count = 0
    output_lines = []
    consensus_lines = []
    for line in genome_lines:
        output_lines.append(line)
        line_count += 1
        if line_count == 1:  # The first line contains the position in the reference where the reads start.
            raw_index = line.split(':')[1]
            line_index = int(raw_index)
        else:
            consensus_lines.append(line[6:])
    ref = consensus_lines[0]
    reads = consensus_lines[1:]
    consensus_string = consensus(ref, reads)
    diff_string = diff(ref, consensus_string)
    snps = snp_calls(ref, consensus_string, line_index)
    output_lines[2:2] = ['Cons: ' + consensus_string]
    output_lines[2:2] = ['Diff: ' + diff_string]
    return snps, output_lines


def consensus(ref, reads):
    """
    :param ref: reference string
    :param reads: the list of reads.
    :return: The most common base found at each position in the reads (i.e. the consensus string)
    """
    consensus_string = ''
    snp_string = ''
    snp_list = []
    line_length = len(ref)
    padded_reads = [read + ' '*(len(ref) - len(read)) for read in reads]
        # The reads are padded with spaces so they are equal in length to the reference
    for i in range(len(ref)):
        base_count = defaultdict(float)
        ref_base = ref[i]
        base_count[ref_base] += 1.1  # If we only have a single read covering a region, we favor the reference.
        read_bases = [read[i] for read in padded_reads if read[i] not in '. ']
            # Spaces and dots (representing the distance between paired ends) do not count as DNA bases
        for base in read_bases:
            base_count[base] += 1
        consensus_base = max(base_count.iterkeys(), key=(lambda key: base_count[key]))
            # The above line chooses (a) key with maximum value in the read_bases dictionary.
        consensus_string += consensus_base
    return consensus_string


def diff(s1, s2):
    chars = [' ' if s1[i] == s2[i] else '*' for i in range(len(s1))]
    return ''.join(chars)


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


if __name__ == "__main__":

    data_folder = 'hw1_W_2'

    input_folder = '/Users/douglasyao/Downloads/' + data_folder
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    alignments, reads, indels = alignment_algorithm(input_reads, reference, 3, 2)
    total_grouped_indels = collapse_indels(indels)

 #   pickle.dump(total_grouped_indels, open(join(input_folder, 'total_grouped_indels.p'),'wb'))

    snps, ins, dels = check_indels(total_grouped_indels, reference)

    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)


    new_snps, lines = generate_consensus(output_fn)

    snps = [x[1:] for x in snps]

    previous_snp_location = new_snps[0][2]
    itersnp = iter(new_snps)
    next(itersnp)
    for snp in itersnp:
        if snp[2] - previous_snp_location <= 3:
            snp[0] = 'TODELETE'
            previous_snp_location = snp[2]
            continue
        for index in total_grouped_indels.keys():
            if abs(snp[2] - index) < 10:
                snp[0] = 'TODELETE'
                previous_snp_location = snp[2]
                continue
        previous_snp_location = snp[2]
    new_snps = [x for x in new_snps if x[0] != 'TODELETE']
    new_snps.extend(snps)
    new_snps.sort(key = lambda x: x[2])

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




