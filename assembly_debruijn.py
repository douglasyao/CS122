### assembly_debruijn.py
# Performs sequence assembly on reads by generating a De Bruijn graph from the read spectrum

from os.path import join
from collections import defaultdict, Counter
import zipfile
from math import ceil
from multiple_sequence_aligner_undergrad_assignment import read_reads


# Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
def simple_de_bruijn(sequence_reads, k, coverage):
    """
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :param coverage: The required converage for a node to be considered 'real'
    :return: A DeBruijn graph where the keys are k-mers and the values are the set of k-mers that form an edge with the previous k-mer
    """

    de_bruijn_counter = defaultdict(Counter)
    for read in sequence_reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = defaultdict(list)
    for key, val in de_bruijn_counter.iteritems():
        if not val:
            continue
        for end, count in val.iteritems():
            if count > 5:
                if count < ceil(coverage/2.0):
                    de_bruijn_graph[key].append([end, 1])
                else:
                    de_bruijn_graph[key].append([end, int(round(float(count)/coverage))])
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    return de_bruijn_graph


# Given a start point, traverses the DeBruijn Graph created by simple_de_bruijn and returns contigs that come from it.
def de_bruijn_reassemble(current_read, de_bruijn_graph, reads):
    '''
    :param current_read: The read that assembly should be started from
    :param de_bruijn_graph: The generated De Bruijn graph
    :param reads: A list of all reads
    :return all_strings: List of assembled contigs
    :return de_bruijn_graph: Remainder of De Bruijn graph after assembly
    '''
    all_strings = []
    assembled_string = current_read
    current_point = current_read
    while True:
        print len(de_bruijn_graph)
        try:
            next_values = de_bruijn_graph[current_point]
            if not next_values:
                break
            if len(next_values) > 1:
                if max([x[1] for x in next_values]) > 5:
                    del de_bruijn_graph[current_point]
                    break
                correct_value, consensi, to_break = resolve_ambiguities(next_values, assembled_string, reads)
                if to_break:
                    all_strings.extend(consensi)
                    del de_bruijn_graph[current_point]
                    break
                else:
                    next_values = correct_value
            if de_bruijn_graph[current_point][0][1] == 1:
                del de_bruijn_graph[current_point]
            else:
                de_bruijn_graph[current_point][0][1] -= 1
            assembled_string += next_values[0][0][-1]
            current_point = next_values[0][0]
        except KeyError:
            break
    all_strings.append(assembled_string)
    return all_strings, de_bruijn_graph


# If the De Bruijn graph contains a fork, attempts to resolve ambiguity by assembling sequence using a greedy method
def resolve_ambiguities(ambiguity, current_read, reads):
    '''
    :param ambiguity: List containing the ambiguities and their relative frequencies in the De Bruijn graph
    :param current_read: Current assembled sequence
    :param reads: List of all reads
    :return correct_path: Which k-mer is better based on greedy assembly
    :return consensi: All assembled sequences generated from greedy assembly
    :return to_break: If no sequence is able to be assembled, return True
    '''
    if len(current_read) < 50:
        tcurrent_read = next(x for x in reads if current_read in x)
    else:
        tcurrent_read = current_read
    forward_matching_reads = [x for x in reads if tcurrent_read[-25:] in x]
    forward_matching_reads.extend([x[::-1] for x in reads if tcurrent_read[-25:] in x[::-1]])
    count = 2
    while len(forward_matching_reads) < 3:
        tcount = 0
        for x in reads:
            if current_read in x:
                tcount += 1
            if tcount == count:
                count += 1
                tcurrent_read = x
                break
        forward_matching_reads = [x for x in reads if tcurrent_read[-25:] in x]
        forward_matching_reads.extend([x[::-1] for x in reads if tcurrent_read[-25:] in x[::-1]])

    matching_locations = map_reads_to_ref(tcurrent_read, forward_matching_reads, 25)
    consensi, to_break = get_consensus(matching_locations, reads)
    consensi = [x for x in consensi if x]
    if not to_break:
        correct_path = [x for x in ambiguity if x[0][-10:] in consensi[0]]
    else:
        correct_path = []
    return correct_path, consensi, to_break


# Maps reads to a reference
def map_reads_to_ref(ref, forward_reads, end):
    '''
    :param ref: The reference sequence
    :param forward_reads: List of reads
    :param end: How far to align
    :return: Alignment locations
    '''
    all_locations = [[ref,0]]
    offset = len(ref) - 50
    if forward_reads:
        for read in forward_reads:
            index = 0
            while (index <= len(ref) - end):
                mismatches = sum([1 if ref[offset + j + index] != read[j] else 0 for j in range(len(read) - index)])
                if mismatches > 2:
                    index += 1
                else:
                    all_locations.append([read, offset+index])
                    break
    all_locations.remove([ref,0])
    return all_locations


# Performs assembly using a greedy algorithm
def get_consensus(matching_locations, reads):
    '''
    :param matching_locations: Read locations relative to partially assembled reference
    :param reads: List of reads to assemble
    :return consensi: List of all assembled sequences
    :return to_break: If no sequence is able to be assembled, return True
    '''
    to_break = False
    consensi = []
    consensus = ''
    bases = defaultdict(Counter)
    for read in matching_locations:
        for i in range(len(read[0])):
            bases[read[1] + i].update([read[0][i]])
    for j in range(0, len(bases)):
        max = bases[j].most_common(1)[0][1]
        to_add = [x for x in bases[j].most_common(2) if x[1] > 10 and x[1] > max/4]
        if len(to_add) > 1:
            for add in to_add:
                temp_string = consensus + add[0]
                forward_matching_reads = [x for x in reads if temp_string[-15:] in x]
                forward_matching_reads.extend([x[::-1] for x in reads if temp_string[-15:] in x[::-1]])
                temp_matching_locations = map_reads_to_ref(temp_string, forward_matching_reads, 25)
                if not temp_matching_locations:
                    to_add = True
                    to_break = True
                    continue
                temp_consensi, temp_to_break = get_consensus(temp_matching_locations, reads)
                consensi.extend(temp_consensi)
                to_break = True
            break
        elif to_add:
            consensus = consensus + to_add[0][0]
    if j == len(bases) - 1:
        consensi.append(consensus)
    return consensi, to_break


# Calls de_bruijn_reassemble repeatedly until entire graph is traversed to generate all contigs
def de_bruijn(de_bruijn_graph, reads):
    '''
    :param de_bruijn_graph: the generated De Bruijn graph
    :param reads: List of all reads
    :return: List of all contigs
    '''
    assembled_strings = []
    ambiguities = {}
    good_starts = []
    all_ins = [x for x in de_bruijn_graph.values()]
    all_ins = [item for sublist in all_ins for item in sublist]

    ins_hash = defaultdict(int)
    for x in all_ins:
        ins_hash[x[0]] += x[1]

    for key, val in de_bruijn_graph.iteritems():
        outs = sum([x[1] for x in val])
        if len(val) > 1:
            ambiguities[key] = val
        if outs == ins_hash[key] + 1:
            good_starts.append(key)

    ind = 0
    current_read = good_starts[0]
    while ind < len(good_starts) - 1:
        strings, de_bruijn_graph = de_bruijn_reassemble(current_read, de_bruijn_graph, reads)
        assembled_strings.extend(strings)
        ind += 1
        current_read = good_starts[ind]
    while len(de_bruijn_graph) > 0:
        current_read = de_bruijn_graph.keys()[0]
        strings, de_bruijn_graph = de_bruijn_reassemble(current_read, de_bruijn_graph, reads)
        assembled_strings.extend(strings)
    return assembled_strings


# For testing cases. Prints out aligned reads with appropriate amount of whitespace to make them line up
def print_aligned_reads_with_ref(read_locations, filename):
    '''
    :param read_locations: List containing reads and corresponding location in sublist
    :param filename: Name of output file to print to
    '''
    all = []
    read_locations.sort(key = lambda x: x[1])
    start = min(x[1] for x in read_locations)
    for x in range(len(read_locations)):
        read_locations[x][1] -= start
    for read, location in read_locations:
        for i in range(location):
            read = ' ' + read
        all.append(read)
    with open(filename, 'w') as f:
        f.write('\n'.join(all))
    return


if __name__ == "__main__":
    chr_name = 'hw3all_A_3_chr_1'
    input_folder = '/Users/douglasyao/Downloads/hw3all_A_3/'
    reads_fn = join(input_folder, 'reads_{}.txt'.format(chr_name))
    reads = read_reads(reads_fn)
    reads = [item for sublist in reads for item in sublist]

    db_graph = simple_de_bruijn(reads, 25, 20)

    output = de_bruijn(db_graph, reads)
    output.sort(key = lambda x: len(x), reverse=True)


    output_fn_end = 'assemble_temp_{}.txt'.format(chr_name)
    output_fn = join(input_folder, output_fn_end)
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))
    zip_fn = join(input_folder, 'assemble_temp_{}.zip'.format(chr_name))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)