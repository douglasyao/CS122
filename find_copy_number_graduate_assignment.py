### find_copy_number_graduate_assignment.py
# Finds copy number of reference sequence only.
#
# Designed to handle much larger datasets than find_copy_number.py


import numpy as np
from os.path import join
from collections import defaultdict, OrderedDict
from itertools import groupby, combinations
import zipfile
import re
from multiple_sequence_aligner_undergrad_assignment import build_index
import cPickle as pickle
from find_strs import read_ans, remove_strs_from_cnv
from multiple_sequence_aligner_graduate_assignment import file_len


# Given a dict of coverage of each section of reference aligned to itself, return start and ends of sections where coverage exceeds 1x
def collapse_dict(all_starts, diff):
    '''
    :param all_starts: Dict of reference coverage
    :return: start and end locations of sections
    '''
    dups = []
    iter_starts = all_starts.iteritems()
    prev, b = next(iter_starts)
    start = prev
    active = False
    for k, v in iter_starts:
        if v != 1:
            if not active:
                start = k
                prev = start
                active = True
            elif k - prev == diff:
                prev = k
        elif active:
            if k - start > diff:
                dups.append(int(np.mean([start,prev])))
            active = False
    return dups

# Find copy number in reference. Reference is divided into several sections and copy number is found independently for each section.
def find_copy_number_reference(ref_fn, size, divisions, fragment_size):
    '''
    :param ref_fn: File name of file containing reference sequence
    :param size: Length of reference sequence
    :param divisions: Number of sections to divide reference into
    :param fragment_size: Size of each 'bin' to calculate coverage
    :return: Reference copy number sequences and locations
    '''
    lines_ref = file_len(ref_fn) - 1
    all_seqs = defaultdict(list)

    # length_ref / divisions * i
    # Divides the reference into i pieces and performs alignment of all reads against each piece separately, then combining the results.
    for num_fragment in range(divisions):
        seqs = []
        print 'Iteration {} out of {}'.format(num_fragment + 1, divisions)
        reference = ''
        with open(ref_fn, 'r') as f:  # Read in the ith piece of the reference
            next(f)
            for j in range(0, lines_ref / divisions * num_fragment):
                next(f)
            count = 0
            for line in f:
                line = line.strip()
                reference += line  # We append each line to the output reference string.
                count += 1
                if count >= lines_ref / divisions:
                    break
        index = build_index(reference, fragment_size)  # Build index from reference subset


        all_starts = OrderedDict([(fragment_size*i, 0) for i in range((size/divisions)/fragment_size)])
        seqs = []
        for i in range((size/divisions)/fragment_size):
            subsection = reference[i * fragment_size:(i + 1) * fragment_size]
            if subsection in index:
                start = round_down(index[subsection][0], fragment_size)
                all_starts[start] += 1
        dups = collapse_dict(all_starts, fragment_size)
        for j in dups:
            print j
            sec = reference[j-fragment_size/2:j+fragment_size/2]
            forward = 0
            backward = 0
            matches = re.findall(sec, reference)
            if len(matches) > 1:
                forward += 1
                sec = reference[j:j+fragment_size+forward]
                new_match = re.findall(sec, reference)
                while len(new_match) == len(matches):
                    forward += 1
                    sec = reference[j:j+fragment_size+forward]
                    new_match = re.findall(sec, reference)
                backward += 1
                sec = reference[j-backward:j+fragment_size+forward-1]
                new_match = re.findall(sec, reference)
                while len(new_match) == len(matches):
                    backward += 1
                    sec = reference[j-backward:j+fragment_size+forward-1]
                    new_match = re.findall(sec, reference)
                sec = sec[1:]
            indices = [x.start() for x in re.finditer(sec, reference)]
            seqs.append([sec, indices])
        seqs = [[sec, indices] for sec, indices in seqs if len(sec) > 15 and len(indices) > 1]
        seqs.sort()
        seqs = [seqs for seqs, x in groupby(seqs)]
        for seq, indices in seqs:
            new = [x+size/divisions*num_fragment for x in indices]
            all_seqs[seq].extend(new)
    return all_seqs


# Round a number down to the nearest multiple of a given number
def round_down(num, divisor):
    return num - (num%divisor)


if __name__ == "__main__":

    data_folder = 'hw2grad_M_1'
    # data_folder = 'practice_E_1'

    input_folder = 'C:/Users/DWYao/Downloads/' + data_folder
    f_base = '{}_chr_1'.format(data_folder)
    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    ans_fn = join(input_folder, 'snps_{}.txt'.format(f_base))


    cn = find_copy_number_reference(reference_fn, 100000000, 20, 25)
    pickle.dump(cn, open("ref_cn.p", "wb"))
    exit(0)

    cn = pickle.load(open("ref_cn.p", "rb"))

    for combo in combinations([[k,v] for k,v in cn.iteritems()], 2):
        combo_seqs = [x[0] for x in combo]
        if combo_seqs[0] in combo_seqs[1] or combo_seqs[1] in combo_seqs[0]:
            to_compare = [[k,v] for k,v in cn.iteritems() if k in combo_seqs]
            if len(to_compare) > 1:
                if abs(len(combo_seqs[0]) - len(combo_seqs[1])) < 5:
                    to_remove = max(to_compare, key=lambda x: len(x[0]))
                    to_keep = min(to_compare, key=lambda x: len(x[0]))
                    cn[to_keep[0]].extend(to_remove[1])
                    cn.pop(to_remove[0])

    cn = [[k, sorted(v)] for k, v in cn.iteritems()]
    cn = sorted(cn, key = lambda x: x[1][0])

    snps, ins, dels, cnv, strs = read_ans(ans_fn)
    cn = remove_strs_from_cnv(cn, strs)

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



