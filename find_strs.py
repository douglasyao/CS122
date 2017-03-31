### find_strs.py
# Detects short tandem repeats (STRs) in a reference sequence

from os.path import join
from itertools import groupby, product
import zipfile
import re
from multiple_sequence_aligner_undergrad_assignment import read_reference
from math import ceil


# Reads in answer file containing donor SNPs, indels, CNVs, and STRs so that they can be appended to
def read_ans(ans_fn):
    '''
    :param ans_fn: Name of answer file
    :return: All SNPs, insertions, deletions, SNVs, and STRs
    '''
    insertions = []
    deletions = []
    snps = []
    cnv = []
    str = []
    with open(ans_fn) as f:
        next(f)
        state = None
        for line in f:
            line = line.strip()
            if line[0] == '>':
                state = line[1:]
            elif state == 'SNP':
                snps.append(line.split(','))
            elif state == 'INS':
                insertions.append(line.split(','))
            elif state == 'DEL':
                deletions.append(line.split(','))
            elif state == 'CNV':
                cnv.append(line.split(','))
            elif state == 'STR':
                str.append(line.split(','))
    return snps, insertions, deletions, cnv, str


# Generates all possible SNP subunits up to a certain length
def get_combinations(length):
    '''
    :param length: Maximum length of SNP subunit
    :return combos: All SNP subunits
    :return final_combos: SNP subunits replicated to form the smallest actual SNP fragment. Length of fragments are
                          smallest multiple of SNP subunit length that is greater than 10.
    '''
    combos = []
    final_combos = []
    for i in range(2, length + 1):
        for combo in product('ACGT',repeat=i):
            to_continue = False
            if len(set(combo)) == 1:
                continue
            combo = ''.join(combo)
            if combo in combos:
                continue
            for j in range(1, i):
                combo = combo[-1] + combo[:i-1]
                if combo in combos:
                    to_continue = True
                    break
            if to_continue:
                continue
            else:
                combos.append(combo)
    for combo in combos:
        num_reps = int(ceil(float(10)/len(combo)))
        to_add = ''.join([combo for i in range(num_reps)])
        final_combos.append(to_add)
    return combos, final_combos


# If multiple STRs detected adjacent to each other, combines them into a single STR
def collapse_list(list, diff):
    '''
    :param list: List of STR locations
    :param diff: Length of STR fragment
    :return:
    '''
    dups = []
    start = list[0]
    prev = start
    for i in list[1:]:
        if i == list[-1]:
            if i - prev == diff:
                dups.append((start, i+diff))
            else:
                dups.append((i, i+diff))
        elif i - prev == diff:
            prev = i
        else:
            dups.append((start, prev+diff))
            start = i
            prev = i
    return dups


# Finds STRs in reference sequence
def identify_strs_in_reference(reference):
    '''
    :param reference: The reference sequence
    :return: List of STRs and locations
    '''
    subset, combos = get_combinations(5)
    strs = []
    for i, combo in enumerate(combos):
        print i, combo
        matches = [x.start() for x in re.finditer(combo, reference)]
        if not matches:
            continue
        matches = collapse_list(matches, len(combo))
        for start, end in matches:
            to_add_back = 0
            to_add_front = 0
            smallest = subset[i]
            match = reference[start:end]
            while reference[start:end+to_add_back+1] == match + smallest[to_add_back%len(smallest)]:
                match += subset[i][to_add_back%len(subset[i])]
                to_add_back += 1
            while reference[start-to_add_front-1:end+to_add_back] == smallest[len(smallest)-to_add_front%len(smallest)-1] + match:
                match = smallest[len(smallest)-to_add_front%len(smallest)-1] + match
                to_add_front += 1
            if len(match) >= 12:
                strs.append([match, start - to_add_front])
    strs.sort()
    strs = [strs for strs, x in groupby(strs)]
    return strs


# Copy number analyses will often incorrectly detect STRs as copy number changes. This function will remove STRs in the CNV section of the answer file
def remove_strs_from_cnv(cnv, strs):
    '''
    :param cnv: List of CNVs from the answer file
    :param strs: List of STRs from identify_strs_in_reference
    :return:
    '''
    str_seqs = [x[0] for x in strs]
    new_cnv = []
    for i in cnv:
        if not in_str(i[0], str_seqs):
            new_cnv.append(i)
    return new_cnv


# Checks whether a sequence is a subset of any STR in a list of STRs
def in_str(seq, strs):
    '''
    :param seq: The sequence
    :param strs: List of STRs
    :return: True if the sequence is in any of the STRs
    '''
    for str in strs:
        if seq in str:
            return True
    return False


if __name__ == '__main__':
    data_folder = 'hw2grad_M_1'
    # data_folder = 'practice_W_1'

    input_folder = 'C:/Users/DWYao/Downloads/hw2grad_M_1'
    f_base = '{}_chr_1'.format(data_folder)
    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    ans_fn = join(input_folder, 'snps_{}.txt'.format(f_base))
    snps, ins, dels, cnv, strs = read_ans(ans_fn)
    strs = identify_strs_in_reference(reference)
    # cnv = remove_strs_from_cnv(cnv, strs)
    strs = sorted(strs, key = lambda x: x[1])
    zip_fn = join(input_folder, 'snps_{}.zip'.format(f_base))

    with open(ans_fn, 'w') as output_file:
        output_file.write('>' + f_base + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>CNV\n')
        for x in cnv:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>STR\n')
        for x in strs:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        tails = ('>' + x for x in ('ALU', 'INV\n'))
        output_file.write('\n'.join(tails))
        output_file.write('>INS\n')
        for x in ins:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in dels:
            output_file.write(','.join([str(u) for u in x]) + '\n')

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(ans_fn)

