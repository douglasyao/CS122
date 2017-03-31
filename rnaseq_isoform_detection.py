### rnaseq_isoform_detection.py
# Given a set of RNA isoforms and RNA reads, calculates the relative frequencies of each isoform using least squares

import numpy as np
import pandas as pd
from multiple_sequence_aligner_undergrad_assignment import read_reference, build_index
from virus_assembly import read_data, compare_reads
from collections import Counter
import re


# Calculates the coverage of each exon. Coverage is defined as the number of reads whose alignment starts within the exon's location.
def calculate_coverage(reads, reference, exons_fn, max_mismatches, num_pieces):
    '''
    :param reads: List of single-end reads
    :param reference: Reference sequence
    :param exons_fn: Filename of file containing exon locations
    :param max_mismatches: Maximum number of mismatches to allow when aligning reads
    :param num_pieces: Number of pieces to divide each read into
    :return: Matrix of exon coverage. Two columns: First column is exon name. Second column is number of reads that align to exon.
    '''
    length_read = len(reads[0])
    length_index = int(length_read / num_pieces)
    exon_locations = []
    coverage = Counter()
    with open(exons_fn, 'r') as f:
        next(f)
        for line in f:
            exon_locations.append(re.split(':|,', line.strip()))

    index = build_index(ref, length_index)
    for read in reads:
        ind = compare_reads(reference, read, index, length_index, max_mismatches, num_pieces)
        if ind is not None:
            exon = find_match(ind, exon_locations)
            if exon is not None:
                coverage.update([exon])
    exon_counts = coverage.most_common()
    exon_counts_df = pd.DataFrame(exon_counts).rename(columns={0: 'exon', 1: 'exon_count'})

    return exon_counts_df


# Detects whether a location falls within an exon
def find_match(ind, exon_locations):
    '''
    :param ind: The location
    :param exon_locations: List of exons names and their start and end locations
    :return: The exon name if the location matches
    '''
    for exon, start, end in exon_locations:
        if ind >= int(start) and ind <= int(end):
            return exon
    return None


# Gets length of each exon from file containing exon start and end locations
def read_exon_lengths(exons_fn):
    '''
    :param exons_fn: File name
    :return: Matrix of exon lengths. Two columns: First column has name of exon. Second column has length of exon.
    '''
    exon_length_dict = {}
    with open(exons_fn) as exons_file:
        exons_file.readline()
        for line in exons_file:
            exon, start_end = line.strip().split(':')
            start, end = (int(x) for x in start_end.split(','))
            exon_length = end - start
            exon_length_dict[exon] = exon_length
    exon_length_df = pd.DataFrame.from_dict(exon_length_dict, orient='index').reset_index().rename(columns={'index': 'exon', 0: 'length'})
    return exon_length_df


# Reads isoform information (which isoforms contain which exons) from file.
def read_isoforms(isoforms_fn):
    '''
    :param isoforms_fn: File name
    :return: Matrix of isoform information. Two columns: First column contains name of isoform. Second column contains
             name of a single exon found in the isoform
    '''
    isoform_exons = []
    with open(isoforms_fn) as isoforms_file:
        isoforms_file.readline()
        for line in isoforms_file:
            isoform, exons = line.strip().split(':')
            exons = exons.split(',')
            gene_id = isoform.split('_ISO')[0]
            full_exons = ['{}_{}'.format(gene_id, exon) for exon in exons]
            these_exons = [(isoform, exon) for exon in full_exons]
            isoform_exons += these_exons

    isoforms_df = pd.DataFrame(isoform_exons).rename(columns={0: 'isoform', 1: 'exon'})
    return isoforms_df


# Reads exon counts from file
def read_exon_counts(exon_counts_fn):
    '''
    :param exon_counts_fn: File name
    :return: Matrix of exon counts. Two columns: First column is exon name. Second column is number of reads that align to exon.
    '''
    exon_counts = []
    with open(exon_counts_fn) as exon_counts_file:
        exon_counts_file.readline()
        for line in exon_counts_file:
            exon, count = line.strip().split(':')
            exon_counts.append((exon, int(count)))
    exon_counts_df = pd.DataFrame(exon_counts).rename(columns={0: 'exon', 1: 'exon_count'})
    return exon_counts_df



if __name__ == "__main__":
    input_folder = 'hw5_M_2'
    input_fn_start = '/Users/douglasyao/Downloads/{0}/{0}'.format(input_folder)
    exons_fn = '{}_exons.txt'.format(input_fn_start)
    isoforms_fn = '{}_isoforms.txt'.format(input_fn_start)
    exon_counts_fn = '{}_exon_counts.txt'.format(input_fn_start)
    reads_fn = '{}_reads.txt'.format(input_fn_start)
    ref_fn = '/Users/douglasyao/Downloads/ref_hw5.txt'

    reads = read_data(reads_fn)
    ref = read_reference(ref_fn)

    exon_length_df = read_exon_lengths(exons_fn)
    isoforms_df = read_isoforms(isoforms_fn)
    isoform_exon_df = pd.merge(isoforms_df, exon_length_df, on='exon')

    combined_df = pd.pivot_table(isoform_exon_df, columns='isoform', index='exon', values='length').fillna(0)

    isoform_columns = combined_df.columns
    combined_df.reset_index(inplace=True)

    exon_counts_df = calculate_coverage(reads, ref, exons_fn, 2, 2)
    complete_df = pd.merge(combined_df, exon_counts_df, on='exon')

    exon_counts_array = complete_df.exon_count.values
    exon_isoform_matrix = complete_df.loc[:, isoform_columns].values

    implied_frequencies = np.linalg.lstsq(exon_isoform_matrix, exon_counts_array)[0]
    implied_frequencies /= implied_frequencies.sum()

    output_fn = '{}_output.txt'.format(input_fn_start)
    with open(output_fn, 'w') as output_file:
        output_file.write('>{}\n'.format(input_folder))
        for isoform, freq in sorted([_ for _ in zip(isoform_columns, implied_frequencies)], key=lambda x: (int(x[0].split('_')[1]), int(x[0].split('_')[3]))):
            output_file.write('{}:{}\n'.format(isoform, freq))
