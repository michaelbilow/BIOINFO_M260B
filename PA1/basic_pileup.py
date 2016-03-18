from collections import defaultdict
import sys
from os.path import join
import zipfile

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
    input_folder = join('../data', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    input_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
    snps, lines = generate_consensus(input_fn)
    output_fn = join(input_folder, 'snps_{}.txt'.format(f_base))
    zip_fn = join(input_folder, 'snps_{}.zip'.format(f_base))
    with open(output_fn, 'w') as output_file:
        header = '>{}\n>{}\n'.format(f_base, 'SNP')
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

    output_fn2 = join(input_folder, 'consensus_{}.txt'.format(f_base))
    output_file2 = open(output_fn2, 'w')
    output_file2.write('\n'.join(lines))
