"""
Created on June 20, 2014

@author: Joseph Korpela
"""

from __future__ import absolute_import
import mmap
import random
import sys
import zipfile
import zlib
import re
import argparse
import pickle
import shutil
import unittest
import os
from heapq import merge
import helpful_scripts
import hash_based_sequencer
import Eval

class chromosome_builder():
    def __init__(self, args=None):
        if not args:
            args = self.parse_system_args()

        self._genome_id = args.id + "_chr_" + str(args.chr_id)
        self._chromosome_id = args.chr_id
        self._chromosome_size = args.chr_size

        if args.verbose == 'y':
            self._verbose = True
        else:
            self._verbose = False

        if args.simple == 'y':
            self._simple = True
        else:
            self._simple = False

        if args.scale == 'k':
            self._chromosome_size *= 1000
        elif args.scale == 'm':
            self._chromosome_size *= 1000000
        elif args.scale == 'b':
            self._chromosome_size *= 1000000000

        if args.alu == 'y':
            self._use_alu = True
            self._base_alu = args.base_alu
        else:
            self._use_alu = False

        if args.assembly == 'y':
            self._use_assembly = True
        else:
            self._use_assembly = False

        if args.spectrum == 'y':
            self._use_spectrum = True
        else:
            self._use_spectrum = False

        self._allele_base_list = ["C", "T", "G", "A"]

        self._working_dir = "TMP_" + str(self._genome_id)
        self._ref_genome_file = "ref_" + str(self._genome_id) + ".txt"
        self._priv_genome_file = "donor_" + str(self._genome_id) + ".txt"
        self._reads_file = "reads_" + str(self._genome_id) + ".txt"
        self._answer_file = "ans_" + str(self._genome_id) + ".txt"
        self._base_alu_file = "alu_" + str(self._genome_id) + ".txt"

        self._overlap_buffer = 5
        self._long_variant_rate = .1

        self._snp_rate = 0.003

        self._ref_str_rate = 0.000075
        self._denovo_str_rate = 0.000025
        self._str_min_copies = 5
        self._str_max_copies = 50
        self._str_min_length = 2
        self._str_max_length = 5
        self._str_mutation_amount = 2

        self._ref_cnv_rate = 0.0001
        self._denovo_cnv_rate = 0.00001
        self._cnv_min_length = 20
        self._cnv_max_length = 500
        self._cnv_min_copies = 2
        self._cnv_max_copies = 10
        self._cnv_mutation_amount = 2

        self._inv_rate = 0.00001
        self._inv_short_min_length = 20
        self._inv_short_max_length = 50
        self._inv_long_min_length = 50
        self._inv_long_max_length = 500

        self._ins_rate = 0.0005
        self._ins_short_min_length = 1
        self._ins_short_max_length = 5
        self._ins_long_min_length = 5
        self._ins_long_max_length = 200

        self._del_rate = 0.0005
        self._del_short_min_length = 1
        self._del_short_max_length = 5
        self._del_long_min_length = 5
        self._del_long_max_length = 200

        self._alu_mutation_rate = 0.3
        self._alu_min_length = 300
        self._alu_max_length = 300

        if self._use_alu:
            self._ref_alu_rate = 0.075
            self._denovo_alu_rate = 0.025
        else:
            self._ref_alu_rate = 0
            self._denovo_alu_rate = 0

        #reduce the max length of mutations for smaller chromosome sizes
        if self._chromosome_size < 500000:
            self._nbr_long_inv = 0
            self._nbr_long_ins = 0
            self._nbr_long_del = 0

            self._str_max_copies = 20 #from 50
            self._cnv_max_length = 50 #from 500
            self._cnv_max_copies = 4 #from 10
        else:
            self._nbr_long_inv = max(4, int(self._inv_rate * self._chromosome_size * self._long_variant_rate))
            self._nbr_long_ins = max(10, int(self._ins_rate * self._chromosome_size * self._long_variant_rate))
            self._nbr_long_del = max(10, int(self._ins_rate * self._chromosome_size * self._long_variant_rate))

        self._nbr_snp = int(self._snp_rate * self._chromosome_size)
        self._nbr_denovo_str = int(self._denovo_str_rate * self._chromosome_size)
        self._nbr_denovo_cnv = int(self._denovo_cnv_rate * self._chromosome_size)
        self._nbr_ref_alu = int(self._ref_alu_rate * self._chromosome_size / float(self._alu_max_length))
        self._nbr_denovo_alu = int(self._denovo_alu_rate * self._chromosome_size / float(self._alu_max_length))

        self._nbr_ref_str = max(4, int(self._ref_str_rate * self._chromosome_size))
        self._nbr_ref_cnv = max(4, int(self._ref_cnv_rate * self._chromosome_size))
        self._nbr_short_inv = max(4, int(self._inv_rate * self._chromosome_size * (1 - self._long_variant_rate)))
        self._nbr_short_ins = max(10, int(self._ins_rate * self._chromosome_size * (1 - self._long_variant_rate)))
        self._nbr_short_del = max(10, int(self._ins_rate * self._chromosome_size * (1 - self._long_variant_rate)))

        if self._simple:
            self._nbr_long_inv = 0
            self._nbr_long_ins = 0
            self._nbr_long_del = 0
            self._nbr_denovo_str = 0
            self._nbr_denovo_cnv = 0
            self._nbr_ref_alu = 0
            self._nbr_denovo_alu = 0
            self._nbr_short_inv = 0
            self._nbr_ref_str = 0

        if self._use_assembly:
            self._nbr_long_inv = 0
            self._nbr_long_ins = 0
            self._nbr_long_del = 0
            self._nbr_snp = 0
            self._nbr_denovo_str = 0
            self._nbr_denovo_cnv = 0
            self._nbr_ref_alu = 0
            self._nbr_denovo_alu = 0
            self._nbr_short_inv = 0
            self._nbr_short_ins = 0
            self._nbr_short_del = 0
            self._cnv_mutation_amount = 0

            if self._chromosome_size < 50000:
                self._nbr_ref_str = max(4, int(self._ref_str_rate * self._chromosome_size))
                self._nbr_ref_cnv = max(4, int(self._ref_cnv_rate * self._chromosome_size))

                self._str_max_copies = 20 #from 50
                self._cnv_max_length = 300 #from 500
                self._cnv_max_copies = 4 #from 10
            else:
                self._nbr_ref_str = int(self._ref_str_rate * self._chromosome_size)
                self._nbr_ref_cnv = int(self._ref_cnv_rate * self._chromosome_size)

                self._str_max_copies = 50
                self._cnv_max_length = 500
                self._cnv_max_copies = 6

        #mutation_list is used when generating the various mutations for the genomes
        self._mutation_list = []
        self._str_list = [] #used when mutating STRs in donor genome
        self._cnv_list = [] #used when mutating CNVs in donor genome
        self._cnv_dict = {}

        if self._use_assembly and not self._use_spectrum:
            self._sequencer_coverage = 100
        else:
            self._sequencer_coverage = 30

        self._sequencer_error_rate = 0.01
        self._sequencer_garbage_rate = 0.1
        self._sequencer_read_length = 50
        self._sequencer_gap_min = 90
        self._sequencer_gap_max = 110

    def insert_newlines(self, sequence, line_size=80):
        return '\n'.join(sequence[i:i+line_size] for i in range(0, len(sequence), line_size)) + '\n'

    def write_genome_lines_to_file(self, genome, file_object):
        genome = self.insert_newlines(genome, 80)
        file_object.write(genome)

    def parse_fasta(self, file_name, buffer_size=100000):
        """Gives buffered access to large fasta files so that the entire file doesn't need to be loaded into
        memory all at once.  Works as a generator, yielding a block of up to buffer_size with each call. For
        general use, use:
        for sequence in parse_fasta(file_name, buffer_size)
        This yield sequences until the end of file or a '>' character is found, at which point it yields None
        Since None is yielded for '>', this can be used with multiple chromosomes separated by '>chr#' in a single
        file.  To do so, the generator should be initialized before iterating through chromosomes, then as each
        chromosome is processed you can anticipate None will be yielded one time to mark the end of the current
        chromoeome
        :param file_name: the file to read in
        :param buffer_size: the number of characters to return for each iteration
        :returns: Sequences of up to size buffer_size, or None if EOF or '>' is encountered
        """
        with open(file_name) as fasta_file:
            #start_of_file = True
            buffer = ""
            while True:
                for line in fasta_file:
                    #skip initial documentation lines
                    if '>' in line:
                        pass
                    else:
                        buffer += line.strip()
                        if len(buffer) >= buffer_size:
                            yield buffer[:buffer_size]
                            buffer = buffer[buffer_size:]
                #clear out any remaining buffer when the file is done
                if len(buffer) > 0:
                    yield buffer
                    buffer = ''
                else:
                    yield None

    #this version may give a slight performance boost, but need to work out the bugs before it can be used
    def parse_fasta_mmap(self, file_name, buffer_size=100000):
        with open(file_name, encoding='utf-8') as fasta:
            fasta_map = mmap.mmap(fasta.fileno(), 0, access=mmap.ACCESS_READ)
            start_of_file = True
            buffer = ""
            for line in fasta_map:
                line = line.decode('utf-8')
                #skip initial documentation lines
                if start_of_file and '>' in line:
                    pass
                #each chromosome is marked by a > line, so need to catch this switch
                elif not start_of_file and '>' in line:
                    if len(buffer) == 0:
                        yield None
                    else:
                        #first yield the buffer, then yeild None to flag the end of the chromosome
                        yield buffer
                        buffer = ''
                        yield None
                else:
                    if start_of_file:
                        start_of_file = False
                    buffer += line.strip()
                    if len(buffer) >= buffer_size:
                        yield buffer[:buffer_size]
                        buffer = buffer[buffer_size:]
            #clear out any remaining buffer when the file is done
            yield buffer

    def compare_intervals(self, stt1, end1, stt2, end2, buffer_space=0):
        """
        Compares two intervals represented by their start and end posns to check which precedes the other,
        or if they overlap.  Adds a buffer space around each interval that increases the region in which they
        are considered to overlap.
        Returns: -1 If interval 1 (stt1 and end1) precedes interval 2, 0 if they overlap, or 1 if interval 2
                 precedes interval 1
        """
        stt1 -= buffer_space
        end1 += buffer_space
        stt2 -= buffer_space
        end2 += buffer_space

        if end1 < stt2:
            return -1
        elif end2 < stt1:
            return 1
        else:
            return 0

    def find_empty_ranges(self, range_size, nbr_posns, buffer_size):
        """
        Searches for ranges of unused posns in the genome for use when introducing new structural variants.
        Finds the number of posns given as a parameter and returns them as a list.
        """
        posn_list = []
        max_posn = self._chromosome_size - range_size - 1

        #Will repeat until enough positions have been found
        while len(posn_list) < nbr_posns:
            raw_posn_list = []
            for posn in posn_list:
                raw_posn_list.append(posn)
            posn_list = []

            #1. Generate 150% needed number of random positions
            for i in range(int(nbr_posns)):
                raw_posn_list.append(random.randint(0, max_posn))

            #2. Sort those positions and then check each to find whether they will overlap a preexisting
            #structural variant.
            raw_posn_list.sort()
            overlap_idx = 0

            #first check that there is no overlap among the generated positions
            last_end = raw_posn_list[0] + range_size + (2 * buffer_size)
            tmp_posn_list = []
            tmp_posn_list.append(raw_posn_list[0])
            for i in range(1, len(raw_posn_list)):
                raw_posn = raw_posn_list[i]
                if raw_posn > last_end:
                    tmp_posn_list.append(raw_posn)
                    last_end = raw_posn + range_size + (2 * buffer_size)
                else:
                    new_posn = last_end + 1
                    new_end = new_posn + range_size + (2 * buffer_size)
                    if new_posn < max_posn:
                        if i == len(raw_posn_list) - 1 or new_end < raw_posn_list[i+1]:
                            tmp_posn_list.append(new_posn)
                            last_end = new_end

            raw_posn_list = tmp_posn_list
            tmp_posn_list = None

            if len(self._mutation_list) == 0:
                posn_list = raw_posn_list
            else:
                #then check that the remaining positions do not overlap existing structural variants
                for i in range(len(raw_posn_list)):
                    raw_posn = raw_posn_list[i]
                    while overlap_idx < len(self._mutation_list):

                        ovlp_stt = self._mutation_list[overlap_idx][0]
                        ovlp_end = self._mutation_list[overlap_idx][1]
                        compare_result = self.compare_intervals(raw_posn, raw_posn+range_size,
                                                                ovlp_stt, ovlp_end,
                                                                buffer_size)

                        #no overlap
                        if compare_result == -1:
                            posn_list.append(raw_posn)
                            break
                        #attempt to shift this interval down, if that doesn't work, then
                        #ignore this position as it overlaps a preexisting position
                        elif compare_result == 0:
                            if overlap_idx > 0:
                                prev_end1 = self._mutation_list[overlap_idx-1][1]
                                prev_end2 = raw_posn_list[i-1]+range_size
                                prev_end = max(prev_end1, prev_end2)

                                new_posn = prev_end + (2*buffer_size)
                                if new_posn + range_size + (2 * buffer_size) < ovlp_stt:
                                    posn_list.append(new_posn)
                            break
                        #no overlap was found, move to the next position in the mutation list to check for overlap
                        elif compare_result == 1:
                            if overlap_idx < len(self._mutation_list) - 1:
                                overlap_idx += 1
                            else:
                                posn_list.append(raw_posn)
                                break

        #3. If there are too many positions, then randomly removes some to reduce list to proper size
        while len(posn_list) > nbr_posns:
            del posn_list[random.randint(0, len(posn_list)-1)]
        return posn_list

    def random_sequence(self, seq_len):
        return "".join(random.choice(self._allele_base_list) for i in range(seq_len))

    def delete_block(self, sequence, index, size):
        """deletes a block of items from a given sequence
        :param sequence: sequence from which to delete items
        :param index: the first position to delete
        :param size: the total number of positions to delete, may extend beyond the end of the sequence
        :returns: modified sequence with block deleted
        """
        if index < 0 and index + size > -1:
            return sequence[:index]
        else:
            return sequence[:index] + sequence[index + size:]

    def insert_block(self, sequence, index, new_block):
        """inserts a block of items into a given sequence
        :param sequence: sequence into which to insert items
        :param index: the position before which to begin the insertion, to append to end use index = len(sequence)
        :param new_block: the items to be inserted
        :returns: modified sequence with block inserted
        """
        return sequence[:index] + new_block + sequence[index:]

    def overwrite_block(self, sequence, index, new_block):
        """overwrites a block of items in a given sequence
        :param sequence: sequence in which to overwrite items
        :param index: the position at which to begin overwriting, to append to end use index = len(sequence)
        :param new_block: the items which will be written, may extend beyond end of original sequence
        :returns: modified sequence with block overwritten
        """
        if (index < 0 and index + len(new_block) > -1) or (index + len(new_block) > len(sequence) - 1):
            return sequence[:index] + new_block
        else:
            return sequence[:index] + new_block + sequence[index + len(new_block):]

    def invert_block(self, sequence, index, size):
        """inverts a block of items in a given sequence
        :param sequence: sequence in which to invert items
        :param index: the position at which to begin inversion
        :param size: the number of items which will be inverted
        :returns: modified sequence with block overwritten, the original block, and the inverted block
        """
        if index < 0:
            stt_idx = len(sequence) + index
        else:
            stt_idx = index
        end_idx = min(stt_idx + size, len(sequence))
        original_block = sequence[stt_idx:end_idx]
        inverted_block = original_block[::-1]
        sequence_with_inversion = self.overwrite_block(sequence, stt_idx, inverted_block)
        return sequence_with_inversion, original_block, inverted_block

    def generate_snp_allele(self, orig_allele):
        allele_list = ["A", "C", "G", "T"]
        allele_list.remove(orig_allele)
        return random.choice(allele_list)

    def generate_str_base(self, seq_len):
        str_seq = ''
        while len(str_seq) == 0:
            str_seq = self.random_sequence(seq_len)

            invalid = True
            #ensure the sequence is not all the same allele
            for idx in range(1, len(str_seq)):
                if str_seq[idx - 1] != str_seq[idx]:
                    invalid = False
            #if the first half of the sequence matches the second half, then consider it invalid
            if not invalid and str_seq[:int(len(str_seq)/2)] == str_seq[int(len(str_seq)/2):]:
                invalid = True
            #if the sequence was invalid, then clear it out and try again
            if invalid:
                str_seq = ''
        return str_seq

    def generate_alu_sequence(self, length):
        if not self._base_alu or len(self._base_alu) == 0:
            raise Exception("No base Alu defined")
        new_alu = self._base_alu
        len_diff = abs(length - len(new_alu))
        for i in range(length - len(new_alu)):
            new_alu = self.insert_block(new_alu,
                                        random.randint(0, len(new_alu)-1),
                                        random.choice(self._allele_base_list))
        for i in range(len(new_alu) - length):
            new_alu = self.delete_block(new_alu,
                                        random.randint(0, len(new_alu)-1),
                                        1)
        for k in range(int(len(new_alu)*self._alu_mutation_rate)-len_diff):
            alu_posn = random.randint(0, len(new_alu)-1)
            snp = self.generate_snp_allele(new_alu[alu_posn])
            new_alu = self.overwrite_block(new_alu, alu_posn, snp)
        return new_alu

    def ranged_length_list(self, min_len, max_len, nbr_items):
        """generates a list of lengths that vary in size between min_len and max_len
        :param min_len: smallest value to return
        :param max_len: largest value to return, at least 2 items will have this value, must be >= min_len
        :param nbr_items: the number of items which will be returned, must be > 0
        :returns: a list of lengths with nbr_items items that vary from min_len to max_len
        """
        if nbr_items < 1:
            raise Exception("Minimum length for the list is 1")
        if max_len < min_len:
            raise Exception("max_len must be greater than or equal to min_len")
        length_list = []

        if nbr_items > 1:
            max_items = max(2, int(nbr_items/10.0))
        else:
            max_items = 1

        for i in range(max_items):
            length_list.append(max_len)

        if nbr_items > 2:
            below_max = nbr_items - max_items
            length_range = max_len - min_len

            for i in range(below_max):
                adj_value = random.randint(0, i) / float(below_max)
                length_list.append(int(min_len + (length_range * adj_value)))

        return length_list

    def mutate_str(self):
        temp_mut_list = []
        for j in range(len(self._str_list)):
            mutation_amount = random.randint(-self._str_mutation_amount, self._str_mutation_amount)
            orig_str = self._str_list[j][1] * self._str_list[j][2]
            new_str = self._str_list[j][1] * (self._str_list[j][2] + mutation_amount)
            str_stt = self._str_list[j][0]
            str_end = self._str_list[j][0] + len(orig_str)
            #self._mutation_list.append([str_stt, str_end, 'MUT_STR', new_str])
        #self._mutation_list.sort()
            temp_mut_list.append([str_stt, str_end, 'MUT_STR', new_str])
        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def mutate_cnv(self):
        new_posn_list = self.find_empty_ranges(self._cnv_max_length,
                                               self._nbr_ref_cnv * self._cnv_mutation_amount,
                                               self._overlap_buffer)
        temp_mut_list = []
        for i in range(self._nbr_ref_cnv):
            cnv_id = self._cnv_list[i][0]
            cnv_len = self._cnv_list[i][1]
            cnv_posn_list = self._cnv_list[i][2]

            mutation_amount = random.randint(-self._cnv_mutation_amount, self._cnv_mutation_amount)

            if mutation_amount > 0:
                for i in range(mutation_amount):
                    cnv_posn = new_posn_list.pop(random.randint(0, len(new_posn_list)-1))
                    cnv_posn_list.append(cnv_posn)
                cnv_posn_list.sort()

            if mutation_amount < 0:
                for j in range(mutation_amount):
                    next_posn = cnv_posn_list.pop(random.randint(0, len(cnv_posn_list)-1))
                    cnv_stt = next_posn
                    cnv_end = next_posn + cnv_len
                    self._mutation_list.append([cnv_stt, cnv_end, 'DEL_CNV', cnv_id])

            for cnv_posn in cnv_posn_list:
                #self._mutation_list.append([cnv_posn, cnv_posn + cnv_len, "DONOR_CNV", cnv_id])
        #self._mutation_list.sort()
                temp_mut_list.append([cnv_posn, cnv_posn + cnv_len, "DONOR_CNV", cnv_id])
        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))


    def allocate_cnv(self, nbr_cnv, variant_tag):
        if nbr_cnv < 1:
            return
        cnv_length_list = self.ranged_length_list(self._cnv_min_length, self._cnv_max_length, nbr_cnv)
        cnv_posn_list = self.find_empty_ranges(self._cnv_max_length,
                                               nbr_cnv * self._cnv_max_copies,
                                               self._overlap_buffer)
        temp_mut_list = []
        for i in range(nbr_cnv):
            posn_list = []
            seq_len = cnv_length_list[i]
            nbr_copies = random.randint(self._cnv_min_copies, self._cnv_max_copies)

            if 'REF' in variant_tag:
                cnv_id = i
            else:
                cnv_id = i + self._nbr_ref_cnv

            for j in range(nbr_copies):
                cnv_posn = cnv_posn_list.pop(random.randint(0, len(cnv_posn_list)-1))
                posn_list.append(cnv_posn)
                #self._mutation_list.append([cnv_posn, cnv_posn+seq_len, variant_tag, cnv_id])
                temp_mut_list.append([cnv_posn, cnv_posn+seq_len, variant_tag, cnv_id])
            self._cnv_list.append([cnv_id, seq_len, posn_list])
        #self._mutation_list.sort()
        temp_mut_list.sort()
        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def allocate_alu(self, nbr_alu, variant_tag):
        if nbr_alu < 1:
            return
        alu_length_list = self.ranged_length_list(self._alu_min_length, self._alu_max_length, nbr_alu)
        alu_posn_list = self.find_empty_ranges(self._alu_max_length,
                                               nbr_alu,
                                               self._overlap_buffer)
        temp_mut_list = []
        for j in range(nbr_alu):
            alu_stt = alu_posn_list[j]
            alu_len = alu_length_list.pop(random.randint(0, len(alu_length_list)-1))
            #donor alus are inserted into the genome, so their end_posn in reference to ref genome is their start
            if variant_tag == 'DONOR_ALU':
                alu_end = alu_stt
            else:
                alu_end = alu_posn_list[j] + alu_len
            #self._mutation_list.append([alu_stt, alu_end, variant_tag, alu_len])
        #self._mutation_list.sort()
            temp_mut_list.append([alu_stt, alu_end, variant_tag, alu_len])
        temp_mut_list.sort()
        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def allocate_str(self, nbr_str, variant_tag):
        if nbr_str < 1:
            return
        str_posn_list = self.find_empty_ranges(self._str_max_copies * self._str_max_length,
                                               nbr_str,
                                               self._overlap_buffer)

        temp_mut_list = []
        for i in range(nbr_str):
            seq_len = random.randint(self._str_min_length, self._str_max_length)
            nbr_copies = random.randint(self._str_min_copies, self._str_max_copies)
            str_seq = self.generate_str_base(seq_len)

            str_posn = str_posn_list.pop(random.randint(0, len(str_posn_list)-1))
            #self._mutation_list.append([str_posn, str_posn + (seq_len*nbr_copies), variant_tag, str_seq, nbr_copies])
            self._str_list.append([str_posn, str_seq, nbr_copies])
        #self._mutation_list.sort()
            temp_mut_list.append([str_posn, str_posn + (seq_len*nbr_copies), variant_tag, str_seq, nbr_copies])

        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def allocate_inversions(self, nbr_inv, min_len, max_len):
        if nbr_inv < 1:
            return
        inv_length_list = self.ranged_length_list(min_len, max_len, nbr_inv)
        inv_posn_list = self.find_empty_ranges(max_len,
                                               nbr_inv,
                                               self._overlap_buffer)
        temp_mut_list = []
        for i in range(nbr_inv):
            inv_stt = inv_posn_list[i]
            inv_len = inv_length_list.pop(random.randint(0, len(inv_length_list)-1))
            inv_end = inv_stt + inv_len
            #self._mutation_list.append([inv_stt, inv_end, 'INV', inv_len])
        #self._mutation_list.sort()
            temp_mut_list.append([inv_stt, inv_end, 'INV', inv_len])
        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def allocate_insertions(self, nbr_ins, min_len, max_len):
        if nbr_ins < 1:
            return
        ins_length_list = self.ranged_length_list(min_len, max_len, nbr_ins)
        ins_posn_list = self.find_empty_ranges(max_len,
                                               nbr_ins,
                                               self._overlap_buffer)
        temp_mut_list = []
        for i in range(nbr_ins):
            ins_stt = ins_posn_list[i]
            ins_len = ins_length_list.pop(random.randint(0, len(ins_length_list)-1))
            ins_end = ins_stt
            #self._mutation_list.append([ins_stt, ins_end, 'INS', ins_len])
        #self._mutation_list.sort()
            temp_mut_list.append([ins_stt, ins_end, 'INS', ins_len])
        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def allocate_deletions(self, nbr_del, min_len, max_len):
        if nbr_del < 1:
            return
        del_length_list = self.ranged_length_list(min_len, max_len, nbr_del)
        del_posn_list = self.find_empty_ranges(max_len,
                                               nbr_del,
                                               self._overlap_buffer)

        temp_mut_list = []
        for i in range(nbr_del):
            del_stt = del_posn_list[i]
            del_len = del_length_list.pop(random.randint(0, len(del_length_list)-1))
            del_end = del_stt + del_len
            #self._mutation_list.append([del_stt, del_end, 'DEL', del_len])
        #self._mutation_list.sort()
            temp_mut_list.append([del_stt, del_end, 'DEL', del_len])
        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))

    def allocate_snps(self):
        if self._nbr_snp < 1:
            return
        snp_posn_list = self.find_empty_ranges(1, self._nbr_snp, 0)
        temp_mut_list = []

        for i in range(self._nbr_snp):
            snp_stt = snp_posn_list[i]
            snp_end = snp_stt + 1
            temp_mut_list.append([snp_stt, snp_end, 'SNP', 1])
            #self._mutation_list.append([snp_stt, snp_end, 'SNP', 1])

        temp_mut_list.sort()

        self._mutation_list = list(merge(self._mutation_list, temp_mut_list))
        #self._mutation_list.sort()


    def generate_ref_genome(self):
        """
        Generates a random reference genome with the specified number of chromosomes,
        each of length length_chromosome
        """
        if self._use_alu:
            if not self._base_alu or len(self._base_alu) == 0:
                raise Exception("No base Alu defined")
            if not os.path.exists(self._base_alu_file):
                with open(self._base_alu_file, "w") as alu_file:
                    alu_file.write(">" + str(self._genome_id) + "\n")
                    self.write_genome_lines_to_file(self._base_alu, alu_file)
        with open(self._ref_genome_file, "w") as ref_file:
            if not os.path.exists(self._working_dir):
                os.makedirs(self._working_dir)
            ref_file.write(">" + str(self._genome_id) + "\n")

            self._mutation_list = []

            if self._use_alu:
                if self._verbose:
                    print('REF GENOME: Allocating ' + str(self._nbr_ref_alu) + ' Alus')
                self.allocate_alu(self._nbr_ref_alu, "REF_ALU")
            if self._verbose:
                print('REF GENOME: Allocating ' + str(self._nbr_ref_cnv) + ' CNVs')
            self.allocate_cnv(self._nbr_ref_cnv, "REF_CNV")
            if self._verbose:
                print('REF GENOME: Allocating ' + str(self._nbr_ref_str) + ' STRs')
            self.allocate_str(self._nbr_ref_str, "REF_STR")
            if self._verbose:
                print('REF GENOME: Writing to File')

            buffer_adj = 0
            buffer = ''
            mut_idx = 0
            mut_max_idx = len(self._mutation_list) - 1
            buffer_size = 80
            count = 0
            if len(self._mutation_list) == 0:
                mut_idx = -1

            while count < self._chromosome_size:
                if len(buffer) > buffer_size or mut_idx == -1:
                    if mut_idx == -1:
                        skip_distance = self._chromosome_size - count
                        buffer += self.random_sequence(skip_distance)
                        count += skip_distance
                        buffer_size = len(buffer)
                        self.write_genome_lines_to_file(buffer, ref_file)
                    else:
                        self.write_genome_lines_to_file(buffer[:buffer_size], ref_file)
                        buffer = buffer[buffer_size:]
                        buffer_adj += buffer_size
                elif len(self._mutation_list) > 0 and count < self._mutation_list[mut_idx][0]:
                    skip_distance = self._mutation_list[mut_idx][0] - count
                    buffer += self.random_sequence(skip_distance)
                    count += skip_distance
                elif mut_idx != -1:
                    mut_type = self._mutation_list[mut_idx][2]

                    if mut_type == 'REF_STR':
                        str_seq = self._mutation_list[mut_idx][3]
                        nbr_copies = self._mutation_list[mut_idx][4]
                        str_seq = str_seq * nbr_copies
                        #pads either side of str with non matching allele to remove ambiguity
                        if buffer[-1] == str_seq[-1]:
                            buffer = buffer[:-1] + self.generate_snp_allele(buffer[-1])
                        right_padding = self.generate_snp_allele(str_seq[0])
                        buffer += str_seq + right_padding
                        count += len(str_seq) + 1

                    elif mut_type == 'REF_CNV':
                        cnv_stt = self._mutation_list[mut_idx][0]
                        cnv_end = self._mutation_list[mut_idx][1]
                        cnv_len = cnv_end - cnv_stt
                        cnv_id = self._mutation_list[mut_idx][3]
                        if cnv_id in self._cnv_dict:
                            cnv_seq = self._cnv_dict[cnv_id]
                        else:
                            cnv_seq = self.random_sequence(cnv_len)
                            self._cnv_dict[cnv_id] = cnv_seq
                        buffer += cnv_seq
                        count += cnv_len

                    elif mut_type == 'REF_ALU':
                        alu_len = self._mutation_list[mut_idx][3]
                        alu_seq = self.generate_alu_sequence(alu_len)
                        buffer += alu_seq
                        count += alu_len

                    if mut_idx < mut_max_idx:
                        mut_idx += 1
                    else:
                        mut_idx = -1 #flags when all mutations have been seen

    def generate_donor_genome(self):
        with open(self._priv_genome_file, "w") as donor_genome_file:
            donor_genome_file.write(">" + str(self._genome_id) + "\n")

        buffer_size = 1000
        fasta_parser = self.parse_fasta(self._ref_genome_file, buffer_size=buffer_size)
        buffer_size = 80

        #plan out all mutation ranges in reference to the ref genome, storing them in the mutation_list
        if self._verbose:
            print('DONOR GENOME: Mutating existing STRs')
        self.mutate_str()
        if self._verbose:
            print('DONOR GENOME: Mutating existing CNVs')
        self.mutate_cnv()

        if self._use_alu:
            if self._verbose:
                print('DONOR GENOME: Allocating ' + str(self._nbr_denovo_alu) + ' Alus')
            self.allocate_alu(self._nbr_denovo_alu, "DONOR_ALU")
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_denovo_cnv) + ' CNVs')
        self.allocate_cnv(self._nbr_denovo_cnv, "DONOR_CNV")
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_denovo_str) + ' STRs')
        self.allocate_str(self._nbr_denovo_str, "DONOR_STR")
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_long_inv) + ' long inversions')
        self.allocate_inversions(self._nbr_long_inv, self._inv_long_min_length, self._inv_long_max_length)
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_long_ins) + ' long insertions')
        self.allocate_insertions(self._nbr_long_ins, self._ins_long_min_length, self._ins_long_max_length)
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_long_del) + ' long deletions')
        self.allocate_deletions(self._nbr_long_del, self._del_long_min_length, self._del_long_max_length)
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_short_inv) + ' short inversions')
        self.allocate_inversions(self._nbr_short_inv, self._inv_short_min_length, self._inv_short_max_length)
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_short_ins) + ' short insertions')
        self.allocate_insertions(self._nbr_short_ins, self._ins_short_min_length, self._ins_short_max_length)
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_short_del) + ' short deletions')
        self.allocate_deletions(self._nbr_short_del, self._del_short_min_length, self._del_short_max_length)
        if self._verbose:
            print('DONOR GENOME: Allocating ' + str(self._nbr_snp) + ' SNPs')
        self.allocate_snps()
        if self._verbose:
            print('DONOR GENOME: Writing to File')

        variant_types = ['STR','CNV','ALU','INV','INS','DEL','SNP']
        answer_files = {}
        for variant in variant_types:
            answer_files[variant] = open(os.path.join(self._working_dir, variant + '_ANS_FILE'), 'w')
        #read in the reference genome, writing out the donor genome out to file using
        #the mutations from the mutation list
        with open(self._priv_genome_file, "a") as donor_genome_file:
            ref_genome_idx = 0
            buffer_adjust = 0

            donor_genome = ''
            ref_genome = ''

            if len(self._mutation_list) > 0:
                mut_idx = 0
            else:
                mut_idx = -1

            mut_max_idx = len(self._mutation_list) - 1

            while ref_genome_idx + buffer_adjust < self._chromosome_size:

                ref_genome = ref_genome[ref_genome_idx:]
                buffer_adjust += ref_genome_idx
                ref_genome_idx = 0

                next_segment = next(fasta_parser)

                if next_segment:
                    ref_genome += next_segment

                if mut_idx == -1:
                    donor_genome += ref_genome[ref_genome_idx:]
                    ref_genome_idx += len(ref_genome) - ref_genome_idx
                elif ref_genome_idx + buffer_adjust != self._mutation_list[mut_idx][0]:
                    donor_genome += ref_genome[ref_genome_idx:self._mutation_list[mut_idx][0] - buffer_adjust]
                    ref_genome_idx = self._mutation_list[mut_idx][0] - buffer_adjust
                else:
                    while len(donor_genome) > buffer_size:
                        self.write_genome_lines_to_file(donor_genome[:buffer_size], donor_genome_file)
                        donor_genome = donor_genome[buffer_size:]
                    mut_type = self._mutation_list[mut_idx][2]
                    ref_genome_stt = self._mutation_list[mut_idx][0] - buffer_adjust
                    ref_genome_end = self._mutation_list[mut_idx][1] - buffer_adjust
                    ref_genome_idx = ref_genome_end

                    if mut_type == 'SNP':
                        orig_allele = ref_genome[ref_genome_stt]
                        snp_allele = self.generate_snp_allele(orig_allele)
                        donor_genome += snp_allele
                        answer_files['SNP'].write('\n' + orig_allele +
                            ',' + snp_allele + ',' + str(self._mutation_list[mut_idx][0]))

                    #the mutation list contains both the original str and the mutated str, so when
                    #one is encountered the other needs to be pulled and dealt with at the same
                    #time
                    elif mut_type == 'MUT_STR':
                        new_str = self._mutation_list[mut_idx][3]
                        mut_idx += 1
                        donor_genome += new_str
                        answer_files['STR'].write('\n' + new_str +
                                                    ',' + str(self._mutation_list[mut_idx][0]))

                    elif mut_type == 'DONOR_STR':
                        str_seq = self._mutation_list[mut_idx][3]
                        nbr_copies = self._mutation_list[mut_idx][4]
                        str_seq = str_seq * nbr_copies

                        #pads either side of str with non matching allele to remove ambiguity
                        left_padding = self.generate_snp_allele(str_seq[-1])
                        right_padding = self.generate_snp_allele(str_seq[0])
                        padded_str_seq = left_padding + str_seq + right_padding
                        donor_genome += padded_str_seq
                        answer_files['STR'].write('\n' + str_seq +
                                                    ',' + str(self._mutation_list[mut_idx][0] + 1))
                        answer_files['INS'].write('\n' + padded_str_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    elif mut_type == 'REF_CNV':
                        cnv_id = self._mutation_list[mut_idx][3]
                        cnv_seq = ref_genome[ref_genome_stt:ref_genome_end]
                        donor_genome += cnv_seq
                        answer_files['CNV'].write('\n' + str(cnv_id) +
                                ',' + str(self._mutation_list[mut_idx][0]) + ',' + cnv_seq)

                    #assumes DEL_CNV is always followed by an entry for the REF_CNV, so the mut_idx is incremented
                    elif mut_type == 'DEL_CNV':
                        mut_idx += 1
                        del_len = self._mutation_list[mut_idx][1] - self._mutation_list[mut_idx][0]
                        del_seq = ref_genome[ref_genome_stt:ref_genome_end]
                        answer_files['DEL'].write('\n' + del_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    #every non deleted CNV will have a DONOR_CNV entry (some will only have DONOR_CNV, no REF_CNV)
                    elif mut_type == 'DONOR_CNV':
                        cnv_stt = self._mutation_list[mut_idx][0]
                        cnv_end = self._mutation_list[mut_idx][1]
                        cnv_len = cnv_end - cnv_stt
                        cnv_id = self._mutation_list[mut_idx][3]
                        if cnv_id in self._cnv_dict:
                            cnv_seq = self._cnv_dict[cnv_id]
                        else:
                            cnv_seq = ref_genome[ref_genome_stt:ref_genome_end]
                            self._cnv_dict[cnv_id] = cnv_seq
                        donor_genome += cnv_seq

                        if mut_idx < mut_max_idx and self._mutation_list[mut_idx+1][2] == 'REF_CNV' and \
                                self._mutation_list[mut_idx][0] == self._mutation_list[mut_idx+1][0]:
                            mut_idx += 1
                        else:
                            answer_files['INS'].write('\n' + cnv_seq + ',' +
                                    str(self._mutation_list[mut_idx][0]))
                        answer_files['CNV'].write('\n' + str(cnv_id) +
                                ',' + str(self._mutation_list[mut_idx][0]) + ',' + cnv_seq)

                    elif mut_type == 'REF_ALU':
                        alu_stt = self._mutation_list[mut_idx][0]
                        alu_end = self._mutation_list[mut_idx][1]
                        alu_len = alu_end - alu_stt
                        alu_seq = ref_genome[ref_genome_stt:ref_genome_end]
                        donor_genome += alu_seq
                        answer_files['ALU'].write('\n' + alu_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))
                        answer_files['INS'].write('\n' + alu_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    elif mut_type == 'DONOR_ALU':
                        alu_len = self._mutation_list[mut_idx][3]
                        alu_seq = self.generate_alu_sequence(alu_len)
                        donor_genome += alu_seq
                        answer_files['ALU'].write('\n' + alu_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    elif mut_type == 'INV':
                        orig_block = ref_genome[ref_genome_stt:ref_genome_end]
                        inv_block = orig_block[::-1]
                        donor_genome += inv_block
                        answer_files['INV'].write('\n' + orig_block + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    elif mut_type == 'INS':
                        ins_len = self._mutation_list[mut_idx][3]
                        ins_seq = self.random_sequence(ins_len)
                        donor_genome += ins_seq
                        answer_files['INS'].write('\n' + ins_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    elif mut_type == 'DEL':
                        del_seq = ref_genome[ref_genome_stt:ref_genome_end]
                        answer_files['DEL'].write('\n' + del_seq + ',' +
                                str(self._mutation_list[mut_idx][0]))

                    if mut_idx < mut_max_idx:
                        mut_idx += 1
                    else:
                        mut_idx = -1 #flags when all mutations have been seen

                writeable = int(len(donor_genome) / 80)
                if writeable >= 1:
                    self.write_genome_lines_to_file(donor_genome[:writeable*80], donor_genome_file)
                    donor_genome = donor_genome[writeable*80:]
            self.write_genome_lines_to_file(donor_genome, donor_genome_file)

        for key in answer_files:
            answer_files[key].close()

        with open(self._answer_file, 'w') as main_ans:
            main_ans.write(">" + str(self._genome_id))
            for variant in variant_types:
                with open(os.path.join(self._working_dir, variant + '_ANS_FILE'), 'r') as temp_file:
                    if variant == 'CNV':
                        main_ans.write("\n>CNV")
                        cnv_dict = {}
                        for line in temp_file:
                            line = line.strip()
                            if line:
                                line_array = line.split(',')
                                cnv_id = line_array[0]
                                cnv_posn = int(line_array[1])
                                cnv_seq = line_array[2]
                                if cnv_id in cnv_dict:
                                    cnv_seq, cnv_posn_list = cnv_dict[cnv_id]
                                else:
                                    cnv_posn_list = []
                                cnv_posn_list.append(cnv_posn)
                                cnv_dict[cnv_id] = (cnv_seq, cnv_posn_list)
                        cnv_list = []
                        for key in cnv_dict:
                            cnv_seq, cnv_posn_list = cnv_dict[key]
                            cnv_posn_list.sort()
                            cnv_list.append([cnv_seq, cnv_posn_list])
                        cnv_list.sort(key = lambda l: l[:][1])
                        for cnv_seq, cnv_posn_list in cnv_list:
                            main_ans.write('\n' + cnv_seq)
                            for posn in cnv_posn_list:
                                main_ans.write(',' + str(posn))
                    else:
                        main_ans.write("\n>" + variant)
                        for line in temp_file:
                            line = line.strip()
                            if line:
                                main_ans.write('\n' + line)
                os.remove(os.path.join(self._working_dir, variant + '_ANS_FILE'))
        shutil.rmtree(self._working_dir)

    def add_sequencer_errors(self, read_sequence):
        error_list = []
        for i in range(len(read_sequence)):
            if random.random() < self._sequencer_error_rate:
                error_list.append(i)
        for i in error_list:
            orig_allele = read_sequence[i]
            error_allele = self.generate_snp_allele(orig_allele)
            read_sequence = self.overwrite_block(read_sequence, i, error_allele)
        return read_sequence

    def create_read_pair(self, donor_genome, left_stt, right_stt):
        left_read = donor_genome[left_stt:left_stt+self._sequencer_read_length]
        right_read = donor_genome[right_stt:right_stt+self._sequencer_read_length]

        #only one is flipped so they are always in opposing directions with their
        #overall direction
        if random.random() > .5:
            right_read = right_read[::-1]
        else:
            left_read = left_read[::-1]
        return left_read, right_read

    def generate_reads(self):
        with open(self._reads_file, "w") as reads_file:
            reads_file.write(">" + str(self._genome_id))
            with open(self._priv_genome_file, "r") as donor_genome_file:
                # skip the first two '>' labels in the donor genome file
                donor_genome_file.readline()
                temp_file_name_list = []
                donor_genome = ""
                for line in donor_genome_file:
                    if ">" in line:
                        break
                    donor_genome += str(line).strip()

                nbr_reads = int(self._chromosome_size * self._sequencer_coverage / self._sequencer_read_length)
                nbr_pairs = int(nbr_reads / 2)
                write_list = []
                for i in range(nbr_pairs):
                    if random.random() < self._sequencer_garbage_rate:
                        left_read = self.random_sequence(self._sequencer_read_length)
                        right_read = self.random_sequence(self._sequencer_read_length)
                    else:
                        gap_len = random.randint(self._sequencer_gap_min, self._sequencer_gap_max)
                        total_len = 2 * self._sequencer_read_length + gap_len
                        left_stt = random.randint(0, self._chromosome_size - total_len - 1)
                        right_stt = left_stt + self._sequencer_read_length + gap_len

                        left_read, right_read = self.create_read_pair(donor_genome, left_stt, right_stt)

                        left_read = self.add_sequencer_errors(left_read)
                        right_read = self.add_sequencer_errors(right_read)

                    reads_file.write('\n' + left_read + ',' + right_read)

    def generate_spectrum_reads(self):
        with open(self._reads_file, "w") as reads_file:
            reads_file.write(">" + str(self._genome_id))
            with open(self._priv_genome_file, "r") as donor_genome_file:
                # skip the first two '>' labels in the donor genome file
                donor_genome_file.readline()
                temp_file_name_list = []
                donor_genome = ""
                for line in donor_genome_file:
                    if ">" in line:
                        break
                    donor_genome += str(line).strip()

                gap_len = 100
                total_len = 2 * self._sequencer_read_length + gap_len
                for i in range(self._chromosome_size - total_len):
                    left_stt = i
                    right_stt = left_stt + self._sequencer_read_length + gap_len
                    left_read = donor_genome[left_stt:left_stt+self._sequencer_read_length]
                    right_read = donor_genome[right_stt:right_stt+self._sequencer_read_length]

                    reads_file.write('\n' + left_read + ',' + right_read)

    def parse_system_args(self):
        parser = argparse.ArgumentParser(
            description="This script generates a reference and donor genome as a set "
                        "of files. The files can be used for various computational "
                        "genetics purposes. The following files are created: 1) "
                        "reference genome \'ref_*.txt\' 2) mutated donor genome "
                        "\'private_*.txt\' 3) paired-end reads \'reads_*.txt\'"
                        "from donor genome 4) mutation answer key \'ans_*.txt\'"
        )
        parser.add_argument(
            "--id",
            type=str,
            default='test',
            help="The name or ID of this genome for identification purposes. The "
                 "genome id will be reflected in the generated file names."
        )
        parser.add_argument(
            "--chr_id",
            type=int,
            default='1',
            help="The id number for this chromosome, defaults to 1."
        )
        parser.add_argument(
            "--chr_size",
            type=int,
            default='10',
            help="The size of each chromosome, multiplied by -s (scaling factor). Change "
                 "scale with -s option"
        )
        parser.add_argument(
            "-s", "--scale",
            type=str,
            choices=["k", "m", "b"],
            default="k",
            help="the amount to scale chromosome-size by. k: thousands, m: millions,"
                 " b: billions. By default, scale is k (thousands)."
        )
        parser.add_argument(
            "--alu",
            type=str,
            choices=["y", "n"],
            default="n",
            help="whether to include Alus in the genome."
        )
        parser.add_argument(
            "--assembly",
            type=str,
            choices=["y", "n"],
            default="n",
            help="whether to generate output for assembly (no reference genome)."
        )
        parser.add_argument(
            "--spectrum",
            type=str,
            choices=["y", "n"],
            default="n",
            help="whether to generate reads as spectrum (in order from position 1 up)."
        )
        parser.add_argument(
            "--simple",
            type=str,
            choices=["y", "n"],
            default="n",
            help="whether to suppress variants other than SNP, INDEL, and CNV, used for early hw assignments"
        )
        parser.add_argument(
            "--verbose",
            type=str,
            choices=["y", "n"],
            default="n",
            help="whether to print status updates while building the genome."
        )
        return parser.parse_args()

class TestClass(unittest.TestCase):

    def setUp(self):
        args = TestSetting()
        args.id = 'test'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'y'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.base_alu = random_sequence(300)
        self.gen = chromosome_builder(args)

        self.scripts = helpful_scripts.helpful_scripts()
        self.seq_ref_file = 'ref_check_reference_file.txt'
        self.seq_reads_file = 'ref_check_reads_file.txt'
        self.seq_ans_file =  'ref_check_ans_file.txt'
        self.seq_temp_ans_file = 'tmp_ans_file.txt'

        self.data_directory = 'test_data'

    def test_using_ref_check(self):
        args = TestSetting()
        args.id = 'test'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'y'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.base_alu = random_sequence(300)
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        ans_key = gen._answer_file
        ref_file = gen._ref_genome_file

        result = self.scripts.ref_check(ans_key, ref_file)
        self.assertEqual(result, "No mismatches")

        os.remove(gen._ref_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._reads_file)

    def test_using_hash_based_assembly(self):
        args = TestSetting()
        args.id = 'test'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'y'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.base_alu = random_sequence(300)
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_ans = self.seq_temp_ans_file
        ans_key = gen._answer_file
        ref_file = gen._ref_genome_file
        reads_file = gen._reads_file
        with open(ref_file, 'r') as reference:
            with open(reads_file, 'r') as reads:
                with open(student_ans, 'w') as student:
                    sequencer = hash_based_sequencer.MultiReSequencer(2, 10, reference, reads, student)
                    sequencer.ref_genome_map_str(6)
                    sequencer.process_reads()
                    sequencer.create_answer_file()

        eval = Eval.Eval()
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = eval.eval(answer, student)
                self.assertTrue(grades['SNP'] > .50)
                self.assertTrue(grades['SNP'] < 1)
                self.assertTrue(grades['STR'] < .95)
                self.assertTrue(grades['STR'] > .50)

                self.assertEqual(grades['INDEL'], 0)
                self.assertEqual(grades['CNV'], 0)
                self.assertEqual(grades['INV'], 0)
                self.assertEqual(grades['ALU'], 0)

                self.assertEqual(grades['A_COV'], -1)
                self.assertEqual(grades['A_ACC'], -1)
                self.assertEqual(grades['A_CON'], -1)
        os.remove(student_ans)
        os.remove(gen._ref_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._reads_file)

    def test_compare_intervals(self):
        self.assertEqual(0, self.gen.compare_intervals(0, 1, 1, 2, buffer_space=0))
        self.assertEqual(-1, self.gen.compare_intervals(0, 0, 1, 10, buffer_space=0))
        self.assertEqual(1, self.gen.compare_intervals(11, 12, 1, 10, buffer_space=0))
        self.assertEqual(0, self.gen.compare_intervals(0, 0, 1, 10, buffer_space=1))
        self.assertEqual(0, self.gen.compare_intervals(1, 10, 0, 0, buffer_space=1))
        self.assertEqual(0, self.gen.compare_intervals(11, 12, 1, 10, buffer_space=1))
        self.assertEqual(0, self.gen.compare_intervals(0, 4, 1, 2, buffer_space=0))
        self.assertEqual(0, self.gen.compare_intervals(0, 4, 1, 2, buffer_space=5))
        self.assertEqual(0, self.gen.compare_intervals(1, 4, 0, 2, buffer_space=0))
        self.assertEqual(0, self.gen.compare_intervals(1, 4, 1, 4, buffer_space=0))

    def test_find_empty_ranges(self):
        self.gen._mutation_list.append([1000, 10000, 'MUT_STR', 9000])
        posn_list = self.gen.find_empty_ranges(10, 25, 5)
        self.assertEqual(len(posn_list), 25)
        for posn in posn_list:
            self.assertTrue(posn < 995 and posn >= 0)

        self.gen._mutation_list = []
        self.gen._mutation_list.append([0, 9000, 'MUT_STR', 9000])
        posn_list = self.gen.find_empty_ranges(10, 25, 5)
        self.assertEqual(len(posn_list), 25)
        for posn in posn_list:
            self.assertTrue(posn <= 10000 and posn > 9005)

        self.gen._mutation_list = []
        self.gen._mutation_list.append([9000, 10000, 'MUT_STR', 1000])
        self.gen._mutation_list.append([0, 1000, 'MUT_STR', 1000])
        self.gen._mutation_list.sort()
        posn_list = self.gen.find_empty_ranges(10, 100, 5)
        self.assertEqual(len(posn_list), 100)
        for posn in posn_list:
            self.assertTrue(posn < 8995 and posn > 1005)

        self.gen._mutation_list = []
        self.gen._mutation_list.append([0, 1000, 'MUT_STR', 1000])
        self.gen._mutation_list.append([1100, 2000, 'MUT_STR', 900])
        self.gen._mutation_list.append([4000, 5000, 'MUT_STR', 1000])
        self.gen._mutation_list.append([7000, 8000, 'MUT_STR', 1000])
        self.gen._mutation_list.append([9000, 10000, 'MUT_STR', 1000])
        self.gen._mutation_list.sort()
        posn_list = self.gen.find_empty_ranges(10, 25, 5)
        self.assertEqual(len(posn_list), 25)
        for posn in posn_list:
            self.assertTrue( (posn > 1005 and posn < 1095) or
                             (posn > 2005 and posn < 3995) or
                             (posn > 5005 and posn < 6995) or
                             (posn > 8005 and posn < 8995) )

        self.gen._mutation_list = []
        posn_list = self.gen.find_empty_ranges(50, 4, 5)
        for posn in posn_list:
            self.assertTrue(posn > 0 and posn < 10000)

    def test_random_sequence(self):
        rand_seq = self.gen.random_sequence(10)
        self.assertEqual(len(rand_seq), 10)
        for allele in rand_seq:
            self.assertTrue(allele in ['A','C','G','T'])

    def test_delete_block(self):
        sequence = 'THIS IS A TEST SEQUENCE'
        sequence = self.gen.delete_block(sequence, 5, 3)
        self.assertEqual(sequence, 'THIS A TEST SEQUENCE')
        sequence = self.gen.delete_block(sequence, 0, 1)
        self.assertEqual(sequence, 'HIS A TEST SEQUENCE')
        sequence = self.gen.delete_block(sequence, -1, 1)
        self.assertEqual(sequence, 'HIS A TEST SEQUENC')
        sequence = self.gen.delete_block(sequence, -2, 1)
        self.assertEqual(sequence, 'HIS A TEST SEQUEC')
        sequence = self.gen.delete_block(sequence, -2, 3)
        self.assertEqual(sequence, 'HIS A TEST SEQU')

    def test_insert_block(self):
        sequence = 'HIS A TEST SEQUE'
        sequence = self.gen.insert_block(sequence, -1, 'EC')
        self.assertEqual(sequence, 'HIS A TEST SEQUECE')
        sequence = self.gen.insert_block(sequence, -2, 'N')
        self.assertEqual(sequence, 'HIS A TEST SEQUENCE')
        sequence = self.gen.insert_block(sequence, 0, 'T')
        self.assertEqual(sequence, 'THIS A TEST SEQUENCE')
        sequence = self.gen.insert_block(sequence, 5, 'IS ')
        self.assertEqual(sequence, 'THIS IS A TEST SEQUENCE')
        sequence = self.gen.insert_block(sequence, len(sequence), '!')
        self.assertEqual(sequence, 'THIS IS A TEST SEQUENCE!')

    def test_overwrite_block(self):
        sequence = 'THIS IS A TEST SEQUENCE!'
        sequence = self.gen.overwrite_block(sequence, 2, 'AT')
        self.assertEqual(sequence, 'THAT IS A TEST SEQUENCE!')
        sequence = self.gen.overwrite_block(sequence, 0, 'W')
        self.assertEqual(sequence, 'WHAT IS A TEST SEQUENCE!')
        sequence = self.gen.overwrite_block(sequence, -1, '?')
        self.assertEqual(sequence, 'WHAT IS A TEST SEQUENCE?')
        sequence = self.gen.overwrite_block(sequence, -1, '?!?')
        self.assertEqual(sequence, 'WHAT IS A TEST SEQUENCE?!?')
        sequence = self.gen.overwrite_block(sequence, len(sequence), '!')
        self.assertEqual(sequence, 'WHAT IS A TEST SEQUENCE?!?!')

    def test_invert_block(self):
        sequence = 'THIS IS A TEST SEQUENCE'
        sequence, orig_block, inverted_block = self.gen.invert_block(sequence, 0, 1)
        self.assertEqual(orig_block, 'T')
        self.assertEqual(inverted_block, 'T')
        self.assertEqual(sequence, 'THIS IS A TEST SEQUENCE')
        sequence, orig_block, inverted_block = self.gen.invert_block(sequence, 0, 2)
        self.assertEqual(orig_block, 'TH')
        self.assertEqual(inverted_block, 'HT')
        self.assertEqual(sequence, 'HTIS IS A TEST SEQUENCE')
        sequence, orig_block, inverted_block = self.gen.invert_block(sequence, -1, 1)
        self.assertEqual(orig_block, 'E')
        self.assertEqual(inverted_block, 'E')
        self.assertEqual(sequence, 'HTIS IS A TEST SEQUENCE')
        sequence, orig_block, inverted_block = self.gen.invert_block(sequence, -2, 2)
        self.assertEqual(orig_block, 'CE')
        self.assertEqual(inverted_block, 'EC')
        self.assertEqual(sequence, 'HTIS IS A TEST SEQUENEC')
        sequence, orig_block, inverted_block = self.gen.invert_block(sequence, 5, 4)
        self.assertEqual(orig_block, 'IS A')
        self.assertEqual(inverted_block, 'A SI')
        self.assertEqual(sequence, 'HTIS A SI TEST SEQUENEC')
        sequence, orig_block, inverted_block = self.gen.invert_block(sequence, len(sequence) - 2, 2)
        self.assertEqual(orig_block, 'EC')
        self.assertEqual(inverted_block, 'CE')
        self.assertEqual(sequence, 'HTIS A SI TEST SEQUENCE')

    def test_generate_snp_allele(self):
        snp_allele = self.gen.generate_snp_allele('A')
        self.assertTrue(snp_allele in ['T','G','C'])
        snp_allele = self.gen.generate_snp_allele('T')
        self.assertTrue(snp_allele in ['A','G','C'])
        snp_allele = self.gen.generate_snp_allele('G')
        self.assertTrue(snp_allele in ['T','A','C'])
        snp_allele = self.gen.generate_snp_allele('C')
        self.assertTrue(snp_allele in ['T','G','A'])

    def test_generate_str_base(self):
        str_base = self.gen.generate_str_base(2)
        self.assertTrue(str_base[0] != str_base[1])
        for allele in str_base:
            self.assertTrue(allele in ['T','A','C','G'])

        str_base = self.gen.generate_str_base(3)
        self.assertTrue((str_base[0] != str_base[1]) or
                        (str_base[1] != str_base[2]) or
                        (str_base[0] != str_base[2]) )
        for allele in str_base:
            self.assertTrue(allele in ['T','A','C','G'])

        str_base = self.gen.generate_str_base(4)
        self.assertEqual(len(str_base), 4)
        self.assertTrue(str_base[:2] != str_base[2:])
        for allele in str_base:
            self.assertTrue(allele in ['T','A','C','G'])

        str_base = self.gen.generate_str_base(5)
        self.assertEqual(len(str_base), 5)
        for allele in str_base:
            self.assertTrue(allele in ['T','A','C','G'])

    def test_generate_alu_sequence(self):
        alu_seq = self.gen.generate_alu_sequence(300)
        self.assertTrue(len(alu_seq) == 300)
        self.assertTrue(alu_seq != self.gen._base_alu)
        for allele in alu_seq:
            self.assertTrue(allele in ['T','A','C','G'])
        alu_seq = self.gen.generate_alu_sequence(295)
        self.assertTrue(len(alu_seq) == 295)
        self.assertTrue(alu_seq != self.gen._base_alu[:295])
        self.assertTrue(alu_seq != self.gen._base_alu[4:])
        for allele in alu_seq:
            self.assertTrue(allele in ['T','A','C','G'])
        alu_seq = self.gen.generate_alu_sequence(305)
        self.assertTrue(len(alu_seq) == 305)
        self.assertTrue(alu_seq != self.gen._base_alu)
        for allele in alu_seq:
            self.assertTrue(allele in ['T','A','C','G'])

    def test_ranged_length_list(self):
        self.assertRaises(Exception, self.gen.ranged_length_list, 10, 10, 0)
        self.assertRaises(Exception, self.gen.ranged_length_list, 1, 0, 1)

        length_list = self.gen.ranged_length_list(10, 10, 10)
        self.assertEqual(10, len(length_list))
        for item in length_list:
            self.assertEqual(item, 10)

        length_list = self.gen.ranged_length_list(10, 20, 10)
        self.assertEqual(10, len(length_list))
        max_count = 0
        for item in length_list:
            if item == 20:
                max_count += 1
            self.assertTrue(item <= 20 and item >= 10)
        self.assertEqual(max_count, 2)

        length_list = self.gen.ranged_length_list(10, 20, 1)
        self.assertEqual(1, len(length_list))
        self.assertTrue(length_list[0] == 20)

        length_list = self.gen.ranged_length_list(10, 20, 2)
        self.assertEqual(2, len(length_list))
        self.assertTrue(length_list[0] == 20 and length_list[1] == 20)

        length_list = self.gen.ranged_length_list(10, 20, 3)
        self.assertEqual(3, len(length_list))
        self.assertTrue(length_list[0] == 20 and length_list[1] == 20 and length_list[2] == 10)

    def test_write_genome_lines_to_file(self):
        length_list = [0,1,79,80,81]
        for i in length_list:
            with open('test_file', 'w') as test_file:
                self.gen.write_genome_lines_to_file(self.gen._base_alu, test_file)
            with open('test_file', 'r') as test_file:
                base_alu = ''
                for line in test_file:
                    base_alu += str(line).strip()
            self.assertEqual(base_alu, self.gen._base_alu)
            os.remove('test_file')

    def test_parse_fasta(self):
        length_list = [1,79,80,81,300,10000]
        for next_len in length_list:
            original_sequence = self.gen.random_sequence(next_len)
            file_name = 'test_file_' + str(1) + '_' + str(next_len)
            with open(file_name, 'w') as test_file:
                test_file.write('>test\n')
                self.gen.write_genome_lines_to_file(original_sequence, test_file)
            for sequence in self.gen.parse_fasta(file_name):
                if sequence:
                    self.assertEqual(sequence, original_sequence)
                else:
                    break
            os.remove(file_name)

    def test_create_read_pair(self):
        donor_genome = self.gen.random_sequence(200)
        left_read, right_read = self.gen.create_read_pair(donor_genome, 0, 0)
        self.assertEqual(left_read, right_read[::-1])
        self.assertTrue( (left_read in donor_genome and right_read[::-1] in donor_genome)
                                    or
                        (right_read in donor_genome and left_read[::-1] in donor_genome) )

        left_stt = random.randint(0, 200 - self.gen._sequencer_read_length - 1)
        right_stt = random.randint(0, 200 - self.gen._sequencer_read_length - 1)
        left_read, right_read = self.gen.create_read_pair(donor_genome, left_stt, right_stt)
        self.assertEqual(len(left_read), len(right_read))
        self.assertTrue( (left_read in donor_genome and right_read[::-1] in donor_genome)
                                    or
                        (right_read in donor_genome and left_read[::-1] in donor_genome) )

        left_stt = 0
        right_stt = 200 - self.gen._sequencer_read_length - 1
        left_read, right_read = self.gen.create_read_pair(donor_genome, left_stt, right_stt)
        self.assertEqual(len(left_read), len(right_read))
        self.assertTrue( (left_read in donor_genome and right_read[::-1] in donor_genome)
                                    or
                        (right_read in donor_genome and left_read[::-1] in donor_genome) )

    def test_add_sequencer_errors(self):
        error_found = False
        donor_genome = self.gen.random_sequence(200)
        left_read, right_read = self.gen.create_read_pair(donor_genome, 0, 0)
        self.assertEqual(left_read, right_read[::-1])
        self.assertTrue( (left_read in donor_genome and right_read[::-1] in donor_genome)
                                    or
                        (right_read in donor_genome and left_read[::-1] in donor_genome) )
        for i in range(25):
            left_read_w_error = self.gen.add_sequencer_errors(left_read)
            right_read_w_error = self.gen.add_sequencer_errors(right_read)
            if left_read != left_read_w_error or right_read != right_read_w_error:
                error_found = True
        self.assertTrue(error_found)

    def test_generate_ref_genome(self):
        for alu in ['y', 'n']:
            for test_args in [[10, 'k'], [100, 'k']]:
                args = TestSetting()
                args.id = 'test'
                args.chr_id = 1
                args.chr_size = test_args[0]
                args.scale = test_args[1]
                args.alu = alu
                args.assembly = 'n'
                args.spectrum = 'n'
                args.base_alu = random_sequence(300)

                self.gen = chromosome_builder(args)
                self.gen._alu_min_length = 300
                self.gen._alu_max_length = 300
                self.gen._alu_mutation_rate = 0.3

                self.gen.generate_ref_genome()
                ref_genome = next(self.gen.parse_fasta(self.gen._ref_genome_file))
                cnv_dict = {}
                cnv_count = 0
                str_count = 0
                for mutation in self.gen._mutation_list:
                    if mutation[2] == 'REF_STR':
                        str_count += 1
                        self.assertTrue(ref_genome[mutation[0]:mutation[1]] == mutation[3]*mutation[4])
                    elif mutation[2] == 'REF_CNV':
                        if mutation[3] in cnv_dict:
                            self.assertTrue(cnv_dict[mutation[3]] == ref_genome[mutation[0]:mutation[1]])
                        else:
                            cnv_count += 1
                            cnv_dict[mutation[3]] = ref_genome[mutation[0]:mutation[1]]
                    elif mutation[2] == 'REF_ALU':
                        base_alu = self.gen._base_alu
                        match_count = 0
                        for i in range(len(base_alu)):
                            if self.gen._base_alu[i] == ref_genome[mutation[0]+i]:
                                match_count += 1
                        self.assertTrue((match_count / float(len(base_alu))) > (.99 - self.gen._alu_mutation_rate))
                self.assertEqual(cnv_count, self.gen._nbr_ref_cnv)
                self.assertEqual(str_count, self.gen._nbr_ref_str)
        os.remove(self.gen._ref_genome_file)
        os.remove(self.gen._base_alu_file)

    def test_generate_donor_genome(self):
        args = TestSetting()
        args.id = 'test'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.base_alu = random_sequence(300)

        self.gen = chromosome_builder(args)
        self.gen._nbr_snp = 0
        self.gen._nbr_denovo_str = 0
        self.gen._nbr_denovo_cnv = 0
        self.gen._nbr_long_inv = 0
        self.gen._nbr_long_ins = 0
        self.gen._nbr_long_del = 0
        self.gen._nbr_ref_alu = 0
        self.gen._nbr_denovo_alu = 0

        self.gen._nbr_ref_str = 0
        self.gen._nbr_ref_cnv = 0
        self.gen._nbr_short_inv = 0
        self.gen._nbr_short_ins = 0
        self.gen._nbr_short_del = 0
        self.gen._cnv_mutation_amount = 0
        self.gen._str_mutation_amount = 0

        self.gen.generate_ref_genome()
        self.gen.generate_donor_genome()
        ref_genome = ''
        for sequence in self.gen.parse_fasta(self.gen._ref_genome_file):
            if sequence:
                ref_genome += sequence
            else:
                break
        donor_genome = ''
        for sequence in self.gen.parse_fasta(self.gen._priv_genome_file):
            if sequence:
                donor_genome += sequence
            else:
                break
        self.assertEqual(ref_genome, donor_genome)
        self.assertEqual(len(ref_genome), 10000)

        for alu in ['y', 'n']:
            for test_args in [[10, 'k'], [100, 'k'], [150, 'k']]:
                args = TestSetting()
                args.id = 'test'
                args.chr_id = 1
                args.chr_size = test_args[0]
                args.scale = test_args[1]
                args.alu = alu
                args.assembly = 'n'
                args.spectrum = 'n'
                args.base_alu = random_sequence(300)

                if args.scale == 'k':
                    expected_size = test_args[0] * 1000
                elif args.scale == 'm':
                    expected_size = test_args[0] * 1000000

                self.gen = chromosome_builder(args)
                self.gen._alu_min_length = 300
                self.gen._alu_max_length = 300
                self.gen._alu_mutation_rate = 0.3

                self.gen.generate_ref_genome()
                self.gen.generate_donor_genome()
                ref_genome = ''
                for sequence in self.gen.parse_fasta(self.gen._ref_genome_file):
                    if sequence:
                        ref_genome += sequence
                    else:
                        break
                donor_genome = ''
                for sequence in self.gen.parse_fasta(self.gen._priv_genome_file):
                    if sequence:
                        donor_genome += sequence
                    else:
                        break

                last_end = 0
                self.assertEqual(expected_size, len(ref_genome))

                for i in range(len(self.gen._mutation_list)):
                    mutation = self.gen._mutation_list[i]
                    self.assertTrue(ref_genome[last_end:mutation[0]] in donor_genome)
                    last_end = mutation[1]
                    mut_type = mutation[2]
                    range_stt = max(0, mutation[0]-20)
                    range_end = min(len(ref_genome)-1, mutation[0]+20)
                    gapped_range_end = min(len(ref_genome)-1, mutation[1]+20)

                    if mut_type == 'SNP':
                        self.assertTrue(ref_genome[range_stt:range_end] not in donor_genome, msg='SNP ' + str(mutation[0]))
                    elif mut_type == 'MUT_STR':
                        new_str = mutation[3]
                        self.assertTrue(new_str in donor_genome, msg='MUT_STR ' + str(mutation[0]))
                    elif mut_type == 'DONOR_STR':
                        str_seq = mutation[3]
                        nbr_copies = mutation[4]
                        new_str = str_seq * nbr_copies
                        self.assertTrue(new_str in donor_genome, msg='DONOR_STR ' + str(mutation[0]))
                    elif mut_type == 'REF_ALU':
                        self.assertTrue(ref_genome[mutation[0]:mutation[1]] in donor_genome, msg='REF_ALU ' + str(mutation[0]))
                    elif mut_type == 'REF_CNV':
                        self.assertTrue(ref_genome[mutation[0]:mutation[1]] in donor_genome, msg='REF_CNV ' + str(mutation[0]))
                    elif mut_type == 'DONOR_ALU':
                        self.assertTrue(ref_genome[range_stt:gapped_range_end] not in donor_genome, msg='DONOR_ALU ' + str(mutation[0]))
                    elif mut_type == 'INV':
                        inv_seq = ref_genome[mutation[0]:mutation[1]]
                        inv_seq = inv_seq[::-1]
                        self.assertTrue(inv_seq in donor_genome, msg='INV ' + str(mutation[0]))
                    elif mut_type == 'INS':
                        self.assertTrue(ref_genome[range_stt:range_end] not in donor_genome, msg='INS ' + str(mutation[0]))
                    elif mut_type == 'DEL':
                        self.assertTrue(ref_genome[range_stt:gapped_range_end] not in donor_genome, msg='DEL ' + str(mutation[0]))
        os.remove(self.gen._priv_genome_file)
        os.remove(self.gen._answer_file)
        os.remove(self.gen._ref_genome_file)
        os.remove(self.gen._base_alu_file)

    def test_generate_spectrum_reads(self):
        args = TestSetting()
        args.id = 'test'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'y'
        args.spectrum = 'y'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_spectrum_reads()

        with open(self.gen._reads_file, "r") as reads_file:
            reads_file.readline()
            with open(self.gen._priv_genome_file, "r") as donor_genome_file:
                donor_genome_file.readline()
                donor_genome = ""
                for line in donor_genome_file:
                    if ">" in line:
                        break
                    donor_genome += str(line).strip()

                i = 0
                for line in reads_file:
                    line = line.strip()
                    left_read = line.split(',')[0]
                    right_read = line.split(',')[1]
                    self.assertEqual(left_read, donor_genome[i:i+self.gen._sequencer_read_length])
                    right_start = i + self.gen._sequencer_read_length + 100
                    self.assertEqual(right_read, donor_genome[right_start:right_start+self.gen._sequencer_read_length])
                    i += 1
        os.remove(self.gen._priv_genome_file)
        os.remove(self.gen._answer_file)
        os.remove(self.gen._ref_genome_file)
        os.remove(self.gen._reads_file)

def random_sequence(seq_len):
    return "".join(random.choice(['A','C','G','T']) for i in range(seq_len))

class TestSetting():
    def __init__(self):
        self.id = None
        self.num_chr = None
        self.chr_size = None
        self.scale = None
        self.alu = None
        self.assembly = None
        self.simple = None
        self.spectrum = None
        self.base_alu = None
        self.verbose = None

class AlgInBioinf():
    def __init__(self):

        self._ans_key_folder = 'answer_keys'
        self._distribution_folder = 'for_distribution'

        if not os.path.exists(self._ans_key_folder):
            os.makedirs(self._ans_key_folder)
        if not os.path.exists(self._distribution_folder):
            os.makedirs(self._distribution_folder)

    """
    Generates data files and places them into appropriate zip files for distribution.
    """
    def hw0(self):
        args = TestSetting()
        args.id = 'hw0_W_0'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'y'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = random_sequence(300)
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        os.remove(gen._base_alu_file)

        student_zip_file_name = 'sample_HW0.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)

    """
    Generates data files and places them into appropriate zip files for distribution.
    """
    def hw1(self):
        #Set 1: Includes all data files for students to use when developing their algorithm
        args = TestSetting()
        args.id = 'practice_W_1'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'y'  #HW1 only uses SNPs, INDELs, and CNVs
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)
            zip_file.write(gen._reads_file)
            zip_file.write(gen._ref_genome_file)
            zip_file.write(gen._priv_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

        #Set 2: For credit, students given only the reads and ref genome
        args = TestSetting()
        args.id = 'hw1_W_2'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'y'  #HW1 only uses SNPs, INDELs, and CNVs
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._reads_file)
            zip_file.write(gen._ref_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

    """
    Generates data files and places them into appropriate zip files for distribution.
    """
    def hw2_undergrad(self):
        #Set 1: Includes all data files for students to use when developing their algorithm
        args = TestSetting()
        args.id = 'practice_W_3'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)
            zip_file.write(gen._reads_file)
            zip_file.write(gen._ref_genome_file)
            zip_file.write(gen._priv_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

        #Set 2: Includes all data files for students to use when developing their algorithm
        args = TestSetting()
        args.id = 'practice_E_1'
        args.chr_id = 1
        args.chr_size = 1
        args.scale = 'm'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)
            zip_file.write(gen._reads_file)
            zip_file.write(gen._ref_genome_file)
            zip_file.write(gen._priv_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

        #Set 3: For credit, students given only the reads and ref genome
        args = TestSetting()
        args.id = 'hw2undergrad_E_2'
        args.chr_id = 1
        args.chr_size = 1
        args.scale = 'm'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._reads_file)
            zip_file.write(gen._ref_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

    def hw2_grad(self):
                #Set 3: For credit, students given only the reads and reference genome
        args = TestSetting()
        args.id = 'hw2grad_M_1'
        args.chr_id = 1
        args.chr_size = 100
        args.scale = 'm'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._reads_file)
            zip_file.write(gen._ref_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._answer_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

    """
    Generates data files and places them into appropriate zip files for distribution.
    """
    def hw3(self):
        #Set 1: Includes all data files for students to use when developing their algorithm
        args = TestSetting()
        args.id = 'spectrum_A_1'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'y'
        args.spectrum = 'y'
        args.simple = 'n'
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_spectrum_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._reads_file)
            zip_file.write(gen._priv_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._priv_genome_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

        #Set 2: Includes all data files for students to use when developing their algorithm
        args = TestSetting()
        args.id = 'practice_A_2'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'y'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = ''
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._reads_file)
            zip_file.write(gen._priv_genome_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._priv_genome_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

        #Set 3: For credit, students given only the reads
        args = TestSetting()
        args.id = 'hw3all_A_3'
        args.chr_id = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'y'
        args.spectrum = 'n'
        args.simple = 'n'
        args.base_alu = 'n'
        args.verbose = 'n'

        gen = chromosome_builder(args)
        gen.generate_ref_genome()
        gen.generate_donor_genome()
        gen.generate_reads()

        student_zip_file_name = args.id + '.zip'
        student_zip_file_name = os.path.join(self._distribution_folder, student_zip_file_name)
        ans_key_zip_file_name = gen._answer_file.split('.')[0] + '.zip'
        ans_key_zip_file_name = os.path.join(self._ans_key_folder, ans_key_zip_file_name)

        with zipfile.ZipFile(student_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._reads_file)
        with zipfile.ZipFile(ans_key_zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            zip_file.write(gen._priv_genome_file)

        os.remove(gen._ref_genome_file)
        os.remove(gen._priv_genome_file)
        os.remove(gen._answer_file)
        os.remove(gen._reads_file)

if __name__ == '__main__':
    gen = AlgInBioinf()
    gen.hw0()
    gen.hw1()
    gen.hw2_undergrad()
    gen.hw2_grad()
    gen.hw3()
