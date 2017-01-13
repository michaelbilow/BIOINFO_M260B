"""
Created on June 20, 2014

@author: Joseph Korpela
"""

import time
import random
import sys
import re
import argparse
import pickle
import shutil
import unittest
import os
import chromosome_builder
import reads_mixer

class genome_builder():
    def __init__(self, args=None):
        if not args:
            args = self.parse_system_args()

        self._genome_id = args.id
        self._nbr_chromosome = args.nbr_chr
        self._chromosome_size = args.chr_size
        self._scale = args.scale
        self._use_alu = args.alu
        self._use_assembly = args.assembly
        self._use_spectrum = args.spectrum
        self._verbose = False

        self._allele_base_list = ["C", "T", "G", "A"]

        self._alu_min_length = 300
        self._alu_max_length = 300

        self._reads_file = "reads_" + str(self._genome_id) + ".txt"
        self._chromosome_reads_file = "reads_" + str(self._genome_id) + "_chr_"
        self._base_alu_file = "alu_" + str(self._genome_id) + ".txt"

        if self._use_alu:
            self._base_alu = self.generate_base_alu()

    def insert_newlines(self, sequence, line_size=80):
        return '\n'.join(sequence[i:i+line_size] for i in range(0, len(sequence), line_size))

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
            start_of_file = True
            buffer = ""
            while True:
                for line in fasta_file:
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
                if len(buffer) > 0:
                    yield buffer
                    buffer = ''
                else:
                    yield None

    def random_sequence(self, seq_len):
        return "".join(random.choice(self._allele_base_list) for i in range(seq_len))

    def generate_base_alu(self):
        alu_len = random.randint(self._alu_min_length, self._alu_max_length)
        self._base_alu = self.random_sequence(alu_len)

    def worker(self, chromosome):
        args = TestSetting()
        args.id = self._genome_id
        args.chr_id = chromosome
        args.chr_size = self._chromosome_size
        args.scale = self._scale
        args.alu = self._use_alu
        args.assembly = self._use_assembly
        args.spectrum = self._use_spectrum
        args.base_alu = self.random_sequence(300)

        gen = chromosome_builder(args)
        if self._verbose:
            print('chromosome ' + str(chromosome) + ' generating ref genome')
        gen.generate_ref_genome()
        if self._verbose:
            print('chromosome ' + str(chromosome) + ' generating donor genome')
        gen.generate_donor_genome()
        if self._verbose:
            print('chromosome ' + str(chromosome) + ' generating reads')
        if self._use_spectrum == 'y':
            gen.generate_spectrum_reads()
        else:
            gen.generate_reads()

    def generate_chromosome(self):
        """
        Generates a random reference genome with the specified number of chromosomes,
        each of length length_chromosome
        """
        if self._use_alu:
            with open(self._base_alu_file, "w") as alu_file:
                self.write_genome_lines_to_file(self._base_alu, alu_file)
        jobs = []
        for chromosome in range(1, self._num_chromosomes + 1):
            chr_builder = multiprocessing.Process(target=self.worker, name='worker')
            jobs.append(chr_builder)
            chr_builder.start()
        for job in jobs:
            job.join()

        mixer = reads_mixer(self._genome_id)

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
            "--num_chr",
            type=int,
            default='1',
            help="The number of chromosomes to generate for the genome."
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
            help="whether to generate output for assembly as the spectrum (from position 1 up in order)."
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
        self.gen = genome_builder()

    def test_generate_base_alu(self):
        self.gen._alu_min_length = 300
        self.gen._alu_max_length = 301
        self.gen._base_alu = ''
        self.gen.generate_base_alu()
        self.assertTrue(len(self.gen._base_alu) >= self.gen._alu_min_length)
        self.assertTrue(len(self.gen._base_alu) <= self.gen._alu_max_length)
        for allele in self.gen._base_alu:
            self.assertTrue(allele in ['T','A','C','G'])

    def test_write_genome_lines_to_file(self):
        length_list = [0,1,79,80,81]
        for i in length_list:
            self.gen._alu_min_length = i
            self.gen._alu_max_length = i
            self.gen.generate_base_alu()

            with open('test_file', 'w') as test_file:
                self.gen.write_genome_lines_to_file(self.gen._base_alu, test_file)
            with open('test_file', 'r') as test_file:
                base_alu = ''
                for line in test_file:
                    base_alu += str(line).strip()
            self.assertEqual(base_alu, self.gen._base_alu)
            os.remove('test_file')

    def test_parse_fasta(self):
        nbr_chr_list = [1,2,3]
        length_list = [1,79,80,81]
        for nbr_chr in nbr_chr_list:
            for next_len in length_list:
                self.gen._alu_min_length = next_len
                self.gen._alu_max_length = next_len
                self.gen.generate_base_alu()
                file_name = 'test_file_' + str(nbr_chr) + '_' + str(next_len)
                with open(file_name, 'w') as test_file:
                    test_file.write('>test')
                    for chr in range(1, nbr_chr + 1):
                        test_file.write('\n>chr' + str(chr) + '\n')
                        self.gen.write_genome_lines_to_file(self.gen._base_alu, test_file)
                for sequence in self.gen.parse_fasta(file_name):
                    if sequence:
                        base_alu = sequence
                        self.assertEqual(base_alu, self.gen._base_alu)
                    else:
                        break
                os.remove(file_name)


    def test_generate_chromosome(self):
        args = TestSetting()
        args.id = 'test'
        args.num_chr = 1
        args.chr_size = 10
        args.scale = 'k'
        args.alu = 'n'
        args.assembly = 'n'
        args.spectrum = 'n'

        self.gen = genome_gen(args)
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

        for alu in ['y', 'n']:
            for test_args in [[10, 'k'], [100, 'k'], [150, 'k']]:
                args = TestSetting()
                args.id = 'test'
                args.num_chr = 1
                args.chr_size = test_args[0]
                args.scale = test_args[1]
                args.alu = alu
                args.assembly = 'n'
                args.spectrum = 'n'

                self.gen = genome_gen(args)
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

class TestSetting():
    def __init__(self):
        self.id = None
        self.num_chr = None
        self.chr_size = None
        self.scale = None
        self.alu = None
        self.assembly = None
        self.spectrum = None
        self.base_alu = None
        self.verbose = None

if __name__ == '__main__':
    unittest.main()
    test_results = []
    for alu in ['y', 'n']:
        for test in [[100, 'k'], [200, 'k'], [300, 'k'], [400, 'k']]:#, [500, 'k'], [600, 'k'], [700, 'k'], [800, 'k'], [900, 'k'], [1, 'm']]:
            args = TestSetting()
            args.id = 'test'
            args.num_chr = 1
            args.chr_size = test[0]
            args.scale = test[1]
            args.alu = alu
            args.assembly = 'n'
            args.spectrum = 'n'

            start = time.clock()
            gen = genome_gen(args)
            print('generating ref genome')
            gen.generate_ref_genome()
            print('generating donor genome')
            gen.generate_donor_genome()
            print('generating reads')
            gen.generate_reads()
            test_results.append('Test: ' + str(test[0]) + test[1] + ' time: ' + str(time.clock() - start))
    for res in test_results:
        print(res)