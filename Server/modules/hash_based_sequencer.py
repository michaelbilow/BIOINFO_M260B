__author__ = 'Joseph'

import unittest
import Eval
import os
import cPickle as pickle

class MultiReSequencer():
    def __init__(self, allowable_error, key_size, ref_file, read_file, ans_file):
        self._snp_dict = {} #dict with key=posn, value=list of alleles for that posn
        self._var_dict = {} #dict with key=type of var
        self._var_dict['STR'] = []
        self._allowable_error = allowable_error #max number of alleles that can differ from reference
        self._genome_id = ''
        self._index = {}
        self._key_size = key_size
        self._read_file = read_file #file object with the reference sequence
        self._ans_file = ans_file #file object which will hold the answers
        self._ref_genome = self.read_ref_file(ref_file)
        self.build_index(store_index=False)
        self._mapped_reads = {}
        self._min_str = 3
        self._max_str = 5

    def read_ref_file(self, ref_file):
        genome = ''
        for line in ref_file:
            if line.startswith('>'):
                if self._genome_id == '':
                    self._genome_id = line
            else:
                line = line.rstrip()
                genome += line
        return genome

    def build_index(self, store_index):
        file_name = self._genome_id[1:].rstrip() + '.p'
        if store_index and os.path.isfile(file_name):
            self._index = pickle.load(open(file_name, 'rb'))
        else:
            for i in range(len(self._ref_genome)):
                if i + self._key_size <=  len(self._ref_genome):
                    next_key = self._ref_genome[i:i+self._key_size]
                    self.add_posn(next_key, i)
            if store_index:
                pickle.dump(self._index, open(file_name, 'wb'))

    def add_posn(self, key, posn):
        if self._index.has_key(key):
            self._index[key].add(posn)
        else:
            self._index[key] = set()
            self._index[key].add(posn)

    def add_snp(self, posn, allele):
        if self._snp_dict.has_key(posn):
            self._snp_dict[posn].append(allele)
        else:
            self._snp_dict[posn] = []
            self._snp_dict[posn].append(allele)

    def add_var(self, type, stt, stop):
        if self._var_dict.has_key(type):
            for i in range(stt, stop + 1):
                self._var_dict[type].append(i)
        else:
            self._var_dict[type] = []
            for i in range(stt, stop + 1):
                self._var_dict[type].append(i)

    def str_check(self, sequence, min_str, max_str):
        while len(sequence) / max_str < 3:
            if max_str > min_str:
                max_str -= 1
            else:
                return None #since this sequence is too short to id an str
        for i in range(len(sequence)):
            for str_size in reversed(range(min_str, max_str+1)):
                if i + (3 * str_size) < len(sequence):
                    if sequence[i:i+str_size] == sequence[i+str_size:i+2*str_size] and \
                        sequence[i+str_size:i+2*str_size] == sequence[i+2*str_size:i+3*str_size]:
                        return (sequence[i:i+str_size], i)
        return None

    def str_expand(self, str_seq, stt, sequence_to_search=None):
        if not sequence_to_search:
            sequence_to_search = self._ref_genome
        new_stt = stt
        repeats = 1
        #back up to find real start
        while new_stt == stt and (stt - len(str_seq)) >= 0:
            new_stt = stt - len(str_seq)
            if sequence_to_search[new_stt:stt] == str_seq:
                stt = new_stt
                repeats += 1
        new_repeats = repeats
        while new_repeats == repeats and len(sequence_to_search) > (stt + len(str_seq) * repeats + 1):
            new_repeats += 1
            if sequence_to_search[stt+len(str_seq)*repeats:stt+len(str_seq)*new_repeats] == str_seq:
                repeats = new_repeats
        return (str_seq, repeats, stt)

    def ref_genome_map_str(self, min_repeats):
        last_end = -1
        for i in range(len(self._ref_genome)):
            if i >= last_end:
                if i+4*self._max_str < len(self._ref_genome):
                    next_str = self.str_check(self._ref_genome[i:i+4*self._max_str], self._min_str, self._max_str)
                    if next_str:
                        next_str = self.str_expand(next_str[0], i+next_str[1])
                        if next_str[1] > 5:
                            self._var_dict['STR'].append((next_str))
                        last_end = next_str[1] * len(next_str[0]) + next_str[2]

    def process_reads(self):
        for line in self._read_file:
            if line.startswith('>'):
                if self._genome_id == '':
                    self._genome_id = line
            else:
                line = line.rstrip()
                for one_read in line.split(','):
                    self.process_one_read(one_read)

    def process_one_read(self, one_read):
        match_sets = []
        max_count = 0
        best_start = -1
        inverted = False
        number_splits = len(one_read) / self._key_size
        for i in range(number_splits):
            stt = i * self._key_size
            stop = stt + self._key_size
            this_set = self._index.get(one_read[stt:stop])
            if this_set:
                match_sets.append(this_set)
            else:
                match_sets.append(set())
        for invert in [False, True]:
            if invert:
                sequence = one_read[::-1]
            else:
                sequence = one_read
            for i in range(number_splits): #using each split of the read as the starting posn
                for j in match_sets[i]: #using each posn found for that split
                    this_count = 0
                    this_start = -1
                    for k in range(number_splits): #check how well all the splits match around that split
                        diff = (i - k) * self._key_size #how to determine posn adjustment for other set
                        if k == 0: #used when mapping snps to align start of one_read with ref
                            this_start = j - diff
                        if ((j - diff) in match_sets[k]):
                            this_count += 1
                    #don't use this count if its not better than the previous max, also don't use if all parts match
                    if (this_count > max_count) and (this_count < number_splits):
                        best_start = this_start
                        max_count = this_count
                        if invert:
                            inverted = True
        if max_count > max(2, number_splits - self._allowable_error):
            if inverted:
                self.map_snps(best_start, one_read[::-1])
            else:
                self.map_snps(best_start, one_read)
            return True
        else:
            return False

    def map_snps(self, posn, read_seq):
        mismatch_count = 0
        last_was_mismatch = False
        snp_list = []
        var_list = []
        for i in range(len(read_seq)):
            if self._ref_genome[i+posn] != read_seq[i]:
                mismatch_count += 1
                if last_was_mismatch:
                    var_list.append(i-1)
                else:
                    last_was_mismatch = True
            else:
                if last_was_mismatch:
                    snp_list.append(i-1)
                last_was_mismatch = False
        if last_was_mismatch and (len(read_seq)-2) not in snp_list:
            snp_list.append(len(read_seq)-1)
        if len(snp_list) < 3:
            for i in snp_list:
                self.add_snp((i+posn),read_seq[i])

    def clean_snps(self, min_snps):
        for key in self._snp_dict:
            if len(self._snp_dict[key]) > 0:
                SNP = ''
                a_count = 0
                c_count = 0
                g_count = 0
                t_count = 0
                for snp in self._snp_dict[key]:
                    if snp == 'A':
                        a_count += 1
                    elif snp == 'C':
                        c_count += 1
                    elif snp == 'G':
                        g_count += 1
                    elif snp == 'T':
                        t_count += 1
                if max(a_count, c_count, g_count, t_count) > min_snps:
                    if a_count == max(a_count, c_count, g_count, t_count):
                        SNP = 'A'
                    elif c_count == max(a_count, c_count, g_count, t_count):
                        SNP = 'C'
                    elif g_count == max(a_count, c_count, g_count, t_count):
                        SNP = 'G'
                    elif t_count == max(a_count, c_count, g_count, t_count):
                        SNP = 'T'
                    self._snp_dict[key] = SNP
                else:
                    self._snp_dict[key] = None
            else:
                self._snp_dict[key] = None

    def create_answer_file(self):
        self.clean_snps(3)
        self._ans_file.write(self._genome_id)
        self._ans_file.write('>SNP')
        for key in self._snp_dict:
            if self._snp_dict[key] != None:
                self._ans_file.write('\n' + self._ref_genome[key] + ',' + self._snp_dict[key] + ',' + str(key))
        if len(self._var_dict['STR']) > 0:
            self._ans_file.write('\n>STR')
            for seq, rpt, posn in self._var_dict['STR']:
                self._ans_file.write('\n' + (seq*rpt) + ',' + str(posn))

def main(): #, ref_file_name, read_file_name, answer_file_name, allowable_errors):
    ref_file_path = 'ref_genomeWALU.txt'
    reads_file_path = 'reads_genomeWALU.txt'
    ans_file_path = 'naive_ans_genomeWALU.txt'
    ans_key_path = 'ans_genomeWALU.txt'
    #with open(ref_file_path, 'r') as ref_file:
    #    with open(reads_file_path, 'r') as read_file:
    #        with open(ans_file_path, 'w') as answer_file:
    #            start = time.time()
    #            naive = MultiReSequencer(2, 10, ref_file, read_file, answer_file)
    #            naive.ref_genome_map_str(6)
    #            naive.process_reads()
    #            naive.create_answer_file()
    #            print 'Sequencing time: ' + str(time.time() - start)
    #sort_file(ans_file_path)
    #ref_check(ans_file_path, ref_file_path)
    with open(ans_file_path, 'r') as answer_file:
        with open(ans_key_path, 'r') as answer_key:
            ans_dict = Eval.Eval(answer_key, answer_file)
            for key in ans_dict:
                print key + ' ' + str(ans_dict[key])

if __name__ == '__main__':
    main()