__author__ = 'Joseph'

import unittest
import os

class helpful_scripts():

    def ref_check(self, answer_file, ref_genome_file):
        """
        Compares the indices used in the answer file with the reference genome to check whether the
        correct index values are being used.
        @params:
        answer_file: filename for the answer file that will be submitted for evaluation
        ref_genome_file: filename for the reference genome file provided with the reads
        @returns:
        a string indicating whether or not mismatches were found
        """
        genome = ""
        mismatches = 0
        with open(ref_genome_file, "r") as gen_file:
            for line in gen_file:
                if line.startswith(">"):
                    pass
                else:
                    line = line.rstrip()
                    genome += line
        with open(answer_file, "r") as ans_file:
            in_snp = False
            for line in ans_file:
                if line.startswith(">SNP"):
                    in_snp = True
                elif line.startswith(">"):
                    in_snp = False
                elif in_snp:
                    line = line.strip()
                    if line:
                        split_line = line.split(",")
                        allele = split_line[0]
                        posn = int(split_line[2])
                        if allele != genome[posn]:
                            mismatches += 1
        if mismatches == 0:
            return "No mismatches"
        else:
            return str(mismatches) + " mismatches found."

    def sort_file(self, file_name):
        unsorted = open(file_name, "r")
        lines = [line for line in unsorted if line.strip()]
        unsorted.close()
        if not lines[-1].endswith("\n"):
            lines[-1] = lines[-1] + "\n"
        lines[-1] = lines[-1] + '\n'
        with open(file_name, "w") as sorted:
            ans_lines = []
            current_type = None
            for file_line in lines:
                if file_line.startswith(">"):
                    if len(ans_lines) > 0:
                        ans_lines = [line for line in ans_lines if line.strip()]
                        if "SNP" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[2]))
                        elif "STR" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[1]))
                        elif "CNV" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[1]))
                        elif "ALU" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[1]))
                        elif "INV" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[1]))
                        elif "INS" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[1]))
                        elif "DEL" in current_type:
                            ans_lines.sort(key=lambda l: int(l.split(",")[1]))
                        sorted.writelines(ans_lines)
                        ans_lines = []
                    if len(file_line.strip()) > 1:
                        sorted.write(file_line)
                        current_type = file_line[1:].rstrip()
                else:
                    #sort all copy posns on a single line
                    if "CNV" in current_type:
                        seq = file_line.split(",")[0]
                        posn_list = file_line.split(",")[1:]
                        posn_list.sort(key=int)
                        new_line = seq
                        for posn in posn_list:
                            new_line += "," + str(posn.rstrip())
                        new_line += "\n"
                        ans_lines.append(new_line)
                    else:
                        ans_lines.append(file_line)

class TestClass(unittest.TestCase):

    def setUp(self):
        self.script = helpful_scripts()

        self.data_directory = 'test_data'

        self.sorted_file = 'ref_check_ans_file.txt'
        self.temp_sorted_file = 'temp_sorted_file.txt'
        self.unsorted_file = 'unsorted_answer_file.txt'

        self.reference = 'ref_check_reference_file.txt'
        self.ref_check_ans = 'ref_check_ans_file.txt'

    def test_ref_check(self):
        ref_filename = os.path.join(self.data_directory, self.reference)
        ans_filename = os.path.join(self.data_directory, self.ref_check_ans)
        result = self.script.ref_check(ans_filename, ref_filename)
        self.assertEqual(result, "No mismatches")

    def test_sort_file(self):
        unsorted_filename = os.path.join(self.data_directory, self.unsorted_file)
        temp_unsorted_filename = os.path.join(self.data_directory, self.temp_sorted_file)
        sorted_filename = os.path.join(self.data_directory, self.sorted_file)
        with open(unsorted_filename, 'r') as unsorted:
            with open(temp_unsorted_filename, 'w') as temp_unsorted:
                for line in unsorted:
                    temp_unsorted.write(line)

        self.script.sort_file(temp_unsorted_filename)

        with open(temp_unsorted_filename, 'r') as temp_sorted:
            with open(sorted_filename, 'r') as sorted:
                for line1 in temp_sorted:
                    line2 = sorted.readline()
                    self.assertEqual(line1.strip(), line2.strip())

        os.remove(temp_unsorted_filename)