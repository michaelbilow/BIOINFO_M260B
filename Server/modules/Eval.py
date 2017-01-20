"""
Created on Apr 22, 2014

@author: Jorge Munoz
@author: Joseph Korpela
@author: Kevin Takeshita
@author: James Gomez
@author: Gabriel Alsheikh
"""

# ########################
#
# To run: python Eval.py studentAnswers.txt genomeX
#
# #####################
import unittest
import os

class Eval():
    def __init__(self):
        self._genome_idx = {}
        self._key_size = 50
        self._key_idx = {}
        self._valid_types = {}

    #Using F1-Score
    def grade(self, studTot, corr, tot):
        if studTot == 0 or corr == 0:
            return 0
        p= float(corr)/max(studTot,tot)
        r= float(corr)/tot
        return 2*p*r/(p+r)

    def create_key_idx(self, ans_key):
        for idxType in [">CNV",">INV",">INS",">DEL",">SNP",">STR",">ALU",">ASSEMBLY" ]:
            self._key_idx[idxType] = -1
            self._valid_types[idxType] = False

        for i in range(0,len(ans_key)):
            for idxType in self._key_idx:
                if ( ans_key[i][0:len(idxType)]==idxType):
                    self._key_idx[idxType] = i
                    self._valid_types[idxType] = True

        for idxType in self._valid_types: #check that the answer key actually had entries for each type
            for idxType2 in self._key_idx:
                if idxType != idxType2:
                    if self._key_idx[idxType2] - 1 == self._key_idx[idxType]:
                        self._valid_types[idxType] = False

        #assembly answer keys are currently just saved as the donor genome file, so there will be no idxType >ASSEMBLY for now
        validTypeFound = False
        for idxType in self._valid_types:
            if self._valid_types[idxType]:
                validTypeFound = True
        if not validTypeFound:
            self._valid_types[">ASSEMBLY"] = True

    #def findIndex(self, arr, temp):
    #    for i in range(0,len(arr)):
    #        if ( arr[i][0:len(temp)]==temp):
    #            return i

    def hash_match(self, n, m):
        key_length = 10
        if len(n) < (3*key_length) or len(m) < (3*key_length):
            return self.needleman_wunsch(n, m)
        else:
            if len(n) > len(m):
                base_seq = n
                srch_seq = m
            else:
                base_seq = m
                srch_seq = n
            idx = set()
            for i in range(len(base_seq)-key_length+1):
                idx.add(base_seq[i:i+key_length])
            end = len(srch_seq) - key_length

            i = 0
            matched = 0
            while i <= end:
                if srch_seq[i:i+key_length] in idx:
                    i += 10
                    matched += 10
                else:
                    i += 1

            if i != len(srch_seq):
                if srch_seq[-key_length:] in idx:
                    matched += len(srch_seq) - i
            return matched

    def needleman_wunsch(self, n, m):
        rows = len(n) + 1
        cols = len(m) + 1
        matrix = []
        for row in range(rows):
            matrix.append(list())
            for col in range(cols):
                matrix[row].append(0)
        #initialize first column
        for row in range(rows):
          matrix[row][0] = -row
        #initalize first row
        for col in range(cols):
          matrix[0][col] = -col
        for row in range(1, rows):
            for col in range(1, cols):
                top = matrix[row-1][col]
                topleft = matrix[row-1][col-1]
                left = matrix[row][col-1]
                if n[row-1] == m[col-1]:
                    charscore = 1
                else:
                    charscore = -1
                max_val = max(top-1, topleft+charscore, left-1)
                matrix[row][col] = max_val
        return matrix[rows-1][cols-1]

    #criteria to grade the copy portion
    #stud contains the student answers
    #key contains the answer key
    # i is the index of the array for which the answers begin on
    # These structs will determine the number of false positives and negatives and
    # use the grade struct to determine an accurate ranking/grade
    #
    # TEST AGAIN FOR MULTIPLE COPIES
    #
    #SNP score function
    def SNPgrade(self, stud, key, stud_index):
        #used to index each line of the answer files
        orig_allele_idx = 0
        snp_idx = 1
        posn_idx = 2

        keyIndex = self._key_idx[">SNP"]
        keyIndex += 1 # index of first answer
        # find end of answers in key
        i = keyIndex
        while (i < len(key) and key[i][0] != '>'):
            i += 1
        key_answers = []
        for ans in key[keyIndex:i]:
            split_ans = ans.split(',')
            key_answers.append(split_ans)
        answer_key_length = len(key_answers)
        i = stud_index
        # find end of student answers
        while (i < len(stud) and stud[i][0]!='>'):
            i += 1
        new_i = i

        #if the student answer file has 400% the number of SNPs as really exists, just assign 0 and don't waste time scoring
        if (new_i - stud_index) > 4 * answer_key_length:
            return new_i, 0
        stud[stud_index:new_i] = sorted(stud[stud_index:new_i], key=lambda l: int(l.split(',')[posn_idx]))

        score = 0
        max_posn_diff = 5
        for student_ans in stud[stud_index:new_i]:
            student_ans = student_ans.split(',')
            remove_list = []

            #loop through remaining answer key entries to find the best match, if any
            for key_ans in key_answers:
                #if this key has already been passed, then mark it for removal
                if int(student_ans[posn_idx]) > int(key_ans[posn_idx])+max_posn_diff:
                    remove_list.append(key_ans)

                elif int(student_ans[posn_idx]) >= int(key_ans[posn_idx])-max_posn_diff and \
                   int(student_ans[posn_idx]) <= int(key_ans[posn_idx])+max_posn_diff and \
                        student_ans[snp_idx] == key_ans[snp_idx]:
                    score += 1
                    remove_list.append(key_ans)
                    break

                #once answers keys that are out of range are reached, then break the for loop
                elif int(student_ans[posn_idx]) < int(key_ans[posn_idx])-max_posn_diff:
                    break

            #get rid of the unneeded answer key entries
            for old_key in remove_list:
                key_answers.remove(old_key)

        return new_i, self.grade(len(stud[stud_index:new_i]), score, answer_key_length)

    def STRgrade(self, stud, key, stud_index):
        #used to index into entries in the answer files
        str_idx = 0
        posn_idx = 1

        keyIndex = self._key_idx[">STR"]
        keyIndex += 1 # index of first answer
        # find end of answers in key
        i = keyIndex
        while (i < len(key) and key[i][0] != '>'):
            i += 1
        key_answers = []
        for ans in key[keyIndex:i]:
            split_ans = ans.split(',')
            key_answers.append(split_ans)
        answer_key_length = len(key_answers)
        i = stud_index
        # find end of student answers
        while (i < len(stud) and stud[i][0]!='>'):
            i += 1
        new_i = i

        stud[stud_index:new_i] = sorted(stud[stud_index:new_i], key=lambda l: int(l.split(',')[posn_idx]))

        score = 0
        max_posn_diff = 20
        for student_ans in stud[stud_index:new_i]:
            student_ans = student_ans.split(',')
            remove_list = []

            #loop through remaining answer key entries to find the best match, if any
            for key_ans in key_answers:
                #if this key has already been passed, then mark it for removal
                if int(student_ans[posn_idx]) > int(key_ans[posn_idx])+max_posn_diff:
                    remove_list.append(key_ans)

                elif int(student_ans[posn_idx]) >= int(key_ans[posn_idx])-max_posn_diff and \
                   int(student_ans[posn_idx]) <= int(key_ans[posn_idx])+max_posn_diff:
                    score += 0.5
                    ans_sequence = key_ans[str_idx]
                    student_sequence = student_ans[str_idx]
                    max_align_score = max(len(ans_sequence), len(student_sequence)) #either sequence can be longer
                    #align_score = self.needleman_wunsch(ans_sequence, student_sequence)
                    align_score = self.hash_match(ans_sequence, student_sequence)
                    if align_score < 0:
                      align_score = 0
                    adj_score = (float(align_score)/(max_align_score)) #normalize to 0 to 1
                    adj_score = (adj_score**4) #since repeats only vary -2 to 2, its fairly easy to get close, so increasing penalty for error
                    score += adj_score / 2
                    remove_list.append(key_ans)
                    break

                #once answers keys that are out of range are reached, then break the for loop
                elif int(student_ans[posn_idx]) < int(key_ans[posn_idx])-max_posn_diff:
                    break

            #get rid of the unneeded answer key entries
            for old_key in remove_list:
                key_answers.remove(old_key)

        return new_i, self.grade(len(stud[stud_index:new_i]), score, answer_key_length)

    #criteria to grade the inversion portion
    def INVgrade (self, stud, key, stud_index):
        #used to index into entries in the answer files
        inv_idx = 0
        posn_idx = 1

        keyIndex=self._key_idx[">INV"]
        keyIndex += 1 # index of first answer
        # find end of answers in key
        i = keyIndex
        while (i < len(key) and key[i][0] != '>'):
            i += 1
        key_answers = []
        for ans in key[keyIndex:i]:
            split_ans = ans.split(',')
            key_answers.append(split_ans)
        answer_key_length = len(key_answers)
        i = stud_index
        # find end of student answers
        while (i < len(stud) and stud[i][0]!='>'):
            i += 1
        new_i = i

        stud[stud_index:new_i] = sorted(stud[stud_index:new_i], key=lambda l: int(l.split(',')[posn_idx]))

        score = 0
        max_posn_diff = 5
        for student_ans in stud[stud_index:new_i]:
            student_ans = student_ans.split(',')
            remove_list = []

            #loop through remaining answer key entries to find the best match, if any
            for key_ans in key_answers:
                #if this key has already been passed, then mark it for removal
                if int(student_ans[posn_idx]) > int(key_ans[posn_idx])+max_posn_diff:
                    remove_list.append(key_ans)

                elif int(student_ans[posn_idx]) >= int(key_ans[posn_idx])-max_posn_diff and \
                   int(student_ans[posn_idx]) <= int(key_ans[posn_idx])+max_posn_diff:
                    score += 0.5
                    ans_sequence = key_ans[inv_idx]
                    student_sequence = student_ans[inv_idx]
                    max_align_score = max(len(ans_sequence), len(student_sequence)) #either sequence can be longer
                    #align_score = self.needleman_wunsch(ans_sequence, student_sequence)
                    align_score = self.hash_match(ans_sequence, student_sequence)
                    if align_score < 0:
                        align_score = 0
                    adj_score = (float(align_score)/(max_align_score)) #normalize to 0 to 1
                    score += adj_score / 2.0
                    remove_list.append(key_ans)
                    break

                #once answers keys that are out of range are reached, then break the for loop
                elif int(student_ans[posn_idx]) < int(key_ans[posn_idx])-max_posn_diff:
                    break

            #get rid of the unneeded answer key entries
            for old_key in remove_list:
                key_answers.remove(old_key)

        return new_i, self.grade(len(stud[stud_index:new_i]), score, answer_key_length)

    #criteria to grade the insert portion
    def INDELgrade (self, stud, key, stud_index, insert_or_delete):
        #used to index into entries in the answer files
        indel_idx = 0
        posn_idx = 1

        keyIndex=self._key_idx[insert_or_delete]
        keyIndex += 1 # index of first answer
        # find end of answers in key
        i = keyIndex
        while (i < len(key) and key[i][0] != '>'):
            i += 1
        key_answers = []
        for ans in key[keyIndex:i]:
            split_ans = ans.split(',')
            key_answers.append(split_ans)
        answer_key_length = len(key_answers)
        i = stud_index
        # find end of student answers
        while (i < len(stud) and stud[i][0]!='>'):
            i += 1
        new_i = i

        stud[stud_index:new_i] = sorted(stud[stud_index:new_i], key=lambda l: int(l.split(',')[posn_idx]))

        score = 0
        max_posn_diff = 5
        for student_ans in stud[stud_index:new_i]:
            student_ans = student_ans.split(',')
            remove_list = []

            #loop through remaining answer key entries to find the best match, if any
            for key_ans in key_answers:
                #if this key has already been passed, then mark it for removal
                if int(student_ans[posn_idx]) > int(key_ans[posn_idx])+max_posn_diff:
                    remove_list.append(key_ans)

                elif int(student_ans[posn_idx]) >= int(key_ans[posn_idx])-max_posn_diff and \
                   int(student_ans[posn_idx]) <= int(key_ans[posn_idx])+max_posn_diff:
                    score += 0.5
                    ans_sequence = key_ans[indel_idx]
                    student_sequence = student_ans[indel_idx]
                    max_align_score = max(len(ans_sequence), len(student_sequence)) #either sequence can be longer
                    #align_score = self.needleman_wunsch(ans_sequence, student_sequence)
                    align_score = self.hash_match(ans_sequence, student_sequence)
                    if align_score < 0:
                        align_score = 0
                    adj_score = (float(align_score)/(max_align_score)) #normalize to 0 to 1
                    score += adj_score / 2.0
                    remove_list.append(key_ans)
                    break

                #once answers keys that are out of range are reached, then break the for loop
                elif int(student_ans[posn_idx]) < int(key_ans[posn_idx])-max_posn_diff:
                    break

            #get rid of the unneeded answer key entries
            for old_key in remove_list:
                key_answers.remove(old_key)

        return new_i, self.grade(len(stud[stud_index:new_i]), score, answer_key_length)

    def COPYgrade (self, stud, key, index):
        #used to index into entries in the answer files
        str_idx = 0
        posn_idx = 1
        max_posn_diff = 5

        ans_key_stt_line = self._key_idx[">CNV"]+1
        tmp_ans_key_line_idx = ans_key_stt_line
        correct=0
        total=0
        copynums=0
        m = index
        keynums=0
        studTot=0

        ans_key_posn_list = []
        ans_cnv_id = 0
        ans_cnv_list = []
        while (tmp_ans_key_line_idx < len(key) and key[tmp_ans_key_line_idx][0]!='>'):
            this_line = key[tmp_ans_key_line_idx].split(',')
            for i in range(1,len(this_line)):
                total += 1
                ans_key_posn_list.append((int(this_line[i]), ans_cnv_id))
            ans_cnv_list.append(tmp_ans_key_line_idx)
            ans_cnv_id += 1
            keynums+=1
            tmp_ans_key_line_idx += 1

        #create a 2x2 dictionary that stores the match score for each CNV entry in the answer key and student key
        #each time a pair is matched, check for the stored score first before running hash_match on the sequences
        #match_scores = {}
        student_posn_list = []
        student_cnv_id = 0
        student_cnv_list = []
        while (m < len(stud) and stud[m][0]!='>'):
            this_line = stud[m].split(',')
            for i in range(1,len(this_line)):
                studTot += 1
                student_posn_list.append((int(this_line[i]), student_cnv_id))
            student_cnv_list.append(m)
            #match_scores[student_cnv_id] = {}
            #for id in range(len(ans_cnv_list)):
            #    match_scores[student_cnv_id][id] = -1
            student_cnv_id += 1
            copynums+=1
            m += 1
        new_i = m
        ans_key_posn_list.sort(key=lambda t: t[0])
        student_posn_list.sort(key=lambda t: t[0])

        correct = 0
        student_ans_idx = 0
        ans_key_idx = 0
        while student_ans_idx < len(student_posn_list):
            #if we have reached the end of the ans key, then no more entries need to be checked
            if ans_key_idx == len(ans_key_posn_list):
                break
            #if the student's position is too low, then we can skip it
            elif student_posn_list[student_ans_idx][0] < ans_key_posn_list[ans_key_idx][0] - max_posn_diff:
                student_ans_idx += 1
            #if the student's position is too high, then we can skip this answer key entry
            elif student_posn_list[student_ans_idx][0] > ans_key_posn_list[ans_key_idx][0] + max_posn_diff:
                ans_key_idx += 1
            #else they match, in which case we should score the two entries
            else:
                student_posn = student_cnv_list[student_posn_list[student_ans_idx][1]]
                student_sequence = stud[student_posn].split(',')[0]
                ans_posn = ans_cnv_list[ans_key_posn_list[ans_key_idx][1]]
                ans_sequence = key[ans_posn].split(',')[0]
                #if this score was not previously computed, then compute and store it now
                #if match_scores[curr_stud_cnv_id][curr_ans_cnv_id] == -1:
                this_score = 0.5
                max_align_score = max(len(ans_sequence), len(student_sequence)) #either sequence can be longer
                align_score = 1.0
                align_score = self.hash_match(ans_sequence, student_sequence)
                if align_score < 0:
                    align_score = 0
                adj_score = float(align_score)/max_align_score #normalize to 0 to 1
                this_score += adj_score / 2.0
                #match_scores[curr_stud_cnv_id][curr_ans_cnv_id] = this_score
                correct += this_score
                #else:
                #    correct += match_scores[curr_stud_cnv_id][curr_ans_cnv_id]
                student_ans_idx += 1
                ans_key_idx += 1
        return new_i, self.grade(studTot, correct, total)

    def longest_increasing_subsequence(self, sequence):
        seq_len = len(sequence)
        M = []
        predecessor = []
        for i in xrange(seq_len):
            M.append(0)
            predecessor.append(0)
        M.append(0)
        M[0] = 0 # not really used
        longest_subseq_len = 0
        for i in xrange(0, seq_len):
            lo = 1
            hi = longest_subseq_len
            while lo <= hi:
              mid = (lo+hi)/2
              if sequence[M[mid]] < sequence[i]:
                  lo = mid+1
              else:
                hi = mid-1
            newL = lo
            predecessor[i] = M[newL-1]
            if newL > longest_subseq_len:
                M[newL] = i
                longest_subseq_len = newL
            elif sequence[i] < sequence[M[newL]]:
                M[newL] = i
        subseq = []
        for i in xrange(longest_subseq_len):
          subseq.append(0)
        k = M[longest_subseq_len]
        for i in xrange(longest_subseq_len):
            subseq[longest_subseq_len-i-1] = sequence[k]
            k = predecessor[k]
        return subseq

    def findCoverage(self, fullRange, length):
        newRange=list()
        count =1
        newRange.append(fullRange[0])
        overlap=0
        for i in range(1,len(fullRange)):
            if(newRange[count-1][1]>=fullRange[i][0] and newRange[count-1][1]<=fullRange[i][1]):
                #count+=1
                overlap+=newRange[count-1][1]-fullRange[i][0]
                newRange[count-1]=[newRange[count-1][0],fullRange[i][1]]

            else:
                newRange.append(fullRange[i])
                count+=1

        sums=0
        for i in range(len(newRange)):
            sums+=newRange[i][1]-newRange[i][0]+1
        return float(sums)/length, overlap

    def merge_intervals(self, unmerged_list):
        unmerged_list.sort(key=lambda t: (t[0], t[1]))
        merged_list = []
        stt = 0
        while stt < len(unmerged_list):
            this_stt = stt
            while stt + 1 < len(unmerged_list) - 1 and unmerged_list[this_stt][1] > unmerged_list[stt+1][0]:
                stt += 1
            merged_list.append((unmerged_list[this_stt][0], unmerged_list[stt][1]))
            stt += 1
        return merged_list

    '''
    Assembly answers should be given with contigs separated by newline characters.  Any contigs that are smaller than
    or equal to the read length (for now set to 50) are automatically thrown out.
    '''
    def ASSEMBLYgrade(self, stud, key, index):
        cov_grade = 0
        contig_grade = 0
        acc_grade = 0

        keyIndex=1
        i=keyIndex
        while (i < len(key) and key[i][0] != '>'):
            i += 1

        key_genome=''
        stud_contigs=[]

        for ans in key[keyIndex:i+1]:
            key_genome+=str(ans).strip()

        corrected_submit_length = 0
        i = index
        # find end of student answers
        while (i < len(stud) and len(stud[i])>0 and stud[i][0]!='>'):
            if len(stud[i].strip()) > self._key_size:
                corrected_submit_length += len(stud[i])
                stud_contigs.append(stud[i].strip())
            i += 1
        new_i = i

        #default to a grade of zero for students that spam with excessively large answers
        if corrected_submit_length > (3*len(key_genome)) or corrected_submit_length == 0:
            cov_grade = 0
            contig_grade = 0
            acc_grade = 0
            return new_i, cov_grade, contig_grade, acc_grade

        #sort the contigs in reverse order of length
        stud_contigs.sort(key=len, reverse=True)

        #create an index of the answer genome to allow contigs to be matched to the key
        for i in range(len(key_genome)-self._key_size+1):
            key = key_genome[i:i+self._key_size]

            if self._genome_idx.has_key(key):
                self._genome_idx[key].append(i)
            else:
                self._genome_idx[key] = list()
                self._genome_idx[key].append(i)
        #GRADE:
        #for each contig, do the following:
        #1) Determine accuracy: the longest perfect match / total length of contig
        #2) Determine coverage intervals: (Start,Stop) for each interval that matched to the genome
        #   These coverage intervals will allow us to both calculate the true coverage of the genome and
        #   the average corrected (not allowing for errors) contig size
        #contig_intervals will match each section of the contig with its longest match (not all matches, so in the case of repetitive structures it may not give all credit possible)
        contig_intervals = []
        #coverage_intervals = [] #using this would make the grading more lenient, it matches each section of the contigs to all possible matched
        contig_accuracies = []

        for contig in stud_contigs:
            temp_interval_list = []
            #temp_coverage_list = []

            next_stt = 0
            #repeat this for each sub contig found
            while next_stt < (len(contig) - self._key_size):

                genome_posn_list = self._genome_idx.get(contig[next_stt:next_stt+self._key_size])
                if genome_posn_list:
                    best_length = 0
                    best_stt = 0
                    best_end = 0

                    for genome_posn in genome_posn_list:
                        this_stt = next_stt
                        next_stt = next_stt + self._key_size
                        while next_stt <= (len(contig) - self._key_size) and self._genome_idx.get(contig[next_stt:next_stt+self._key_size]):
                            next_stt = next_stt + self._key_size
                        gen_idx = genome_posn + (next_stt - this_stt)
                        while next_stt < len(contig) and gen_idx < len(key_genome):
                            if key_genome[gen_idx] == contig[next_stt]:
                                gen_idx += 1
                                next_stt += 1
                            else:
                                break
                        if (next_stt - this_stt) > best_length:
                            best_length = next_stt - this_stt
                            best_stt = genome_posn
                            best_end = gen_idx

                    #only include portions that haven't been mapped by previous contigs
                    #overlap is possible with multiple intervals, so first we check for all possible overlaps, then we add the non-overlapping portions only
                    overlap_list = []
                    for interval in contig_intervals:
                        #if there is overlap
                        if interval[0] < best_end and interval[1] > best_stt:
                            #if it is contained, then entire interval is overlap
                            if best_stt >= interval[0] and best_end <= interval[1]:
                                overlap_list.append((best_stt, best_end))
                                break
                            if best_stt <= interval[0] and best_end >= interval[1]:
                                overlap_list.append((interval[0], interval[1]))
                            elif best_stt > interval[0]:
                                overlap_list.append((best_stt, interval[1]))
                            elif best_end < interval[1]:
                                overlap_list.append((best_end, interval[1]))

                    if len(overlap_list) == 0:
                        temp_interval_list.append((best_stt, best_end))
                    else:
                        overlap_list = self.merge_intervals(overlap_list)
                        curr_stt = best_stt
                        for i in range(len(overlap_list)):
                            if curr_stt < overlap_list[i][0]:
                                temp_interval_list.append((curr_stt, overlap_list[i][0]))
                            if i == len(overlap_list) - 1 and best_end > overlap_list[i][1]:
                                temp_interval_list.append((overlap_list[i][1], best_end))
                            curr_stt = overlap_list[i][1]
                        del overlap_list
                else:
                    next_stt = next_stt + 1

            #merge any overlapping subcontigs
            #temp_interval_list.sort(key=lambda t: (t[0], t[1]))
            #temp_coverage_list.sort(key=lambda t: (t[0], t[1]))
            #merged_temp_coverage_list = []
            #merged_temp_interval_list = []
            #stt = 0
            #while stt < len(temp_interval_list):
            #    this_stt = stt
            #    while stt + 1 < len(temp_interval_list) - 1 and temp_interval_list[this_stt][1] > temp_interval_list[stt+1][0]:
            #        stt += 1
            #    merged_temp_interval_list.append((temp_interval_list[this_stt][0], temp_interval_list[stt][1]))
            #    stt += 1
            #temp_interval_list = merged_temp_interval_list
            #del merged_temp_interval_list
            temp_interval_list = self.merge_intervals(temp_interval_list)

            #stt = 0
            #while stt < len(temp_coverage_list):
            #    this_stt = stt
            #    while stt + 1 < len(temp_coverage_list) - 1 and temp_coverage_list[this_stt][1] > temp_coverage_list[stt+1][0]:
            #        stt += 1
            #    merged_temp_coverage_list.append((temp_coverage_list[this_stt][0], temp_coverage_list[stt][1]))
            #    stt += 1
            #temp_coverage_list = merged_temp_coverage_list
            #del merged_temp_coverage_list

            #find the largest sub contig and use it to compute accuracy
            temp_interval_list.sort(key=lambda t: t[1] - t[0], reverse=True)
            #if len(temp_interval_list) > 1:
            #    length1 = temp_interval_list[0][1] - temp_interval_list[0][0]
            #    length2 = temp_interval_list[1][1] - temp_interval_list[1][0]
            if len(temp_interval_list) > 0:
                length1 = temp_interval_list[0][1] - temp_interval_list[0][0]
                length2 = 0
            else:
                length1 = 0
                length2 = 0

            #record accuracy and store the corrected contigs
            contig_accuracies.append((length1+length2)/float(len(contig)))
            for interval in temp_interval_list:
                if interval[1] - interval[0] > self._key_size:
                    contig_intervals.append((interval[0], interval[1]))
            #for interval in temp_coverage_list:
            #    coverage_intervals.append((interval[0], interval[1]))

        #merge overlapping student contigs
        #contig_intervals.sort(key=lambda t: (t[0], t[1]))
        #coverage_intervals.sort(key=lambda t: (t[0], t[1]))

        #merged_interval_list = []
        #merged_coverage_list = []
        #stt = 0
        #while stt < len(contig_intervals):
        #    this_stt = stt
        #    while stt + 1 < len(contig_intervals) and contig_intervals[this_stt][1] > contig_intervals[stt+1][0]:
        #        stt += 1
        #    merged_interval_list.append((contig_intervals[this_stt][0], contig_intervals[stt][1]))
        #    stt += 1
        #contig_intervals = merged_interval_list
        #del merged_interval_list
        contig_intervals = self.merge_intervals(contig_intervals)
        contig_intervals.sort(key=lambda t: t[1] - t[0], reverse=True)
        #stt = 0
        #while stt < len(coverage_intervals):
        #    this_stt = stt
        #    while stt + 1 < len(coverage_intervals) and coverage_intervals[this_stt][1] > coverage_intervals[stt+1][0]:
        #        stt += 1
        #    merged_coverage_list.append((coverage_intervals[this_stt][0], coverage_intervals[stt][1]))
        #    stt += 1
        #coverage_intervals = merged_coverage_list
        #del merged_coverage_list

        contig_grade = 0

        #determine coverage and contig grades
        #student_coverage = 0
        #for interval in coverage_intervals:
        #    int_length = interval[1] - interval[0]
        #    student_coverage += int_length

        student_coverage = 0
        total_interval_length = 0
        for interval in contig_intervals:
            int_length = interval[1] - interval[0]
            #int_weight = int_length / float(len(key_genome))
            #contig_grade += int_length * int_weight
            #total_interval_length += int_length
            student_coverage += int_length
        cov_grade = min(student_coverage / float(len(key_genome)), 1)
        #contig_grade = contig_grade / float(total_interval_length)

        #an error in the above contig grade caused a score over 100%, so for now we are falling back to N50 as our contig grading until more data can be collected
        #compute N50 grade
        assembly_size = float(len(key_genome))
        n50_contig_size = 0
        accumulated_size = 0
        for interval in contig_intervals:
            int_length = interval[1] - interval[0]
            accumulated_size += int_length
            if n50_contig_size == 0 and accumulated_size >= student_coverage / 2.0:
                n50_contig_size = int_length
        contig_grade = min(float(n50_contig_size) / assembly_size, 1)
        acc_grade = round(sum(contig_accuracies) / float(len(contig_accuracies)),2)

        return new_i, cov_grade, contig_grade, acc_grade
        '''

        fullRange.sort(key= lambda fullRange:(int(fullRange[1])-int(fullRange[0])), reverse=True)

        sums=0
        N50=0
        avgN50=0

        for i in range(len(fullRange)):
            N50=fullRange[i][1]-fullRange[i][0]+1
            sums+=N50
            avgN50=sums/i
            if(sums>=len(key_answers)/2):
                break

        score= min(avgN50/(len(key_answers)/4.0),1)*0.5'''

    def eval(self, answerKey, studentAns):

        studAns = [line.rstrip() for line in studentAns]
        studentAns.close()

        for i in range(0,len(studAns)-1):
            if (studAns[i][0]==">"):
                filename = studAns[i+1]
                filename=filename.translate(None,'\n>')
                break
        ansKey = [line.rstrip() for line in answerKey]
        answerKey.close()

        self.create_key_idx(ansKey)

        copyGrade=0
        invGrade=0
        insertGrade=0
        deleteGrade=0
        snpGrade=0
        strGrade=0
        aluGrade=0
        aCoverage=0
        aAccuracy=0
        aContig=0

        for i in range(len(studAns)):
            new_i = -1
            if (studAns[i][0:5]==">CNV") and self._key_idx[">CNV"] != -1:
                new_i, copyGrade=self.COPYgrade(studAns,ansKey,i+1)
            if (studAns[i][0:10]==">INV"):
                new_i, invGrade=self.INVgrade(studAns,ansKey,i+1)
            if (studAns[i][0:7]==">INS"):
                new_i, insertGrade =self.INDELgrade(studAns,ansKey,i+1, ">INS")
            if (studAns[i][0:7]==">DEL"):
                new_i, deleteGrade=self.INDELgrade(studAns,ansKey,i+1, ">DEL")
            if (studAns[i][0:4]==">SNP"):
                new_i, snpGrade=self.SNPgrade(studAns,ansKey,i+1)
            if (studAns[i][0:4]==">STR"):
                new_i, strGrade=self.STRgrade(studAns,ansKey,i+1)
            if (studAns[i][0:4]==">ALU"):
                new_i, aluGrade=self.INDELgrade(studAns,ansKey,i+1, ">ALU")
            if (studAns[i][0:9]==">ASSEMBLY"):
                new_i, aCoverage, aContig, aAccuracy=self.ASSEMBLYgrade(studAns, ansKey, i+1)
            #allows loop to skip over the portions that were already evaluated
            if new_i != -1:
                i = new_i - 1

        #if the answer key didn't have entries for a given type, then set the value to -1 so that we can correctly compute the total score
        if not self._valid_types[">CNV"]:
            copyGrade=-1
        if not self._valid_types[">INV"]:
            invGrade=-1
        if not self._valid_types[">INS"]:
            insertGrade=-1
        if not self._valid_types[">DEL"]:
            deleteGrade=-1
        if not self._valid_types[">SNP"]:
            snpGrade=-1
        if not self._valid_types[">STR"]:
            strGrade=-1
        if not self._valid_types[">ALU"]:
            aluGrade=-1
        if not self._valid_types[">ASSEMBLY"]:
            aCoverage=-1
            aContig=-1
            aAccuracy=-1

        grades = {'SNP': snpGrade,'INDEL':(insertGrade+deleteGrade)/2,'CNV': copyGrade, 'INV': invGrade,
                  'STR': strGrade, 'ALU': aluGrade, 'A_COV':aCoverage, 'A_ACC':aAccuracy, 'A_CON':aContig }
        return grades

class TestClass(unittest.TestCase):

    def setUp(self):
        self.gen = Eval()
        self.data_directory = 'test_data'

        self.nonassembly_10k_key = 'nonassembly_10k_key_wALU.txt'
        self.assembly_10k_key = 'assembly_10k_key.txt'
        self.nonassembly_1m_key = 'nonassembly_1m_key_wALU.txt'
        self.nonassembly_100m_key = 'nonassembly_100m_key.txt'

        self.assembly_10k_100percent = 'assembly_10k_100percent.txt'

        self.empty_answer_file = 'empty.txt'
        self.headers_only_answer_file = 'headers_only.txt'

        self.duplicates_answer_file = 'nonassembly_10k_duplicates_wALU.txt'

        self.nonassembly_10k_50percent = 'nonassembly_10k_50percent_wALU.txt'

        self.assembly_headers_only = 'assembly_headers_only.txt'

        self.assembly_10k_reads = 'assembly_10k_reads.txt'
        self.assembly_2contigs = 'assembly_10k_2contigs.txt'
        self.assembly_5contigs = 'assembly_10k_5contigs.txt'
        self.assembly_10contigs = 'assembly_10k_10contigs.txt'

        self.assembly_10k_50percent_accuracy = 'assembly_10k_50percent_accuracy.txt'

        self.assembly_undersize_50percent_5contigs = 'assembly_undersize_50percent_5contigs.txt'
        self.assembly_oversize_150percent_5contigs = 'assembly_oversize_150percent_5contigs.txt'
        self.assembly_oversize_200percent_5contigs = 'assembly_oversize_200percent_5contigs.txt'
        self.assembly_oversize_300percent_5contigs = 'assembly_oversize_300percent_5contigs.txt'

        self.assembly_results_padding = 'assembly_results_padding.txt'
        self.assembly_previous_error1_key = 'assembly_studentSubmissionsKey_1.txt'
        self.assembly_previous_error1_1 = '37_39'

        self.non_assembly_additional_10k_100 = 'ans_hw0_W_0_chr_1.txt'

    def test_100_percent_10k_answer(self):
        student_ans = os.path.join(self.data_directory, self.nonassembly_10k_key)
        ans_key = os.path.join(self.data_directory, self.nonassembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], 1)
                self.assertEqual(grades['INDEL'], 1)
                self.assertEqual(grades['CNV'], 1)
                self.assertEqual(grades['INV'], 1)
                self.assertEqual(grades['STR'], 1)
                self.assertEqual(grades['ALU'], 1)
                self.assertEqual(grades['A_COV'], -1)
                self.assertEqual(grades['A_ACC'], -1)
                self.assertEqual(grades['A_CON'], -1)

    def test_100_percent_10k_answer2(self):
        student_ans = os.path.join(self.data_directory,  self.non_assembly_additional_10k_100)
        ans_key = os.path.join(self.data_directory,  self.non_assembly_additional_10k_100)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], 1)
                self.assertEqual(grades['INDEL'], 1)
                self.assertEqual(grades['CNV'], 1)
                self.assertEqual(grades['INV'], 1)
                self.assertEqual(grades['STR'], 1)
                self.assertEqual(grades['ALU'], 1)
                self.assertEqual(grades['A_COV'], -1)
                self.assertEqual(grades['A_ACC'], -1)
                self.assertEqual(grades['A_CON'], -1)

    def test_100_percent_1m_answer(self):
        student_ans = os.path.join(self.data_directory, self.nonassembly_1m_key)
        ans_key = os.path.join(self.data_directory, self.nonassembly_1m_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], 1)
                self.assertEqual(grades['INDEL'], 1)
                self.assertEqual(grades['CNV'], 1)
                self.assertEqual(grades['INV'], 1)
                self.assertEqual(grades['STR'], 1)
                self.assertEqual(grades['ALU'], 1)
                self.assertEqual(grades['A_COV'], -1)
                self.assertEqual(grades['A_ACC'], -1)
                self.assertEqual(grades['A_CON'], -1)

    def test_empty_answer(self):
        student_ans = os.path.join(self.data_directory, self.empty_answer_file)
        ans_key = os.path.join(self.data_directory, self.nonassembly_1m_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], 0)
                self.assertEqual(grades['INDEL'], 0)
                self.assertEqual(grades['CNV'], 0)
                self.assertEqual(grades['INV'], 0)
                self.assertEqual(grades['STR'], 0)
                self.assertEqual(grades['ALU'], 0)
                self.assertEqual(grades['A_COV'], -1)
                self.assertEqual(grades['A_ACC'], -1)
                self.assertEqual(grades['A_CON'], -1)

    def test_answer_with_headers_but_no_data(self):
        student_ans = os.path.join(self.data_directory, self.headers_only_answer_file)
        ans_key = os.path.join(self.data_directory, self.nonassembly_1m_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], 0)
                self.assertEqual(grades['INDEL'], 0)
                self.assertEqual(grades['CNV'], 0)
                self.assertEqual(grades['INV'], 0)
                self.assertEqual(grades['STR'], 0)
                self.assertEqual(grades['ALU'], 0)
                self.assertEqual(grades['A_COV'], -1)
                self.assertEqual(grades['A_ACC'], -1)
                self.assertEqual(grades['A_CON'], -1)

    def test_duplicate_answers(self):
        student_ans = os.path.join(self.data_directory, self.duplicates_answer_file)
        ans_key = os.path.join(self.data_directory, self.nonassembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertTrue(grades['SNP'] < .7)
                self.assertTrue(grades['INDEL'] < .7)
                self.assertTrue(grades['CNV'] < .7)
                self.assertTrue(grades['INV'] < .7)
                self.assertTrue(grades['STR'] < .7)
                self.assertTrue(grades['ALU'] < .7)
                self.assertTrue(grades['A_COV'] == -1)
                self.assertTrue(grades['A_ACC'] == -1)
                self.assertTrue(grades['A_CON'] == -1)

                self.assertTrue(grades['SNP'] > .6)
                self.assertTrue(grades['INDEL'] > .6)
                self.assertTrue(grades['CNV'] > .5)
                self.assertTrue(grades['INV'] > .6)
                self.assertTrue(grades['STR'] > .6)
                self.assertTrue(grades['ALU'] > .6)
                self.assertTrue(grades['A_COV'] == -1)
                self.assertTrue(grades['A_ACC'] == -1)
                self.assertTrue(grades['A_CON'] == -1)

    def test_50_percent_answers(self):
        student_ans_10k = os.path.join(self.data_directory, self.nonassembly_10k_50percent)
        ans_key = os.path.join(self.data_directory, self.nonassembly_10k_key)
        with open(student_ans_10k,'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertTrue(grades['SNP'] < .55)
                self.assertTrue(grades['INDEL'] < .55)
                self.assertTrue(grades['CNV'] < .55)
                self.assertTrue(grades['INV'] < .55)
                self.assertTrue(grades['STR'] < .55)
                self.assertTrue(grades['ALU'] < .55)
                self.assertTrue(grades['A_COV'] == -1)
                self.assertTrue(grades['A_ACC'] == -1)
                self.assertTrue(grades['A_CON'] == -1)

                self.assertTrue(grades['SNP'] > .45)
                self.assertTrue(grades['INDEL'] > .45)
                self.assertTrue(grades['CNV'] > .45)
                self.assertTrue(grades['INV'] > .45)
                self.assertTrue(grades['STR'] > .45)
                self.assertTrue(grades['ALU'] > .45)
                self.assertTrue(grades['A_COV'] == -1)
                self.assertTrue(grades['A_ACC'] == -1)
                self.assertTrue(grades['A_CON'] == -1)

    def test_assembly_100percent(self):
        student_ans = os.path.join(self.data_directory, self.assembly_10k_100percent)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertEqual(grades['A_COV'], 1)
                self.assertEqual(grades['A_ACC'], 1)
                self.assertEqual(grades['A_CON'], 1)

    def test_assembly_vary_contigs(self):
        student_ans2 = os.path.join(self.data_directory, self.assembly_2contigs)
        student_ans5 = os.path.join(self.data_directory, self.assembly_5contigs)
        student_ans10 = os.path.join(self.data_directory, self.assembly_10contigs)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans2, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades2 = self.gen.eval(answer, student)
                self.assertEqual(grades2['SNP'], -1)
                self.assertEqual(grades2['INDEL'], -1)
                self.assertEqual(grades2['CNV'], -1)
                self.assertEqual(grades2['INV'], -1)
                self.assertEqual(grades2['STR'], -1)
                self.assertEqual(grades2['ALU'], -1)
                self.assertTrue(grades2['A_COV'] == 1)
                self.assertTrue(grades2['A_ACC'] == 1)
                self.assertTrue(grades2['A_CON'] < .6)
                self.assertTrue(grades2['A_CON'] > .4)
        with open(student_ans5, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades5 = self.gen.eval(answer, student)
                self.assertEqual(grades5['SNP'], -1)
                self.assertEqual(grades5['INDEL'], -1)
                self.assertEqual(grades5['CNV'], -1)
                self.assertEqual(grades5['INV'], -1)
                self.assertEqual(grades5['STR'], -1)
                self.assertEqual(grades5['ALU'], -1)
                self.assertTrue(grades5['A_COV'] == 1)
                self.assertTrue(grades5['A_ACC'] == 1)
                self.assertTrue(grades5['A_CON'] < .25)
                self.assertTrue(grades5['A_CON'] > .2)
        with open(student_ans10, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades10 = self.gen.eval(answer, student)
                self.assertEqual(grades10['SNP'], -1)
                self.assertEqual(grades10['INDEL'], -1)
                self.assertEqual(grades10['CNV'], -1)
                self.assertEqual(grades10['INV'], -1)
                self.assertEqual(grades10['STR'], -1)
                self.assertEqual(grades10['ALU'], -1)
                self.assertTrue(grades10['A_COV'] == 1)
                self.assertTrue(grades10['A_ACC'] == 1)
                self.assertTrue(grades10['A_CON'] < .2)
                self.assertTrue(grades10['A_CON'] > 0)
        self.assertTrue(grades2['A_CON'] > grades5['A_CON'])
        self.assertTrue(grades5['A_CON'] > grades10['A_CON'])

    def test_assembly_50percent_accuracy(self):
        student_ans = os.path.join(self.data_directory, self.assembly_10k_50percent_accuracy)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertTrue(grades['A_COV'] > .95) #the coverage is actually 100% but the contigs in this test are so small that there are some aligned improperly due to repetitive structures in the test data
                self.assertTrue(grades['A_COV'] <= 1)
                self.assertTrue(grades['A_CON'] > 0)
                self.assertTrue(grades['A_CON'] <= 1)
                self.assertTrue(grades['A_ACC'] >= .45)
                self.assertTrue(grades['A_ACC'] <= .55)

    def test_assembly_results_padding(self):
        student_ans = os.path.join(self.data_directory, self.assembly_results_padding)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertTrue(grades['A_COV'] < .75)
                self.assertTrue(grades['A_COV'] > 0)
                self.assertTrue(grades['A_CON'] > .2)
                self.assertTrue(grades['A_CON'] < .3)
                self.assertTrue(grades['A_ACC'] < 1)
                self.assertTrue(grades['A_ACC'] > 0)

    def test_assembly_coverage(self):
        student_ans50 = os.path.join(self.data_directory, self.assembly_undersize_50percent_5contigs)
        student_ans150 = os.path.join(self.data_directory, self.assembly_oversize_150percent_5contigs)
        student_ans200 = os.path.join(self.data_directory, self.assembly_oversize_200percent_5contigs)
        student_ans300 = os.path.join(self.data_directory, self.assembly_oversize_300percent_5contigs)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans50, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades50 = self.gen.eval(answer, student)
                self.assertEqual(grades50['SNP'], -1)
                self.assertEqual(grades50['INDEL'], -1)
                self.assertEqual(grades50['CNV'], -1)
                self.assertEqual(grades50['INV'], -1)
                self.assertEqual(grades50['STR'], -1)
                self.assertEqual(grades50['ALU'], -1)
                self.assertTrue(grades50['A_COV'] >= .40)
                self.assertTrue(grades50['A_COV'] <= .65)
                self.assertTrue(grades50['A_ACC'] >= .85)
                self.assertTrue(grades50['A_ACC'] <= 1)
                self.assertTrue(grades50['A_CON'] > .1)
                self.assertTrue(grades50['A_CON'] < .2)

        with open(student_ans150, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades150 = self.gen.eval(answer, student)
                self.assertEqual(grades150['SNP'], -1)
                self.assertEqual(grades150['INDEL'], -1)
                self.assertEqual(grades150['CNV'], -1)
                self.assertEqual(grades150['INV'], -1)
                self.assertEqual(grades150['STR'], -1)
                self.assertEqual(grades150['ALU'], -1)
                self.assertTrue(grades150['A_COV'], 1)
                self.assertTrue(grades150['A_ACC'] > .30)
                self.assertTrue(grades150['A_ACC'] <= 1)
                self.assertTrue(grades150['A_CON'] > 0)
                self.assertTrue(grades150['A_COV'] <= 1)

        with open(student_ans200, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades200 = self.gen.eval(answer, student)
                self.assertEqual(grades200['SNP'], -1)
                self.assertEqual(grades200['INDEL'], -1)
                self.assertEqual(grades200['CNV'], -1)
                self.assertEqual(grades200['INV'], -1)
                self.assertEqual(grades200['STR'], -1)
                self.assertEqual(grades200['ALU'], -1)
                self.assertTrue(grades200['A_COV'], 1)
                self.assertTrue(grades200['A_ACC'] > .20)
                self.assertTrue(grades200['A_ACC'] <= 1)
                self.assertTrue(grades200['A_CON'] <= 1)
                self.assertTrue(grades200['A_CON'] > 0)
        with open(student_ans300, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades300 = self.gen.eval(answer, student)
                self.assertEqual(grades300['SNP'], -1)
                self.assertEqual(grades300['INDEL'], -1)
                self.assertEqual(grades300['CNV'], -1)
                self.assertEqual(grades300['INV'], -1)
                self.assertEqual(grades300['STR'], -1)
                self.assertEqual(grades300['ALU'], -1)
                self.assertEqual(grades300['A_COV'], 0)
                self.assertEqual(grades300['A_ACC'], 0)
                self.assertEqual(grades300['A_CON'], 0)

    def test_assembly_empty(self):
        student_ans = os.path.join(self.data_directory, self.empty_answer_file)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertEqual(grades['A_COV'], 0)
                self.assertEqual(grades['A_ACC'], 0)
                self.assertEqual(grades['A_CON'], 0)

    def test_assembly_headers_only(self):
        student_ans = os.path.join(self.data_directory, self.assembly_headers_only)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertEqual(grades['A_COV'], 0)
                self.assertEqual(grades['A_ACC'], 0)
                self.assertEqual(grades['A_CON'], 0)

    def test_assembly_reads_only(self):
        student_ans = os.path.join(self.data_directory, self.assembly_10k_reads)
        ans_key = os.path.join(self.data_directory, self.assembly_10k_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertEqual(grades['A_COV'], 0)
                self.assertEqual(grades['A_ACC'], 0)
                self.assertEqual(grades['A_CON'], 0)

    def test_assembly_previous_error(self):
        student_ans = os.path.join(self.data_directory, self.assembly_previous_error1_1)
        ans_key = os.path.join(self.data_directory, self.assembly_previous_error1_key)
        with open(student_ans, 'r') as student:
            with open(ans_key, 'r') as answer:
                grades = self.gen.eval(answer, student)
                self.assertEqual(grades['SNP'], -1)
                self.assertEqual(grades['INDEL'], -1)
                self.assertEqual(grades['CNV'], -1)
                self.assertEqual(grades['INV'], -1)
                self.assertEqual(grades['STR'], -1)
                self.assertEqual(grades['ALU'], -1)
                self.assertTrue(grades['A_COV'] < .05)
                self.assertTrue(grades['A_CON'] > 0)
                self.assertTrue(grades['A_CON'] < .01)
                self.assertTrue(grades['A_ACC'] == .5)