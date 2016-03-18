import sys
from os.path import join
import zipfile
import time
from collections import defaultdict
from CM124.helpers import read_reference, read_reads


def make_STR_dicts(ref):
    STR_3_dict = defaultdict(list)
    STR_4_dict = defaultdict(list)
    STR_5_dict = defaultdict(list)

    for i in range(len(ref)-5):
        trinuc = ref[i:i+3]
        STR_3_dict[trinuc].append(i)

        fournuc = ref[i:i+4]
        STR_4_dict[fournuc].append(i)

        fivenuc = ref[i:i+5]
        STR_5_dict[fivenuc].append(i)

    return STR_3_dict, STR_4_dict, STR_5_dict


def identify_potential_STRs(seq_dict):
    potential_strs_dict = {}
    str_length = len(seq_dict.keys()[0])
    for short_sequence in seq_dict:
        position_list = seq_dict[short_sequence]

        possible_strs = find_repetititive_sequences(position_list, str_length)
        potential_strs_dict[short_sequence] = possible_strs

    return potential_strs_dict


def find_repetititive_sequences(position_list, short_sequence_length):
    all_sequences = []
    current_sequence = []
    for i in range(len(position_list)-1):
        current_position = position_list[i]
        next_position = position_list[i+1]
        if current_position + short_sequence_length != next_position:
            if current_sequence:
                current_sequence.append(current_position)
                if len(current_sequence) > 3:
                    repeat_begin = current_sequence[0]
                    repeat_end = current_sequence[-1] + short_sequence_length
                    all_sequences.append((repeat_begin,repeat_end))
                current_sequence=[]
            else:
                pass
        else:
            current_sequence.append(current_position)

    if current_sequence:
        current_sequence.append(next_position)
        if len(current_sequence) > 3:
            repeat_begin = current_sequence[0]
            repeat_end = current_sequence[-1] + short_sequence_length
            all_sequences.append((repeat_begin, repeat_end))
    else:
        pass
    return all_sequences


def STR_sequence_output(ref, potential_strs_dict):
    output_strs = []
    for k in potential_strs_dict:
        bounds_list = potential_strs_dict[k]
        for bound in bounds_list:
            lower_bound = bound[0]
            upper_bound = bound[1]
            ref_piece = ref[lower_bound: upper_bound]
            output_str = '{},{}'.format(ref_piece, lower_bound)
            output_strs.append(output_str)
    output = '\n'.join(output_strs)
    return output


def run_reference_STRs(ref):
    output = []
    STR_dict_list = make_STR_dicts(ref)
    for d in STR_dict_list:
        potential_STRs_dict = identify_potential_STRs(d)
        dict_output = STR_sequence_output(ref, potential_STRs_dict)
        print dict_output
        output.append(dict_output)
    return output


def main(input_folder):
    f_base = '{}_chr_1'.format(input_folder)
    reference_fn = join(folder, 'ref_{}.txt'.format(f_base))
    ref = read_reference(reference_fn)
    output_fn = join(folder, 'STRs_{}.txt'.format(f_base))

    STRs_from_ref = run_reference_STRs(ref)
    filter_output(STRs_from_ref)
    final_output = ['>{}'.format(f_base), '>STR'] + STRs_from_ref

    final_STRs = '\n'.join(final_output)
    with open(output_fn, 'w') as f:
        f.write(final_STRs)
        # print final_STRs


def filter_output(results):
    print len(results)
    print results[0]
    # re_results = [line.strip().split(',') for line in results]
    # print len(re_results)
    # print re_results[:3]
    # sorted_results = sorted(results, key = lambda x: )

if __name__ == "__main__":
    fake_ref = 'GCA'*20
    fake_ref2 = 'GGCTACGG' + 'ACG'*10 + 'TTGATCTA'
    fake_ref3 = 'GCA'*20 + 'CTAG'*13 + 'ACGTG'*8
    for r in (fake_ref, fake_ref2, fake_ref3):
        run_reference_STRs(r)

    folder = "hw2undergrad_E_2"
    main(folder)
