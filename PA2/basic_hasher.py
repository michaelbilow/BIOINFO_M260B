from collections import defaultdict, Counter
from helpers.helpers import *
import cPickle as pickle
from os.path import join, exists, splitext

def hash_read(read, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param read: A single read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """
    key_length = len(genome_ht.keys()[0])
    read_pieces = [read[i * key_length: (i + 1) * key_length]
                   for i in range(len(read) / key_length)]

    hashed_read_locations = [genome_ht[read_piece]
                             for read_piece in read_pieces]
    start_positions = [[x - i * key_length for x in hashed_read_locations[i]]
                       for i in range(len(hashed_read_locations))]
    start_counter = Counter()

    for position_list in start_positions:
        start_counter.update(position_list)

    if not start_counter:
        return -1
    else:
        best_alignment_location, best_alignment_count = \
            start_counter.most_common(1)[0]

    if best_alignment_count < 2:
        return -1
    else:
        return best_alignment_location


def hash_genome(reference, key_length):
    """

    :param reference: The reference as a string stored
    :param key_length: The length of keys to use.
    :return:
    """
    genome_hash = defaultdict(list)
    for i in range(len(reference) - key_length):
        ref_piece = reference[i: i + key_length]
        genome_hash[ref_piece].append(i)
    return genome_hash


def build_hash_and_pickle(ref_fn, key_length, force_rebuild=False):
    reference_hash_pkl_fn = '{}_hash.pkl'.format(splitext(ref_fn)[0])
    if exists(reference_hash_pkl_fn) and not force_rebuild:
        ref_genome_hash = pickle.load(open(reference_hash_pkl_fn, 'rb'))
        if len(ref_genome_hash.keys()[0]) == key_length:
            return ref_genome_hash
        else:
            pass
    else:
        pass
    reference = read_reference(ref_fn)
    ref_genome_hash = hash_genome(reference, key_length)
    pickle.dump(ref_genome_hash, open(reference_hash_pkl_fn, 'wb'))
    return ref_genome_hash


if __name__ == "__main__":
    folder = 'hw1_W_2'
    f_base = '{}_chr_1'.format(folder)
    reads_fn = join(folder, 'reads_{}.txt'.format(f_base))
    import time
    start_time = time.clock()
    reference_fn = join(folder, 'ref_{}.txt')
    genome_hash = build_hash_and_pickle(reference_fn, key_length=4)
    # Pickle allows python objects (like this hash table)
    # to be saved to disk. This allows them to be reloaded,
    # rather than rebuilt every time you run the program.
    # If you want to reload your file, COMMENT OUT THE hash_genome line
    # and use the following line:
    # reference_hash = pickle.load(open(reference_pkl_fn,'rb'))

    print time.clock() - start_time

    print sum([len(genome_hash[k]) for k in genome_hash])

    # for k in genome_hash:
    #     print k, genome_hash[k]
    ## Read this output--you can start identifying STRs using this data.
