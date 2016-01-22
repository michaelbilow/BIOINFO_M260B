def read_reads(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
        all_reads.append(paired_end_reads)
    return all_reads


def read_reference(ref_fn):
    f = open(ref_fn, 'r')
    first_line = True
    output_reference = ''
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        output_reference += line  # We append each line to the output reference string.
    return output_reference