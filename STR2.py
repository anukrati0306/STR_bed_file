from Bio import SeqIO
from motif_finder import find_repeating_unit

for record in SeqIO.parse('/home/ccmb/imperfact_repeats/human_ref/chr1.fa', 'fasta'):
    seq = record.seq

    shift_list = [1, 2, 3, 4, 5, 6]
    cutoffs = {
        1: range(0, 3),
        2: range(0, 4),
        3: range(0, 3),
        4: range(0, 5),
        5: range(0, 5),
        6: range(0, 7)
    }  
    lengths = {
        1: range(12, 200),
        2: range(12, 200),
        3: range(12, 200),
        4: range(12, 200),
        5: range(15, 250),
        6: range(6, 360)
    }

    with open('/home/ccmb/imperfact_repeats/bed_files/code2_chr1.bed', 'w') as bed_file:
        for shift in shift_list:
            shift_seq = seq[shift:] + 'N' * shift
            match = ''
            for i in range(len(seq)):
                if seq[i] == shift_seq[i]:
                    match += '1'
                else:
                    match += '0'

            i = 0
            n = len(match)
            while i < n:
                if match[i] == '1':
                    start = i
                    count_zeros = 0
                    cutoff_range = cutoffs[shift]
                    while i < n and (match[i] == '1' or (match[i] == '0' and count_zeros in cutoff_range)):
                        if match[i] == '0':
                            count_zeros += 1
                        i += 1
                    end = i - 1
                    length = (end - start) + 1
                    length_range = lengths[shift]
                    if length in length_range:
                        bed_entry_start = start
                        bed_entry_end = end + 1
                        bed_entry_string = seq[start:end + 1]
                        motifs = find_repeating_unit(bed_entry_string, shift)
                        if motifs:
                            motif = motifs[0]
                            bed_file.write(f"chr1\t{bed_entry_start}\t{bed_entry_end}\tShift={shift}\t{bed_entry_string}\t{match[start:end]}\t{motif}\n")
                else:
                    i += 1
