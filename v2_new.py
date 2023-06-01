from Bio import SeqIO
from motif_finder import find_repeating_unit
for record in SeqIO.parse('/home/ccmb/imperfact_repeats/human_ref/chr1.fa', 'fasta'):
    seq = record.seq

    bed_entries = []
    covered_positions = set()
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
        results = []
        while i < n:
            if match[i] == '1':
                start = i
                count_zeros = 0
                cutoff_range = cutoffs[shift]  
                while i < n and (match[i] == '1'or (match[i] == '0' and count_zeros in cutoff_range)):
                    if match[i] == '0':
                        count_zeros += 1
                    i += 1
                end = i - 1
                length = (end - start) + 1
                length_range = lengths[shift]
                if length in length_range:#and (match[start:start + shift] == '1' * (shift)):              
                        if not any(pos in covered_positions for pos in range(start, end + 1)):    
                            results.append((start, end, match[start: end], seq[start:end + 1]))
                            covered_positions.update(range(start, end + 1))

            else:
                i += 1
                
        
        if results:
            for result in results:
                start, end, match_string, string = result
                motifs = find_repeating_unit(string, shift)
                if motifs:
                    motif = motifs[0]  
                    bed_entry = (shift, start, end + 1, match_string, string, motif)
                    bed_entries.append(bed_entry)



    if bed_entries:
        with open('/home/ccmb/imperfact_repeats/bed_files/code2_chr1.bed', 'w') as bed_file:
            for entry in bed_entries:
                end_pos = entry[2]
                shift, start, end_pos, match_string, string, motifs = entry
                bed_file.write(f"chr1\t{start}\t{end_pos}\tShift={shift}\t{string}\t{match_string}\t{motifs}\n")
