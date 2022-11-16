
# ----------------------------------------------------------------------------------------------------------------------
# all the util functions
# ----------------------------------------------------------------------------------------------------------------------

def mutation_iducer(spike_seq, mut_list):
    """
    this function will take in a og seq and a mutation list
    and will return a spike seq with the mutations inserted
    :param spike_seq: og spike sequence
    :param mut_list: a list of mutations
    :return: the mutated seq
    """
    spike_list = list(spike_seq)
    for m in mut_list:
        pos = int(m[1:-1])
        if m[0] != spike_list[pos - 1]:
            print(m[0])
            print(spike_list[pos - 1])
            print('index not correct with mutation ', m)
            break

        mut = m[-1]
        spike_list[pos - 1] = mut
    seq_str = ''.join(spike_list)
    return seq_str


def phisical_change(row_data, phis, ref_df_f):
    Hydropathicity = {'Ala': 1.800, 'Arg': -4.500, 'Asn': -3.500, 'Asp': -3.500, 'Cys': 2.500, 'Gln': -3.500,
                      'Glu': -3.500, 'Gly': -0.400, 'His': -3.200, 'Ile': 4.500, 'Leu': 3.800, 'Lys': -3.900,
                      'Met': 1.900, 'Phe': 2.800, 'Pro': -1.600, 'Ser': -0.800, 'Thr': -0.700, 'Trp': -0.900,
                      'Tyr': -1.300, 'Val': 4.200}
    Bulkiness = {'Ala': 11.500, 'Arg': 14.280, 'Asn': 12.820, 'Asp': 11.680, 'Cys': 13.460, 'Gln': 14.450,
                 'Glu': 13.570, 'Gly': 3.400, 'His': 13.690, 'Ile': 21.400, 'Leu': 21.400, 'Lys': 15.710, 'Met': 16.250,
                 'Phe': 19.800, 'Pro': 17.430, 'Ser': 9.470, 'Thr': 15.770, 'Trp': 21.670, 'Tyr': 18.030, 'Val': 21.570}
    aa_three_to_one = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly',
                       'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
                       'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}

    pos = row_data[1]
    mut_aa = row_data[2]
    prot_name = row_data[0]
    # print(mut_aa)
    prot_seq = ref_df_f.loc[prot_name, 'sequence']
    wt_aa = prot_seq[pos - 1]

    if phis == "B":
        delta_phis = Bulkiness[aa_three_to_one[mut_aa]] - Bulkiness[aa_three_to_one[wt_aa]]
        return abs(delta_phis)
    if phis == "H":
        delta_phis = Hydropathicity[aa_three_to_one[mut_aa]] - Hydropathicity[aa_three_to_one[wt_aa]]
        return abs(delta_phis)


def row_anchor_paprams(row_data):
    mut_pos = row_data[0]
    statr_pos = row_data[2]
    end_pos = row_data[3]
    peptide = row_data[1]
    if mut_pos == end_pos or mut_pos == statr_pos + 1:
        return True
    if mut_pos == statr_pos:
        return 1
    return peptide[2:-1]


def count_number_of_binders(score_list):
    score_list = list(score_list)
    count = 0
    # print(score_list)
    for score in score_list:
        if score >= 2:
            count += 1
    return count


def count_number_of_binders(score_list):
    score_list = list(score_list)
    count = 0
    # print(score_list)
    for score in score_list:
        if score >= 2:
            count += 1
    return count


def get_allele_type(allel_str):
    # a = 'HLA-A02:02'
    return allel_str[4:5]


def read_in_allels(prealent_allel_file_and_path):
    """
    this function will read in the file and return all the allelels in it in the following format:
    'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    so the netHmcpan predicter can read it in
    :param prealent_allel_file_and_path:
    :return: string
    """
    with open(prealent_allel_file_and_path) as file_in:
        lines = ''
        for line in file_in:
            temp_line = line.rstrip() + ','
            if len(temp_line) < 2:
                continue
            else:
                lines += temp_line
    return lines[:-1]
