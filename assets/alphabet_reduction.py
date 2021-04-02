# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  get_distances - returns a pandas dataframe distance matrix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd



######-------Alphabet #4---------###########
def AlphaD_reduce(sequence):
    B = ["M", "H", "V", "Y", "N", "D", "I"]
    J = ["Q", "L", "E", "K", "F"]
    U = ["W", "P", "R", "G", "S", "A", "T", "C"]

    lenSeq = len(sequence)
    let = 0
    reducedSeq1 = ""
    while let < lenSeq:
        if sequence[let] in B:
            reducedSeq1 = reducedSeq1 + "B"
        if sequence[let] in J:
            reducedSeq1 = reducedSeq1 + "J"
        if sequence[let] in U:
            reducedSeq1 = reducedSeq1 + "U"
        let = let + 1
    return reducedSeq1


def reduce_dataset(df, output_folder):
    df_red = pd.DataFrame(columns=['Label', 'PDBID', 'Sequence',
                                   'Structure', 'redA', 'redB',
                                   'redC', 'redD', 'redE'])
    for index, row in df.iterrows():
        redB = AlphaB_reduce(row['Sequence'])

        df_red = df_red.append({'Label': row['Label'],
                                'PDBID': row['PDBID'],
                                'Sequence': row['Sequence'],
                                'Structure': row['Structure'],
                                'redB': redB}, ignore_index=True)

    return df_red


def reduce_dataset_st(orig_seq):
    redD = AlphaD_reduce(orig_seq)
    return redD
