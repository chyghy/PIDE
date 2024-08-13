# Change sequences into tokens
from Bio import SeqIO
import esm
import pandas as pd


def get_token(input_data):
    alphabet = esm.pretrained.esm2_t33_650M_UR50D()[1]
    batch_converter = alphabet.get_batch_converter()

    data = []
    for seq_record in SeqIO.parse(input_data, "fasta"):
        data.append((seq_record.id, str(seq_record.seq).replace('J', 'X').rstrip('*')))
    batch_tokens = pd.DataFrame(batch_converter(data)[2]).iloc[:, :1024]
    return(batch_tokens)
