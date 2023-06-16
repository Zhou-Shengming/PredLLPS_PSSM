from __future__ import division
import re
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
from gensim.models.doc2vec import TaggedDocument
from gensim.models.doc2vec import Doc2Vec
import random as python_random




np.random.seed(42)
python_random.seed(42)


def extract_seq(fasta_fname):
    """
    param fasta_fname: Fasta file address.

    """
    seq_list = []
    for seq_record in SeqIO.parse(fasta_fname, 'fasta'):

        bool=re.search(r'[^ACDEFGHIKLMNPQRSTVWY]',(str(seq_record.seq)).upper())
        if bool:continue
        seq=(str(seq_record.seq)).upper()
        seq_list.append(seq)
    return seq_list


def excel_read(excel_fname, tag):
    """
    param excel_fname: Excel file address.
    param tag: The excel of the tag.

    """
    df = pd.read_excel(excel_fname, sheet_name=0)
    seq_list = df.loc[:, tag].values
    # y = df.iloc[:, 0]
    return seq_list


def get_documents(seq_list, k, s):
    """
    param seq_list: The sequence to be processed.
    param k: The length of the word.
    param s:Step length.

    Example sequence: MALFFFNNN
        k = 3
        s = 1
        ['MAL', 'ALF', 'LFF', 'FFF', 'FFN', 'FNN', 'NNN']

    """
    documents = []

    for seq, seq_id in zip(seq_list, seq_list):
        codes = seq
        words = [codes[i: i + k] for i in range(0, len(codes) - (k - 1), s)]
        documents.append(TaggedDocument(words, tags=[seq_id]))
    return documents


def training_doc2vec(documents, vec_size, alpha, min_alpha, dm, min_count, dm_mean, dm_concat, epochs, fname, seed):

    model = Doc2Vec(documents,
                    size=vec_size,
                    alpha=alpha,
                    min_alpha=min_alpha,
                    dm=dm,
                    min_count=min_count,
                    dm_mean=dm_mean,
                    dm_concat=dm_concat,
                    epochs=epochs,
                    seed=seed)

    model.save(fname)
    return


def save_csv(data, file_name):
    data1 = pd.DataFrame(data=data)

    data1.to_csv(file_name, encoding='gbk')
    return


fasta_fname = 'input/example.fasta'

k = 3
s = 1
alpha = 0.025
min_alpha = 0.00025
dm = 1
min_count = 1
dm_mean = 1
dm_concat = 0
epochs = 10
seed = 42


seq_list = extract_seq(fasta_fname)
documents = get_documents(seq_list, k, s)


# training
for i in [2**3, 2**4, 2**5, 2**6, 2**7, 2**8, 2**9, 2**10]:
    vec_size = i
    fname = 'doc2vec_'+str(i)+'.model'

    start = time.time()

    training_doc2vec(documents, vec_size, alpha, min_alpha, dm, min_count, dm_mean, dm_concat, epochs, fname, seed)

    end = time.time()

    print('doc2vec'+str(i)+'.model Running time: %s Seconds'%(end-start))