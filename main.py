from __future__ import division
import os
import re
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from tensorflow.keras import models


def extract_seq(fasta_fname, min_len, max_len):
    """
    param fasta_fname: Fasta file address.

    """
    seq_list, seq_ids, seq_description = [], [], []
    with open(fasta_fname, "r", encoding='utf-8') as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):

            bool = re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', (str(seq_record.seq)).upper())
            if bool: continue
            seq = (str(seq_record.seq)).upper()
            if len(seq) >= min_len and len(seq) <= max_len:
                seq_list.append(seq)
                seq_ids.append(seq_record.id)
                seq_description.append(seq_record.description)
            else:
                print(seq_record.id, 'Not encoding')
                continue
    return seq_list, seq_ids, seq_description


def extract_seq1(fasta_fname, min_len):
    """
    param fasta_fname: Fasta file address.

    """
    seq_list, seq_ids, seq_description = [], [], []
    with open(fasta_fname, "r", encoding='utf-8') as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):

            bool = re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', (str(seq_record.seq)).upper())
            if bool: continue
            seq = (str(seq_record.seq)).upper()
            if len(seq) > min_len:
                seq_list.append(seq)
                seq_ids.append(seq_record.id)
                seq_description.append(seq_record.description)
            else:
                print(seq_record.id, 'Not encoding')
                continue
    return seq_list, seq_ids, seq_description


def read_feature(filepath):
    feature = pd.read_csv(filepath)
    x_feature = np.array(feature)
#    y_train_feature = feature[:, 0]
    return x_feature


def judging_sample_label(y_predict):
    labels = []

    for i in range(len(y_predict)):
        if y_predict[i] >= 0.5:
            label = [1]
            labels.append(label)
        else:
            label = [0]
            labels.append(label)
    return labels


def model_training1(x_1, filepath1):

    model_1_ = models.load_model(filepath1)
    y_predict_1 = model_1_.predict(x_1)

    return y_predict_1


def main():
    base_path = os.path.dirname(__file__)

    parser = argparse.ArgumentParser(
        description='Relevant forecast of LLPS.')
    parser.add_argument('--input_fasta_file', dest='inputfile1', type=str, required=True,
                        help='Input sequence fasta files.')
    parser.add_argument('--input_AB_PSSM_csv_file', dest='inputfile2', type=str, required=True,
                        help='Input the AB PSSM feature file of the sequence obtained in PUSSUM. The file format should be csv.')
    args = parser.parse_args()

    inputfile1 = args.inputfile1
    inputfile2 = args.inputfile2

    file_path1, file_name1 = os.path.split(inputfile1)
    file_path2, file_name2 = os.path.split(inputfile2)

    new_file_name = file_name1.split('.')[0] + '_result.csv'

    fasta_fname = inputfile1
    seq_list, seq_ids, seq_description = extract_seq(fasta_fname, 50, 5000)
    seq_list1, seq_ids1, seq_description1 = extract_seq1(fasta_fname, 0)

    x = read_feature(inputfile2)

    if len(seq_list) == len(seq_list1):
        if len(seq_list) == len(x) and len(x[0]) == 400:
            x = x.reshape(len(x), 20, 20)

            filepath1 = os.path.join(base_path, 'model', 'All model.h5')
            filepath2 = os.path.join(base_path, 'model', 'Self model.h5')
            filepath3 = os.path.join(base_path, 'model', 'Part model.h5')

            y_pre = model_training1(x, filepath1)
            y_pre = np.around(y_pre, 2)
            y_pre_score = np.array(y_pre)
            y_pre_label = judging_sample_label(y_pre)
            y_pre_score = np.array(y_pre_score).T
            y_pre_score = list(y_pre_score[0])
            y_pre_label = np.array(y_pre_label).T
            y_pre_label = list(y_pre_label[0])

            y_pre_self = model_training1(x, filepath2)
            y_pre_self = np.around(y_pre_self, 2)
            y_pre_self_score = np.array(y_pre_self)
            y_pre_self_label = judging_sample_label(y_pre_self)
            y_pre_self_score = np.array(y_pre_self_score).T
            y_pre_self_score = list(y_pre_self_score[0])
            y_pre_self_label = np.array(y_pre_self_label).T
            y_pre_self_label = list(y_pre_self_label[0])

            y_pre_part = model_training1(x, filepath3)
            y_pre_part = np.around(y_pre_part, 2)
            y_pre_part_score = np.array(y_pre_part)
            y_pre_part_label = judging_sample_label(y_pre_part)
            y_pre_part_score = np.array(y_pre_part_score).T
            y_pre_part_score = list(y_pre_part_score[0])
            y_pre_part_label = np.array(y_pre_part_label).T
            y_pre_part_label = list(y_pre_part_label[0])

            results = [seq_description, seq_list, y_pre_score, y_pre_label, y_pre_self_score, y_pre_self_label, y_pre_part_score, y_pre_part_label]
            results = np.array(results)
            results = results.T
            results = pd.DataFrame(results, columns=['Description', 'Sequence', 'LLPS score', 'LLPS label', 'LLPS self score', 'LLPS self label', 'LLPS partner score', 'LLPS partner label'])
            results.to_csv(os.path.join(file_path1, new_file_name), index=False, header=True, escapechar=',')
        else:print('The character of the sequence is problematic. It may not come from possum or the content format may have changed.')
    else:
        print('The sequence length less than 50 or more than 5000 cannot be predicted.')


if __name__ == "__main__":
    main()
