import re
import os
import sys
import math
import pandas as pd
from collections import Counter




def read_protein_sequences(file):
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '-', ''.join(array[1:]).upper())
        fasta_sequences.append([header, sequence])
    return fasta_sequences


def AAC(fastas):
    # AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
    AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    header = ['name', 'sequence']
    for i in AA:
        header.append(i)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name, sequence]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings


def DPC(fastas):
    # AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
    AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    triPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
    header = ['name', 'sequence'] + triPeptides
    encodings.append(header)

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name, sequence]
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]]*20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]]*20 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings


def TPC(fastas):
    # AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
    AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
    header = ['name', 'sequence'] + triPeptides
    encodings.append(header)

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name, sequence]
        tmpCode = [0] * 8000
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings


def GAAC(fastas):
    group = {
        'AA_aci': 'DE',
        'AA_bas': 'KHR',
        'AA_neu': 'ACFGILMNPQSTVWY',
        'AA_int': 'FILMV',
        'AA_ext': 'DEHKNQR',
        'AA_amb': 'ACGPSTWY',
        'AA_aro': 'FYW',
        'AA_naro': 'ARNDCQEGHILKMPSTV'
    }

    groupKey = group.keys()

    encodings = []
    header = ['name', 'sequence']
    for key in groupKey:
        header.append(key)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name, sequence]
        count = Counter(sequence)
        myDict = {}
        for key in groupKey:
            for aa in group[key]:
                myDict[key] = myDict.get(key, 0) + count[aa]

        for key in groupKey:
            code.append(myDict[key]/len(sequence))
        encodings.append(code)

    return encodings


def GDPC_acid_base(fastas):
    group = {
        'AA_aci': 'DE',
        'AA_bas': 'KHR',
        'AA_neu': 'ACFGILMNPQSTVWY'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['name', 'sequence'] + dipeptide
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name, sequence]
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 2 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                sequence[j + 1]]] + 1
            sum = sum + 1

        if sum == 0:
            for t in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)

    return encodings


def GDPC_hydrophobicity(fastas):
    group = {
        'AA_int': 'FILMV',
        'AA_ext': 'DEHKNQR',
        'AA_amb': 'ACGPSTWY',
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['name', 'sequence'] + dipeptide
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name, sequence]
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 2 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                sequence[j + 1]]] + 1
            sum = sum + 1

        if sum == 0:
            for t in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)

    return encodings


def GDPC_aromatic(fastas):
    group = {
        'AA_aro': 'FYW',
        'AA_naro': 'ARNDCQEGHILKMPSTV'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['name', 'sequence'] + dipeptide
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name, sequence]
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 2 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                sequence[j + 1]]] + 1
            sum = sum + 1

        if sum == 0:
            for t in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)

    return encodings


def GDPC(fastas):
    GDPC_1 = GDPC_acid_base(fastas)
    GDPC_2 = GDPC_hydrophobicity(fastas)
    GDPC_3 = GDPC_aromatic(fastas)

    GDPC_all = []

    for i in range(len(GDPC_1)):
        GDPC_all_1 = GDPC_1[i] + GDPC_2[i][2:] + GDPC_3[i][2:]
        GDPC_all.append(GDPC_all_1)

    return GDPC_all


def GTPC_acid_base(fastas):
    group = {
        'AA_aci': 'DE',
        'AA_bas': 'KHR',
        'AA_neu': 'ACFGILMNPQSTVWY'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['name', 'sequence'] + triple
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name, sequence]
        myDict = {}
        for t in triple:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 3 + 1):
            myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
            sum = sum +1

        if sum == 0:
            for t in triple:
                code.append(0)
        else:
            for t in triple:
                code.append(myDict[t]/sum)
        encodings.append(code)

    return encodings


def GTPC_hydrophobicity(fastas):
    group = {
        'AA_int': 'FILMV',
        'AA_ext': 'DEHKNQR',
        'AA_amb': 'ACGPSTWY',
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['name', 'sequence'] + triple
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name, sequence]
        myDict = {}
        for t in triple:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 3 + 1):
            myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
            sum = sum +1

        if sum == 0:
            for t in triple:
                code.append(0)
        else:
            for t in triple:
                code.append(myDict[t]/sum)
        encodings.append(code)

    return encodings


def GTPC_aromatic(fastas):
    group = {
        'AA_aro': 'FYW',
        'AA_naro': 'ARNDCQEGHILKMPSTV'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['name', 'sequence'] + triple
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name, sequence]
        myDict = {}
        for t in triple:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 3 + 1):
            myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
            sum = sum +1

        if sum == 0:
            for t in triple:
                code.append(0)
        else:
            for t in triple:
                code.append(myDict[t]/sum)
        encodings.append(code)

    return encodings


def GTPC(fastas):
    GTPC_1 = GTPC_acid_base(fastas)
    GTPC_2 = GTPC_hydrophobicity(fastas)
    GTPC_3 = GTPC_aromatic(fastas)

    GTPC_all = []

    for i in range(len(GTPC_1)):
        GTPC_all_1 = GTPC_1[i] + GTPC_2[i][2:] + GTPC_3[i][2:]
        GTPC_all.append(GTPC_all_1)

    return GTPC_all


def CTDC_Count(seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum


def CTDC(fastas):
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity':        'HQRKNED',
        'polarizability':  'KMHFRYW',
        'charge':          'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess':   'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
    'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
    'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
    'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = ['name', 'sequence']
    for p in property:
        for g in range(1, len(groups) + 1):
            header.append(p + '.G' + str(g))
    encodings.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name, sequence]
        for p in property:
            c1 = CTDC_Count(group1[p], sequence) / len(sequence)
            c2 = CTDC_Count(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        encodings.append(code)
    return encodings


def CTDD_Count(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence) * 100)
                    break
        if myCount == 0:
            code.append(0)
    return code


def CTDD(fastas):
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity':        'HQRKNED',
        'polarizability':  'KMHFRYW',
        'charge':          'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess':   'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
    'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
    'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
    'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')


    encodings = []
    header = ['name', 'sequence']
    for p in property:
        for g in ('1', '2', '3'):
            for d in ['0', '25', '50', '75', '100']:
                header.append(p + '.' + g + '.residue' + d)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name, sequence]
        for p in property:
            code = code + CTDD_Count(group1[p], sequence) + CTDD_Count(group2[p], sequence) + CTDD_Count(group3[p], sequence)
        encodings.append(code)
    return encodings


def CTDT(fastas):
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity':        'HQRKNED',
        'polarizability':  'KMHFRYW',
        'charge':          'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess':   'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
    'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
    'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
    'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = ['name', 'sequence']
    for p in property:
        for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
            header.append(p + '.' + tr)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name, sequence]
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        for p in property:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code = code + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
        encodings.append(code)
    return encodings


def CTD(fastas):
    CTDC_1 = CTDC(fastas)
    CTDT_1 = CTDT(fastas)
    CTDD_1 = CTDD(fastas)

    CTD_all = []

    for i in range(len(CTDC_1)):
        CTD_all_1 = CTDC_1[i] + CTDT_1[i][2:] + CTDD_1[i][2:]
        CTD_all.append(CTD_all_1)

    return CTD_all


def save_csv(data, file_name):
    name = data[0]
    data1 = pd.DataFrame(columns=name, data=data[1:])

    data1.to_csv(file_name, encoding='gbk')
    return


file_1 = 'data/LLPS+ (high LLPS-propensity sequences).fasta'
file_2 = 'data/LLPS- (low LLPS-propensity sequences).fasta'
file_3 = 'data/Subset of PDB used for training.fasta'


fastas_1 = read_protein_sequences(file_1)
fastas_2 = read_protein_sequences(file_2)
fastas_3 = read_protein_sequences(file_3)


file_name_1 = '七种序列特征/LLPS+_AAC.csv'
file_name_2 = '七种序列特征/LLPS+_DPC.csv'
file_name_3 = '七种序列特征/LLPS+_TPC.csv'
file_name_4 = '七种序列特征/LLPS+_GAAC.csv'
file_name_5 = '七种序列特征/LLPS+_GDPC.csv'
file_name_6 = '七种序列特征/LLPS+_GTPC.csv'
file_name_7 = '七种序列特征/LLPS+_CTD.csv'


x1 = AAC(fastas_1)
x2 = DPC(fastas_1)
x3 = TPC(fastas_1)
x4 = GAAC(fastas_1)
x5 = GDPC(fastas_1)
x6 = GTPC(fastas_1)
x7 = CTD(fastas_1)


save_csv(x1, file_name_1)
save_csv(x2, file_name_2)
save_csv(x3, file_name_3)
save_csv(x4, file_name_4)
save_csv(x5, file_name_5)
save_csv(x6, file_name_6)
save_csv(x7, file_name_7)


file_name_1 = '七种序列特征/LLPS-_AAC.csv'
file_name_2 = '七种序列特征/LLPS-_DPC.csv'
file_name_3 = '七种序列特征/LLPS-_TPC.csv'
file_name_4 = '七种序列特征/LLPS-_GAAC.csv'
file_name_5 = '七种序列特征/LLPS-_GDPC.csv'
file_name_6 = '七种序列特征/LLPS-_GTPC.csv'
file_name_7 = '七种序列特征/LLPS-_CTD.csv'


x1 = AAC(fastas_2)
x2 = DPC(fastas_2)
x3 = TPC(fastas_2)
x4 = GAAC(fastas_2)
x5 = GDPC(fastas_2)
x6 = GTPC(fastas_2)
x7 = CTD(fastas_2)


save_csv(x1, file_name_1)
save_csv(x2, file_name_2)
save_csv(x3, file_name_3)
save_csv(x4, file_name_4)
save_csv(x5, file_name_5)
save_csv(x6, file_name_6)
save_csv(x7, file_name_7)


file_name_1 = '七种序列特征/PDB_training_AAC.csv'
file_name_2 = '七种序列特征/PDB_training_DPC.csv'
file_name_3 = '七种序列特征/PDB_training_TPC.csv'
file_name_4 = '七种序列特征/PDB_training_GAAC.csv'
file_name_5 = '七种序列特征/PDB_training_GDPC.csv'
file_name_6 = '七种序列特征/PDB_training_GTPC.csv'
file_name_7 = '七种序列特征/PDB_training_CTD.csv'


x1 = AAC(fastas_3)
x2 = DPC(fastas_3)
x3 = TPC(fastas_3)
x4 = GAAC(fastas_3)
x5 = GDPC(fastas_3)
x6 = GTPC(fastas_3)
x7 = CTD(fastas_3)


save_csv(x1, file_name_1)
save_csv(x2, file_name_2)
save_csv(x3, file_name_3)
save_csv(x4, file_name_4)
save_csv(x5, file_name_5)
save_csv(x6, file_name_6)
save_csv(x7, file_name_7)