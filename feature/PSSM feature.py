import sys
import numpy as np
import pandas as pd
import math
import re
import fileinput
from os import listdir
from os.path import isfile, join


def readToMatrix(input_matrix):
    #print "start to read PSSM matrix"
    PSSM = []
    p = re.compile(r'-*[0-9]+')
    for line, strin in enumerate(input_matrix):
        if line > 2:
            str_vec = []
            overall_vec = strin.split()
            #print len(overall_vec)
            if len(overall_vec) == 0:
                break
            str_vec.extend(overall_vec[1])
            if(len(overall_vec) < 44):
                print ("There is a mistake in the pssm file")
                print ("Try to correct it")
                for cur_str in overall_vec[2:]:
                    str_vec.extend(p.findall(cur_str))
                    if(len(str_vec) >= 21):
                        if(len(str_vec)) >21:
                            #print len(str_vec)
                            #print str_vec
                            #print overall_vec
                            #print "Exit with an error"
                            exit(1)
                        break;
                print ("Done")
            else:
                str_vec = strin.split()[1:42]
            if len(str_vec) == 0:
                break
            #str_vec_positive=map(int, str_vec[1:])
            PSSM.append(str_vec)
    fileinput.close()
    #print "finish to read PSSM matrix"
    PSSM = np.array(PSSM)
    return PSSM

def average(matrixSum, seqLen):
    # average the summary of rows
    matrix_array = np.array(matrixSum)
    matrix_array = np.divide(matrix_array, seqLen)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [(np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1], )))]
    return matrix_average

def handleRows(PSSM, SWITCH, COUNT):
    '''
    if SWITCH=0, we filter no element.
    if SWITCH=1, we filter all the negative elements.
    if SWITCH=2, we filter all the negative and positive elements greater than expected.
    '''
    '''
    if COUNT=20, we generate a 20-dimension vector.
    if COUNT=400, we generate a 400-dimension vector.
    '''
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"

    matrix_final = [ [0.0] * 20 ] * int((COUNT/20))
    matrix_final=np.array(matrix_final)
    seq_cn = 0

    PSSM_shape=np.shape(PSSM)
    for i in range(PSSM_shape[0]):
        seq_cn += 1
        str_vec=PSSM[i]
        str_vec_positive=map(int, str_vec[1:21])
        str_vec_positive=np.array(str_vec_positive)
        if SWITCH==1:
            str_vec_positive[str_vec_positive<0]=0
        elif SWITCH==2:
            str_vec_positive[str_vec_positive<0]=0
            str_vec_positive[str_vec_positive>7]=0
        #print "str_vec_positive="
        #print str_vec_positive
        if COUNT==20:
            matrix_final[0]=map(sum, zip(str_vec_positive, matrix_final[0]))
        elif COUNT==400:
            matrix_final[Amino_vec.index(str_vec[0])] = map(sum, zip(str_vec_positive, matrix_final[Amino_vec.index(str_vec[0])]))

        #print "matrix_final="
        #print matrix_final

    return matrix_final

def aac_pssm(input_matrix):
    #print "start aac_pssm function"
    SWITCH = 0
    COUNT = 20
    seq_cn=float(np.shape(input_matrix)[0])
    aac_pssm_matrix=handleRows(input_matrix,SWITCH,COUNT)
    aac_pssm_matrix=np.array(aac_pssm_matrix)
    aac_pssm_vector=average(aac_pssm_matrix,seq_cn)
    #print "end aac_pssm function"
    return aac_pssm_vector

def ab_pssm(input_matrix):
    #print "start ab_pssm function"
    seq_cn=np.shape(input_matrix)[0]
    BLOCK=seq_cn/20
    #print BLOCK
    matrix_final=[]
    for i in range(19):
        tmp=input_matrix[i*BLOCK:(i+1)*BLOCK]
        #print tmp
        matrix_final.append(aac_pssm(tmp)[0])
    tmp=input_matrix[19*BLOCK:]
    #print tmp
    matrix_final.append(aac_pssm(tmp)[0])
    ab_pssm_matrix=average(matrix_final,1.0)
    #print "finish ab_pssm function"
    return ab_pssm_matrix

def pssm_composition(input_matrix):
    #print "start pssm_composition function"
    SWITCH = 0
    COUNT = 400
    seq_cn=float(np.shape(input_matrix)[0])
    pssm_composition_matrix=handleRows(input_matrix, SWITCH, COUNT)
    pssm_composition_vector=average(pssm_composition_matrix,seq_cn)
    #print "end pssm_composition function"
    return pssm_composition_vector

def rpm_pssm(input_matrix):
    #print "start rpm_pssm function"
    SWITCH = 1
    COUNT = 400
    seq_cn=float(np.shape(input_matrix)[0])
    rpm_pssm_matrix=handleRows(input_matrix,SWITCH,COUNT)
    rpm_pssm_vector=average(rpm_pssm_matrix,seq_cn)
    #print "end rpm_pssm function"
    return rpm_pssm_vector

def s_fpssm(input_matrix):
    #print "start s_fpssm function"
    SWITCH = 2
    COUNT = 400
    seq_cn = 1
    s_fpssm_matrix=handleRows(input_matrix, SWITCH, COUNT)
    s_fpssm_matrix = np.array(s_fpssm_matrix)
    s_fpssm_matrix_shape = np.shape(s_fpssm_matrix)
    matrix_average = [(np.reshape(s_fpssm_matrix, (s_fpssm_matrix_shape[0] * s_fpssm_matrix_shape[1], )))]
    #print "end s_fpssm function"
    return matrix_average


inputFile = 'input/example.fasta' # input a file in fasta format
check_head = re.compile(r'\>')

smplist = []
smpcnt = 0
a = fileinput.input(inputFile)
for line, strin in enumerate(fileinput.input(inputFile)):
    if not check_head.match(strin):
        smplist.append(strin.strip())
        smpcnt += 1

pssmdir = 'input/pssm_files' # specify the directory of pssm files

onlyfiles = [ f for f in listdir(pssmdir) if isfile(join(pssmdir,f)) ]

fastaDict = {}
for fi in onlyfiles:
    cntnt = ''
    pssmContentMatrix=readToMatrix(fileinput.input(pssmdir+'/'+fi))
    pssmContentMatrix=np.array(pssmContentMatrix)
    sequence=pssmContentMatrix[:,0]
    seqLength=len(sequence)
    for i in range(seqLength):
        cntnt+=sequence[i]
    if cntnt in fastaDict:
        #print strin
        continue
    fastaDict[cntnt] = fi


finalist = []
for smp in smplist:
    #print "smp="+smp
    #print "fastaDict[smp]="+fastaDict[smp]
    finalist.append(pssmdir+'/'+fastaDict[smp])

x_ab_pssm = []
x_pssm_composition = []
x_rpm_pssm = []
x_s_fpssm = []

for fi in finalist:
    input_matrix=fileinput.input(fi)
    input_matrix=readToMatrix(input_matrix)

    x1 = ab_pssm(input_matrix)
    x2 = pssm_composition(input_matrix)
    x3 = rpm_pssm(input_matrix)
    x4 = s_fpssm(input_matrix)

    x_ab_pssm.append(x1)
    x_pssm_composition.append(x2)
    x_rpm_pssm.append(x3)
    x_s_fpssm.append(x4)

x_ab_pssm = np.array(x_ab_pssm)
x_pssm_composition = np.array(x_pssm_composition)
x_rpm_pssm = np.array(x_rpm_pssm)
x_s_fpssm = np.array(x_s_fpssm)

x_ab_pssm = x_ab_pssm.reshape([smpcnt,400])
x_pssm_composition = x_pssm_composition.reshape([smpcnt,400])
x_rpm_pssm = x_rpm_pssm.reshape([smpcnt,400])
x_s_fpssm = x_s_fpssm.reshape([smpcnt,400])

x_ab_pssm = pd.DataFrame(x_ab_pssm)
x_ab_pssm.to_csv('ab_pssm.csv', index=False, header=False, escapechar=',')
x_pssm_composition = pd.DataFrame(x_pssm_composition)
x_pssm_composition.to_csv('pssm_composition.csv', index=False, header=False, escapechar=',')
x_rpm_pssm = pd.DataFrame(x_rpm_pssm)
x_rpm_pssm.to_csv('rpm_pssm.csv', index=False, header=False, escapechar=',')
x_s_fpssm = pd.DataFrame(x_s_fpssm)
x_s_fpssm.to_csv('s_fpssm.csv', index=False, header=False, escapechar=',')

print('good')