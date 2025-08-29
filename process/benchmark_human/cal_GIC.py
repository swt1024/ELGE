import os
import sys
import re
import csv
import time
import math
import string
import pandas as pd
from tqdm import tqdm


# record all triplet combinations
BASES = ['a', 't', 'c', 'g']
TRIPLETS = []
for B1 in BASES:
    for B2 in BASES:
        for B3 in BASES:
            TRIPLETS.append(B1+B2+B3)

mers = {} #construct a list of triplets that are initialized to 0

# obtain the corresponding relationship between lncRNA and transcript
def readlnc_transID(filename):
    lnc_trans = pd.read_csv(filename, sep='\t', header=None, names=['lncRNA_ID', 'transcript_id'])
    #Create the mapping relationship between genes and transcripts
    lnc_trans_mapping = lnc_trans.groupby('lncRNA_ID')['transcript_id'].aggregate(list).to_dict()

    return lnc_trans_mapping

# Get the secondary structure minimum free energy of each transcript of lncRNA
def readlnc_mfe(filename):
    trans_mfe = pd.read_csv(filename)
    #Create the mapping relationship between transcripts and MFE
    trans_mfe_mapping = pd.Series(trans_mfe['MFE'].values, index=trans_mfe['transcript_id']).to_dict()

    return trans_mfe_mapping

# Get transcript sequence
def readfasta(filename):
    f = open(filename,'r')
    res = {} #Preserve transcripts and their sequence
    for line in f:
        if line.startswith('>'):
            line = line.strip()
            ID = line.split('>',2)[1]
            res[ID] = ''
        else:
            res[ID] += line
    return res

#sliding window
def slidingWindow(seq, l, win, step=1): 
    length = l
    mod = divmod((length-win), step)[1]
    if (win >= length):
        return seq
    else:
        start = 0
        end = win
        fragments = []
        while (len(seq[start:end]) == win):
            fragments.append(seq[start:end])
            start += step
            end += step
        if (mod > 0):
            fragments.append(seq[(length-win):])
        return fragments
 
#calculate the frequency of the triplet
def stat3mer(seq, l):   
    freq = {}
    for item in TRIPLETS:
        mers[item] = 0
    num3mer = float(l-2)
    all3mer = slidingWindow(seq, l, win=3, step=1)
    for i in set(TRIPLETS):
        mers[i] = all3mer.count(i)
    for triplet in TRIPLETS:
        freq[triplet] = mers[triplet]/num3mer
    return freq


if __name__ == '__main__':

    lnc_trans_path = './filtered_lnc_trans.txt'  
    MFE_path = './trans_MFE.csv'
    seq_path = './transcript_sequences.fasta'  

    lnc_transID = readlnc_transID(lnc_trans_path)
    trans_mfe = readlnc_mfe(MFE_path)
    trans_seq = readfasta(seq_path)

    # parameter Settings of GIC
    LRMODEL_FEATURES_7 = ['intercept', 'length', 'mfe/L', 'cga', 'gcg', 'tcg', 'acg', 'tca']
    LRMODEL_COEFS_7 = [0.7417, 2.612e-04, 4.295, 48.66, 15.64, 76.23, -1.113, -60.29]
    LRMODEL_7 = dict(zip(LRMODEL_FEATURES_7, LRMODEL_COEFS_7))

    #calculated the eigenvalues of each transcript
    features_trans = {}  #save the features of each transcript
    for id, seq in tqdm(trans_seq.items(), desc="Processing Transcripts"):
        if id not in trans_mfe:
            continue
        feature = {}
        seq = seq.replace('\n','')
        seq = seq.lower()
        length = len(seq)
        feature['trans'] = id
        feature['length'] = length
        feature['mfe/L'] = trans_mfe[id]/length
        freq = stat3mer(seq, length)
        for item in TRIPLETS:
            feature[item] = freq[item]
        # tmp = LRMODEL_7['intercept']+sum([feature[_]*LRMODEL_7[_] for _ in LRMODEL_FEATURES_7[1:]])
        # gic = math.exp(tmp)/(math.exp(tmp)+1)
        # feature['GIC_7'] = gic
        features_trans[id] = feature
    
    #calculated the eigenvalues of each lncRNA
    lncRNA_GIC_score = {}  ##save the GIC score of each lncRNA
    for lnc,transcripts in lnc_transID.items():
        feature = {} #save the features of lncRNA gene
        feature['length'] = 0 # initialize the value of each feature to 0
        feature['mfe/L'] = 0
        for item in TRIPLETS:
            feature[item] = 0
    
        for trans in transcripts:
            if trans not in features_trans:
                feature = {}
                break
            else:
                tranfeature = features_trans[trans]
                feature['length'] += tranfeature['length']
                feature['mfe/L']  += tranfeature['mfe/L']
                for item in TRIPLETS:
                    feature[item] += tranfeature[item]
        if len(feature)>0:
            feature['lncRNA_ID'] = lnc
            feature['length'] = feature['length']/len(transcripts)
            feature['mfe/L'] = feature['mfe/L']/len(transcripts)
            for item in TRIPLETS:
                feature[item] =  feature[item]/len(transcripts)
            tmp = LRMODEL_7['intercept']+sum([feature[_]*LRMODEL_7[_] for _ in LRMODEL_FEATURES_7[1:]])
            gic = math.exp(tmp)/(math.exp(tmp)+1)
            feature['GIC_score'] = gic
            lncRNA_GIC_score[lnc] = feature['GIC_score']

    # sort all LncRNAs according to GIC scores
    lncRNA_GIC_score = sorted(lncRNA_GIC_score.items(),key=lambda d:d[1],reverse=False)
    
    # Convert the lncRNA GIC scores to a DataFrame
    lncRNA_scores_df_sorted = pd.DataFrame(lncRNA_GIC_score, columns=['lncRNA_ID', 'GIC_score'])

    # Write the sorted results to a CSV file
    sorted_GIC_path = 'sorted_GIC_score.csv'
    lncRNA_scores_df_sorted.to_csv(sorted_GIC_path, index=False)

    print("GIC scores have been calculated and saved successfully.")
    print(f"Output file is located at: {sorted_GIC_path}")