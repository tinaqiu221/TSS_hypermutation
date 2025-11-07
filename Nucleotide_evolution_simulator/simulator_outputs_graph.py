import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from sequences_generations_mr import *


def loading_sequences(filename,bins):
    """ Load the sequences created into a dictionary
    Keys are the gene name
    Content is the sequence
    """
    input_file = open(filename+".txt", 'r')
    dict_seq = dict()
    gene_name = input_file.readline().strip()
    seq_line = input_file.readline().strip()
    while gene_name[0] != '*':
        dict_seq[gene_name] = seq_line     
        gene_name = input_file.readline().strip()
        seq_line = input_file.readline() .strip()
    input_file.close()
    return dict_seq

def extract_SS(name):
    pattern = r'SS:(\d+)'
    match = re.search(pattern, name)
    if match:
        number = int(match.group(1))
        return number
    else:
        return None

def freq_per_bin_exon(seq,gene_name, bins, nuc,TSS):
    """ Given a sequence and a number of bin, compute some freq of a nuc
    seq : the sequence that is from dict_seq
    bins : number of bins
    nuc: which nucleotide I want to compute the frequency
    """   
    freq_bin=[]
    SS = extract_SS(gene_name)
    # you extract this info from the artificial gene that evolves
    #tells you the position of the SS that you keep in the gene name of the results files like gBGC_evo_5000bp_300rounds_500seqs.txt
    #TSS = len(seq)//2
    start_analysis = TSS-(2*SS)
    end_analysis = TSS+(2*SS)
    bin_size = ((end_analysis - start_analysis) / ((bins*4)+1)) 
    if bin_size>2:
        beginning = start_analysis
        for a in range(bins*4+1):
            end = int(round(((a+1)*bin_size),0)) + start_analysis
            next_beginning = end+1
            if nuc =="GC":
                freq_bin.append((seq[beginning:end].count('G')+seq[beginning:end].count('C'))/(end-beginning))
            elif nuc =="CpG":
                freq_bin.append(seq[beginning:end].count('CG')/(end-beginning-1))
            else:
                freq_bin.append(seq[beginning:end].count(nuc)/(end-beginning))
            beginning = next_beginning
    return freq_bin
     

def freq_per_bin(seq, bins, nuc,seq_length):
    """ Given a sequence and a number of bin, compute some freq of a nuc
    seq : the sequence that is from dict_seq
    bins : number of bins
    nuc: which nucleotide I want to compute the frequency
    """
    freq_bin=[]
    bin_size = seq_length/bins
    if bin_size>2:
        beginning = 0  
        for a in range(bins):
            end = int(round(((a+1)*bin_size),0))
            next_beginning = end+1
            if nuc =="GC":
                freq_bin.append((seq[beginning:end].count('G')+seq[beginning:end].count('C'))/(end-beginning))
            elif nuc =="CpG":
                freq_bin.append(seq[beginning:end].count('CG')/(end-beginning-1))
            else:
                freq_bin.append(seq[beginning:end].count(nuc)/(end-beginning))
            beginning = next_beginning
    return freq_bin

""" Analyze nucleotide content """

def analyze_nucleotide_content(analysis_type, T_iter, nrep,add_B_iter,bg_B_iter, bins,seq_length):
    filename = f"{analysis_type}_B_TSS_Tcurrent{T_iter}_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_iter{nrep}"
    # Load the input file
    dict_seq = loading_sequences(filename, bins)
    dict_nuc_summaries = dict()
    dict_exon_summaries = dict()

    for nuc in ["A", "C", "G", "T", "CpG","GC"]:
        dict_nuc_summaries[nuc] = dict()
        dict_exon_summaries[nuc] = dict()

    for n_seq in dict_seq.keys():
        cur_seq = dict_seq[n_seq]
        for nuc in ["A", "C", "G", "T", "CpG","GC"]:
            dict_nuc_summaries[nuc][n_seq] = freq_per_bin(cur_seq, bins, nuc, seq_length)
            dict_exon_summaries[nuc][n_seq] = freq_per_bin_exon(cur_seq, n_seq, 10, nuc, TSS=2500) #exon_bins = 10

    for nuc in ["A", "C", "G", "T", "CpG","GC"]:
        # Write nucleotide summaries
        with open(f"All_{nuc}_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt", 'w') as file:
            for key, value in dict_nuc_summaries[nuc].items():
                file.write('%s\t%s\n' % (key, ", ".join(map(str, value))))
        with open(f"Exon_All_{nuc}_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt", 'w') as file:
            for key, value in dict_exon_summaries[nuc].items():
                file.write('%s\t%s\n' % (key, ", ".join(map(str, value)))) 

    """ Write sequences summaries """
def output_nuc_content(summary_nuc, nrep, seq_length, ngenes,cpg, current_rounds, background_B, longterm_rounds, additional_B): #, GtoA_cpg
    # Compute mean per position
    nuc_axis = summary_nuc.mean(axis=0)
    
    # Prepare the data row for appending
    data_row = [f"{value}" for value in nuc_axis] #faire un dico plutot car ceci est une chaine de caractères
    return data_row

def summary_analysis(analysis_type,evo_stats,current_rounds,longterm_rounds,T_iter,add_B_iter,bg_B_iter,background_B,additional_B,cpg,nrep,bins): #GtoA_cpg,
    input_files = {
        'A': f"All_A_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'T': f"All_T_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'G': f"All_G_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'C': f"All_C_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'GC': f"All_GC_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'CpG': f"All_CpG_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt"
    }

    nuc_summaries = {key: [] for key in input_files.keys()}

    seq_length = evo_stats['seq_length']
    ngenes = evo_stats['ngenes']


    for key, filename in input_files.items():
        with open(filename, 'r') as input_file:
            lines = input_file.readlines()

        data = {}
        for line in lines:
            index, values = line.split('\t')
            values = list(map(float, values.split(',')))
            data[index] = values
        
        df = pd.DataFrame.from_dict(data, orient='index')
        nuc_summaries[key] = df

    for key, summary in nuc_summaries.items():
        file_name = f"{key}_summary_{analysis_type}_analysis.csv"

        # Generate the data row using the refactored function
        data_row = output_nuc_content(summary, nrep, seq_length, ngenes, cpg, current_rounds, background_B, longterm_rounds, additional_B) #, GtoA_cpg,
        with open(file_name,mode='w', newline="") as file:
            file.write(",".join(map(str, data_row)) + "\n")
        # Append the row to the summary file safely
        #safe_append_summary_ACGT(file_name, data_row)
    

    A_real= [0.27472,0.27410,0.27545,0.27451,0.27511,0.27469,0.27375,0.27257,0.27145,0.27066,0.26971,0.26751,0.26425,0.26042,0.25554,0.24944,0.24176,0.23111,0.21499,0.18239,0.16343,0.16753,0.17012,0.18278,0.19788,0.21220,0.22431,0.23319,0.23953,0.24443,0.24885,0.25133,0.25347,0.25535,0.25743,0.25876,0.25948,0.25955,0.25987,0.25956,0.26076]
    T_real = [0.27489,0.27457,0.27350,0.27412,0.27297,0.27350,0.27272,0.27168,0.27012,0.26852,0.26677,0.26371,0.26004,0.25736,0.25115,0.24365,0.23339,0.22278,0.20583,0.17985,0.17286,0.17649,0.18658,0.20624,0.22409,0.23859,0.25021,0.25950,0.26670,0.27095,0.27548,0.27620,0.27878,0.27953,0.28162,0.28277,0.28391,0.28480,0.28494,0.28485,0.28554]
    G_real = [0.22446,0.22439,0.22470,0.22495,0.22541,0.22497,0.22560,0.22662,0.22789,0.22949,0.23054,0.23415,0.23731,0.24031,0.24568,0.25237,0.26065,0.27011,0.28464,0.32015,0.34684,0.34566,0.34003,0.31928,0.29992,0.28467,0.27319,0.26341,0.25624,0.25133,0.24670,0.24518,0.24260,0.24105,0.23838,0.23763,0.23640,0.23531,0.23387,0.23520,0.23353]
    C_real = [0.22593,0.22694,0.22635,0.22642,0.22652,0.22684,0.22793,0.22913,0.23055,0.23133,0.23299,0.23464,0.23839,0.24191,0.24763,0.25455,0.26421,0.27600,0.29455,0.31760,0.31687,0.31032,0.30326,0.29170,0.27811,0.26454,0.25229,0.24389,0.23754,0.23330,0.22897,0.22728,0.22515,0.22407,0.22257,0.22083,0.22021,0.22034,0.22131,0.22039,0.22017]
    GC_real = [0.45039,0.45133,0.45105,0.45137,0.45193,0.45181,0.45353,0.45575,0.45844,0.46082,0.46352,0.46878,0.47571,0.48223,0.49330,0.50691,0.52486,0.54611,0.57919,0.63775,0.66371,0.65598,0.64330,0.61098,0.57803,0.54921,0.52548,0.50730,0.49377,0.48462,0.47567,0.47247,0.46775,0.46512,0.46095,0.45846,0.45661,0.45566,0.45519,0.45559,0.45371]
    CpG_real = [0.01450,0.01455,0.01470,0.01487,0.01518,0.01537,0.01594,0.01658,0.01738,0.01842,0.01936,0.02120,0.02325,0.02568,0.02936,0.03377,0.03959,0.04640,0.05753,0.07824,0.09273,0.08261,0.07612,0.06482,0.05429,0.04437,0.03722,0.03186,0.02776,0.02490,0.02235,0.02102,0.01953,0.01838,0.01773,0.01669,0.01610,0.01564,0.01531,0.01515,0.01491]
    X_real = [-2500, -2375, -2250, -2125, -2000, -1875, -1750, -1625,-1500, -1375, -1250, -1125, -1000,  -875,  -750,  -625,-500,  -375,  -250,  -125,     0,   125,   250,   375,500,   625,   750,   875,  1000,  1125,  1250,  1375,1500,  1625,  1750,  1875,  2000,  2125,  2250,  2375,2500]

    # Load the data from the CSV file
    filename= {
    "A_file" : f"A_summary_{analysis_type}_analysis.csv",
    "C_file" : f"C_summary_{analysis_type}_analysis.csv",
    "G_file" : f"G_summary_{analysis_type}_analysis.csv",
    "T_file" : f"T_summary_{analysis_type}_analysis.csv",
    "GC_file" : f"GC_summary_{analysis_type}_analysis.csv",
    "CpG_file" : f"CpG_summary_{analysis_type}_analysis.csv"
    }

    data_A = pd.read_csv(filename["A_file"], delimiter=",", header = None)
    data_C = pd.read_csv(filename["C_file"], delimiter="," , header = None)
    data_G = pd.read_csv(filename["G_file"], delimiter="," , header = None)
    data_T = pd.read_csv(filename["T_file"], delimiter="," , header = None)
    data_GC = pd.read_csv(filename["GC_file"], delimiter=",", header = None)
    data_CpG = pd.read_csv(filename["CpG_file"], delimiter=",", header = None)

    X_axis = range(0, bins)

    # Calculate X_mapped
    bin_width = (max(X_real) - min(X_real)) / bins
    X_mapped = [min(X_real) + bin_width * (i + 0.5) for i in X_axis] 


    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(X_mapped, data_A.iloc[0,:], linestyle='dashed', marker='', color='black' , label='A Sim')
    ax.plot(X_mapped, data_T.iloc[0,:], linestyle='dashed', marker='', color='green' , label='T Sim')
    ax.plot(X_mapped, data_G.iloc[0,:], linestyle='dashed', marker='', color='blue' , label='G Sim')
    ax.plot(X_mapped, data_C.iloc[0,:], linestyle='dashed', marker='', color='red' , label='C Sim')

    ax.plot(X_mapped, A_real, linestyle='-', marker='', color='black' , label='A Real')
    ax.plot(X_mapped, T_real, linestyle='-', marker='', color='green' , label='T Real')
    ax.plot(X_mapped, G_real, linestyle='-', marker='', color='blue' , label='G Real')
    ax.plot(X_mapped, C_real, linestyle='-', marker='', color='red' , label='C Real')


    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Position",fontsize=18 )
    ax.set_ylabel("Nucleotide content",fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim(0.10,0.4)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f"{analysis_type}_Simulation_vs_Observation_Tcurrent{T_iter}_additional_B{add_B_iter}.png")  # Save the figure
    plt.close(fig)

    # Second plot: CpG
    fig, ax = plt.subplots(figsize=(10, 8))  # Create a new figure
    ax.plot(X_mapped, data_CpG.iloc[0,:], linestyle='dashed', marker='', color='grey', label='CpG Sim')
    ax.plot(X_mapped, CpG_real, linestyle='-', marker='', color='grey', label='CpG Real')
        
    # Add vertical lines

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Position",fontsize=18 )
    ax.set_ylabel("Nucleotide content",fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim(0, 0.15)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f"{analysis_type}_CpG_Simulation_vs_Observation_Tcurrent{T_iter}_additional_B{add_B_iter}.png")  # Save the figure
    plt.close(fig)

    # third plot: GC
    fig, ax = plt.subplots(figsize=(10, 8)) # Create a new figure
    ax.plot(X_mapped, data_GC.iloc[0,:], linestyle='dashed', marker='', color='purple' , label='GC Sim')
    ax.plot(X_mapped, GC_real, linestyle='-', marker='', color='purple' , label='GC Real')
        
    # Add vertical lines

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Position",fontsize=18 )
    ax.set_ylabel("Nucleotide content",fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim(0.4, 0.75)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f"{analysis_type}_GC_Simulation_vs_Observation_Tcurrent{T_iter}_additional_B{add_B_iter}.png")  # Save the figure
    plt.close(fig)

    for filename in input_files.values():
        if os.path.exists(filename):  # Avoid FileNotFoundError
            os.remove(filename)


def summary_analysis_exon(analysis_type,evo_stats,current_rounds,longterm_rounds,T_iter,add_B_iter,bg_B_iter,background_B,additional_B,cpg,nrep):
    input_files = {
        'A': f"Exon_All_A_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'T': f"Exon_All_T_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'G': f"Exon_All_G_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'C': f"Exon_All_C_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'GC': f"Exon_All_GC_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt",
        'CpG': f"Exon_All_CpG_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_Tcurrent{T_iter}_iter{nrep}_{analysis_type}analysis.txt"
    }

    nuc_summaries = {key: [] for key in input_files.keys()}

    seq_length = evo_stats['seq_length']
    ngenes = evo_stats['ngenes']


    for key, filename in input_files.items():
        with open(filename, 'r') as input_file:
            lines = input_file.readlines()

        data = {}
        for line in lines:
            index, values = line.split('\t')
            values = list(map(float, values.split(',')))
            data[index] = values
        
        df = pd.DataFrame.from_dict(data, orient='index')
        nuc_summaries[key] = df

    for key, summary in nuc_summaries.items():
        file_name = f"Exon_{key}_summary_{analysis_type}_analysis.csv"

        # Generate the data row using the refactored function
        data_row = output_nuc_content(summary, nrep, seq_length, ngenes, cpg, current_rounds, background_B, longterm_rounds, additional_B) #, GtoA_cpg,
        with open(file_name,mode='w', newline="") as file:
            file.write(",".join(map(str, data_row)) + "\n")
        # Append the row to the summary file safely
        #safe_append_summary_ACGT(file_name, data_row)
    

    A_real= [0.22827,0.22563,0.22550,0.22253,0.22009,0.21596,0.21444,0.21437,0.21001,0.20579,0.20261,0.19870,0.19534,0.19169,0.18703,0.18405,0.17538,0.17413,0.16871,0.15715,0.17526,0.16143,0.16321,0.16341,0.16285,0.16683,0.16804,0.17143,0.17688,0.18578,0.21383,0.19140,0.16389,0.16661,0.17077,0.17686,0.17964,0.18267,0.18586,0.18890,0.19124]
    T_real = [0.22362,0.22222,0.21612,0.21850,0.21511,0.21331,0.21015,0.20753,0.20559,0.20312,0.19918,0.19653,0.19230,0.19010,0.18679,0.18138,0.17754,0.17467,0.17199,0.17099,0.19368,0.17638,0.16768,0.16935,0.17134,0.17188,0.17332,0.17472,0.17679,0.17679,0.17839,0.19553,0.19379,0.19893,0.20391,0.20566,0.21232,0.21484,0.21508,0.21872,0.22199]
    G_real = [0.27180,0.27482,0.27768,0.27724,0.28169,0.28418,0.28661,0.29074,0.29325,0.29786,0.29993,0.30450,0.30903,0.31313,0.31703,0.32558,0.33144,0.33383,0.33917,0.34622,0.31277,0.34433,0.35402,0.34862,0.34473,0.34056,0.33585,0.33384,0.32819,0.32182,0.37773,0.37205,0.35955,0.34902,0.34417,0.33790,0.33174,0.32556,0.32293,0.31672,0.31510]
    C_real = [0.27632,0.27732,0.28070,0.28173,0.28311,0.28655,0.28881,0.28736,0.29116,0.29323,0.29828,0.30027,0.30333,0.30507,0.30915,0.30900,0.31563,0.31737,0.32014,0.32563,0.31828,0.31785,0.31509,0.31862,0.32108,0.32073,0.32280,0.32001,0.31814,0.31562,0.23005,0.24101,0.28276,0.28543,0.28115,0.27959,0.27630,0.27693,0.27614,0.27566,0.27167]
    GC_real = [0.54812,0.55215,0.55838,0.55897,0.56481,0.57073,0.57542,0.57810,0.58440,0.59108,0.59822,0.60477,0.61236,0.61821,0.62618,0.63458,0.64707,0.65120,0.65931,0.67185,0.63105,0.66219,0.66911,0.66723,0.66580,0.66129,0.65864,0.65385,0.64633,0.63744,0.60778,0.61306,0.64232,0.63445,0.62532,0.61748,0.60804,0.60249,0.59906,0.59238,0.58677]
    CpG_real = [0.04794,0.04906,0.05055,0.05166,0.05305,0.05483,0.05621,0.05691,0.05956,0.06162,0.06401,0.06607,0.07012,0.07210,0.07487,0.07843,0.08433,0.08686,0.09110,0.09695,0.08233,0.09481,0.09348,0.09115,0.08940,0.08786,0.08604,0.08402,0.08349,0.08017,0.06250,0.06674,0.07244,0.07033,0.06555,0.06433,0.06071,0.05941,0.05863,0.05683,0.05461]
    X_real = np.arange(0,41).tolist()  # Assuming X_real is a list of positions from 0 to 40

    # Load the data from the CSV file
    filename= {
    "A_file" : f"Exon_A_summary_{analysis_type}_analysis.csv",
    "C_file" : f"Exon_C_summary_{analysis_type}_analysis.csv",
    "G_file" : f"Exon_G_summary_{analysis_type}_analysis.csv",
    "T_file" : f"Exon_T_summary_{analysis_type}_analysis.csv",
    "GC_file" : f"Exon_GC_summary_{analysis_type}_analysis.csv",
    "CpG_file" : f"Exon_CpG_summary_{analysis_type}_analysis.csv"
    }

    data_A = pd.read_csv(filename["A_file"], delimiter=",", header = None)
    data_C = pd.read_csv(filename["C_file"], delimiter="," , header = None)
    data_G = pd.read_csv(filename["G_file"], delimiter="," , header = None)
    data_T = pd.read_csv(filename["T_file"], delimiter="," , header = None)
    data_GC = pd.read_csv(filename["GC_file"], delimiter=",", header = None)
    data_CpG = pd.read_csv(filename["CpG_file"], delimiter=",", header = None)


    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(X_real, data_A.iloc[0,:], linestyle='dashed', marker='', color='black' , label='A Sim')
    ax.plot(X_real, data_T.iloc[0,:], linestyle='dashed', marker='', color='green' , label='T Sim')
    ax.plot(X_real, data_G.iloc[0,:], linestyle='dashed', marker='', color='blue' , label='G Sim')
    ax.plot(X_real, data_C.iloc[0,:], linestyle='dashed', marker='', color='red' , label='C Sim')

    ax.plot(X_real, A_real, linestyle='-', marker='', color='black' , label='A Real')
    ax.plot(X_real, T_real, linestyle='-', marker='', color='green' , label='T Real')
    ax.plot(X_real, G_real, linestyle='-', marker='', color='blue' , label='G Real')
    ax.plot(X_real, C_real, linestyle='-', marker='', color='red' , label='C Real')


    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Bins",fontsize=18 )
    ax.set_ylabel("Nucleotide content",fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim(0.10,0.4)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f"Exon_{analysis_type}_Simulation_vs_Observation_Tcurrent{T_iter}_additional_B{add_B_iter}.png")  # Save the figure
    plt.close(fig)

    # Second plot: CpG
    fig, ax = plt.subplots(figsize=(10, 8))  # Create a new figure
    ax.plot(X_real, data_CpG.iloc[0,:], linestyle='dashed', marker='', color='grey', label='CpG Sim')
    ax.plot(X_real, CpG_real, linestyle='-', marker='', color='grey', label='CpG Real')
        
    # Add vertical lines

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Bins",fontsize=18 )
    ax.set_ylabel("Nucleotide content",fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim(0, 0.15)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f"Exon_{analysis_type}_CpG_Simulation_vs_Observation_Tcurrent{T_iter}_additional_B{add_B_iter}.png")  # Save the figure
    plt.close(fig)

    # third plot: GC
    fig, ax = plt.subplots(figsize=(10, 8)) # Create a new figure
    ax.plot(X_real, data_GC.iloc[0,:], linestyle='dashed', marker='', color='purple' , label='GC Sim')
    ax.plot(X_real, GC_real, linestyle='-', marker='', color='purple' , label='GC Real')
        
    # Add vertical lines

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Bins",fontsize=18 )
    ax.set_ylabel("Nucleotide content",fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim(0.4, 0.75)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f"Exon_{analysis_type}_GC_Simulation_vs_Observation_Tcurrent{T_iter}_additional_B{add_B_iter}.png")  # Save the figure
    plt.close(fig)

    for filename in input_files.values():
        if os.path.exists(filename):  # Avoid FileNotFoundError
            os.remove(filename)
        # Remove specific files
    file_path1 = f"active_B_TSS_Tcurrent{T_iter}_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_iter{nrep}.txt"
    file_path2 = f"background_B_TSS_Tcurrent{T_iter}_additional_B{add_B_iter}_bg_B_iter{bg_B_iter}_iter{nrep}.txt"
    os.remove(file_path1)
    os.remove(file_path2)

"""  Study of the initialisation of the sequences  """
def create_data_frame(file_path):
    
    #create a data frame from csv
    df = pd.read_csv(file_path)
    
    # -- filter columns --
    df_intergenic = df[df["region"] == "intergenic"]
    df_exonic = df[df["region"] == "exon"]
    df_intronic = df[df["region"] == "intron"]
    
    # -- write output files --
    #df_intergenic.to_csv("C:/Users/eh263/OneDrive/Documents/Emilie/Research/GC_Simulator/Results/intergenic_init.csv")
    #df_exonic.to_csv("C:/Users/eh263/OneDrive/Documents/Emilie/Research/GC_Simulator/Results/exonic_init.csv")
    #df_intronic.to_csv("C:/Users/eh263/OneDrive/Documents/Emilie/Research/GC_Simulator/Results/intronic_init.csv")
    
    return df_intergenic, df_exonic, df_intronic

def mean_df(df): 
    #calculate the mean value for A,C,G,T content
    A_mean = df.loc[:,"A"].mean()
    C_mean = df.loc[:,"C"].mean() 
    G_mean = df.loc[:,"G"].mean() 
    T_mean = df.loc[:,"T"].mean()
    # build dataframe
    df_mean = pd.DataFrame({
        "A": [A_mean],
        "C": [C_mean],
        "G": [G_mean],
        "T": [T_mean]
    })
    return df_mean

"""def plot_nuc_content(df, region, ngenes):
    #Plot the mean nucleotide content for each region 
    plt.figure(figsize=(10,8))
    ax = sns.barplot(data=df, palette="Set2")
    plt.title(f"Nucleotide content at initialisation of sequences \n(over {ngenes} sequences) in {region} region")
    plt.xlabel("Nucleotide")
    plt.ylabel("Mean content")

    # Annotate each bar with its mean value
    for i, p in enumerate(ax.patches):
        height = p.get_height()
        ax.text(
            p.get_x() + p.get_width() / 2,  # x position
            height,                         # y position
            f"{height:.3f}",                # text (rounded to 3 decimals)
            ha="center", va="bottom", fontsize=12
        )

    plt.savefig(f"nucleotide_content_initialisation_{region}.png")
    plt.close()"""

""" def plot_nuc_content_init(file_path, ngenes):
    # -- Creation of the data frame per region --
    df_intergenic, df_exonic, df_intronic = create_data_frame(file_path)
    region_dfs = {
        "intergenic": df_intergenic,
        "exonic": df_exonic,
        "intronic": df_intronic
    }
    for region in ["intergenic", "exonic", "intronic"]:
        df = region_dfs[region].drop(columns="region")
        plot_nuc_content(mean_df(df), region, ngenes)
    os.remove(file_path)  # Remove the input file after processing """