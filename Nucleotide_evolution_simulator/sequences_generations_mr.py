import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
from mutation_rate import *
import csv
import os
#import seaborn as sns
#import matplotlib.pyplot as plt



# ajouter background_gBGC sur l'ensemble de la séquence
def newsequence_equilibrium(seq_length,SS,TSS, mr_TSS, mr_EIB,mut_rate_matrix,nuc_combNr,residual_tss,residual_eib,TSS_args, EIB_args):
    """
    *** GOAL : Generate a random nucleotide sequence of given length based on equilibrium nucleotide probabilities per position.

    *** Input :
        -   seq_length (int): Length of the sequence to generate.
        -   ATGC_mut (dict): Base mutation rates for transitions between nucleotides (intergenic baseline), hypothesis : strand symmetry in intergenic region.
        -   base_mut_rate (float): Base mutation rate multiplier.
        -   .
        -   GtoA, CtoT (list), parametres to generate the CpG distribution (strand assymetry).
        -   TSS (int): Start position for TSS region.
        -   SS (int): Start position for SS region.

    *** Output : bytes: Byte array representing the sequence where each position is a nucleotide ('A', 'C', 'G', or 'T').
    """
    def compute_equilibrium_probs(seq_length, mr_TSS, mr_EIB,mut_rate_matrix,SS,TSS, nuc_combNr,residual_tss,residual_eib,TSS_args, EIB_args):
        # --- Implement specific mutation rates corrections ---

        """ # CpG effect distribution
        CtoT_dist = CpG_mutation_rate(CtoT, half_length)
        GtoA_dist = CpG_mutation_rate(GtoA, half_length)"""
        
        """overall_rates = [mut_rate_matrix[i,0] for i in range(12)] #mutation rates are constant (for now)
        mut_types = nuc_combNr.keys()
        # build DataFrame
        df = pd.DataFrame({
            'Mutation_Type': mut_types,
            'Mutation_Rate': overall_rates
        })

        # plot
        plt.figure(figsize=(8,5))
        sns.barplot(data=df, x='Mutation_Type', y='Mutation_Rate', hue='Mutation_Type', palette="Set2")
        plt.ylim(0, 0.0014)
        plt.title('Overall intergenic baseline mutation rates')
        plt.show()"""
        if TSS_args : 
            #TSS effect
            f_mean_tss = np.array(TSS_mutation_rate(mr_TSS))
        if EIB_args :
            #EIB effect
            f_mean_eib =np.array(EIB_mutation_rate(mr_EIB))
        #print(len(scaled_profile)

        for mut_name in nuc_combNr.keys():
            dist = np.ones(seq_length)
            mut = nuc_combNr[mut_name]
            ## -- Exon / Intron specific effects --
            if TSS_args :
                if mut_name in residual_tss.keys():
                    mut_residual_tss = residual_tss[mut_name]
                    scaled_profile = residual_TSS(f_mean_tss,mut_name,mut_residual_tss)
                    dist[TSS-500:TSS+SS] =  scaled_profile[0:500+SS] #tss effect
                else:
                    dist[TSS-500:TSS+SS] =  f_mean_tss[TSS-500:TSS+SS] #tss effect
                
                mut_rate_matrix[mut, TSS-500:TSS+SS] *= dist[TSS-500:TSS+SS]
                """plt.plot(np.arange(0,5000),mut_rate_matrix[mut,:], label=f"{mut_name}")
                plt.title(f"{mut_name} initialisation of sequence after TSS")
                plt.show()"""
            if EIB_args : 
                if mut_name in residual_eib.keys():
                    mut_residual_eib = residual_eib[mut_name]
                    scaled_profile = residual_EIB(f_mean_eib,mut_name,mut_residual_eib)
                    dist[TSS+SS:] = scaled_profile[0:seq_length-(TSS+SS)] #eib effect
                else:
                    dist[TSS+SS:] = f_mean_eib[0:seq_length-(TSS+SS)] #eib effect
                #print("shape mut_rate_matrix:",mut_rate_matrix.shape)
                mut_rate_matrix[mut, TSS+SS:seq_length] *= dist[TSS+SS:seq_length]
                """plt.plot(np.arange(0,5000),mut_rate_matrix[mut,:], label=f"{mut_name}")
                plt.title(f"{mut_name} initialisation of sequence after TSS and EIB modifications")
                plt.show()"""

        """plt.figure(figsize=(12, 6))  # create one figure

        for mut_name in nuc_combNr:
            mut_idx = nuc_combNr[mut_name]
            plt.plot(np.arange(seq_length), mut_rate_matrix[mut_idx, :], label=f"{mut_name}")

        plt.title("mutation rate per mutation type hotspot gBGC before the rounds")
        plt.legend()
        plt.show()
        plt.close()"""

        # --- BUILD Q matrices and compute equilibrium per position ---

        equilibrium_probs = []

        for i in range(0,5000):
            # Build transition rates at position i
            QAG =  mut_rate_matrix[nuc_combNr['AG'],i]
            QAC =  mut_rate_matrix[nuc_combNr['AC'],i]
            QAT =  mut_rate_matrix[nuc_combNr['AT'],i]

            QCT =  mut_rate_matrix[nuc_combNr['CT'],i]
            QCA =  mut_rate_matrix[nuc_combNr['CA'],i]
            QCG =  mut_rate_matrix[nuc_combNr['CG'],i]

            QGA =  mut_rate_matrix[nuc_combNr['GA'],i]
            QGT =  mut_rate_matrix[nuc_combNr['GT'],i]
            QGC =  mut_rate_matrix[nuc_combNr['GC'],i]

            QTC =  mut_rate_matrix[nuc_combNr['TC'],i]
            QTG =  mut_rate_matrix[nuc_combNr['TG'],i]
            QTA =  mut_rate_matrix[nuc_combNr['TA'],i]

            QAA = -(QAG + QAC + QAT)
            QCC = -(QCT + QCA + QCG)
            QGG = -(QGA + QGT + QGC)
            QTT = -(QTC + QTG + QTA)

            Q = np.array([
                [QAA, QAC, QAG, QAT],
                [QCA, QCC, QCG, QCT],
                [QGA, QGC, QGG, QGT],
                [QTA, QTC, QTG, QTT]
            ])

            # Solve for equilibrium distribution
            Q_aug = np.vstack([Q.T, np.ones(Q.shape[0])])
            b = np.zeros(Q.shape[0]+1)
            b[-1] = 1

            pi = np.linalg.lstsq(Q_aug, b, rcond=None)[0]
            equilibrium_probs.append(pi)

        equilibrium_probs = np.array(equilibrium_probs)
        return equilibrium_probs

    # Step 1: Compute equilibrium probabilities for each region
    eq_probs = compute_equilibrium_probs(seq_length, mr_TSS, mr_EIB,mut_rate_matrix,SS,TSS,nuc_combNr,residual_tss, residual_eib,TSS_args, EIB_args)
    print(eq_probs.shape)

    # Step 2: Generate the sequence
    nucleotides = np.frombuffer(b'ACGT', dtype=np.uint8)
    sequence = np.zeros(seq_length, dtype=np.uint8)

    # Apply default equilibrium probabilities to each position of the sequence
    #print("example of eq_probs:",eq_probs[5])
    sequence = np.zeros(seq_length, dtype=np.uint8)
    for i in range(seq_length):
        sequence[i] = np.random.choice(nucleotides, p=eq_probs[i])
    result = bytearray(sequence) #.tobytes()

    #plot the proportion of ACGT in every region 
    # Count A, C, G, T
    #region_names = ['intergenic', 'exon', 'intron']
    #regions = [sequence[0:TSS], sequence[TSS:TSS+SS],sequence[TSS+SS:]]
    
    #print("length intergenic",len(sequence[0:TSS]))
    #print("length exon",len(sequence[TSS:TSS+SS]))
    #print("length intronic",len(sequence[TSS+SS:]))

    """for name, seq in zip(region_names, regions):
        #plt.figure(figsize=(10,8))
        composition={}
        tab=[]

        counts = np.bincount(seq, minlength=256)  # since sequence is uint8
        counts = counts[ord('A')], counts[ord('C')], counts[ord('G')], counts[ord('T')]
        total = sum(counts)
        freqs = [c/total for c in counts]

        nuc_dico = {'A':0,'C':1,'G':2,'T':3}
        for nuc in nuc_dico.keys():
            composition[nuc] = tab
            tab.append(freqs[nuc_dico[nuc]])
            #df = pd.DataFrame.from_dict(composition,orient='index')
            #sns.barplot(df, palette="Set2")
        
        #print(f"frequency of nuc in {name}:", tab)
        path_freq = os.path.join(absolute_path, "nuc_content_freq_initialisation.csv")
        with open(path_freq, 'a', newline='') as file:
            writer = csv.writer(file)
            header = ["region", "A", "C", "G", "T"]
            #print(header)
            
            if  os.path.getsize(path_freq) == 0:
                writer.writerow(header)
                writer.writerow([f"{name}", tab[0] ,tab[1],tab[2],tab[3]])
            else: 
                writer.writerow([f"{name}", tab[0] ,tab[1],tab[2],tab[3]])
    
        plt.title(f'{name} nucleotide frquency')
        plt.show()
        plt.close()"""

    return result



def choose_base_with_weight(w):
    return random.choices(b'ACGT', weights=w, k=1)[0]


def consensus_splice_bytes():
    motif = bytearray()
    gtat_splicesite = random.random()
    if (gtat_splicesite <= (53.58/(53.58+45.10))):
        # if we fall into the first consensus sequence
        motif.append(choose_base_with_weight([32, 33, 18, 17]))
        motif.append(choose_base_with_weight([50, 15, 15, 20]))
        motif.append(choose_base_with_weight([15, 5, 70, 10]))
        motif += b"GT"
        motif.append(choose_base_with_weight([1, 0, 1, 0]))
        motif += b"AG"
        motif.append(choose_base_with_weight([16, 15, 20, 50]))
        motif.append(choose_base_with_weight([32, 22, 28, 18]))
        motif.append(choose_base_with_weight([22, 27, 23, 28]))
        motif.append(choose_base_with_weight([22, 28, 23, 27]))
        motif.append(choose_base_with_weight([22, 23, 27, 28]))
    else:
        # otherwise we consider the other one
        motif.append(choose_base_with_weight([38, 40, 15, 7]))
        motif.append(choose_base_with_weight([80, 5, 7, 10]))
        motif.append(choose_base_with_weight([4, 0, 94, 2]))
        motif += b"GT"
        motif.append(choose_base_with_weight([72, 6, 16, 6]))
        motif.append(choose_base_with_weight([25, 17, 25, 23]))
        motif.append(choose_base_with_weight([20, 12, 50, 18]))
        motif.append(choose_base_with_weight([20, 15, 20, 45]))
        motif.append(choose_base_with_weight([27, 18, 30, 25]))
        motif.append(choose_base_with_weight([23, 25, 22, 30]))
        motif.append(choose_base_with_weight([22, 25, 23, 30]))
        motif.append(choose_base_with_weight([22, 23, 25, 30]))
    return motif


def conserve_inplace(sequence, TSS, SS, CDS_start, splice_motif, TSSnuc):
    """KN: in place version of conserve above which "fixes" the input sequence given
           as argument. We only overwrite the ranges that are modified in place.
    """
    sequence[TSS-1] = ord(TSSnuc)  # assumes len = 1
    idx = TSS + SS - 4
    sequence[idx: idx + len(splice_motif)] = splice_motif
    # If CDS_start is provided, replace the codon at that position with the start motif
    if CDS_start is not None:
        sequence[TSS+ CDS_start:TSS+ CDS_start + 3] = b'ATG'



def complement(nucleotide):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complement_dict.get(nucleotide, nucleotide)


def pdf_trans(x, epsilon, delta):
    """  mathematical transformation and evaluation of a probability density function (PDF)
    Parameters:
    x: The variable for which the PDF is being calculated.
    epsilon: A parameter that shifts the transformed variable.
    delta: A scaling parameter that affects both the transformation and the PDF.
    The output is a scaled and normalized version of the standard normal PDF applied to the transformed variable.
    """
    transformed = delta * np.arcsinh(x) - epsilon
    pdf_value = stats.norm.pdf(np.sinh(transformed))
    scaling_factor = delta * np.cosh(transformed)
    normalization = np.sqrt(1 + np.power(x, 2))

    return pdf_value * scaling_factor / normalization


def SS_generate(a=56.1, c=0.1966, loc=1.876, scale=0.0000003479):
    """Generate random samples from a generalized gamma distribution.

    Parameters:
    a (float): Shape parameter.
    c (float): Shape parameter.
    loc (float): Location parameter.
    scale (float): Scale parameter.

    Returns:
    float: Random sample from the distribution.
    see Movassat et al., RNA 2019, curve fitted above (exon 1 size)
    SS = 5'SS = splice site, or exon/intron boundary
    """

    # Generate random sample from the distribution
    return stats.gengamma.rvs(a, c, loc, scale)


def ATG_generate(a=1.795,  c=1.048, loc=1.792,  scale=136.3):
    """Generate random samples from a generalized gamma distribution.

    Parameters:
    a (float): Shape parameter.
    c (float): Shape parameter.
    loc (float): Location parameter.
    scale (float): Scale parameter.
    Values: see Leppek et al., Nat Rev MCB 2018, curve fitted above (5'UTR size)

    Returns:
    float: Random sample from the distribution.
    """
    return stats.gengamma.rvs(a, c, loc, scale)


def prep_mutation_rate(mutation_rate, N,CpG_distribution): #, GtoA_CpG_distribution
    """ *** GOAL : gather per position the probability of each nuc to mutate to another nuc and implement the CpG effect
    by adding C in CpG and Gin CpG as two new nucleotides
    *** Input : 
        -   mutation_rate : numpy array 12 lines * 5000 columns that contains per position the probability of mutation
        -   N : dictionnary containing the 12 mutations possible, and associates to each mr a number
        -   CtoT_CpG_distribution : a list, that is 5000 long that contains the probability of C>T per position
        -   GtoA_CpG_distribution : a list that is 5000 long that contains the probability of G>A per position

    *** Output : 
        -   output : a list, that contains for each position for each nuc the probability of that nuc to change for another
        for position 1 for nuc A for instance output = [ {65: [np.float64(0.00018842111999999999), np.float64(0.0009188912354145894), np.float64(0.00014402303999999996)], ...}]"""
    
    output = []
    bases_num = [65, 67, 71, 84]
    for idx in range(len(mutation_rate)): # parcourir les positions de la sequence
        #preparation of mutation_rate to 
        d = {}
        output.append(d)
        for b1 in bases_num: #parcourir les 4 nt
            tab = []
            d[b1] = tab #ex print 1er tour : d = {65 : []}
            for b2 in bases_num: # parcourir 2e nt
                if b1 == b2:
                    continue 
                tab.append(mutation_rate[idx, N[chr(b1) + chr(b2)]]) 


        fCpG = CpG_distribution[idx] # creation de la correction
        d[68] = [d[67][0], d[67][1], d[67][2] * fCpG] #ajouter au dico 2 nouvelles clees
        d[72] = [d[71][0] * fCpG, d[71][1], d[71][2]]
        
        
        dmax = 0
        for tab in d.values(): #ex tab = [0.1, 0.2, 0.3]
            tab[1] += tab[0] #tab [1] = 0.2 + 0.1 = 0.3
            tab[2] += tab[1] #tab [2] = 0.3 + 0.3 = 0.6
            #tab = [0.1, 0.3, 0.6]
            if tab[2] > dmax:
                dmax = tab[2]
        d['max'] = dmax 

    return output #liste de dictionnaires 



TARGET_BASES = {65: b'CGT', 67: b'AGT', 71: b'ACT', 84: b'ACG'}
# KN: a dictionnary to convert an index 0 1 2 to the corresponding
#    code. For instance for base C (67) the lowest target base is A, then G, then T


#11 realiser les mutations sur la sequence directement
def mutate_inplace(byte_sequence: bytearray, mutation_rate): #
    """KN: mutate the given byte_sequence in place. We remark that modifications
           are local, so we can just overwrite the sequence with new bases or
           leave them unchanged. mutation_rate is an array returned by pre_mutation_rate
    """
    sequence_length = len(byte_sequence)

    # We temporarily add a fake "T" at the end. having a last T
    # will not influence the probability of the last real base and
    # allows us to not test idx < sequence_length - 1 every time to
    # check the next letter
    # we carefully compute sequence_length BEFORE, so the for loop
    # will stop before this fake last base.

    byte_sequence.append(84) #ajout de 1 T a la fin de la seq pour faciliter l'indexation

    # rolling window to access the array once per loop
    prev_base = None 
    cur_base = None
    next_base = byte_sequence[0] # element de la sequence ie : code ASCII
    idx = 0 #initialisation a la position 0 ie : 1ere base de la seq
    for idx in range(sequence_length):  # won't see the fake T
        prev_base = cur_base
        cur_base = next_base 
        next_base = byte_sequence[idx+1]

        mutation_chance = random.random()
        d = mutation_rate[idx]

        if mutation_chance >= d['max']:
            # shortcut
            continue

        # The second part of the sum is true if and only if:
        # the current base is G and the previous was C, yielding an index of 71 + int(True) = 72
        # the current base is C and the next is G, yielding an index of 67 + int(True) = 68

        key = cur_base + ((cur_base == 71 and prev_base == 67)
                          or (cur_base == 67 and next_base == 71)) #((cur_base == 71 and prev_base == 67) or (cur_base == 67 and next_base == 71)) renvoient 1 si elles sont vraies -> on obtient donc : key=68 ou key=72
        target_rates = d[key] #example with C
        if mutation_chance < target_rates[0]: #if mutation_chance < P(C->A) in CpG context (no change)
            byte_sequence[idx] = TARGET_BASES[cur_base][0]
        elif mutation_chance < target_rates[1]: #if mutation_chance < P(C->G) in CpG context (no change of mutation rate)
            byte_sequence[idx] = TARGET_BASES[cur_base][1]
        elif mutation_chance < target_rates[2]:#if mutation_chance < P(C->T) in CpG context -> it increases
            byte_sequence[idx] = TARGET_BASES[cur_base][2]# EX : TARGET_BASES[67]=b'AGT' ; TARGET_BASES[67][2] = 'AGT'[2] = 84 (ASCII) => T
        # else:
        #   don't mutate

    byte_sequence.pop()  # clean-up
