import time
from sequences_generations_mr import *
from simulator_outputs_graph import *
from mutation_rate import *
from main_graph import *
import numpy as np
import argparse
import ast
import os


start = time.time()

def main(background_B,bg_B_iter,additional_B,add_B_iter,seq_length,longterm_rounds,current_rounds,T_iter,ngenes, bins_all,cpg,nrep,mr_TSS, mr_EIB,TSS_args, EIB_args):
    ''' REPRENDRE LE TEXTE EXPLICATIF
    
    
    --nseq NSEQ           an integer for the number of ngenes (required)
  --bg_B           a float for the background B (default : 0.8)
  --add_B   a float for B at hotspot (previously called additional_B, default : 1.3)
  --additional_B
                        an integer for number of rounds with active BGC (default : 300)
  --current_rounds
                        an integer for the number of round with the background BGC (default: 100)
  --length LENGTH       an integer for the sequence length (default : 5000)
  --bins BINS [BINS ...]
                        an integer for the number of bins  of the region around TSS and for only the exonic region (default: [40, 10])
  --cpg CPG [CPG ...]   a list of integers for CpG (default: . These numbers correspond to a CpG peak 160bp before TSS (the curve has a sigma of 1300, a delta of 1). The peak CpG rate is half and the baseline CpG rate is 13 times the baseline mutation rate
'''
    #np.random.seed(42) #generation de nb aleatoires
    #random.seed(42)

    control_sum = 0
    
    #seq_length = 5000  # 5000
    half_length = seq_length // 2
    TSS = half_length 
    base_mut_rate = 0.001 #taux de mutation de base
    #ngenes = 500
    # valueFanny
    #longterm_rounds = 300
    #current_rounds = 100
    UI5 = 0

    bins = bins_all[0]
    exon_bins = bins_all[1]

    """ basic mutation base rates """
    ATGC_mut = {
        'AT_rate': 0.386, 'GC_rate': 0.61447, 'AC': 0.186, 'AG': 0.635,
        'AT': 0.179, 'GA': 0.560, 'GC': 0.205, 'GT': 0.236}

    # Interconversion tools to manipulate basic mutation rates
    nuc_combN = {'AG': 0, 'AC': 1, 'AT': 2, 'GA': 3, 'GT': 4, 'GC': 5} #dico
    nuc_combL = ['AG', 'AC', 'AT', 'GA', 'GT', 'GC'] #liste

    """ Stats """
    evo_stats = {
        'longterm_rounds': longterm_rounds ,
        'current_rounds': current_rounds,
        'seq_length': seq_length,
        'ngenes': ngenes,
        'base_mut_rate': base_mut_rate,
        'time': time
    }

    # Initialize mut_rate_init using numpy for efficiency
    mut_rate_init = np.zeros((6, seq_length), dtype=float)
    # Mapping for the corresponding rate keys
    rate_keys = {'A': 'AT_rate', 'G': 'GC_rate'}
    # Precompute the rates for each combination
    rates = np.array([
        [base_mut_rate *
         ATGC_mut[rate_keys[mut[0]]] * ATGC_mut[mut]] for mut in nuc_combL])#calcul taux mutation pour chaque paire en ajustant
    # Use broadcasting to fill the array
    mut_rate_zeros = np.zeros((len(rates), seq_length), dtype=float)
    mut_rate_init = rates + mut_rate_zeros
    
    #overall_rates_baseline = [mut_rate_init[i,0] for i in range(6)] #mutation rates are constant (for now)

    if background_B != 0 :
        """ Add background longterm_rounds, formula from Glemin et al 2015 Genome Research """
        mut_rate_init = background_gBGC(background_B,nuc_combN, mut_rate_init) #bg_gBGC modification is a modification by a constant proper to each mutation type

        #overall_rates_gBGC = [mut_rate_init[i,0] for i in range(6)] #mutation rates are constant (for now)
        #graph_baseline(overall_rates_baseline,overall_rates_gBGC,nuc_combN)

    """ Evolution of the sequence under additional_B at TSS """
    """ peak center: peak; Increase in fixation of weak(AT) to strong(GC): WS_factor; Decrease in fixation of strong(GC) to weak(AT): SW_factor; peak width (one sigma): sigma; distribution: calculation of normal curve """

    # Initialize mut_rate_init_gBGC as a numpy array
    mut_rate_init_gBGC = np.copy(mut_rate_init)

    #Compute additional gBGC at TSS (hotspot)
    if additional_B!=0 :
        gBGC_TSS = {'GC_peak_rel_TSS': peak,
        'GC_sigma': width,
        'GC_delta': 1,
        'GC_epsilon': 0} #best aprox place for GC_peak_rel_TSS + GC_sigma
        
        mut_rate_init_gBGC = hotspot_gBGC(additional_B, gBGC_TSS, nuc_combN,length=half_length,mut_rate_gBGC=mut_rate_init_gBGC)
    
        """for mut in range(mut_rate_init.shape[0]):
            mut_name = nuc_combL[mut]
            positions = np.arange(len(mut_rate_init[mut]))
            plt.plot(positions, mut_rate_init[mut,:], label=f'{mut_name} baseline', color = "green")
            plt.plot(positions, mut_rate_init_gBGC[mut,:], label=f'{mut_name} hotspot gBGC', color = "blue")
            plt.legend()
            plt.title(f"mut rate with and without gbgc {mut_name}")
            plt.show()"""

    """ CpG rate - gaussian curves  """

    CpG_distribution = CpG_mutation_rate(cpg, half_length)
    #GtoA_CpG_distribution = CpG_mutation_rate(GtoA_cpg, half_length)

    """plt.plot(np.linspace(1,5000,5000),CtoT_CpG_distribution[0:], label = "ctot distribution")
    plt.legend()
    plt.title("CtoT distribution")
    plt.show()"""


    """TSS correction factors"""
    """residual_tss = { "AG":[0.409,174,507,-0.37,-0.453,532,6.77], "CG":[0.256,42.9,495,-0.230,-0.299,527,25.1],
                    "CT":[0.431,244,264,-0.347,-0.270,543,35], "GA":[-0.164,147,120,0.159,0.317,511,33.9],
                    "GC":[-0.197,333,246,0.146,0.194,462,95.6], "TC":[-0.239,121,182,0.234,0.153,528,34.2], "TG":[-0.136,380,78.1,0.0160,-0.231,479,33.9]}"""
    residual_tss = ast.literal_eval(args.res_tss)



    """EIB correction factors"""
    """residual_eib = { "CA":[-0.137,0.00136,0.175], "CG":[-0.336,0.000675,0.0496],"CT":[-0.375,0.0025,-0.221], 
                    "GA":[0.415,0.0006,-0.0853],"GC":[-0.140,0.000931,0.22],"GT":[-0.366,0.00162,0.00367],
                    "TA":[0,0.0853], "TC":[0.305,0.001,0.103], "TG":[0,0.0188]}"""
    residual_eib = ast.literal_eval(args.res_eib)

    """ conservation of 5'UTR and CDS """
    # wherever you are at UTR or CDS, for the moment we reduce mut rate by 30 or 80%
    # instead of reducing the fixation rate we want a constant mut rate but sln to occur and decrease the fiation rate (ie. not exacttly a fix 80% reduction, simplifying assumtion)
    conserve_UTR = 0.7
    conserve_CDS = 0.2
    """ define mut_rate_rep with interconversion tools """
    # mut_rate_rep = np.ndarray(shape=(12, seq_length), dtype=float)
    # mut_rate_rep_gBGC = np.ndarray(shape=(12, seq_length), dtype=float)

    # KN: reshape before the loop (matrix is doubled so that it goes from 6 to 12, gBGC was applied before doubling so gBGC is symmetric)
    mut_rate_init = np.concatenate((mut_rate_init,) * 2)
    mut_rate_init_gBGC = np.concatenate((mut_rate_init_gBGC,) * 2)

    nuc_combNr = dict([('AG', 0), ('AC', 1), ('AT', 2), ('GA', 3), ('GT', 4), (
        'GC', 5), ('TC', 6), ('TG', 7), ('TA', 8), ('CT', 9), ('CA', 10), ('CG', 11)])
    nuc_combLr = ['AG', 'AC', 'AT', 'GA', 'GT',
                  'GC', 'TC', 'TG', 'TA', 'CT', 'CA', 'CG']


    """ set directory """
    # Relative path to the directory
    relative_path = "../Results"
    # Get the absolute path
    absolute_path = os.path.abspath(relative_path)
    os.chdir(absolute_path)

    """ evolve each gene as a replicate """
    longterm_rounds_output_file = open("active_B_TSS_Tcurrent"+str(T_iter)+"_additional_B"+str(add_B_iter)+"_bg_B_iter"+str(bg_B_iter)+"_iter"+str(nrep)+".txt", 'w')
    current_rounds_output_file = open("background_B_TSS_Tcurrent"+str(T_iter)+"_additional_B"+str(add_B_iter)+"_bg_B_iter"+str(bg_B_iter)+"_iter"+str(nrep)+".txt", 'w')


    for replicate in range(ngenes):
        """ reset mutation rate for each replicate """
        mut_rate_rep = np.copy(mut_rate_init)
        mut_rate_rep_gBGC = np.copy(mut_rate_init_gBGC)
        
        """ define 5'Splice Site (SS) position based on human genomic data Movassat et al., 2019 """
        SS = int(SS_generate())
        while (SS > (seq_length/4 - 10)) or (SS < 40):
            SS = int(SS_generate())
    
        
        #GC_content = 0.42
        #sequence = newsequence(seq_length, GC_content)
        sequence = newsequence_equilibrium(seq_length,SS,TSS, mr_TSS, mr_EIB,mut_rate_rep,nuc_combNr,residual_tss,residual_eib,TSS_args, EIB_args)

        """ define 5'Splice Site (SS) motif based on Sibley et al., 2016 """
        splice_motif = consensus_splice_bytes()

        """ define the nucleotide at the TSS"""  # to check  do the analysis with ensembl
        TSSnuc = chr(choose_base_with_weight([0.474044195,0.122530106,0.377586812, 0.025838887]))

        """ define ATG position based on Leppek et al., Nat Rev MCB 2018 """
        CDS_start = int(ATG_generate())

        if (CDS_start > (SS-4)) or (CDS_start < 3):
            CDS_start = None
            UI5 = UI5+1
           # based on the SS and the position of ATG we modify the mutation rate accordingly, for each replicate. bc of polak paper
        
        """ decrease/increase mutations around TSS """
        # Create final scaling profile: elementwise product of mr_TSS_distribution and base factor
        f_mean = np.array(TSS_mutation_rate(mr_TSS))
        
        """plt.plot(np.linspace(2000,3000,1000), scaled_profile[2000:3000], label= "average curve fit", color='red' )
        plt.title("average curve fit mr_TSS")
        plt.legend()
        plt.show()"""
        if TSS_args :
            for mut_name in nuc_combLr:
                mut_idx = nuc_combNr[mut_name]
                if mut_name in residual_tss.keys():
                    #specific corrections for certain mutation types
                    mut_residual = residual_tss[mut_name]
                    scaled_profile = residual_TSS(f_mean,mut_name,mut_residual)
                    # Apply per-position mutation rate adjustments (preserving the shape) 
                    #mut_rate_rep[mut_idx, TSS-500:TSS+SS] *= scaled_profile[0:500+SS]  #already implemented in newsequence_equilibrium()
                    mut_rate_rep_gBGC[mut_idx, TSS-500:TSS+SS] *= scaled_profile[0:500+SS] #not implemented in newsequence_equilibrium()
                    """plt.plot(np.arange(TSS-500,TSS+SS), scaled_profile[0:500+SS])
                    plt.title(f"{mut_name} residual TSS")
                    plt.show()"""
                else:
                    mut_rate_rep_gBGC[mut_idx, TSS-500:TSS+SS] *= f_mean[TSS-500:TSS+SS] #not implemented in newsequence_equilibrium()
                    #mut_rate_rep[mut_idx, TSS-500:TSS+SS] *= f_mean[TSS-500:TSS+SS] #already implemented in newsequence_equilibrium()
            
        if EIB_args :
            """ EIB specific mutation rate implementation """
            f_mean = np.array(EIB_mutation_rate(mr_EIB))
            """x= np.arange(TSS+SS,seq_length)
            plt.plot(x, scaled_profile[0:seq_length-(TSS+SS)], label= "EIB average curve fit", color='purple' )
            plt.title("average curve fit EIB")
            plt.legend()
            plt.show()"""
            
            #plt.figure(figsize=(12, 6))
            for mut_name in nuc_combLr:
                mut_idx = nuc_combNr[mut_name]
                if mut_name in residual_eib.keys():
                    #print("mut_name correction eib: ", mut_name)
                    mut_residual = residual_eib[mut_name]
                    #print("mut_residual:", mut_residual)

                    scaled_profile = residual_EIB(f_mean, mut_name, mut_residual)

                    #mut_rate_rep[mut_idx, TSS+SS :seq_length] *= scaled_profile[0:seq_length-(TSS+SS)]
                    mut_rate_rep_gBGC[mut_idx, TSS+SS :seq_length] *= scaled_profile[0:seq_length-(TSS+SS)] 
                    '''plt.figure(figsize=(12, 6))
                    plt.plot(np.arange(seq_length), mut_rate_rep[mut_idx, :], label=f"{mut_name}", color= "blue")
                    plt.plot(np.arange(seq_length), mut_rate_rep_gBGC[mut_idx, :], label=f"{mut_name} hotspot B", color='green')
                    plt.plot(np.arange(TSS+SS, seq_length), scaled_profile[0:seq_length-(TSS+SS)], label=f"{mut_name} EIB correction", color='purple')
                    plt.legend()
                    plt.title(f"eib specific correction {mut_name}")
                    plt.show()'''
                
                else: 
                    mut_rate_rep_gBGC[mut_idx, TSS+SS :seq_length] *= f_mean[0:seq_length-(TSS+SS)]
                    #mut_rate_rep[mut_idx, TSS+SS :seq_length] *= f_mean[0:seq_length-(TSS+SS)]
            """overall_rates_gBGC = [mut_rate_rep[i,1700] for i in range(12)] #mutation rates are constant (for now)
            mut_types = nuc_combNr.keys()
            # build DataFrame
            df = pd.DataFrame({
                'Mutation_Type': mut_types,
                'Mutation_Rate': overall_rates_gBGC
            })

            # plot
            plt.figure(figsize=(8,5))
            sns.barplot(data=df, x='Mutation_Type', y='Mutation_Rate', hue='Mutation_Type', palette="Set2")
            plt.ylim(0, 0.0014)
            plt.title('Overall intergenic baseline mutation rates after hotspot gBGC correction +TSS+EIB')
            plt.show()

            #check before the longterm rounds whether there are changes in the matrix
            plt.figure(figsize=(12, 6))  # create one figure

            for mut_name in nuc_combLr:
                mut_idx = nuc_combNr[mut_name]
                plt.plot(np.arange(seq_length), mut_rate_rep[mut_idx, :], label=f"{mut_name}")

            plt.title("mutation rate per mutation type hotspot gBGC before the rounds")
            plt.legend()
            plt.show()"""

        """ decrease mutations in 5'UTR """
        if CDS_start == None:
            for nuc_change in range(12):
                mut_rate_rep_gBGC[nuc_change, TSS:TSS+SS] *= conserve_UTR
                mut_rate_rep[nuc_change, TSS:TSS+SS] *= conserve_UTR
        else:
            for nuc_change in range(12):
                mut_rate_rep_gBGC[nuc_change, TSS:TSS+CDS_start] *= conserve_UTR
                mut_rate_rep[nuc_change, TSS:TSS+CDS_start] *= conserve_UTR

            """ decrease mutations in CDS """
            for nuc_change in range(12):
                mut_rate_rep_gBGC[nuc_change, TSS+CDS_start:TSS+SS] *= conserve_CDS
                mut_rate_rep[nuc_change, TSS+CDS_start:TSS+SS] *= conserve_CDS

        """#check before the current rounds whether there are changes in the matrix
        plt.figure(figsize=(12, 6))  # create one figure

        for mut_name in nuc_combLr:
            mut_idx = nuc_combNr[mut_name]
            plt.plot(np.arange(seq_length), mut_rate_rep[mut_idx, :], label=f"{mut_name}")

        plt.title("mutation rate per mutation type no hotspot gBGC before current rounds")
        plt.legend()
        plt.show()

        #check before the longterm rounds whether there are changes in the matrix
        plt.figure(figsize=(12, 6))  # create one figure

        for mut_name in nuc_combLr:
            mut_idx = nuc_combNr[mut_name]
            plt.plot(np.arange(seq_length), mut_rate_rep_gBGC[mut_idx, :], label=f"{mut_name}")

        plt.title("mutation rate per mutation type with hotspot gBGC before longterm rounds")
        plt.legend()
        plt.show()"""

        """ rounds of gBGC mutations with additional_B at TSS and background_B"""

        # mut_sequence = np.array(list(sequence))
        #consideration of CpG effect and additional_B effect on the mutation rates of the sequence (bytearray)
        if additional_B == 0:
            mut_rate_rep_noadd_opt = prep_mutation_rate(
                np.transpose(mut_rate_rep), nuc_combNr, CpG_distribution) #GtoA_CpG_distribution

            for generation in range(longterm_rounds):
                mutate_inplace(sequence, mut_rate_rep_noadd_opt) #evolution of the sequence under additional_B at TSS and background_B on the whole sequence

            conserve_inplace(sequence, TSS, SS,
                            CDS_start, splice_motif, TSSnuc)

            """ Save sequences """
            longterm_rounds_output_file.write(">Sequence_")
            longterm_rounds_output_file.write(str(replicate))
            longterm_rounds_output_file.write("_SS:")
            longterm_rounds_output_file.write(str(SS))
            longterm_rounds_output_file.write("\n")
            longterm_rounds_output_file.write(sequence.decode())
            longterm_rounds_output_file.write("\n")

        else:
            mut_rate_rep_gBGC_opt = prep_mutation_rate(
                np.transpose(mut_rate_rep_gBGC), nuc_combNr,CpG_distribution )#,GtoA_CpG_distribution
        
            for generation in range(longterm_rounds):
                mutate_inplace(sequence, mut_rate_rep_gBGC_opt) #evolution of the sequence under additional_B at TSS and background_B on the whole sequence
            conserve_inplace(sequence, TSS, SS,
                            CDS_start, splice_motif, TSSnuc)
            #print("CDS_start:",CDS_start)

            """ Save sequences """
            longterm_rounds_output_file.write(">Sequence_")
            longterm_rounds_output_file.write(str(replicate))
            longterm_rounds_output_file.write("_SS:")
            longterm_rounds_output_file.write(str(SS))
            longterm_rounds_output_file.write("\n")
            longterm_rounds_output_file.write(sequence.decode())
            longterm_rounds_output_file.write("\n")

        """ switch to the new phase: no longer additional_B at TSS """

        # Update mutation rate for the current phase
        mut_rate_rep_opt_current = prep_mutation_rate(
            np.transpose(mut_rate_rep), nuc_combNr,CpG_distribution)#GtoA_CpG_distribution #only background_B and CpG effect
        # Current evolution
        for generation in range(current_rounds):
            mutate_inplace(sequence, mut_rate_rep_opt_current) #modification of the sequence (that already evolved under additional_B) but continues to evolve under background_B

        # Conserve regions and save sequences
        conserve_inplace(sequence, TSS, SS, CDS_start, splice_motif, TSSnuc)

        """ Save sequences """
        current_rounds_output_file.write(">Sequence_")
        current_rounds_output_file.write(str(replicate))
        current_rounds_output_file.write("_SS:")
        current_rounds_output_file.write(str(SS))
        current_rounds_output_file.write("\n")
        current_rounds_output_file.write(sequence.decode())
        current_rounds_output_file.write("\n")

        """#check after the longterm rounds whether there are changes in the matrix
        plt.figure(figsize=(12, 6))  # create one figure

        for mut_name in nuc_combLr:
            mut_idx = nuc_combNr[mut_name]
            plt.plot(np.arange(seq_length), mut_rate_rep[mut_idx, :], label=f"{mut_name}")

        plt.title("mutation rate per mutation type no hotspot gBGC after current rounds")
        plt.legend()
        plt.show()

        #check after the longterm rounds whether there are changes in the matrix
        plt.figure(figsize=(12, 6))  # create one figure

        for mut_name in nuc_combLr:
            mut_idx = nuc_combNr[mut_name]
            plt.plot(np.arange(seq_length), mut_rate_rep_gBGC[mut_idx, :], label=f"{mut_name}")

        plt.title("mutation rate per mutation type with hotspot gBGC after longterm rounds")
        plt.legend()
        plt.show()"""

        """ Counter for every 50 """
        if replicate % 50 == 0:
            print(str(replicate)+" sequences done")
    """ Finalize sequence files """
    current_rounds_output_file.write("*")
    current_rounds_output_file.close()
    longterm_rounds_output_file.write("*")
    longterm_rounds_output_file.close()

    #analyze_nucleotide_content("active",T_iter, nrep,add_B_iter,bg_B_iter, bins,seq_length)
    analyze_nucleotide_content("background", T_iter, nrep,add_B_iter,bg_B_iter, bins,seq_length)


    print("Starting Summary analysis")

    #summary_analysis("active",evo_stats,current_rounds,longterm_rounds,T_iter,add_B_iter,bg_B_iter,background_B,additional_B,cpg,nrep,bins)
    summary_analysis("background",evo_stats,current_rounds,longterm_rounds,T_iter,add_B_iter,bg_B_iter,background_B,additional_B,cpg,nrep,bins)
    summary_analysis_exon("background",evo_stats,current_rounds,longterm_rounds,T_iter,add_B_iter,bg_B_iter,background_B,additional_B,cpg,nrep)
    # plot_nuc_content_init("nuc_content_freq_initialisation.csv")
    #file_path = os.path.join(absolute_path, "nuc_content_freq_initialisation.csv")
    #plot_nuc_content_init(file_path, ngenes)


    os.chdir("../..")

    end = time.time()
    print("Durée du calcul : " + str(end-start))
    print("Contrôle:", control_sum % (1 << 128))
    pass



if __name__ == "__main__":
    """ basic stats """
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser = argparse.ArgumentParser(
        description=(
            "A script to process nucleotides evolution around the transcription start sites (TSS).\n\n"
            "Example usage with default values and 100 sequences simulated :\n"
            "  python simulator_gc_tss.py --nseq 100 \n"
            "Example usage with all parameters :\n"
            "  python simulator_gc_tss.py --nseq 50 --background_B 0.8 --additional_B 1.3 --longterm_rounds 300 --current_rounds 100 --length 5000 --bins 40 10 --CtoT_cpg160 0.5 13 1300   --GtoA_cpg160 0.5 13 1300 1\n\n\n"
            "It creates a folder Results/Bval[additional_B]+[bc_B]_Bstep[longterm_rounds]+[current_rounds]_L[LENGTH]_nseq[NSEQ] with figures\n\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('--nseq', type=int, required=True, help='an integer for the number of ngenes (required)')
        # Optional argument with default value
    parser.add_argument('--bg_B', type=float, default=0.5, help='a float for background gBGC value (current human gBGC ie : non active gBGC at TSS)')
    parser.add_argument('--bg_B_iter', type =int, help = 'an integer for the index of bg_B values')

    parser.add_argument('--add_B', type=float, default=1.75, help='a float for B in the past, active gBGC at TSS')
    parser.add_argument('--add_B_iter', type =int, help = 'an integer for the index of additional_B values')
    parser.add_argument("--width", type=int,default= 300, help="sigma of the GC_distribution curve")
    parser.add_argument('--W_iter', type =int, help = 'an integer for the index of width values')
    parser.add_argument("--peak", type=int, default=80 , help="X position relative to TSS of the peak of the GC_distribution curve")
    parser.add_argument('--P_iter', type =int, help = 'an integer for the index of peak values')


    parser.add_argument('--longterm_rounds', type=int, default=500, help='an integer for number of rounds with active BGC (default : 300)')
    parser.add_argument('--current_rounds', type=int, default=100, help='an integer for the number of round with the background BGC (default: 100)')
    parser.add_argument('--T_iter', type= int, help='an integer for the index of current_rounds values')
    parser.add_argument('--length', type=int, default=5000, help='an integer for the sequence length (default : 5000)')
   # Optional list argument with default value
    parser.add_argument('--bins', type=int, nargs='+', default=[41,10], help='an integer for the number of bins  of the region around TSS and for only the exonic region (default: [41, 10])')
    parser.add_argument('--cpg', type=float, nargs='+', default=[35,1.526055912,15.123163619,1320.4637893284423,1], help='a list of float for CtoT_CpG. These numbers correspond to a C>T in CpG mutation rate that peaks 35 bp after TSS (the curve has a sigma of 1320.45, a delta of 1). The peak CpG rate is 1.52 and the baseline CpG rate is 15 times the baseline mutation rate')
    #parser.add_argument('--GtoA_cpg', type=float, nargs='+', default=[160,1.392120774,15.490597877,1368.8945044038244,1],help= 'These numbers correspond to a G>A in CpG mutation rate that peaks 160 bp after TSS (the curve has a sigma of 1368.89, a delta of 1). The peak CpG rate is 1.39 and the baseline CpG rate is 15.49 times the baseline mutation rate')
    parser.add_argument('--EIB', action='store_true')
    parser.add_argument('--mr_EIB', type=float, nargs='+', default=[0.49,0.0014,-0.157],help= 'These numbers correspond to mutation rate that peaks ... bp after TSS (the curve has a sigma of ..., a delta of 1). The peak rate is ... and the baseline rate is ... times the baseline mutation rate')
    parser.add_argument('--res_eib', type=str, default='{"CA":[-0.137,0.00136,0.175],"GC":[-0.140,0.000931,0.220],"TA": [0, 0.0853]}')
    parser.add_argument('--res_eib_iter', type= int, help='an integer for the index of res_eib ')
    parser.add_argument('--TSS', action='store_true')
    parser.add_argument('--mr_TSS', type=float, nargs='+', default=[-112,1.35,1,81.8,1,22,1.55,31.4,1],help= 'TThese numbers correspond to mutation rate that peaks ... bp after TSS (the curve has a sigma of ..., a delta of 1). The peak CpG rate is ... and the baseline rate is ... times the baseline mutation rate')
    parser.add_argument('--res_tss', type=str, default='{"GA":[-0.164,147,120,0.159,0.317,511,33.9]}')
    parser.add_argument('--res_tss_iter', type= int, help='an integer for the index of res_tss ')
    parser.add_argument("--nrep", type=int, help="Number of replicates")


    args = parser.parse_args()
    
    background_B = args.bg_B
    bg_B_iter = args.bg_B_iter
    additional_B = args.add_B
    add_B_iter = args.add_B_iter
    width= args.width
    peak= args.peak
    P_iter= args.P_iter
    W_iter= args.W_iter
    seq_length = args.length
    longterm_rounds = args.longterm_rounds
    current_rounds = args.current_rounds
    T_iter = args.T_iter
    ngenes = args.nseq
    bins_all= args.bins
    cpg= args.cpg
    TSS_args = args.TSS
    mr_TSS = args.mr_TSS
    res_tss = args.res_tss
    res_tss_iter = args.res_tss_iter
    EIB_args = args.EIB
    mr_EIB = args.mr_EIB
    res_eib = args.res_eib
    res_eib_iter = args.res_eib_iter
    nrep= args.nrep


    


    main(background_B,bg_B_iter,additional_B,add_B_iter,seq_length,longterm_rounds,current_rounds,T_iter,ngenes,bins_all,cpg,nrep, mr_TSS, mr_EIB, TSS_args, EIB_args)