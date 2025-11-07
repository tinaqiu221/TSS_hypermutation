import numpy as np
import pandas as pd
import os
import ast
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

###### Preparation of the calculation of the AIC score ######

def filter_df(df,factor1_name, factor2_name, factor1, factor2):
    """ ***GOAL : filter the lines of the df to find the ones that match factor 1 and factor 2"""
    factor1_name = factor1_name.strip()
    factor2_name = factor2_name.strip()
    df.columns = df.columns.str.strip()
    if factor1_name and factor2_name in ['res_eib', 'res_tss']:
        df[factor1_name] = df[factor1_name].apply(ast.literal_eval)
        df[factor2_name] = df[factor2_name].apply(ast.literal_eval)

    df_filtered = df[df[factor1_name].apply(lambda x: x == factor1) &
                     df[factor2_name].apply(lambda x: x == factor2)]
    if df_filtered.empty:
        print(f"[Warning] No match for {factor1_name}={factor1}, {factor2_name}={factor2}")
        return None
    return df_filtered



#compute the mean and variance between observed and simulated values for a given model replicate set

def prep_score(A_path, C_path, G_path, out_dir, factor1_name, factor2_name, factor1, factor2):
    simu_data = {'A': A_path, 'C': C_path, 'G': G_path}
    dfs_mean = []
    dfs_var = []
    
    for base, simu_choice in simu_data.items():
        df = pd.read_csv(simu_choice, sep=";")
        
        df_filtered = filter_df(df, factor1_name, factor2_name, factor1, factor2)
        if df_filtered is None:
            print(f"No data found for model {factor1_name} in file {simu_choice}. Skipping.")
            continue
        
        base_cols = ['nrep','background_B','additional_B',
                     'res_tss','res_tss_iter','res_eib','res_eib_iter']
        
        df_mean = df_filtered[base_cols].copy()
        df_var  = df_filtered[base_cols].copy()

        for col in ['res_eib', 'res_tss']:
            df_mean[col] = df_mean[col].astype(str)
            df_var[col] = df_var[col].astype(str)
        
        if factor1_name == "res_eib" :
            df_mean[base] = df_filtered.iloc[:, -12:].mean(axis=1)
            df_var[base]  = df_filtered.iloc[:, -12:].var(axis=1)
            #print("shape df_filtered: ",df_filtered.iloc[:, -12:].shape)
        dfs_mean.append(df_mean)
        dfs_var.append(df_var)
    
    if not dfs_mean:
        print("No valid data found in any file.")
        return None, None
    
    df_score_mean = dfs_mean[0]
    df_score_var  = dfs_var[0]
    for dfm, dfv in zip(dfs_mean[1:], dfs_var[1:]):
        df_score_mean = pd.merge(df_score_mean, dfm, on=base_cols, how="outer")
        df_score_var  = pd.merge(df_score_var, dfv, on=base_cols, how="outer")
    
    out_mean_path = os.path.join(out_dir, "score_mean.csv")
    out_var_path  = os.path.join(out_dir, "score_var.csv")
    df_score_mean.to_csv(out_mean_path, index=False, sep=";", header= not os.path.exists(out_mean_path),mode="a")
    df_score_var.to_csv(out_var_path,  index=False, sep=";", header = not os.path.exists(out_var_path), mode="a")
    
    print(f"Processed files saved to: {out_dir}")
    return df_score_mean, df_score_var

#compute AIC for a given model replicate set
def compute_aic(df_mean, df_var, factor1,factor1_name,factor2, factor2_name,k, out_dir,EIB):
    """
    *** GOAL : 
        - Compute AIC for a given model replicate set.
    *** INPUT :
        - res_eib : str, Identifier of the model (column "res_eib" in the CSV)
    *** OUTPUT :
        - aic : float, Computed AIC value
    *** METHOD :
        - Isolate replicates for the specified model.
        - Compute mean and variance across replicates.
        - Calculate Gaussian log-likelihood.
        - Compute AIC using the formula: AIC = 2k - 2ln(L)
    *** NOTE :
        - Ensure variance is not zero to avoid division errors.
    """

    if factor1_name in ["res_eib","res_tss"] and factor2_name in ["res_eib","res_tss"]:
        if EIB:
            #from EIB to end of sequence (EIB analysis)
            A_real= [0.24443,0.24885,0.25133,0.25347,0.25535,0.25743,0.25876,0.25948,0.25955,0.25987,0.25956,0.26076]
            G_real = [0.25133,0.24670,0.24518,0.24260,0.24105,0.23838,0.23763,0.23640,0.23531,0.23387,0.23520,0.23353]
            C_real = [0.23330,0.22897,0.22728,0.22515,0.22407,0.22257,0.22083,0.22021,0.22034,0.22131,0.22039,0.22017]
        else : 
            #around TSS (-875, +875 around TSS) (TSS analysis)
            A_real= [0.26425,0.26042,0.25554,0.24944,0.24176,0.23111,0.21499,0.18239,0.16343,0.16753,0.17012,0.18278,0.19788,0.2122,0.22431]
            G_real = [0.23731,0.24031,0.24568,0.25237,0.26065,0.27011,0.28464,0.32015,0.34684,0.34566,0.34003,0.31928,0.29992,0.28467,0.27319]
            C_real = [0.23839,0.24191,0.24763,0.25455,0.26421,0.276,0.29455,0.3176,0.31687,0.31032,0.30326,0.2917,0.27811,0.26454,0.25229]
        
    real_data = {'A': A_real, 'C': C_real, 'G': G_real}

    logL_dic = {"A" : [],"C" : [],"G" : []}
    
    for base in real_data.keys():
        obs = np.array(real_data[base])
        sim_mean = np.mean(df_mean[base])
        print("sim_mean:",sim_mean)
        #print(type(sim_mean))
        sim_var = np.mean(df_var[base])

        
        # Gaussian log-likelihood
        logL_dic[base] = -0.5 * np.sum(((obs - sim_mean)**2 / sim_var) + np.log(2 * np.pi * sim_var))
        print("sim_var=",sim_var)
        print("(obs-sim_mean **2 /var ) = ",((obs - sim_mean)**2 / sim_var))
        print("np.log(2 * np.pi * sim_var)=", np.log(2 * np.pi * sim_var))   
        print(f"log_L for {base}, {factor1}, {factor2}:", logL_dic[base])
        #Somme des score logL

    
    df_result = pd.DataFrame([logL_dic])
    df_result["aic"] = 2*k - (logL_dic['A'] + logL_dic['C'] + logL_dic['G'])
    df_result[f"{factor1_name}"] = f"{factor1}"
    df_result[f"{factor2_name}"] = f"{factor2}"
    df_result["k"] = k
    
    path = os.path.join(out_dir, f"aic_score_{factor1_name}_{factor2_name}.csv")
    df_result.to_csv(path, mode="a", sep=";",header=not os.path.exists(path), index=False)


def nb_parametrers(model):
    """
    Return the number of parameters (param) for a given model based on its identifier.
    """
    nb_param = 0
    for base in model.keys():
        nb_param += len(model[base])
    return nb_param

###### Selection of the minimal AIC score and Plot #######

#pick the lowest score
def min_aic(aic):
    minimum = aic["aic"].min()
    row = aic[aic["aic"] == minimum]

    print(f"the row which has the lowest AIC is/are")
    print(row)

    return minimum
    
def aic_score_comp(aic, dir, minimum):
    aic["min_diff"] = np.abs(aic["aic"] - minimum)
    aic.to_csv(os.path.join(dir, "aic_score_diff.csv"),sep = ";" ,index=False)

    return aic

def plot_aic_scores(aic, dir):
    plt.figure(figsize=(10,6))

    # Crée un axe x si tu n’as pas de colonne "index"
    if "index" not in aic.columns:
        aic = aic.reset_index()
    
    x = aic["index"]
    y = aic["min_diff"]

    # Couleurs selon les intervalles
    colors = []
    for val in y:
        if val <= 2:
            colors.append("green")
        elif 4 <= val <= 7:
            colors.append("orange")
        elif val > 10:
            colors.append("red")
        else:
            colors.append("blue")  # optionnel: pour les autres valeurs

    # Scatter plot avec couleurs
    plt.scatter(x, y, c=colors, s=50)
    plt.title("AIC scores (min_diff)")
    plt.xlabel("Index")
    plt.ylabel("Difference from min AIC")
    plt.grid(True)
    
    legend_patches = [
        mpatches.Patch(color='green', label='≤ 2'),
        mpatches.Patch(color='blue', label='3–4'),
        mpatches.Patch(color='orange', label='4–7'),
        mpatches.Patch(color='red', label='> 10')
    ]
    plt.legend(handles=legend_patches)

    # Sauvegarde
    plt.tight_layout()
    plt.savefig(os.path.join(dir, "AIC_score.png"))
    plt.close()

def main(factor1_name, factor2_name,dir,EIB,factor1_file=None,factor2_file=None):
    
    A_path = os.path.join(dir, "A_summary.csv")
    C_path = os.path.join(dir, "C_summary.csv")
    G_path = os.path.join(dir, "G_summary.csv")

    out_dir = os.path.join(dir, "AIC")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    factors_values={}

    if factor1_file and factor2_file :
        files = {f"{factor1_name}" : factor1_file,
                 f"{factor2_name}" : factor2_file }
        
        # Initialize the container for the parsed factor dictionaries
        factors_values = {}

        for factor in files.keys():
            # 1. Read the raw file content
            with open(files[factor], "r") as f:
                factor_data = f.read().strip()
            
            factor_values = []
            
            # 2. Iterate through lines and convert them to Python dict objects
            for line in factor_data.splitlines():
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                    
                try:
                    # CRITICAL: Convert the string representation to a Python dictionary
                    row_dic = ast.literal_eval(line)
                    if not isinstance(row_dic, dict):
                        print(f"[Error] Line in {files[factor]} was not a valid dictionary: {line}")
                        continue
                    factor_values.append(row_dic)
                except (ValueError, SyntaxError) as e:
                    print(f"[Error] Failed to evaluate dictionary from line: {line}. Error: {e}")
                    continue
            
            # 3. Store the list of dictionary objects
            factors_values[factor] = factor_values
    

    for factor1 in factors_values[factor1_name]:
        print(f"Processing factor1: {factor1}")
        print("type factor1", type(factor1))
        for factor2 in factors_values[factor2_name]:
            factor
            print(f"  Processing factor2: {factor2}")
            print("type factor2", type(factor2))

            df_mean, df_var = prep_score(A_path,C_path,G_path,out_dir,factor1_name, factor2_name, factor1, factor2)
            if df_mean is None and df_var is None:
                print(f"[Warning] No data for factor1={factor1}, factor2={factor2}. Skipping.")
                continue   # saute cette combinaison et passe à la suivante
            else:
                model = filter_df(df_mean,factor1_name, factor2_name, factor1, factor2)
                if model is None or model.empty:
                    print(f"[Warning] No rows match factor1={factor1}, factor2={factor2} in df_mean. Skipping.")
                    continue
            # Assurez-vous que les variables eib_param et tss_param sont initialisées à 0
            eib_param = 0
            tss_param = 0

            if factor1_name in ["res_eib", "res_tss"] and factor2_name in ["res_eib", "res_tss"]:
                
                # --- Coping with Factor 1 (res_eib or res_tss) ---
                factor1_data = model[f"{factor1_name}"].tolist()[0]
                
                if factor1_name == "res_eib":
                    eib_param = nb_parametrers(factor1_data)
                    print(f"eib_param (from {factor1_name}): {eib_param}")
                elif factor1_name == "res_tss":
                    tss_param = nb_parametrers(factor1_data)
                    print(f"tss_param (from {factor1_name}): {tss_param}")

                # --- Coping with Factor 2  (res_eib or res_tss) ---
                factor2_data = model[f"{factor2_name}"].tolist()[0]

                if factor2_name == "res_eib":
                    # S'il est différent, on ajoute. Si factor1_name était déjà res_eib, on additionne les contributions (moins probable dans ce contexte)
                    eib_param += nb_parametrers(factor2_data)
                    print(f"eib_param (from {factor2_name}): {eib_param}")
                elif factor2_name == "res_tss":
                    tss_param += nb_parametrers(factor2_data)
                    print(f"tss_param (from {factor2_name}): {tss_param}")

                # --- Calculation of k ---
                # k is the total number of parameters in the simulation
                
                k =  28 + eib_param + tss_param 
                
                print(f"Total k = 28 + {eib_param} (EIB) + {tss_param} (TSS) = {k}")
                print(f"Computing AIC for {factor1_name}={factor1}, {factor2_name}={factor2}, k={k}")
                compute_aic(df_mean, df_var, factor1, factor1_name, factor2, factor2_name, k, out_dir,EIB)

    #calculate the min score
    final_df = pd.read_csv(os.path.join(out_dir,f"aic_score_{factor1_name}_{factor2_name}.csv"), sep=";")
    minimum = min_aic(final_df)

    #compute the difference between min and all the other scores
    aic_comp = aic_score_comp(final_df, out_dir, minimum)
    print("score comp was computed")
    #plot the results
    plot_aic_scores(aic_comp, out_dir)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Simulation results analysis script, selecting factors and paths via the command line.",formatter_class=argparse.RawTextHelpFormatter)

    # Mandatory parameters
    parser.add_argument('-d', '--dir', required=True, type =str , help="Base directory where the simulation results are located.")
    parser.add_argument('--factor1_name', required=True,type=str,help="Name of the first factor (e.g., 'res_eib', 'background_B', 'width').")
    parser.add_argument('--factor2_name',required=True,type=str, help="Name of the second factor (e.g., 'res_tss', 'additional_B', 'peak').")
    parser.add_argument('--EIB',action='store_true', help="if equal to 1, the AIC score will be computed from EIB to end of sequence")
    
    # Optional parameters (for EIB/TSS combinations)
    parser.add_argument('--factor1_file',default=None,type=str,help="Full path to the combinations file for Factor 1 (Optional, only needed for res_eib).")
    parser.add_argument('--factor2_file', default=None,type=str, help="Full path to the combinations file for Factor 2 (Optional, only needed for res_tss).")
    args = parser.parse_args()

    main(args.factor1_name, 
         args.factor2_name, 
         args.dir,
         args.EIB, 
         args.factor1_file, 
         args.factor2_file)

#Exemple:
#python AIC_score.py --factor1_name "res_eib" --factor2_name "res_tss" --dir "results/res_eib" --factor1_file "results/res_eib/" --factor2_file "results/res_eib/"

#if res_TSS is of interest (the parameter that varies)

"""python AIC_score.py --factor1_name "res_eib" 
--factor2_name "res_tss" 
--dir "C:/Results/res_TSS" 
--factor1_file "C:/Results/residual/combinations_EIB_2.txt"  
--factor2_file "C:/Results/residual/combinations_TSS.txt" """

#if res_eib is of interest (the parameter that varies)

"""python AIC_score.py --factor1_name "res_eib" 
--factor2_name "res_tss" 
--dir "C:/Users/eh263/OneDrive/Documents/Emilie/Research/Cluste/Results/res_eib" 
--factor1_file "C:Results/residual/combinations_EIB.txt"  
--factor2_file "C:Results/residual/combinations_TSS.txt" 
--EIB """