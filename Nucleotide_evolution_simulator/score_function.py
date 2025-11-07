import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse

def score_function(A_path,C_path,G_path, output_file):
    A_real= [0.27472,0.27410,0.27545,0.27451,0.27511,0.27469,0.27375,0.27257,0.27145,0.27066,0.26971,0.26751,0.26425,0.26042,0.25554,0.24944,0.24176,0.23111,0.21499,0.18239,0.16343,0.16753,0.17012,0.18278,0.19788,0.21220,0.22431,0.23319,0.23953,0.24443,0.24885,0.25133,0.25347,0.25535,0.25743,0.25876,0.25948,0.25955,0.25987,0.25956,0.26076]
    #T_real = [0.27489,0.27457,0.27350,0.27412,0.27297,0.27350,0.27272,0.27168,0.27012,0.26852,0.26677,0.26371,0.26004,0.25736,0.25115,0.24365,0.23339,0.22278,0.20583,0.17985,0.17286,0.17649,0.18658,0.20624,0.22409,0.23859,0.25021,0.25950,0.26670,0.27095,0.27548,0.27620,0.27878,0.27953,0.28162,0.28277,0.28391,0.28480,0.28494,0.28485,0.28554]
    G_real = [0.22446,0.22439,0.22470,0.22495,0.22541,0.22497,0.22560,0.22662,0.22789,0.22949,0.23054,0.23415,0.23731,0.24031,0.24568,0.25237,0.26065,0.27011,0.28464,0.32015,0.34684,0.34566,0.34003,0.31928,0.29992,0.28467,0.27319,0.26341,0.25624,0.25133,0.24670,0.24518,0.24260,0.24105,0.23838,0.23763,0.23640,0.23531,0.23387,0.23520,0.23353]
    C_real = [0.22593,0.22694,0.22635,0.22642,0.22652,0.22684,0.22793,0.22913,0.23055,0.23133,0.23299,0.23464,0.23839,0.24191,0.24763,0.25455,0.26421,0.27600,0.29455,0.31760,0.31687,0.31032,0.30326,0.29170,0.27811,0.26454,0.25229,0.24389,0.23754,0.23330,0.22897,0.22728,0.22515,0.22407,0.22257,0.22083,0.22021,0.22034,0.22131,0.22039,0.22017]
    #GC_real = [0.45039,0.45133,0.45105,0.45137,0.45193,0.45181,0.45353,0.45575,0.45844,0.46082,0.46352,0.46878,0.47571,0.48223,0.49330,0.50691,0.52486,0.54611,0.57919,0.63775,0.66371,0.65598,0.64330,0.61098,0.57803,0.54921,0.52548,0.50730,0.49377,0.48462,0.47567,0.47247,0.46775,0.46512,0.46095,0.45846,0.45661,0.45566,0.45519,0.45559,0.45371]

    real_data = {'A': A_real, 'C': C_real, 'G': G_real}
    simu_data = {'A': A_path,'C': C_path, 'G': G_path}
    score_dict = {'A': [], 'C': [], 'G': []}
    df_score = pd.DataFrame.from_dict(score_dict)

    for base in ['A','C', 'G']:
        real = real_data[base]
        simu = simu_data[base]
        # Load the CSV file
        df = pd.read_csv(simu, delimiter=";" )
        if df_score.empty:
            df_score = df[['nrep','background_B', 'additional_B','res_tss','res_tss_iter','res_eib','res_eib_iter', 'width','W_iter','peak','P_iter',"current_rounds"]].copy()
            print(df_score.head(5))
            print(df_score.shape)

        # Extract columns (dynamic columns, numeric keys as strings)
        sim_columns = df.iloc[:,18:]
        
        # Convert the columns values to float
        sim_columns = sim_columns.apply(pd.to_numeric, errors='coerce')
        
        # Ensure real is a Series and has the same length as the columns in sim_columns
        real_series = pd.Series(real)

        sumdelta_list = []
        # Loop through each row in the DataFrame
        for _, row in sim_columns.iterrows():
            # Calculate sumdelta_GC for the current row
            sumdelta = (((row - real_series.values) ** 2)/((real_series.values)**2)).sum()
            sumdelta_list.append(sumdelta)
        print(len(sumdelta_list))
        # Assign the calculated sum values to the DataFrame
        df_score[base] = sumdelta_list
    df_score['score'] = df_score[['A','C','G']].sum(axis=1)

    # Save the updated DataFrame to a new CSV file
    df_score.to_csv(output_file, index=False,sep=";")

    return df_score


def plot_graph(df,path, gBGC, B_time):
    
    # Define plots
    if gBGC:
        plots = [
            ("additional_B", "background_B", "Score as a Function of bg B and add B ", "score_summary.png", "Additional B", "Background B")
    ] 
    elif B_time:
        plots = [
            ("additional_B","current_rounds","Score as a Function add B and Time off ", "score_summary.png", "Additional B", "Time OFF (rounds of Evolution)")
    ] 
    else: 
        plots = [
            ("width", "peak", "Score as a Function of the width and peak", "score_summary.png", "B width", "B peak")
    ]

    print(plots)
    for x_col,y_col, title, filename, x_name, y_name in plots:
        plt.figure(figsize=(10, 8))

        # Scatter plot with colors based on background_B
        print(x_col)
        scatter = plt.scatter(
        df[f'{x_col}'],
        df["score"],
        c=df[f'{y_col}'],   
        cmap="magma",
        edgecolor='k',
        alpha=0.8,
        s=80
    )
        # Add a colorbar to explain the mapping
        colorbar = plt.colorbar(scatter)
        colorbar.set_label(y_name, fontsize=20)  # larger colorbar label
        colorbar.ax.tick_params(labelsize=20)   # larger colorbar tick labels

        # Add plot labels and title
        plt.xlabel(x_name, fontsize=20)
        plt.ylabel("Score", fontsize=20)
        #plt.title(title, fontsize=18)

        # Increase tick label sizes on both axes
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        if B_time:
            array= np.array([0.51,0.36,0.26,0.19,0.15,0.14,0.17,0.21,0.28,0.37])
            plt.plot(np.linspace(0.75,3,10),array, linestyle="-", color="red")

        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(path, filename))
        plt.close()

def main(dir,gBGC,B_time):
    #score
    A_path = os.path.join(dir, "A_summary.csv")
    C_path = os.path.join(dir, "C_summary.csv")
    G_path = os.path.join(dir, "G_summary.csv")
    score_file = os.path.join(dir, "Score.csv")
    score = score_function(A_path,C_path,G_path,score_file)

    
    score.dropna(subset=["score"], inplace=True) 
    
    if not score.empty:
        plot_graph(score, dir, gBGC=gBGC, B_time=B_time)

if __name__ == "__main__":

    # Set up the parser for reusability
    parser = argparse.ArgumentParser(description="Analysis and visualization of simulation scores from summary files.")
    
    # Mandatory parameter
    parser.add_argument('-d', '--dir',required=True, help="Base directory containing the A_summary.csv, C_summary.csv, G_summary.csv files.")
    
    # Optional parameters for plotting (use action='store_true' for boolean flags)
    parser.add_argument('--gBGC', action='store_true', help="Activates gBGC-related plotting options (default: False).")
    parser.add_argument('--B_time', action='store_true',help="Activates B_time-related plotting options (default: False).")
    
    args = parser.parse_args()

    main(args.dir, args.gBGC, args.B_time)

#EX 
"""
python score_function.py --dir "C:\Results\addB_bgB" --gBGC

"""
