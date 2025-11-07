import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse

def main(file_path,dir, B_time=0,gBGC=0 ):
    # Load the data into a pandas DataFrame
    df = pd.read_csv(file_path, delimiter=';')


    
    # Define plots
    if gBGC:
        plots = [
            ("additional_B", "background_B", "heatmap_summary_gBGC.png", "Additional B", "Background B")
    ] 
    elif B_time:
        plots = [
            ("additional_B","current_rounds","heatmap_summary_Btime.png", "Additional B", "Time OFF (rounds of Evolution)")
    ] 
    else: 
        plots = [
            ( "peak","width", "heatmap_summary_Width_peak.png", "Hotspot B peak","Hotspot B width")
    ]


    for x_col,y_col, filename, x_name, y_name in plots:
        if not gBGC and not B_time:
            df[x_col] = df[x_col].astype(int)
            df[y_col] = df[y_col].astype(int)
        # Extract relevant columns for the heatmap
        heatmap_data = df[[x_col, y_col, 'score']]

        # Pivot the data to structure it for the heatmap
        heatmap_data_pivot = heatmap_data.pivot_table(
            index=y_col, 
            columns=x_col, 
            values='score', 
            aggfunc='mean'
        )
        if B_time or gBGC:
        # Create the heatmap
            plt.figure(figsize=(20, 18))
            ax = sns.heatmap(
                heatmap_data_pivot,
                annot=True,
                cmap='viridis',
                linewidths=0.7,
                annot_kws={"size": 30}
            )
                    # Add axis labels
            plt.ylabel(y_name, fontsize=37)
            plt.xlabel(x_name, fontsize=37)
            plt.xticks(fontsize=35)
            plt.yticks(fontsize=35,rotation=0)

            
            colorbar = ax.collections[0].colorbar
            colorbar.set_label("Score", fontsize=37)
            colorbar.ax.tick_params(labelsize=35)    

            plt.tight_layout()
            # Save the plot 
            plt.savefig(os.path.join(dir,filename))
        else :
            
            plt.figure(figsize=(20, 16))
            ax = sns.heatmap(
                heatmap_data_pivot,
                annot=True,
                cmap='viridis',
                linewidths=0.7,
                annot_kws={"size": 25}
            )

            # Add axis labels
            plt.ylabel(y_name, fontsize=26)
            plt.xlabel(x_name, fontsize=26)
            plt.xticks(fontsize=30)
            plt.yticks(fontsize=30, rotation=0)

            
            colorbar = ax.collections[0].colorbar
            colorbar.set_label("Score", fontsize=30)
            colorbar.ax.tick_params(labelsize=30)    

            plt.tight_layout()
            # Save the plot 
            plt.savefig(os.path.join(dir,filename))

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="plot a heatmap from Score.csv.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Argument pour le chemin du fichier de score (remplace 'path')
    parser.add_argument('-p', '--file', 
                        required=True, 
                        help="path of the Score.csv file (ex: .../width_peak/Score.csv).")
    
    # Argument pour le répertoire de sortie (remplace 'dir')
    parser.add_argument('-d', '--dir', 
                        required=True, 
                        help="directory in which the heatmap is stored (ex: .../width_peak).")
    parser.add_argument('--gBGC', action='store_true',
                        help="1 if the score is computed to select best Background_B and best Additional_B")
    parser.add_argument('--B_Time',action='store_true', 
                        help="1 if the score is computed to select best Additional B and Best Time off")

    args = parser.parse_args()

    main(args.file, args.dir)

#Ex for computing a heatmap 
"""
python heatmap_score.py --file Results/addB_bgB/Score_addB_bgB.csv --dir Results/addB_bgB --gBGC """