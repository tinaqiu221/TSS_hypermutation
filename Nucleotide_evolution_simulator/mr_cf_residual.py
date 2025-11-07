import numpy as np 
import pandas as pd
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

### plot the residual from curve fit of mutation rate and curve fit of the average mutation rate ###

def gaussian(x, a1, b1, c1, d1, a2, b2, c2):
    gaussian_1 = d1 + a1 * np.exp(-((x - b1) ** 2) / (2 * c1 ** 2))
    gaussian_2 = a2 * np.exp(-((x - b2) ** 2) / (2 * c2 ** 2))
    return gaussian_1 + gaussian_2

def gaussian_2(x, param):
    gaussian_1 = param[3] + param[0] * np.exp(-((x - param[1]) ** 2) / (2 * param[2] ** 2))
    gaussian_2 = param[4] * np.exp(-((x - param[5]) ** 2) / (2 * param[6] ** 2))
    return gaussian_1 + gaussian_2

def exponential_inv(x, a, b, c) : 
    return 1/(1- (a * np.exp(-b*x)+c))

def exponential_decay_2(x, param) : 
    return param[0]* np.exp(-param[1]*x)+param[2]

def exponential(x,a,b,c):
    return a * np.exp(-b*x)+c

def log(x, a,b,c) : 
    return a * np.log(b*x) + c

def linear(x,a,b):
    return a*x + b

def linear_2(x,param):
    return x*param[0] + param[1]

def from_param(directory, file_mean, region):

    """ *** GOAL : get best parametres of the fitting curves for each region of the genome 
        *** Input : 
            - file_mean : str
            - region :
            - directory : str, path in which the files are
        *** Output : 
            - result : csv file, containing the data for each curve"""
        

    results = os.path.join(directory, "Results", file_mean)
    df = pd.read_csv(results, delimiter=",") #parametres of mean fitted curve
    
    
    #Extract data from results.csv
    for index, row in df.iterrows():
        
        if row["mr_type"] != "mean":
            continue
        else:
            x = np.linspace(1,999,999)
            
            if region == "TSS":
                #Fitted curve parametres (mean)
                a1, b1, c1, d1, a2, b2, c2 = row["a1"], row["b1"], row["c1"], row["d1"], row["a2"], row["b2"], row["c2"]
                f = gaussian(x,a1, b1, c1, d1, a2, b2, c2)
                
                # Initialisation of parametres for mutation rates
            
            elif region=="EIB":
                #Fitted curve parametres (mean)
                a1, b1, c1 = row["a1"], row["b1"], row["c1"]
                
                x= np.linspace(500,3499,3000)
                
                f = exponential_inv(x,a1, b1, c1)
                
            elif region== "Exon_1_forward":
                #Fitted curve parametres (mean)
                a,b,c = [row["a1"], row["b1"], row["c1"]]
                f= exponential(x,a,b,c)


            elif region=="Exon_1_backward":
                #Fitted curve parametres (mean)
                x = np.linspace(850,999,150)
                a,b = [row["a1"], row["b1"]]
                f = linear(x,a,b) 
    #print("type f", type(f))
    return f

import matplotlib.pyplot as plt
import pandas as pd

def plot_residual_score(file_path, output_dir, region):
    df = pd.read_csv(file_path, sep=";")
    plt.figure(figsize=(10,6))
    plt.bar(df['mr_type'], df['residual_score'], color='skyblue')
    plt.xlabel('Mutation Rate Type')
    plt.ylabel('Residual Score')
    plt.title('Residual Scores by Mutation Rate Type')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"score_residual_{region}.png"))


def get_diff(file_mr, file_mean, directory, region):

    output_dir = os.path.join(directory, "Results","residual")
    os.makedirs(output_dir, exist_ok=True)

    region_residual = os.path.join(output_dir, f'{region}_diff')
    os.makedirs(region_residual, exist_ok=True)
    f_mean= from_param(directory, file_mean, region)
    df_mr = pd.read_csv(os.path.join(directory, "Data",file_mr), header=0, sep=";" )
    df_mr.dropna(how='any', inplace=True)

    residual_df = pd.DataFrame(columns= df_mr.columns[:])
    residual_score = {}
    if region=='EIB':

        df_mr_masked = df_mr.mask(((df_mr['Position'] <500 ) & (df_mr['Position'] >0)), None)
        df_mr_masked.dropna(how='any', inplace=True)
        residual_df["Position"]= df_mr_masked["Position"]
        x_res_fit = np.linspace(500,3499,3000)
        x_plot= np.linspace(0,3000,3000)
        x_raw = np.linspace(-500,3000,3499)
        ymin = -1.1
        ymax = 2.6
        xmin=-550
        xmax=3050
    else: 
        residual_df["Position"]= df_mr["Position"]
        x_res_fit = np.linspace(1,999,999)
        x_plot = np.linspace(-500,500,999)
        ymin = -1.1
        ymax = 2.6
        xmin = -550
        xmax=550

    for mr_type in df_mr.columns[1:]:
        print(mr_type)
        #determine r(x) for each mutation rate in the file r(x) = f(x) - s(x)
        if region =="EIB":
            residual_df[mr_type] = f_mean - df_mr_masked[f'{mr_type}']
            residual_score[mr_type] = residual_df[mr_type].abs().sum()
        
        else: 
            residual_df[mr_type] = f_mean - df_mr[mr_type]
            residual_score[mr_type] = residual_df[mr_type].abs().sum()
    
        residual_cf_param = {"C to A per C":[-0.137,0.00136,0.175],"G to C per G":[-0.1402775, 0.0009305, 0.2201743], "T to A per T":[0,0.0853], "G to A per G":[-0.164,147,120,0.159,0.317,511,33.9]}
        cf_type= {"C to A per C":exponential_decay_2, "G to C per G": exponential_decay_2, "T to A per T": linear_2, "G to A per G": gaussian_2}
        
        if (region=="EIB" and mr_type in ["C to A per C", "G to C per G", "T to A per T"]) or (region=="TSS" and mr_type in ["G to A per G"]):
            print("ok")
            param = residual_cf_param[mr_type]
            curve = cf_type[mr_type]
            print(curve)
            print(param)

            plt.figure(figsize=(10,8))
            plt.plot(x_plot, residual_df[mr_type], label = 'residual curve', color='maroon')
            if region=="EIB":
                plt.plot(x_raw, df_mr[mr_type], label = f'{mr_type} raw data', color= "black",alpha=0.7)
            else: 
                plt.plot(x_plot, df_mr[mr_type], label = f'{mr_type} raw data',color= "black",alpha=0.7)
            plt.plot(x_plot, f_mean, label = 'mean fit', color="red",linewidth=5)
            plt.plot(x_plot,f_mean - curve(x_res_fit, param),color="green",linewidth=5)
            print(np.mean(residual_df[mr_type]), np.mean(curve(x_res_fit, param)), np.mean(f_mean))
            if region=="TSS":
                plt.xticks([-500, -250, 0, 250, 500], fontsize=20)
            else:
                plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            #plt.xlabel("Position", fontsize=24)
            #plt.ylabel("Normalised Mutation rate", fontsize=24)
            #plt.legend()
            plt.title(f"residual {mr_type}")
            plt.savefig(os.path.join(region_residual, f"residual_{mr_type}.png"))
            plt.close()
        else:
            plt.figure(figsize=(10,8))
            plt.plot(x_plot, residual_df[mr_type], label = 'residual curve', color='maroon')
            if region=="EIB":
                plt.plot(x_raw, df_mr[mr_type], label = f'{mr_type} raw data', color= "black",alpha=0.7)
            else:
                plt.plot(x_plot, df_mr[mr_type], label = f'{mr_type} raw data', color= "black",alpha=0.7)
            plt.plot(x_plot, f_mean, label = 'mean fit', color="red",linewidth=5)
            if region=="TSS":
                plt.xticks([-500, -250, 0, 250, 500], fontsize=20)
            else:
                plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            #plt.legend()
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            plt.title(f"residual {mr_type}")
            plt.savefig(os.path.join(region_residual, f"residual_{mr_type}.png"))
            plt.close()

directory = "C:\\Mutation_rate" 
region_init = ["TSS", "EIB"]
for elem in region_init:
    file_mean = f"cf_mr_{elem}.csv"
    file_mr = f"{elem}.csv"
    get_diff(file_mr, file_mean, directory, elem)

#EX
#python mr_cf_residual.py






















