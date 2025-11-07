import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import csv

### Curve fitting for average fit and residual ###

def gaussian(x, a1, b1, c1, d1, a2, b2, c2):
    gaussian_1 = d1 + a1 * np.exp(-((x - b1) ** 2) / (2 * c1 ** 2))
    gaussian_2 = a2 * np.exp(-((x - b2) ** 2) / (2 * c2 ** 2))
    return gaussian_1 + gaussian_2

def simple_gaussian(x, a1, b1, c1, d1):
    gaussian = d1 + a1 * np.exp(-((x - b1) ** 2) / (2 * c1 ** 2))
    return gaussian

def exponential_inv(x, a, b, c) : 
    return 1/(1- (a * np.exp(-b*x)+c))

def exponential_decay(x, a, b, c) : 
    return a * np.exp(-b*x)+c

def log(x, a, b, c) : 
    return a * np.log(b*x) + c

def linear(x,a,b):
    return a*x + b

def get_parametres(region, directory, mean=True):
    """ *** GOAL : get best parametres of the fitting curves for each region of the genome 
        *** Input : 
            - region, string, a region of interest (TSS, EIB, etc)
            - directory : str, path in which the files are
        *** Output : 
            - result : csv file, containing the data of each curve
            - df, data frame containing for a given region """

    
    filepath = os.path.join(directory,"Source", f"{region}.csv")
    df = pd.read_csv(filepath, delimiter=";", header=0 )
    
    if mean:
        df['mean'] = df.iloc[:, 1:].mean(axis=1)
        df.dropna(how='any', inplace=True)


    #Initialisation of the curve parametres
    for column in df.columns[1:]:
        a1,b1,c1,d1,a2,b2,c2 = [None for i in range(7)]
        
        if region == "TSS": 
                params, _ = curve_fit(gaussian, df['Position'], df[f'{column}'], p0= [1.35, 380, 100, 1, 1.55, 500, 40])
                fitted_curve = gaussian(df['Position'], *params)
                a1,b1,c1,d1,a2,b2,c2 = params
                #print("params TSS:", params)
        
        elif region=="Exon_1_forward":
                params, _ = curve_fit(exponential_inv, df['Position'], df[f'{column}'])
                fitted_curve = exponential_inv(df['Position'], *params)
                a1, b1, c1 = params

        elif region=="EIB":
                     
                masked_df = df.mask(((df['Position'] <500 ) & (df['Position'] >0)), None)
                masked_df.dropna(how='any', inplace= True)
                #print("masked_df[Position] for EIB :",masked_df['Position'])

                params, _ = curve_fit(exponential_inv, masked_df['Position'], masked_df[f'{column}'],p0 = [0.49,0.005,0],maxfev=100000)
                fitted_curve = exponential_inv(masked_df['Position'], *params)
                a1,b1,c1= params

        elif region=="Exon_1_backward":
                masked_df = df.mask((df['Position'] < 850), None)
                masked_df.dropna(how='any', inplace= True)

                params, _ = curve_fit(linear, masked_df['Position'], masked_df[f'{column}'])
                fitted_curve = linear(masked_df['Position'], *params)
                a1, b1 = params

        # Analyze the Fitted curve 
        ssr = np.sum((df[f'{column}'] - fitted_curve)**2)
        rmse = np.sqrt(ssr / len(df[f'{column}']))

        ##stock the values -> in a file

        result = os.path.join(directory,"Results", f"cf_mr_{region}.csv")

        with open(result, 'a', newline='') as file:
            writer = csv.writer(file)
            header = ["mr_type", "a1", "b1", "c1", "d1", "a2", "b2", "c2" ,"SSR", "RMSE"]
            #print(header)
            
            if  os.path.getsize(result) == 0:
                writer.writerow(header)
                writer.writerow([f"{column}",f"{a1}" ,f"{b1}",f"{c1}",f"{d1}",f"{a2}" ,f"{b2}",f"{c2}", f"{ssr}", f"{rmse}"])
            else: 
                writer.writerow([f"{column}",f"{a1}" ,f"{b1}",f"{c1}",f"{d1}",f"{a2}" ,f"{b2}",f"{c2}", f"{ssr}", f"{rmse}"])
        
        
    return df

def get_parametres_residual(region, directory):
    """ *** GOAL : get best parametres of the fitting curves for each region of the genome 
        *** Input : 
            - region, string, a region of interest (TSS, EIB, etc)
            - directory : str, path in which the files are
        *** Output : 
            - result : csv file, containing the data of each curve
            - df, data frame containing for a given region """

    
    filepath = os.path.join(directory,"Results", "residual", f"residual_{region}.csv")
    df = pd.read_csv(filepath, delimiter=";")


    #Initialisation of the curve parametres
    for column in df.columns[1:]:
        a1,b1,c1,d1,a2,b2,c2 = [None for i in range(7)]
        
        if region == "TSS": 
                params, _ = curve_fit(gaussian, df['Position'], df[f'{column}'], p0= [0.3, 300, 100, 0,0.4,500,50],maxfev = 100000)
                fitted_curve = gaussian(df['Position'], *params)
                a1,b1,c1,d1,a2,b2,c2 = params
                #print("params TSS:", params)
        

        elif region=="EIB":
                if column in ["T to A per T", "T to G per T","A to T per A","A to C per A"]:
                     
                    masked_df = df.mask(((df['Position'] <500 ) & (df['Position'] >0)), None)
                    masked_df.dropna(how='any', inplace= True)
                    #print("masked_df[Position] for EIB :",masked_df['Position'])

                    params, _ = curve_fit(linear, masked_df['Position'], masked_df[f'{column}'],maxfev=100000)
                    fitted_curve = linear(masked_df['Position'], *params)
                    a1,b1= params
                else:
                    masked_df = df.mask(((df['Position'] <500 ) & (df['Position'] >0)), None)
                    masked_df.dropna(how='any', inplace= True)
                    #print("masked_df[Position] for EIB :",masked_df['Position'])

                    params, _ = curve_fit(exponential_decay, masked_df['Position'], masked_df[f'{column}'],p0 = [5000,0.01,-130],maxfev=100000)
                    fitted_curve = exponential_decay(masked_df['Position'], *params)
                    a1,b1,c1= params

        # Analyze the Fitted curve 
        ssr = np.sum((df[f'{column}'] - fitted_curve)**2)
        rmse = np.sqrt(ssr / len(df[f'{column}']))

        ##stock the values -> in a file

        result = os.path.join(directory,"Results","residual", f"cf_mr_residual_{region}.csv")

        with open(result, 'a', newline='') as file:
            writer = csv.writer(file)
            header = ["mr_type", "a1", "b1", "c1", "d1", "a2", "b2", "c2" ,"SSR", "RMSE"]
            #print(header)
            
            if  os.path.getsize(result) == 0:
                writer.writerow(header)
                writer.writerow([f"{column}",f"{a1}" ,f"{b1}",f"{c1}",f"{d1}",f"{a2}" ,f"{b2}",f"{c2}", f"{ssr}", f"{rmse}"])
            else: 
                writer.writerow([f"{column}",f"{a1}" ,f"{b1}",f"{c1}",f"{d1}",f"{a2}" ,f"{b2}",f"{c2}", f"{ssr}", f"{rmse}"])
        
        
    return df

def plot_fitted_curve(x_raw, y_raw, x_fit, y_fit, analysis, mr_type, save_dir, residual=False):
    plt.figure(figsize=(10, 6))
    if residual:
        plt.plot(x_raw, y_raw, label=f"{mr_type}", linestyle='dotted', color='maroon')
        plt.plot(x_fit, y_fit, label="Fitted curve", color="green")
    else:
        plt.plot(x_raw, y_raw, label=f"{mr_type}", linestyle='dotted')
        plt.plot(x_fit, y_fit, label="Fitted curve", color="red")
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    plt.xlabel("Position")
    plt.ylabel("Mutation Rate")
    plt.title(f"Curve fit for {analysis} {mr_type}")
    plt.savefig(os.path.join(save_dir, f"{mr_type}.png"))
    plt.close()
                    
def graph_fit_curve(directory,filename_cf,df_mr,analysis):
        
        """*** GOAL : plot best parametres of the fitting curves for mutation type 
        *** Input : 
            - directory : str, path in which the files are
            - filename_cf
            - df_realdata : data frame, that contains mutation rates per region or per mr_type
            - analysis : string, that is either : summary, TSS, EIB, Exon_1_forward,Exon_1_backward
        *** Output : 
            - graphs : png files"""


        df_cf = pd.read_csv(filename_cf) #file with curve parametres

        #creation of the folder for the analysis needed (TSS, EIB, etc)
        output_dir = os.path.join(directory, "Results", analysis)
        os.makedirs(output_dir, exist_ok=True)    

        for index, row in df_cf.iterrows():
            mr_type = row["mr_type"]

            if analysis=="TSS":
                a1,b1,c1,d1,a2,b2,c2 = row["a1"], row["b1"], row["c1"], row["d1"], row["a2"], row["b2"], row["c2"]
                x_raw = np.linspace(1,999,999)
                
                y_raw = df_mr[mr_type]
                y_fit = gaussian(x_raw,a1,b1,c1,d1,a2,b2,c2)
                
                plot_fitted_curve(x_raw, y_raw, x_raw, y_fit, analysis, mr_type,output_dir)

            elif analysis == "Exon_1_backward":
                a1,b1 = row["a1"], row["b1"]

                x_raw = np.linspace(1,999,999)
                y_raw = df_mr[mr_type]
                x_fit = np.linspace(850,999,150)
                y_fit = linear(x_fit, a1,b1)

                plot_fitted_curve(x_raw, y_raw, x_fit, y_fit, analysis, mr_type, output_dir)

            elif analysis in ["Exon_1_forward", "EIB"]:
                a1,b1,c1 = row["a1"], row["b1"], row["c1"]
                if analysis == "Exon1_forward":
                    x_raw = np.linspace(1,999,999)
                    x_fit = x_raw
                else:
                    x_fit = np.linspace(500,3499,3000)
                    x_raw = np.linspace(1,3499,3499)
                    #print("length x_raw: ", len(x_raw))
                
                y_raw = df_mr[mr_type]
                y_fit = exponential_inv(x_fit,a1,b1,c1)

                plot_fitted_curve(x_raw, y_raw, x_fit, y_fit, analysis, mr_type, output_dir)

            else:
                continue  # Skip unknown analysis types

def graph_fit_curve_residual(directory,filename_cf,df_mr,analysis):
        
        """*** GOAL : plot best parametres of the fitting curves for mutation type 
        *** Input : 
            - directory : str, path in which the files are
            - filename_cf
            - df_realdata : data frame, that contains mutation rates per region or per mr_type
            - analysis : string, that is either : summary, TSS, EIB, Exon_1_forward,Exon_1_backward
        *** Output : 
            - graphs : png files"""


        df_cf = pd.read_csv(filename_cf) #file with curve parametres

        #creation of the folder for the analysis needed (TSS, EIB, etc)
        output_dir = os.path.join(directory, "Results","residual", f"{analysis}_cf_diff")
        os.makedirs(output_dir, exist_ok=True)    

        for index, row in df_cf.iterrows():
            mr_type = row["mr_type"]

            if analysis=="TSS":
                a1,b1,c1,d1,a2,b2,c2 = row["a1"], row["b1"], row["c1"], row["d1"], row["a2"], row["b2"], row["c2"]
                x_raw = np.linspace(1,999,999)
                
                y_raw = df_mr[mr_type]
                y_fit = gaussian(x_raw,a1,b1,c1,d1,a2,b2,c2)
                
                plot_fitted_curve(x_raw, y_raw, x_raw, y_fit, analysis, mr_type, output_dir, residual=True)

            else:
                x_fit = np.linspace(500,3499,3000)
                #print("length x_raw: ", len(x_raw))
                
                y_raw = df_mr[mr_type]
                if mr_type in ["T to A per T", "T to G per T","A to T per A","A to C per A"] :
                    a1,b1 = row["a1"], row["b1"]
                    y_fit = linear(x_fit,a1,b1)
                else: 
                     a1,b1,c1 = row["a1"], row["b1"], row["c1"]
                     y_fit = exponential_decay(x_fit,a1,b1,c1)

                plot_fitted_curve(x_fit, y_raw, x_fit, y_fit, analysis, mr_type, output_dir, residual=True)


def main(directory, residual=False):

    #list = ["TSS","EIB","Exon_1_forward","Exon_1_backward"]
    list = ["TSS","EIB"]
    for elem in list :
        if residual: 
            df_final = get_parametres_residual(elem,directory)
            filename_cf = os.path.join(directory,"Results","residual", f"cf_mr_residual_{elem}.csv")
            graph_fit_curve_residual(directory,filename_cf,df_final,elem)
        
        else:
            df_final = get_parametres(elem,directory)
            #print("df_final:",df_final.head(n=10))
            filename_cf = os.path.join(directory,"Results",f"cf_mr_{elem}.csv")
            #print("filename_cf:", filename_cf)
            graph_fit_curve(directory,filename_cf,df_final,elem)
            

if __name__ == "__main__":
    path = "C:\\Users\\eh263\\OneDrive\\Documents\\Emilie\\Research\\Mutation_rate"
    main(path, residual=True)
