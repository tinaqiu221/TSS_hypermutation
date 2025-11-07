import pandas as pd
import os

def set_baseline(df):
    #reads the file TSS -> set the baseline
    dico = {}
    for col in df.columns[1:]:
        #print("col: ",col)
        #set the baseline for a given mr_type based on intergenic region
        baseline = df[col].iloc[0:200].median()
        #print("baseline:",baseline)
        dico[col] = baseline
    return dico

def normalisation(df, baseline):

    #read input file
    df_normalised= pd.DataFrame(columns=df.columns)
    df_normalised[df.columns[0]] = df[df.columns[0]]
    for mr_type in baseline.keys():
        #print("mr_type:", mr_type)
        df_normalised[f"{mr_type}"] = df[f"{mr_type}"]/baseline[mr_type]
            
    return df_normalised

def main(list, directory):
    
    #read the csv file as a df to set the baseline for each mr_type
    df_baseline = pd.read_csv(os.path.join(directory,f"{list[0]}.csv"), sep=";")
    #print("df_baseline:", df_baseline.head(n=5))
    
    #set the baseline value for each mr_type
    dico_baseline = set_baseline(df_baseline)
    print("dico_baseline",dico_baseline )

    for file in list:
        input_file =  f"{file}.csv"
        output_file = f"{file}_normalised.csv"

        #input df that is to be normalised
        df_mr = pd.read_csv(os.path.join(directory,input_file), sep=";")
        
        #conduct the normalisation
        df_normalised = normalisation(df_mr, dico_baseline)
        #print("df_normalised:",df_normalised.head(n=5))
        
        #write the output file
        df_normalised.to_csv(os.path.join(directory, output_file), sep=";",index=False)

list = ["TSS_raw", "EIB_raw"]
directory = "C:/Users/eh263/OneDrive/Documents/Emilie/Research/Mutation_rate/Source"
main(list, directory)

