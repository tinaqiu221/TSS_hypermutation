import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt


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

def CpG_mutation_rate(cpg, half_length):

    """ *** GOAL : generate a distribution (gaussian) of the CpG mutation rates in a region of interest
        *** Input : 
                    - cpg : list, containing CtoT or GtoA gaussian information
                    - half_length : integer corresponding to the position of TSS
        *** Output : 
                    - distribution : list, containing CpG mutation rates in the region of interest"""

    peak_rel_TSS = cpg[0]
    peak_mut = cpg[1]
    base_mut = cpg[2]
    sigma = cpg[3]
    delta = cpg[4]
    epsilon = 0


    distribution = [
        (base_mut-peak_mut) * (1-(pdf_trans(i/sigma,
                                                    epsilon, delta)/pdf_trans(0, epsilon, delta)))+peak_mut
        for i in range(- half_length - peak_rel_TSS, half_length - peak_rel_TSS)
    ]
    #print("distribution cpg length after function was computed:",len(distribution))
    return distribution

def exponential_inv(x,a,b,c):
    return 1/(1-(a*np.exp(-b*x)+c))

def exponential_decay(x,a,b,c):
    return a*np.exp(-b*x)+c

def linear(x,a,b):
    return a * x + b

def gaussian(x, a1, b1, c1, d1, a2, b2, c2):
    gaussian_1 = d1 + a1 * np.exp(-((x - b1) ** 2) / (2 * c1 ** 2))
    gaussian_2 = a2 * np.exp(-((x - b2) ** 2) / (2 * c2 ** 2))
    return gaussian_1 + gaussian_2

def TSS_mutation_rate(mr_region):

    """ *** GOAL : generate a distribution of specific mutation rates in a region of interest
        *** Input : 
                    - mr_region : list, containing the information related to the two gaussian curves to model the region of interest
                    - length : integer corresponding to the position of the region og interest (TSS, EIB)
        *** Output : 
                    - distribution : list, containing CpG mutation rates in the region of interest"""
    #print("position of TSS:",length)
    ## Parametres
    #1st gaussian
    length= 2500
    peak_rel_region = mr_region[0]
    peak_mut= mr_region[1]
    base_mut= mr_region[2]
    sigma = mr_region[3]
    delta = mr_region[4]
    epsilon = 0

    #2nd gaussian
    peak_rel_region_2 = mr_region[5]
    peak_mut_2 = mr_region[6]
    sigma_2 = mr_region[7]
    delta_2 = mr_region[8]
    epsilon_2 = 0


    ### region distribution
    distribution = [
        (
    base_mut +
    #1st Gaussian component
    (peak_mut - base_mut) * (pdf_trans((i - peak_rel_region)/sigma,epsilon, delta)
                                                /pdf_trans(0, epsilon, delta)))
    +

    (
    #2nd Gaussian component
    (peak_mut_2 - base_mut) * (pdf_trans((i - peak_rel_region_2)/sigma_2,
                                                    epsilon_2, delta_2)/pdf_trans(0, epsilon_2, delta_2)))
        for i in range(-length, length)
    ]
    #print("distribution length after function was computed:",len(distribution))
    return distribution

def EIB_mutation_rate(mr_region):

    x = np.arange(500,3501)# this range is set based on real_data from Tina (x data ranging from 500 to 3500)
    a = mr_region[0]
    b = mr_region[1]
    c = mr_region[2]
    distribution = exponential_inv(x,a,b,c) #the X values matter because the curve fit was done on x data ranging from 500 to 3500

    """#plot
    plt.plot(np.arange(TSS+SS, (TSS+SS)+3000), distribution, color='purple')
    plt.title("EIB distribution")
    plt.show()"""
    #print("distribution length after function was computed:",len(distribution))
    return distribution

"""***GOAL : the two following functions enable to correct the average_curve that model the region effect of TSS, or EIB 
            on the baseline mutation rate
    *** Input :
                - f_mean : np array, average curve data
                - mut_type : str, type of mutation (only for EIB)
                - mut_residual : list, parameters of the residual curve 
    *** Output : f_mean - residual : np array, corrected TSS/EIB effect depending on the mutation type"""

def residual_EIB(f_mean, mut_type, mut_residual):
    x = np.arange(500,3501)
    if mut_type in ["TA", "TG","AT","AC"]:
        a,b  = mut_residual[0], mut_residual[1]
        residual = linear(x,a,b)
    else: 
        a,b,c = mut_residual[0], mut_residual[1], mut_residual[2]
        residual = exponential_decay(x,a,b,c)
    distribution = f_mean - residual
    """plt.plot(x, distribution)
    plt.title(f"{mut_type} eib residual correction")
    plt.show()"""
    
    return distribution

def residual_TSS(f_mean,mut_type,mut_residual):
    TSS=2500
    x = np.arange(1,3001)
    a1,b1,c1,d1,a2,b2,c2 = mut_residual[0], mut_residual[1], mut_residual[2], mut_residual[3],mut_residual[4],mut_residual[5],mut_residual[6]
    residual = gaussian(x,a1,b1,c1,d1,a2,b2,c2)
    distribution = f_mean[TSS-500:] - residual
    """plt.plot(x, distribution)
    plt.title(f"{mut_type} residual TSS")
    plt.show()"""
    return distribution


#gBGC 
def calculate_factor(B):
    """ using B = 4*Ne*b """
    Ne = 10000
    b = B/(4*Ne)
    """ calculate pi (probablity of fixation) for WS mutations at gBGC sites = (1-e^2b)/(1-e^4Neb) """
    piWSgBGC = (1-np.exp(-2*b))/(1-np.exp(-4*Ne*b))
    """ finding the relative probability of pi WS mutations at gBGC vs other sites = piWSgBGC/piWS """
    """ fixation of neutral mutation is 1/2Ne """
    piWS = 1/(2*Ne)
    return (piWSgBGC/piWS)

#bg gBGC
def background_gBGC(background_B,nuc_combN,mut_rate_matrix):
    #background_B = 0.8 # 0.8 or 0.23 # instead of 0.8
    background_WS = calculate_factor(background_B) #seq simulation
    background_SW = calculate_factor(-background_B)
    # Modify the values of mut_rate_init directly using numpy indexing and operations
    #Transitions
    mut_rate_matrix[nuc_combN['AG']] *= background_WS
    mut_rate_matrix[nuc_combN['GA']] *= background_SW

    #Transversions
    mut_rate_matrix[nuc_combN['AC']] *= background_WS
    mut_rate_matrix[nuc_combN['GT']] *= background_SW

    return mut_rate_matrix


def hotspot_gBGC(additional_B, gBGC_TSS, nuc_combN,length,mut_rate_gBGC):
    #gbgc at TSS
    WS_factor = calculate_factor(additional_B)
    SW_factor = calculate_factor(-additional_B)

    GC_peak_rel_TSS = gBGC_TSS["GC_peak_rel_TSS"]
    GC_sigma = gBGC_TSS["GC_sigma"]
    GC_delta = gBGC_TSS["GC_delta"]
    GC_epsilon = gBGC_TSS["GC_epsilon"]
    GC_range = np.arange( - length - GC_peak_rel_TSS,
                        length - GC_peak_rel_TSS)
    GC_distribution = (pdf_trans(GC_range / GC_sigma, GC_epsilon, GC_delta) /
                    pdf_trans(0, GC_epsilon, GC_delta))

    # Compute mut_rate_init_gBGC using vectorized operations and modify the mutation rates for the nuc affected by gBGC
    #Transitions
    mut_rate_gBGC[nuc_combN['AG']] *=(GC_distribution * (WS_factor - 1) + 1) 
    mut_rate_gBGC[nuc_combN['GA']] *= (GC_distribution * (SW_factor - 1) + 1)
    #Transversions
    mut_rate_gBGC[nuc_combN['AC']] *= (GC_distribution * (WS_factor - 1) + 1)
    mut_rate_gBGC[nuc_combN['GT']] *= (GC_distribution * (SW_factor - 1) + 1)
    
    return mut_rate_gBGC
