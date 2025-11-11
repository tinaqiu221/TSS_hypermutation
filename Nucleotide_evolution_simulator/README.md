# Data processing

## normalization.py
    Normalizes mutation rate data by dividing values by a baseline calculated from intergenic regions.

## curve_fitting.py
    Fits mathematical functions to mutation rate patterns in specific genomic regions to model and analyze mutation distribution trends.
    TSS (Transcription Start Sites) - fitted with double Gaussian
    EIB (Exon-Intron Boundaries) - fitted with exponential/linear
    Exon regions - fitted with exponential/linear

## mr_cf_residual.py
    Compares actual mutation rate data with fitted mathematical models and analyzes the residuals (differences) between them.

# Simulation

## sequences_generations_mr.py
	Generates the sequences to be mutated in the simulations, by implementing the equilibrium probability for each nucleotide per position in the sequence. (TSS, EIB, background B, hotspot B factors taken into consideration to initialize the nucleotide sequences)

## mutation_rate.py
	Gathers functions that make the mutation rate matrix evolve (hotspot_B, background_B, etc)

## simulator_outputs_graph.py
	Creates an output file per combination of factors (additional_B, background_B, current_T), gathering data from each generated sequence. One summary file per nucleotide is generated, containing only one row (the average mutation rates along the sequences).

## simulator_gc_tss_mr_graph.py
    Simulates the evolution of sequences by varying the factors: additional_B, background_B, and current_T.
	Implementation of TSS correction, EIB correction as a distribution.
	Creates output graphs

To run the codes to get the data + graphs : 
Example in terminal : python simulator_gc_tss_mr_graph.py --nseq 50 --bg_B 0.8 --add_B 1.3

# Performance analysis

## score_function.py
    Compares simulated nucleotide sequence data against real nucleotide distribution patterns and scores how closely they match.
    Calculates a similarity score using a weighted squared difference formula.
    Creates scatter plots visualizing the scores

## AIC_score.py
    Compares different simulation models (varying in parameters like EIB and TSS residuals) to determine which best matches observed nucleotide content (A, C, G bases).
    AIC Calculation:
        Computes mean and variance between observed and simulated values
        Calculates Gaussian log-likelihood for each nucleotide base
        Computes AIC score using: AIC = 2k - 2ln(L) where k is number of parameters
    Model Comparison:
        Finds the model with minimum AIC score
        Computes differences from the best model
        Visualizes results with colored scatter plots (green = best models, red = poor models)

## heatmap_score.py
    Creates customizable heatmap visualizations from score data.
    Default mode: Shows hotspot parameters (peak vs width)
    gBGC mode: Shows additional B vs background B parameters
    B_Time mode: Shows additional B vs time off parameters
