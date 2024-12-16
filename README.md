This repository contains a set of scripts designed to benchmark DeepFRI, a deep learning-based tool for functional annotation, and to evaluate its performance alongside other methods such as Blast2GO. The pipeline has been tailored specifically for small peptides, with a focus on removing biases and ensuring high-quality evaluations.

**Tools and Scripts**

**1. SimRel.py**
This script implements the SimRel scoring to evaluate the similarity between predicted and known Gene Ontology (GO) terms. It calculates scores for three GO subontologies (biological process, molecular function, cellular component) . Note that the input file for predicted go-terms is named go_terms_DeepFRI.tsv in this instance. This file can be replaced by any other file containing predicted GO-terms as we also did when running SimRel.py on other prediction tools. We would like to highlight the importance of the formatting of the input file: one row per peptide, with its protein id listed in the first column and all its predicted GO-terms listed in the second column, separated by commas. To see an example: al_go_terms.tsv is also fromatted in this way


**2. DeepFRI_RUN.py**
Automates the process of running DeepFRI. It downloads AlphaFold PDB files, modifies headers for compatibility, and predicts GO annotations for peptides across all three GO subontologies. Note that the results file must be reformatted before it can be given as input to SimRel.py

**3. make_datatfile.py**
Reformats the results.tsv obtained from DeepFRI_RUN.py and generates go_terms_DeepFRI.tsv, which can be used as input for SimRel.py

**4. filtering_DeepFRI_train_val.sh**
Filters out peptides from the test dataset that are highly similar to the training and validation data used in DeepFRI. This ensures unbiased evaluation by creating a final dataset free of overlapping sequences.

**5. boxplot.py**
Generates boxplots to visualize the distribution of SimRel scores across GO subontologies, providing a summary of prediction performance.

**Workflow**

**Data Preparation**: Use filtering_DeepFRI_train_val.sh to remove peptides with high similarity to DeepFRI’s training/validation data.

**Prediction**: Run DeepFRI_RUN.py to predict GO annotations for the filtered peptide dataset.

**Reformatting**: use make_datafile.py to reformat DeepFRI_RUN.py's output

**Scoring**: Evaluate the predictions using SimRel.py, which generates similarity scores.

**Visualization**: Use boxplot.py to create graphical summaries of the scoring results.

**Reuse SimRel.py**: score and visualize predictions from other tools for comparative analysis. Mind the formatting of the input file

**Requirements**

Python (3.7+)

BLAST+ tools

Required Python libraries: pandas, numpy, matplotlib, requests, goatools

DeepFRI setup with predict.py, can be obtained from the DeepFRI github: https://github.com/flatironinstitute/DeepFRI
