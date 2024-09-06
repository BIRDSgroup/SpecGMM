# Supplementary Data Files

## File D1: Dataset Details

The file contains the following details of the datasets used in the study.

- Number and names of categories for each dataset
- Number of sequences per category
- Proportion of each category and the analytically computed chance accuracies using the proportions
- Descriptive statistics for lengths of sequences of each category across all the datasets — maximum, minimum, mean, mode, std. dev., mean absolute deviation, and median absolute deviation

## File D2: 16S-ITGDB Taxonomy and HVR Information

The file contains the taxonomy information (Kingdom to Species labels) and the start and end positions of the 16S rRNA HVRs
obtained using the QIIME2 analysis.

## File D3: Performance — Baseline vs. SpecGMM

The file contains detailed classification results for all the datasets analysed in the study for Linear Discriminant (LD), Linear SVM (LSVM), Quadratic SVM (QSVM), FineKNN, Subspace Discriminant, and Subspace KNN classifiers comparing baseline and SpecGMM methods. The following performance metrics, computed over four folds, are available in the file:

- average accuracy
- standard deviation
- average weighted precision
- average weighted recall
- average weighted specificity
- average weighted F1-score

The file also contains results for the baseline and the SpecGMM methods after performing homology reduction on our datasets using the GraphPart algorithm for different threshold values and different numbers of partitions. The LD and LSVM classifiers were used in the analysis.

Lastly, the file contains a comparative analysis of SpecGMM and DNABERT-S models.

## File D4: Comparative analysis of numerical representations

The file reports SpecGMM's performance on various datasets from the baseline study using different numerical representations.
