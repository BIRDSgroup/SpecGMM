%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Main Wrappper Script -- DL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clearvars -except svc seedVals;
clc;

%% Dataset
% Mention the name of the dataset folders stored in the "SpecGMM/Database" folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Randhawa et al., 2019 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Primates', 'Protists', 'Fungi', 'Plants', 'Amphibians', 'Insects', ...
% 'threeClasses', 'Vertebrates', 'BacteriaTest', 'Birds_Fish_Mammals', ...
% 'Dengue', 'Mammalia', ...
% 'DomainToKingdom_Eukaryota', 'DomainToKingdom_Eukaryota_noProtists', 'KingdomToPhylum_Animalia', ...
% 'PhylumToSubphylum_Chordata', 'SubphylumToClass_Vertebrata', 'ClassToSubclass_Actinopterygii', ...
% 'SubclassToSuperorder_Neopterygii', 'SuperorderToOrder_Ostariophysi', 'OrderToFamily_Cypriniformes', ...
% 'FamilyToGenus_Cyprinidae', 'SubfamilyToGenus_Acheilognathinae'

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Randhawa et al., 2019 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Test1', 'Test2', 'Test3a', 'Test3b', 'Test4', 'Test5', 'Test6', ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16S rRNA datasets %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'KingdomToPhylum_Bacteria_16SrRNA', ...
% 'PhylumToClass_Firmicutes_16SrRNA', ...
% 'ClassToOrder_Bacilli_16SrRNA', ...
% 'OrderToFamily_Bacillales_16SrRNA', ...
% 'FamilyToGenus_Bacillaceae_16SrRNA', ...
% 'GenusToSpecies_Bacillus_16SrRNA'

% Datasets
dataSets = {
    'Primates'
};

% Classifiers
classifiers = {
	'LinearDiscriminant'
	'LinearSVM'
    'QuadraticSVM'
	'FineKNN'
	'SubspaceDiscriminant'
	'SubspaceKNN'
};

%% Arguments
args.seedClassifier = 15;   % Seed for classifiers
args.computeEmbeddings = 1; % Compute DNA embeddings

%% Parameters
parameters.folds = 4;       % Number of folds

%% To store the accuracies for each class
dataSetAccuracies.DL = struct();

%% Getting the current directory -- ~/path_to_folder/SpecGMM/
currentDirectory = pwd;

%% Add folders and subfolders in the current directory to the MATLAB path
% addpath(genpath('./'));

% Loop over datasets
for ds = 1:length(dataSets)

    % Display the dataset name
    disp(dataSets{ds});

    close all;
    clearvars -except parameters classifiers dataSets ds currentDirectory taxonomyTable dataSetAccuracies dataSetMetrics args featureNames svc seedVals;
    set(0, 'DefaultFigureRenderer', 'painters');

    % Random seed
    rng('default');

    % Change the directory to SpecGMM/
    cd(currentDirectory)

    % Dataset
    dataSet = dataSets{ds};

    %% Reading Data
    disp("Reading data...");
    [data.AcNmb, data.Sequences, data.numberOfClusters, data.clusterNames, data.pointsPerCluster, data.Sequences_withNs] = ReadData(dataSet);

    %% Getting the starting point of each class as data is read in 1D cell array
    StartingPoint = [];
    counter = 1;
    for cls = 1:data.numberOfClusters
        StartingPoint = [StartingPoint, counter];
        counter = counter + data.pointsPerCluster{cls};
    end
    data.StartingPoint = StartingPoint;

    % DNA Embedding
    DLEmbFileName = strcat(dataSet, '-', '-DLEmbPooled.mat');

    if args.computeEmbeddings == 1
        % Compute Embeddings using DNA-BERT
        fprintf('\nComputing Embeddings using DNA-BERT .... \n');
        SequenceDataPooled = ComputeEmbeddingsDNABERTS(data, dataSet);

        % Saving extracted embeddings
        save(DLEmbFileName, 'SequenceDataPooled', '-v7.3')
    else
        % Load pre-computed embeddings
        if isfile(DLEmbFileName)
            fprintf('\nLoading Pre-computed DL Embeddings .... \n');
            load(DLEmbFileName, 'SequenceDataPooled');
        else
            fprintf('\nPre-computed DL Embeddings not found. Computing DL Embeddings .... \n');
            % Compute Embeddings using DNA-BERT
            fprintf('\nComputing Embeddings using DNA-BERT .... \n');
            SequenceDataPooled = ComputeEmbeddingsDNABERTS(data, dataSet);

            % Saving extracted embeddings
            save(DLEmbFileName, 'SequenceDataPooled', '-v7.3')
        end
    end

    %% Train-test split for K-fold cross-validation
    counter = 1;
    for cls = 1:data.numberOfClusters

        % Starting point
        ppc = data.pointsPerCluster{cls};

        % Class sequences
        classSequences = data.Sequences(counter:counter+ppc-1);

        % Get the magnitude spectra data for the current class
        classData = SequenceDataPooled{cls};

        % Get Accession number
        classAccessionNumbers = data.AcNmb(counter:counter+ppc-1);

        % Sequene Lengths
		seqLengths = cellfun('length', data.Sequences(counter:counter+ppc-1));

        % Number of sequences in this class
        numDataPoints = numel(classSequences);

        % Create a cvpartition object for k-fold cross-validation
        cvp = cvpartition(numDataPoints, 'KFold', parameters.folds);

        % Loop through each fold
        for fold = 1:parameters.folds

            % Get the training indices for the current fold
            trainIndices = training(cvp, fold);
            
            % Get the testing indices for the current fold
            testIndices = test(cvp, fold);

            % Points per mixture component
			pointsPerClusterTrain{fold}{cls} = sum(trainIndices);
			pointsPerClusterTest{fold}{cls} = sum(testIndices);

            % AcNmb
            classAccessionNumbersTrain{fold}{cls} = classAccessionNumbers(trainIndices);
            classAccessionNumbersTest{fold}{cls} = classAccessionNumbers(testIndices);

            % Extract sequence lengths for current fold
            SequenceLengthsTrain{fold}{cls} = seqLengths(trainIndices);
            SequenceLengthsTest{fold}{cls} = seqLengths(testIndices);
            
            % Extract the magnitude spectra data for the current fold
            trainEmbeddings{fold}{cls} = classData(trainIndices,:);
            testEmbeddings{fold}{cls} = classData(testIndices,:);
        end
        counter = counter + ppc;
    end

    %% K-fold Cross Validation
    foldAccuracies.DLEmbeddings = zeros(parameters.folds, length(classifiers)); % DL Embeddings

    % To store fold metrics (5 metrics -- accuracy, precision, recall, specificity, f1_score)
    numMetrics=5;
    foldMetrics.DLEmbeddings = zeros(parameters.folds, numMetrics*length(classifiers)); % Supervectors

    % Loop over each fold
	for fold = 1:parameters.folds

        % Prepare training and test data
        trainData = zeros(sum([pointsPerClusterTrain{fold}{:}]), 768);
        testData = zeros(sum([pointsPerClusterTest{fold}{:}]), 768);

        % Generate labels for train and test sets
        trainLabels = [];
        testLabels = [];

        c1=1;
        for i1 = 1:data.numberOfClusters
            for j1 = 1:pointsPerClusterTrain{fold}{i1}
                trainData(c1,:) = trainEmbeddings{fold}{i1}(j1,:);
                trainLabels(c1,1) = i1;
                c1=c1+1;
            end
        end

        c1=1;
        for i1 = 1:data.numberOfClusters
            for j1 = 1:pointsPerClusterTest{fold}{i1}
                testData(c1,:) = testEmbeddings{fold}{i1}(j1,:);
                testLabels(c1,1) = i1;
                c1=c1+1;
            end
        end

        %% Classification Code
        fprintf('\nPerforming classification .... \n');

        % Direct classification using embeddings
        tic
            [accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, predLabels, conf_Matrix] = classificationCode(trainData, trainLabels, testData, testLabels, classifiers, args.seedClassifier);
            testLabelsFold{fold} = testLabels;
            predLabelsFold{fold} = predLabels;
            foldAccuracies.DLEmbeddings(fold,:) = accuracy';
            metrics = [];
            for clsfr = 1:length(classifiers)
                metrics = [metrics, accuracy(clsfr), precision(clsfr), recall(clsfr), specificity(clsfr), f1_score(clsfr)];
            end
            foldMetrics.DLEmbeddings(fold, :) = metrics;
        toc
    end

    disp(dataSet)

    % Average accuracy
    avgAcc = mean(foldAccuracies.DLEmbeddings);
    disp('Average Accuracies:')
    dataSetAccuracies.(dataSet) = avgAcc;
    disp(avgAcc)

    % Standard deviation of accuracies across folds
    stdAcc = std(foldAccuracies.DLEmbeddings, 1, 1);

    % Average metrics
    avgMetrics = mean(foldMetrics.DLEmbeddings);
    dataSetMetrics.(dataSet) = [avgMetrics, stdAcc];

end % End of datasets loop

% Writing the results to an Excel file
metricsConsolidated = zeros(length(dataSets), 6 * length(classifiers));
for ds = 1:length(dataSets)
    metricsConsolidated(ds, :) = dataSetMetrics.(dataSets{ds});
end
metricsConsolidated = [dataSets', num2cell(round(metricsConsolidated, 2))];
writecell(metricsConsolidated, "ConsolidatedMetrics-DNABERT-S.xlsx");

disp("Processing finished!")