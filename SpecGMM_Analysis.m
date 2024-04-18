%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Main Wrappper Script -- SpecGMM Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clearvars -except svc seedVals;
clc;

%% Dataset
% Mention the name of the dataset folders stored in the "SpecGMM/Database" folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Randhawa et al., 2019 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Primates', 'Protists', 'Fungi', 'Plants', 'Amphibians', 'Insects', ...
% 'threeClasses', 'Vertebrates', 'BacteriaTest', 'Birds_Fish_Mammals', ...
% 'Dengue', 'Mammalia', 'NewVertSequences', ...
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

dataSets = {'GenusToSpecies_Bacillus'};

%% Arguments
% If the dataset is one of the above mentioned 16S rRNA datasets and you want to do HVR analysis, set the following variable to 1
args.analysis16SrRNA = 1;

if args.analysis16SrRNA == 1
    % Load Taxonomy Table which contains the taxonomy and the HVR start-end information obtained using QIIME2
    load('DataBase/taxonomyTable.mat');
end

args.magSpectra = 1;        % Window-based magnitude spectra extraction: 1--extract, 0--load pre-computed spectra
args.assignHVR = 1;         % Assign HVR to each window based on maximum overlap: 1--assign, 0--load pre-assigned HVR
args.trainUBM = 1;          % Build UBM-GMM on the training data: 1--train, 0--load pre-built UBM-GMM
args.seedUBM = 15;          % Seed for training UBM-GMM
args.seedClassifier = 15;   % Seed for classifiers

%% Parameters
parameters.winLen = 63;                                 % Window length
parameters.fftOrder= 2^nextpow2(parameters.winLen);     % FFT Order
parameters.winShift = 9;                                % Window shift
parameters.K = 5;                          % Number of mixture components in UBM-GMM
parameters.folds = 4;                      % Number of parameters.folds
parameters.magSpectrumStart = 1;           % 1--DC
parameters.magSpectrumEnd = parameters.fftOrder/2;
parameters.magSpectrumDimension = parameters.magSpectrumEnd-parameters.magSpectrumStart+1;

%% Define the classifiers to be used ('LinearDiscriminant', 'LinearSVM')
classifiers = {
    'LinearDiscriminant';
    'LinearSVM';
};

%% 16S rRNA analysis
featureNames = {'SV', 'MW', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7'};
numFeatures = length(featureNames);
acVals = struct();
for clsfr = 1:length(classifiers)
    classifierName = classifiers{clsfr};
    acVals.(classifierName) = struct();
end
for fold = 1:parameters.folds
        foldName = ['fold', num2str(fold)]; % Convert fold number to a string for field name
        acVals.(classifierName).(foldName) = zeros(length(dataSets), numFeatures);
end

%% To store the accuracies for each class
dataSetAccuracies.SV = struct();
dataSetMetrics.SV = struct();
if args.analysis16SrRNA == 1
    dataSetAccuracies.MW = struct();
    dataSetAccuracies.V2 = struct();
    dataSetAccuracies.V3 = struct();
    dataSetAccuracies.V4 = struct();
    dataSetAccuracies.V5 = struct();
    dataSetAccuracies.V6 = struct();
    dataSetAccuracies.V7 = struct();
end

%% Getting the current directory -- ~/path_to_folder/SpecGMM/
currentDirectory = pwd;

%% Add folders and subfolders in the current directory to the MATLAB path
addpath(genpath('./'))

%% Loop over datasets
for ds = 1:length(dataSets)

    % Display the dataset name
    disp(dataSets{ds});

    close all;
    clearvars -except parameters classifiers dataSets ds currentDirectory taxonomyTable dataSetAccuracies dataSetMetrics args featureNames acVals svc seedVals;
    set(0, 'DefaultFigureRenderer', 'painters');

    % Random seed
    rng('default');

    % Change the directory to SpecGMM/
    cd(currentDirectory)

    % Dataset
    dataSet = dataSets{ds};

    %% Reading Data
    [data.AcNmb, data.Sequences, data.numberOfClusters, data.clusterNames, data.pointsPerCluster] = ReadData(dataSet);

    %% Getting the starting point of each class as data is read in 1D cell array
    StartingPoint = [];
    counter = 1;
    for cls = 1:data.numberOfClusters
        StartingPoint = [StartingPoint, counter];
        counter = counter + data.pointsPerCluster{cls};
    end
    data.StartingPoint = StartingPoint;

    %% Compute length statistics -- Each row: clusterName, pointsPerCluster, maxLen, minLen, meanLen, medLen, modLen, stdDev, meanAbsDev, medAbsDev
    stats = LengthStatistics(data);

    %% Window-based magnitude spectra
    if args.magSpectra == 1
        fprintf('\nComputing Magnitude Spectra .... \n');
        % Extract magnitude spectra
        SequenceDataPooled = ComputeMagnitudeSpectra(parameters, data);

        % Saving the extracted magnitude spectra
        save(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-MagSpecPooled.mat'), 'SequenceDataPooled', '-v7.3');
    else
        fprintf('\nLoading Pre-computed Magnitude Spectra .... \n');
        % Loading pre-computed magnitude spectra
        load(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-MagSpecPooled.mat'), 'SequenceDataPooled');
    end

    if args.analysis16SrRNA == 1
        %% Get the HVR-assignment to each window
        if args.assignHVR == 1
            fprintf('\nAssigning HVRs to windows .... \n');
            % Assigning HVRs
            assignedHVRs_DataPooled = AssignHVR(parameters, data, taxonomyTable);
            
            % Saving the assigned HVRs
            save(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-assignedHVRs_DataPooled.mat'), 'assignedHVRs_DataPooled', '-v7.3');
        else
            fprintf('\nLoading pre-assigned HVRs for windows .... \n');
            % Loading the pre-assigned HVRs
            load(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-assignedHVRs_DataPooled.mat'), 'assignedHVRs_DataPooled');
        end 
    end  

    %% Train-test split for K-fold cross validation
    for cls = 1:data.numberOfClusters

        % Get the magnitude spectra data for the current class
        classData = SequenceDataPooled{cls};

        % Get the assigned HVR data for the current class
        if args.analysis16SrRNA == 1
            classHVR = assignedHVRs_DataPooled{cls};
        end

        % Number of sequences in this class
        numDataPoints = numel(classData);

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
            
            % Extract the magnitude spectra data for the current fold
            SequenceDataPooledTrain{fold}{cls} = classData(trainIndices);
            SequenceDataPooledTest{fold}{cls} = classData(testIndices);

            if args.analysis16SrRNA == 1
                % Extract the HVR assignment data for the current fold
                assignedHVRs_DataPooledTrain{fold}{cls} = classHVR(trainIndices);
                assignedHVRs_DataPooledTest{fold}{cls} = classHVR(testIndices);
            end
        end
    end


    %% K-fold Cross Validation

    % To store fold accuracies
    foldAccuracies.SV = zeros(parameters.folds, length(classifiers)); % Supervectors
    if args.analysis16SrRNA == 1
        foldAccuracies.MW = zeros(parameters.folds, length(classifiers)); % Mixture weights
        foldAccuracies.V2 = zeros(parameters.folds, length(classifiers)); % V2 Activations
        foldAccuracies.V3 = zeros(parameters.folds, length(classifiers)); % V3 Activations
        foldAccuracies.V4 = zeros(parameters.folds, length(classifiers)); % V4 Activations
        foldAccuracies.V5 = zeros(parameters.folds, length(classifiers)); % V5 Activations
        foldAccuracies.V6 = zeros(parameters.folds, length(classifiers)); % V6 Activations
        foldAccuracies.V7 = zeros(parameters.folds, length(classifiers)); % V7 Activations
    end

    % To store fold metrics (5 metrics -- accuracy, precision, recall, specificity, f1_score)
    numMetrics=5;
    foldMetrics.SV = zeros(parameters.folds, numMetrics*length(classifiers)); % Supervectors
    if args.analysis16SrRNA == 1
        foldMetrics.MW = zeros(parameters.folds, numMetrics*length(classifiers)); % Mixture weights
        foldMetrics.V2 = zeros(parameters.folds, numMetrics*length(classifiers)); % V2 Activations
        foldMetrics.V3 = zeros(parameters.folds, numMetrics*length(classifiers)); % V3 Activations
        foldMetrics.V4 = zeros(parameters.folds, numMetrics*length(classifiers)); % V4 Activations
        foldMetrics.V5 = zeros(parameters.folds, numMetrics*length(classifiers)); % V5 Activations
        foldMetrics.V6 = zeros(parameters.folds, numMetrics*length(classifiers)); % V6 Activations
        foldMetrics.V7 = zeros(parameters.folds, numMetrics*length(classifiers)); % V7 Activations
    end

    % Loop over each fold
    for fold = 1:parameters.folds
    
        fprintf('\nFold %d................................................. \n', fold);
	    disp(datetime);

        % Random seed
        rng(args.seedUBM,'twister');

        % Height of matrix to store all magnitude spectra in training sequences
        matHeight=0;
        for cls =1:data.numberOfClusters
            for a = 1:pointsPerClusterTrain{fold}{cls}
                nWin = size(SequenceDataPooledTrain{fold}{cls}{a}, 1);
                matHeight=matHeight+nWin;
            end
        end
        
        % Stacking data from all the classes to build a UBM-GMM
        StackedMagnitudeSpectra = zeros(matHeight, parameters.magSpectrumDimension); % To stack the data
        counter = 1;
        for cls =1:data.numberOfClusters
            for a = 1:pointsPerClusterTrain{fold}{cls}
                seqData = SequenceDataPooledTrain{fold}{cls}{a};
                nWin = size(seqData,1);
                StackedMagnitudeSpectra(counter:counter+nWin-1, :) = SequenceDataPooledTrain{fold}{cls}{a};
                counter=counter+nWin;
            end
        end
        
        % UBM-GMM Training
        if args.trainUBM == 1
            fprintf('\nUBM-GMM Training .... \n');
            tic
            options = statset('MaxIter',1000);
            UBMGMM = fitgmdist(StackedMagnitudeSpectra,parameters.K,'CovarianceType','diagonal','Options',options);
            toc
            
            save(strcat(dataSet, '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(fold), '-', string(args.seedUBM), '-UBMGMM.mat'), 'UBMGMM', '-v7.3');
        else
            fprintf('\nLoading Pre-trained UBM-GMM .... \n');
            load(strcat(dataSet, '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(fold), '-', string(args.seedUBM), '-UBMGMM.mat'), 'UBMGMM');
        end

        % MAP Adaptation and Average Posteriors for each HVR
        fprintf('\nMAP Adaptation and Average Posteriors for HVRs .... \n');

        if args.analysis16SrRNA == 1
            % Call the function for train data
            [SuperVectorsMatrixTrain, trainLabels, mixtureActivationsPerHVRPerSequencePerClusterTrain] = MAPAdaptationAndAveragePosteriorsHVR(SequenceDataPooledTrain, assignedHVRs_DataPooledTrain, pointsPerClusterTrain, UBMGMM, data, parameters, fold);

            % Call the function for test data
            [SuperVectorsMatrixTest, testLabels, mixtureActivationsPerHVRPerSequencePerClusterTest] = MAPAdaptationAndAveragePosteriorsHVR(SequenceDataPooledTest, assignedHVRs_DataPooledTest, pointsPerClusterTest, UBMGMM, data, parameters, fold);
        else
             % Call the function for train data
            [SuperVectorsMatrixTrain, trainLabels, mixtureComponentWeightsTrain] = MAPAdaptation(SequenceDataPooledTrain, pointsPerClusterTrain, UBMGMM, data, parameters, fold);

            % Call the function for test data
            [SuperVectorsMatrixTest, testLabels, mixtureComponentWeightsTest] = MAPAdaptation(SequenceDataPooledTest, pointsPerClusterTest, UBMGMM, data, parameters, fold);
        end

        if args.analysis16SrRNA == 1
            % Activations for each HVR across all the train seqeunces
            hvrActivationsTrain = {};
            for hvrIndex = 1:8
                vrActivations = [];
                for i = 1:data.numberOfClusters
                    vrActivations = [vrActivations; mixtureActivationsPerHVRPerSequencePerClusterTrain{i}{hvrIndex}];
                end
                hvrActivationsTrain{hvrIndex} = vrActivations;
            end
    
            % Activations for each HVR across all the test seqeunces
            hvrActivationsTest = {};
            for hvrIndex = 1:8
                vrActivations = [];
                for i = 1:data.numberOfClusters
                    vrActivations = [vrActivations; mixtureActivationsPerHVRPerSequencePerClusterTest{i}{hvrIndex}];
                end
                hvrActivationsTest{hvrIndex} = vrActivations;
            end        
        end
        

        %% Classification Code
        cd(pwd)
        
        fprintf('\nPerforming classification .... \n');

        % Classification based on Supervectors
        fprintf('\nFeatures: Supervectors .... \n');
        tic
            [accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, conf_Matrix] = classificationCode(SuperVectorsMatrixTrain, trainLabels, SuperVectorsMatrixTest, testLabels, classifiers, args.seedClassifier);
            foldAccuracies.SV(fold,:) = accuracy';
            metrics = [];
            for clsfr=1:length(classifiers)
                metrics = [metrics, accuracy(clsfr), precision(clsfr), recall(clsfr), specificity(clsfr), f1_score(clsfr)];
            end
            foldMetrics.SV(fold,:) = metrics;
        toc

        % Classification based on avaerage posteriors from different HVRs
        if args.analysis16SrRNA == 1
            featureNames = {'MW', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7'};
            numFeatures = length(featureNames);

            for featureIndex = 1:numFeatures
                featureName = featureNames{featureIndex};
                
                fprintf('\nFeatures: %s .... \n', featureName);
                tic
                    [accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, conf_Matrix] = classificationCode(hvrActivationsTrain{featureIndex}, trainLabels, hvrActivationsTest{featureIndex}, testLabels, classifiers, args.seedClassifier);
                    foldAccuracies.(featureName)(fold, :) = accuracy';
                    metrics = [];
                    for clsfr=1:length(classifiers)
                        metrics = [metrics, accuracy(clsfr), precision(clsfr), recall(clsfr), specificity(clsfr), f1_score(clsfr)];
                    end
                    foldMetrics.(featureName)(fold,:) = metrics;
                toc
            end
        end
    end

    disp(dataSet)

    % If 16S rRNA Analysis
    if args.analysis16SrRNA == 1
        featureNames = {'SV', 'MW', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7'};
        numFeatures = length(featureNames);

        for featureIndex = 1:numFeatures
            featureName = featureNames{featureIndex};
            
            % Calculate average accuracy
            avgAcc = mean(foldAccuracies.(featureName));
            disp(['Average Accuracies: ', featureName, ':'])
            dataSetAccuracies.(featureName).(dataSet) = avgAcc;
            disp(avgAcc)

            % Average metrics based on feature
            avgMetrics = mean(foldMetrics.(featureName));
            dataSetMetrics.(featureName).(dataSet) = avgMetrics;
        end

        % for 16S rRNA plot
        if args.analysis16SrRNA == 1
            for clsfr = 1:length(classifiers)
                classifierName = classifiers{clsfr};        % Assuming classifiers is a cell array of strings
                for fold = 1:parameters.folds
                    foldName = ['fold', num2str(fold)];     % Convert fold number to a string for field name            
                    for featureIndex = 1:numFeatures
                        featureName = featureNames{featureIndex};              
                        acVals.(classifierName).(foldName)(ds, featureIndex) = foldAccuracies.(featureName)(fold, clsfr);
                    end
                end
            end
        end
    % If not 16S rRNA Analysis
    else
        % Average accuracy based on Supervectors
        avgAcc = mean(foldAccuracies.SV);
        disp('Average Accuracies: SuperVectors:')
        dataSetAccuracies.SV.(dataSet) = avgAcc;
        disp(avgAcc)

        % Standard deviation of accuracies based on Supervectors across folds
        stdAcc = std(foldAccuracies.SV, 1, 1);

        % Average metrics based on Supervectors
        avgMetrics = mean(foldMetrics.SV);
        dataSetMetrics.SV.(dataSet) = [avgMetrics, stdAcc];

    end
end % End of dataSets loop

% Writing the results to an excel file for 16S rRNA analysis
if args.analysis16SrRNA == 1
    excelFileName = strcat('16SrRNA', '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(args.seedUBM), '-', string(args.seedClassifier), '.xlsx')

    classifiersNames = fieldnames(acVals); % Get all classifier names

    for clsfrIdx = 1:length(classifiersNames)
        classifierName = classifiersNames{clsfrIdx};
        foldNames = fieldnames(acVals.(classifierName)); % Get all fold names for this classifier
        dataForExcel = {}; % Initialize cell array to hold data to be written to Excel
        header = ['DataSet', featureNames]; % Assuming your first column will be the DataSet index, followed by feature names
        
        for foldIdx = 1:length(foldNames)
            foldName = foldNames{foldIdx};
            foldData = acVals.(classifierName).(foldName); % Extract fold data
            foldDataWithIndex = [(1:size(foldData, 1))', foldData]; % Add an index column if needed
            % Create a header for this fold
            foldHeader = [{['Fold ' num2str(foldIdx)]}, cell(1, length(featureNames))]; % Adjust as per your layout preferences
            % Append fold header and data to the dataForExcel
            dataForExcel = [dataForExcel; foldHeader; header; num2cell(foldDataWithIndex)]; % Convert numeric data to cell for writecell
            if foldIdx < length(foldNames) % Add an empty row between folds for clarity
                dataForExcel = [dataForExcel; cell(1, size(dataForExcel, 2))];
            end
        end
        
        % Write to a separate sheet for each classifier
        writecell(dataForExcel, excelFileName, 'Sheet', classifierName);
    end
end

% Consolidate all the metrics -- (5 metrics -- accuracy, precision, recall, specificity, f1_score) + std. dev
metricsConsolidated = zeros(length(dataSets), 6*length(classifiers));
for ds=1:length(dataSets)
    metricsConsolidated(ds,:) = dataSetMetrics.SV.(dataSets{ds});
end
metricsConsolidated = [dataSets', num2cell(round(metricsConsolidated,2))];
writecell(metricsConsolidated, "ConsolidatedMetrics.xlsx");

disp("Processing finished!")
