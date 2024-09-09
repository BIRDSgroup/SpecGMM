%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Main Wrappper Script -- Baseline Method
% Adapted from Randhawa et al. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;
set(0, 'DefaultFigureRenderer', 'painters');

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

dataSets = {
'Primates'
};

% To store the accuracies for each class
dataSetAccuracies = struct();
dataSetMetrics = struct();

%% Define the classifiers to be used ('LinearDiscriminant', 'LinearSVM')
classifiers = {
	'LinearDiscriminant'
	'LinearSVM'
    'QuadraticSVM'
	'FineKNN'
	'SubspaceDiscriminant'
	'SubspaceKNN'
};

% Getting the current directory
currentDirectory = pwd;

%% Getting parallel pool
LASTN = maxNumCompThreads(12); % Max number of computational threads
p = gcp('nocreate');
if isempty(p)           % There is no parallel pool
    parpool(12, 'IdleTimeout', Inf);
else                    % There is a parallel pool of <p.NumWorkers> workers
    disp(strcat("Parpool: ", string(p.NumWorkers)))
end

%% Add folders and subfolders in the current directory to the path
% addpath(genpath('./'))

%% Loop over datasets
for ds = 1:length(dataSets)

    close all;
    clearvars -except classifiers dataSets dataSetAccuracies dataSetMetrics ds currentDirectory p taxonomyTable;
    set(0, 'DefaultFigureRenderer', 'painters');

    % Random seed
    rng('default');

    % Change the directory to SpecGMM/
    cd(currentDirectory)

    %% Dataset
    dataSet = dataSets{ds};         

    lenNormType = 'med';            % length-normalization
    distType = 'cor';               % Distance type, Euclidean -- 'euc', PCC -- 'cor'
    folds = 4;                      % Number of folds

    disp(strcat('------- ', dataSet, ' ------- ', lenNormType, ' ------- ', distType, '------- '))

    %% Reading Data
    if isfile(strcat(dataSet,'.mat'))
        disp('Retreiving the preloaded dataset...')
        load(strcat(dataSet,'.mat'))
        if exist('Seq', 'var') == 1
            Sequences = Seq;
        end
    else
        % Reading data from fasta files and saving the same
        fprintf('Reading sequences .... \n');
        [AcNmb, Sequences, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet);
        save(strcat(dataSet,'.mat'), 'AcNmb', 'Sequences', 'numberOfClusters', 'clusterNames', 'pointsPerCluster');
    end

    totalSeq = length(Sequences);

    % Calculate length stats
    [maxLen, minLen, meanLen, medLen] = lengthCalc(Sequences);

    %% Length normalization
    % numerical sequences, length normalized by median length
    % code can be modified to use other length stats for length normalization
    if lenNormType == 'max'
        mLen = maxLen;
    elseif lenNormType == 'med'
        mLen = medLen;
    end

    % Preallocate cell arrays
    nmValSH=cell(1,totalSeq);
    f=cell(1,totalSeq);
    lg=cell(1,totalSeq);

    tic
        fprintf('Generating numerical sequences, applying DFT, computing magnitude spectra .... \n');
   
        parfor a = 1:totalSeq
            % using Purin/Pyramidine representation by default
            % change method call for other representations
            
            ns = numMapping(Sequences{a}, 'PP');            % Mapping
    
            I = mLen-length(ns);
    
            % Padding
            if(I>0)
                nsTemp = wextend('1','asym',ns,I);
                nsNew = nsTemp((I+1):length(nsTemp));
            elseif(I<0)
                nsNew=ns(1:mLen);
            else
                nsNew = ns;
            end
    
            nmValSH{a} = nsNew;
            f{a} = fft(nsNew);              % Fourier transform
            lg{a} = abs(f{a});              % Magnitude spectra
        end
    toc

    %% Classwise data
    SequenceDataPooled={};
    counter=1;
    for cls = 1:numberOfClusters
        SequenceDataPooled{cls} = lg(counter:counter+pointsPerCluster{cls}-1);
        counter = counter + pointsPerCluster{cls};
    end

    %% Train-test split for K-fold cross validation
    for cls = 1:numberOfClusters

        % Get the data for the current class
        classData = SequenceDataPooled{cls};

        % Calculate the number of data points in this class
        numDataPoints = numel(classData);

        % Create a cvpartition object for k-fold cross-validation
        cvp = cvpartition(numDataPoints, 'KFold', folds);

        % Loop through each fold
        for fold = 1:folds
            % Get the training indices for the current fold
            trainIndices = training(cvp, fold);
            
            % Get the testing indices for the current fold
            testIndices = test(cvp, fold);
            
            % Extract the data for the current fold
            SequenceDataPooledTrain{fold}{cls} = classData(trainIndices);
            SequenceDataPooledTest{fold}{cls} = classData(testIndices);
        end
    end

    % To store fold accuracies
    foldAccuracies = [];

    %% K-fold Cross Validation
    for fold = 1:folds
        fprintf('\nFold %d................................................. \n', fold);
        disp(datetime);
        % Random seed
        rng(15,'twister');

        %% Points per cluster in the train and test data in eac fold
        for cls=1:numberOfClusters
            pointsPerClusterTrain{cls} = size(SequenceDataPooledTrain{fold}{cls}, 2);
        end

        for cls=1:numberOfClusters
            pointsPerClusterTest{cls} = size(SequenceDataPooledTest{fold}{cls}, 2);
        end

        %% Stacking train and test
        
        VectorsMatrixTrain = zeros(sum([pointsPerClusterTrain{:}]), mLen); % Train Supervectors
        VectorsMatrixTest = zeros(sum([pointsPerClusterTest{:}]), mLen); % Test Supervectors
        
        trainLabels = [];
        testLabels = [];
        
        counterTrain=1;
        counterTest=1;
        
        for i=1:numberOfClusters
            for j=1:pointsPerClusterTrain{i}
                VectorsMatrixTrain(counterTrain, :) = SequenceDataPooledTrain{fold}{i}{j};
                trainLabels=[trainLabels; i];
                counterTrain = counterTrain + 1;
            end
        
            for j=1:pointsPerClusterTest{i}
                VectorsMatrixTest(counterTest, :) = SequenceDataPooledTest{fold}{i}{j};
                testLabels=[testLabels; i];
                counterTest = counterTest + 1;
            end
        end
        
        %% distance calculation by Pearson correlation coefficient
        % change distType from 'cor' to 'euc' for Euclidean
        totalTrainSeq = sum([pointsPerClusterTrain{:}]);
        lg=cell(1,totalTrainSeq);
        
        for i = 1:totalTrainSeq
            lg{i} = VectorsMatrixTrain(i, :);
        end
        
        tic
            fprintf('\nComputing Distance matrix .... \n');
            fm=cell2mat(lg(:));
            disMat = f_dis(fm,distType,0,1);
        toc
        
        %% Test
        tic
            fprintf('\nGenerating test vectors .... \n');

            batchSize = 10; % Adjust based on memory constraints and performance
            numTrainSamples = size(VectorsMatrixTrain, 1);
            numTestSamples = size(VectorsMatrixTest, 1);
            numBatches = ceil(numTestSamples / batchSize);
            
            % Preallocate cell array to store results from each batch
            testMatrixBatches = cell(1, numBatches);
            
            for batchIdx = 1:numBatches
                batchStart = (batchIdx - 1) * batchSize + 1;
                batchEnd = min(batchIdx * batchSize, numTestSamples);
                testBatch = VectorsMatrixTest(batchStart:batchEnd, :);
                
                % Compute correlations for the current batch
                correlationsBatch = corr(testBatch', VectorsMatrixTrain');
                
                % Calculate the dissimilarity for the current batch and store in the cell array
                testMatrixBatches{batchIdx} = 0.5 * (1 - correlationsBatch');
            end
            
            % Combine the results from each batch into the final testMatrix
            testMatrix = zeros(numTrainSamples, numTestSamples);
            for batchIdx = 1:numBatches
                batchStart = (batchIdx - 1) * batchSize + 1;
                batchEnd = min(batchIdx * batchSize, numTestSamples);
                
                % Assign the results from each batch to the corresponding columns in testMatrix
                testMatrix(:, batchStart:batchEnd) = testMatrixBatches{batchIdx};
            end
        toc
        
        
        %% Classification Code
        cd(pwd)
        
        fprintf('\nPerforming classification .... \n');

        % Classification based on pairwise distance matrix columns (random seed = 15)
        tic
            [accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, conf_Matrix] = classificationCode(disMat, trainLabels, testMatrix', testLabels, classifiers, 15);
        toc
        foldAccuracies(fold,:) = accuracy';
        metrics = [];
        for clsfr=1:length(classifiers)
            metrics = [metrics, accuracy(clsfr), precision(clsfr), recall(clsfr), specificity(clsfr), f1_score(clsfr)];
        end
        foldMetrics(fold,:) = metrics;
    end


    % Average accuracy based on pairwise distance matrix
    avgAcc = mean(foldAccuracies);
    fprintf('\nAverage Accuracies for %s .... \n', dataSet);
    dataSetAccuracies.(dataSet) = avgAcc;
    disp(avgAcc)

    % Standard deviation of accuracies based on Supervectors across folds
    stdAcc = std(foldAccuracies, 1, 1);

    % Average metrics based on Supervectors
    avgMetrics = mean(foldMetrics);
    dataSetMetrics.(dataSet) = [avgMetrics, stdAcc];
end

% Initialize the metricsConsolidated matrix with zeros
metricsConsolidated = zeros(length(dataSets), 6 * length(classifiers));

% Consolidate all the metrics in the desired order
for ds = 1:length(dataSets)
    metrics = dataSetMetrics.(dataSets{ds});
    
    % Number of metrics excluding stdDev
    numBasicMetrics = 5;
    
    for i = 1:length(classifiers)
        % Extract indices for the basic metrics
        basicMetricStartIdx = (i - 1) * numBasicMetrics + 1;
        basicMetricEndIdx = basicMetricStartIdx + numBasicMetrics - 1;
        
        % Extract the basic metrics (accuracy, precision, recall, specificity, f1_score)
        basicMetrics = metrics(basicMetricStartIdx:basicMetricEndIdx);
        
        % Extract the standard deviation, which is placed after all basic metrics
        stdDevIdx = length(classifiers) * numBasicMetrics + i;
        stdDev = metrics(stdDevIdx);
        
        % Rearrange and store metrics in the desired order
        startIdx = (i - 1) * 6 + 1;
        metricsConsolidated(ds, startIdx:startIdx+5) = [basicMetrics(1), stdDev, basicMetrics(2:end)];
    end
end

% Generate headers in the new order
headers = {'Dataset_name'};
for i = 1:length(classifiers)
    classifierName = classifiers{i};
    headers = [headers, ...
            strcat(classifierName, '_Accuracy'), ...
            strcat(classifierName, '_StdDev'), ...
            strcat(classifierName, '_Precision'), ...
            strcat(classifierName, '_Recall'), ...
            strcat(classifierName, '_Specificity'), ...
            strcat(classifierName, '_F1Score')];
end

% Combine headers with the data
metricsConsolidated = [headers; dataSets', num2cell(round(metricsConsolidated, 2))];

% Write to Excel file
writecell(metricsConsolidated, "PP-ConsolidatedMetrics-Baseline-LargeDatasets.xlsx");



disp("Processing finished!")
