%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Main Wrappper Script -- Baseline GraphPart
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
    'LinearDiscriminant';
    'LinearSVM';
};

% Getting the current directory
currentDirectory = pwd;

%% Getting parallel pool
LASTN = maxNumCompThreads(4); % Max number of computational threads
p = gcp('nocreate');
if isempty(p)           % There is no parallel pool
    parpool(4, 'IdleTimeout', Inf);
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
    folds = 2;                      % Number of folds

    disp(strcat('------- ', dataSet, ' ------- ', lenNormType, ' ------- ', distType, '------- '))

    % Define the path to your CSV file
    csvFile = strcat("graphpart_", dataSet, "-Merged-Labeled_part4_thresh0.9.csv");

    % Load the cluster data from the CSV
    opts = detectImportOptions(csvFile, 'VariableNamingRule', 'preserve');

    % Set the first two columns as strings
    opts = setvartype(opts, {'AC', 'priority'}, 'string');

    % Check if the 'between_connectivity' column exists in the file
    fileContents = readtable(csvFile, 'ReadVariableNames', true);  % Read the CSV without specifying types initially

    % Adjust import options based on the presence of 'between_connectivity'
    if ismember('between_connectivity', fileContents.Properties.VariableNames)
        % If the column exists, set its type as double along with other columns
        opts = setvartype(opts, {'label-val', 'between_connectivity', 'cluster'}, 'double');
    else
        % If the column doesn't exist, only set 'label-val' and 'cluster' as double
        opts = setvartype(opts, {'label-val', 'cluster'}, 'double');
    end

    % Read the table with the adjusted import options
    clusterTable = readtable(csvFile, opts);
    
    
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
        [AcNmb, Sequences, numberOfClusters, clusterNames, pointsPerCluster, Sequences_withNs] = readFasta(dataSet);
        save(strcat(dataSet,'.mat'), 'AcNmb', 'Sequences', 'numberOfClusters', 'clusterNames', 'pointsPerCluster');
    end

    %% Map cluster IDs from the CSV to the sequences
    [isMatched, idx] = ismember(AcNmb, clusterTable.AC);

    % Filter out sequences that are not present in the clusterTable
    AcNmb = AcNmb(isMatched);
    Sequences = Sequences(isMatched);
    
    % Update ClassIDs based on filtered data
    ClassIDs = clusterTable.('label-val')(idx(isMatched));

    % Update the number of points per cluster (class)
    uniqueClasses = unique(ClassIDs);
    pointsPerCluster = cell(numel(uniqueClasses), 1);
    
    % Populate the cell array with the number of sequences per cluster
    for i = 1:numel(uniqueClasses)
        classID = uniqueClasses(i);
        % Store the count of sequences that belong to this cluster/class in the cell array
        pointsPerCluster{i} = sum(ClassIDs == classID);
    end

    % Calculate length stats
    [maxLen, minLen, meanLen, medLen] = lengthCalc(Sequences);

    %% Getting the starting point of each class as data is read in 1D cell array
    StartingPoint = [];
    counter = 1;
    for cls = 1:numberOfClusters
		StartingPoint = [StartingPoint, counter];
		counter = counter + pointsPerCluster{cls};
    end

    %% Determine the unique clusters
    ClusterIDs = clusterTable.cluster(idx(isMatched));
    uniqueClusters = unique(ClusterIDs);

    totalSeq = length(Sequences);


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

        % Getting the starting point for the class
        sp = StartingPoint(cls);

        % Starting point
        ppc = pointsPerCluster{cls};

        % Get the data for the current class
        classData = SequenceDataPooled{cls};

        % Get cluster IDs for the class data
        classClusterIDs = ClusterIDs(sp:sp+ppc-1);

        % Loop through each fold
		for testClusterIndex = 1:length(uniqueClusters)

            fold = testClusterIndex;

            testCluster = uniqueClusters(testClusterIndex);

            % train and test indices
            % Test indices are where the class's cluster ID matches the current testCluster
            testIndices = (classClusterIDs == testCluster);
            % Train indices are where the class's cluster ID does not match the testCluster
            trainIndices = ~testIndices;
            
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


metricsConsolidated = zeros(length(dataSets), 6*length(classifiers));
for ds=1:length(dataSets)
    metricsConsolidated(ds,:) = dataSetMetrics.(dataSets{ds});
end

metricsConsolidated = [dataSets', num2cell(round(metricsConsolidated,2))];
writecell(metricsConsolidated, "PP-ConsolidatedMetrics-Baseline-GP-0.9-4.xlsx");

disp("Processing finished!")
