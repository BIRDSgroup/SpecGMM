%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Main Wrappper Script -- SpecGMM Method with Homology Reduction
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

dataSets = {
'Primates'
};

%% Arguments
% If the dataset is one of the above mentioned 16S rRNA datasets and you want to do HVR analysis, set the following variable to 1
args.analysis16SrRNA = 0;

if args.analysis16SrRNA == 1
    % Load Taxonomy Table which contains the taxonomy and the HVR start-end information obtained using QIIME2
    load('DataBase/taxonomyTable.mat');
end

args.magSpectra = 1;        % Window-based magnitude spectra extraction: 1--extract, 0--load pre-computed spectra
args.assignHVR = 0;         % Assign HVR to each window based on maximum overlap: 1--assign, 0--load pre-assigned HVR
args.trainUBM = 1;          % Build UBM-GMM on the training data: 1--train, 0--load pre-built UBM-GMM
args.seedUBM = 15;          % Seed for training UBM-GMM
args.seedClassifier = 15;   % Seed for classifiers
args.numRepresentations = {"PP"};

%% Parameters
parameters.winLen = 351;                                 % Window length
parameters.fftOrder= 2^nextpow2(parameters.winLen);     % FFT Order
parameters.winShift = 99;                                % Window shift
parameters.K = 5;                          % Number of mixture components in UBM-GMM
parameters.folds = 4;                      % Number of parameters.folds / partitions
parameters.magSpectrumStart = 1;           % 1--DC
parameters.magSpectrumEnd = parameters.fftOrder/2;
parameters.magSpectrumDimension = parameters.magSpectrumEnd-parameters.magSpectrumStart+1;
parameters.threshold=0.9;

%% Define the classifiers to be used ('LinearDiscriminant', 'LinearSVM')
classifiers = {
    'LinearDiscriminant';
    'LinearSVM';
};

%% Add folders and subfolders in the current directory to the MATLAB path
% addpath(genpath('./'));

for numRep = 1:length(args.numRepresentations)
	
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

	%% Loop over datasets
	for ds = 1:length(dataSets)

		    % Display the dataset name
		    disp(dataSets{ds});
		    
		    % Display the numerical representation
		    disp(args.numRepresentations{numRep});

		    close all;
		    clearvars -except parameters classifiers dataSets ds currentDirectory taxonomyTable dataSetAccuracies dataSetMetrics args featureNames acVals svc seedVals numRep;
		    set(0, 'DefaultFigureRenderer', 'painters');

		    % Random seed
		    rng('default');

		    % Change the directory to SpecGMM/
		    cd(currentDirectory)

		    % Dataset
		    dataSet = dataSets{ds};

            % Define the path to your CSV file
			csvFile = strcat("graphpart_", dataSet, "-Merged-Labeled_part", string(parameters.folds), "_thresh", string(parameters.threshold),".csv");

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
		    [data.AcNmb, data.Sequences, data.numberOfClusters, data.clusterNames, data.pointsPerCluster, data.Sequences_withNs] = ReadData(dataSet);

            %% Map cluster IDs from the CSV to the sequences
            [isMatched, idx] = ismember(data.AcNmb, clusterTable.AC);

            % Filter out sequences that are not present in the clusterTable
            data.AcNmb = data.AcNmb(isMatched);
            data.Sequences = data.Sequences(isMatched);
            if ~isempty(data.Sequences_withNs)
                data.Sequences_withNs = data.Sequences_withNs(isMatched);
            end
            
            % Update ClassIDs based on filtered data
            data.ClassIDs = clusterTable.('label-val')(idx(isMatched));

            % Update the number of points per cluster (class)
            uniqueClasses = unique(data.ClassIDs);
            data.pointsPerCluster = cell(numel(uniqueClasses), 1);
            
            % Populate the cell array with the number of sequences per cluster
            for i = 1:numel(uniqueClasses)
                classID = uniqueClasses(i);
                % Store the count of sequences that belong to this cluster/class in the cell array
                data.pointsPerCluster{i} = sum(data.ClassIDs == classID);
            end


		    %% Getting the starting point of each class as data is read in 1D cell array
		    StartingPoint = [];
		    counter = 1;
		    for cls = 1:data.numberOfClusters
				StartingPoint = [StartingPoint, counter];
				counter = counter + data.pointsPerCluster{cls};
		    end
		    data.StartingPoint = StartingPoint;

            %% Determine the unique clusters
            data.ClusterIDs = clusterTable.cluster(idx(isMatched));
            uniqueClusters = unique(data.ClusterIDs);

		    %% Window-based magnitude spectra
		    if args.magSpectra == 1
				fprintf('\nComputing Magnitude Spectra .... \n');
				% Extract magnitude spectra
				SequenceDataPooled = ComputeMagnitudeSpectra(parameters, data, args.numRepresentations{numRep});

				% Saving the extracted magnitude spectra
				save(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-', args.numRepresentations{numRep}, '-', string(parameters.folds), '-', string(parameters.threshold), '-MagSpecPooled-GP.mat'), 'SequenceDataPooled', '-v7.3');
			else
				fprintf('\nLoading Pre-computed Magnitude Spectra .... \n');
				% Loading pre-computed magnitude spectra
				load(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-', args.numRepresentations{numRep}, '-', string(parameters.folds), '-', string(parameters.threshold), '-MagSpecPooled-GP.mat'));
            end

			%% Train-test split for K-fold cross validation
			counter=1;
			for cls = 1:data.numberOfClusters

                % Getting the starting point for the class
                sp = data.StartingPoint(cls);

				% Starting point
				ppc = data.pointsPerCluster{cls};

				% Get the magnitude spectra data for the current class
				classData = SequenceDataPooled{cls};

				% Number of sequences in this class
				numDataPoints = numel(classData);

                classClusterIDs = data.ClusterIDs(sp:sp+ppc-1);

				% Loop through each fold
				for testClusterIndex = 1:length(uniqueClusters)

                    fold = testClusterIndex;

                    testCluster = uniqueClusters(testClusterIndex);

                    % train and test indices
                    % Test indices are where the class's cluster ID matches the current testCluster
                    testIndices = (classClusterIDs == testCluster);
                    % Train indices are where the class's cluster ID does not match the testCluster
                    trainIndices = ~testIndices;

					% Points per mixture component
					pointsPerClusterTrain{fold}{cls} = sum(trainIndices);
					pointsPerClusterTest{fold}{cls} = sum(testIndices);

                    % Extract the magnitude spectra data for the current fold
					SequenceDataPooledTrain{fold}{cls} = classData(trainIndices);
					SequenceDataPooledTest{fold}{cls} = classData(testIndices);

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

		    correctLengths = cell(1,length(classifiers));
		    incorrectLengths = cell(1,length(classifiers));

		    correctCountNs = cell(1,length(classifiers));
		    incorrectCountNs = cell(1,length(classifiers));

		    correctFractionNs = cell(1,length(classifiers));
		    incorrectFractionNs = cell(1,length(classifiers));

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
					
					save(strcat(dataSet, '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(fold), '-', string(args.seedUBM), '-', args.numRepresentations{numRep}, '-UBMGMM-GP.mat'), 'UBMGMM', '-v7.3');
				else
					fprintf('\nLoading Pre-trained UBM-GMM .... \n');
					load(strcat(dataSet, '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(fold), '-', string(args.seedUBM), '-', args.numRepresentations{numRep}, '-UBMGMM-GP.mat'), 'UBMGMM');
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
					[accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, predLabels, conf_Matrix] = classificationCode(SuperVectorsMatrixTrain, trainLabels, SuperVectorsMatrixTest, testLabels, classifiers, args.seedClassifier);
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
						[accuracy, clNames, sensitivity, specificity, precision, recall, f1_score, predLabels, conf_Matrix] = classificationCode(hvrActivationsTrain{featureIndex}, trainLabels, hvrActivationsTest{featureIndex}, testLabels, classifiers, args.seedClassifier);
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

	% Consolidate all the metrics -- (5 metrics -- accuracy, precision, recall, specificity, f1_score) + std. dev
	metricsConsolidated = zeros(length(dataSets), 6*length(classifiers));
	for ds=1:length(dataSets)
	    metricsConsolidated(ds,:) = dataSetMetrics.SV.(dataSets{ds});
	end
	metricsConsolidated = [dataSets', num2cell(round(metricsConsolidated,2))];
	writecell(metricsConsolidated, strcat(args.numRepresentations{numRep}, "-", "ConsolidatedMetrics-GP.xlsx"));

end % end of numRepresentation loop

disp("Processing finished!")