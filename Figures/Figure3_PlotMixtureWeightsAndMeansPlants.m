close all;
clear;
clc;

% Dataset
dataSet = 'Plants';

% Getting the current directory
currentDirectory = pwd;

% Arguments
args.magSpectra = 0;
args.trainUBM = 0;
args.visualizeTrain = 1;

% Parameters
parameters.distType = 'cor';               % Distance type, Euclidean -- 'euc', PCC -- 'cor'
parameters.winLen = 351;                   	% Window length
parameters.fftOrder= 2^nextpow2(parameters.winLen);   % FFT Order
parameters.winShift = 99;                   % Window shift
parameters.K = 5;                          % Number of mixture components in UBM-GMM
parameters.folds = 4;                      % Number of parameters.folds
parameters.magSpectrumStart = 1;           % 1--DC
parameters.magSpectrumEnd = parameters.fftOrder/2;
parameters.magSpectrumDimension = parameters.magSpectrumEnd-parameters.magSpectrumStart+1;

% Define the number of subplots
K = parameters.K;
N = parameters.magSpectrumDimension;

set(0, 'DefaultFigureRenderer', 'painters');

% Random seed
rng('default');

cd(currentDirectory)

% Reading Data
[data.AcNmb, data.Sequences, data.numberOfClusters, data.clusterNames, data.pointsPerCluster] = ReadData(dataSet);

% Getting the starting point of each class
StartingPoint = [];
counter = 1;
for cls = 1:data.numberOfClusters
    StartingPoint = [StartingPoint, counter];
    counter = counter + data.pointsPerCluster{cls};
end
data.StartingPoint = StartingPoint;

% Compute magnitude spectra
if args.magSpectra == 1
    fprintf('\nComputing Magnitude Spectra .... \n');
    SequenceDataPooled = ComputeMagnitudeSpectra(parameters, data);

    % Saving the extracted magnitude spectra
    save(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-MagSpecPooled.mat'), 'SequenceDataPooled', '-v7.3');
else
    fprintf('\nLoading Pre-computed Magnitude Spectra .... \n');
    load(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-MagSpecPooled.mat'), 'SequenceDataPooled');
end

% Train-test split
for cls = 1:data.numberOfClusters

    % Get the data for the current class
    classData = SequenceDataPooled{cls};

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
    end
end

% Mention the fold
fold = 1;

% Random seed
seed2 = 15;
rng(seed2,'twister');

fprintf("\n----- Fold: %d -----\n", fold);
disp(datetime)

% Height of matrix to store all magnitude spectra
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
    
    save(strcat(dataSet, '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(fold), '-', string(seed2), '-UBMGMM.mat'), 'UBMGMM', '-v7.3');
else
    fprintf('\nLoading Pre-trained UBM-GMM .... \n');
    load(strcat(dataSet, '-', string(parameters.K), '-', string(parameters.winLen), '-', string(parameters.winShift), '-', string(fold), '-', string(seed2), '-UBMGMM.mat'), 'UBMGMM');
end

% MAP Adaptation and Average Posteriors for each HVR
fprintf('\nMAP Adaptation and Average Posteriors for HVRs .... \n');

% Train data
% Initialize variables
superVectorsMatrixTrain = zeros(sum([pointsPerClusterTrain{fold}{:}]), parameters.K*parameters.magSpectrumDimension);
trainLabels = zeros(sum([pointsPerClusterTrain{fold}{:}]), 1);
meanVectorsTrain = cell(sum([pointsPerClusterTrain{fold}{:}]), 1);
mixtureComponentWeightsTrain = zeros(sum([pointsPerClusterTrain{fold}{:}]), parameters.K);

counter = 1;
% Loop over each class
for i = 1:data.numberOfClusters
    % Loop over points per cluster
    for j = progress(1:pointsPerClusterTrain{fold}{i})
        % MAP adaptation for data
        [superVector, pointsPerMixtureComponent, mixtureWeights, meanVectors, posteriors] = AdaptGMM(UBMGMM, SequenceDataPooledTrain{fold}{i}{j});

        % Store supervectors and labels
        meanVectorsTrain{counter} = meanVectors;
        superVectorsMatrixTrain(counter, :) = superVector;
        mixtureComponentWeightsTrain(counter, :) = mixtureWeights;
        trainLabels(counter,1) = i;
        counter = counter + 1;
    end
end

% Test data
% Initialize variables
superVectorsMatrixTest = zeros(sum([pointsPerClusterTest{fold}{:}]), parameters.K*parameters.magSpectrumDimension);
testLabels = zeros(sum([pointsPerClusterTest{fold}{:}]), 1);
meanVectorsTest = cell(sum([pointsPerClusterTest{fold}{:}]), 1);
mixtureComponentWeightsTest = zeros(sum([pointsPerClusterTest{fold}{:}]), parameters.K);

counter = 1;
% Loop over each class
for i = 1:data.numberOfClusters
    % Loop over points per cluster
    for j = progress(1:pointsPerClusterTest{fold}{i})
        % MAP adaptation for data
        [superVector, pointsPerMixtureComponent, mixtureWeights, meanVectors, posteriors] = AdaptGMM(UBMGMM, SequenceDataPooledTest{fold}{i}{j});

        % Store supervectors and labels
        meanVectorsTest{counter} = meanVectors;
        superVectorsMatrixTest(counter, :) = superVector;
        mixtureComponentWeightsTest(counter, :) = mixtureWeights;
        testLabels(counter,1) = i;
        counter = counter + 1;
    end
end

%% Visualizing mixture component activations
if args.visualizeTrain == 1
    X = mixtureComponentWeightsTrain;
    X = diag(1./sum(X,2))*X;
else
    X = mixtureComponentWeightsTest;
    X = diag(1./sum(X,2))*X;
end

% Define the figure size, adjust this to your article's requirements
figure('Units', 'inches', 'Position', [0, 0, 18, 3]); % Increase width to accommodate 6 subplots

% Plot the mixture weights as the first subplot
subplot(1, 6, 1);
imagesc(X);
colorbar;
xlabel('Mixture Component');
ylabel('Sequence');
title('Mixture Weights');

% Add the lines separating clusters
if args.visualizeTrain == 1
    x = pointsPerClusterTrain{fold}{1};
    for i = 1:data.numberOfClusters-1
        yline(x+0.5,'LineWidth',2,'Color','r');
        x = x + pointsPerClusterTrain{fold}{i+1};
    end
else
    x = pointsPerClusterTest{fold}{1};
    for i = 1:data.numberOfClusters-1
        yline(x+0.5,'LineWidth',2,'Color','r');
        x = x + pointsPerClusterTest{fold}{i+1};
    end
end
xticks([1 2 3 4 5]);

set(gca, 'FontSize', 12); % Match the font size with other subplots

% Select the sequence from each species of plants
if args.visualizeTrain == 1
    % Fold-1: Chlorophyta: 1–33 ; Streptophyta: 34–131 ;
    % Fold-2: Chlorophyta: 1–33 ; Streptophyta: 34–130 ;
    % Fold-3: Chlorophyta: 1–33 ; Streptophyta: 34–130 ;
    % Fold-4: Chlorophyta: 1–33 ; Streptophyta: 34–131 ;
    % 1 and 35 -- used for manuscript
    Species1Means = meanVectorsTrain{1};
    Species2Means = meanVectorsTrain{35};
else
    % Fold-1: Chlorophyta: 1–11 ; Streptophyta: 12–43 ;
    % Fold-2: Chlorophyta: 1–11 ; Streptophyta: 12–44 ;
    % Fold-3: Chlorophyta: 1–11 ; Streptophyta: 12–44 ;
    % Fold-4: Chlorophyta: 1–11 ; Streptophyta: 12–43 ;
    Species1Means = meanVectorsTest{6};
    Species2Means = meanVectorsTest{13};
end

% Loop through each of the remaining subplots to plot the data
for i = 1:K
    subplot(1, 6, i+1); % Adjust the index for subplots to start from 2
    f = linspace(0, pi, N); 
    plot(f, Species1Means(i, :), 'LineWidth', 1.5, 'MarkerFaceColor', '#2c85c5'); 
    hold on;
    plot(f, Species2Means(i, :), 'LineWidth', 1.5, 'MarkerFaceColor', '#db6839');
    hold off;
    
    % Set the axis limits
    xlim([0 pi]);
    ylim([10 30]);
    
    % Improve readability by setting ticks and labels only on the edges
    if i == 1
        ylabel('Magnitude', 'FontSize', 14);
    else
        set(gca, 'YTickLabel', []);
    end
    if i == K
        xlabel('Frequency (rad/sample)', 'FontSize', 14);
    else
        set(gca, 'XTickLabel', []);
    end
    set(gca, 'XTick', [0 pi], 'XTickLabel', {'0', '\pi'}, 'FontSize', 12);
    
    % Add a title to each subplot
    title(sprintf('Component %d', i), 'FontSize', 14);
end

% Add a single legend for the entire figure
legend('Chlorophyta', 'Streptophyta', 'FontSize', 14, 'Position', [0.9, 0.9, 0.1, 0.1]);

% Make sure the figure is displayed correctly
set(gcf, 'PaperPositionMode', 'auto');