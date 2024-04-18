close all;
clear;
clc;

% Dataset
dataSet = 'Primates';

% Getting the current directory
currentDirectory = pwd;

% Arguments
args.magSpectra = 0;
args.trainUBM = 0;

% Parameters
parameters.distType = 'cor';               % Distance type, Euclidean -- 'euc', PCC -- 'cor'
parameters.winLen = 63;                   	% Window length
parameters.fftOrder= 2^nextpow2(parameters.winLen);   % FFT Order
parameters.winShift = 9;                   % Window shift
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

% % Plotting average spectrum for all the seuquences of given categories
% figure
% hold on
% frequencies = linspace(0, pi, parameters.magSpectrumDimension);
% overallAvgSpec = zeros(1,parameters.magSpectrumDimension);
% for cls =1:data.numberOfClusters
%     avgSpec = zeros(1,parameters.magSpectrumDimension);
%     for a = 1:data.pointsPerCluster{cls}
%         avgSpec = avgSpec + mean(SequenceDataPooled{1,cls}{1,a});
%     end
% 
%     avgSpec = avgSpec./data.pointsPerCluster{cls};
%     overallAvgSpec = overallAvgSpec + avgSpec;
%     %avgSpec = avgSpec./max(avgSpec);
%     plot(frequencies, avgSpec, 'LineWidth',2);
% end
% xlim([0,pi]);

% Plotting average spectrum for one sequence per category
figure
hold on
frequencies = linspace(0, pi, parameters.magSpectrumDimension);
for cls =1:data.numberOfClusters
    avgSpec = zeros(1,parameters.magSpectrumDimension);
    for a = 1:1
        avgSpec = avgSpec + mean(SequenceDataPooled{1,cls}{1,a});
    end
    plot(frequencies, avgSpec, 'LineWidth',2);
end
xlim([0,pi]);
xticks([0 1.05 2.1 pi]);
xticklabels({'0', '\pi/3', '2\pi/3', '\pi'});