%% If feature extraction, HVR assignment, and UBM training needs to be done, then set this variable to 1
args.magSpectra = 0;

%% Parameters
parameters.distType = 'cor';               % Distance type, Euclidean -- 'euc', PCC -- 'cor'
parameters.winLen = 63;                   	% Window length
parameters.fftOrder= 2^nextpow2(parameters.winLen);   % FFT Order
parameters.winShift = 9;                   % Window shift
parameters.K = 5;                          % Number of mixture components in UBM-GMM
parameters.folds = 4;                      % Number of parameters.folds
parameters.magSpectrumStart = 1;           % 1--DC
parameters.magSpectrumEnd = parameters.fftOrder/2;
parameters.magSpectrumDimension = parameters.magSpectrumEnd-parameters.magSpectrumStart+1;

%% Dataset
dataSet = 'New_HVR_16SrRNA_GenusToSpecies_Bacillus';

%% Reading Data
[data.AcNmb, data.Sequences, data.numberOfClusters, data.clusterNames, data.pointsPerCluster] = ReadData(dataSet);

%% Getting the starting point of each class
StartingPoint = [];
counter = 1;
for cls = 1:data.numberOfClusters
    StartingPoint = [StartingPoint, counter];
    counter = counter + data.pointsPerCluster{cls};
end
data.StartingPoint = StartingPoint;

%% Compute magnitude spectra
if args.magSpectra == 1
    fprintf('\nComputing Magnitude Spectra .... \n');
    SequenceDataPooled = ComputeMagnitudeSpectra(parameters, data);

    % Saving the extracted magnitude spectra
    save(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-MagSpecPooled.mat'), 'SequenceDataPooled', '-v7.3');
else
    fprintf('\nLoading Pre-computed Magnitude Spectra .... \n');
    load(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-MagSpecPooled.mat'), 'SequenceDataPooled');
end

% %% Load Assigned HVR
load(strcat(dataSet, '-', string(parameters.winLen), '-', string(parameters.winShift), '-assignedHVRs_DataPooled.mat'), 'assignedHVRs_DataPooled');


%{
%% Create spectrograms for each class' numSeq sequences
numSeq = 200;
for cls = 1:data.numberOfClusters
    sp = data.StartingPoint(cls);
    for a = progress(1:numSeq)
        specGram = SequenceDataPooled{1,cls}{1,a};
        h = figure('Position', [100, 100, 1800, 300], 'visible','off');
        imagesc(specGram');

        % Set the y-axis limits from 0 to pi and label it
        yticks([1, parameters.magSpectrumDimension]); % Assuming parameters.magSpectrumDimension corresponds to pi
        yticklabels({'0', '\pi'}); % Label the y-ticks accordingly

        % Set the x-axis label
        xlabel('Window', 'FontSize', 14, 'FontWeight', 'bold');

        % Set the y-axis label
        ylabel('Frequency (rad/sample)', 'FontSize', 14, 'FontWeight', 'bold');

        % Generate the title with appropriate formatting
        titleStr = strcat(data.clusterNames{cls}, '-', data.AcNmb{sp + a - 1}, '-', string(cls), '-', string(a));
        titleStr = replace(titleStr, '_', '\_');
        title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'YDir', 'normal', 'FontSize', 14);

        % Create the directory if it does not exist
        plotDir = fullfile('Plots', dataSet, data.clusterNames{cls});
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end

        % Construct the file path
        filePath = fullfile(plotDir, strcat(data.AcNmb{sp + a - 1}, '-', string(cls), '-', string(a), '.png'));
        
        % Save the figure
        exportgraphics(h, filePath, 'Resolution', 300)
        close(h);
    end
end

%% Make a video for each class
% Set up parameters
for cls = 1:data.numberOfClusters
    outputVideoFile = strcat('Plots/', dataSet, '/', data.clusterNames{cls}, '.mp4');   % Choose an MP4 video file name
    imageFolder = strcat('Plots/', dataSet, '/', data.clusterNames{cls});               % Specify the path to the folder containing PNG images
    imageFiles = dir(fullfile(imageFolder, '*.png'));  % Use '*.png' for PNG images

    % Create VideoWriter object with MPEG-4 profile
    outputVideo = VideoWriter(outputVideoFile, 'MPEG-4');
    outputVideo.FrameRate = 16;  % Adjust the frame rate as needed
    open(outputVideo);

    % Read images, resize, and write to video
    for i = 1:length(imageFiles)
        % Read the PNG image
        img = imread(fullfile(imageFolder, imageFiles(i).name));

        % Resize the image to the expected dimensions
        resizedImg = imresize(img, [1212, 5948]);

        % Write the resized image to the video
        writeVideo(outputVideo, resizedImg);
    end

    % Close the VideoWriter object
    close(outputVideo);
end
%}

%% Getting magnitude spectra computed over sliding windows for a sequence in a class
% In manuscript Bacillus Subtilis: cls=7, a=29
cls=7;  % class
a=29;   % data point
sp = data.StartingPoint(cls);
disp(data.AcNmb{sp + a - 1})

specGram = SequenceDataPooled{1,cls}{1,a};  % spectrogram
h = figure('Position', [100, 100, 1800, 300]);
imagesc(specGram');
set(gca, 'YDir', 'normal', 'FontSize', 14);

% Set the y-axis limits from 0 to pi and label it
yticks([1, parameters.magSpectrumDimension]); % Create 2 y-ticks from 0 to pi
yticklabels({'0', '\pi'}); % Label the y-ticks accordingly

% Set the x-axis label
xlabel('Window', 'FontSize', 14, 'FontWeight', 'bold');

% Set the y-axis label
ylabel('Frequency (rad/sample)', 'FontSize', 14, 'FontWeight', 'bold');

% Set the title
titleStr = strcat(data.clusterNames{cls}, '-', data.AcNmb{sp + a - 1}, '-', string(cls), '-', string(a));
titleStr = replace(titleStr, '_', '\_');
title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');

%% Get sequence corresponding to some consecutive windows
sp = data.StartingPoint(cls);
seq = data.Sequences{sp + a - 1};

w1=13; % starting window of motif
w2=16; % ending window of motif

startIndex = (w1-1)*parameters.winShift + 1;
endIndex = (w2-1)*parameters.winShift + parameters.winLen;

subSeq = seq(startIndex:endIndex);
disp(subSeq)
%% Find the number of sequences in each class containing the subseq
motifPresence = zeros(1,data.numberOfClusters);
for cls = 1:data.numberOfClusters
    sp = data.StartingPoint(cls);
    count=0;
    for a = progress(1:data.pointsPerCluster{cls})
        if contains(data.Sequences{sp + a - 1}, subSeq)
            motifPresence(cls) = motifPresence(cls) + 1;
        end
    end
end

%% Print the motif presence
for cls = 1:data.numberOfClusters
    fprintf("%-20s:", data.clusterNames{cls});
    fprintf(" %d/%d", motifPresence(cls), data.pointsPerCluster{cls});
    fprintf("\n");
end