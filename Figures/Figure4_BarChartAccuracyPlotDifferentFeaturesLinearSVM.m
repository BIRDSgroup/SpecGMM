% Specify the directory where your Excel files are stored
excelFilesDir = 'Seeds-16SrRNA/';

% List all Excel files in the directory
excelFiles = dir(fullfile(excelFilesDir, '*.xlsx'));

% Classifier name you're interested in
targetClassifier = 'LinearSVM'; % Use the exact sheet name as in your Excel files

% Initialize the accuracyCell array
accuracyCell = {};

% 'chanceAccuracy' is an array containing chance accuracy for each level
chanceAccuracy = [8.72, 33.4, 20.89, 28.86, 34.66, 11.83]; % Replace with your actual data

for fileIdx = 1:length(excelFiles)
    % Construct full file path
    filePath = fullfile(excelFilesDir, excelFiles(fileIdx).name);
    
    % Attempt to read the data from the target sheet
    try
        rawData = readcell(filePath, 'Sheet', targetClassifier);
    
        % Find the rows where the 'Fold' starts
        foldStartRows = find(cellfun(@(c) ischar(c) && startsWith(c, 'Fold'), rawData(:, 1)));
    
        % Append an extra element to handle the last fold
        foldStartRows(end+1) = size(rawData, 1) + 1;
    
        for foldIdx = 1:length(foldStartRows) - 1
            % Define the range for this fold's data
            startRow = foldStartRows(foldIdx) + 2; % Skip the fold title and header row
            endRow = startRow + 6 - 1;
        
            % Extract the accuracy data for the fold, assuming the first column is 'DataSet'
            foldData = rawData(startRow:endRow, 3:9); % Exclude the first column
        
            % Convert the cell array to a numeric matrix and store it in accuracyCell
            accuracyCell{end+1} = cell2mat(foldData);
        end
    catch ME
        warning('Failed to read classifier %s from file %s. Error: %s', targetClassifier, filePath, ME.message);
    end
end


% Load the accuracy cell
% load('accuracyCell_New.mat');

% Calculate the mean and standard deviation of accuracies across different seeds and folds
meanAccuracies = mean(cat(3, accuracyCell{:}), 3);  % Mean across the third dimension (folds)
stdAccuracies = std(cat(3, accuracyCell{:}), 0, 3);  % Standard deviation across the third dimension (folds)

% Define the taxonomic ranks for labeling
taxonomicRanks = {'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'};

% Define the features for labeling
features = {'All Windows', 'V2 Windows', 'V3 Windows', 'V4 Windows', 'V5 Windows', 'V6 Windows', 'V7 Windows'};

% Define a pastel color palette
colors = [0.7010, 0.7804, 0.8902;  % Pastel blue
          0.9882, 0.6902, 0.7529;  % Pastel red
          0.6431, 0.8588, 0.6667;  % Pastel green
          1.0000, 0.9333, 0.6353;  % Pastel yellow
          0.8392, 0.7059, 0.9098;  % Pastel purple
          0.9686, 0.7137, 0.8235;  % Pastel pink
          0.7333, 0.8706, 0.9843]; % Pastel light blue

% % Original colors
% originalColors = [
%     [203, 187, 136] / 255;
%     [201, 167, 123] / 255;
%     [168, 140, 108] / 255;
%     [101, 138, 110] / 255;
%     [145, 186, 141] / 255;
%     [195, 209, 140] / 255;
%     [247, 236, 168] / 255;
% ];
% 
% % Create pastel colors by averaging with white
% colors = (originalColors + 0.7) / 2;


% Create the figure
figure('Units', 'inches', 'Position', [0, 0, 12, 6]);

% Number of groups and number of bars in each group
numGroups = size(meanAccuracies, 1);
numBars = size(meanAccuracies, 2); 

% Set up the bar graph
barHandles = bar(meanAccuracies, 'grouped');
set(barHandles, {'FaceColor'}, num2cell(colors, 2));  % Assign colors to each bar group

% Add error bars
hold on;
for i = 1:numBars
    % Calculate the position for each set of feature error bars
    xPos = barHandles(i).XEndPoints;
    
    % Plot error bars for each feature
    errorbar(xPos, meanAccuracies(:, i), stdAccuracies(:, i), 'k', 'linestyle', 'none', 'HandleVisibility','off');
end

% Set the x-axis labels and their position
set(gca, 'xticklabel', taxonomicRanks);

% Improve the axes and add labels and title
xlabel('Taxonomic Rank', 'FontSize', 14);
ylabel('Accuracy', 'FontSize', 14);
title('Classification Accuracy by Feature and Taxonomic Rank', 'FontSize', 16);
set(gca, 'FontSize', 12);

% Add a legend
legend(features, 'FontSize', 10, 'Location', 'bestoutside');

% Adjust the axis limits if necessary
ylim([0, max(meanAccuracies(:) + stdAccuracies(:)) + 5]); % Adjust the limits based on your data range


% Loop through each group to place a horizontal line for chance accuracy
for i = 1:numGroups
    % Calculate the x position for the current group
    xPos = barHandles(i).XEndPoints;
    
    % Draw the horizontal line for chance accuracy
    line([i-0.4, i+0.4], [chanceAccuracy(i), chanceAccuracy(i)], 'Color', [200 100 16]/256, 'LineStyle', '--', 'LineWidth', 2);
end
hold off;