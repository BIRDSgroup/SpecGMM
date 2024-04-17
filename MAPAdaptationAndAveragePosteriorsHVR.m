%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: MAP Adaptation and Getting Average Posteriors for each HVR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SuperVectorsMatrix, Labels, mixtureActivationsPerHVRPerSequencePerCluster] = MAPAdaptationAndAveragePosteriorsHVR(SequenceDataPooled, assignedHVRs_DataPooled, pointsPerCluster, UBMGMM, data, parameters, fold)
    % Initialize variables
    SuperVectorsMatrix = zeros(sum([pointsPerCluster{fold}{:}]), parameters.K*parameters.magSpectrumDimension);;
    Labels = [];
    mixtureActivationsPerHVRPerSequencePerCluster = cell(1, data.numberOfClusters);
    counter = 1;

    % Loop over each class
    for i = 1:data.numberOfClusters
        % Initialize -- looping over HVRs (1 to be ignored)
        averagePosteriorsPerHVRPerSequence = cell(1, 8);
        for hvrIndex = 1:8
            averagePosteriorsPerHVRPerSequence{hvrIndex} = [];
        end

        % Loop over points per cluster
        for j = progress(1:pointsPerCluster{fold}{i})
            % MAP adaptation for data
            [SuperVector, pointsPerMixtureComponent, mixtureWeights, meanVectors, posteriors] = AdaptGMM(UBMGMM, SequenceDataPooled{fold}{i}{j});

            % Store supervectors and labels
            SuperVectorsMatrix(counter, :) = SuperVector;
            Labels = [Labels; i];
            counter = counter + 1;

            % Getting assigned HVR for each window
            assignedHVRs = assignedHVRs_DataPooled{fold}{i}{j};

            % Process assigned HVRs
            [averagePosteriorsPerHVRPerSequence] = processAssignedHVRs(assignedHVRs, mixtureWeights, posteriors, averagePosteriorsPerHVRPerSequence, parameters);
        end

        % Average posteriors for HVRs in the class
        mixtureActivationsPerHVRPerSequencePerCluster{i} = averagePosteriorsPerHVRPerSequence;
    end
end

function [averagePosteriorsPerHVRPerSequence] = processAssignedHVRs(assignedHVRs, mixtureWeights, posteriors, averagePosteriorsPerHVRPerSequence, parameters)
    % Initialize cell arrays to store average posteriors and counts for each HVR
    averagePosteriorsPerHVR = cell(1, 8);
    windowCountPerHVR = zeros(1, 8);

    % Initialize counters for each HVR
    for hvrIndex = 1:8
        averagePosteriorsPerHVR{hvrIndex} = zeros(1, parameters.K);
    end

    % Iterate through each window
    for windowIndex = 1:length(assignedHVRs)
        % Retrieve the assigned HVR for the current window
        assignedHVR = assignedHVRs(windowIndex);

        % Check if the assigned HVR is valid (greater than 0)
        if assignedHVR > 0
            % sum(gamma_nk for HVR)
            averagePosteriorsPerHVR{assignedHVR} = averagePosteriorsPerHVR{assignedHVR} + posteriors(windowIndex, :);
            windowCountPerHVR(assignedHVR) = windowCountPerHVR(assignedHVR) + 1;
        end
    end

    % Normalization -- average posteriors = sum(gamma_nk for HVR)/No. of win for HVR
    averagePosteriorsPerHVR{1} = mixtureWeights; % First one stores the mixture weights [all windows]
    for hvrIndex = 2:8
        if windowCountPerHVR(hvrIndex) > 0
            averagePosteriorsPerHVR{hvrIndex} = averagePosteriorsPerHVR{hvrIndex} / windowCountPerHVR(hvrIndex);
        end
    end

    % Store average posteriors per cluster per sequence
    for hvrIndex = 1:8
        averagePosteriorsPerHVRPerSequence{hvrIndex} = [averagePosteriorsPerHVRPerSequence{hvrIndex}; averagePosteriorsPerHVR{hvrIndex}];
    end
end

