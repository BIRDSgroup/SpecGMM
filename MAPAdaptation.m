%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: MAP Adaptation — Sequence-specific GMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [superVectorsMatrix, labels, mixtureComponentWeights] = MAPAdaptation(SequenceDataPooled, pointsPerCluster, UBMGMM, data, parameters, fold)
    % Initialize variables
    superVectorsMatrix = zeros(sum([pointsPerCluster{fold}{:}]), parameters.K*parameters.magSpectrumDimension);
    labels = zeros(sum([pointsPerCluster{fold}{:}]), 1);
    counter = 1;
    mixtureComponentWeights = zeros(sum([pointsPerCluster{fold}{:}]), parameters.K);

    % Loop over each class
    for i = 1:data.numberOfClusters
        % Loop over points per cluster
        for j = progress(1:pointsPerCluster{fold}{i})
            % MAP adaptation for data
            [superVector, pointsPerMixtureComponent, mixtureWeights, meanVectors, posteriors] = AdaptGMM(UBMGMM, SequenceDataPooled{fold}{i}{j});

            % Store supervectors and labels
            superVectorsMatrix(counter, :) = superVector;
            mixtureComponentWeights(counter, :) = mixtureWeights;
            labels(counter,1) = i;
            counter = counter + 1;
        end
    end
end
