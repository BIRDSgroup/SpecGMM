%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AcNmb, Sequences, numberOfClusters, clusterNames, pointsPerCluster, Sequences_withNs] = ReadData(dataSet)
    if isfile(strcat(dataSet, '.mat'))
        disp('Retrieving the preloaded dataset...');
        load(strcat(dataSet, '.mat'));
        
        % Check if Sequences_withNs exists in the loaded file
        if exist('Sequences_withNs', 'var') == 0
            disp('Sequences_withNs not found in the preloaded dataset, assigning empty value.');
            Sequences_withNs = {}; % Assign an empty cell array or any default value
        end
        
    else
        % Reading data from fasta files and saving the same
        disp('Reading sequences ....');
        [AcNmb, Sequences, numberOfClusters, clusterNames, pointsPerCluster, Sequences_withNs] = readFasta(dataSet);
        save(strcat(dataSet, '.mat'), 'AcNmb', 'Sequences', 'numberOfClusters', 'clusterNames', 'pointsPerCluster', 'Sequences_withNs');
    end
end
