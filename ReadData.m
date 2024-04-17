%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AcNmb, Sequences, numberOfClusters, clusterNames, pointsPerCluster] = ReadData(dataSet)
	if isfile(strcat(dataSet,'.mat'))
        disp('Retreiving the preloaded dataset...')
        load(strcat(dataSet,'.mat'));
    else
        % Reading data from fasta files and saving the same
        [AcNmb, Sequences, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet);
        save(strcat(dataSet,'.mat'), 'AcNmb', 'Sequences', 'numberOfClusters', 'clusterNames', 'pointsPerCluster') 
    end
end