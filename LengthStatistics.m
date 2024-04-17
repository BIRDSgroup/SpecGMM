%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Compute Length Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stats] = LengthStatistics(data)
    counter = 1;
    for cls = 1:data.numberOfClusters
        % Starting point
        ppc = data.pointsPerCluster{cls};

        % Length statistics
        [maxLen, minLen, meanLen, medLen, modLen, stdDev, meanAbsDev, medAbsDev] = lengthCalc(data.Sequences(counter:counter+ppc-1));
        statistics(cls,:) = {maxLen, minLen, meanLen, medLen, modLen, stdDev, meanAbsDev, medAbsDev};

        counter = counter + ppc;
    end
    stats = [data.clusterNames', data.pointsPerCluster', statistics];
end