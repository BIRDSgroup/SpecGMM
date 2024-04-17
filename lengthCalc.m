%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: To calculate descriptive statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ maxLen, minLen, meanLen, medLen, modLen, stdDev, meanAbsDev, medAbsDev ] = lengthCalc( Seq )
    %function calculates length stats
    len = cellfun('length',Seq);
    maxLen = max(len);
    minLen = min(len);
    meanLen = round(mean(len));
    medLen = round(median(len));
    modLen = mode(len);
    stdDev = round(std(len));
    meanAbsDev = round(mad(len, 0));
    medAbsDev = round(mad(len, 1));
end
