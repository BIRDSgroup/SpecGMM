%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Compute Magnitude Spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SequenceDataPooled] = ComputeMagnitudeSpectra(parameters, data, numRepresentation)
	% Looping over classes
    for cls = 1:data.numberOfClusters

        fprintf('\n%s .... \n', data.clusterNames{cls});

        % Getting the starting point for the class
        sp = data.StartingPoint(cls);

        % Empty cell for storing magnitude spectra for class data
        SequenceData = cell(1,data.pointsPerCluster{cls});
        
        % Looping over class sequences
        for a = progress(1:data.pointsPerCluster{cls})

            % Mapping sequence to numerical representation
            ns = numMapping(data.Sequences{sp + a - 1}, numRepresentation);

            % Number of windows
            nWin = floor((length(ns)-parameters.winLen)/parameters.winShift);

            % Sequence length
            sequence_length = length(ns);

            % To store the magnitude spectra for a sequence
            seqData = zeros(nWin, parameters.magSpectrumDimension);

            % Looping over windows in a sequence
            for w = 1:nWin

                window_start = (w-1)*parameters.winShift + 1;                                   % Window start
                window_end = (w-1)*parameters.winShift + parameters.winLen;                     % Window end

                subSeq = ns(window_start: window_end);                                          % Getting subsequence

                sft = fft(subSeq, parameters.fftOrder);                                         % Computing FFT
                sft = sft(parameters.magSpectrumStart:parameters.magSpectrumEnd);               % Considering 0 to pi
                seqData(w,:) = abs(sft);                                                        % Computing magnitude spectrum
            end

            % Storing magnitude spectra for each sequence windows in a cell
            SequenceData{1, a} = seqData;
        end

        % Magnitude spectra for each class sequencces
        SequenceDataPooled{cls} = SequenceData;
    end
end
