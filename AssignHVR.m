%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Assign HVR to each window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Assigned_HVRs_DataPooled] = AssignHVR(parameters, data, taxonomyTable)

	% Loop over the classes
	for cls = 1:data.numberOfClusters

        % Getting the starting point for the class
        sp = data.StartingPoint(cls);

        % Empty cell for storing assigned HVR each sequence of the class
        Assigned_HVRs = cell(1,data.pointsPerCluster{cls});

        % Looping over class sequences
        for a = progress(1:data.pointsPerCluster{cls})

            % Get the sequence header name
            header_name = data.AcNmb{sp + a - 1};

            % Find the row in taxonomyTable that matches the header name
            match_row = find(strcmp(taxonomyTable.Header, header_name));

            % Sequence length
            sequence_length = length(data.Sequences{sp + a - 1});

            % Number of windows
            nWin = floor((sequence_length-parameters.winLen)/parameters.winShift);

            % To store HVR assignment for the sequence
            assigned_HVRs_seq = zeros(1,nWin);

            % Looping over windows in a sequence
            for w = 1:nWin

                window_start = (w-1)*parameters.winShift + 1;      					% Window start
                window_end = (w-1)*parameters.winShift + parameters.winLen;   	% Window end

                max_overlap = 0;                        						% Maximum overlap
                max_overlap_HVR = 0;                    						% Maximum overlap HVR

                % Iterate through all HVRs (V2 to V7)
                for hvrIndex = 2:7

                    % Start and end positions of HVR header in CSV
                    hvr_start_col = sprintf('V%dStart', hvrIndex);
                    hvr_end_col = sprintf('V%dEnd', hvrIndex);

                    % Extract the HVR start and end positions for the matching row
                    hvr_start = taxonomyTable{match_row, hvr_start_col};
                    hvr_end = taxonomyTable{match_row, hvr_end_col};

                    % Calculate overlap with the extracted HVR
                    overlap = max(0, min(window_end, hvr_end) - max(window_start, hvr_start));

                    % If the overlap is more than the previous one, update
                    if overlap > max_overlap
                        max_overlap = overlap;
                        max_overlap_HVR = hvrIndex;
                    end

                end

                % Assign the window to the HVR with maximum overlap
                assigned_HVRs_seq(1, w) = max_overlap_HVR;
            end

            % Storing assigned HVR for each sequence windows in a cell
            Assigned_HVRs{1,a} = assigned_HVRs_seq;
        end

        % Assigned HVRs for each class
        Assigned_HVRs_DataPooled{cls} = Assigned_HVRs;
    end
    
end
