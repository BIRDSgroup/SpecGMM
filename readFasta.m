%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Read Fasta File
% Adapted from Randhawa et al., 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster, Seq_withNs] = readFasta(dataSet)
    % Read fasta files in folder DataBase/'dataSet'
    % Make subFolders for each cluster and place respective fasta sequences inside

    % Base paths
    basePath = fullfile(pwd, 'DataBase');
    dataSetPath = fullfile(basePath, dataSet);
    matFilePath = fullfile(basePath, strcat(dataSet, '.mat'));

    % Initialize output variables
    Seq = {};
    Seq_withNs = {};
    AcNmb = {};
    clusterNames = {};
    pointsPerCluster = {};

    % Check if the .mat file already exists
    if isfile(matFilePath)
        disp('Retrieving the preloaded dataset...');
        load(matFilePath, 'AcNmb', 'Seq', 'numberOfClusters', 'clusterNames', 'pointsPerCluster', 'Seq_withNs');
        return;
    end

    % Get folder info
    folderInfo = dir(dataSetPath);
    if isempty(folderInfo)
        error('Error: DataSet does not exist.');
    end

    % Process each subfolder
    numFolders = length(folderInfo) - 2;
    for i = 3:numFolders + 2
        subFolderPath = fullfile(dataSetPath, folderInfo(i).name);
        if ~isfolder(subFolderPath)
            continue;
        end
        
        clusterNames{end+1} = folderInfo(i).name; % Add cluster name
        subFolderInfo = dir(fullfile(subFolderPath, '*.fasta')); % Only process fasta files
        pts = length(subFolderInfo);

        % Preallocate cell arrays
        seqTemp = cell(1, pts);
        seqTemp_withNs = cell(1, pts);
        acTemp = cell(1, pts);

        for j = 1:pts
            fastaFile = fullfile(subFolderPath, subFolderInfo(j).name);
            [Header, Sequence] = fastaread(fastaFile);

            % Check and store sequences containing Ns
            if contains(Sequence, 'N', 'IgnoreCase', true)
                seqTemp_withNs{j} = Sequence;
            end

            % Remove other characters
            Sequence = regexprep(Sequence, '[^A,^C,^G,^T]', '', 'ignorecase');

            % Extract only the part before the first '|'
            Header = strtok(Header, '|');

            % Store the processed header and sequence
            seqTemp{j} = Sequence;
            acTemp{j} = Header;
        end

        % Update output variables
        Seq = [Seq, seqTemp];
        Seq_withNs = [Seq_withNs, seqTemp_withNs(~cellfun('isempty', seqTemp_withNs))];
        AcNmb = [AcNmb, acTemp];
        pointsPerCluster{end+1} = pts;
    end

    numberOfClusters = length(clusterNames);

    % Save the results in a .mat file for future use
    save(matFilePath, 'AcNmb', 'Seq', 'numberOfClusters', 'clusterNames', 'pointsPerCluster', 'Seq_withNs');

end



% function [AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster, Seq_withNs] = readFasta(dataSet)
%     % read fasta files in folder DataBase/'dataSet'
%     % make subFolders for each cluster and place respective fasta sequences insides
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Author: Gurjit Singh Randhawa  %
%     % Department of Computer Science,%
%     % Western University, Canada     %
%     % email: grandha8@uwo.ca         %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         path = pwd;
%         dbPath = strcat(path,'/','DataBase');
%         cd(dbPath);
%         tStr = strcat(dataSet,'.mat');
%         if ~isempty(dir(tStr))
%             load(tStr);
%         else
%             dbPath = strcat(dbPath,'/',dataSet);
%             folderInfo = dir (dbPath);
%             Seq = {};
%             Seq_withNs = {};
%             AcNmb = {};    
    
%             if(isempty(folderInfo))
%                  error('Error : DataSet does not exist.')
%             else
%     %             numberOfClusters = length(folderInfo)-2;
%     %             clusterNames = cell(1,numberOfClusters);
%     %             pointsPerCluster = cell(1,numberOfClusters);
%                 cd(dbPath); 
%                 for i=3:length(folderInfo)
%                    fnm=folderInfo(i).name; 
%                    subFolderPath = strcat(dbPath,'/',fnm);
%                    subFolderInfo = dir(subFolderPath);
%                    result = subFolderInfo.isdir;
%                    if ~(result)
%                         flag = downloadFasta( fnm );
%                    end
%                 end
                
%                 numberOfClusters = 0;
%                 index = 1;
%                 folderInfo = dir (dbPath);
%                 for i=3:length(folderInfo)
%                     cd(dbPath);                
%                     subFolderPath = strcat(dbPath,'/',folderInfo(i).name);
%                     subFolderInfo = dir(subFolderPath);
%                     result = subFolderInfo.isdir;
%                     if ~(result)
%                         continue;
%                     end
%                     cd(subFolderPath);
%                     clusterNames{index} =  folderInfo(i).name;
%                     pts = length(subFolderInfo)-2;
%                     pointsPerCluster{index} = pts;
%                     seqTemp = cell(1,pts);
%                     seqTemp_withNs = cell(1,pts);
%                     acTemp = cell(1,pts);
%                     for j=progress(3:length(subFolderInfo))
%                         [Header, Sequence] = fastaread(subFolderInfo(j).name);
                        
%                         % If the sequence contains Ns
%                         if contains(Sequence, 'N', 'IgnoreCase', true)
%                             seqTemp_withNs{j-2} = Sequence;
%                         end

%                         % Remove other characters
%                         Sequence = regexprep(Sequence,'[^A,^C, ^G, ^T]','','ignorecase');
    
%                         % Extract only the part before the first '|'
%                         Header = strtok(Header, '|');
                        
%                         % Store the processed header and sequence
%                         seqTemp{j-2} = Sequence;
%                         acTemp{j-2} = Header; 
%                     end
%                     Seq = [Seq seqTemp];
%                     Seq_withNs = [Seq_withNs seqTemp_withNs];
%                     AcNmb = [AcNmb acTemp]; 
%                     numberOfClusters = numberOfClusters+1;
%                     index = index+1;
%                 end 
%                 if(numberOfClusters<2)
%                     error('Error : DataSet should have atleast two clusters.')            
%                 end
%             end
%         end
%         cd(path);
%     end
    