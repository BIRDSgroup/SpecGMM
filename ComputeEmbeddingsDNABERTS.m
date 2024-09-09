%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Compute Embeddings using DNA-BERT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SequenceDataPooled] = ComputeEmbeddingsDNABERTS(data, dataset_name)

    % Define paths for the input and output files
    inputFilePath = 'sequences.txt';  % File with sequences
    outputFilePath = strcat('embeddings-DNABERTS-', dataset_name, '.npy'); % File to store embeddings

    % Write the sequences to a text file for the Python script
    fid = fopen(inputFilePath, 'w');
    for i = 1:length(data.Sequences)
        fprintf(fid, '%s\n', data.Sequences{i});
    end
    fclose(fid);

    % Call the Python script to extract embeddings
    disp("Calling Python function to extract embeddings...");
    pyenv('Version', 'path-to-venv/bin/python'); % Set up Python environment
    %pyenv('Version', '/tts3/saish/Research/SpecGMM/DNABERT-S/bin/python'); % Set up Python environment
    system(sprintf('python extract_embeddings_DNABERT_S.py %s %s', inputFilePath, outputFilePath));

    % Load the embeddings from the output file
    disp("Reading embeddings from pickle file...");
    embeddings = readNPY(outputFilePath); % Make sure to have npy-matlab package folder

    % Split the embeddings into clusters based on the data.StartingPoint
    SequenceDataPooled = cell(1, data.numberOfClusters);
    for cls = 1:data.numberOfClusters
        startIdx = data.StartingPoint(cls);
        endIdx = startIdx + data.pointsPerCluster{cls} - 1;
        SequenceDataPooled{cls} = embeddings(startIdx:endIdx, :);
    end
end