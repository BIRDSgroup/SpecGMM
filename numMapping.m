%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saish Jaiswal
% Indian Institute of Technology Madras
% Function: Numerical Mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for different nucleotide sequence mapping methods
% Input:
%   - seq: input nucleotide sequence
%   - mappingMethod: string specifying the mapping method
% Output:
%   - numSeq: numeric sequence based on the selected mapping method

function numSeq = numMapping(seq, mappingMethod)

    switch mappingMethod
        case 'AT_CG'
            numSeq = numMappingAT_CG(seq);
        case 'Atomic'
            numSeq = numMappingAtomic(seq);
        case 'Codons'
            numSeq = numMappingCodons(seq);
        case 'Doublet'
            numSeq = numMappingDoublet(seq);
        case 'EIIP'
            numSeq = numMappingEIIP(seq);
        case 'Int'
            numSeq = numMappingInt(seq);
        case 'IntN'
            numSeq = numMappingIntN(seq);
        case 'JustA'
            numSeq = numMappingJustA(seq);
        case 'JustC'
            numSeq = numMappingJustC(seq);
        case 'JustG'
            numSeq = numMappingJustG(seq);
        case 'JustT'
            numSeq = numMappingJustT(seq);
        case 'PP'
            numSeq = numMappingPP(seq);
        case 'Random3'
            numSeq = numMappingRandom3(seq);
        case 'Random13'
            numSeq = numMappingRandom13(seq);
        case 'Real'
            numSeq = numMappingReal(seq);
        otherwise
            error('Invalid mapping method. Please choose a valid method.');
    end
end

function numSeq = numMappingAT_CG(seq)
    % PairedNumeric representation
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if t == 'A'
            numSeq(K) = 1;
        elseif t == 'C'
            numSeq(K) = -1;
        elseif t == 'G'
            numSeq(K) = -1;
        elseif t == 'T'
            numSeq(K) = 1;
        end
    end
end

function numSeq = numMappingAtomic(seq)
    % Atomic representation
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if t == 'A'
            numSeq(K) = 70;
        elseif t == 'C'
            numSeq(K) = 58;
        elseif t == 'G'
            numSeq(K) = 78;
        elseif t == 'T'
            numSeq(K) = 66;
        end
    end
end

function numSeq = numMappingCodons(seq)
    % Codon representation
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    codons = {'TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', ...
              'TAA', 'TAG', 'TGA', 'TGT', 'TGC', 'TGG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', ...
              'CGA', 'CGG', 'AGA', 'AGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', ...
              'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'};

    for K = 1:len
        if K < len - 1
            t = seq(K:K + 2);
        elseif K == len - 1
            t = strcat(seq(K:K + 1), seq(1:1));
        else
            t = strcat(seq(K:K), seq(1:2));
        end
        tp = find(strcmpi(codons, t)) - 1;
        numSeq(K) = tp;
    end
end

function numSeq = numMappingDoublet(seq)
    % Nearest-neighbor based doublet representation
    len = length(seq);
    doublet = {'AA', 'AT', 'TA', 'AG', 'TT', 'TG', 'AC', 'TC', 'GA', 'CA', 'GT', 'GG', 'CT', 'GC', 'CG', 'CC'};
    alpha = 0;
    numSeq = zeros(1, len, 'double');
    kStrings = (2 * alpha) + 1;
    for K = 1:len
        if alpha == 0
            if K < len
                t = seq(K:K + 1);
            else
                t = strcat(seq(K:K), seq(1:1));
            end
            tp = find(strcmpi(doublet, t)) - 1;
            numSeq(K) = tp;
        else
            loc = 0;
            for index = K - alpha:K + alpha
                sPos = index;
                ePos = sPos + 1;
                if sPos < 1
                    sPos = len + sPos;
                elseif sPos > len
                    sPos = sPos - len;
                end
                if ePos > len
                    ePos = ePos - len;
                elseif ePos < 1
                    ePos = ePos + len;
                end
                if sPos == len && ePos == 1
                    t = strcat(seq(len:len), seq(1:1));
                else
                    t = seq(sPos:ePos);
                end
                loc = loc + (find(strcmpi(doublet, t)) - 1);
            end
            tp = loc / kStrings;
            numSeq(K) = tp;
        end
    end
end

function numSeq = numMappingEIIP(seq)
    % EIIP representation
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if t == 'A'
            numSeq(K) = 0.1260;
        elseif t == 'C'
            numSeq(K) = 0.1340;
        elseif t == 'G'
            numSeq(K) = 0.0806;
        elseif t == 'T'
            numSeq(K) = 0.1335;
        end
    end
end

function numSeq = numMappingInt(seq)
    % Integer representation
    dob = {'T', 'C', 'A', 'G'};
    len = length(seq);
    numSeq = zeros(1, len, 'double');

    for K = 1:len
        t = seq(K);
        tp = find(strcmpi(dob, t)) - 1;
        numSeq(K) = tp;
    end
end

function numSeq = numMappingIntN(seq)
    % Integer (other variant) representation
    dob = {'T', 'C', 'A', 'G'};
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        tp = find(strcmpi(dob, t));
        numSeq(K) = tp;
    end
end

function numSeq = numMappingJustA(seq)
    % JustA representation
    a = 'A';
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if strcmpi(t, a)
            numSeq(K) = 1;
        else
            numSeq(K) = 0;
        end
    end
end

function numSeq = numMappingJustC(seq)
    % JustC representation
    a = 'C';
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if strcmpi(t, a)
            numSeq(K) = 1;
        else
            numSeq(K) = 0;
        end
    end
end

function numSeq = numMappingJustG(seq)
    % JustG representation
    a = 'G';
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if strcmpi(t, a)
            numSeq(K) = 1;
        else
            numSeq(K) = 0;
        end
    end
end

function numSeq = numMappingJustT(seq)
    % JustT representation
    a = 'T';
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if strcmpi(t, a)
            numSeq(K) = 1;
        else
            numSeq(K) = 0;
        end
    end
end

function numSeq = numMappingPP(seq)
    % Purine/Pyramidine representation
    len = length(seq);
    numSeq = zeros(1, len, 'double');
    for K = 1:len
        t = seq(K);
        if strcmpi(t, 'A')
            numSeq(K) = -1;
        elseif strcmpi(t, 'C')
            numSeq(K) = 1;
        elseif strcmpi(t, 'G')
            numSeq(K) = -1;
        elseif strcmpi(t, 'T')
            numSeq(K) = 1;
        end
    end
end

function result = numMappingRandom3(seq)
    % Random among 3 representations
    a = 3;
    r = randi([1 a], 1, 1);
    if r == 1
        result = numMappingPP(seq);
    elseif r == 2
        result = numMappingReal(seq);
    elseif r == 3
        result = numMappingJustA(seq);
    end
end

function numSeq = numMappingReal(seq)
    % Real number representation
    len = length(seq);  
    numSeq = zeros(1,len,'double');
    for K = 1:len
       t = seq(K);
       if(strcmpi(t,'A'))
           numSeq(K) = -1.5;
       elseif(strcmpi(t,'C'))
           numSeq(K) = 0.5;
       elseif(strcmpi(t,'G'))
           numSeq(K) = -0.5; 
       elseif(strcmpi(t,'T'))
           numSeq(K) = 1.5;
       end           
   end   
end

function result = numMappingRandom13(seq)
    % Random among 13 representations
    a = 13;
    r = randi([1 a], 1, 1);
    if r == 1
        result = numMappingPP(seq);
    elseif r == 2
        result = numMappingReal(seq);
    elseif r == 3
        result = numMappingJustA(seq);
    elseif r == 4
        result = numMappingAtomic(seq);
    elseif r == 5
        result = numMappingEIIP(seq);
    elseif r == 6
        result = numMappingInt(seq);
    elseif r == 7
        result = numMappingAT_CG(seq);
    elseif r == 8
        result = numMappingDoublet(seq);
    elseif r == 9
        result = numMappingCodons(seq);
    elseif r == 10
       result = numMappingIntN(seq);  
    elseif r == 11
       result = numMappingJustC(seq); 
    elseif r == 12
       result = numMappingJustG(seq);  
    elseif r == 13
       result = numMappingJustT(seq); 
    end 
end