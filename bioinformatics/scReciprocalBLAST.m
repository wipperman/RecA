% scReciprocalBLAST ---   script to start working from the reciprocal blast
% data of recA (MSMEG_2723) and ku (MSMEG_5580).

clear, close all; clc;
%     clears all workspaces

basedir         = '/Users/Dave/Desktop/recA_ku_orthologSearch/';
%     defines the base directory this script uses

[MSMEG_5580, locus5580, raw5580] = ...
    xlsread([basedir 'MSMEG_5580_orthologs.xls']);
[MSMEG_2723, locus2723, raw2723] = ...
    xlsread([basedir 'MSMEG_2723_orthologs.xls']);
clear locus2723 locus5580 raw2723 raw5580;
%     import the data. the two options denoted with locus and raw are
%     superfluous and therefore deleted

MSMEG_2723 = MSMEG_2723(2:end, [9 11]);
MSMEG_5580 = MSMEG_5580(2:end, [9 11]);
%     column 9 contains the taxID for the organisms; column 11 contains the
%     information about the NCBI gene ID. The first row is
%     header-information and is superfluous.

missingTaxID2723 = find(isnan(MSMEG_2723(:, 1)));
missingTaxID5580 = find(isnan(MSMEG_5580(:, 1)));
%     there are a few genes there there is no taxonomy ID
%     specified in the file. 

MSMEG_2723(missingTaxID2723, :) = [];
MSMEG_5580(missingTaxID5580, :) = [];
%     removing all the lines that have no mention of a taxonomy ID
    clear missingTaxID2723 missingTaxID5580

[C, ia, ib]   = intersect(MSMEG_2723, MSMEG_5580);
%     intersect the data sets to allow identification of only those strains
%     that have both ku and recA. This is merely to exclude any strains
%     that, through this way of searching their genomic data, are supposed
%     to only have ku.

ku          = MSMEG_5580(ib, :);
recAwithKu  = MSMEG_2723(ia, :);
recAnoKu    = MSMEG_2723;
    recAnoKu(ia, :) = [];
recA        = MSMEG_2723;
%     parse out the data. ku maintains only those records that intersect
%     with recA, so those organisms have recA. Then there are three sets of
%     recA data: the full set (recA), the intersect with ku (recAwithKu)
%     and any recA that does not intersect with ku (recAnoKu).
    clear C ia ib MSMEG_2723 MSMEG_5580;

%% QUERY KEGG FOR GENE IDS --- TIME CONSUMING!!
base            = 'http://rest.kegg.jp/';
operation       = 'conv/';
database        = 'genes/';
dbentry         = 'ncbi-gi:';

%%
kegg_ku = cell(length(ku), 1);
for i = 1:length(ku)
    dbentry1 = strcat(dbentry, num2str(ku(i, 2)));
    kegg_ku{i} = regexpi(urlread(strcat(base, operation, database,...
        dbentry1)), '(?<=(??@dbentry1)\s+)\w+\W+\w*','match');
end; clear i

%%
kegg_recAnoKu = cell(length(recAnoKu), 1);
for i = 1:length(recAnoKu)
    dbentry1 = strcat(dbentry, num2str(recAnoKu(i, 2)));
    kegg_recAnoKu{i} = regexpi(urlread(strcat(base, operation, database,...
        dbentry1)), '(?<=(??@dbentry1)\s+)\w+\W+\w*','match');
end; clear i

%%
kegg_recAwithKu = cell(length(recAwithKu), 1);
for i = 1:length(recAwithKu)
    dbentry1 = strcat(dbentry, num2str(recAwithKu(i, 2)));
    kegg_recAwithKu{i} = regexpi(urlread(strcat(base, operation, database,...
        dbentry1)), '(?<=(??@dbentry1)\s+)\w+\W+\w*','match');
end; clear i

%% QUERY KEGG FOR KU GENE SEQUENCES
base            = 'http://rest.kegg.jp/';
operation       = 'get/';
option          = 'ntseq';

%%
emptyKu = find(cellfun(@isempty,kegg_ku));
kegg_ku(emptyKu) = [];
%     remove all the genes that did not give any hits in KEGG. 
for i = 1:length(kegg_ku)
    kuSeq = urlread(char(strcat(base, operation, kegg_ku{i}, '/', option)));
    headerFind = findstr(') atg', kuSeq);
        if numel(headerFind) == 0
            headerFind = findstr(') gtg', kuSeq);
        elseif numel(headerFind) == 0
            headerFind = findstr(') ctg', kuSeq);
        elseif numel(headerFind) == 0
            headerFind = findstr(') ttg', kuSeq);
        end
    sequenceData = kuSeq(headerFind(1):end);
    fastawrite([basedir 'kuSeq.fasta'], kegg_ku{i},...
        sequenceData);
end; clear i

%%
emptyRecAnoKu = find(cellfun(@isempty,kegg_recAnoKu));
kegg_recAnoKu(emptyRecAnoKu) = [];
%     remove all the genes that did not give any hits in KEGG. 
for i = 1:length(kegg_recAnoKu)
    recAnoKuSeq = urlread(char(strcat(base, operation,...
        kegg_recAnoKu{i}, '/', option)));
    headerFind = findstr('atg', recAnoKuSeq);
        if numel(headerFind) == 0
            headerFind = findstr('gtg', recAnoKuSeq);
        end
    sequenceData = recAnoKuSeq(headerFind(1):end);
    fastawrite([basedir 'recAnoKuSeq.fasta'], kegg_recAnoKu{i},...
        sequenceData);
end; clear i

%%
emptyRecAwithKu = find(cellfun(@isempty,kegg_recAwithKu));
kegg_recAwithKu(emptyRecAwithKu) = [];
%     remove all the genes that did not give any hits in KEGG. 
for i = 1:length(kegg_recAwithKu)
    recAwithKuSeq = urlread(char(strcat(base, operation,...
        kegg_recAwithKu{i}, '/', option)));
    headerFind = findstr('atg', recAwithKuSeq);
        if numel(headerFind) == 0
            headerFind = findstr('gtg', recAwithKuSeq);
        end
    sequenceData = recAwithKuSeq(headerFind(1):end);
    fastawrite([basedir 'recAwithKuSeq.fasta'], kegg_recAwithKu{i},...
        sequenceData);
end; clear i

%%

% 
%     %======================================================================
%     % ALIGNMENT AND CLEAN UP FOR ALL THREE PROTEINS
%     %======================================================================
% fleNseqs        = fastaread([basedir 'fleNorthologs_140120.fasta']);
% flhFseqs        = fastaread([basedir 'flhForthologs_140120.fasta']);
% %         Import the sequence dump files.

% fleNalign       = clearPaGaps(multialign(fleNseqs, 'ScoringMatrix', 'BLOSUM62'));
% fleQalign       = clearPaGaps(multialign(fleQseqs, 'ScoringMatrix', 'BLOSUM62'));
% flhFalign       = clearPaGaps(multialign(flhFseqs, 'ScoringMatrix', 'BLOSUM62'));
% %         Do the multiple alignment.
% %         In the assessment of which scoring matrix to use, various have
% %         been used and visually compared using the seqalignviewer. This
% %         was done on the full sequences of the proteins in PA14, while
% %         excluding any other multiple sequence alignment gaps.
% %         
% %         For all three proteins, BLOSUM62 visually seemed to be a
% %         consensus between either lower or higher divergence BLOSUM
% %         scoring matrices, as well as the Dayhoff and Gonnet matrices.
% %         Furthermore, BLOSUM62 is used relatively standard when doing
% %         multiple sequence alignments. 
% 
% %%
% 
%     %======================================================================
%     % CONCATENATION OF THE THREE SEQUENCES FOR FURTHER ANALYSIS
%     %======================================================================
% fleNfleQflhFconc= fleNalign;
% for i           = 1:length(fleNalign)
%     fleNfleQflhFconc(i).Sequence ...
%                 =   [fleNalign(i).Sequence fleQalign(i).Sequence...
%                         flhFalign(i).Sequence];
% end; clear i
% %         This for loop creates a horizontal concatenation of the different
% %         aligned cleaned up sequences. This should be useable afterwards
% %         for calculations of mutual information and potentially EVfold
% %         purposes. 
% 
% c = vertcat(fleNfleQflhFconc.Sequence);
% %         Vertically concatenate said sequences.
% 
% %% TIME-CONSUMING CALCULATION!!
% 
% I = [];
% for i = 1:size(c, 2)
%     fprintf('%d of %d\n', i, size(c, 2));
%     for j = 1:size(c, 2)
%         I(i, j) = MutualInformation(c(:, i), c(:, j));
%     end
% end
% %         Here a very large xy-matrix is created containing the mutual
% %         information for all the sequences. This is cordoned off because
% %         this calculation is very long. Furthermore, it seems not to work
% %         if it is in the same section as the previous operations.
% 
% %%
% 
% fastawrite([basedir 'micatoTest.fasta'], fleNfleQflhFconc);
% 
% %%
% load([basedir...
%     'scReadAlignmentFleN/micatoMatrixVis/micato/micato-0.8.9/micatoOutput.txt'])
% 
% %
% miNormalized = zeros(1199);
% miNormalized(tril(ones(1199), -1)==1) = I(tril(ones(1199), -1)==1);
% for i = 1:size(micatoOutput, 1)
%     miNormalized(micatoOutput(i, 1)+1, micatoOutput(i, 2)+1) =...
%         micatoOutput(i, 3);
% end