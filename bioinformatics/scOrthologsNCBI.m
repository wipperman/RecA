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
% recA        = MSMEG_2723;
%     parse out the data. ku maintains only those records that intersect
%     with recA, so those organisms have recA. Then there are three sets of
%     recA data: the full set (recA), the intersect with ku (recAwithKu)
%     and any recA that does not intersect with ku (recAnoKu).
    clear C ia ib MSMEG_2723 MSMEG_5580;

%% QUERY NC
baseURL         = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=';
restURL         = '&rettype=fasta&retmode=text';
%     part FOUR of the URL to query, returning the accession codes
%         NOTE: the URL to submit for the query will consist of the base
%         URL, followed by the accID part, then the GI numbers with commas
%         in between them, and finished by rettype.
% samples2search = 1:length(recAnoKu);
GInumbers   = sprintf('%d,', recAnoKu(1:400, 2));
GInumbers   = char(GInumbers(1:end-1));
searchValue = strcat(baseURL, GInumbers, restURL);
% 
k = websave([basedir 'recAnoKu_protein.fasta'], searchValue); clear k