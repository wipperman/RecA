% scReciprocalBLASTdb ---   this script is devised to attempt making one
%                           database out of the sequences coming from
%                           reciprocal BLAST of fleN, fleQ and flhF. 

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                              --------
%                              !!NOTE!!
%                              --------
% 
% Please run this script section per section. There are a few
% time-consuming calculations in separate sections, of which it would be
% wise to avoid running them time and time again.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% IMPORTING THE DATA
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
clear, close all; clc;

    %======================================================================
    % DEFINE THE BASE DIRECTORY
    %======================================================================
basedir         = '/Users/Dave/Desktop/recA_ku_orthologSearch/';

    %======================================================================
    % IMPORT THE DIFFERENT FILES
    %======================================================================
recA            = readOrthologFile([basedir 'MSMEG_2723_orthologs.txt']);
ku              = readOrthologFile([basedir 'MSMEG_5580_orthologs.txt']);

%         The function that is called upon here will automatically parse
%         the entire datafile and just output the high quality calls from
%         the different Ortholuge search results. The output will only have
%         the project IDs for the bugs and their different protein IDs.

%     %======================================================================
%     % INTERSECTING THE DATA TO MAINTAIN SPECIES PRESENT IN ALL DATA SETS
%     %======================================================================
% [C, ia, ib]     = intersect(fleN(:, 1), flhF(:, 1));
% fleNflhF        = [fleN(ia, 1:2) flhF(ib, 2)];
% %                 clear C ia ib fleN flhF
% %         Intersect fleN and flhF
% 
%     %======================================================================
%     % QUERYING THE NCBI
%     %======================================================================
% baseURL         = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=';
% restURL         = '&rettype=fasta&retmode=text';
% %         Define the general part of the URL to query the NCBI using their
% %         own E-Utilities.
% 
% concFleN        = sprintf('%d,',fleNflhF(:,2));
% concFleN        = ['116049399,' concFleN(1:end-1)];
% concFlhF        = sprintf('%d,',fleNflhF(:,3));
% concFlhF        = ['116049398,' concFlhF(1:end-1)];
% %         Variable part for fleN, fleQ and flhF, respectively
% 
% searchFleN      = [baseURL concFleN restURL];
%                 urlwrite(searchFleN, [basedir 'fleNorthologs_140120.fasta']);
% searchFlhF      = [baseURL concFlhF restURL];
%                 urlwrite(searchFlhF, [basedir 'flhForthologs_140120.fasta']);
%                 clear concFleN concFlhF baseURL restURL...
%                     fleNflhF searchFleN searchFlhF
% %         The arguments above here use the E-Utilities of the NCBI to
% %         create a big URL querying their database for the sequences of the
% %         proteins within each data set. The returned item will be one
% %         sequence dump of fasta files.
% %         Said fasta files will include the full headers, which not only
% %         have the bug name, but also the protein description and the
% %         unique identifier, as well as, obviously, the sequence. 
% %         urlwrite is chosen to output the files immediately, as this
% %         maintains the formatting of the fasta file, which would otherwise
% %         be lost. 
% 
% % %% TIME-CONSUMING CALCULATION!!
% % 
% %     %======================================================================
% %     % ALIGNMENT AND CLEAN UP FOR ALL THREE PROTEINS
% %     %======================================================================
% % fleNseqs        = fastaread([basedir 'fleNorthologs_140120.fasta']);
% % flhFseqs        = fastaread([basedir 'flhForthologs_140120.fasta']);
% % %         Import the sequence dump files.
% 
% % fleNalign       = clearPaGaps(multialign(fleNseqs, 'ScoringMatrix', 'BLOSUM62'));
% % fleQalign       = clearPaGaps(multialign(fleQseqs, 'ScoringMatrix', 'BLOSUM62'));
% % flhFalign       = clearPaGaps(multialign(flhFseqs, 'ScoringMatrix', 'BLOSUM62'));
% % %         Do the multiple alignment.
% % %         In the assessment of which scoring matrix to use, various have
% % %         been used and visually compared using the seqalignviewer. This
% % %         was done on the full sequences of the proteins in PA14, while
% % %         excluding any other multiple sequence alignment gaps.
% % %         
% % %         For all three proteins, BLOSUM62 visually seemed to be a
% % %         consensus between either lower or higher divergence BLOSUM
% % %         scoring matrices, as well as the Dayhoff and Gonnet matrices.
% % %         Furthermore, BLOSUM62 is used relatively standard when doing
% % %         multiple sequence alignments. 
% % 
% % %%
% % 
% %     %======================================================================
% %     % CONCATENATION OF THE THREE SEQUENCES FOR FURTHER ANALYSIS
% %     %======================================================================
% % fleNfleQflhFconc= fleNalign;
% % for i           = 1:length(fleNalign)
% %     fleNfleQflhFconc(i).Sequence ...
% %                 =   [fleNalign(i).Sequence fleQalign(i).Sequence...
% %                         flhFalign(i).Sequence];
% % end; clear i
% % %         This for loop creates a horizontal concatenation of the different
% % %         aligned cleaned up sequences. This should be useable afterwards
% % %         for calculations of mutual information and potentially EVfold
% % %         purposes. 
% % 
% % c = vertcat(fleNfleQflhFconc.Sequence);
% % %         Vertically concatenate said sequences.
% % 
% % %% TIME-CONSUMING CALCULATION!!
% % 
% % I = [];
% % for i = 1:size(c, 2)
% %     fprintf('%d of %d\n', i, size(c, 2));
% %     for j = 1:size(c, 2)
% %         I(i, j) = MutualInformation(c(:, i), c(:, j));
% %     end
% % end
% % %         Here a very large xy-matrix is created containing the mutual
% % %         information for all the sequences. This is cordoned off because
% % %         this calculation is very long. Furthermore, it seems not to work
% % %         if it is in the same section as the previous operations.
% % 
% % %%
% % 
% % fastawrite([basedir 'micatoTest.fasta'], fleNfleQflhFconc);
% % 
% % %%
% % load([basedir 'scReadAlignmentFleN/micatoMatrixVis/micato/micato-0.8.9/micatoOutput.txt'])
% % 
% % %
% % miNormalized = zeros(1199);
% % miNormalized(tril(ones(1199), -1)==1) = I(tril(ones(1199), -1)==1);
% % for i = 1:size(micatoOutput, 1)
% %     miNormalized(micatoOutput(i, 1)+1, micatoOutput(i, 2)+1) =...
% %         micatoOutput(i, 3);
% % end