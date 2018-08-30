% scClusteringPairwiseDistances --- this script will plot the dendrogram
% for the pairwise distance information that Matt obtained for the recA
% sequences.

clc; clear, close all;

basedir = '/Users/Dave/Desktop/recA_ku_orthologSearch/pairwiseDistance/';
%     define the base directory

distanceMatrix = importdata([basedir '151202.csv']);
    distanceMatrix.textdata = distanceMatrix.textdata(2:end, 1:2);
%     import the data

attributes = importdata([basedir 'attributeTable.xlsx']);
    attributes.AttributeFile = attributes.AttributeFile(2:end, :);
%     import the attribute file. Only 'AttributeFile' is of importance.

% %% THIS IS TO CREATE THE CLUSTERGRAM FROM THE DATA
% distanceMatrix.data = flipud(distanceMatrix.data);
% distanceMatrix.textdata = flipud(distanceMatrix.textdata);
% %     data flipped upside down
% 
% data4clustering = squareform(distanceMatrix.data(:, 1));
% %     make the data into a square matrix for clustering
% CGobj = clustergram(data4clustering);
% %     standard clustergram, making an object that is not easily manipulated

%% ATEMPT TO RECREATE THE DENDROGRAM USING LINKAGE AND DENDROGRAM COMMANDS
distanceMatrix.data = flipud(distanceMatrix.data);
distanceMatrix.textdata = flipud(distanceMatrix.textdata);
%     data flipped upside down; this is done to get out the proper format
%     for creating the squareform matrix downstream. it also means that all
%     subsequent data will need to be reorganized. 

labels = flipud(unique(distanceMatrix.textdata));
%     try to get the labels out for the tree and flip for compatibility

data4clustering = squareform(distanceMatrix.data(:, 1));
%     make the data into a square matrix for clustering 

attributes.AttributeFile = flipud(sortrows(attributes.AttributeFile, 1));
%     alphabetically order the attribute file and flip upside down to be
%     compatible with the rest of the data

NHEJ = strfind(attributes.AttributeFile(:, 2), 'NHEJ');
% 	  denote all NHEJ organisms with 1
NHEJ = ~cellfun(@isempty, NHEJ);
%     replace all empty cells with 0
NHEJ = repmat(NHEJ, 1, 50);

res207 = NaN(length(attributes.AttributeFile), 1);
for i = 1:length(attributes.AttributeFile)
    if strcmp('S', attributes.AttributeFile(i, 3)) == 1;
        res207(i, 1) = 0;
    elseif strcmp('N', attributes.AttributeFile(i, 3)) == 1;
        res207(i, 1) = 1;
    else
        res207(i, 1) = 2;
    end
end; clear i
%     create a vector with different values depending on what is considered
%     residue 207 (S, N, or other).
res207 = repmat(res207, 1, 50);

tree = linkage(data4clustering, 'average', 'euclidean');
%     make the data for the tree. The tree will be made in the exact same
%     fashion as the standard clustergram dendrogram is made. This means
%     'average' linkage and 'euclidian' distance metric.

leafOrder = optimalleaforder(tree, data4clustering);
%     optimize the leaf order.

%%
figure;
    subplot (20, 1, 1:8 )
        dendrogram(tree, 0, 'Labels', [], 'Reorder', leafOrder,...
            'colorthreshold', 75);
            ax = gca;
%             ax.XTickLabelRotation = 90; % make sure the labels are rotated
           set(gca, 'YTick', [], 'XTick', [], 'Box', 'off', 'Color', 'none');
           axis off;
    h = subplot (20, 1, 9:18);
        imagesc(data4clustering(leafOrder, leafOrder));
            set(gca, 'YTick', [], 'XTick', [], 'Color', 'none');
            colormap(h, jet);
    i = subplot (20, 1, 19);
        imagesc(NHEJ(leafOrder))
            set(gca, 'YTick', [], 'XTick', [], 'Color', 'none');
            colormap(i, gray);
    j = subplot (20, 1, 20);
        imagesc(res207(leafOrder))
            set(gca, 'YTick', [], 'XTick', [], 'Color', 'none');
            colormap(j, parula);
clear h i j

%%
N = 0; NKu = 0; S = 0; SKu = 0;
for i = 1:length(attributes.AttributeFile)
    if strcmp(attributes.AttributeFile(i, 3), 'N') &&...
            strcmp(attributes.AttributeFile(i, 2), 'NHEJ')
        NKu = NKu + 1;
    end
    if strcmp(attributes.AttributeFile(i, 3), 'N')
        N = N + 1;
    end
    if strcmp(attributes.AttributeFile(i, 3), 'S') &&...
            strcmp(attributes.AttributeFile(i, 2), 'NHEJ')
        SKu = SKu + 1;
    end
    if strcmp(attributes.AttributeFile(i, 3), 'S')
        S = S + 1;
    end
end; clear i
            
%% SAMPLE FROM JOAO
% tree = linkage(Y, 'complete');
% leafOrder = optimalleaforder(tree, Y);
%  
% figure(1)
% 
% % this plots the tree
% subplot(2, 1, 1)
% dendrogram(tree, 0, 'Labels', labels, 'Reorder?, leafOrder) 
% ax = gca;
% ax.XTickLabelRotation = 90; % make sure the labels are rotated
% set(gca, 'YTick', []) % in my case I didn?t want any text on the y axis
% 
% % this plots the data using 'imagesc'
% subplot(2, 1, 2)
% imagesc(correlationMatrix(leafOrder, leafOrder))
% title('correlation between samples')
% set(gca, 'XTick', [], 'YTick', [])
% title('clustering before sample correction')
