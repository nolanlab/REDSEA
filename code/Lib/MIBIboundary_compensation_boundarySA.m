%% Yunhao Bai for MIBIsubgroup, 20200421

%assumption, for adjacent cell A and BCD:
%(cell B channel X counts) = (cell B channel X counts) - (cell A channel X counts, only around the boundary) * BoundaryAB/perimeterA - cell C/D...
% + ((cell B channel X counts, only around the boundary)*BoundaryAB/perimeterB + (cell B channel X counts, only around the boundary)*BoundaryBC/perimterB +...

%% main function
function MIBIdataNorm2 = MIBIboundary_compensation_boundarySA(newLmod,data,countsNoNoiseRec,channelNormIdentity,elementShape,elementSize,REDSEAChecker)
%newLmod: the cell label image, boundary pixels are 0 while each cell has different label from 1 to cell number.
%MIBIdata: matrix of cellNum by channelNum, contains all signal counts of each cells.
%countsNoNoiseRec: imageWidth X imageLength X channelNum, signal with background removed and denoised, can be read from denoised Tiff
%channelNormIdentity: passed from MibiExtractSingleCellDataFromSegmentation3.m, this vector labels which
%channels need to be normalized as 1, while others are 0.

[rowNum, colNum] = size(newLmod);
cellNum = max(max(newLmod));
channelNum = size(data,2);

%this code use a pad to avoid boundary issues padarray(Image,[niche niche],0,'both')
niche = 1;
newLmod_padded = padarray(newLmod,[niche niche],0,'both');

%% calculate the boundary length of each cell pair
cellPairMap = zeros(cellNum,cellNum);
for x = 1+niche:rowNum+niche
    for y = 1+niche:colNum+niche
        if newLmod_padded(x,y) == 0
            tempMatrix = reshape(newLmod_padded(x-1:x+1,y-1:y+1),[9 1]);
            tempFactors = unique(tempMatrix);
            %the theoretical extreme situation is 5, in a Four Courners
            %style, and this do appear quite frequently, so consider 3-5

            %for this part, we triple/quadruple count one pixel 3/4 times
            %for each pair.
            if length(tempFactors) >= 3 
                %only assign 1 times, as unique() sort the sequence,
                %remember that the tempFactors(2) < tempFactors(3)
                cellPairMap = cellPairMapUpdater(cellPairMap,tempFactors);
            end
        end
    end
end

%flip the matrix to make it double-direction
cellPairMap = cellPairMap+cellPairMap';
%sum the matrix to get the total boundary length of each cells
cellBoundaryTotal = sum(cellPairMap,1);
%divide by the total number of cell boundary to get the fractions
cellBoundaryTotalMatrix = repmat(cellBoundaryTotal',[1 cellNum]);
cellPairNorm = REDSEAChecker*eye(cellNum) - cellPairMap./cellBoundaryTotalMatrix;


%% calculate the signals from pixels along the boundary of the cells
MIBIdataNearEdge1 = zeros(cellNum,channelNum);

%define the structural element
if elementShape == 1
    strElement = strel('square',2*elementSize+1);
elseif elementShape == 2
    strElement = strel('diamond',elementSize);
end

%bw image of the newLmod for morphological operation
newLmod_BW = imbinarize(newLmod,1/cellNum);
newLmod_BW = 1 - newLmod_BW; %reverse the direction

%subtract the boundary regions, from morphological operation
boundaryRegions = imdilate(newLmod_BW,strElement) - newLmod_BW;

%loop through each cells to get the boundary signals
for i = 1:cellNum
    [tempRow,tempCol] = find(newLmod==i);
    for j = 1:length(tempRow)
        if boundaryRegions(tempRow(j),tempCol(j)) == 1
            MIBIdataNearEdge1(i,:) = MIBIdataNearEdge1(i,:) + reshape(countsNoNoiseRec(tempRow(j),tempCol(j),:),[1 channelNum]);
        end
    end
end
%% combine that two matrix with data to give final compensated data
%original data matrix is a cellNum by channelNum matrix, so 
MIBIdataNorm2 = (MIBIdataNearEdge1'*cellPairNorm)';
MIBIdataNorm2 = MIBIdataNorm2 + data; %reinforce with boundary signal of itself
MIBIdataNorm2(MIBIdataNorm2<0) = 0;

%flip the channelNormIdentity for calculation
rev_channelNormIdentity = ones(length(channelNormIdentity),1) - channelNormIdentity;
%composite the normalized channels with non-normalized channels
MIBIdataNorm2 = data.*repmat(rev_channelNormIdentity,[1 cellNum])' + MIBIdataNorm2.*repmat(channelNormIdentity,[1 cellNum])';
end

%% auxiliary function 
function cellPairMap = cellPairMapUpdater(cellPairMap,tempFactors)
%functionalize the cellPair map updating process
    cellPairs = nchoosek(tempFactors(2:end),2); 
    %the function will order the pairs with elements consequentially
    for i = 1:size(cellPairs,1)
        cellPairMap(cellPairs(i,1),cellPairs(i,2)) = cellPairMap(cellPairs(i,1),cellPairs(i,2)) + 1;
    end
end