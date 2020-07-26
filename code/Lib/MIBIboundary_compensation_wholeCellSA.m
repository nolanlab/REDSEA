%% Yunhao Bai for MIBIsubgroup, 20200420
%compensation the boundary leakover with CODEX-like logic
%(cell B channel X counts) = (cell B channel X counts) 
%- (cell A channel X counts) * BoundaryAB/perimeterA - cell C/D
%+ (cell B channel X counts)

%assumption and approximation: 
%1. this code leaveout the very boundary pixels, which is just 1 or 2 for each cells.
%2. at corner of 3/4 different cells boundary, this code will triple/quadruple count that pixel.

%% 
function MIBIdataNorm = MIBIboundary_compensation_wholeCellSA(newLmod, MIBIdata, channelNormIdentity,REDSEAChecker)
%newLmod: the cell label image, boundary pixels are 0 while each cell has
%different label from 1 to cell number.
%MIBIdata: matrix of cellNum by channelNum, contains all signal counts of
%each cells.
%channelNormIdentity: passed from
%MibiExtractSingleCellDataFromSegmentation3.m, this vector labels which
%channels need to be normalized as 1, while others are 0.

[rowNum, colNum] = size(newLmod);
cellNum = max(max(newLmod));
cellPairMap = zeros(cellNum,cellNum);
%factors = -1;
for x = 2:rowNum-1 % we just leave out the boundary cases
    for y = 2:colNum-1 % we just leave out the boundary cases
        if newLmod(x,y) == 0
            tempM = reshape(newLmod(x-1:x+1,y-1:y+1),[9 1]);
            tempFactors = unique(tempM);
            %the theoretical extreme situation is 5, in a cross style, and
            %this do appear quite frequently, so consider 3-5

            %for this part, we triple/quadruple count one pixel 3/4 times
            %for each pair.
            if length(tempFactors) == 3 
                %only assign 1 times, as unique() sort the sequence,
                %remember that the tempFactors(2) < tempFactors(3)
                cellPairMap(tempFactors(2),tempFactors(3)) = cellPairMap(tempFactors(2),tempFactors(3))+1;
            elseif length(tempFactors) == 4
                cellPairMap(tempFactors(2),tempFactors(3)) = cellPairMap(tempFactors(2),tempFactors(3))+1;
                cellPairMap(tempFactors(2),tempFactors(4)) = cellPairMap(tempFactors(2),tempFactors(4))+1;
                cellPairMap(tempFactors(3),tempFactors(4)) = cellPairMap(tempFactors(3),tempFactors(4))+1;
            elseif length(tempFactors) == 5
                cellPairMap(tempFactors(2),tempFactors(3)) = cellPairMap(tempFactors(2),tempFactors(3))+1;
                cellPairMap(tempFactors(2),tempFactors(4)) = cellPairMap(tempFactors(2),tempFactors(4))+1;
                cellPairMap(tempFactors(2),tempFactors(5)) = cellPairMap(tempFactors(2),tempFactors(5))+1;
                
                cellPairMap(tempFactors(3),tempFactors(4)) = cellPairMap(tempFactors(3),tempFactors(4))+1;
                cellPairMap(tempFactors(3),tempFactors(5)) = cellPairMap(tempFactors(3),tempFactors(5))+1;
                
                cellPairMap(tempFactors(4),tempFactors(5)) = cellPairMap(tempFactors(4),tempFactors(5))+1;
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
cellPairNorm = (REDSEAChecker+1)*eye(cellNum) - cellPairMap./cellBoundaryTotalMatrix;

%flip the channelNormIdentity for calculation
channelNum = length(channelNormIdentity);
rev_channelNormIdentity = ones(channelNum,1) - channelNormIdentity;

%original data matrix is a cellNum by channelNum matrix, so 
MIBIdataNorm = (MIBIdata'*cellPairNorm)';
MIBIdataNorm(MIBIdataNorm<0) = 0;

%composite the normalized channels with non-normalized channels
MIBIdataNorm = MIBIdata.*repmat(rev_channelNormIdentity,[1 cellNum])' + MIBIdataNorm.*repmat(channelNormIdentity,[1 cellNum])';
end
