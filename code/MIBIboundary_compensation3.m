%% Yunhao Bai for MIBIsubgroup, 20200421
%compensation the boundary leakover with CODEX-like logic, but modified

%for cell AB contact
%(cell B channel X counts) = (cell B channel X counts) - (cell A channel X counts, only around the boundary) * BoundaryAB/perimeterA - cell C/D...
% + ((cell B channel X counts, only around the boundary)/perimeterB)*BoundaryAB + cell B * BoundaryBC/perimterB +...
%(cell A chennel X counts) = (cell A channel X counts) - (cell B channel X counts, only around the boundary) * BoundaryAB/perimeterB - cell C/D...


%assumption and approximation: 
%1. this code leaveout the very boundary pixels, which is just 2 pixels for boundary cells.

%% main function
function MIBIdataNorm2 = MIBIboundary_compensation3(newLmod,MIBIdata,countsNoNoiseRec,channelNormIdentity)

%newLmod: the cell label image, boundary pixels are 0 while each cell has different label from 1 to cell number.
%MIBIdata: matrix of cellNum by channelNum, contains all signal counts of each cells.
%countsNoNoiseRec: 
%channelNormIdentity: passed from MibiExtractSingleCellDataFromSegmentation3.m, this vector labels which
%channels need to be normalized as 1, while others are 0.

[rowNum, colNum] = size(newLmod);
cellNum = max(max(newLmod));

%% calculate the boundary length of each cell pair
cellPairMap = zeros(cellNum,cellNum);
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
cellPairNorm = eye(cellNum) - cellPairMap./cellBoundaryTotalMatrix;
%flip the channelNormIdentity for calculation
channelNum = length(channelNormIdentity);
rev_channelNormIdentity = ones(channelNum,1) - channelNormIdentity;

%% calculate the signals that n=1 pixel from the boundary of the cells
MIBIdataNearEdge1 = zeros(cellNum,channelNum);

for i = 1:cellNum
    [tempRow,tempCol] = find(newLmod==i);
    for j = 1:length(tempRow)
        if (2<tempRow(j))&&(tempRow(j)<rowNum-1)&&(2<tempCol(j))&&(tempCol(j)<colNum-1)
            %Sudoku style, 3x3, 8 pixels considered
            %tempM = reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-1:tempCol(j)+1),[9 1]);
            
            %cross style, 3x3 4 pixels considered
            %tempM = [reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)),[1 3]), newLmod(tempRow(j),tempCol(j)+1), newLmod(tempRow(j),tempCol(j)-1)];
            
            %Sudoku style, 5x5, 24 pixels considered
            %need to change the boundary conditions
            %tempM = reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)-2:tempCol(j)+2),[25 1]);
            
            %cross style, 5x5 12 pixels considered
            %need to change the boundary conditions
            %tempM = reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-1:tempCol(j)+1),[13 1]);
            
            tempM = [reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)),[1 5]),...
                reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)+1),[1 3]),...
                reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-1),[1 3]),...
                newLmod(tempRow(j),tempCol(j)+2), newLmod(tempRow(j),tempCol(j)-2)];
            
            tempFactors = unique(tempM);
            if tempFactors(1) == 0
                MIBIdataNearEdge1(i,:) = MIBIdataNearEdge1(i,:) + reshape(countsNoNoiseRec(tempRow(j),tempCol(j),:),[1 channelNum]);
            end
        end
    end
end
%% combine that two matrix with data to give final compensated data
%original data matrix is a cellNum by channelNum matrix, so 
MIBIdataNorm2 = (MIBIdataNearEdge1'*cellPairNorm)';
MIBIdataNorm2 = MIBIdataNorm2 + MIBIdata;
MIBIdataNorm2(MIBIdataNorm2<0) = 0;

%composite the normalized channels with non-normalized channels
MIBIdataNorm2 = MIBIdata.*repmat(rev_channelNormIdentity,[1 cellNum])' + MIBIdataNorm2.*repmat(channelNormIdentity,[1 cellNum])';
end
