%% Yunhao Bai for MIBIsubgroup, 20200421

%assumption, for adjacent cell A and BCD:
%(cell B channel X counts) = (cell B channel X counts) - (cell A channel X counts, only around the boundary) * BoundaryAB/perimeterA - cell C/D...
% + ((cell B channel X counts, only around the boundary)*BoundaryAB/perimeterB + (cell B channel X counts, only around the boundary)*BoundaryBC/perimterB +...

%approximation: 
%1. this code leaves out the very boundary pixels, as the stitched MIBI 
% images have boundary gaps while 

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

%% calculate the boundary length of each cell pair
cellPairMap = zeros(cellNum,cellNum);
for x = 2:rowNum-1 % leave out the boundary cases
    for y = 2:colNum-1 % leave out the boundary cases
        if newLmod(x,y) == 0
            tempMatrix = reshape(newLmod(x-1:x+1,y-1:y+1),[9 1]);
            tempFactors = unique(tempMatrix);
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
cellPairNorm = REDSEAChecker*eye(cellNum) - cellPairMap./cellBoundaryTotalMatrix;


%% calculate the signals from pixels along the boundary of the cells
MIBIdataNearEdge1 = zeros(cellNum,channelNum);

for i = 1:cellNum
    [tempRow,tempCol] = find(newLmod==i);
    for j = 1:length(tempRow)
        tempM = [-1]; %ignore the very boundary signals by initiating tempM that can never be reached
        if elementSize == 1
            if (1<tempRow(j))&&(tempRow(j)<rowNum)&&(1<tempCol(j))&&(tempCol(j)<colNum)
                if elementShape == 1
                    %M1, Sudoku style, 3x3, 8 pixels considered
                    tempM = reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-1:tempCol(j)+1),[9 1]);
                elseif elementShape == 2
                    %M2, star style, 3x3 4 pixels considered
                    tempM = [reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)),[1 3]),...
                        newLmod(tempRow(j),tempCol(j)+1), newLmod(tempRow(j),tempCol(j)-1)];
                end
            end
        elseif elementSize == 2
            if (2<tempRow(j))&&(tempRow(j)<rowNum-1)&&(2<tempCol(j))&&(tempCol(j)<colNum-1)
                if elementShape == 1
                    %M3, Sudoku style, 5x5, 24 pixels considered
                    tempM = reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)-2:tempCol(j)+2),[25 1]);
                elseif elementShape == 2
                    %M4, star style, 5x5 12 pixels considered
                    tempM = [reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)),[1 5]),...
                        reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)+1),[1 3]),...
                        reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-1),[1 3]),...
                        newLmod(tempRow(j),tempCol(j)+2), newLmod(tempRow(j),tempCol(j)-2)];
                end
            end
        elseif elementSize == 3
            if (3<tempRow(j))&&(tempRow(j)<rowNum-2)&&(3<tempCol(j))&&(tempCol(j)<colNum-2)
                if elementShape == 1
                    %M5, Sudoku style, 7x7, 48 pixels considered
                    tempM = reshape(newLmod(tempRow(j)-3:tempRow(j)+3,tempCol(j)-3:tempCol(j)+3),[49 1]);
                elseif elementShape == 2
                    %M6, star style, 7x7, 24 pixels considered
                    tempM = [reshape(newLmod(tempRow(j)-3:tempRow(j)+3,tempCol(j)),[1 7]),...
                        reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)+1),[1 5]),...
                        reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)-1),[1 5]),...
                        reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)+2),[1 3]),...
                        reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-2),[1 3]),...
                        newLmod(tempRow(j),tempCol(j)+3),newLmod(tempRow(j),tempCol(j)-3)];
                end
            end
        elseif elementSize == 4
            if (4<tempRow(j))&&(tempRow(j)<rowNum-3)&&(4<tempCol(j))&&(tempCol(j)<colNum-3)
                if elementShape == 1
                    %M7, Sudoku style, 9x9, 80 pixels considered
                    tempM = reshape(newLmod(tempRow(j)-4:tempRow(j)+4,tempCol(j)-4:tempCol(j)+4),[81 1]);
                elseif elementShape == 2
                    %M8, star style, 9x9, 40 pixels considered
                    tempM = [reshape(newLmod(tempRow(j)-4:tempRow(j)+4,tempCol(j)),[1 9]),...
                        reshape(newLmod(tempRow(j)-3:tempRow(j)+3,tempCol(j)+1),[1 7]),...
                        reshape(newLmod(tempRow(j)-3:tempRow(j)+3,tempCol(j)-1),[1 7]),...
                        reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)+2),[1 5]),...
                        reshape(newLmod(tempRow(j)-2:tempRow(j)+2,tempCol(j)-2),[1 5]),...
                        reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)+3),[1 3]),...
                        reshape(newLmod(tempRow(j)-1:tempRow(j)+1,tempCol(j)-3),[1 3]),...
                        newLmod(tempRow(j),tempCol(j)+4),newLmod(tempRow(j),tempCol(j)-4)];
                end
            end    
        end
        tempFactors = unique(tempM);
        if tempFactors(1) == 0
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
