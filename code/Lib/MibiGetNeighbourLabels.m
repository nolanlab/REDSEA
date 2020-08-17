function [neighbourLabels , neighbouringRegionSize] = MibiGetNeighbourLabels (L, labelID)
% function [neighbourLabels , neighbouringRegionSize] = MibiGetNeighbourLabels (L, labelID)
% The function receives a matrix of labeled objects L and a specific label
% ID labelID. It returns neighbourLabels, a vector of indexes of labels
% adjacent to L

object = L == labelID;
se = ones(3);   % 8-connectivity for neighbours - could be changed
objectWborder = imdilate(object, se);
neighbours = imdilate(objectWborder, se) & ~objectWborder;
neighbourLabelsAll = L(neighbours);
a=histogram(neighbourLabelsAll,[unique(neighbourLabelsAll); max(neighbourLabelsAll)+1]);
neighbourLabels = a.BinEdges([1:end-1]);
neighbouringRegionSize = a.Values;