function [newLabel] = MibiMergeLabels (L, labelID, labelTargetID)
% function [newLabel] = MibiMergeLabels (L, labelID, labelTargetID)
% function receives a label matrix L. It merges the pixels with label
% labelID into labelTargetID

object1 = L == labelID;
object2 = L == labelTargetID;
newLabel = L;

%% 1. Change label pixels to the new values
newLabel(object1) = labelTargetID;

%% 2. Erase border between labels
se = ones(3);   % 8-connectivity for neighbours - could be changed

% get border of label 1
object1Wborder = imdilate(object1, se);

% get border of label 2
object2Wborder = imdilate(object2, se);

% get joint border 
jointBorder = object1Wborder & object2Wborder;
% if the joint border also partitions another label leave it as is
jointBorderDilate = imdilate(jointBorder, se);
jointBorderDilate((newLabel ~= labelTargetID) & (newLabel ~=0)) = 0;
jointBorderFix = imerode(jointBorderDilate, se);

borderToDelete = jointBorder & jointBorderFix;

newLabel(borderToDelete) = labelTargetID;

