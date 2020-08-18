%% Nuclear segmentation using pixel probabilities from deepCell
% Modified from original script from Leeat Keren
% Changed the pipeline to use the probablities from the deepcell instead of 
% the raw nuclear intensities.
% 25July2020, by Yunhao Bai and Sizun Jiang

% Main path for the all the data
mainPath = 'sampleData_MIBI'; %for MIBI
%mainPath = 'sampleData_cycIF'; %for CyCIF
resultsPath = [mainPath,'/segmentResults/']; % for the next step, put the 
% output segmentationParams.mat and PureSegmentation.tif to the
% originalTiff/Point? folder

for p=1:1
for segmentThres=0.01 % change from 0.01 to 0.1, segmentThres defines the 
% local maximum, higher segmentThres leads to fewer cells and more merged
% regions
for probNucThres=0.35 % change from 0.05 to 0.5, probNucThres defines how 
% board the cells will expand
    pointNumber=p;disp(['point',num2str(p)]);
    
    % read .tiff image of nucleus marker to matrix for image  
    pathNucleusMarker = [mainPath,'/originalTIFF/Point',num2str(p),'/dsDNA.tiff']; %for MIBI
    %pathNucleusMarker = [mainPath,'/originalTIFF/Point',num2str(p),'/x7500y3500_1700_DAPI.tif']; %for CyCIF
    t = Tiff(pathNucleusMarker,'r');
    nucIm = read(t);
    %the max value of nucleus channel .tiff
    maxv=25; %for MIBI
    %maxv=3000; %for CyCIF
    rgb_image = MibiGetRGBimageFromMat(nucIm,maxv);

    % %% Get maxima from deepCell probabilities
    % read possibility map from deepCell/ilastik/other segmentation methods
    probNuc = double(imread([mainPath,'/deepCell/feature_1_frame_1_p',num2str(p),'_dsDNA.tif'])); %for MIBI
    %probNuc = double(imread([mainPath,'/deepCell/feature_1_frame_1_p',num2str(p),'_DAPI.tif'])); %for CyCIF
    figure;imagesc(probNuc);

    % find local maxima in probability map
    maxs = imextendedmax(probNuc,segmentThres); % change from 0.1 to 0.02
    rgb_image_perim_extMax = imoverlay(rgb_image , maxs, [1 0 0]);   
    figure;
    imagesc(rgb_image_perim_extMax);
    
    %% watershed over the deep results
    bw1 = zeros(size(probNuc));
    bw1(probNuc>probNucThres) = 1; %change from 0.05 to 0.5
    figure; imagesc(bw1);
    bw2 = bwareaopen(bw1,40);
    SE = strel('disk',4);
    bw=imdilate(bw2,SE);
    
    [B,L] = bwboundaries(maxs,4,'noholes');
    maxsFix = bw & maxs;

    % modify the image so that the background pixels and the extended maxima 
    % pixels are forced to be the only local minima in the image.
    Jc = imcomplement(probNuc);
    I_mod = imimposemin(Jc, ~bw | maxsFix);
    L = watershed(I_mod);
    labeledImage = label2rgb(L);
    
    cellPerimNewMod= L;
    cellPerimNewMod(L>0) = 100;
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    
    rgb_image_cellPerim = imoverlay(rgb_image, cellPerimNewMod, [1 0 0]);
    figure;
    imagesc(rgb_image_cellPerim);

    %% 1. For each label, decide whether it is of a nucleus/ background
    t=40;
    labelNum = length(unique(L(:)));
    labelIdentity = zeros (labelNum,1);
    labelPixelsPercentInNucleiMask = zeros(labelNum,1);
    labelSize = zeros(labelNum,1);
    
    %spend a lot of time during this loop 
    for i=1:labelNum
        [r c] = find(L==i);
        labelSize(i) = length(r); 
        labelMask = (L==i);
        labelPixelsNumInNucleiMask = sum(bw(labelMask));
        labelPixelsPercentInNucleiMask(i) = labelPixelsNumInNucleiMask / labelSize(i);
        if (labelPixelsPercentInNucleiMask(i) > 0.7)
            labelIdentity(i) = 1;
        end
    end

    % 2. Merge small regions within the nuclei mask with their neighbours
    keepVec = ones(labelNum,1);
    newL = L;
    for i=1:labelNum
        if (labelIdentity(i) == 1) && (labelSize(i) < t)
            disp(['Removing label ',num2str(i),'. Size: ',num2str(labelSize(i))]);
            % get neighbour with largest border that is also in nuclear region
            [neighbourLabels , neighbouringRegionSize] = MibiGetNeighbourLabels (newL, i);
            found = 0;
            [neighbouringRegionSizeS , neighbouringRegionSizeSInd] = sort(neighbouringRegionSize,'descend');
            neighbourLabelsS = neighbourLabels(neighbouringRegionSizeSInd);
            maxInd = 1;
            while ~found
                mergeLabelId = neighbourLabelsS(maxInd);
                if (~(mergeLabelId == 0) && (labelIdentity(mergeLabelId) == 1))
                    found = 1;
                else
                    maxInd = maxInd+1;
                end
                if (maxInd >length(neighbourLabelsS)) % reached end of neighbours with no good merging candidate
                    disp (['Warning: no good merging target found for label', num2str(i), '. Keeping it.']);
                    break;
                end
            end
            % update
            if (maxInd <= length(neighbourLabelsS))
                [newL] = MibiMergeLabels (newL, i, mergeLabelId);
                keepVec(i) = 0;
            end
        end
    end

    % Update label numbers to account for deleted labels
    allLabels = [1:labelNum];
    currLabels =  allLabels(keepVec == 1);
    labelIdentityNew = zeros(length(currLabels),1); 
    newLmod = newL;
    for i = 1:length(currLabels)
        newLmod(newLmod == currLabels(i)) = i;
        labelIdentityNew(i) = labelIdentity(currLabels(i));
    end

    cellPerimNewMod = bwperim(newLmod);
    cellPerimNewMod = newLmod;
    cellPerimNewMod(newL>0) = 100;
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    rgb_image_cellPerimNewMod = imoverlay(rgb_image, cellPerimNewMod, [1 0 0]);
    % add ilastik segmentation
    % ilsegParams = load([path,'/SegmentPerim/Point',num2str(pointNumber),'/segmentationParams.mat']);
    % rgb_image_cellPerimNewModIl = imoverlay(rgb_image, ilsegParams.cellPerimNewMod, [0 0 0]);
    % rgb_image_cellPerimNewModBoth = imoverlay(rgb_image_cellPerimNewModIl, cellPerimNewMod, [1 0 0]);
    
    figure;
    imagesc(rgb_image_cellPerimNewMod);
    mkdir([resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres)]);
    imwrite(rgb_image_cellPerimNewMod,[resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres),'/compareIlastikDeep.tif'],'tif');
    imwrite(cellPerimNewMod,[resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres),'/PureSegmentation.tif'],'tif');
    save([resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres),'/segmentationParams.mat'],'newLmod','cellPerimNewMod','labelIdentityNew');
end
end
end

disp('Please put the PureSegmentation.tif and segmentationParams.mat to the corresponding originalTiff folder for next step.');

