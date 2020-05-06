

<p align="left"><img width=20%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/redsea.jpg"></p>

# REDSEA
We present **RE**inforcement **D**ynamic **S**pillover **E**limin**A**tion (REDSEA) as a solution for spillover compensation without loss of signal. This is the step-by-step guidance of how to use REDSEA to produced a FCS file with comepensated channel information.

## Table of content

- [Required Inputs](#required-inputs)
    - [MATLAB scripts](#matlab-scripts)
    - [tiff images and channel information](#tiff-images-and-channel-information)
    - [Segmentation Mask](#segmentation-mask)

- [REDSEA](#probe-designing-showcase)
    - [Part 1 arb](#part-1-arb)
    - [Part 2 R](#part-2-r)
    - [Part 3 Optional multiple probe design](#part-3-optional-multiple-probe-design)
- [F&Q](#f&q)    

## Required Inputs



### MATLAB scripts
The MATLAB scripts include the REDSEA algorithm, along with other scripts created by [Leeat Keren](https://github.com/lkeren/MIBIAnalysis).

### tiff images and channel information

### Segmentation Mask
This method requires you to provide a cell segmentation mask. A cell nuclei probability matrix can be produced by your own chose (popular options include [ilastik](https://www.ilastik.org/) or [deepcell](https://github.com/vanvalenlab/deepcell-tf)). In our case we implemented a trained-in-house deepcell CNN model.

<p align="left"><img width=50%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/probability_matrix.png"></p>

After producing a nuclei probablity mask, we will then use the script ```MibiSegmentByDeepProbWithPerim3.m``` to implement a watershed algorithm for whole cell segmentation. This will produce something like this:

<p align="left"><img width=50%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/watershed-result.png"></p>

<details><summary>See MibiSegmentByDeepProbWithPerim3.m Script</summary>
<p>



```MATLAB
%% Pipeline for nuclear segmentation using pixel probabilities from deepCell
% Changed the pipeline to use the boundaries from the deepcell as boudaries
% instead of the raw nuclear intensities.

t1 = clock;
% maxs = imextendedmax(probNuc,0.015);
% probNuc>0.05

path = 'inputs/dsDNA.tiff';
%path = '1FOV/Histone H3.tiff';
resultsPath = 'watershed_result';
%resultsPath = 'segmentated_H3_1FOV_resize_back';
%deepPath = 'deepcell_1FOV';
deepPath = 'deepcell';

%segmentThres = 0.015; %change from 0.1 to 0.02
%probNucThres = 0.5; %change from 0.05 to 0.5
for segmentThres=0.01:0.005:0.02
for probNucThres=0.05:0.15:0.5
for p=1:1
    disp(['point',num2str(p)]);
    pointNumber=p;
    
    %{
    %% Get perimiters of nuclei for plotting
    % load data and get nuclear markers
    load([path,'/Point',num2str(pointNumber),'/dataDeNoiseCohort.mat']);
    imSize = (size(countsNoNoise,1));
    % sum nuclear markers to increase contrast
    nucleiChannels = {'Histone H3'};
    [tf loc] = ismember(nucleiChannels,massDS.Label);
    nucIm = sum(countsNoNoise(:,:,loc),3);
    %}
    
    % read tiff image back to matrix for prediction
    t = Tiff(path,'r');
    nucIm = read(t);
    maxv=50;
    rgb_image = MibiGetRGBimageFromMat(nucIm,maxv);

    % %% Get maxima from deep learning probabilities
    % read nuclear segmentation from Ilastik/deepcell
    probNuc = double(imread([deepPath,'/nuclei-probability-matrix.tif']));
    figure;
    imagesc(probNuc);

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

    % modify the image so that the background pixels and the extended maxima pixels are forced to be the only local minima in the image.
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

    cellPerimNewMod= bwperim(newLmod);
    cellPerimNewMod= newLmod;
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
    %close all;
end
end
end

t2=clock;
t = etime(t2,t1);
disp(['elapsed time: ',num2str(t)]);

```

</p>
</details>

