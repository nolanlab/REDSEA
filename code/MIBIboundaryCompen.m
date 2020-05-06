%% Yunhao Bai for MIBIsubgroup, 20200430
%main function of boundary compensation 

%{
%add this part to the MibiExtractSingleCellDataFromSegmentation3.m before
the for loop to make the channelNormIdentity() vector

%selected channels to do the boundary compensation
normChannels = {'dsDNA';'Histone H3';'CD16';'SMA';'CD209 (DC-SIGN)';'NFkB-p100 (pS865)';'CD4';'CD11c';'CD56';'FoxP3';'CD39';'Granzyme B';'Biotin';'CD21 (CR2)';'Ki-67';'PD-1';'Pax-5';'CCR7';'CD163';'CD68';'FoxO1';'CD8';'CD3';'CD45-RA';'Dinitrophenyl (DNP)';'CD86';'CTLA-4';'CD20';'Lamin AC';'MPO';'HLA-DR';'IL10';'CD169 (Sialoadhesin)';'CD8a';'Pan-Keratin';'CD11b';'CD36';'Digoxigenin (DIG)';'CD25';'CD45'};
[~, normChannelsInds] = ismember(normChannels,massDS.Label);
channelNormIdentity = zeros(length(clusterChannels),1);
for i = 1:length(normChannelsInds)
    channelNormIdentity(normChannelsInds(i)) = 1;
end
%}

%besides, add 'Centroid' to the regionprops() parameters if we want to map
%stats = regionprops(newLmod,'Area','PixelIdxList','Centroid');


%% main function
dataNorm = MIBIboundary_compensation3(newLmod,data,countsNoNoiseRec,channelNormIdentity);

%add additional 1 line at the beginning as the cell label
[cellNum,channelNum] = size(data);
dataCopy = zeros(cellNum,channelNum+1);
dataNormCopy = zeros(cellNum,channelNum+1);

dataScaledCopy = zeros(cellNum,channelNum+1);
dataScaledNormCopy = zeros(cellNum,channelNum+1);

%add a cell label at the begining column
for i = 1:cellNum
    if labelIdentityNew2(i) == 1
        dataCopy(i,:) = [i,data(i,:)];
        dataNormCopy(i,:) = [i,dataNorm(i,:)];
        
        dataScaledCopy(i,:) = [i,data(i,:)./cellSizes(i)];
        dataScaledNormCopy( = [i,dataNorm(i,:)./cellSizes(i)];
    end
end

%remove 0 values (labeled by labelIdentityNew2 == 0
dataCopy = dataCopy(dataCopy(:,1) ~= 0,:);
dataNormCopy = dataNormCopy(dataNormCopy(:,1) ~= 0,:);
dataScaledCopy = dataScaledCopy(dataScaledCopy(:,1) ~= 0,:);
dataScaledNormCopy = dataScaledNormCopy(dataScaledNormCopy(:,1) ~= 0,:);

%write the csv files
csvwrite('20200430_Slide99_SIV_LNs.csv',dataCopy);
csvwrite('20200430_Slide99_SIV_LNs_Norm.csv',dataNormCopy);
csvwrite('20200430_Slide99_SIV_LNs_Scaled.csv',dataScaledCopy);
csvwrite('20200430_Slide99_SIV_LNs_ScaledNorm.csv',dataScaledNormCopy);
