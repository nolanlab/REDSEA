function rgb_data = MibiGetRGBimageFromMat(data,maxv)
% rgb_data = MibiGetRGBimageFromMat(data,maxv)
% function gets a data matrix and returns an rgb matrix where colors are
% scaled according to the intensity in the data.
% maxv - an optional parameter for capping the data matrix. Otherwise the
% max value in the matrix is used

if ~exist('maxv') || isempty(maxv)
    maxv = max(nucIm(:));
end

map = colormap;
minv = min(data(:));

ncol = size(map,1);
s = round(1+(ncol-1)*(data-minv)/(maxv-minv));
rgb_data = ind2rgb(s,map);