function [minVal, maxVal] = recommendContrastRange(I, lowerPercentile, upperPercentile)
% -------------------------------------------------------------------------
% recommendContrastRange: Suggests optimal min and max intensity values 
%                         for maximizing image contrast, based on percentiles.
%
% Inputs:
%   I                : A 2D or 3D grayscale image (single channel)
%   lowerPercentile  : Lower cutoff percentile (e.g., 1)
%   upperPercentile  : Upper cutoff percentile (e.g., 99)
%
% Outputs:
%   minVal, maxVal   : Recommended intensity range [minVal, maxVal] derived 
%                      from the specified percentiles, to optimally enhance 
%                      the image contrast when displaying or scaling 'I'.
%
% Example:
%   [minVal, maxVal] = recommendContrastRange(I, 1, 99);
%   figure;
%   imshow(I, [minVal, maxVal]); % Display image with the recommended contrast range
%
% -------------------------------------------------------------------------


    data = I(:);
    data(~isfinite(data)) = [];

    minVal = prctile(data, lowerPercentile);
    maxVal = prctile(data, upperPercentile);

    if minVal >= maxVal
        warning('minVal >= maxVal 상태 발생! global [min, max]로 대체합니다.');
        minVal = min(data);
        maxVal = max(data);
    end
end
