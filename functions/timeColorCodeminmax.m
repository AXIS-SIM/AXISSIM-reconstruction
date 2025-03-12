function [colorStack, MIP] = timeColorCode(stack, colormapUsed, intensityMin, intensityMax)
    % ----------------------------------------------------
    % Time-lapse color coding with a user-defined intensity range.
    %
    % Inputs:
    %   stack         : A 3D array (H x W x N) representing single-channel frames/slices
    %   colormapUsed  : A colormap of size (N x 3), where N is the number of frames/slices
    %   intensityMin  : Desired minimum intensity value
    %   intensityMax  : Desired maximum intensity value
    %
    % Outputs:
    %   colorStack : A 4D array (H x W x 3 x N) containing RGB color-coded frames
    %   MIP        : A 3D array (H x W x 3) for the maximum intensity projection of colorStack
    % ----------------------------------------------------


    [H, W, nFrames] = size(stack);


    colorStack = zeros(H, W, 3, nFrames);

    for i = 1:nFrames

        frame = stack(:,:,i);

        frameNorm = (frame - intensityMin) / (intensityMax - intensityMin);
        frameNorm(frameNorm < 0) = 0;
        frameNorm(frameNorm > 1) = 1;


        color = colormapUsed(i,:);  % [R_i, G_i, B_i]
        R = frameNorm * color(1);
        G = frameNorm * color(2);
        B = frameNorm * color(3);


        colorStack(:,:,:,i) = cat(3, R, G, B);
    end

    MIP = squeeze(max(colorStack, [], 4)); 
end
