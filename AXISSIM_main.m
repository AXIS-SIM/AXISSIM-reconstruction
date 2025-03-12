%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example Code for:
% "Near-isotropic Super-Resolution Microscopy with Axial Interference Speckle Illumination"
%
% This script demonstrates the AXIS-SIM reconstruction pipeline, including:
%   - Data loading and pre-processing
%   - Layer extraction, summation, and RMS computation
%   - PSF loading and pre-deconvolution
%   - SACD reconstruction using cumulant analysis and deconvolution
%   - Optional 3D linearization and Denoising 
%
% Parameter descriptions are provided inline for user modification.
%
% This algorithm is based on W. Zhao et al.'s work with significant modifications.
% Please also cite:
%  "Enhanced detection of fluorescence fluctuations for high-throughput 
%   super-resolution imaging," Nat. Photon. 17, 806–813 (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INITIALIZATION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% Add necessary paths (modify these paths as needed)
addpath('./AXIS-SIM data');
addpath('./functions');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USER-DEFINED PARAMETERS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Data processing parameters ---
layer_num   = 100;    % Number of frames per layer (acquisition parameter)
total_layers = 26;    % Total number of layers to process
save_images = false;  % Set to 'true' to save output images

% --- Imaging parameters ---
pixel_size_nm  = 104;    % Camera pixel size (in nm)
wavelength_nm  = 515;    % emission wavelength (in nm)
magnification  = 3;      % Magnification factor for upscaling
dz             = 40;     % Axial step size (in nm)
RI             = 1.338;  % Refractive index of medium

% --- Deconvolution parameters ---
iter_preRL     = 7;  % Number of iterations for pre-deconvolution
iter_postRL    = 8;  % Number of iterations for post-deconvolution
cumulant_order = 3;  % Cumulant order for SACD reconstruction
scale_factor   = 1;  % Scaling factor for PSF exponentiation

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOAD DATA & PRE-PROCESSING
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Loading data...");
load("AXIS-SIM data/Fig3a_MicroTubule_AXIS_SIM.mat"); % Load dataset

% Compute the SUM projection for each layer (sum along the frame dimension)
WF = zeros(size(Rawdata,1), size(Rawdata,2), total_layers);
for i = 1:total_layers
    WF(:,:,i) = sum( Rawdata(:,:,:,i), 3 );
end

% Compute RMS projection for each layer using rms across frames
RMS = zeros(size(Rawdata,1), size(Rawdata,2), total_layers);
for i = 1:total_layers
    RMS(:,:,i) = rms( Rawdata(:,:,:,i), 3 );
end

% Normalize images
normalize = @(q) (q - min(q(:))) / (max(q(:)) - min(q(:))) * 65535;
WF_norm = normalize(WF);
RMS_norm = normalize(RMS);

% Save processed images (if enabled)
if save_images
    for i = 1:total_layers
        imwrite(uint16(WF_norm(:,:,i)), 'AXIS-SIM results/Fig3a_MicroTubule_WF.tif', 'WriteMode', 'append');
        imwrite(uint16(RMS_norm(:,:,i)), 'AXIS-SIM results/Fig3a_MicroTubule_RMS.tif', 'WriteMode', 'append');
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PSF LOADING & PRE-DECONVOLUTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Loading PSF...");
psf_filename = sprintf('./AXIS-SIM data/PSF BW_r%.3f_w%d_xy%d_z%d.tif', RI, wavelength_nm, pixel_size_nm, dz);
psf_filename_rescaled = sprintf('./AXIS-SIM data/PSF BW_r%.3f_w%d_xy%d_z%d.tif', RI, wavelength_nm, round(pixel_size_nm/magnification), dz);
psf_3DBW = double(tiffreadVolume(psf_filename));
psf_3DBW_rescaled = double(tiffreadVolume(psf_filename_rescaled));

% Select the central slice as PSF for deconvolution
psf = psf_3DBW(:,:,65);
psf_rescaled = psf_3DBW_rescaled(:,:,65);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PRE-DECONVOLUTION (LUCY-RICHARDSON) & DSI reweighting
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Performing pre-deconvolution...");

for z = 1:total_layers
    current_layer = Rawdata(:,:,:,z);
    deconvData = zeros(size(current_layer));

    for i = 1:size(current_layer,3)
        deconvData(:,:,i) = deconvlucy(current_layer(:,:,i), psf, iter_preRL);
    end

    RMS_deconv = rms(deconvData, 3);
    for i = 1:layer_num
        deconvDSI(:,:,i) = RMS_deconv .* deconvData(:,:,i);
    end

    datadeconDSI3D(:,:,z,:) = deconvDSI;
    disp(['Layer ', num2str(z), ' pre-deconvolution complete. ']);
end

datadeconDSI3D_linear = sqrt(datadeconDSI3D);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SACD RECONSTRUCTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subfactor = 0.8;  % Subtraction factor for background correction
disp("Performing SACD reconstruction...");

for z = 1:total_layers
    disp(['Processing cumulant analysis for layer ', num2str(z)]);
    
    interpDSI = fourierInterpolation(squeeze(datadeconDSI3D_linear(:,:,z,:)), [magnification, magnification, 1], 'lateral');
    interpDSI = abs(interpDSI - mean(interpDSI, 3) * subfactor);
    cumDataDSI(:,:,z) = abs(cumulant(interpDSI, cumulant_order));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  POST-DECONVOLUTION (LUCY-RICHARDSON)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Performing 3D deconvolution...");
% Center PSF slice selection
z_trim = round(size(psf_3DBW_rescaled,3) / 2) - round(size(cumDataDSI,3)/2);  
psf_central = psf_3DBW(:, :, (z_trim + 1):(z_trim + floor(size(cumDataDSI,3)/2)*2-1));
psf_central_rescaled = psf_3DBW_rescaled(:, :, (z_trim + 1):(z_trim + floor(size(cumDataDSI,3)/2)*2-1));
RMSresult_decon = deconvlucy(RMS_norm,psf_central, iter_postRL); % RL-Deconvolved RMS
SACDDSI3dresult = deconvlucy(cumDataDSI, psf_central_rescaled.^cumulant_order, iter_postRL); % Raw AXIS-SIM images
SACDDSI3d_norm = normalize(SACDDSI3dresult);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3D Linearization and Denoising 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SACDDSI3dresult_decon = deconvlucy(SACDDSI3d_norm,psf_central_rescaled, iter_postRL);
SACDDSI3dresult_decon_linear = (SACDDSI3dresult_decon).^(1/cumulant_order);

% -------------------------------
% User-defined parameters (Adjustable)
do_bg_subtraction = false;  
gaussian_filter_sigma = 100;  % Gaussian filter sigma value for background estimation
bg_threshold = -0.0;          % Background subtraction threshold
% -------------------------------

SACDDSI3dresult_decon_linear_bg = SACDDSI3dresult_decon_linear;

if do_bg_subtraction
    % --- Background subtraction using Gaussian filtering ---
    for z = 1:total_layers
        SACDDSI3dresult_decon_linear_bg(:,:,z) = ...
            SACDDSI3dresult_decon_linear_bg(:,:,z) ...
            - imgaussfilt(SACDDSI3dresult_decon_linear(:,:,z), gaussian_filter_sigma);
    end
    
    SACDDSI3dresult_decon_linear_bg(SACDDSI3dresult_decon_linear_bg < bg_threshold) = bg_threshold;
    SACDDSI3dresult_decon_linear_bg = SACDDSI3dresult_decon_linear_bg - bg_threshold;
    
    disp('Background subtraction performed.');
else
    disp('Background subtraction skipped.');    
end

% Convolution for final reconstruction
for z = 1:total_layers
    SACDDSI3dresult_decon_linear_conv2d(:,:,z) = convnfft(SACDDSI3dresult_decon_linear_bg(:,:,z), psf_rescaled.^cumulant_order, 'same');
end

RMSresult_decon_norm = normalize(RMSresult_decon);
SACDDSI3dresult_decon_linear_conv2d_norm = normalize(SACDDSI3dresult_decon_linear_conv2d);

if save_images
    for i = 1:total_layers
        imwrite(uint16(RMSresult_decon_norm(:,:,i)), ['AXIS-SIM results\Fig3a_MicroTubule_RMS_decon.tif'],'WriteMode','append' );
    end 
    for i = 1:total_layers
        imwrite(uint16(SACDDSI3dresult_decon_linear_conv2d_norm(:,:,i)), ['AXIS-SIM results\Fig3a_MicroTubule_AXISSIM_linearization.tif'],'WriteMode','append' );
    end 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DISPLAY RESULTS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Preparing for visualization...");

display_layer = 19;
contrast_min_per = 25;
contrast_max_per = 100;
colormap_name = 'darkjet';  
colormap_used = eval(sprintf('%s(%d)', colormap_name, display_layer));

stack1 = []; stack2 = []; stack3 = [];
for i = 1:display_layer
    display_layer0 = 3;
    stack1(:,:,i) = RMS_norm(:,:, i+display_layer0);
    stack2(:,:,i) = RMSresult_decon_norm(:,:, i+display_layer0);
    stack3(:,:,i) = SACDDSI3dresult_decon_linear_conv2d_norm(:,:, i+display_layer0);
end
[minVal, maxVal] = recommendContrastRange(stack1, contrast_min_per, contrast_max_per);
[colorStack1, MIP1] = timeColorCodeminmax(stack1, colormap_used, minVal, maxVal);
[minVal, maxVal] = recommendContrastRange(stack2, contrast_min_per, contrast_max_per);
[colorStack2, MIP2] = timeColorCodeminmax(stack2, colormap_used, minVal, maxVal);
[minVal, maxVal] = recommendContrastRange(stack3, contrast_min_per, contrast_max_per);
[colorStack3, MIP3] = timeColorCodeminmax(stack3, colormap_used, minVal, maxVal);

% Display results
disp("Displaying visualization...");
figure('Position', [200, 200, 1200, 500]);
% Full-field images
subplot(2,3,1); imshow(MIP1); title('Color-coded Image (DL)');
subplot(2,3,2); imshow(MIP2); title('Color-coded Image (Deconvolved DL)');
subplot(2,3,3); imshow(MIP3); title('Color-coded Image (AXIS-SIM)');

% Zoomed-in images
xm = 1050; ym = 1110; sm = 420; sym = 420;
subplot(2,3,4); imshow(MIP1(xm/magnification:xm/magnification+sm/magnification, ym/magnification:ym/magnification+sym/magnification, :));
title('Zoomed-in Image (DL)');
subplot(2,3,5); imshow(MIP2(xm/magnification:xm/magnification+sm/magnification, ym/magnification:ym/magnification+sym/magnification, :));
title('Zoomed-in Image (Deconvolved DL)');
subplot(2,3,6); imshow(MIP3(xm:xm+sm, ym:ym+sym, :));
title('Zoomed-in Image (AXIS-SIM)');
disp("Visualization complete.");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
