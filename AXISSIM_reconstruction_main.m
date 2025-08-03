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
% This algorithm is based on W. Zhao et al.'s work with modifications.
% Please also cite:
%  "Enhanced detection of fluorescence fluctuations for high-throughput 
%   super-resolution imaging," Nat. Photon. 17, 806–813 (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INITIALIZATION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% Add necessary paths (modify these paths as needed)
addpath('./functions');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USER-DEFINED PARAMETERS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Data processing parameters ---
layer_num   = 100;    % Number of frames per layer 
total_layers = 26;    % Total number of layers to process

% --- Imaging parameters ---
pixel_size_nm  = 104;    % Camera pixel size (nm)
wavelength_nm  = 515;    % emission wavelength (nm)
magnification  = 3;      % Magnification factor for upscaling
dz             = 40;     % Axial step size (nm)
RI             = 1.338;  % Refractive index of medium

% --- Deconvolution parameters ---
iter_preRL     = 7;  % Number of iterations for pre-deconvolution
iter_postRL    = 8;  % Number of iterations for post-deconvolution (recommendation: 8-9)
cumulant_order = 3;  % Cumulant order for SACD reconstruction
scale_factor   = 1;  % Scaling factor for PSF exponentiation

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOAD DATA & PRE-PROCESSING
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Loading TIFF files...");

firstFile = sprintf("./AXIS-SIM data/U2OS/U2OS_MicroTubule_AXIS_SIM_01.tif");
info = imfinfo(firstFile);
imgWidth = info(1).Width;  
imgHeight = info(1).Height; 

% Rawdata (512x512x100x26, uint16)
Rawdata = zeros(imgHeight, imgWidth, layer_num, total_layers, 'uint16');

for i = 1:total_layers
    fileName = sprintf("./AXIS-SIM data/U2OS/U2OS_Microtubule_AXIS_SIM_%02d.tif", i);
    fprintf("Reading: %s\n", fileName);
    for j = 1:layer_num
        Rawdata(:,:,j,i) = imread(fileName, j);
    end
end

Rawdata = double(Rawdata);

disp("Data successfully loaded into Rawdata array.");

% Compute the SUM projection for each layer (sum along the frame dimension)
WF = zeros(size(Rawdata,1), size(Rawdata,2), total_layers);
for i = 1:total_layers
    WF(:,:,i) = sum( Rawdata(:,:,:,i), 3 );
end

% Compute RMS projection for each layer using rms across frames
RMS = zeros(size(Rawdata,1), size(Rawdata,2), total_layers);
for i = 1:total_layers
    RMS(:,:,i) = std( Rawdata(:,:,:,i), 0,3 ); 
end

% Normalize images
normalize = @(q) (q - min(q(:))) / (max(q(:)) - min(q(:))) * 65535;
WF_norm = normalize(WF);
RMS_norm = normalize(RMS);

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

    RMS_deconv = std(deconvData, 0,3); 
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
    
    interpDSI = fourierInterpolation(squeeze(datadeconDSI3D_linear(:,:,z,:)), [magnification, magnification, 1], 'lateral'); interpDSI(interpDSI < 0) = 0;    
    interpDSI = abs(interpDSI - mean(interpDSI, 3) * subfactor); interpDSI(interpDSI < 0) = 0;    
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
RMS_decon = deconvlucy(RMS_norm,psf_central, iter_postRL); % RL-Deconvolved RMS
AXISSIMresult = deconvlucy(cumDataDSI, psf_central_rescaled.^cumulant_order, iter_postRL); % Raw AXIS-SIM images
AXISSIM_norm = normalize(AXISSIMresult);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3D Linearization and Denoising 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AXISSIM_decon = deconvlucy(AXISSIM_norm,psf_central_rescaled, iter_postRL);
AXISSIM_decon_linear = (AXISSIM_decon).^(1/cumulant_order);

% -------------------------------
% User-defined parameters (Adjustable)
do_bg_subtraction = true;  
gaussian_filter_sigma = 300;  % Gaussian filter sigma value for background estimation
bg_threshold = -0.00;          % Background subtraction threshold
% -------------------------------

AXISSIM_decon_linear_bg = AXISSIM_decon_linear;

if do_bg_subtraction
    % --- Background subtraction using Gaussian filtering ---
    for z = 1:total_layers
        AXISSIM_decon_linear_bg(:,:,z) = ...
            AXISSIM_decon_linear_bg(:,:,z) ...
            - imgaussfilt(AXISSIM_decon_linear(:,:,z), gaussian_filter_sigma);
    end
    
    AXISSIM_decon_linear_bg(AXISSIM_decon_linear_bg < bg_threshold) = bg_threshold;
    AXISSIM_decon_linear_bg = AXISSIM_decon_linear_bg - bg_threshold;
    
    disp('Background subtraction performed.');
else
    disp('Background subtraction skipped.');    
end

% Convolution for final reconstruction
for z = 1:total_layers
    AXISSIM_decon_linearization(:,:,z) = convnfft(AXISSIM_decon_linear_bg(:,:,z), psf_rescaled.^cumulant_order, 'same');
end

RMS_decon_norm = normalize(RMS_decon);
AXISSIM_decon_linearization_norm = normalize(AXISSIM_decon_linearization);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DISPLAY RESULTS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Preparing for visualization...");

% Parameters
px = pixel_size_nm;         % original XY pixel size (nm)
px_corr = px / 3;           % corrected pixel size after upsampling
scale_len_nm = 5000;        % scale bar in nm (5 µm)
bar_px = scale_len_nm / px;      
bar_px_corr = scale_len_nm / px_corr;
acf = tan(asin(0.5*1.49/1.518))/tan(asin(0.5*1.49/RI));
dz_corr = dz/acf;

% Cropping coordinates (normalized and AXIS)
xm = 1050; ym = 1110; sm = 420; z_idx = 6;
x0n = floor(xm/3)+1; y0n = floor(ym/3)+1; w = floor(sm/3);
x1n = x0n+w-1; y1n = y0n+w-1;
x0a = xm+1; y0a = ym+1; x1a = xm+sm; y1a = ym+sm;

% Cropped images
rms_crop  = RMS_norm(x0n:x1n, y0n:y1n, z_idx);
rms_cropd = RMS_decon_norm(x0n:x1n, y0n:y1n, z_idx);
axis_crop = AXISSIM_decon_linearization_norm(x0a:x1a, y0a:y1a, z_idx);

% Line for plotting
x1 = 295; y1 = 138; x2 = 353; y2 = 303;
x1a = 1405; y1a = 1188; x2a = 1463; y2a = 1353;

% Visualization
figure('Color','w','Position',[100 100 1200 800]);

subplot(2,3,1); imagesc(rms_crop); axis image off; colormap hot;
title('DSI (layer 6)'); hold on;
xbar = 10; ybar = size(rms_crop,1) - 10;
line([xbar, xbar+bar_px], [ybar, ybar], 'Color','w', 'LineWidth',5);
text(xbar+13, ybar-15, '5 µm', 'Color','w', 'FontSize',12, 'FontWeight','bold');
plot([round(x1/3), round(x2/3)], [round(y1/3), round(y2/3)], 'w--', 'LineWidth',2);

subplot(2,3,2); imagesc(rms_cropd); axis image off; colormap hot;
title('DSI decon (layer 6)'); hold on;
plot([round(x1/3), round(x2/3)], [round(y1/3), round(y2/3)], 'w--', 'LineWidth',2);

subplot(2,3,3); imagesc(axis_crop); axis image off; colormap hot;
title('AXIS-SIM (layer 6)'); hold on;
plot([x1, x2], [y1, y2], 'w--', 'LineWidth',2);

fprintf('XY pixel: %.2f nm, Z pixel: %.2f nm\n', px_corr, dz_corr);

% Generate line profile coordinates
n_xy = max(round(abs(x2a - x1a)/3), round(abs(y2a - y1a)/3)) + 1;
n_ax = max(abs(x2a - x1a), abs(y2a - y1a)) + 1;
xLine = linspace(round(x1a/3), round(x2a/3), n_xy);
yLine = linspace(round(y1a/3), round(y2a/3), n_xy);
xLine2 = linspace(x1a, x2a, n_ax);
yLine2 = linspace(y1a, y2a, n_ax);

% Extract axial cross-section
nz = size(RMS_norm,3);
line_rms  = zeros(n_xy, nz);
line_rmsd = zeros(n_xy, nz);
line_axis = zeros(n_ax, nz);

for z = 1:nz
    line_rms(:,z)  = interp2(double(RMS_norm(:,:,z)), xLine,  yLine);
    line_rmsd(:,z) = interp2(double(RMS_decon_norm(:,:,z)), xLine,  yLine);
    line_axis(:,z) = interp2(double(AXISSIM_decon_linearization_norm(:,:,z)), xLine2, yLine2);
end

% Axial view plots
subplot(2,3,4);
imshow(imrotate(imresize(line_rms, [round(size(line_rms,1)/(dz_corr/px)), total_layers], 'nearest'), 90), []);
axis equal tight; title('Axial View - DSI'); xlabel('Z-slice');

subplot(2,3,5);
imshow(imrotate(imresize(line_rmsd, [round(size(line_rmsd,1)/(dz_corr/px)), total_layers], 'nearest'), 90), []);
axis equal tight; title('Axial View - DSI Decon'); xlabel('Z-slice');

subplot(2,3,6);
imshow(imrotate(imresize(line_axis, [round(size(line_axis,1)/(dz_corr/px)), total_layers*3], 'nearest'), 90), []);
axis equal tight; title('Axial View - AXIS-SIM'); xlabel('Z-slice');
colormap hot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
