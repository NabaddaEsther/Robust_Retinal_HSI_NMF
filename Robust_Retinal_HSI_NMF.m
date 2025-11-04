%% QE_MNF(adaptiveK)_NMF_450_770_SGpresmooth_displayOnly.m
clear all
close all
clc

% -----------------------------------------------------------
% Global plotting style (bigger labels, ticks, titles)
% -----------------------------------------------------------
set(groot, 'defaultAxesFontSize',                16);   % ticks + default label size baseline
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.25); % label > tick size
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.35); % title > label size
set(groot, 'defaultAxesLineWidth',               1.2);  % thicker axes box
set(groot, 'defaultLineLineWidth',               1.6);  % thicker plotted lines by default
set(groot, 'defaultLegendFontSize',              14);
set(groot, 'defaultColorbarFontSize',            14);

% -----------------------------------------------------------
% User Configuration
% -----------------------------------------------------------
cube_path = fullfile('data','example_cube.mat'); % expects spectr (X×Y×B), wavelengths (1×B)
analysis_window = [450 770];         % analyze ONLY this wavelength range (nm)

% Savitzky–Golay pre-smoothing (BEFORE QE, BEFORE MNF)
use_sg          = false;             % toggle SG smoothing on/off
sg_frame        = 11;                % odd window length (bands). Try 9–15.
sg_order        = 3;                 % polynomial degree (<= frame-1), usually 2–3

% QE file (CSV or XLSX/XLS)
qe_file         = 'QE_Hamamatsu_BackThinned_450_700nm.csv';
apply_qe        = true;              % toggle QE correction
qe_low_thresh   = 0.05;              % mask bands where QE < 5% to avoid noise blow-up

% NMF / QC
num_components  = 4;                 % desired NMF components
binsize         = 4;                 % QC binning (for binned QC plot)
nmf_binsize     = 8;                 % spectral bin size for smoothing NMF H (viz only) & binned RGB at end

% MNF options (Adaptive K, hard truncation)
use_mnf         = true;              % toggle MNF denoise
retain_frac     = 0.99999999;          % keep 99.9999% of MNF cumulative eigenvalues
minK_factor     = 2;                 % at least 2× NMF components retained (safety floor)

% -----------------------------------------------------------
% Load Data
% -----------------------------------------------------------
data = load(cube_path);              % expects: spectr (X×Y×B), wavelengths (1×B)
cube = data.spectr;
w    = data.wavelengths(:)';         % row vector
orig_w = w;
[x, y, ~] = size(cube);

% -----------------------------------------------------------
% Optional: Show SLO Image with ROI and crop cube to ROI (if available)
% -----------------------------------------------------------
if isfield(data,'slo_image') && isfield(data,'slo_roi')
    figure; imshow(data.slo_image, []);
    drawrectangle('Position', data.slo_roi, 'StripeColor','r','LineWidth',1.5,'HandleVisibility','off');
    title('SLO Image with ROI','FontSize',18,'FontWeight','bold');

    roi = data.slo_roi; % [x0, y0, width, height]
    rowStart = max(1, floor(roi(2))); colStart = max(1, floor(roi(1)));
    rowEnd   = min(x, ceil(roi(2) + roi(4) - 1));
    colEnd   = min(y, ceil(roi(1) + roi(3) - 1));
    if rowStart < rowEnd && colStart < colEnd
        cube = cube(rowStart:rowEnd, colStart:colEnd, :);
        [x, y, ~] = size(cube);
    else
        warning('ROI appears invalid or outside bounds. Proceeding without cropping.');
    end
elseif isfield(data,'slo_roi')
    roi = data.slo_roi;
    rowStart = max(1, floor(roi(2))); colStart = max(1, floor(roi(1)));
    rowEnd   = min(x, ceil(roi(2) + roi(4) - 1));
    colEnd   = min(y, ceil(roi(1) + roi(3) - 1));
    if rowStart < rowEnd && colStart < colEnd
        cube = cube(rowStart:rowEnd, colStart:colEnd, :);
        [x, y, ~] = size(cube);
    else
        warning('ROI appears invalid or outside bounds. Proceeding without cropping.');
    end
end

% -----------------------------------------------------------
% Remove hot pixels (per-band median filtering + outlier replace)
% -----------------------------------------------------------
sd  = std(cube, 0, [1 2]);
thr = 5 * median(sd);
cube_clean = zeros(size(cube), 'like', cube);
for i = 1:size(cube,3)
    img = cube(:,:,i);
    img_median = medfilt2(img, [3 3]);
    noisemap = imabsdiff(img, img_median) > thr;
    img(noisemap) = img_median(noisemap);
    cube_clean(:,:,i) = img;
end
cube = cube_clean;

% -----------------------------------------------------------
% Spectral filtering: ANALYZE ONLY analysis_window
% -----------------------------------------------------------
in_window = (orig_w >= analysis_window(1)) & (orig_w <= analysis_window(2));
if ~any(in_window)
    error('No wavelengths fall inside the requested analysis window [%g, %g] nm.', analysis_window(1), analysis_window(2));
end
cube = cube(:,:,in_window);
w    = orig_w(in_window);
B    = numel(w);
fprintf('Analyzing ONLY [%g, %g] nm -> kept %d bands.\n', analysis_window(1), analysis_window(2), B);

% -----------------------------------------------------------
% QC: mean spectrum BEFORE SG (for overlay)
% -----------------------------------------------------------
mean_preSG = squeeze(mean(rescale(cube), [1 2]));

% -----------------------------------------------------------
% Savitzky–Golay smoothing (BEFORE QE)
% -----------------------------------------------------------
if use_sg
    if mod(sg_frame,2)==0, sg_frame = sg_frame + 1; end
    sg_order = min(sg_order, sg_frame-1);
    if sg_frame > B
        sg_frame = max(3, 2*floor(B/2)-1); % largest odd <= B
        sg_order = min(sg_order, sg_frame-1);
        warning('SG frame reduced to %d to fit available bands.', sg_frame);
    end
    D = reshape(cube, [], B);           % [Npix x Bands]
    D = sgolayfilt(D, sg_order, sg_frame, [], 2);
    cube = reshape(D, x, y, B);
end

% -----------------------------------------------------------
% Overlay plot: BEFORE vs AFTER SG
% -----------------------------------------------------------
mean_postSG = squeeze(mean(rescale(cube), [1 2]));
figure; plot(w, mean_preSG, 'LineWidth',1.0); hold on;
plot(w, mean_postSG, 'LineWidth',1.6);
ylabel('Mean intensity'); xlabel('Wavelength (nm)');
title(sprintf('Spectral QC — BEFORE vs AFTER SG (frame=%d, order=%d)', sg_frame, sg_order), ...
      'FontSize',18,'FontWeight','bold');
legend({'Before SG','After SG'}, 'Location','best','FontSize',14);
xline(analysis_window(1),'--'); xline(analysis_window(2),'--');

% -----------------------------------------------------------
% Detector QE correction (AFTER SG)
% -----------------------------------------------------------
if apply_qe
    [qe_wl, qe_val] = load_qe_curve(qe_file);  % qe_val in fraction (0–1)
    qe_interp = interp1(qe_wl, qe_val, w, 'linear', 'extrap');

    if any(w < min(qe_wl)) || any(w > max(qe_wl))
        warning('QE: wavelengths outside QE table range; extrapolating.');
    end

    lowQE = qe_interp < qe_low_thresh;         % low QE mask
    for b = 1:B
        if ~lowQE(b)
            cube(:,:,b) = cube(:,:,b) ./ max(qe_interp(b), eps);
        else
            cube(:,:,b) = 0;                   % hard mask (avoid noise blow-up)
        end
    end

    mean_postQE = squeeze(mean(rescale(cube), [1 2]));
    figure;
    yyaxis left;  plot(w, mean_postQE, 'LineWidth',1.4); ylabel('Mean intensity (post-QE)');
    yyaxis right; plot(w, qe_interp*100, '--', 'LineWidth',1.0); ylabel('QE (%)');
    xlabel('Wavelength (nm)');
    title('QE-corrected mean spectrum and detector QE','FontSize',18,'FontWeight','bold');
    legend({'Mean spectrum (post-QE)','Detector QE'}, 'Location','best','FontSize',14);
else
    mean_postQE = []; % not used if QE off
end

% -----------------------------------------------------------
% MNF denoise (Adaptive K, hard truncation) BEFORE NMF
% -----------------------------------------------------------
if use_mnf
    minK = min(max(minK_factor*num_components, 1), B);
    cube = mnf_denoise_cube_adaptiveK(cube, retain_frac, minK);
    cube(cube < 0) = 0;

    mean_postMNF = squeeze(mean(rescale(cube), [1 2]));
    ref_curve = mean_postSG;
    if apply_qe && ~isempty(mean_postQE), ref_curve = mean_postQE; end

    figure; plot(w, ref_curve, 'LineWidth',1.0); hold on;
    plot(w, mean_postMNF, 'LineWidth',1.6);
    ylabel('Mean intensity'); xlabel('Wavelength (nm)');
    title(sprintf('Spectral QC — %s vs Post MNF (retain=%.6f, minK=%d)', ...
        ternStr(apply_qe,'Post-QE','Post-SG'), retain_frac, minK), ...
        'FontSize',18,'FontWeight','bold');
    legend({ternStr(apply_qe,'Post-QE','Post-SG'),'Post-MNF'}, 'Location','best','FontSize',14);
end

% -----------------------------------------------------------
% Create binned version for QC (binned plot)
% -----------------------------------------------------------
total_bands = size(cube,3);
nBands = floor(total_bands / binsize) * binsize;
if (nBands / binsize) < num_components
    warning('QC-binned bands (%d) < components (%d); disabling binning.', nBands/binsize, num_components);
    binsize = 1; nBands = total_bands;
end

cube_qc = cube(:,:,1:nBands);
w_qc    = w(1:nBands);
cube_temp = reshape(cube_qc, x, y, binsize, []);
cube_qc   = sum(cube_temp, 3);                 % bin spectrally (sum)
w_qc      = mean(reshape(w_qc, binsize, []), 1);
cube_qc   = rescale(cube_qc);                  % normalize [0,1] for visualization

% QC plot (binned)
spectr_mean_binned = squeeze(mean(cube_qc, [1 2]));
figure; plot(w_qc, spectr_mean_binned, 'LineWidth',1.2);
ylabel('Mean intensity (binned)'); xlabel('Wavelength (nm)');
title('Spectral QC (Binned)','FontSize',18,'FontWeight','bold');
xline(analysis_window(1),'--'); xline(analysis_window(2),'--');

% -----------------------------------------------------------
% NMF on full (unbinned) cube
% -----------------------------------------------------------
cube_nmf = rescale(cube);
reshaped  = reshape(cube_nmf, [], size(cube_nmf,3));
[~, nBandsFinal] = size(reshaped);
if nBandsFinal < num_components
    warning('Only %d spectral band(s) available; reducing components from %d to %d.', nBandsFinal, num_components, nBandsFinal);
    num_components = nBandsFinal;
end

[W, H] = nnmf(reshaped, num_components);

% Optional smoothing of H for nicer viz (only for plotting)
nBands_H = size(H,2);
nSmooth  = floor(nBands_H / nmf_binsize) * nmf_binsize;
if nSmooth >= nmf_binsize
    H_smooth = reshape(H(:,1:nSmooth), num_components, nmf_binsize, []);
    H_smooth = squeeze(mean(H_smooth,2));
    w_smooth = mean(reshape(w(1:nSmooth), nmf_binsize, []),1);
else
    H_smooth = H; w_smooth = w;
end

% -----------------------------------------------------------
% Abundance Maps
% -----------------------------------------------------------
abundance_maps = reshape(W, x, y, []);
nC   = num_components;
nCol = ceil(sqrt(nC)); nRow = ceil(nC / nCol);

figure; tiledlayout(nRow, nCol, 'Padding','compact','TileSpacing','compact');
for k = 1:nC
    nexttile; imagesc(mat2gray(abundance_maps(:,:,k)));
    axis image off; colormap gray; title(sprintf('Comp %d', k),'FontSize',16,'FontWeight','bold');
end
sgtitle('Abundance Maps (All Components)','FontSize',20,'FontWeight','bold');

% -----------------------------------------------------------
% NMF Spectral Components
% -----------------------------------------------------------
figure; plot(w_smooth, H_smooth', 'LineWidth',1.2);
xlabel('Wavelength (nm)'); ylabel('Component intensity');
title(sprintf('NMF Spectral Components (Window %d–%d nm)%s%s', ...
    analysis_window(1), analysis_window(2), ternStr(apply_qe,' + QE',''), ternStr(use_mnf,' + MNF(adaptiveK)','')), ...
    'FontSize',18,'FontWeight','bold');
legend(arrayfun(@(i) sprintf('Comp %d', i), 1:num_components, 'UniformOutput', false), ...
       'Location','best','FontSize',14);

% -----------------------------------------------------------
% Optional RGB Visualization (Unbinned cube)
% -----------------------------------------------------------
try
    if isfield(data, 'metadata'); hcube = hypercube(cube_nmf, w, data.metadata);
    else;                          hcube = hypercube(cube_nmf, w);
    end
    [rgbImg, bands] = colorize(hcube, Method="falsecolored", ContrastStretching=true);
    figure; imshow(rgbImg);
    title({'False-colored RGB Composite (Unbinned)', ...
           sprintf('Bands: %s', strjoin(string(round(w(bands))), ', '))}, ...
           'FontSize',18,'FontWeight','bold');
catch
    warning('Hypercube/colorize not available; skipping unbinned RGB.');
end

% -----------------------------------------------------------
% Binned RGB AT THE VERY END (binned-for-smoothing by nmf_binsize)
% -----------------------------------------------------------
if exist('nSmooth','var') && nSmooth >= nmf_binsize
    cube_nmf_binned = reshape(cube_nmf(:,:,1:nSmooth), x, y, nmf_binsize, []);
    cube_nmf_binned = squeeze(mean(cube_nmf_binned, 3));
    w_binned = mean(reshape(w(1:nSmooth), nmf_binsize, []), 1);
else
    nb = floor(size(cube_nmf,3)/nmf_binsize) * nmf_binsize;
    if nb >= nmf_binsize
        cube_nmf_binned = reshape(cube_nmf(:,:,1:nb), x, y, nmf_binsize, []);
        cube_nmf_binned = squeeze(mean(cube_nmf_binned, 3));
        w_binned = mean(reshape(w(1:nb), nmf_binsize, []), 1);
    else
        cube_nmf_binned = [];
    end
end

if exist('cube_nmf_binned','var') && ~isempty(cube_nmf_binned)
    try
        if isfield(data, 'metadata'); hcube_binned = hypercube(cube_nmf_binned, w_binned, data.metadata);
        else;                          hcube_binned = hypercube(cube_nmf_binned, w_binned);
        end
        [rgbImg_binned, bands_binned] = colorize(hcube_binned, Method="falsecolored", ContrastStretching=true);
        figure; imshow(rgbImg_binned);
        title({'False-colored RGB Composite (Binned for Smoothing)', ...
               sprintf('Bands: %s', strjoin(string(round(w_binned(bands_binned))), ', '))}, ...
               'FontSize',18,'FontWeight','bold');
    catch
        % Toolbox fallback: simple 3-band pseudo-RGB
        idx = max(1, round(linspace(1, size(cube_nmf_binned,3), 3)));
        rgb_simple = cat(3, mat2gray(cube_nmf_binned(:,:,idx(1))), ...
                            mat2gray(cube_nmf_binned(:,:,idx(2))), ...
                            mat2gray(cube_nmf_binned(:,:,idx(3))));
        figure; imshow(rgb_simple);
        title(sprintf('Pseudo RGB (Binned for Smoothing) ~ [%g, %g, %g] nm', ...
              w_binned(idx(1)), w_binned(idx(2)), w_binned(idx(3))), ...
              'FontSize',18,'FontWeight','bold');
    end
end

disp('Done: ROI → hot-pixel → window → SG (pre) → QE → MNF(adaptiveK) → QC-bin → NMF → RGB (unbinned) → RGB (binned-for-smoothing).');

%% ----------------------- Helpers -----------------------
function out = ternStr(cond, a, b)
% Return a or b depending on cond (for concise titles/labels)
    if cond; out = a; else; out = b; end
end

function [qe_wl, qe_val] = load_qe_curve(qe_file)
% Load QE curve with CSV/XLSX support. Returns wavelength (nm) and QE (0–1).
    if ~isfile(qe_file), error('QE file not found: %s', qe_file); end
    [~,~,ext] = fileparts(qe_file);
    switch lower(ext)
        case '.csv'
            T = readtable(qe_file);
        case {'.xlsx', '.xls'}
            T = readtable(qe_file, 'FileType','spreadsheet');
        otherwise
            % Try CSV first, then spreadsheet as fallback
            try
                T = readtable(qe_file);
            catch
                T = readtable(qe_file,'FileType','spreadsheet');
            end
    end
    varnames = lower(string(T.Properties.VariableNames));
    wl_idx = find(contains(varnames,"wavelength") | contains(varnames,"nm"), 1, 'first');
    qe_idx = find(contains(varnames,"qe") | contains(varnames,"quantum") | contains(varnames,"efficiency"), 1, 'first');
    if isempty(wl_idx), wl_idx = 1; end
    if isempty(qe_idx), qe_idx = 2; end
    qe_wl  = T{:, wl_idx};
    qe_val = T{:, qe_idx};
    mask = isfinite(qe_wl) & isfinite(qe_val);
    qe_wl  = qe_wl(mask); 
    qe_val = qe_val(mask);
    [qe_wl, order] = sort(qe_wl(:)); 
    qe_val = qe_val(order);
    if max(qe_val) > 1.5, qe_val = qe_val / 100; end   % % → fraction
    qe_val = max(0, min(1, qe_val));                   % clamp to [0,1]
end

function cube_out = mnf_denoise_cube_adaptiveK(cube_in, retain_frac, minK)
% MNF with Adaptive K: keep enough modes to reach retain_frac of cumulative λ.
% Reconstruct using ONLY the top-K MNF modes (hard truncation).
    [x, y, b] = size(cube_in);
    N = x * y;

    % Reshape and center
    D = reshape(double(cube_in), N, b);
    mu = mean(D, 1);
    X  = D - mu;

    % Noise estimate via spatial differences (horizontal & vertical)
    I = reshape(X, x, y, b);
    Ih = []; Iv = [];
    if y > 1, Ih = reshape(I(:,1:end-1,:) - I(:,2:end,:), [], b); end
    if x > 1, Iv = reshape(I(1:end-1,:,:) - I(2:end,:,:), [], b); end
    if isempty(Ih) && isempty(Iv)
        Nmat = randn(N, b) * 1e-6;
    else
        Nmat = [Ih; Iv];
    end

    % Covariances
    Sx = (X.' * X)   / max(N,1);
    Sn = (Nmat.'*Nmat)/ max(size(Nmat,1),1);

    % Regularize Sn slightly for numerical stability
    eps_reg = 1e-6 * trace(Sn) / max(b,1);
    Sn = Sn + eps_reg * eye(b);

    % Generalized eigen: Sx v = λ Sn v
    [V, Lambda] = eig(Sx, Sn);
    lam = real(diag(Lambda));
    [lam, idx] = sort(lam, 'descend');
    V = real(V(:, idx));

    % --- Adaptive K selection ---
    lam_pos = max(lam, 0);
    if sum(lam_pos) <= 0
        K = b;
    else
        cum = cumsum(lam_pos) / sum(lam_pos);
        K = find(cum >= retain_frac, 1, 'first');
        if isempty(K), K = b; end
    end
    K = max(minK, K); K = min(K, b);

    % Transform to MNF space
    Y = X * V;

    % Hard truncate to K modes
    Yd = zeros(size(Y)); Yd(:,1:K) = Y(:,1:K);

    % Reconstruct in spectral domain (stable right-division)
    Xhat = Yd / V;

    % Re-add mean and reshape
    Dhat = Xhat + mu;
    cube_out = reshape(Dhat, x, y, b);
end
