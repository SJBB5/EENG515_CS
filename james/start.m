function k = fft2c(x)
    x = ifftshift(x); % move image center (low frequency) to [1,1]
    k = fftshift(fft2(x)) / sqrt(numel(x)); % take 2D DFT and normalize
end

function x = ifft2c(k)
    k = ifftshift(k);
    x = fftshift(ifft2(k)) * sqrt(numel(k)); % inverse 2D DFT and normalize
end

function [img,k_full,scale] = prep_slice(k_raw, crop_size)
    img_full = ifft2c(k_raw); % get from k space to image space
    [num_rows, num_cols] = size(img_full); % find image size
    % crop image
    rowStart = floor((num_rows - crop_size) / 2) + 1;
    colStart = floor((num_cols - crop_size) / 2) + 1;
    img = img_full(rowStart:rowStart + crop_size - 1, colStart:colStart + crop_size - 1);

    scale = max(abs(img(:))); % find highest intensity
    img = img / scale; % normalize image intensity
    k_full = fft2c(img); % retain full-resolution k-space
end

function [mask,line_mask] = make_mask(num_rows,num_cols,accel_factor,center_frac)
    target_total = round(num_cols/accel_factor); % defines how many total columns we want
    num_center = round(num_cols * center_frac); % how many are in fully sampled center

    % define center
    center_start = floor((num_cols - num_center) / 2) + 1;
    center_end = center_start + num_center - 1;
    
    % get probabilities for each non-center column
    num_outside = num_cols - num_center;
    prob_outside = (target_total - num_center) / num_outside;
    
    % get random columns and make sure we have our center
    line_mask = rand(1, num_cols) < prob_outside;
    line_mask(center_start:center_end) = true;

    mask = repmat(line_mask, num_rows, 1); % expand column vector into full mask 
end

function [zf, y, mask, info] = zero_filled_baseline(k_full, accel, center_frac,seed)

    [num_rows,num_cols] = size(k_full); 

    if nargin >= 4 && ~isempty(seed)
        rng(seed);          % reset the random generator so the mask is repeatable
    end

    [mask,line_mask] = make_mask(num_rows,num_cols,accel,center_frac); % build our sampling mask
    y = k_full .* mask; % retain sampled k-space measurements
    zf = ifft2c(y); % reconstruct zero-filled image
    
    % reporitng
    info.num_sampled_cols = nnz(line_mask);
    info.actual_accel = num_cols / nnz(line_mask);
    info.sampled_fraction = nnz(line_mask) / num_cols;
end


% load a slice 
kspace_raw = double(load_mri_data_alt());
kspace_raw = kspace_raw.'; % get dimensions right

% find our ground truth
[img, k_full] = prep_slice(kspace_raw,320);

% undersample and reconstruct
[zf, y, mask, info] = zero_filled_baseline(k_full, 4, 0.08,0);
fprintf('Requested R = 4, realized R = %.2f (%d columns)\n',info.actual_accel, info.num_sampled_cols);

% error metrics (compare magnitudes)
a = abs(zf);   b = abs(img);
nmse    = norm(a - b, 'fro')^2 / norm(b, 'fro')^2;
psnr_db = psnr(a, b, max(b(:)));
ssim_v  = ssim(a, b, 'DynamicRange', max(b(:)));
fprintf('Zero-filled: NMSE = %.4f, PSNR = %.2f dB, SSIM = %.4f\n',nmse, psnr_db, ssim_v);

% graphs
figure;
subplot(2,2,1); imshow(abs(img), []);
title('Ground truth |x|');
subplot(2,2,2); imshow(mask, []);                     
title('Sampling mask M');
subplot(2,2,3); imshow(abs(zf), []);                  
title('Zero-filled |A^H y|');
subplot(2,2,4); imshow(abs(zf) - abs(img), []);       
title('Error (zero-filled - truth)');
