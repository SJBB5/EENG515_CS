%% MRI Reconstruction Using Haar Wavelets and ISTA 
% Sophia Mimlitz, 11/18/26
% This function expects two inputs from the sampling functions: K and mask

function x = reconstruct_wavelet_ista(y, mask, lambda, nIter)
% y:        complex, centered k-space, zero-filled at missing locations
% mask:     logical array, true where a measurement was kept
% lambda:   wavelet threshold with the unit ISTA step used here

% Check for errors in function call
assert(ismatrix(y) && isequal(size(y), size(mask)), ...
    'y and mask must be 2D arrays of the same size.');
assert(all(mod(size(y), 2) == 0), ...
    'This requires even dimensions.');

% Rearrange and scale
scale = sqrt(numel(y));
F  = @(z) fftshift(fft2(ifftshift(z))) / scale;
FH = @(z) fftshift(ifft2(ifftshift(z))) * scale;

mask = logical(mask);
y = mask .* y;
x = FH(y);  % Start from the zero-filled reconstruction

% Shrink the magnitude of each coefficient while preserving its phase.
shrink = @(c) c .* max(0, 1 - lambda ./ max(abs(c), eps));

for k = 1:nIter
    % 1. Compare predicted and actual data at measured locations.
    residual = mask .* (F(x) - y);

    % 2. Bring that error back into the image and correct the estimate.
    z = x - FH(residual);

    % 3. Split each 2-by-2 block into one average and three details.
    a = z(1:2:end, 1:2:end);
    b = z(1:2:end, 2:2:end);
    c = z(2:2:end, 1:2:end);
    d = z(2:2:end, 2:2:end);

    avg = (a + b + c + d) / 2;
    h = shrink((a - b + c - d) / 2);
    v = shrink((a + b - c - d) / 2);
    q = shrink((a - b - c + d) / 2);

    % 4. Turn the average and updated details back into pixels.
    x(1:2:end, 1:2:end) = (avg + h + v + q) / 2;
    x(1:2:end, 2:2:end) = (avg - h + v - q) / 2;
    x(2:2:end, 1:2:end) = (avg + h - v - q) / 2;
    x(2:2:end, 2:2:end) = (avg - h - v + q) / 2;
end
end
