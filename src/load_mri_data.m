function kspace = load_mri_data(file_idx, slice_idx, multicoil)
% load_mri_data  Load a k-space slice from the fastMRI knee dataset.
%                Data source: https://fastmri.med.nyu.edu/
%
%   kspace = load_mri_data(file_idx, slice_idx, multicoil)
%
%   file_idx  - which .h5 scan to load (default: 1)
%   slice_idx - which slice to extract (default: middle slice)
%   multicoil - true to return all coils (default: false)
%
%   Returns:
%     kspace  - single-coil: (width x height)
%               multi-coil:  (width x height x num_coils)

if nargin < 1, file_idx = 1; end
if nargin < 3, multicoil = false; end

raw_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'raw');

% Extract tar if not yet done
tar_files = dir(fullfile(raw_dir, '*.tar'));
if ~isempty(tar_files) && isempty(dir(fullfile(raw_dir, '**', '*.h5')))
    fprintf('Extracting %s...\n', tar_files(1).name);
    untar(fullfile(raw_dir, tar_files(1).name), raw_dir);
end

% Locate .h5 files
h5_files = dir(fullfile(raw_dir, '**', '*.h5'));
if isempty(h5_files)
    error('No .h5 files found in %s. Check that the tar extracted correctly.', raw_dir);
end

h5_path = fullfile(h5_files(file_idx).folder, h5_files(file_idx).name);
fprintf('Loading: %s\n', h5_files(file_idx).name);

% Read k-space
% fastMRI stores kspace as complex64. Python shape: (slices, coils, height, width).
% MATLAB h5read reverses dims:                      (width, height, coils, slices).
% For single-coil knee: (368, 640, 1, num_slices).
raw = h5read(h5_path, '/kspace');

% Handle complex layout for older MATLAB versions
if isstruct(raw)
    kspace_all = raw.r + 1j * raw.i;
else
    kspace_all = raw;
end

% kspace_all: (width, height, coils, slices)
num_slices = size(kspace_all, 4);
if nargin < 2, slice_idx = round(num_slices / 2); end

kspace = squeeze(kspace_all(:, :, :, slice_idx));  % (width, height, coils)

if ~multicoil
    kspace = kspace(:, :, 1);  % (width, height)
end

fprintf('k-space ready: %dx%d  coils: %d  (slice %d/%d)\n', ...
    size(kspace,1), size(kspace,2), size(kspace_all,3), slice_idx, num_slices);
