% Volumetric Fit Test
% Demonstrate Matlab fitting of parametric volume through volumetric data
% Daniel Cox
%
% Note: Requires utilities repository to be in path
% Note: Load the 3D data that you want to fit
% close all; clc; clear
% load('/home/daniel/Downloads/200310_OoC_collagenimaging/stack_3D.mat')


%% Create basis polynomial functions
% Create coordinate arrays for x,y,z
[nx, ny, nz] = size(stack_3D);
x = linspace(-1, 1, nx)';
y = linspace(-1, 1, ny)';
z = linspace(-1, 1, nz)';

% Create 2D arrays containing 1, x, x�, x�..., 1, y, y�, y�..., 1, z, z�, z�...
npowersx = 7;                                    % Polynomial powers (including 0)
npowersy = 7;                                    % Polynomial powers (including 0)
npowersz = 7;                                    % Polynomial powers (including 0)

npowers = npowersx*npowersy*npowersz;

% powers = 0:(npowers-1);                         % Array containing all powers
xpowers = x.^(0:(npowersx-1));
ypowers = y.^(0:(npowersy-1));
zpowers = z.^(0:(npowersz-1));
xyzpow = zeros(nx, ny, nz, npowers);     % Initialize basis

% Loop over all powers in x,y,z
m = 1;
for xpow = xpowers
    for ypow = ypowers
        for zpow = zpowers
            % Add 3D polynomial to set of basis functions
            xyzpow(:, :, :, m) = xpow .* ypow' .* permute(zpow, [3 2 1]);
            m = m+1;
        end
    end
end

xyzpowlin = reshape(xyzpow, nx*ny*nz, npowers);    % Reshape 4D basis set to matrix

% %% Create test data
% % Generate random coefficients for creating test data
% a = randn(npowers.^3, 1) ./ linspace(1, 30, npowers.^3)';
% noisefactor = 0.15;
% 
% % Create test data from coefficients
% Plin = xyzpowlin * a;                           % Ground truth, linear array
% Pmlin = Plin + noisefactor * randn(n^3,1);      % 'Measured' = ground truth + noise, linear array
% P = reshape(Plin, n, n, n);                     % Ground truth
% Pm = reshape(Pmlin, n, n, n);                   % 'Measured' = ground truth + noise

Pm = stack_3D;
clipvalue = 4000;
Pm(Pm>clipvalue) = clipvalue;     % Clip values above
Pmlin = Pm(:);

%% Plot test data
ramp = linspace(0,1,256)';
transparency = ((ramp-0.4).*(ramp>0.4))*0.15;
figure; plot(transparency); title('Transparency')
% figure; ims(Pm, 0.1);
figure; volshow(Pm, 'Alphamap', transparency)

%% Fit data
cf = xyzpowlin \ Pmlin;                         % Compute coefficients
Pfitlin = xyzpowlin * cf;                       % Compute fit
Pfit = reshape(Pfitlin, nx, ny, nz);               % Reshape to volumetric array

%% Compare fit to data in plot
% Pcompare = cat(2, Pm, Pfit);
Pmthresh = Pmmedfilt;
Pmthresh(Pmmedfilt > 3500) = 1.3;
Pmthresh(Pmmedfilt <= 3500) = 1.4;
Pcompare = cat(2, permute(Pm./max(Pm(:)), [2 3 1]), permute(Pmthresh, [2 3 1]));
figure; ims(Pcompare, 0.1)
figure; volshow(Pfit, 'Alphamap', transparency)
