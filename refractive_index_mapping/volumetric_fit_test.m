% Volumetric Fit Test
% Demonstrate Matlab fitting of parametric volume through volumetric data
% Daniel Cox
%
% Note: Requires utilities repository to be in path
close all; clc; clear

%% Create basis polynomial functions
% Create coordinate arrays for x,y,z
n = 32;
x = linspace(-1, 1, n)';
y = linspace(-1, 1, n)';
z = linspace(-1, 1, n)';

% Create 2D arrays containing 1, x, x², x³..., 1, y, y², y³..., 1, z, z², z³...
npowers = 7;                                    % Polynomial powers (including 0)

powers = 0:(npowers-1);                         % Array containing all powers
xpowers = x.^powers;
ypowers = y.^powers;
zpowers = z.^powers;
xyzpow = zeros(n, n, n, length(powers).^3);     % Initialize basis

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

xyzpowlin = reshape(xyzpow, n^3, npowers^3);    % Reshape 4D basis set to matrix

%% Create test data
% Generate random coefficients for creating test data
a = randn(npowers.^3, 1) ./ linspace(1, 30, npowers.^3)';
noisefactor = 0.15;

% Create test data from coefficients
Plin = xyzpowlin * a;                           % Ground truth, linear array
Pmlin = Plin + noisefactor * randn(n^3,1);      % 'Measured' = ground truth + noise, linear array
P = reshape(Plin, n, n, n);                     % Ground truth
Pm = reshape(Pmlin, n, n, n);                   % 'Measured' = ground truth + noise


%% Plot test data
ramp = linspace(0,1,256)';
transparency = ((ramp-0.4).*(ramp>0.4))*0.15;
figure; plot(transparency); title('Transparency')
figure; ims(Pm, 0.1);
figure; volshow(Pm, 'Alphamap', transparency)

%% Fit data
cf = xyzpowlin \ Pmlin;                         % Compute coefficients
Pfitlin = xyzpowlin * cf;                       % Compute fit
Pfit = reshape(Pfitlin, n, n, n);               % Reshape to volumetric array

%% Compare fit to data in plot
Pcompare = cat(2, Pm, Pfit);
figure; ims(Pcompare, 0.1)
figure; volshow(Pfit, 'Alphamap', transparency)
