function [X_d_T, X_r_T, Z_r, d_ind] = transform(X_d, X_r, image_size, sf, dorsal)
% Centre and rotate coordinates, and express X_r wrt distance from X_d.
% Inputs 
% 'X_d', X_r': coordinate pairs from 'parameterise.m', 
% 'image_size': image dimensions
% 'dorsal': flag for dorsal orientation (True/False)
% Outputs 
% 'X_d', 'X_r': matrices of transformed coordinate pairs.
% 'Z_r': distance from X_r to closeest point on X_d.
% 'd_ind': index of cut site

% Centre and rotate traces
xr_c = X_r(1, :) - X_d(1, 1);
yr_c = X_d(2, 1) - X_r(2, :);
X_r = [xr_c; yr_c];
xd_c = X_d(1, :) - X_d(1, 1);
yd_c = image_size(1) - X_d(2, :);
yd_c = yd_c - yd_c(1);
X_d = [xd_c; yd_c];

[~, d_ind] = mindist(X_r(:, 1), X_d, 0); % index of cut site
init_Ang = atan2(yd_c(d_ind+5)-yd_c(d_ind), xd_c(d_ind+5)-xd_c(d_ind));
Rot = [cos(-init_Ang+pi/2), -sin(-init_Ang+pi/2);
    sin(-init_Ang+pi/2), cos(-init_Ang+pi/2)];
X_d = Rot*X_d;
X_r = Rot*X_r;

% Parameterise regen with respect to displacement from degen.
xr_min = zeros(1, size(X_r, 2));         % min distance to distal degen axon
xr_min_coord = zeros(1, size(X_r, 2));   % assocaited arclength coordinate
x_hat = zeros(1, size(X_r, 2));         % min distance to entire orginal axon 
y_hat = zeros(1, size(X_r, 2));         % assocaited arclength coordinate
for m = 1:length(xr_min)
    [xr_min(m), xr_min_coord(m)] = mindist(X_r(:, m), X_d(:, d_ind:end), 0);
    [x_hat(m), y_hat(m)] = mindist(X_r(:, m), X_d, 1);
end
y_hat = y_hat - d_ind;
if dorsal == 1
    x_hat = -x_hat; % flip x-coordinate if dorsal orientation 
end

X_d_T = X_d/sf;
X_r_T = X_r/sf;
Z_r = [x_hat; y_hat]/sf;

%% Sub-functions
function [md, d_ind] = mindist(p, r, signed)
% Find minimum distance between point p and curve r.
d = zeros(1, length(r));
for k=1:length(r)
    d(k) = sqrt((p(1) - r(1, k)).^2 + (p(2) - r(2, k)).^2);
end
[~, d_ind] = min(abs(d));
md = d(d_ind);
if signed
   T = r(:, d_ind+1) - r(:, d_ind);
   nor = [-T(2); T(1)];
   v = p - r(:, d_ind);
   md = -sign(v'*nor)*(md);
end
