function [X_d, X_r, X_b] = parameterise(coords_d, coords_r)
% Parameterise coordinates by arclength and create background trace.
% Inputs: extracted xy coordsinate pairs from 'get_coords.m'.
% Outputs: X_d, X_r, X_b are matrices of coordinate pairs for 
% degenerating, regenerateing and backgound traces respectively. 

% Parameters
filt_span = 50; % moving average filter span for background trace 

% Arclength parmeterisation of degenerating and regenerating axons.
[xd, yd, ~] = arclength(coords_d(1, :), coords_d(2, :));
[xr, yr, ~] = arclength(coords_r(1, :), coords_r(2, :));

% Create a curve alongside degen axon for local background subtraction.
xb = xd;
yb = yd;
xb = smooth(xb, filt_span, 'moving');
yb = smooth(yb, filt_span, 'moving');
nor = [-diff(yb), diff(xb)];
nor = nor./repmat(sqrt(sum(nor.^2,2)), 1, 2);
xb(end) = [];
yb(end) = [];
xb = xb + 10*nor(:, 1); % Shift 10 pixels normal to (x_d, y_d)
yb = yb + 10*nor(:, 2);

X_d = [xd; yd];
X_r = [xr; yr];
X_b = [xb'; yb'];

%% Sub-functions
function [xi, yi, ti] = arclength(x,y)
% Parameterise coordinates by arclength.
[~,ii] = unique([x', y'], 'rows', 'stable');
x = x(ii); y = y(ii);
t = cumsum(sqrt([0, diff(x)].^2 + [0, diff(y)].^2));
ti = linspace(t(1), t(end), round(t(end)));
xi = interp1(t, x, ti, 'pchip');
yi = interp1(t, y, ti, 'pchip');
if isempty(xi)
    xi = x;
    yi = y;
end
