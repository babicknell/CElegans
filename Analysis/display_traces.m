function display_traces(filenum)
% Display traces overlaid on images, and plot background-subtracted 
% fluorsecence along degenerating axon. 
% Input: optional integer argument 'filenum' to display specified image 
% (otherwise loops through all images in directory). 
%% Setup
sf = 3.08; % Pixels per micron scale factor

fprintf('Select directory containing images to process...\n')
imageDir = uigetdir('','Select directory containing images to process...');
if ~imageDir
    disp('Cancelled.')
    return
end
list = dir([imageDir,'/2*.txt']);
fnames = {list.name};

if nargin == 0
    ki = 1;
    kf = length(fnames);
else
    ki = filenum;
    kf = filenum;
end

imtype = dir([imageDir,'/*.jpg']);  %check if files are jpg (tif otheriwse) 
if isempty(imtype)
    imtype = '.tif';
else
    imtype = '.jpg';
end

%% Parameterise traces and extract pixel values 
for k = ki:kf
    fprintf([fnames{k}, '\n'])
    I = imread([imageDir, '/', fnames{k}(1:end-4),imtype]);
    dims = size(I);
    fileID = fopen([imageDir, '/', fnames{k}]);
    coords=textscan(fileID, '%s %s');
    N = find(strcmp(coords{1}, 'Tracing'));
    if length(N) < 2 || length(N) > 3
        fprintf('Tracing error.\n')
        continue
    elseif length(N) == 2
        N = [N, length(coords{1}+1)];
    end
    % Arclength parmeterisation of degenerating and regenerating axons.
    [xd, yd] = arclength(str2double(coords{1}(2:N(2)-1)),...
               str2double(coords{2}(2:N(2)-1)));
    [xr, yr] = arclength(str2double(coords{1}(N(2)+1:N(3)-1)),...
               str2double(coords{2}(N(2)+1:N(3)-1)));
    % Create a curve alongside degen axon for local background subtraction.
    xb = xd';
    yb = yd';
    filtN = 100;
    xb = smooth(xb,floor(filtN/2), 'moving');
    yb = smooth(yb,floor(filtN/2), 'moving');
    nor = [-diff(yb), diff(xb)];
    nor = nor./repmat(sqrt(sum(nor.^2,2)), 1, 2);
    xb(end) = [];
    yb(end) = [];
    xb = xb + 10*nor(:, 1); % Shift 10 pixels normal to (x_d, y_d)
    yb = yb + 10*nor(:, 2);
    
    f_radb = 2;            % Clip curve to allow smoothing
    xb(yb>(dims(1)-f_radb)) = [];
    yb(yb>(dims(1)-f_radb)) = [];
    yb(xb>(dims(2)-f_radb)) = [];
    xb(xb>(dims(2)-f_radb)) = [];
    
    % Extract pixel values long degenerating axon
    V = zeros(1, length(xd));
    f_rad = 2;             % Pixel radius for smoothing
    for j = 1:length(xd)
        Vpix = [round(yd(j)), round(xd(j))];
        Vav = I(Vpix(1)-f_rad:Vpix(1)+f_rad, Vpix(2)-f_rad:Vpix(2)+f_rad);    
        V(j) = mean(Vav(:));
    end
    
    % Extract local background values
    B = zeros(1, length(xb));
    for l = 1:length(xb)
        Bpix = [round(yb(l)), round(xb(l))];
        if Bpix(1)<(f_radb+1) || Bpix(2)<(f_radb+1) ||...
           Bpix(1)>(dims(1)-f_radb-1) || Bpix(2)>(dims(2)-f_radb-1)
           Bav = I(max(Bpix(1),1), max(Bpix(2),1));
        else
           Bav = I(Bpix(1)-f_radb:Bpix(1)+f_radb, Bpix(2)-f_radb:Bpix(2)+...
                   f_radb);      
        end
        B(l) = mean(Bav(:));
    end
    pad = length(V)-length(B);
    if pad
        B = [B, ones(1, pad)*B(end)];
    end
    
    % Centre and rotate traces
    xr_c = xr-xd(1);
    yr_c = dims(1) - yr;
    yr_c = yr_c -(dims(1) - yd(1));
    Xr = [xr_c; yr_c];
    xd_c = xd - xd(1);
    yd_c = dims(1) - yd;
    yd_c = yd_c - yd_c(1);
    Xd = [xd_c; yd_c];
    
    [~,d_ind] = mindist(Xr(:,1), Xd, 0); % index of cut site
    init_Ang = atan2(yd_c(d_ind+5)-yd_c(d_ind), xd_c(d_ind+5)-xd_c(d_ind));
    Rot = [cos(-init_Ang+pi/2), -sin(-init_Ang+pi/2);
        sin(-init_Ang+pi/2), cos(-init_Ang+pi/2)];
    Xd = Rot*Xd;
    Xr = Rot*Xr;

    % Parameterise regen with respect to signed distance from degen.
    x_min = zeros(1, size(Xr,2));      % min distance to distal degen axon
    x_min_coord = zeros(1, size(Xr,2));% assocaited arclength coordinate
    x_hat = zeros(1, size(Xr,2));      % min distance to entire orginal axon 
    y_hat = zeros(1, size(Xr,2));      % assocaited arclength coordinate
    for m = 1:length(x_min)
        [x_min(m), x_min_coord(m)] = mindist(Xr(:,m), Xd(:,d_ind:end), 0);
        [x_hat(m), y_hat(m)] = mindist(Xr(:,m), Xd, 1);
    end
    y_hat = y_hat - d_ind;
    if strcmp(fnames{k}(end-4), 'D')
        x_hat = -x_hat; % flip x-coordinate if dorsal orientation 
    end
    
    distal_s = max(y_hat) + d_ind; % first index of degen distal to regen tip 
    % to measure degen over 100 microns:
    % distal_f = min(max(y_hat)+d_ind+round(sf*100),length(V));
    % to meaure over whole distal section:
    distal_f = length(V);
    V_sub = V - B;  %background subtraction
    V_sub(V_sub<0) = 0;
    V_rel = V_sub./mean(V_sub(1:d_ind)); %relative to initial segment
    % to print average relative fluorescence:
    % mean(Vrel(distal_s:distal_f))
    
    %% Plotting %%
    figure(1)
    subplot(2, 2, 1)
    imagesc(I)
    colormap gray
    brighten(0.5)
    hold on
    plot(xd, yd, 'b', 'LineWidth', 1);
    plot(xr, yr, 'r', 'LineWidth', 1);
    plot(xb, yb, 'c' ,'LineWidth', 1)
    legend('Orginal', 'Regen', 'Background')
    set(gca,'Xtick', [], 'Ytick', []);
    axis equal
    hold off
    
    subplot(2,2,2)
    plot(Xd(1,:)/sf, Xd(2,:)/sf, 'b', 'LineWidth', 2);
    hold on
    plot(Xr(1,:)/sf, Xr(2,:)/sf, 'r', 'LineWidth', 2);
    axis equal
    axis([-150, 150, 0, 400])
    hold off
    title('Rotated')
    xlabel('$X$ ($\mu{m}$)', 'Interpreter', 'Latex','FontSize', 20)
    ylabel('$Y$ ($\mu{m}$)', 'Interpreter', 'Latex','FontSize', 20)
    
    subplot(2,2,3)
    plot([0,0], [0,250], 'b--', 'LineWidth', 2)
    hold on
    plot([0,0], [-50,0], 'b-', 'LineWidth', 2)
    plot(x_hat/sf, y_hat/sf, 'r', 'LineWidth', 2)
    axis equal
    axis([-50, 50, -50, 150])
    xlabel('$\hat{X}$ ($\mu{m}$)', 'Interpreter', 'Latex','FontSize', 20)
    ylabel('$\hat{Y}$ ($\mu{m}$)', 'Interpreter', 'Latex','FontSize', 20)
    title('Transformed')
    hold off

    subplot(2,2,4)
    plot((1:length(V_rel))/sf, V_rel)
    hold on
    plot([d_ind, d_ind]/sf, [0,1.5], 'r--')
    axis([0, length(V_rel)/sf, 0, 1.5])
    title('Relative Fluorescence')
    hold off
    
    % To view raw, background and background subtracted fluoresence values
%     figure(2)
%     plot((1:length(V))/sf, V, 'k')
%     hold on
%     plot((1:length(B))/sf, B, 'c')
%     plot((1:length(B))/sf, V_sub, 'b')
%     plot([d_ind,d_ind]/sf, [0,100], 'r--')
%     legend('Raw', 'Backgound', 'Subtracted', 'Cut')
%     hold off
    
    pause
end

%% Sub-functions
function [xi, yi, ti] = arclength(x,y)
% Parameterise coordinates by arclength.
[~,ii] = unique([x,y], 'rows', 'stable');
x = x(ii); y = y(ii);
t = cumsum(sqrt([0, diff(x')].^2 + [0, diff(y')].^2));
ti = linspace(t(1), t(end), round(t(end)));
xi = interp1(t, x', ti, 'pchip');
yi = interp1(t, y', ti, 'pchip');
if isempty(xi)
    xi = x;
    yi = y;
end

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



