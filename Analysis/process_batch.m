function process_batch()
% Process batch of images and tracing files, and save results.
%% Setup
sf = 3.08; %pixels per micron

fprintf('Select directory containing images to process...\n')
imageDir = uigetdir('','Select directory containing images to process...');
if ~imageDir
	disp('Cancelled.')
	return
end
% read axotomy score file
AS = dir([imageDir, '/AS*']);
fid = fopen([imageDir, '/', AS.name]);
Recon = textscan(fid, '%s %s %s %s %s %s'); 
fclose(fid);

imtype = dir([imageDir, '/*.jpg']);  %check if files are jpg (tif otheriwse) 
if isempty(imtype)
    imtype = '.tif';
else
    imtype = '.jpg';
end

list = dir([imageDir, '/2*.txt']);
fnames = {list.name};

name{length(fnames)} = [];
len{length(fnames)} = [];
len_p{length(fnames)} = [];
len_d{length(fnames)} = [];
tort{length(fnames)} = [];
fDist{length(fnames)} = [];
mDist{length(fnames)} = [];
fAng{length(fnames)} = [];
mAng{length(fnames)} = [];
degen{length(fnames)} = [];
recon{length(fnames)} = [];
regrew{length(fnames)} = [];
fused{length(fnames)} = [];
%% Process
for k=1:length(fnames)
    fprintf([fnames{k}, '\n'])
    parts = textscan(fnames{k}, '%s %s %s %s %s %s', 'Delimiter', '_');
    for ind = 2:length(Recon{2})
        if strcmp(parts{1}, Recon{1}(ind)) && strcmp(parts{3}, Recon{2}(ind))
            regrew{k} = str2double(Recon{4}(ind));
            recon{k} = str2double(Recon{5}(ind));
            fused{k} = str2double(Recon{6}(ind));
            break
        end
    end
    I = imread([imageDir, '/', fnames{k}(1:end-4), imtype]);
    dims = size(I);    
    fileID = fopen([imageDir, '/', fnames{k}]);
    coords=textscan(fileID, '%s %s');
    fclose(fileID);
    N = find(strcmp(coords{1}, 'Tracing'));
    if length(N) < 2 || length(N) > 3
       fprintf('Tracing error.\n')
    elseif length(N)==2
        N = [N, length(coords{1})+1];
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
    V_rel_av = mean(V_rel(distal_s:distal_f));
  
    name{k} = fnames{k}(1:end-4);
    len{k} = length(Xr)/sf;
    len_p{k} = y_hat(end)/sf;
    len_d{k} = length(Xd)/sf;
    tort{k} = len{k}/(sqrt((Xr(1, end) - Xr(1, 1)).^2 + (Xr(2, end) - Xr(2, 1)).^2)/sf);
    fDist{k} = x_min(end)/sf;
    mDist{k} = mean(abs((x_min/sf)));
    l_thr = 20;
    if length(x_hat) > l_thr
        fAng{k} = atan2d(-x_hat(end), y_hat(end));
        mAng{k} = mean(abs(atan2d(-x_hat(l_thr:end), y_hat(l_thr:end))));
    else
        fAng{k}=nan;
        mAng{k}=nan;
        tort{k}=nan;
    end
    degen{k} = V_rel_av;

    % Plot
    figure(1);
    plot([0, 0], [0, 250],'b--','LineWidth',2);
    hold on
    plot([0, 0], [-50, 0],'b-','LineWidth',2)
    plot(x_hat/sf, y_hat/sf, 'r', 'LineWidth', 0.5)
    plot(x_hat(end)/sf, y_hat(end)/sf, 'ko', 'MarkerSize', 4)
    axis equal
    axis([-50, 50, -50, 200])
    xlabel('displacement ($\mu{m}$)', 'Interpreter','Latex', 'FontSize',18)
    ylabel({'original axon arclength ($\mu{m}$)'}, 'Interpreter', 'Latex', 'FontSize',18)
    title(parts{2})
    set(gca,'FontSize',20)
end
    
figname = ['../Results/',parts{2}{1},'.pdf'];
figure(1);
print('-dpdf',figname)
close(figure(1));
S = struct('Name', name, 'Length', len, 'Length_proj', len_p, 'Length_distal', len_d, 'Tortuosity', tort, ...
    'Final_dist', fDist, 'Mean_dist', mDist, 'Final_ang', fAng, 'Mean_ang', mAng,...
    'Degen', degen, 'Regrew', regrew, 'Reconnected', recon, 'Fused', fused);
save(['../Results/', parts{2}{1}], 'S')


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


