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
n = length(fnames);
name{n} = [];
len{n} = [];
len_d{n} = [];
cut_site{n} = [];
tort{n} = [];
fDist{n} = [];
mDist{n} = [];
fAng{n} = [];
mAng{n} = [];
degen{n} = [];
recon{n} = [];
regrew{n} = [];
fused{n} = [];
xy_coords{n} = [];

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
    I = imread([imageDir, '/', fnames{k}(1:end-4),imtype]);
    dims = size(I);
    file = [imageDir, '/', fnames{k}];
    dorsal_flag = strcmp(fnames{k}(end-4), 'D');
    [coords_d, coords_r] = get_coords(file);

    % Parameterise traces and extract pixel values 
    [X_d, X_r, X_b] = parameterise(coords_d, coords_r);
    [X_d_T, ~, Z_r, d_ind] = transform(X_d, X_r, dims, sf, dorsal_flag);
    [V, B] = get_fluorescence(X_d, X_r, X_b, I);
    
    % measurements
    distal_s = d_ind + round(100*sf); % average from 100 um distal to cut 
    V_sub = V - B;  % background subtraction
    V_sub(V_sub<0) = 0;
    V_rel = V_sub./mean(V_sub(1:d_ind)); %relative to initial segment
    V_rel_av = mean(V_rel(distal_s:end));
    
    r_dist = abs(Z_r(1, :));
    r_dist(Z_r(2, :) < 0) = sqrt(sum(X_d_T(:, Z_r(2, :) < 0).^2)); 
    
    name{k} = fnames{k}(1:end-4);
    len{k} = size(X_r, 2)/sf;
    len_d{k} = size(X_d, 2)/sf;
    cut_site{k} = d_ind/sf;
    tort{k} = len{k}/(sqrt((X_r(1, end) - X_r(1, 1)).^2 + (X_r(2, end) - X_r(2, 1)).^2)/sf);
    fDist{k} = r_dist(end);
    mDist{k} = mean(r_dist);
    l_thr = 20;
    if size(Z_r, 2) > l_thr
        fAng{k} = atan2d(-Z_r(1, end), Z_r(2, end));
        mAng{k} = mean(abs(atan2d(-Z_r(1, l_thr:end), Z_r(2, l_thr:end))));
    else
        fAng{k}=nan;
        mAng{k}=nan;
        tort{k}=nan;
    end
    degen{k} = V_rel_av;
    xy_coords{k} = Z_r;

    % Plot
    figure(1);
    plot([0, 0], [0, 250],'b--','LineWidth',2);
    hold on
    plot([0, 0], [-50, 0],'b-','LineWidth',2)
    plot(Z_r(1, :), Z_r(2, :), 'r', 'LineWidth', 0.5)
    plot(Z_r(1, end), Z_r(2, end), 'ko', 'MarkerSize', 4)
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
S = struct('Name', name, 'Coords', xy_coords, 'Length', len, 'Length_distal', len_d, 'Cut_site', cut_site, 'Tortuosity', tort, ...
    'Final_dist', fDist, 'Mean_dist', mDist, 'Final_ang', fAng, 'Mean_ang', mAng,...
    'Degen', degen, 'Regrew', regrew, 'Reconnected', recon, 'Fused', fused);
save(['../Results/', parts{2}{1}], 'S')

