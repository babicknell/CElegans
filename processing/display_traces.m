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

%% Main loop
for k = ki:kf
    fprintf([fnames{k}, '\n'])
    I = imread([imageDir, '/', fnames{k}(1:end-4),imtype]);
    dims = size(I);
    file = [imageDir, '/', fnames{k}];
    dorsal_flag = strcmp(fnames{k}(end-4), 'D');
    [coords_d, coords_r] = get_coords(file);
    
    % Parameterise traces and extract pixel values 
    [X_d, X_r, X_b] = parameterise(coords_d, coords_r);
    [X_d_T, X_r_T, Z_r, d_ind] = transform(X_d, X_r, dims, sf, dorsal_flag);
    [V, B] = get_fluorescence(X_d, X_r, X_b, I);
    
    V_sub = V - B;  %background subtraction
    V_sub(V_sub<0) = 0;
    V_rel = V_sub./mean(V_sub(1:d_ind)); %relative to initial segment

    % Plotting
    figure(1)
    subplot(2, 2, 1)
    imagesc(I)
    colormap gray
    brighten(0.5)
    hold on
    plot(X_d(1, :), X_d(2, :), 'b', 'LineWidth', 1);
    plot(X_r(1, :), X_r(2, :), 'r', 'LineWidth', 1);
    plot(X_b(1, :), X_b(2, :), 'c' ,'LineWidth', 1)
    legend('Orginal', 'Regen', 'Background')
    set(gca,'Xtick', [], 'Ytick', []);
    axis equal
    hold off
    
    subplot(2,2,2)
    plot(X_d_T(1, :), X_d_T(2, :), 'b', 'LineWidth', 2);
    hold on
    plot(X_r_T(1, :), X_r_T(2, :), 'r', 'LineWidth', 2);
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
    plot(Z_r(1, :), Z_r(2, :), 'r', 'LineWidth', 2)
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

