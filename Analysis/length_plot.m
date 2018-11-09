function L = length_plot(strains)

%Input: strains is a cell of strings, eg. strains = {'CZ10175','QH6084'}
%Output: L is a matrix of lengths for each strain
%        padded with NaNs (necessary for using boxplot). 

if strcmp(strains,'all')
    strains = {'CZ10175', 'QH6084', 'QH6095', 'QH6101', 'QH6106', 'QH6108',...
        'QH6162', 'QH6166', 'QH6200', 'QH6314', 'QH6338', 'QH6342', 'QH6367',...
        'QH6396', 'QH6607'};
end

L_max = 300;
N=length(strains);
Mean = zeros(1 ,N);
STD = zeros(1, N);
L = nan(L_max, N);

for k=1:length(strains)
    S = load(['../Results/', strains{k}, '.mat']);
    S = S.S;
    len = [S.Length]';
    reg_ind = find([S.Regrew]==0);
    fuse_ind = find([S.Fused]);
    del = unique([reg_ind, fuse_ind]); % ignore non-regrowing or fused
    len(del)=[];
    L(1:length(len),k) = len;
    Mean(k) = mean(len);
    STD(k) = std(len);
end

boxplot(L, 'labels', strains, 'labelorientation', 'inline')
ylabel('length ($\mu{m}$)', 'Interpreter', 'Latex', 'FontSize', 20)
set(gca,'FontSize',20)
hold on
plot(1:length(strains), Mean, 'ks')
hold off
