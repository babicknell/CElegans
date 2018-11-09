function [strains,D,P,T] = degen_plot(strains)

%Input: strains is a cell of strings, eg. strains = {'CZ10175','QH6084'}
%Output: D is a matrix of relative fluorescence values for each strain
%        padded with NaNs (necessary for using boxplot).
%        P, T are aresults of ranksum and t-tests.        

if strcmp(strains,'all')
    strains = {'CZ10175', 'QH6084', 'QH6095', 'QH6101', 'QH6106', 'QH6108',...
        'QH6162', 'QH6166', 'QH6200', 'QH6314', 'QH6338', 'QH6342', 'QH6367',...
        'QH6396', 'QH6607'};
end

D_max = 300;
N=length(strains);
Mean = zeros(1, N);
STD = zeros(1, N);
D = nan(D_max ,N);

for k=1:length(strains)
    S = load(['../Results/', strains{k}, '.mat']);
    S = S.S;
    deg = [S.Degen]';
    rec_ind = find([S.Reconnected]);
    deg(rec_ind)=[]; % ignore reconnected
    D(1:length(deg), k) = deg;
    Mean(k) = mean(deg);
    STD(k) = std(deg);
end

boxplot(D, 'labels', strains, 'labelorientation', 'inline')
ylabel('relative fluorescence', 'Interpreter', 'Latex', 'FontSize', 20)
set(gca, 'FontSize', 20)
hold on
plot(1:length(strains), Mean, 'ks')
hold off

for k=1:N
    P(k) = ranksum(D(:,1), D(:,k));
    [~, p] = ttest2(D(:,1), D(:,k));
    T(k) = p;
end

