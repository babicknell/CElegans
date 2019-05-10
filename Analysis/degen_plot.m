function [strains,D,P,T] = degen_plot(strains, len_flag)

%Input: strains is a cell of strings, eg. strains = {'CZ10175','QH6084'},
%       len_flag (True/False) scales by distal/av_original lengths if the 
%       measurement exists.
%Output: D is a matrix of relative fluorescence values for each strain
%        padded with NaNs (necessary for using boxplot).
%        P, T are aresults of ranksum and t-tests.        

if strcmp(strains,'all')
    strains = {'CZ10175', 'QH6084', 'QH6095', 'QH6101', 'QH6106', 'QH6108',...
        'QH6162', 'QH6166', 'QH6200', 'QH6314', 'QH6338', 'QH6342', 'QH6367',...
        'QH6396', 'QH6607'};
end

D_max = 300;
N = length(strains);
Mean = zeros(1, N);
STD = zeros(1, N);
D = nan(D_max ,N);

for k=1:length(strains)
    S = load(['../Results/', strains{k}, '.mat']);
    S = S.S;
    deg = [S.Degen]';
    distal_length = S(k).Length_distal;
    if len_flag && isfield(S, 'Length_orig')
        av_orig_length = mean(S(k).Length_orig);
        prox_length = 100 + S(k).Cut_site;
        deg_index = 1 - deg*(distal_length - prox_length)/(av_orig_length - prox_length);
    else
        deg_index = 1 - deg;
    end
    rec_ind = find([S.Reconnected]);
    deg_index(rec_ind)=[]; % ignore reconnected
    D(1:length(deg_index), k) = deg_index;
    Mean(k) = mean(deg_index);
    STD(k) = std(deg_index);
end

boxplot(D, 'labels', strains, 'labelorientation', 'inline')
ylabel('degeneration index', 'Interpreter', 'Latex', 'FontSize', 20)

set(gca, 'FontSize', 20, 'ylim', [0, 1])
hold on
plot(1:length(strains), Mean, 'ks')
hold off

for k=1:N
    P(k) = ranksum(D(:,1), D(:,k));
    [~, p] = ttest2(D(:,1), D(:,k));
    T(k) = p;
end

