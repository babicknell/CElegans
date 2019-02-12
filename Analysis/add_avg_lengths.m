function add_avg_lengths(strains, vals)
% Add measurementns of original axon length to results data structures.
%Input: 'strains' is a cell of strings eg. {'CZ10175', 'QH6084'}
%       'vals' is a cell of length vectors for each strain eg. {[100, 45,
%       90], [100,100,123]}

if strcmp(strains,'all')
    strains = {'CZ10175', 'QH6084', 'QH6095', 'QH6101', 'QH6106', 'QH6108',...
        'QH6162', 'QH6166', 'QH6200', 'QH6314', 'QH6338', 'QH6342', 'QH6367',...
        'QH6396', 'QH6607'};
end

for j = 1:length(strains)
    S = load(['../Results/', strains{j}, '.mat']);
    S = S.S;
    for k = 1:length(S)
        S(k).Length_orig = vals(j);
    end
    save(['../Results/', strains{j}, '.mat'], 'S')
    clear S
end
