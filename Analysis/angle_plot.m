function A = angle_plot(strains)

%Input: strains is a cell of strings, eg. strains = {'CZ10175','QH6084'}
%Output: A is a matrix of final angle values for each strain padded with 
%        NaNs (necessary for using boxplot). 

if strcmp(strains,'all')
    strains = {'CZ10175', 'QH6084', 'QH6095', 'QH6101', 'QH6106', 'QH6108',...
        'QH6162', 'QH6166', 'QH6200', 'QH6314', 'QH6338', 'QH6342', 'QH6367',...
        'QH6396', 'QH6607'};
end

A_max = 300;
N = length(strains);
A = nan(A_max,N);
for k=1:length(strains)
    S = load(['../Results/', strains{k}, '.mat']);
    S = S.S;
    ang = [S.Final_ang]';
    reg_ind = find([S.Regrew]==0);  % ignore non-regrowing
%     r_ind = find([S.Reconnected]);
    ang(reg_ind)=[];
    A(1:length(ang),k) = abs(ang);    
end

boxplot(A, 'labels', strains, 'labelorientation', 'inline')
ylabel('final angle (degrees)', 'Interpreter', 'Latex', 'FontSize', 20)
set(gca, 'FontSize', 20)
