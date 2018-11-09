function A = angle_dist(strains)
%Input: 'strains' is a cell of strings eg. {'CZ10175', 'QH6084'}
%Ouput: 'A' is a cell of angel measurements

if strcmp(strains,'all')
    strains = {'CZ10175', 'QH6084', 'QH6095', 'QH6101', 'QH6106', 'QH6108',...
        'QH6162', 'QH6166', 'QH6200', 'QH6314', 'QH6338', 'QH6342', 'QH6367',...
        'QH6396', 'QH6607'};
end

N = length(strains);
A{N} = [];

for j = 1:length(strains)
    S = load(['../Results/', strains{j}, '.mat']);
    S = S.S;
    %exclude non-regrow  
    reg_ind = find([S.Regrew]==0);
    rec_ind = find([S.Reconnected]);
    del = unique([reg_ind, rec_ind]); % ignore non-regrowing and reconnected
    S(del)=[];

    theta = 0:180;
    Pf = zeros(1, length(theta));
    Pavg = zeros(1, length(theta));

    fang = abs([S.Final_ang]);
    ang = [S.Mean_ang];
    ang(isnan(ang)) = [];
    fang(isnan(fang)) = [];

    for k = 1:length(Pf)
        Pf(k) = sum((fang)<theta(k))/length(fang);
        Pavg(k) = sum((ang)<theta(k))/length(ang);    
    end

    %final angle%
    plot(theta,Pf)
    xlabel('$\theta$ (degrees)', 'Interpreter', 'latex', 'FontSize',18)
    ylabel('$P(|\theta_{final}|>\theta)$', 'Interpreter', 'latex', 'FontSize',18)
    hold on

    %mean angle%
    % plot(theta,Pavg)
    % xlabel('$\theta$ (degrees)', 'Interpreter', 'latex', 'FontSize', 18)
    % ylabel('$P(|\theta|_{avg}>\theta)$', 'Interpreter', 'latex', 'FontSize', 18)
    % hold on
     
    A{j} = fang;
    %A{j} = ang;
end
legend(strains)