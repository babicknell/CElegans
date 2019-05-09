function [X_d, X_r] = get_coords(text_file)
% Extract coordinates of traced axons from saved NeuronJ .txt file
% Input: absolute path to tracing .txt file
% Outputs: X_d, X_r are matrices of xy coordinate pairs for 
% degenerating and regenerateing traces respectively. 

fileID = fopen(text_file);
coords=textscan(fileID, '%s %s');
N = find(strcmp(coords{1}, 'Tracing'));
if length(N) < 2 || length(N) > 3
    error('Tracing error.')
elseif length(N) == 2
    N = [N; length(coords{1}) + 1];
end

X_d = [str2double(coords{1}(2:N(2)-1)), str2double(coords{2}(2:N(2)-1))]';
X_r = [str2double(coords{1}(N(2)+1:N(3)-1)), ... 
       str2double(coords{2}(N(2)+1:N(3)-1))]';
