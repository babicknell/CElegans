function [V, B] = get_fluorescence(X_d, X_b, I)
% Extract pixel values for background and degenerating traces.
% Inputs
% 'X_d', X_b': coordiantes for degen and background 
% 'I': image
% Outputs 
% 'V': degen pixel values (filtered)
% 'B': background pixel values (filtered)

% Parameters
filt_rad_b = 2; % Pixel radius for filtering
filt_rad_d = 2;
dims = size(I);

% Extract pixel values along degenerating axon
V = zeros(1, size(X_d, 2));    
for j = 1:length(V)
    Vpix = [round(X_d(2, j)), round(X_d(1, j))];
    Vav = I(Vpix(1)-filt_rad_d:Vpix(1)+filt_rad_d, Vpix(2)-filt_rad_d:Vpix(2)+filt_rad_d);    
    V(j) = mean(Vav(:));
end

% Clip background curve to allow smoothing
X_b(:, X_b(2, :)>(dims(1)-filt_rad_b)) = [];
X_b(:, X_b(1, :)>(dims(2)-filt_rad_b)) = [];

% Extract local background values
B = zeros(1, size(X_b, 2));
for l = 1:length(B)
    Bpix = [round(X_b(2, l)), round(X_b(1, l))];
    if Bpix(1)<(filt_rad_b+1) || Bpix(2)<(filt_rad_b+1) ||...
       Bpix(1)>(dims(1)-filt_rad_b-1) || Bpix(2)>(dims(2)-filt_rad_b-1)
       Bav = I(max(Bpix(1),1), max(Bpix(2),1));
    else
       Bav = I(Bpix(1)-filt_rad_b:Bpix(1)+filt_rad_b, Bpix(2)-filt_rad_b:Bpix(2)+...
               filt_rad_b);      
    end
    B(l) = mean(Bav(:));
end
pad = length(V)-length(B);
if pad
    B = [B, ones(1, pad)*B(end)];
end
