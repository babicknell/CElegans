function [V, B] = get_fluorescence(X_d, X_r, X_b, I)
% Extract pixel values for background and degenerating traces.
% Inputs
% 'X_d', X_b': coordinates for degen and background 
% 'I': image
% Outputs 
% 'V': degen pixel values (filtered)
% 'B': background pixel values (filtered)

% Parameters
filt_rad_b = 4; % Pixel radius for filtering background signal
filt_rad_d = 2; % Pixel radius for filtering degen signal
ov_thr = 5; % distance threshold to define overlapping traces 
dims = size(I);

% Extract pixel values along degenerating axon
V = zeros(1, size(X_d, 2));    
for j = 1:length(V)
    Vpix = [round(X_d(2, j)), round(X_d(1, j))];
    Vav = I(Vpix(1)-filt_rad_d:Vpix(1)+filt_rad_d, Vpix(2)-...
            filt_rad_d:Vpix(2)+filt_rad_d);    
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
       Bav = I(Bpix(1)-filt_rad_b:Bpix(1) + filt_rad_b, ...
               Bpix(2)-filt_rad_b:Bpix(2) + filt_rad_b);      
    end
    B(l) = mean(Bav(:));
end
pad = length(V)-length(B);
if pad
    B = [B, ones(1, pad)*B(end)];
end

V = remove_overlaps(X_r, X_d, V, ov_thr);
B = remove_overlaps(X_r, X_b, B, ov_thr);

%% Sub-functions
function F = remove_overlaps(X_r, X, F, ov_thr)
% Remove artefacts in fluorescence signal from overlapping regnerating
% axon, and interpolat across removed interval.
dd = zeros(1, size(X, 2));
for k = 1:length(dd)
    dd(k) = min(sqrt(sum((X(:, k) - X_r).^2, 1)));
end
ov_ind = find(dd < ov_thr);
if length(ov_ind) > 2
    sec_ind = find(diff(ov_ind) > 1);
    if ~isempty(sec_ind)
        ov = [ov_ind(1), ov_ind(sec_ind+1); ov_ind(sec_ind), ov_ind(end)];
    else
        ov = [ov_ind(1); ov_ind(end)];
    end
        for k = 1:size(ov, 2)
            s = ov(1, k):ov(2, k);
            if length(s) > 1
                f_start = mean(F(s(1)-10:s(1)));
                f_end = mean(F(s(end):s(end)+10));
                f_interp = interp1([s(1), s(end)], [f_start, f_end], s);
                F(s) = f_interp;
            end
        end
end
