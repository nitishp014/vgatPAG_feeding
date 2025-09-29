%%%%%%%%%%%%%%%%%%%
%This section is footprints of Mouse F, Week 2, Session 2 (Hungry)
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Hunger_CoReg/Mouse_F/F_CoReg.mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_F/Week 2/Session 2/F_hungry.mat");



A = permute(A, [2, 3, 1]);  % Now A is [280, 300, 43] (needs to be height x width x n_cells)

% === Load A ===
[nx, ny, n_cells] = size(A);

% === Threshold and assign labels ===
cell_labels = zeros(nx, ny);  % label matrix
for i = 1:n_cells
    cell_img = A(:,:,i);
    
    % Apply Gaussian smoothing with sigma=2 (adjust if needed)
    cell_img_smooth = imgaussfilt(cell_img, 2);
    
    % Normalize
    cell_img_smooth = cell_img_smooth / max(cell_img_smooth(:));
    
    % Threshold
    mask = cell_img_smooth > 0.2;
    
    % Assign label
    cell_labels(mask) = i;
end


% === Generate pastel version of 'jet' ===
%base_cmap = jet(n_cells);  % or use any base colormap
%pastel_cmap = 0.6 + 0.4 * base_cmap;  % shift colors toward white

% Now use it in label2rgb
%rgb_img = label2rgb(cell_labels, pastel_cmap, 'k', 'shuffle');
% === Create colored label image ===
rgb_img = label2rgb(cell_labels, 'parula', 'k', 'shuffle');  % 'k' = black background

% === Plot ===
figure('Color','w');
imshow(rgb_img);
axis image off;




% Above is just plotting footprints from a single mouse