%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is for coregistration footprints from hungry vs sated Mouse D
% Assume A1, A2 already loaded and permuted: [height, width, n_cells]
% Assume cell_registered_struct loaded (load the output from Cell Reg, from
% whichever mouse you want to run (MouseD_CoReg.mat)

load("/Users/nitishpatel/Desktop/Miniscope Feeding/Hunger_CoReg/Mouse_D/D_CoReg.mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_D/Week 2/Session 2/D_hungry.mat");
%load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_D/Week 1/Session 3/D_hungry.mat");




A1 = permute(A, [2, 3, 1]);  % Now A is [280, 300, 43] Mouse from Week 2,
%Session 2
A2 = permute(A, [2, 3, 1]);  % Now A is [280, 300, 96] Mouse from Week 1,
%Session 3





% === Parameters ===
sigma = 2;
threshold = 0.2;
alpha_val = 0.6;
[H, W, ~] = size(A1);

% === Colors ===
yellow = [249, 218, 120]/255;
blue   = [65, 111, 178]/255;
white  = [0.3, 0.3, 0.3];

% === CellReg Mapping ===
map = cell_registered_struct.cell_to_index_map;
co_registered_idx = find(map(:,1) > 0 & map(:,2) > 0);  % co-registered only
n_coreg = numel(co_registered_idx);
n_coreg = numel(co_registered_idx);
% === Pastel spectrum: red → violet in soft tones ===
% === Pastel spectrum: red → violet ===
hues = linspace(0, 0.83, n_coreg);   % red to violet (HSV hue)
saturation = 0.4;                    % reduced saturation = pastel
value = 1.0;                         % max brightness allowed by hsv2rgb

% === Helper function ===
normalize = @(x) x / max(x(:));
cmap = hsv2rgb([hues(:), repmat(saturation, n_coreg, 1), repmat(value, n_coreg, 1)]);


% === ----- Plot A1 Figure ----- === THIS IS MOUSE D, WEEK 2, SESSION 2
overlay_rgb_A1 = zeros(H, W, 3);
overlay_alpha_A1 = zeros(H, W);

coreg_counter = 1;  % Index for cmap

for i = 1:size(map,1)
    idx1 = map(i, 1);
    idx2 = map(i, 2);

    if idx1 == 0, continue; end

    img = normalize(imgaussfilt(A1(:,:,idx1), sigma));
    mask = img > threshold;

    if idx2 > 0
        color = cmap(coreg_counter, :);  % Use next color from cmap
        coreg_counter = coreg_counter + 1;
    else
        color = white;  % A1-only
    end

    for ch = 1:3
        channel = overlay_rgb_A1(:,:,ch);
        channel(mask) = channel(mask) + color(ch);
        overlay_rgb_A1(:,:,ch) = channel;
    end
    overlay_alpha_A1(mask) = overlay_alpha_A1(mask) + alpha_val;
end

overlay_alpha_A1 = min(overlay_alpha_A1, 1);

figure('Color', 'w');
imshow(overlay_rgb_A1);
hold on;
h = imshow(overlay_rgb_A1);
set(h, 'AlphaData', overlay_alpha_A1);
axis image off;
exportgraphics(gca,'MouseD_HungryFootprints.eps','ContentType','vector')


% === Parameters ===
sigma = 2;
threshold = 0.2;
alpha_val = 0.6;
[H, W, ~] = size(A2);

% === Colors ===
yellow = [249, 218, 120]/255;
blue   = [65, 111, 178]/255;
white  = [0.3, 0.3, 0.3];

% === CellReg Mapping ===
map = cell_registered_struct.cell_to_index_map;
co_registered_idx = find(map(:,1) > 0 & map(:,2) > 0);  % co-registered only
n_coreg = numel(co_registered_idx);
% === Pastel spectrum: red → violet in soft tones ===
% === Pastel spectrum: red → violet ===
hues = linspace(0, 0.83, n_coreg);   % red to violet (HSV hue)
saturation = 0.4;                    % reduced saturation = pastel
value = 1.0;                         % max brightness allowed by hsv2rgb

% Convert HSV to pastel RGB
cmap = hsv2rgb([hues(:), repmat(saturation, n_coreg, 1), repmat(value, n_coreg, 1)]);

% === Helper function ===
normalize = @(x) x / max(x(:));

% === ----- Plot A2 Figure ----- === THIS IS MOUSE D, WEEK 1, SESSION 3
overlay_rgb_A2 = zeros(H, W, 3);
overlay_alpha_A2 = zeros(H, W);

coreg_counter = 1;  % Index for cmap

for i = 1:size(map,1)
    idx1 = map(i, 1);
    idx2 = map(i, 2);

    if idx2 == 0, continue; end  % skip if not present in A2

    img = normalize(imgaussfilt(A2(:,:,idx2), sigma));
    mask = img > threshold;

    if idx1 > 0
        color = cmap(coreg_counter, :);  % Use next color from cmap
        coreg_counter = coreg_counter + 1;
    else
        color = white ;  % A2-only
    end

    for ch = 1:3
        channel = overlay_rgb_A2(:,:,ch);
        channel(mask) = channel(mask) + color(ch);
        overlay_rgb_A2(:,:,ch) = channel;
    end
    overlay_alpha_A2(mask) = overlay_alpha_A2(mask) + alpha_val;
end

overlay_alpha_A2 = min(overlay_alpha_A2, 1);

figure('Color', 'w');
imshow(overlay_rgb_A2);
hold on;
h = imshow(overlay_rgb_A2);
set(h, 'AlphaData', overlay_alpha_A2);
axis image off;
exportgraphics(gca,'MouseD_SatedFootprints.eps','ContentType','vector')


