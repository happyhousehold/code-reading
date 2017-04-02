load ../../../../Results/PPI/735640.8990/PPI1Linear;
imS = [sense_map(:, :, 1), sense_map(:, :, 2), sense_map(:, :, 3), sense_map(:, :, 4);
    sense_map(:, :, 5), sense_map(:, :, 6), sense_map(:, :, 7), sense_map(:, :, 8)];
h = figure;
imshow(abs(imS), []);
print(h, '-deps', '../../../../Results/PPI/PPI1_smap');
h = figure;
imshow(abs(xTrue), []);
funCropEdge(h);
print(h, '-deps', '../../../../Results/PPI/PPI1_true');
h = figure;
imshow(abs(fftshift(mask(:,:,1))), []);
funCropEdge(h);
print(h, '-deps', '../../../../Results/PPI/PPI1_LinearMask');
load ../../../../Results/PPI/735640.8990/PPI1Poisson;
h = figure;
imshow(abs(fftshift(mask(:,:,1))), []);
funCropEdge(h);
print(h, '-deps', '../../../../Results/PPI/PPI1_PoissonMask');

load ../../../../Results/PPI/735640.8990/PPI2Linear;
imS = [sense_map(:, :, 1), sense_map(:, :, 2), sense_map(:, :, 3), sense_map(:, :, 4);
    sense_map(:, :, 5), sense_map(:, :, 6), sense_map(:, :, 7), sense_map(:, :, 8)];
h = figure;
imshow(abs(imS), []);
print(h, '-deps', '../../../../Results/PPI/PPI2_smap');
h = figure;
imshow(abs(xTrue), []);
funCropEdge(h);
print(h, '-deps', '../../../../Results/PPI/PPI2_true');
h = figure;
imshow(abs(fftshift(mask(:,:,1))), []);
funCropEdge(h);
print(h, '-deps', '../../../../Results/PPI/PPI2_LinearMask');
load ../../../../Results/PPI/735640.8990/PPI2Poisson;
h = figure;
imshow(abs(fftshift(mask(:,:,1))), []);
funCropEdge(h);
print(h, '-deps', '../../../../Results/PPI/PPI2_PoissonMask');
