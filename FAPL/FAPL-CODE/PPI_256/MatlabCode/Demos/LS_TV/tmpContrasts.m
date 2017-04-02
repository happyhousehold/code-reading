load ../../../../Results/PPI/735640.8990/PPI1Linear;
addpath(genpath('../../../../External/exact_alm_rpca'));
addpath(genpath('../../../../External/svt'));

imS = [sense_map(:, :, 1), sense_map(:, :, 2), sense_map(:, :, 3), sense_map(:, :, 4);
    sense_map(:, :, 5), sense_map(:, :, 6), sense_map(:, :, 7), sense_map(:, :, 8)];
xAll = abs(bsxfun(@times, xTrue, sense_map));
xAll(xAll==0) = eps;
xAllImg = [xAll(:, :, 1), xAll(:, :, 2), xAll(:, :, 3), xAll(:, :, 4);
    xAll(:, :, 5), xAll(:, :, 6), xAll(:, :, 7), xAll(:, :, 8)];
% figure;imshow(abs(xAllImg), []);
xAllLnImg = log(xAllImg);
% figure;imshow((xAllLnImg), []);
% figure;imshow(exp(xAllLnImg), []);
% figure;imshow(abs(log(imS)), [0, 10]);

%% Robust PCA
xAllLn = log(xAll);
tmp = (reshape(xAllLn(:), [nRow*nCol, nCh]));
% [A_hat, E_hat, iter] = exact_alm_rpca(tmp);
[tmpA,tmpE,tmpY] = singular_value_rpca(tmp, .002);
sum(sum(double(abs(tmpE)>0)))
% tmpA = tmpA';
% tmpE = tmpE';
%%
% tmpA(tmpA>=10) = 0;
% tmpE(tmpE>=10) = 0;
imA = reshape(tmpA, [nRow, nCol, nCh]);
imE = reshape(tmpE, [nRow, nCol, nCh]);
xRPCA = sqrt(sum(exp(imA).^2, 3)/nCh);
xRPCA = (xRPCA - min(xRPCA(:)))/(max(xRPCA(:)) - min(xRPCA(:)));
norm(xRPCA(:) - xTrue(:)) / norm(xTrue(:))
figure;
subplot(1, 3, 1);
imshow(xTrue, []);
subplot(1, 3, 2);
imshow(xRPCA, []);
subplot(1, 3, 3);
imshow(abs(xTrue-xRPCA), []); colorbar;
figure;
imshow(exp([imA(:, :, 1), imA(:, :, 2), imA(:, :, 3), imA(:, :, 4);...
    imA(:, :, 5), imA(:, :, 6), imA(:, :, 7), imA(:, :, 8)]),...
    [])
figure;
imshow(exp([imE(:, :, 1), imE(:, :, 2), imE(:, :, 3), imE(:, :, 4);...
    imE(:, :, 5), imE(:, :, 6), imE(:, :, 7), imE(:, :, 8)]),...
    [])
figure;
imshow(abs(imS), []);