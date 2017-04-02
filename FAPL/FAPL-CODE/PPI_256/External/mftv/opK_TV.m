function Kx = opK_TV(x)
    [nRow, nCol] = size(x);
    Kx = zeros(nRow, nCol, 2);
    Kx(:, :, 1) = [diff(x, 1, 1); zeros(1, nCol)];
    Kx(:, :, 2) = [diff(x, 1, 2), zeros(nRow, 1)];
end
