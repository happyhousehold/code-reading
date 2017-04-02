function Kty = opKt_TV(y)
    [nRow, nCol, ~] = size(y);
    D1p = diff([zeros(1, nCol); y(1:end-1, :, 1); zeros(1, nCol)], 1, 1);
    D2p = diff([zeros(nRow, 1), y(:, 1:end-1, 2), zeros(nRow, 1)], 1, 2);
    Kty = - D1p - D2p;
end
