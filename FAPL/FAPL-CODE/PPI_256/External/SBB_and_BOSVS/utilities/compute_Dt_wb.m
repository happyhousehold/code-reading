function RHS = compute_Dt_wb(Wx, Wy, bx, by, rho)

RHS = (DxtU(rho*Wx-bx)+DytU(rho*Wy-by));

end

% compute D'_x(U)
function dxtu = DxtU(U)
    dxtu = [U(:,end)-U(:, 1), U(:,1:end-1)-U(:,2:end)];
end

% compute D'_y(U)
function dytu = DytU(U)
    dytu = [U(end,:)-U(1, :); U(1:end-1,:)-U(2:end,:)];
end