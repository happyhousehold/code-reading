function err = relerr(u, u0)
% compute relative error (RMS)

err = norm(u(:)-u0(:),'fro')/norm(u0(:),'fro');

return;