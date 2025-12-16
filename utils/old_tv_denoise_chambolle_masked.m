function U = tv_denoise_chambolle_masked(f, mask, lambda, maxIter)
% TV denoising restricted to mask area
if nargin<4, maxIter = 200; end
tau = 0.125;
[M,N] = size(f);
px = zeros(M,N);
py = zeros(M,N);

for iter = 1:maxIter
    % divergence
    div_p = [px(1,:); diff(px);] + [py(:,1), diff(py,1,2)];
    u = f - lambda*div_p;
    u(~mask) = 0;

    % compute gradients only inside mask
    ux = zeros(M,N);
    uy = zeros(M,N);
    ux(1:end-1,:) = diff(u,1,1);
    uy(:,1:end-1) = diff(u,1,2);
    ux(~mask) = 0; uy(~mask) = 0;

    px_new = px + tau*ux;
    py_new = py + tau*uy;
    norm_p = max(1, sqrt(px_new.^2 + py_new.^2));
    px = px_new ./ norm_p;
    py = py_new ./ norm_p;

    px(~mask) = 0; py(~mask) = 0;
end

div_p = [px(1,:); diff(px);] + [py(:,1), diff(py,1,2)];
U = f - lambda*div_p;
U(~mask) = NaN;
end
