function U = tv_denoise_chambolle(f, lambda, maxIter)
% Simple implementation of Chambolle's TV denoising for 2D images.
% Solves: min_u 0.5*||u - f||^2 + lambda * TV(u)
% f can be MxN image. Returns denoised U.
if nargin<3, maxIter = 200; end
tau = 0.125; % time step (stable)
[M,N] = size(f);
px = zeros(M,N); py = zeros(M,N);
div_p = zeros(M,N);
for iter=1:maxIter
    % compute divergence
    div_p = [px(1,:) ; px(2:end,:) - px(1:end-1,:)] + [py(:,1), py(:,2:end) - py(:,1:end-1)];
    % update dual variable
    u = f - lambda*div_p;
    % gradients of u
    ux = [u(2:end,:) - u(1:end-1,:); zeros(1,N)];
    uy = [u(:,2:end) - u(:,1:end-1), zeros(M,1)];
    px_new = px + tau*ux;
    py_new = py + tau*uy;
    norm_p = max(1, sqrt(px_new.^2 + py_new.^2));
    px = px_new ./ norm_p;
    py = py_new ./ norm_p;
end
div_p = [px(1,:) ; px(2:end,:) - px(1:end-1,:)] + [py(:,1), py(:,2:end) - py(:,1:end-1)];
U = f - lambda*div_p;
end
