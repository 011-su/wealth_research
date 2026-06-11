function grids = grids_build(params)
% GRIDS_BUILD  Construct the (a, y, h) state grids (appendix section 3.1).
%
%   grids = grids_build(params)
%
% Input:   params  struct from params_default()
% Output:  grids   struct with fields
%   .a   (Na x 1)  wealth grid on [0, a_max], power-spaced with curvature
%                  params.a_curv so points concentrate near a = 0
%   .y   (Ny x 1)  log-productivity grid, uniform on
%                  mu_y +/- y_sd_span * stationary OU sd
%   .h   (Nh x 1)  age grid, uniform, dh = 1 year
%   .da  (Na-1 x 1) forward wealth steps (nonuniform)
%   .dy, .dh       scalar grid steps
%   .aa, .yy, .hh  (N x 1) flattened state vectors, column-major over the
%                  (a, y, h) tensor: a varies fastest, h slowest
%   .Na, .Ny, .Nh, .N  dimensions
%
% Grid construction only -- sparse operators are assembled in
% operator_build.m (later session).

% Wealth: power-spaced, refinement near the borrowing constraint a = 0
u = linspace(0, 1, params.Na)';
grids.a  = params.a_max * u.^params.a_curv;
grids.da = diff(grids.a);

% Log productivity: uniform around the OU stationary distribution
% (stationary sd of dy = -theta (y - mu) dt + sigma dW is sigma/sqrt(2 theta))
sd_stat  = params.sigma_y / sqrt(2 * params.theta_y);
grids.y  = linspace(params.mu_y - params.y_sd_span * sd_stat, ...
                    params.mu_y + params.y_sd_span * sd_stat, params.Ny)';
grids.dy = grids.y(2) - grids.y(1);

% Age: uniform, one-year steps
grids.h  = (params.h0:1:params.h_max)';
grids.dh = 1;

% Dimensions
grids.Na = params.Na;
grids.Ny = params.Ny;
grids.Nh = numel(grids.h);
grids.N  = grids.Na * grids.Ny * grids.Nh;

% Flattened state vectors, column-major (Moll convention): reshaping an
% (Na x Ny x Nh) array with (:) gives a fastest, then y, then h
[A, Y, H] = ndgrid(grids.a, grids.y, grids.h);
grids.aa = A(:);
grids.yy = Y(:);
grids.hh = H(:);

end
