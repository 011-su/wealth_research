function test_aiyagari_limit()
% TEST_AIYAGARI_LIMIT  Regression vs Moll's Huggett-diffusion PE script (appendix 5.1.3).
%
% Collapse age (Nh = 1), switch off mortality and bequests: the model is
% Huggett (PE, fixed r0) with log-OU diffusion income. The reference below
% is huggett_diffusion_partialeq.m with its commented log-OU branch
% activated, amin = 0, and our Step 1 parameters, solved on the identical
% grid -- so V and the stationary masses must agree to solver precision.
% (The appendix names the GE Aiyagari script as target; with fixed r0 the
% PE Huggett script is the correct limit, see CLAUDE.md discrepancy 3.)

p = params_default();
p.h_max      = p.h0;     % Nh = 1: age collapsed, infinite horizon
p.lambda_bar = 0;        % no mortality
p.theta_b    = 0;        % no bequest motive
p.G_eur      = 0;        p.G = 0;
p.a_curv     = 1;        % uniform wealth grid, as in the reference
p.a_max      = 100;      % model units, comparable to Moll's ranges
p.Na         = 300;
p.tol        = 1e-10;
p.maxit_hjb  = 200;

grids = grids_build(p);
ops   = operator_build(grids, p);
[V, ~, ~, A] = hjb_solve(grids, ops, p);
m = kfe_solve(grids, ops, p, A);

[V_ref, m_ref] = reference_huggett_logou(p);

err_V = max(abs(V - V_ref)) / max(abs(V_ref));
err_m = max(abs(m - m_ref));
fprintf('Aiyagari/Huggett limit: max rel V err %.2e, max mass err %.2e\n', ...
    err_V, err_m);

assert(err_V < 1e-8, 'value function deviates from Moll reference');
assert(err_m < 1e-10, 'stationary distribution deviates from Moll reference');
disp('PASS test_aiyagari_limit');

end


function [V_out, m_out] = reference_huggett_logou(p)
% huggett_diffusion_partialeq.m (Moll/Ahn), adapted minimally: log-OU branch
% (state y = log z, constant coefficients, exp(y) in the budget), amin = 0,
% parameters from p. Kept structurally as close to the original as possible.

ga = p.gamma; r = p.r0; rho = p.rho; w = p.w;
the = p.theta_y; sig2 = p.sigma_y^2;

I = p.Na;
amin = 0; amax = p.a_max;
a = linspace(amin, amax, I)';
da = (amax - amin) / (I - 1);

J = p.Ny;
sd = sqrt(sig2 / (2 * the));
y = linspace(p.mu_y - p.y_sd_span * sd, p.mu_y + p.y_sd_span * sd, J);
dy = (y(end) - y(1)) / (J - 1);
dy2 = dy^2;

mu = -the * (y - p.mu_y);    % OU drift in logs (no Ito term)
s2 = sig2 .* ones(1, J);

aa = a * ones(1, J);
zz = ones(I, 1) * exp(y);    % income levels

maxit = p.maxit_hjb;
crit  = p.tol;
Delta = p.Delta_hjb;

Vaf = zeros(I, J); Vab = zeros(I, J);

% B_switch: evolution of y (verbatim Moll construction)
chi  = -min(mu, 0) / dy + s2 / (2 * dy2);
yyc  =  min(mu, 0) / dy - max(mu, 0) / dy - s2 / dy2;
zeta =  max(mu, 0) / dy + s2 / (2 * dy2);

updiag = zeros(I, 1);
for j = 1:J
    updiag = [updiag; repmat(zeta(j), I, 1)]; %#ok<AGROW>
end
centdiag = repmat(chi(1) + yyc(1), I, 1);
for j = 2:J-1
    centdiag = [centdiag; repmat(yyc(j), I, 1)]; %#ok<AGROW>
end
centdiag = [centdiag; repmat(yyc(J) + zeta(J), I, 1)];
lowdiag = repmat(chi(2), I, 1);
for j = 3:J
    lowdiag = [lowdiag; repmat(chi(j), I, 1)]; %#ok<AGROW>
end
Aswitch = spdiags(centdiag, 0, I*J, I*J) + spdiags(lowdiag, -I, I*J, I*J) ...
        + spdiags(updiag, I, I*J, I*J);

v = (w * zz + r .* aa).^(1 - ga) / (1 - ga) / rho;

for n = 1:maxit
    V = v;
    Vaf(1:I-1, :) = (V(2:I, :) - V(1:I-1, :)) / da;
    Vaf(I, :) = (w * exp(y) + r .* amax).^(-ga);
    Vab(2:I, :) = (V(2:I, :) - V(1:I-1, :)) / da;
    Vab(1, :) = (w * exp(y) + r .* amin).^(-ga);

    cf = Vaf.^(-1 / ga);  sf = w * zz + r .* aa - cf;
    cb = Vab.^(-1 / ga);  sb = w * zz + r .* aa - cb;
    c0 = w * zz + r .* aa;
    Va0 = c0.^(-ga);

    If = sf > 0; Ib = sb < 0; I0 = (1 - If - Ib);
    Va_Upwind = Vaf .* If + Vab .* Ib + Va0 .* I0;

    c = Va_Upwind.^(-1 / ga);
    u = c.^(1 - ga) / (1 - ga);

    X = -min(sb, 0) / da;
    Y = -max(sf, 0) / da + min(sb, 0) / da;
    Z =  max(sf, 0) / da;

    updiag = 0;
    for j = 1:J
        updiag = [updiag; Z(1:I-1, j); 0]; %#ok<AGROW>
    end
    centdiag = reshape(Y, I*J, 1);
    lowdiag = X(2:I, 1);
    for j = 2:J
        lowdiag = [lowdiag; 0; X(2:I, j)]; %#ok<AGROW>
    end
    AA = spdiags(centdiag, 0, I*J, I*J) + spdiags([updiag; 0], 1, I*J, I*J) ...
       + spdiags([lowdiag; 0], -1, I*J, I*J);

    A = AA + Aswitch;
    B = (1 / Delta + rho) * speye(I*J) - A;

    b = reshape(u, I*J, 1) + reshape(V, I*J, 1) / Delta;
    V_stacked = B \ b;
    V = reshape(V_stacked, I, J);

    Vchange = V - v;
    v = V;
    if max(max(abs(Vchange))) < crit
        break
    end
end

% Fokker-Planck: fix one value, solve, normalise to unit mass
AT = A';
b = zeros(I*J, 1);
i_fix = 1;
b(i_fix) = .1;
row = [zeros(1, i_fix - 1), 1, zeros(1, I*J - i_fix)];
AT(i_fix, :) = row;
gg = AT \ b;

V_out = reshape(v, I*J, 1);
m_out = gg / sum(gg);

end
