function lambda = mortality_step2(grids, params)
% MORTALITY_STEP2  Block C.2: age-specific hazard from the Destatis Sterbetafel.
%
%   lambda = mortality_step2(grids, params)
%
% Returns lambda (Nh x 1): the continuous-time mortality hazard at each model
% age. Reads data/derived/mortality_hazard.csv (built by
% data/make_mortality_hazard.py from the Destatis 2022/2024 single-year,
% all-Germany period life table, sex-pooled by life-table survivors), which
% stores lambda(age) = -ln(1 - q_pooled(age)) indexed by BIOLOGICAL age, and
% maps it onto the model age grid grids.h.
%
% In the current code model age == biological age (h0 = 18), so the map is the
% identity; interp1 keeps it correct if the age convention changes. Replaces
% the constant Step 1 hazard, concentrating deaths at old ages (fixing the
% spurious aging-out pile-up at h_max under constant mortality).

f = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
             'data', 'derived', 'mortality_hazard.csv');
if ~exist(f, 'file')
    error('mortality_step2:nodata', ...
        ['%s not found. Run: python3 data/make_mortality_hazard.py'], f);
end
M = readmatrix(f);                       % cols: age, qx_m, qx_f, qx_pooled, lambda
age_d = M(:, 1);
lam_d = M(:, end);

lambda = interp1(age_d, lam_d, grids.h, 'linear', 'extrap');
lambda = max(lambda(:), 0);

end
