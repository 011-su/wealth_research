function params = params_derive(params)
% PARAMS_DERIVE  Recompute derived model-unit fields from the EUR inputs.
%
%   params = params_derive(params)
%
% Call after changing any *_eur field (or eur_per_unit, h0, h_max) so the
% model-unit quantities stay consistent. params_default() calls this at the
% end; experiment scripts must call it again after overriding EUR inputs.

params.a_max         = params.a_max_eur         / params.eur_per_unit;
params.F             = params.F_eur             / params.eur_per_unit;
params.G             = params.G_eur             / params.eur_per_unit;
params.pension       = params.pension_eur       / params.eur_per_unit;
params.bequest_shift = params.bequest_shift_eur / params.eur_per_unit;
params.Nh            = params.h_max - params.h0 + 1;  % informational; grids_build derives from h0/h_max

end
