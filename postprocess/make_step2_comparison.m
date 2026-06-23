function T = make_step2_comparison()
% MAKE_STEP2_COMPARISON  Compare the wealth distributions across all Step 2
% (Destatis mortality) Grunderbe results (no-exemption surtax) and the Step 2
% status quo.
%
%   T = make_step2_comparison()
%
% Loads results/step1/flatgrant_{20k,100k,200k}_step2.mat, prints a table of
% key metrics, and saves a 2-panel figure overlaying, for each grant size:
%   (1) the wealth-over-age marginal (mean wealth by age), and
%   (2) the non-cumulative wealth distribution (density g(a)).

proj    = fileparts(fileparts(mfilename('fullpath')));
out_dir = fullfile(proj, 'results', 'step2');
files   = {'flatgrant_20k_step2.mat', 'flatgrant_100k_step2.mat', 'flatgrant_200k_step2.mat'};
glabels = {'Grunderbe 20k', 'Grunderbe 100k', 'Grunderbe 200k'};

% Status quo + grids (same grid for every run) from the first file
S0     = load(fullfile(out_dir, files{1}), 'out', 'params');
params = S0.params;  grids = grids_build(params);  eur = params.eur_per_unit;

moms   = {S0.out.mom_sq};
labels = {'status quo'};
G = 0; tau = 0;
b50 = S0.out.mom_sq.bottom50; gini = S0.out.mom_sq.gini; meanw = S0.out.mom_sq.wealth_mean_eur;
for k = 1:numel(files)
    S = load(fullfile(out_dir, files{k}), 'G_eur', 'out');
    moms{end+1}   = S.out.mom;          %#ok<AGROW>
    labels{end+1} = glabels{k};         %#ok<AGROW>
    G(end+1)   = S.G_eur/1e3;           %#ok<AGROW>
    tau(end+1) = S.out.tau_add;         %#ok<AGROW>
    b50(end+1) = S.out.mom.bottom50;    %#ok<AGROW>
    gini(end+1)= S.out.mom.gini;        %#ok<AGROW>
    meanw(end+1)= S.out.mom.wealth_mean_eur; %#ok<AGROW>
end

T = table(G(:), tau(:)*100, b50(:)*100, gini(:), meanw(:), ...
    'VariableNames', {'grant_kEUR','surtax_pct','bottom50_pct','gini','mean_wealth_EUR'});
disp('--- Step 2 (Destatis mortality), no-exemption surtax ---'); disp(T);

% ---- figure: wealth-by-age and non-cumulative wealth distribution ----
f = figure('visible','off','Position',[100 100 1150 470]);

subplot(1,2,1); hold on; grid on;
for i = 1:numel(moms)
    plot(grids.h, moms{i}.wealth_by_age * eur/1e3, 'LineWidth', 1.7);
end
xline(params.hR, '--', 'retirement', 'HandleVisibility','off');
xlabel('Age'); ylabel('mean wealth (1000 EUR)');
legend(labels, 'Location', 'northwest'); title('Wealth over age');

subplot(1,2,2); hold on; grid on;
a_keur = grids.a * eur/1e3;
for i = 1:numel(moms)
    dens = moms{i}.wealth_marginal ./ grids.wa;     % density g(a), non-cumulative
    plot(a_keur(2:end), dens(2:end), 'LineWidth', 1.7);
end
set(gca, 'YScale', 'log'); xlim([0 1500]);
xlabel('Wealth (1000 EUR)'); ylabel('density g(a)');
legend(labels, 'Location', 'northeast'); title('Wealth distribution (non-cumulative)');

sgtitle('Grunderbe under Step 2 mortality — wealth distribution by grant size');
exportgraphics(f, fullfile(out_dir, 'step2_comparison.png'), 'Resolution', 150);
close(f);
fprintf('saved figure: %s\n', fullfile(out_dir, 'step2_comparison.png'));

end
