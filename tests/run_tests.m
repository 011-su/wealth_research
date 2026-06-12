function run_tests()
% RUN_TESTS  Run all validation tests in one session.
%
%   matlab -batch "startup; run_tests"
%
% Add new tests to the list below as they are implemented (appendix §5.1).

tests = {@test_aiyagari_limit, @test_mass_conservation};

n_fail = 0;
for k = 1:numel(tests)
    name = func2str(tests{k});
    fprintf('\n===== %s =====\n', name);
    try
        tests{k}();
    catch err
        n_fail = n_fail + 1;
        fprintf(2, 'FAIL %s: %s\n', name, err.message);
    end
end

if n_fail > 0
    error('run_tests:failures', '%d of %d tests failed', n_fail, numel(tests));
end
fprintf('\nAll %d tests passed.\n', numel(tests));

end
