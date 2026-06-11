# WP0 — Estate-Tax-Financed Capital Grants and the German Wealth Distribution

MATLAB implementation of a continuous-time heterogeneous-agent life-cycle
model with bequests, for the DIW preliminary paper. Specification documents:
`docs/preliminary_paper_outline.md` and `docs/implementation_appendix.md`
(the appendix is the implementation contract).

## Setup

- MATLAB R2021a+ (developed on R2026a).
- Run `startup.m` from the project root at the start of each session.
- `external/` is git-ignored and holds read-only copies of the reference
  libraries. To recreate it:
  - `git clone https://github.com/schaab-lab/SparseEcon external/SparseEcon`
    (developed against commit `9d04a89`, Feb 2023)
  - Moll reference scripts from <https://benjaminmoll.com/codes/> into
    `external/Moll-codes/` (see CLAUDE.md for the list)

## Layout

Per the implementation appendix §6: model code in `src/`, policy experiments
in `experiments/`, validation in `tests/`, paper output in `postprocess/`,
German source data in `data/`, run output in `results/` (git-ignored).

The `params` struct returned by `src/params_default.m` is the single source
of truth for parameters and Step 1 ↔ Step 2 module toggles.

## Status

Session 1 (bootstrap/scaffold) complete. See `CLAUDE.md` for current state
and next steps.
