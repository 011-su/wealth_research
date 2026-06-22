#!/usr/bin/env python3
"""Build data/derived/mortality_hazard.csv from the Destatis Sterbetafel.

Extracts the single-year, all-Germany 2022/2024 period life table (sheets
csv-12613-b01 male, b02 female), pools the sexes weighting each sex's qx by
its life-table survivors lx, and converts the annual death probability q to a
continuous-time hazard  lambda = -ln(1 - q_pooled)  (constant within the year,
matching the annual survival probability). Output indexed by biological age.

Run once from the project root:  python3 data/make_mortality_hazard.py
"""
import csv, math, os
import openpyxl

XLSX = "data/statistischer-bericht-sterbetafeln-5126207247005.xlsx"
OUT  = "data/derived/mortality_hazard.csv"
# column indices (0-based) in the csv-* sheets: Alter=4, qx=5, lx=7
COL_AGE, COL_QX, COL_LX = 4, 5, 7


def parse(wb, sheet):
    out = {}
    for r in list(wb[sheet].iter_rows(values_only=True))[1:]:  # skip header
        if r[COL_AGE] is None:
            continue
        try:
            age = int(r[COL_AGE])
        except (TypeError, ValueError):
            continue
        out[age] = (float(r[COL_QX]), float(r[COL_LX]))
    return out


def main():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    male, female = parse(wb, "csv-12613-b01"), parse(wb, "csv-12613-b02")
    ages = sorted(set(male) & set(female))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["age", "qx_male", "qx_female", "qx_pooled", "lambda_hazard"])
        for a in ages:
            qm, lm = male[a]
            qf, lf = female[a]
            qp = (lm * qm + lf * qf) / (lm + lf)            # lx-weighted pool
            lam = -math.log(max(1e-12, 1.0 - qp))            # continuous hazard
            w.writerow([a, f"{qm:.8g}", f"{qf:.8g}", f"{qp:.8g}", f"{lam:.8g}"])

    print(f"wrote {OUT} for ages {ages[0]}-{ages[-1]} ({len(ages)} rows)")
    for a in (18, 40, 65, 80, 90, 99):
        if a in male:
            qm, lm = male[a]
            qf, lf = female[a]
            qp = (lm * qm + lf * qf) / (lm + lf)
            print(f"  age {a:3d}: qx_pooled={qp:.5f}  lambda={-math.log(1-qp):.5f}/yr")


if __name__ == "__main__":
    main()
