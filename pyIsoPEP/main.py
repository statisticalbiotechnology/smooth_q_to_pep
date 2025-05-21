#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd
from .IsotonicPEP import IsotonicPEP

def parse_args():
    parser = argparse.ArgumentParser(
        prog="IsotonicPEP",
        description="Estimate PEPs from empirical null models using isotonic regression.",
        epilog="""
Example usage:
  1) Q2PEP:
     a) a concatenated target and decoy input file
     ./main.py q2pep --target-file ../example/peptide.target.txt --calc-q-from-pep --output ../example/results
     b) separate target and decoy input files
     ./main.py q2pep --target-file ../example/peptide.target.txt --decoy-file ../example/peptide.decoy.txt --label-col type --target-label 0 --decoy-label 1 --calc-q-from-fdr --calc-q-from-pep --output ../example/results

  2) D2PEP:
     a) a concatenated target and decoy input file
     ./main.py d2pep --cat-file ../example/peptide.cat.txt --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-pep --output ../example/results
     b) separate target and decoy input files
     ./main.py d2pep --target-file ../example/peptide.target.txt --decoy-file ../example/peptide.decoy.txt --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-fdr --calc-q-from-pep --output ../example/results
  """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Subparsers for method selection (q2pep or d2pep)
    subparsers = parser.add_subparsers(dest="method", required=True, help="PEP estimation method: 'q2pep' or 'd2pep'")
    
    # q2pep parser
    parser_q = subparsers.add_parser("q2pep", help="Estimate PEPs from q-values (q2pep method).")
    parser_q.add_argument("--qcol", type=str, default="q-value",
                          help="Column name for q-values in the input file (default: 'q-value').")
    parser_q.add_argument("--ip", action="store_true",
                          help="Apply monotonic interpolation (only used when regression-algo is 'PAVA', default: False).")
    parser_q.add_argument("--ip-algo", type=str, choices=["ispline", "pchip"], default="ispline",
                          help="Interpolation algorithm to use on PAVA-derived block centers when --ip is used (default: ispline).")
    parser_q.add_argument("--center-method", type=str, choices=["mean", "median"], default="mean",
                          help="Method for computing block centers (default: mean).")

    # d2pep parser
    parser_d = subparsers.add_parser("d2pep", help="Estimate PEPs from the target-decoy competition-resulting list (d2pep method).")

    # Input selection: either concatenated or separate files
    for p in (parser_q, parser_d):
        input_group = p.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--cat-file", metavar="FILE", help="Path to a concatenated target and decoy TSV file, containing score and label columns (used in concatenated input mode).")
        input_group.add_argument("--target-file", metavar="FILE", help="Path to the TSV file containing target scores (used in separate input mode).")
        # If in separate mode, decoy is optional for q2pep without --calc-q-from-fdr, required otherwise
        p.add_argument("--decoy-file", metavar="FILE", help="Path to the TSV file containing decoy scores (required in separate input mode).")

        # Shared options
        p.add_argument("--score-col", type=str, default="score",
                            help="Column name for scores (default: 'score').")
        p.add_argument("--label-col", type=str, default="label",
                            help="Column name for target/decoy labels(default: 'label').")
        p.add_argument("--target-label", type=str, default="target",
                            help="Target label in concatenated input file (default: 'target').")
        p.add_argument("--decoy-label", type=str, default="decoy",
                            help="Decoy label in concatenated input file (default: 'decoy').")
        
        p.add_argument("--regression-algo", type=str, choices=["PAVA", "ispline"], default="ispline",
                            help="Regression algorithm to use: 'PAVA' or 'ispline' (default: ispline).")

        p.add_argument("--calc-q-from-fdr", action="store_true", dest="calc_q_from_fdr",
                            help="Estimate FDRs and q-values from target and decoy scores before calculating PEPs. Required if no q-values are provided to q2pep; optional when q-values are provided to q2pep or for d2pep. (default: disabled).")
        p.add_argument("--calc-q-from-pep", action="store_true", dest="calc_q_from_pep",
                            help="Calculate q-values from estimated PEPs (default: disabled).")
        p.add_argument("--output", type=str, required=True, metavar="FILE|DIR", 
                            help="Output file path (or output directory; if directory, a default name 'outputPEP.target.dbased.txt' or 'outputPEP.target.qbased.txt' will be used. Note: only target PEPs are saved.)")
        
        p.add_argument("--verbose", action="store_true",
                            help="Print parameter information.")
    
    return parser.parse_args()

def load_data(args):
    """
    Load input data and return a tuple (df_target, obs):
      - df_target: DataFrame of target entries preserving input order.
      - obs: numpy array shape (n,2) of [score, label_numeric] for all entries (target and decoy), or None.

    In concatenated input mode, interpret --label-col/--target-label/--decoy-label.
    In separate input mode, read target and optionally decoy, assign labels automatically.
    """
    if args.cat_file:   # concatenated input mode
        df = pd.read_csv(args.cat_file, sep="\t")
        for col in (args.score_col, args.label_col):
            if col not in df.columns:
                sys.exit(f"Column '{col}' is missing from --cat-file.")

        df = df.copy()
        df['pyIsoPEP label'] = df[args.label_col].astype(str).str.lower()
        df_target = df[df["pyIsoPEP label"] == args.target_label.lower()].reset_index(drop=True)
        label_map = {args.target_label.lower(): 0.0, args.decoy_label.lower(): 1.0}
        num = df["pyIsoPEP label"].map(label_map)
        if num.isnull().any():
            sys.exit("Unrecognised labels in --cat-file (check --target-label/--decoy-label).")
        scores = df[args.score_col].astype(float).values
        obs = np.column_stack([scores, num.values])

        return df_target, obs
    
    # separate input mode
    df_target = pd.read_csv(args.target_file, sep="\t")
    if args.score_col not in df_target.columns:
        sys.exit(f"Column '{args.score_col}' missing from --target-file.")
    df_target = df_target.copy()
    df_target['pyIsoPEP label'] = 'target'

    if args.decoy_file:
        df_decoy = pd.read_csv(args.decoy_file, sep="\t")
        if args.score_col not in df_decoy.columns:
            sys.exit(f"Column '{args.score_col}' missing from --decoy-file.")
        df_decoy = df_decoy.copy()
        df_decoy['pyIsoPEP label'] = 'decoy'

        df = pd.concat([df_target, df_decoy], ignore_index=True)
        num = df['pyIsoPEP label'].map({'target': 0.0, 'decoy': 1.0}).values
        scores = df[args.score_col].astype(float).values
        obs = np.column_stack([scores, num])

        return df_target, obs
    else:   # target‑only input (allowed for q2pep when q‑values are provided)
        return df_target.reset_index(drop=True), None

def main():
    args = parse_args()

    if args.verbose:
        print("Parameters:")
        for arg, value in sorted(vars(args).items()):
            print(f"  {arg}: {value}")

    # Load data
    df_target, obs = load_data(args)

    # Initialize PEP processor
    pep_regressor = IsotonicPEP(
        regression_algo=args.regression_algo,
        ip_algo=getattr(args, "ip_algo", None),
        center_method=getattr(args, "center_method", None)
    )

    if args.method == 'q2pep':
        if args.calc_q_from_fdr and obs is None:
            sys.exit("--calc-q-from-fdr requires decoy information.")
        if not args.calc_q_from_fdr and args.qcol not in df_target.columns:
            sys.exit(f"Missing q-value column '{args.qcol}'. Either supply it or add --calc-q-from-fdr together with decoy information.")
        fdr_arr, q1_arr, pep_arr, q2_arr = pep_regressor.pep_regression(
            method="q2pep",
            q_values=df_target[args.qcol].astype(float).values if not args.calc_q_from_fdr else None,
            obs=obs,
            regression_algo=args.regression_algo,
            ip=getattr(args, "ip", False),
            ip_algo=getattr(args, "ip_algo", None),
            center_method=getattr(args, "center_method", None),
            calc_q_from_fdr=args.calc_q_from_fdr,
            calc_q_from_pep=args.calc_q_from_pep,
        )

    else:   # d2pep
        if obs is None:
            sys.exit("d2pep requires target AND decoy scores (use --target-file and --decoy-file, or --cat-file).")
        fdr_arr, q1_arr, pep_arr, q2_arr = pep_regressor.pep_regression(
            obs=obs,
            method='d2pep',
            regression_algo=args.regression_algo,
            calc_q_from_fdr=args.calc_q_from_fdr,
            calc_q_from_pep=args.calc_q_from_pep
        )

    if args.calc_q_from_fdr and fdr_arr is not None:
        df_target['pyIsoPEP FDR'] = fdr_arr
        df_target['pyIsoPEP q-value from FDR'] = q1_arr
    df_target['pyIsoPEP PEP'] = pep_arr
    if args.calc_q_from_pep and q2_arr is not None:
        df_target['pyIsoPEP q-value from PEP'] = q2_arr
    df_target = df_target.drop(columns=['pyIsoPEP label'])

    out_path = args.output
    if os.path.isdir(out_path):
        suffix = "qbased" if args.method == "q2pep" else "dbased"
        out_path = os.path.join(out_path, f"outputPEP.target.{suffix}.txt")
    try:
        df_target.to_csv(out_path, sep="\t", index=False)
        print(f"Saved target results to: {out_path}")
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

if __name__ == "__main__":
    main()
