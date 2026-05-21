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
        description="Estimate monotone PEPs.",
        epilog="""
6 supported combinations (pathway, regression x-axis, regressor):

  1) Default: q-value pathway, rank-based regression using I-splines
     pyisopep q2pep --target-file example/peptide.target.txt --qcol q-value --calc-q-from-pep --output results/

  2) q-value pathway, rank-based regression using PAVA
     pyisopep q2pep --pava --target-file example/peptide.target.txt --qcol q-value --calc-q-from-pep --output results/

  3) q-value pathway, score-based regression using I-splines
     pyisopep q2pep --score-based --target-file example/peptide.target.txt --qcol q-value --score-col score --calc-q-from-pep --output results/

  4) q-value pathway, score-based regression using PAVA
     pyisopep q2pep --score-based --pava --target-file example/peptide.target.txt --qcol q-value --score-col score --calc-q-from-pep --output results/

  5) TDC pathway, score-based regression using I-splines
     pyisopep d2pep --target-file example/peptide.target.txt --decoy-file example/peptide.decoy.txt --score-col score --calc-q-from-pep --output results/

  6) TDC pathway, score-based regression using PAVA
     pyisopep d2pep --pava --target-file example/peptide.target.txt --decoy-file example/peptide.decoy.txt --score-col score --calc-q-from-pep --output results/

  Notes:
  - TDC pathway always requires decoy information.
  - For the q-value pathway, decoy information is only needed when computing
    q-values from TDC scores via --calc-q-from-fdr.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="method", required=True, help="PEP estimation method: 'q2pep' or 'd2pep'")

    parser_q = subparsers.add_parser("q2pep", help="Estimate PEPs from target q-values.")
    parser_q.add_argument("--qcol", type=str, default="q-value",
                          help="Column name for q-values in the input file (default: 'q-value').")
    parser_q.add_argument("--score-based", action="store_true", dest="score_based",
                          help="Use raw scores as the independent variable instead of rank.")
    parser_q.add_argument("--trim-plateaus", action="store_true", dest="trim_plateaus",
                          default=False,
                          help="Enable trimming of leading/trailing q-value plateaus before "
                               "isotonic regression (default: trimming disabled).")
    parser_q.add_argument("--no-pseudo-count", action="store_false", dest="pseudo_count",
                          default=True,
                          help="Disable the Jeffreys-style smoothing added to raw PEPs before "
                               "isotonic regression (default: enabled, value 0.5 distributed as "
                               "0.5 / n_mid per rank).")

    parser_d = subparsers.add_parser("d2pep", help="Estimate PEPs from TDC scores.")

    for p in (parser_q, parser_d):
        input_group = p.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--cat-file", metavar="FILE",
                                 help="Concatenated target+decoy tsv file.")
        input_group.add_argument("--target-file", metavar="FILE",
                                 help="Target-only tsv file (separate input mode).")
        p.add_argument("--decoy-file", metavar="FILE",
                       help="Decoy tsv file (separate input mode).")

        p.add_argument("--score-col", type=str, default="score",
                       help="Score column name (default: 'score').")
        p.add_argument("--label-col", type=str, default="label",
                       help="Label column name in concatenated file (default: 'label').")
        p.add_argument("--target-label", type=str, default="target",
                       help="String marking target rows (default: 'target').")
        p.add_argument("--decoy-label", type=str, default="decoy",
                       help="String marking decoy rows (default: 'decoy').")

        p.add_argument("--pava", action="store_true",
                       help="Use PAVA regression instead of the default I-splines.")
        p.add_argument("--calc-q-from-fdr", action="store_true", dest="calc_q_from_fdr",
                       help="Estimate FDRs and q-values from TDC scores before PEP estimation.")
        p.add_argument("--calc-q-from-pep", action="store_true", dest="calc_q_from_pep",
                       help="Derive q-values from estimated PEPs.")
        p.add_argument("-k", "--keep-input-order", action="store_true", dest="keep_input_order",
                       help="Write output in original input order (default: sort by ascending PEP).")
        p.add_argument("--output", type=str, required=True, metavar="FILE|DIR",
                       help="Output path. If a directory, a default filename is used. "
                            "Only targets are written.")
        p.add_argument("--verbose", action="store_true",
                       help="Print parameter information.")

    return parser.parse_args()


def load_data(args, require_score):
    if args.cat_file:
        df = pd.read_csv(args.cat_file, sep="\t")
        if args.label_col not in df.columns:
            sys.exit(f"Column '{args.label_col}' is missing from --cat-file.")
        if require_score and args.score_col not in df.columns:
            sys.exit(f"Column '{args.score_col}' is missing from --cat-file.")
        df = df.copy()
        df["pyIsoPEP label"] = df[args.label_col].astype(str).str.lower()
        label_map = {args.target_label.lower(): 0.0, args.decoy_label.lower(): 1.0}
        if not set(df["pyIsoPEP label"]).issubset(label_map.keys()):
            sys.exit("Unrecognised labels in --cat-file (check --target-label/--decoy-label).")

        df_target = df[df["pyIsoPEP label"] == args.target_label.lower()].reset_index(drop=True)
        if not require_score:
            return df_target, None
        scores = df[args.score_col].astype(float).values
        num = df["pyIsoPEP label"].map(label_map)
        obs = np.column_stack([scores, num.values])
        return df_target, obs

    # Separate input mode
    df_target = pd.read_csv(args.target_file, sep="\t").copy()
    df_target["pyIsoPEP label"] = "target"
    if not require_score:
        return df_target, None

    if args.decoy_file is None:
        sys.exit("--decoy-file is required.")
    df_decoy = pd.read_csv(args.decoy_file, sep="\t").copy()
    df_decoy["pyIsoPEP label"] = "decoy"
    if args.score_col not in df_target.columns:
        sys.exit(f"Column '{args.score_col}' missing from --target-file.")
    if args.score_col not in df_decoy.columns:
        sys.exit(f"Column '{args.score_col}' missing from --decoy-file.")
    df = pd.concat([df_target, df_decoy], ignore_index=True)
    num = df["pyIsoPEP label"].map({"target": 0.0, "decoy": 1.0}).values
    scores = df[args.score_col].astype(float).values
    obs = np.column_stack([scores, num])
    return df_target.reset_index(drop=True), obs


def main():
    args = parse_args()

    if args.verbose:
        print("Parameters:")
        for arg, value in sorted(vars(args).items()):
            print(f"  {arg}: {value}")

    score_based = getattr(args, "score_based", False)

    # Decoy information (obs) is only needed when computing q-values from TDC
    # scores (--calc-q-from-fdr) or running the TDC pathway (d2pep).
    if args.method == "q2pep":
        require_score = args.calc_q_from_fdr
    else:  # d2pep
        require_score = True

    df_target, obs = load_data(args, require_score)

    if score_based and args.score_col not in df_target.columns:
        sys.exit(f"Column '{args.score_col}' missing from target input (required for --score-based).")

    pep_regressor = IsotonicPEP(pava=args.pava)

    if args.method == "q2pep":
        if score_based:
            qcol = args.qcol
            if not args.calc_q_from_fdr and qcol not in df_target.columns:
                sys.exit(f"Missing q-value column '{qcol}'. Either supply it or add --calc-q-from-fdr together with decoy information.")
            fdr_arr, q1_arr, pep_arr, q2_arr = pep_regressor.pep_regression(
                method="qns2pep",
                q_values=df_target[qcol].astype(float).values if not args.calc_q_from_fdr else None,
                obs=obs,
                target_scores=df_target[args.score_col].astype(float).values,
                pava=args.pava,
                calc_q_from_fdr=args.calc_q_from_fdr,
                calc_q_from_pep=args.calc_q_from_pep,
                trim_plateaus=args.trim_plateaus,
                pseudo_count=args.pseudo_count,
            )
        else:
            if args.calc_q_from_fdr and obs is None:
                sys.exit("--calc-q-from-fdr requires decoy information.")
            if not args.calc_q_from_fdr and args.qcol not in df_target.columns:
                sys.exit(f"Missing q-value column '{args.qcol}'. Either supply it or add --calc-q-from-fdr together with decoy information.")
            fdr_arr, q1_arr, pep_arr, q2_arr = pep_regressor.pep_regression(
                method="q2pep",
                q_values=df_target[args.qcol].astype(float).values if not args.calc_q_from_fdr else None,
                obs=obs,
                pava=args.pava,
                calc_q_from_fdr=args.calc_q_from_fdr,
                calc_q_from_pep=args.calc_q_from_pep,
                trim_plateaus=args.trim_plateaus,
                pseudo_count=args.pseudo_count,
            )

    else:  # d2pep subcommand
        if obs is None:
            sys.exit("d2pep requires target AND decoy scores (--target-file + --decoy-file, or --cat-file).")
        fdr_arr, q1_arr, pep_arr, q2_arr = pep_regressor.pep_regression(
            obs=obs,
            method="tdc2pep",
            pava=args.pava,
            calc_q_from_fdr=args.calc_q_from_fdr,
            calc_q_from_pep=args.calc_q_from_pep,
        )

    if args.calc_q_from_fdr and fdr_arr is not None:
        df_target["pyIsoPEP FDR"] = fdr_arr
        df_target["pyIsoPEP q-value from FDR"] = q1_arr
    df_target["pyIsoPEP PEP"] = pep_arr
    if args.calc_q_from_pep and q2_arr is not None:
        df_target["pyIsoPEP q-value from PEP"] = q2_arr
    df_target = df_target.drop(columns=["pyIsoPEP label"])

    if not args.keep_input_order:
        by = ["pyIsoPEP PEP"]
        ascending = [True]
        if args.score_col in df_target.columns:
            by.append(args.score_col)
            ascending.append(False)
        df_target = df_target.sort_values(by=by, ascending=ascending, kind="mergesort").reset_index(drop=True)

    out_path = args.output
    if os.path.isdir(out_path):
        if args.method == "q2pep":
            suffix1 = "qns2pep" if score_based else "q2pep"
        else:
            suffix1 = "tdc2pep"
        if args.pava:
            suffix2 = "pava"
        else:
            suffix2 = "ispline"
        out_path = os.path.join(out_path, f"outputPEP.target.{suffix1}.{suffix2}.txt")
    try:
        df_target.to_csv(out_path, sep="\t", index=False)
        print(f"Saved target results to: {out_path}")
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")


if __name__ == "__main__":
    main()
