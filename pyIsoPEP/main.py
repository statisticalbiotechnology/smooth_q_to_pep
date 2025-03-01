#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
from IsotonicPEP import IsotonicPEP

def parse_args():
    parser = argparse.ArgumentParser(
        description="Estimate PEPs for identifications using isotonic regression.",
        epilog="""
Example usage:
  1) Q2PEP:
     ./main.py q2pep --input example/peptide.target.txt --qcol q-value --pava-method ip --center-method mean --output /data

  2) OBS2PEP (a concatenated target and decoy input file):
     ./main.py obs2pep --cat-file example/peptide.cat.txt --score-col score --type-col type --target_label 0 --decoy_label 1 --output /data

  3) OBS2PEP (separate target and decoy input files):
     ./main.py obs2pep --target-file example/peptide.target.txt --decoy-file example/peptide.decoy.txt --score-col score --output /data
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Subparsers for method selection (q2pep or obs2pep)
    subparsers = parser.add_subparsers(dest="method", required=True,
                                       help="PEP estimation method: 'q2pep' or 'obs2pep'")
    
    # q2pep parser
    parser_q = subparsers.add_parser("q2pep", help="Estimate PEP from q-values (q2pep method).")
    parser_q.add_argument("--input", type=str, required=True,
                          help="Path to the TSV file containing q-values.")
    parser_q.add_argument("--qcol", type=str, default="q-value",
                          help="Column name for q-values in the input file (default: 'q-value').")
    parser_q.add_argument("--smooth", action="store_true",
                          help="Apply block-merge pre-processing (default: False).")
    parser_q.add_argument("--pava-method", type=str, choices=["basic", "ip"], default="basic",
                          help="PAVA method to use: 'basic' for PAVA regression or 'ip' for PAVA interpolation (default: basic).")
    parser_q.add_argument("--center-method", type=str, choices=["mean", "median"], default="mean",
                          help="Center method to use in PAVA interpolation (default: mean).")
    parser_q.add_argument("--output", type=str, required=True,
                          help="Output file path (or output directory; if directory, a default name 'outputPEP.target.txt' will be used).")
    
    # obs2pep parser
    parser_obs = subparsers.add_parser("obs2pep", help="Estimate PEP from target-decoy observations (obs2pep method).")
    group_obs = parser_obs.add_mutually_exclusive_group(required=True)
    group_obs.add_argument("--cat-file", type=str,
                           help="Path to a concatenated TSV file containing score and label columns.")
    group_obs.add_argument("--target-file", type=str,
                           help="Path to the TSV file containing target scores (used in separate input mode).")
    parser_obs.add_argument("--decoy-file", type=str,
                           help="Path to the TSV file containing decoy scores (required in separate input mode).")
    parser_obs.add_argument("--score-col", type=str, default="score",
                           help="Column name for score (default: 'score').")
    parser_obs.add_argument("--type-col", type=str, default="label",
                           help="Column name for target/decoy (default: 'label').")
    parser_obs.add_argument("--target-label", type=str, default="target",
                           help="Target label in concatenated input file (default: 'target').")
    parser_obs.add_argument("--decoy-label", type=str, default="decoy",
                           help="Decoy label in concatenated input file (default: 'decoy').")
    parser_obs.add_argument("--output", type=str, required=True,
                           help="Output file path or directory for target results. (Note: only target PEPs are saved.)")
    
    # Global argument to optionally disable q-value calculation.
    parser.add_argument("--no-calc-q", action="store_false", dest="calc_q",
                        help="Do not estimate q-values from calculated PEPs (default: calculate q-values).")
    
    parser.add_argument("--verbose", action="store_true",
                        help="Print parameter information.")
    
    return parser.parse_args()

def main():
    args = parse_args()

    if args.verbose:
        print("Parameters:")
        for arg, value in sorted(vars(args).items()):
            print(f"  {arg}: {value}")
    
    # For q2pep, pava_method and center_method are defined; for obs2pep, use defaults.
    pava_method = args.pava_method if hasattr(args, "pava_method") else "basic"
    center_method = args.center_method if hasattr(args, "center_method") else "mean"
    pep_processor = IsotonicPEP(pava_method=pava_method, center_method=center_method)
    
    method = args.method.lower()
    if method == "q2pep":
        # Load q-values TSV file.
        try:
            df = pd.read_csv(args.input, sep="\t")
        except Exception as e:
            sys.exit(f"Error reading input file: {e}")

        if args.qcol not in df.columns:
            sys.exit(f"Column '{args.qcol}' not found in input file.")

        q_values = df[args.qcol].values
        if args.calc_q:
            pep_series, q_series = pep_processor.pep_regression(
                q_values=q_values,
                method="q2pep",
                calc_q=args.calc_q,
                smooth=args.smooth,
                pava_method=args.pava_method,
                center_method=args.center_method
            )
            # Add the computed PEP and estimated q-values to the DataFrame.
            df["pep"] = pep_series
            df["est_q"] = q_series
        else:
            pep_series = pep_processor.pep_regression(
                q_values=q_values,
                method="q2pep",
                calc_q=args.calc_q,
                smooth=args.smooth,
                pava_method=args.pava_method,
                center_method=args.center_method
            )
            df["pep"] = pep_series
        # Determine output filename.
        output_path = args.output
        if os.path.isdir(output_path):
            output_path = os.path.join(output_path, "outputPEP.target.txt")
        try:
            df.to_csv(output_path, sep="\t", index=False)
            print(f"Output saved to: {output_path}")
        except Exception as e:
            sys.exit(f"Error writing output file: {e}")
    
    elif method == "obs2pep":
        # Determine input mode: concatenated or separate.
        if args.cat_file:
            # Concatenated mode: load single TSV file.
            try:
                df_obs = pd.read_csv(args.cat_file, sep="\t")
            except Exception as e:
                sys.exit(f"Error reading input concatenated file: {e}")

            # Check if required columns exist.
            if args.score_col not in df_obs.columns or args.type_col not in df_obs.columns:
                sys.exit(f"Columns '{args.score_col}' and/or '{args.type_col}' not found in input file.")
            # Filter targets and decoys.
            df_target = df_obs[df_obs[args.type_col] == args.target_label]
            df_decoy = df_obs[df_obs[args.type_col] == args.decoy_label]
            # Call pep_regression using the concatenated observations.
            if args.calc_q:
                pep_series, q_series = pep_processor.pep_regression(
                    obs=df_obs[[args.score_col, args.type_col]].values,
                    method="obs2pep",
                    calc_q=args.calc_q,
                    target_label=args.target_label,
                    decoy_label=args.decoy_label
                )
                df_target["pep"] = pep_series
                df_target["est_q"] = q_series
            else:
                pep_series = pep_processor.pep_regression(
                    obs=df_obs[[args.score_col, args.type_col]].values,
                    method="obs2pep",
                    calc_q=args.calc_q,
                    target_label=args.target_label,
                    decoy_label=args.decoy_label
                )
                df_target["pep"] = pep_series
            # Determine output filename.
            output_path = args.output
            if os.path.isdir(output_path):
                output_path = os.path.join(output_path, "outputPEP.target.txt")
            try:
                df_target.to_csv(output_path, sep="\t", index=False)
                print(f"Output saved to: {output_path}")
            except Exception as e:
                sys.exit(f"Error writing output file: {e}")
        else:
            # Separate mode: target and decoy files are provided.
            if not (args.target_file and args.decoy_file):
                sys.exit("For separate input mode, both --target-file and --decoy-file must be provided.")
            try:
                df_target = pd.read_csv(args.target_file, sep="\t")
                df_decoy = pd.read_csv(args.decoy_file, sep="\t")
            except Exception as e:
                sys.exit(f"Error reading input target/decoy file(s): {e}")
            
            # Check if the specified score column exists in both files.
            if args.score_col not in df_target.columns:
                sys.exit(f"Column '{args.score_col}' not found in target file.")
            if args.score_col not in df_decoy.columns:
                sys.exit(f"Column '{args.score_col}' not found in decoy file.")

            # Call pep_regression using separate target and decoy inputs.
            if args.calc_q:
                pep_series, q_series = pep_processor.pep_regression(
                    target=df_target[args.score_col].values,
                    decoy=df_decoy[args.score_col].values,
                    method="obs2pep",
                    calc_q=args.calc_q
                )
                df_target["pep"] = pep_series
                df_target["est_q"] = q_series
            else:
                pep_series = pep_processor.pep_regression(
                    target=df_target[args.score_col].values,
                    decoy=df_decoy[args.score_col].values,
                    method="obs2pep",
                    calc_q=args.calc_q
                )
                df_target["pep"] = pep_series
            # Create output directory if it does not exist.
            output_dir = args.output
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            target_out = os.path.join(output_dir, "outputPEP.target.txt")
            try:
                df_target.to_csv(target_out, sep="\t", index=False)
                print(f"Target output saved to: {target_out}")
            except Exception as e:
                sys.exit(f"Error writing output files: {e}")
    else:
        sys.exit("Unknown method. Use 'q2pep' or 'obs2pep'.")

if __name__ == "__main__":
    main()
