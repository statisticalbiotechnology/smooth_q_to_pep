import pandas as pd
import argparse

def revise_files(target_file, decoy_file, revised_target_file, revised_decoy_file, N):
    target_df = pd.read_csv(target_file, sep="\t")
    decoy_df = pd.read_csv(decoy_file, sep="\t")

    # Identify top N rows from the target file by highest 'score'
    top_target = target_df.nlargest(N, 'score').copy()

    # Remove these top entries from the target DataFrame
    target_df = target_df.drop(top_target.index)

    # Update the PSMId of the top entries to change "target" to "decoy"
    # top_target['PSMId'] = top_target['PSMId'].str.replace('^target', 'decoy', regex=True)

    # Append the top target entries to the decoy DataFrame
    decoy_df = pd.concat([top_target, decoy_df], ignore_index=True)

    # Save the revised DataFrames to the provided output file paths
    target_df.to_csv(revised_target_file, sep="\t", index=False)
    decoy_df.to_csv(revised_decoy_file, sep="\t", index=False)

    print(f"[INFO]Revised files saved as {revised_target_file} and {revised_decoy_file} :)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Move top N scoring PSM entries from target to decoy file.")
    parser.add_argument("--target", required=True, help="Path to the peptide.target.txt file.")
    parser.add_argument("--decoy", required=True, help="Path to the peptide.decoy.txt file.")
    parser.add_argument("--revised_target", required=True, help="Path to save the revised peptide.target.txt file.")
    parser.add_argument("--revised_decoy", required=True, help="Path to save the revised peptide.decoy.txt file.")
    parser.add_argument("--N", type=int, required=True, help="Number of top-scoring PSMs to move from target to decoy.")
    
    args = parser.parse_args()
    revise_files(
        target_file=args.target,
        decoy_file=args.decoy,
        revised_target_file=args.revised_target,
        revised_decoy_file=args.revised_decoy,
        N=args.N
    )