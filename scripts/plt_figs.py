#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# module load buildtool-easybuild/4.9.4-hpc71cbb0050 GCCcore/12.3.0
# module load Python/3.11.3
# source /proj/proteoforma_nsc/11analysis/envs/annotate_peptide_env/bin/activate

def estimate_q(df):
    df.index += 1
    df["est_q"] = df["posterior_error_prob"].cumsum() / df.index
    return df

def estimate_and_save(output_dir, N_values):
    """
    For each value in N_values, reads the peptide.target.txt file,
    estimates the q-values and saves the output as q_pep.txt.
    """
    for n in N_values:
        # Build the input file path
        input_file = os.path.join(output_dir, str(n), "peptide.target.txt")
        df = pd.read_csv(input_file, sep="\t")
        df = estimate_q(df)
        # Drop unnecessary columns
        columns_to_drop = ["PSMId", "filename", "peptide", "proteinIds"]
        df = df.drop(columns=columns_to_drop)
        # Build the output file path
        output_file = os.path.join(output_dir, str(n), "q_pep.txt")
        df.to_csv(output_file, sep='\t', index=False)

def cal_max_rel_diff(q_value, q_est):
    """
    Calculates the maximum relative difference (in percent) between q_value and q_est.
    """
    relative_difference = np.abs(q_est - q_value) / np.where(q_value != 0, q_value, np.nan)
    return np.nanmax(relative_difference) * 100

def main():
    ### Part 1: Estimate and Save q-values for each dataset, run, and model
    datasets = [
        "PXD003868", "PXD004325", "PXD004424", "PXD004467", "PXD004536",
        "PXD004565", "PXD004947", "PXD004948", "PXD005025", "PXD013274"
    ]
    num_runs = 10
    models = ["irls", "ispline.rank", "ispline.score", "pava.rank", "pava.score"]
    N_values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 20, 40]
    
    # Modify this base path to match where your data resides on the HPC system
    base_output_dir = "/proj/proteoforma_nsc/smooth_q_to_pep/perc_test"
    
    for dataset in datasets:
        output_dir = os.path.join(base_output_dir, dataset)
        for run in range(1, num_runs + 1):
            for model in models:
                print(f"Processing: run {run}, dataset {dataset}, model {model}")
                estimate_and_save(os.path.join(output_dir, f"run{run}", model), N_values)
    
    ### Part 2: Calculate metrics and generate plots
    # Directory where the input files for plotting are located
    base_dir = "/proj/proteoforma_nsc/smooth_q_to_pep/perc_test"
    # Define two groups of models (for rank-based and score-based comparisons)
    groupA = ["irls", "pava.rank", "ispline.rank"]
    groupB = ["irls", "pava.score", "ispline.score"]
    
    # Mapping of model names to labels for plotting
    label_map = {
        "irls": r"IRLS",
        "pava.rank": r"PAVA.rank",
        "ispline.rank": r"ISpline.rank",
        "pava.score": r"PAVA.score",
        "ispline.score": r"ISpline.score"
    }
    
    # Mapping of model names to colors
    color_map = {
        "irls": "tab:blue",
        "pava.rank": "tab:orange",
        "ispline.rank": "tab:green",
        "pava.score": "tab:orange",
        "ispline.score": "tab:green"
    }
    
    # Directory to save the generated figures; create if it doesn't exist
    save_dir = "/proj/proteoforma_nsc/smooth_q_to_pep/perc_test/figs"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    for dataset in datasets:
        output_dir = os.path.join(base_dir, dataset)
        for N in N_values:
            # Loop through the two groups
            for group_label, suffix, group_models in zip(
                ["rank-based methods", "score-based methods"],
                ["rank", "score"],
                [groupA, groupB]
            ):
                plt.figure(figsize=(8, 8), dpi=600)
                # Add text annotations for N and header info
                plt.text(0.05, 0.9, f"N = {N}", transform=plt.gca().transAxes,
                         fontsize=18, fontweight='bold')
                plt.text(0.05, 0.83, f"max relative difference:",
                         transform=plt.gca().transAxes, fontsize=18, fontweight='bold')
                cnt = 0
                # Process each model in the current group
                for model in group_models:
                    color = color_map.get(model, "black")
                    model_runs = []
                    # Concatenate results from all runs
                    for run in range(1, num_runs + 1):
                        file_name = "q_pep.txt"
                        file_path = os.path.join(output_dir, f"run{run}", model, str(N), file_name)
                        # Read only the required columns
                        df_tmp = pd.read_csv(file_path, sep="\t", usecols=["q-value", "est_q"])
                        df_tmp["run"] = run
                        model_runs.append(df_tmp)
                    df_model = pd.concat(model_runs, ignore_index=True)
                    
                    # Collect all unique q-value points
                    unique_q = np.unique(df_model["q-value"].values)
                    unique_q.sort()
                    aligned_df = pd.DataFrame({"q-value": unique_q})
                    
                    # Align the step function for each run onto the shared grid
                    for run in range(1, num_runs + 1):
                        df_run = df_model[df_model["run"] == run].sort_values("q-value")
                        run_qvals = df_run["q-value"].values
                        run_qests = df_run["est_q"].values
                        aligned_values = []
                        current_index = 0
                        current_qest = run_qests[0]
                        for q in unique_q:
                            while current_index < len(run_qvals) and q >= run_qvals[current_index]:
                                current_qest = run_qests[current_index]
                                current_index += 1
                            aligned_values.append(current_qest)
                        aligned_df[f"run_{run}"] = aligned_values
                    
                    # Compute the mean and standard deviation across runs at each q-value
                    run_cols = [col for col in aligned_df.columns if col.startswith("run_")]
                    aligned_df["mean"] = aligned_df[run_cols].mean(axis=1)
                    aligned_df["std"] = aligned_df[run_cols].std(axis=1)
                    
                    # Compute the maximum relative differences across runs
                    max_diffs = []
                    for run in range(1, num_runs + 1):
                        diff = cal_max_rel_diff(aligned_df["q-value"].values,
                                                aligned_df[f"run_{run}"].values)
                        max_diffs.append(diff)
                    avg_max_diff = np.mean(max_diffs)
                    std_max_diff = np.std(max_diffs)
                    
                    # Print the average max relative difference on the figure
                    plt.text(0.05, 0.76 - cnt * 0.05,
                             f"{label_map.get(model, model)}: {avg_max_diff:.0f} ± {std_max_diff:.0f}%",
                             transform=plt.gca().transAxes, fontsize=18)
                    cnt += 1
                    
                    # Plot the mean step curve and fill the area between mean ± std
                    plt.step(aligned_df["q-value"], aligned_df["mean"],
                             where="post", color=color, linewidth=2,
                             label=label_map.get(model, model))
                    y_lower = aligned_df["mean"] - aligned_df["std"]
                    y_upper = aligned_df["mean"] + aligned_df["std"]
                    plt.fill_between(aligned_df["q-value"], y_lower, y_upper,
                                     step="post", color=color, alpha=0.2)
                
                # Plot reference lines
                a = np.linspace(0, 1, 100)
                b = a / 10 ** 0.25
                c = a * 10 ** 0.25
                plt.plot(a, b, c="k", linewidth=0.5, linestyle="--")
                plt.plot(a, c, c="k", linewidth=0.5, linestyle="--")
                plt.plot(a, a, c="k", linewidth=0.5)
                
                # Set title, axis labels, legend and scales
                plt.title(f"{group_label}: {dataset}", fontsize=24)
                plt.xlabel("FDR-derived $q$ values", fontsize=24)
                plt.ylabel("PEPs-derived $q$ values", fontsize=24)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.legend(loc='lower right', fontsize=24)
                plt.xscale("log")
                plt.yscale("log")
                plt.xlim(1e-5, 1)
                plt.ylim(1e-5, 1)
                
                # Save the figure
                fig_filename = f"{dataset}_N{N}_{suffix}.png"
                plt.savefig(os.path.join(save_dir, fig_filename), bbox_inches="tight")
                plt.close()

if __name__ == "__main__":
    main()
