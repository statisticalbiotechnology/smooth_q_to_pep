main_dir: pyiso_test
    - each dataset has a folder
        - in each dataset's folder, there are 10 runs' folders
            - in each run's folder, there are six folders for six models:
                percolator, qPAVA, qipPAVA, qISpline, dPAVA, dISpline,
                a crux folder,
                and, a interm folder to store files with flipped labels
                - in each model's folder, there is one folder for every knocking-out target PSMs
            - the run_crux.sh script run crux 10 times with 10 seeds, and save each run's result in the crux folder under each run's directory


in each run, execute run_crux.sh first using different seeds
read crux's percolator's output make-pin.pin file
run qIRLS percolator using the make-pin.pin file, and get the weights.pin, features.pin, target.peptide.txt, and decoy.peptide.txt
insert weights to features.pin
revise labels and files
save them into interm directory

run percolator using N=1, 2, 3, ... knocking-out target PSMs and save the results in percolator/1, 2, 3, ... using the same seed

run the other five models using N=0, 1, 2, 3, ... knocking-out target PSMs and save the results in model/0, 1, 2, 3, ...

data_dir: datasets
    - each dataset has a folder
        - in each folder, there is a raw file, an mzml file, a fasta file, and a run_crux.sh script

use sbatch, each run of each dataset is excuted using 1 job array
thus there are in total 10x10 nodes requested
