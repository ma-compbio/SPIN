# SPIN
SPIN - Spatial Position Inference of the Nuclear genome.

## Required Packages
SPIN requires the following Python packages to be installed:
- Python (tested on version 3.6)
- `scikit-learn` (tested on version 0.22.2)
- `NumPy` 
- `SciPy`
- `pandas`
- `pickle`

[Juicer tools](https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start) is required to extract Hi-C data from [.hic files](https://github.com/aidenlab/juicer/wiki/Pre). Requires Java Runtime Engine  installed on your file system.

## Usage

After install all dependencies, run the following python command:

`python main.py -i <input_signal> --hic <hic_interactions> -w <window_size> -n <number_of_states> -o <output_path> -g <genome_bin> [--prev <previous_model>] [--save]`

The options:

- -i \<input_signal\> : 1D genomic measurements of nuclear organization. `input_signal` file should be a tab-separated text file where each line is a vector of 1D genomic measurements from TSA-Seq/DamID. 

- --hic \<hic_interactions\> : List of Hi-C interactions added as edges in the model. `hic_interactions` file should be a tab-separated text file where the first two columns are bin numbers (needs to be consistent with `genome_bin`) and the third column is edge weight (opitional).

- -w \<window_size\> : Window size of each genome bin. Default = 100000

- -n \<number_of_states\> : Number of states to estimate for SPIN. Default = 5

- -o \<output_path\> : Output path.

- -g \<genome_bin\>: Genomic coordinates of each bin. `genome_bin` file should be a tab-separated text file where the first three columns are the genomic coordinates of each bin and the fourth column is the bin number indicator. 

- --prev \<previous_model\>: (opitional) Load previously saved model.

- --save : (opitional) Save curent model to .pkl file.


Example:

`python main.py -i input_chr1.txt --hic Hi-C_chr1.txt -w 25000 -n 5 -o example_chr1 -g bin_chr1.bed --save`
