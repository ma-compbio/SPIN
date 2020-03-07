# SPIN
SPIN - Spatial Position Inference of the Nuclear genome.

## Required Package
- Python (>=3.6)
- scikit-learn
- NumPy 
- SciPy
- pandas
- pickle

## Running command
python main.py [Options]

The options:

- -i : 1D genomic measurements of nuclear organization

- --hic : Hi-C interactions

- -w : Resolution

- -n : Number of states

- -m : Choose mode. Supported: full, hic

- -o : Output dir

- -g : Genome bin file

- --prev : reload existing model
