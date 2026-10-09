# genome-informed-proteomics
A workflow for constructing and optimising genome-informed  protein search databases for bacterial proteomics and antigen discovery.

- ## Installation

Clone the repository:

```bash
git clone https://github.com/hjortswang/genome-informed-proteomics.git
cd genome-informed-proteomics
```

Install the required Python dependencies:

```bash
pip install -r requirements.txt
```


### CD-HIT

[CD-HIT](https://github.com/weizhongli/cdhit) is required for protein sequence clustering and must be installed separately.

On macOS, CD-HIT can be compiled from source without OpenMP:

```bash
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make openmp=no
```

Verify the installation:

```bash
./cd-hit -h
```

To make `cd-hit` available from the command line and to the Python pipeline:

```bash
sudo cp cd-hit /usr/local/bin/
```

Verify that the executable is available on your PATH:

```bash
which cd-hit
cd-hit -h
```

CD-HIT must be available on the system PATH before running the clustering pipeline.
