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

**macOS (Homebrew):**

```bash
brew install cd-hit
```

**Ubuntu/Debian:**

```bash
sudo apt install cd-hit
```

Verify that CD-HIT is installed and available on your PATH:

```bash
cd-hit -h
```
