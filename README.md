# Gyper
Gyper (Graph genotYPER) is a genotyper for aligned DNA sequencing data. The input files are BAM files that have been sorted and indexed (e.g. using samtools). Gyper creates a partial order graph which represent haplotypes. Currently only six HLA genes are supported: HLA-A, HLA-B, HLA-C, HLA-DQA1, HLA-DQB1, and HLA-DRB1. Gyper aligns reads from certain position of the input BAM file to the graphs and determines the individual's haplotype in a fast and accurate manner.

## Dependencies
* Boost>=1.56.0
* SeqAn>=2.1.0 (in development)
* zlib>=1.2.8

Furthermore, Gyper is released with a CMake build system which requires CMake>=3.0.

## Simple installation (unix-like systems)
### Install dependencies
Install CMake, zlib, and Boost with a package manager of choice. E.g.

* apt-get (Ubuntu/Debian): `sudo apt-get install cmake zlib libboost-all-dev`
* pacman (Arch Linux): `sudo pacman -S cmake zlib boost`
* RPM (Fedora/RHEL): `sudo yum install cmake zlib boost`

Installing SeqAn is a little tougher, since we need their development branch on Github. Here, I assume SeqAn is cloned to ~/git/seqan and built to ~/git/seqan-build.

```sh
mkdir -p ~/git
cd ~/git
git clone git@github.com:seqan/seqan.git seqan
cd seqan
git checkout develop
cd ~/git
mkdir seqan-build
cd seqan-build
cmake ../seqan -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY
make dox
sudo make install
```

### Install Gyper
Here we assume Gyper is cloned to ~/git/gyper and built in ~/git/gyper/build. 

```sh
cd ~/git
git clone git@github.com:hannespetur/gyper.git gyper
cd gyper
mkdir build
cd build
cmake ..
make
sudo make install
```

Gyper should now be installed to `/usr/local/bin/gyper`.
