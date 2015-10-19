[![Build Status](https://travis-ci.org/gyper/gyper.svg)](https://travis-ci.org/gyper/gyper)

Gyper (Graph genotYPER) is a genotyper for aligned DNA sequencing data. The input files are BAM files that have been sorted and indexed (e.g. using samtools). Gyper creates a partial order graph which represent haplotypes. Currently only six HLA genes are supported: HLA-A, HLA-B, HLA-C, HLA-DQA1, HLA-DQB1, and HLA-DRB1. Gyper aligns reads from certain position of the input BAM file to the graphs and determines the individual's haplotype in a fast and accurate manner.

## Dependencies
* Boost>=1.56.0
* SeqAn>=2.1.0 (in development)
* zlib>=1.2.8

Furthermore, Gyper is released with a CMake build system which requires CMake>=2.8.

## Simple installation (unix-like systems)
### Install dependencies
In case you don't have CMake or zlib installed, they should be available in any package manager, e.g.

* apt-get (Ubuntu/Debian): `sudo apt-get install cmake zlib`
* pacman (Arch Linux): `sudo pacman -S cmake zlib`
* RPM (Fedora/RHEL): `sudo yum install cmake zlib`

Boost and SeqAn do not need to be installed, Gyper's cmake build system will automatically fetch them for you if they are missing.

### Install Gyper
Here we assume Gyper is cloned to ~/git/gyper and built in ~/git/gyper/build. 

```sh
cd ~/git
git clone git@github.com:gyper/gyper.git gyper
cd gyper
mkdir build && cd build
cmake ..
make
sudo make install
```
The last command installs Gyper to `/usr/local/bin/gyper`. If you don't have root access, you can add `~/git/gyper/build/bin/gyper` to your `PATH` environment variable manually instead.

### Use specific Boost or SeqAn libraries.
You can specify the locations of both libraries when you use cmake, e.g.

```sh
cmake -DBOOST_INCLUDEDIR=/some/path/you/installed/boost_1_XX_0/ -DSEQAN_INCLUDE_PATH=/some/path/you/installed/seqan/include/ ..
make
```
If you're feeling adventurous, you can also set the `BOOST_INCLUDEDIR` and `SEQAN_INCLUDE_PATH` environment variables.
