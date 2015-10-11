[![Build Status](https://travis-ci.org/gyper/gyper.svg)](https://travis-ci.org/gyper/gyper)

Gyper (Graph genotYPER) is a genotyper for aligned DNA sequencing data. The input files are BAM files that have been sorted and indexed (e.g. using samtools). Gyper creates a partial order graph which represent haplotypes. Currently only six HLA genes are supported: HLA-A, HLA-B, HLA-C, HLA-DQA1, HLA-DQB1, and HLA-DRB1. Gyper aligns reads from certain position of the input BAM file to the graphs and determines the individual's haplotype in a fast and accurate manner.

## Dependencies
* Boost>=1.56.0
* SeqAn>=2.1.0 (in development)
* zlib>=1.2.8

Furthermore, Gyper is released with a CMake build system which requires CMake>=2.8.

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
git clone git@github.com:gyper/gyper.git gyper
cd gyper
mkdir build
cd build
cmake ..
make
sudo make install
```
Gyper should now be installed to `/usr/local/bin/gyper`.

### User-only install
Maybe you are a pawn in a big company, and you don't have root access on your computer. In this case, worry not, because you can still compile and use Gyper. I'll assume you have or can install cmake and zlib. Follow the instructions on [how to get Boost](http://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html). Then, [install and cmake SeqAn](http://seqan.readthedocs.org/en/latest/BuildManual/UsingTheSeqAnBuildSystem.html#user-library-installation) using the `-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY` options because we only need the library.

You can specify the locations of both libraries when you use cmake, e.g.

```sh
cmake -DBOOST_INCLUDEDIR=/some/path/you/installed/boost_1_XX_0/ -DSEQAN_INCLUDE_PATH=/some/path/you/installed/seqan/include/ ..
make
```
If you're feeling adventurous, you can also set the `BOOST_INCLUDEDIR` and `SEQAN_INCLUDE_PATH` environment variables. The Gyper binary file will be located in `./bin/gyper`. If you'd like, add it to your `PATH` variable.
