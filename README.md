topscan V2.0
============

`topscan` is a program for generating a topology string for a protein
and for scanning the topology string of a protein against a library of
such strings to find similar folds.

The topology string represents the secondary structure elements
together with directional information.

The method is described in: Martin, A.C.R. (2000) 'Ups and downs of
protein topology; rapid comparison of protein structure.' Protein
Engineering 13:829-837. [PMID: 11239082]

Prerequisites
-------------

`topscan` requires that you have a program installed for secondary
structure calculation. This can be:

- `pdbsecstr` (http://www.bioinf.org.uk/software/bioptools/)
- `STRIDE` (http://webclu.bio.wzw.tum.de/stride/), or
- `DSSP` (https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html)

The default is to use `pdbsecstr`.

Installation
------------

You need to install your chosen secondary structure calculation
program in your path (e.g. in `~/bin`)

To install `topscan` simply enter the source directory and run make:

```
cd src
make
```
to compile.

Once the program has compiled, install by typing:

```
make install
```
to copy the executables to your `~/bin` directory and data files to
`~/data`.

Test the install by typing:

```
make test
```
(Note that this must be done after the install.)

Getting Help
------------

To get help, type:

```
topscan -h
```
