This directory contains the semi-analytic code that has been used in De Lucia
and Blaizot (2007). This code builds on the methodologies originally introduced
in Springel et al. 2001, and De Lucia, Kauffmann and White 2004.

Details about modelling of different physical processes can be found in De
Lucia, Kauffmann and White 2004, Croton et al. 2006, and De Lucia and Blaizot
2007.


Here you have some basic information: 

DIRECTORIES:
============

code: contains the semi-analytic code. 

CoolFunctions: contains the cooling functions.  Currently Sutherland & Dopita.

devel: contains code used in developing and testing L-Galaxies.  Only useful
if you are a developer.

Idl: contains some reading and plotting examples.

input: contains input files that govern the run-time behaviour of
L-Galaxies.  See INPUT DATA below.  Also contains some test data.

output: default location for data produced by the code.  May contain
some sample output.

PhotTables: contains the photometric tables.

run: contains some examples of how to run L-Galaxies via a script, or
in parallel.

INPUT DATA:
===========

The directory ./input/ contains the following files:

--- input.par contains information about the location of the merger tree files
and code output, the recipe and parameter choices, as well as some units needed
by the code.

--- desired_outputsnaps.txt contains the number of the snapshot for which you
desire to save the output. (Remember to change NOUT in the Makefile if you want
to save more than one snapshot!).

--- millennium_zlist.txt contains the expansion factors corresponding to the
snapshots of the Millennium Simulation

MAKEFILE
========

Choose your options and compile the file ("make").
This should create an executable "L-Galaxies" in this directory.

RUNNING L-Galaxies
==================

The simplest way to run L-Galaxies is via the command
> L-Galaxies input/input.par

For examples of executing via a job control script see the entries in
the ./run directory.

Peter Thomas
29-Oct-10
