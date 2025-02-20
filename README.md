# cluster-copoly
Stochastic model for PLGA with random transesterification and output
file of dyads/triads/tetrads/pentads

	The stochastic model of PLG copolymerization models forward
	polymerization, depolymerization, and random transesterification for
	the copolymerization of lactide and glycolide. Before you use it for
	your own system, there are a few things to be aware of.
	1.	Lactide and glycolide are cyclic dimers, so each
	polymerization adds TWO repeat units to the chain. Likewise, each
	depolymerization reaction removes two repeat units to reform lactide
	or glycolide. Depolymerization of LG or GL dimers to form 3-methyl
	glycolide is neglected. Some aspects of the code would need to be
	changed for systems in which reactants are strictly monomers rather
	than cyclic dimers. 2.	The transesterification included in this
	version is random. Neither the identity of the chain end nor the
	ester affect the rate, other than that the selected chain end and
	ester may not be on the same chain. In this way, formation of cycles
	is neglected in this model.

The files included in the package are the following:
1.	PLGA-stochastic_blockiness.cpp 2.	polymer.cpp 3.	polymer.h
4.	7525_008kT.dat – This is an example output file we generated using
the following command line input:

start [exe_name] [run_name] 0.04332 1.9737 5.9211 20.59 4.829 10.25
3.492 0.885 1.302 0.728 1.134 0.008 2880 1000 1

Of course, you should replace the exe name and run_name with the name of
your exe and whatever you want to call your output files. If you use the
same rate constants as our example, you should be able to produce
something quite similar to the example output file. Because this is a
stochastic model, you should not be able to reproduce this file exactly.
However, plots of conversion vs. time from the example output should
look the same as plots generated from your output.

5.	7525_008.seq – This is an example output file we generated using
the same command line input as in (4). Again, you should not be able to
reproduce this file exactly but should get the same format and plots of
dyad/triad/tetrad populations over time should look quite similar to the
example output when the same command line input is used.

A few other notes:

1.	The command line input is ordered as follows:

[Exe_name] [run_name] [initiator_conc.] [A_conc.] [B_conc.] [kAA] [kAB]
[kBA] [kBB] [kA-A] [kA-B] [kB-A] [kB-B] [kT] [termination_time(min)]
[number chains] [catalyst ratio]

All concentrations are in mol/L, and rate constants are in corresponding
units (either s-1 or Lmol-1s-1). These are converted to number of
molecules and corresponding units and rate constants inside the code,
but converted back to mol/L for output. For our case, the difference
between the rate constants that we input and what Gillespie calls
reaction constants (all in units of s-1 is basically a factor of volume
and avogadros number. For some reactions there might be more of a
difference, for example a bimolecular reaction involving two of the same
reactant would also involve a factor of 2. Refer to Gillespie’s work
(Gillespie, D. T. A General Method for Numerically Simulating the
Stochastic Time Evolution of Coupled Chemical Reactions. Journal of
Computational Physics 1976, 22 (4), 403–434.
https://doi.org/10.1016/0021-9991(76)90041-3.) to make sure these are
converted correctly.

2.	There are several possible output files. Each output file will be
named with whatever name you give the run (which we will call “run_name”
for example purposes) and a file extension designating what type of data
is included in the file. These include: a.	run_name.dat – General
kinetic data. Time, concentrations of I, A, B, PA, and PB  (in mol/L),
overall conversion. PA and PB are polymers terminated with A or B,
respectively. Currently, this file contains 1000 data points, which are
output at even intervals throughout the reaction time. The user may
easily modify the number of data points by changing time_report_int in
line 148. b.	run_name.seq – Polymer sequence data. This reports
fractions of each of the possible dyads, triads, tetrads, and pentads.
For example, AA here corresponds to the number of AA dyads divided by
the total number of AA + AB + BA + BB dyads and similar for longer
sequence motifs. Currently, this file contains only 2 data points which
are output after the first step and end of the simulation. The number of
data points may be modified by changing time_report_int2 in line 149.
This data takes some time to generate so I wouldn’t recommend outputting
very often if you don’t need to. c.	run_name.mwd – currently this
output file is commented out. You can add it back in at line 1018. This
report returns a molecular weight distribution for the polymers at the
end of the reaction. d.	run_name.composition – currently this output
file is commented out. You can add it back in at line 1019. This report
returns the composition of all the chains in your model system.
e.	run_name.matrix – currently this output file is commented out. You
can add it back in at line 1020. This report returns the sequence of
every chain in order. Each line represents one polymer and each repeat
unit is separated by a comma. “A” repeat units are designated by a “1”
and “B” repeat units are designated with a “2”. f.	run_name.maldi –
currently this output file is commented out. You can add it back in at
line 1021. This report returns the number of A and B repeat units in
each chain. It is structured as #A’s (tab) #B’s (tab) # polymers. So if
you have 3 polymers which have 1 A repeat unit and 2 B repeat units,
that line would read: 1	2	3 Combinations which do not exist in
your model data are not reported. We used this data in combination with
Isospec (IsoSpec: Hyperfast Fine Structure Calculator Mateusz K. Łącki,
Michał Startek, Dirk Valkenborg, and Anna Gambin Analytical Chemistry
2017 89 (6), 3272-3277 DOI: 10.1021/acs.analchem.6b01459) to generate
model MALDI spectra.

3.	Some data types are more sensitive to the number of chains than
others. For example, conversion vs. time data requires really very few
chains to get data that is qualitatively the same every time. Model
MALDI data on the other hand, requires a lot more chains to be
consistent, potentially more than is really feasible for running the
code on your personal computer.
