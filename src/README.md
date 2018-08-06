analyse.perl            Takes an output file from topscan and rewrites with CATH codes
buildtoplib.perl        Generate a library for topscan using reps from the database
findegs.sql             SQL code to find example for a given CAT ordered by rosolution
runanalyse.sh           Runs analysis on a PDB file for e3h3,e3h4,e4h3,e4h4
topscan.c               The actual topscan program

runworst.sh             Run the worst.sh script on a specified set of .out files
worst.sh                Finds the lowest ranked example of a CAT in a .out file
findline.perl           Finds the last occurrence of a CAT in output from analyse.perl
                        Used by worst.sh

runbest.sh              Run the best.sh script on a specified set of .out files
best.sh                 Find the highest ranked example of a CAT in a .out file ignoring
                        the probe structure
findbestline.perl       Finds the first occurrence of a CAT in output from analyse.perl
                        Used by best.sh
