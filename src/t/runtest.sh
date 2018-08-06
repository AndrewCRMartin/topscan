echo "Checking with pre-built secstr data"
topscan -b 1yqvY.ss >1yqvY.out
diff 1yqvY.out 1yqvY.ss.out.ref

echo "Checking while building secstr data"
topscan -p -b 1yqvY.pdb >1yqvY.out
diff 1yqvY.out 1yqvY.pdb.out.ref

\rm -f 1yqvY.out

