# Script to make a xyz file

python symmetry_polar.py
python nt_to_flat.py

n=10
motif=20

touch atoms
for (( unit=0; unit<$n; unit++))
do
    cat atoms.dat >> atoms
done

((natoms=$n*$motif))

rm flat_cell.xyz
touch flat_cell.xyz
echo "$natoms
" >> flat_cell.xyz

paste atoms flat_cell >> flat_cell.xyz

rm atoms flat_cell cell_polar
