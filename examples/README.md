
# Gmxstrain example: Calculation of elastic constants of a monoclinic diphenyl-sulfone crystal

1. Starting with a crystal structure (usually a cif file from the [CCDC](https://www.ccdc.cam.ac.uk/),
here dps.cif) generate a super-cell (multiples of the unit cell in 3 dimensions) of a size that ensures that
no atom is closer than twice the Lennard Jones cut-off paramter (usually 1-1.2 Angstrom) plus a small
buffer of at least one Angstrom, to any of is periodic replicas. Here this requires a supercell comprising
3x4x3 unit cells. We use chimera (after converting the cif to PDB format), or gdis for this purpose. Save
the super cell in PDB format, Here this file is called dps.pdb. Make sure that each of the molecules in
the structure file has canonical atom ordering (openbabel can be used for this purpose)

2. Extract one of the molecules in dps.pdb, save as PDB file (dps1.pdb), convert to mol2
format (dps.mol2, e.g., with openbabel),
```
obabel dps1.pdb -Odps.mol2
```
and replace the partial charges in the mol2 file by ESPD charges. We generally use
AM1BCC, RESP, or ABCG2 charges (see [Ambertools](http://ambermd.org/AmberTools.php).
In the given example we use [ABCG2](https://clickff.org/amberweb/Antechamber/Create)

3. To generate a Gromacs topology input file for the molecule with GAFF2 parameters, we use acpype, 
as in:
```
acpype -i dps.mol2 -c user -a gaff2
```
This should create a number of files including dps_GMX.itp and dps_GMX.top which contain force field
parameters based on the GAFF2 parameter set, including initial values for LJ parameters for the molecule
in Gromacs topology format. Combine the two files, and, at the end, include a molecules section
with the number of molecules in your super-cell:
```
[ molecules ]
 DPS 144  
```

3. Put all the required input files in a folder and start the run

As input files you need the topology file (dps.top), and the structure (dps.pdb) that were
generated above. Other input files incldue two gromacs mdp parameter files (mdnvt.mdp, mdnpt.mdp),
which should remain un-changed for all cases, apart from the temperate of the thermostate (parameter ref_t)
which you might want to change, according to your preferences. A few text files used to generate STDIN
for various tasks (0, 1, 1.ndx, 4stress), generally remain unchanged, and are also included there.
Now start the script on the command line:
```
./gmxstrain.awk -v id=dps -v nmax=3 -v dd=0.004 -v t0=1000 1 > er-dps 2>&1 &
```
On a state of the art workstation with a decent GPU this should take about a couple of hours to finish.

4. Analyze

The gmxstrain.awk command generates numberous output files, which can be used for de-bugging purposes in
case anything goes wrong (this probably requires some experinece with Gromacs). If the run finishes
gracefully, the main output is in nmax files, here called ec-dps-004000, ec-dps-008000, ec-dps-012000.
Each contains elastic constants as obtained with a given stress. Once you made sure that are in the
linear regime, you can average these results using the accessory script avgec-column.awk:
```
./avgec-column.awk ec-dps-* > avg-ec-dps
```
Finally you can also calculate elastic moduli and anisotropy using another script:
```
./anisotropy3.awk avg-ce-dps
avg-ce-dps: BV     7.97 SV:     3.86 BR     7.47 SR:     2.44 BH     7.72 SH:     3.15 UA:     1.03
```
This provides the bulk modulus (B), the shear modulus (S), each as Voigt (V), Reuss (R), or Hill average,
in units of GPa, as well as the (universal) elastic anisotroy (UA).

