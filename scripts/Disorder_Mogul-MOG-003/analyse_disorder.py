#! /usr/bin/env python
###############################################################################
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-09-29: created by Andrew J. Peel, Cambridge Crystallographic Data Centre
#
###############################################################################

from ccdc.io import EntryReader
from ccdc.conformer import GeometryAnalyser
import matplotlib.pyplot as plt

# Read CSD entry
entry_reader = EntryReader('CSD')
purzio_entry = entry_reader.entry('PURZIO')
print(f'Refcode: {purzio_entry.identifier}')

# Get crystal
purzio_crystal = purzio_entry.crystal

# Examine disorder groups
print(f'Crystal has disorder: {purzio_crystal.has_disorder}')
print(f'Number of disorder assemblies: '
      f'{len(purzio_crystal.disorder.assemblies)}')
for i, assembly in enumerate(purzio_crystal.disorder.assemblies, start=1):
    print(f'Assembly {i} has {len(assembly.groups)} disorder groups.')

# Activate disorder assembly, group 1 of 2
first_assembly = purzio_crystal.disorder.assemblies[0]
first_group = first_assembly.groups[0]
first_group.activate()

print(f'Activated assembly_id: {first_assembly.id}')
print(f'Activated group_id: {first_assembly.active.id}')
print(f'Occupancy of disorder group: {first_group.occupancy:.3f}')

# Examine the molecules
mol = purzio_crystal.molecule
print(f'Number of molecules in crystal: '
      f'{len( mol.components)}')

# Select component by index
comp_index = 0

# Examine the first molecule
comp = mol.components[comp_index]
print(f'Analysing molecule {comp_index + 1} of '
      f'{len(mol.components)}')

# Check if this molecule is disordered
if any([atom.occupancy < 1 for atom in comp.atoms]):
    print('This molecule is disordered.')
else:
    print('This molecule is not disordered.')

# Identify disordered atoms
print([(atom.label, atom.occupancy) for atom in comp.atoms
      if atom.occupancy < 1])

# Create an instance of geometry analyser
engine = GeometryAnalyser()

# Adjust settings - exclude organometallics and powder structures
engine.settings.organometallic_filter = 'organics_only'
engine.settings.powder_filter = True

# Analyse geometry
geom_analysed_mol = engine.analyse_molecule(comp)

# For example, look at torsion measurements
for torsion in geom_analysed_mol.analysed_torsions:
    if torsion.unusual and torsion.enough_hits:
        print(f'Unusual torsion: {",".join(torsion.atom_labels)}  '
              f'Torsion: {torsion.value:.3f} °')

# Accessing histogram data
# Take example of an unusual torsion
unusual_tors = [t for t in geom_analysed_mol.analysed_torsions
                if t.unusual and t.enough_hits]

# Retrieves the first occurance of specified torsion
tor_query = next(t for t in unusual_tors if
                 t.atom_labels == ['C39', 'C28', 'C29', 'C30'])

counts = tor_query.histogram()

# Create histogram
n_bins = len(counts)
bin_min, bin_max = 0, 180
bin_width = (bin_max - bin_min) / n_bins
bin_edges = [bin_min + i * bin_width for i in range(n_bins + 1)]

plt.figure(figsize=(8, 5))
plt.bar(
    bin_edges[:-1],
    counts,
    width=bin_width,
    edgecolor="black",
    align="edge",
    color="limegreen"
)

# Red line for torsion value
torsion_val = tor_query.value
plt.axvline(
    torsion_val,
    color="red",
    linestyle="--",
    linewidth=2
)

# Add text label above the histogram, slightly shifted
y_max = max(counts)
plt.text(
    torsion_val + 2,     # shift slightly right
    y_max * 0.9,         # place near the top
    f"{torsion_val:.2f}°",
    color="red",
    fontsize=10,
    rotation=90,
    va="bottom"
)

plt.xlabel("Torsion Angle (degrees)")
plt.ylabel("Count of Observation Hits")
plt.title(f"Histogram for torsion {tor_query.atom_labels}")
plt.xlim(bin_min, bin_max)

plt.show()

print('\n----------------------------------\n')
print('Now using disorder combinations\n')

# Interate over disorder combinations
for comb_id, combination in enumerate(purzio_crystal.disorder.combinations):
    for comp_id, component in enumerate(purzio_crystal.molecule.components):
        print(f'Combination_id {comb_id}, Component_id {comp_id}')
        analysed_mol = engine.analyse_molecule(component)
        unusual_tors = [t for t in analysed_mol.analysed_torsions
                        if t.enough_hits and t.unusual]
        if len(unusual_tors) > 0:
            for t in unusual_tors:
                print(f'Unusual torsion: {",".join(t.atom_labels)}'
                      f'  Torsion: {t.value:.3f} °')
        else:
            print('No unusual torsions with enough hits.')
        print()
