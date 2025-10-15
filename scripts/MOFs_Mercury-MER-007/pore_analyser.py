from ccdc.io import EntryReader
from ccdc.descriptors import CrystalDescriptors
from pathlib import Path
import pandas as pd


# Read CSD entries
entry_reader = EntryReader('CSD')

# Selection of MOF-76 type structures M = lanthanoid, yttrium
refcodes = ['MARXEK', 'NADZID', 'QOTZEG', 'SADLIU',
            'SEHXIN', 'UPIBOM', 'YIMPAP', 'YIMPIX', 'YIMSAS']

# Headers for table
print(f"{'Refcode':<10} {'Formula':<20} {'He Volume (Å³)':>15} {'System Vol (Å³)':>18}")
print("-" * 65)

# Calculate pore properties
results = []
for refcode in refcodes:
    mof = entry_reader.entry(refcode)
    formula = mof.formula
    crys = mof.crystal  # Need crystal object to calcuate descriptors
    pore_analyser = CrystalDescriptors.PoreAnalyser(crys)
    He_vol_tot = pore_analyser.total_helium_volume
    sys_vol = pore_analyser.system_volume
    results.append({'Refcode': refcode,
                    'Formula': formula,
                    'He Volume (Å³)': He_vol_tot,
                    'System Vol (Å³)': sys_vol})
    print(f"{refcode:<10} {formula:<20} {He_vol_tot:15.3f} {sys_vol:18.3f}")

# Create DataFrame
df = pd.DataFrame(results)
df = df.round(3)

# Write data to csv
workdir = Path.cwd()
df.to_csv(f'{workdir}/results.csv', index=False, encoding='utf-8-sig')
