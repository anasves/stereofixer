# stereofixer

[![Tests](https://github.com/anasves/stereofixer/actions/workflows/tests.yml/badge.svg)](https://github.com/anasves/stereofixer/blob/main/.github/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/anasves/stereofixer/badge.svg?branch=main)](https://coveralls.io/github/anasves/stereofixer?branch=main)

## Installation
```bash
git clone https://github.com/anasves/stereofixer.git 
cd stereofixer
conda create -n stereofixer python=3.9 -y
conda activate stereofixer
pip install -e .
```

## Usage
1. Check reaction for stereo issues
```python
from stereofixer.StereoChecker import StereoChecker
sc = StereoChecker()
input_stereo_string1 = 'C(C)(CC)(CCC)[O-]>>[C@](C)(CC)(CCC)[O-]'
result = sc.smiles_stereo_analysis(input_stereo_string1)
print(result)
# Output:
# {'mapped_rxn': '[CH3:1][CH2:2][CH2:3][C:4]([CH3:5])([O-:6])[CH2:7][CH3:8]>>
#             [CH3:1][CH2:2][CH2:3][C@@:4]([CH3:5])([O-:6])[CH2:7][CH3:8]',
# 'stereo_difference': '4-?;4-S',
# 'reaction_stereoisomer_count_total': '2>>2',
# 'reaction_stereoisomer_count_unenumerated': '2>>1',
# 'current_stereo_atoms': '>>4',
# 'all_theoretically_stereo_atoms': '4>>4',
# 'disappearing_stereo_atoms': [],
# 'unassigned_stereo_atoms': [4],
# 'num_unassigned_stereo_atoms': 1,
# 'case': 'missing stereo, enumerate',
# 'mismatched_atoms': [(4, {'S', '?'})]}
```

2. Enumerate stereo versions of a reaction
```python
from stereofixer.StereoEnumerator import StereoEnumerator
stereo_enumerator = StereoEnumerator()
result = stereo_enumerator.enumerate_stereoisomers_per_reaction('[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]', mismatched_atoms=[3,])
print(result)
# Output:
#(['[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]',
# '[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]',
# '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]',
# '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]'], 4)
```

3. Filter out versions of reaction with stereoissues
```python
from stereofixer.PandasStereoEnumerator import PandasStereoEnumerator
pandas_stereo_enumerator = PandasStereoEnumerator()
reaction_stero_options =
  ['[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]',
   '[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]',
   '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]',
   '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]']

# If executed after previous block of code,
#reaction_stero_options = result[0]

df_filtered = self.pandas_stereo_enumerator.filter_stereo_options_results(reaction_stero_options)
print(df_filtered)
# Output:
#    SMILES  ...      case
# 1  [CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>...  ...  all good
# 2  [CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>...  ...  all good
#
# [2 rows x 3 columns]

```
