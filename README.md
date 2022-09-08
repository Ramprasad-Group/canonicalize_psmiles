# Canonicalize PSMILES 

PSMILES (Polymer SMILES) is a chemical language to represent polymer structures. PSMILES strings have two stars (`[*]` or `*`) symbols that indicate the two endpoints of the polymer repeat unit and otherwise follow the daylight SMILES syntax defined at [OpenSmiles](http://opensmiles.org/opensmiles.html). 

The raw PSMILES syntax is ambiguous and non-unique; i.e., the same polymer may be written using many PSMILES strings:

Polyethylene | Polyethylene oxide | Polypropylene |
|-|-|-|
| `[*]C[*]`   | `[*]CCO[*]` | `[*]CC([*])C` | 
| `[*]CC[*]`  | `[*]COC[*]` | `[*]CC(CC([*])C)C` | 
| `[*]CCC[*]` | `[*]OCC[*]` | `CC([*])C[*]` | 

The canonicalization routine of the `PSMILES` packages finds a canonicalized version of the SMILES string by

1. Finding the shortest representation of a PSMILES string 

`[*]CCOCCO[*]` ->  `[*]CCO[*]`

2. Making the PSMILES string cyclic

`[*]CCO[*]` -> `C1 CCO C1`

3. Applying the canonicalization routine as implemented in RDKit

`C1 CCO C1` -> `C1 COC C1`

4. Breaking the cyclic bond

`C1 COC C1` -> `[*]COC[*]`

## Install 

```bash
pip install git+ssh://git@github.com/Ramprasad-Group/canonicalize_psmiles.git
```

## How to use

See also [test.ipynb](tests/test.ipynb)

```bash
from canonicalize_psmiles.canonicalize import canonicalize

smiles = "[*]NC(C)CC([*])=O"
print(smiles)
print(canonicalize(smiles))
```

