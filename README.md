# Enumeration of Minimum Weight Codewords of Pre&#x2011;Transformed Polar Codes by Tree Intersection
This is a C script for enumerating the minimum weight codewords of pre-transformed polar codes (PTPCs).

## Overview
Pre-transformation can improve the distance properties polar codes as reducing the number of minimum weight codewords 
and/or increasing the minimum distance.
PTPCs can be split into disjoint cosets whereby the number of minimum weight codewords can be determined by counting their number in each coset.

The implemented method uses the explicit construction of minimum weight codewords of "universal" polar cosets introduced in [[1]](#1).
A coset of a PTPC is always contained in the correspoding "universal" polar coset.
Consequently, the minimum weight codewords in a coset of a PTPC can be found by computing the intersection with the set of minimum weight codewords in the corresponding "universal" polar coset.
The algorithm efficiently counts the number of minimum weight codewords by traversing the intersection tree of the two codeword sets.

Please note that this method cannot find the number of minimum weight codewords if the pre-transformation increases the minimum distance of the polar code.

Further information on the algorithm can be found in the following paper: https://arxiv.org/abs/1812.08562

## Usage
``enumerate_minimum_weight_codewords.c`` gives a simple example for counting the minimum weight codewords of a PAC code.
Compile with:
```
gcc -march=native -Ofast -o enumerate_minimum_weight_codewords enumerate_minimum_weight_codewords.c -lm
```

## Citing
Please cite the following reference if you want to use this algorithm in your research:
```bibtex
@misc{zunker2023enumeration,
      title={Enumeration of Minimum Weight Codewords of Pre-Transformed Polar Codes by Tree Intersection}, 
      author={Andreas Zunker and Marvin Geiselhart and Stephan ten Brink},
      year={2023},
}
```

## References
<a id="1">[1]</a>
M. Rowshan, S. H. Dau and E. Viterbo, "On the Formation of Min-Weight Codewords of Polar/PAC Codes and Its Applications," in IEEE Transactions on Information Theory, vol. 69, no. 12, pp. 7627-7649, Dec. 2023, doi: 10.1109/TIT.2023.3319015.
