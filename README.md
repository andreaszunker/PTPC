# Enumeration of Minimum Weight Codewords of Pre&#x2011;Transformed Polar Codes by Tree Intersection
This is a C source file for enumerating the minimum weight codewords of pre-transformed polar codes (PTPCs).

## Overview
Pre-transformation can improve the distance properties polar codes by reducing the number of minimum weight codewords 
and/or increasing the minimum distance.
PTPCs can be split into disjoint cosets whereby the number of minimum weight codewords can be determined by counting their number in each coset.

The implemented method uses the explicit construction of the messages that generate the minimum weight codewords of "universal" polar cosets introduced in [[1]](#1).
A coset of a PTPC is always contained in the corresponding "universal" polar coset.
Consequently, the minimum weight codewords in a coset of a PTPC can be found by computing the intersection with the set of minimum weight codewords in the corresponding "universal" polar coset.
The algorithm efficiently counts the number of minimum weight codewords by traversing the intersection tree of the messages that form the two codeword sets.

Please note that this method cannot find the number of minimum weight codewords if the pre-transformation increases the minimum distance of the polar code.

Further information on the algorithm can be found in the following paper: https://arxiv.org/abs/2311.17774

## Usage
``enumeration.c`` gives a simple example for counting the minimum weight codewords of a PAC code.
Compile with:
```
gcc -march=native -Ofast -o enumeration enumeration.c
```
The code was tested on an Intel® Core™ i7-4790K CPU @ 4.00GHz.

## Citing
Please cite the following reference if you want to use this algorithm in your research:
```bibtex
@INPROCEEDINGS{10480163,
  author={Zunker, Andreas and Geiselhart, Marvin and Ten Brink, Stephan},
  booktitle={2024 58th Annual Conference on Information Sciences and Systems (CISS)}, 
  title={Enumeration of Minimum Weight Codewords of Pre-Transformed Polar Codes by Tree Intersection}, 
  year={2024},
  doi={10.1109/CISS59072.2024.10480163}}
```

## References
<a id="1">[1]</a>
M. Rowshan, S. H. Dau and E. Viterbo, "On the Formation of Min-Weight Codewords of Polar/PAC Codes and Its Applications," in IEEE Transactions on Information Theory, vol. 69, no. 12, pp. 7627-7649, Dec. 2023, doi: 10.1109/TIT.2023.3319015.
