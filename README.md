# A one-mode quantization of the clock bound on a null cone, and its sensitivity to operator ordering

Author: Michaël Perreault

Version 1.0.1.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24575468.svg)](https://doi.org/10.5281/zenodo.24575468)

On a vacuum null cone the affine length of a generator is a clock, bounded below by a convex function of the conformal length. The paper asks whether this bound survives quantization in the null-cone scheme of Chandrasekaran. For a free system with N canonical pairs it proves an operator family bound, the expectation bound in every state, and the operator inequality on a spectral window of the length operator. Ordering matters at order hbar: the bound holds in the Weyl ordering and fails at the vacuum when the zero-point constant is removed from the shear, yet survives, at the vacuum and for small kappa, removal of any fraction below one of it. The full operator inequality is open; the paper states the obstruction and gives small-N numerics.

Companion: the classical bound is in [Affine length bounds the conformal excursion on a null cone: a sharp focusing inequality](https://doi.org/10.5281/zenodo.24572220).

Source: [paper.md](paper.md) and [a-one-mode-quantization-of-the-clock-bound-on-a-null-cone-and-its-sensitivity-to-operator-ordering.pdf](a-one-mode-quantization-of-the-clock-bound-on-a-null-cone-and-its-sensitivity-to-operator-ordering.pdf). Text: CC BY 4.0 ([LICENSE-paper.md](LICENSE-paper.md)). Code in `code/`: MIT ([LICENSE](LICENSE)).

Archive contents: `paper.md` and the PDF; `code/` (Python scripts and their outputs, with the commands in `code/out/README.md`).

## Re-run

Python 3 with numpy, scipy and mpmath. The commands and the output files are listed in [code/out/README.md](code/out/README.md); `code/run_table.sh PYTHON` reruns the numeric table.

## Cite as

Perreault, M. A one-mode quantization of the clock bound on a null cone, and its sensitivity to operator ordering, version 1.0.1. https://doi.org/10.5281/zenodo.24575468 (concept DOI, all versions). The first release, v1.0, is https://doi.org/10.5281/zenodo.24575469. See [CITATION.cff](CITATION.cff).
