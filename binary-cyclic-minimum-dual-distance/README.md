# Binary Cyclic Codes with Simultaneously Near-Linear Minimum and Dual Distances

**Mathematical proof note — Version 1.0 — October 16, 2026**

Author: **Yang Shifu**

ORCID: [0009-0002-2768-8506](https://orcid.org/0009-0002-2768-8506)

Affiliation: National Mobile Communications Research Laboratory and the Frontiers Science Center for Mobile Information Communication and Security, Southeast University, Nanjing 211189, China.

## Result overview

For every odd integer $m\ge3$, put $n=2^m-1$. This note gives a self-contained probabilistic proof that ordinary binary cyclic codes of dimension $(n-1)/2$ can have both minimum distance and dual distance at least $(1/2-o(1))n/(m\log_2m)$. A finite-length sufficient condition is

$$H_2(\delta)+\delta\le\frac1{2m}.$$

The proof counts only low-weight words whose zeroth coordinate is one. Shift symmetry contributes a factor $\delta$ to the count, and one first-moment estimate for random reciprocal-factor choices completes the argument. A probability-tending-to-one version retains the same leading distance constant.

This is a probabilistic existence result. The finite-length condition is sufficient, not necessary. No fixed positive relative distance, deterministic efficient generator-selection algorithm, or decoding algorithm is claimed.

## Full text and source

- [Project page](https://stalin1943.github.io/Coding-Theory-Proof-Notes/binary-cyclic-minimum-dual-distance/)
- [Full proof note (PDF)](proof.pdf)
- [LaTeX source](proof.tex)
- [BibTeX](citation.bib)
- [PDF on GitHub](https://github.com/stalin1943/Coding-Theory-Proof-Notes/blob/main/binary-cyclic-minimum-dual-distance/proof.pdf)

## Suggested citation

S. Yang, "Binary Cyclic Codes with Simultaneously Near-Linear Minimum and Dual Distances," mathematical proof note, version 1.0, Oct. 16, 2026. [Online]. Available: https://github.com/stalin1943/Coding-Theory-Proof-Notes/blob/main/binary-cyclic-minimum-dual-distance/proof.pdf

```bibtex
@misc{Yang2026CyclicMinimumDualDistance,
  author = {Yang, Shifu},
  title = {Binary Cyclic Codes with Simultaneously Near-Linear Minimum and Dual Distances},
  year = {2026},
  month = oct,
  howpublished = {Mathematical proof note},
  note = {Version 1.0, October 16, 2026},
  url = {https://github.com/stalin1943/Coding-Theory-Proof-Notes/blob/main/binary-cyclic-minimum-dual-distance/proof.pdf}
}
```

## Build

```sh
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error proof.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error proof.tex
```

Keywords: binary cyclic codes; minimum distance; dual distance; self-orthogonal codes; cyclic rank; entropy method; probabilistic existence.

中文关键词：二元循环码；最小距离；对偶距离；自正交码；循环秩；熵方法；概率存在性。
