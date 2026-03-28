# Geodesic-Algorithm-for-Gaussian-Integer-Continued-Fractions
This contains a Python implementation of a geodesic algorithm for computing Gaussian integer continued fractions (GICFs), developed in my fourth-year mathematics dissertation.

## Overview

The algorithm extends the correspondence between continued fractions of real numbers and geodesics in hyperbolic two-space to hyperbolic three-space. This implementation provides a computational realisation of the geometric link.

## Contents

- `RealCF.py`
<br> Implementation of the classical continued fraction algorithm for real numbers.

- `ComplexCF.py`
<br> Implementation of a geodesic algorithm for Gaussian integer continued fractions. Also provided in the output of this code are the stages of simplifying a word over translations, inversions and rotations to allow us to read the expansion.

- `Report.pdf`
<br> This contains the full context including the definition of the maps that are used throughout the code and the reasoning behind their consideration.

## Notes
- Numerical behaviour may depend on floating point precision. The code defines a small epsilon error to help accommodate any rounding that may occur, however the continued fraction expansion may still suffer from floating point instability leading to incorrect digits for large iteration counts
- Full details of the theory and construction can be found in the accompanying dissertation
