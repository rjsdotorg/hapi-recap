# README Refactoring Notes

This document summarizes the recent documentation and code-maintenance work in this repository.

## README updates

The main README was expanded to better explain:

- What HAPI-RECAP does
- What inputs the script expects
- What outputs it produces
- What HAPI2 outputs are required for HAPI-RECAP to run
- An example workflow for running HAPI2 first, then HAPI-RECAP

## Script documentation updates

The Python script was updated with:

- A module-level docstring
- SciPy/NumPy-style docstrings for all major functions

## Refactoring and robustness improvements

A few low-risk changes were added to improve speed and data safety:

- Added lightweight input validation helpers
- Added checks for required HAPI2 JSON keys
- Added checks for required IBD Feather columns
- Added a guard for missing IBD Feather files
- Reduced repeated DataFrame filtering in the inner loops by grouping IBD segments per chromosome and relative ID
- Replaced some repeated overlap-coordinate calculations with a helper function
- Disabled unnecessary pandas groupby sorting where order is not needed
- Added a guard against division by zero when no overlaps are found

## Validation

The script was checked after the edits and no syntax errors were reported.
