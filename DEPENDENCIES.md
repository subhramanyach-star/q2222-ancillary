# Dependencies

The stand-alone verifiers that support the computer-assisted claims of the accompanying paper need:

- Python 3
- python-flint 0.9 (`pip install python-flint==0.9.0`), for the detection-loss, CHSH-variant, and ball-arithmetic checks
- sympy, for some exact identities in the tilted-CHSH polygon checks

The published checkers use only the Python standard library together with those packages. They do not call an SDP or conic solver.

No network access is required once the repository is present. Run every verifier from the root of the tree, or from the directory named in `MANIFEST.md`; data are loaded by a path relative to the script.

The modules in `rigorous_rates/` are imported by the rate and variant verifiers. They are not a second, independent certificate collection.
