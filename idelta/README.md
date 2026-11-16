# I_δ supporting quadrilateral

`certify_idelta.py` checks the exact algebraic identities used in the
appendix on the supporting quadrilateral of the I_δ family of
Wooltorton, Brown and Colbeck (0 < δ ≤ π/6): the closed-form
vertices, multipliers, tight local behaviours, the twelve nontight
slacks, and the elementary sign comparisons on
cos δ ∈ [√3/2, 1).

```
python3 idelta/certify_idelta.py
```

The expected output is `All I_δ identities and interval sign checks passed.`
The check uses sympy. It is a verification of the algebra in the paper,
not a computer-assisted existence proof.
