"""Pickle-free text storage of the exact I-family polygon data (sympy objects stored as srepr strings in JSON).
dump(obj, path) / load(path).  Every sympy object round-trips exactly (sympify(srepr(e)) == e)."""
import json, sympy as sp
def _enc(o):
    if isinstance(o, sp.MatrixBase): return {'__matrix__': [[sp.srepr(e) for e in row] for row in o.tolist()]}
    if isinstance(o, sp.Basic): return {'__sympy__': sp.srepr(o)}
    if isinstance(o, dict): return {'__dict__': [[_enc(k), _enc(v)] for k, v in o.items()]}
    if isinstance(o, tuple): return {'__tuple__': [_enc(v) for v in o]}
    if isinstance(o, list): return [_enc(v) for v in o]
    if isinstance(o, (str, int, float, bool)) or o is None: return o
    raise TypeError(type(o))
def _dec(o):
    if isinstance(o, dict):
        if '__matrix__' in o: return sp.Matrix([[sp.sympify(e) for e in row] for row in o['__matrix__']])
        if '__sympy__' in o: return sp.sympify(o['__sympy__'])
        if '__dict__' in o: return {_dec(k): _dec(v) for k, v in o['__dict__']}
        if '__tuple__' in o: return tuple(_dec(v) for v in o['__tuple__'])
    if isinstance(o, list): return [_dec(v) for v in o]
    return o
def dump(obj, path): json.dump({'format': 'polygon_io v1', 'data': _enc(obj)}, open(path, 'w'), indent=0)
def load(path): return _dec(json.load(open(path))['data'])
