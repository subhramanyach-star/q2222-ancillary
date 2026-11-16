"""Canonical 112-point efficiency grid of the continuum rate certificate, as exact rationals.

The list is the intended set: no duplicates, endpoints 0.75 and 0.98, and the six
table efficiencies.  Interval membership is decided with these rationals, not with
float tolerances.
"""

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

from fractions import Fraction as Fr

# Filenames in lossycone/data/t4, in increasing order.
EXPECTED_ETA_STRS = (
    '0.750', '0.751', '0.752', '0.753', '0.754', '0.755', '0.756', '0.757', '0.758', '0.759',
    '0.760', '0.761', '0.762', '0.763', '0.764', '0.765', '0.766', '0.767', '0.768', '0.769',
    '0.770', '0.771', '0.772', '0.773', '0.774', '0.775', '0.776', '0.777', '0.778', '0.779',
    '0.780', '0.782', '0.784', '0.786', '0.788', '0.790', '0.792', '0.794', '0.796', '0.798',
    '0.800', '0.802', '0.804', '0.806', '0.808', '0.810', '0.812', '0.814', '0.816', '0.818',
    '0.820', '0.822', '0.824', '0.826', '0.828', '0.830', '0.832', '0.834', '0.836', '0.838',
    '0.840', '0.842', '0.844', '0.846', '0.848', '0.850', '0.853', '0.856', '0.859', '0.862',
    '0.865', '0.868', '0.871', '0.874', '0.877', '0.880', '0.883', '0.886', '0.889', '0.892',
    '0.895', '0.898', '0.900', '0.901', '0.904', '0.907', '0.910', '0.913', '0.916', '0.919',
    '0.922', '0.925', '0.928', '0.931', '0.934', '0.937', '0.940', '0.943', '0.946', '0.949',
    '0.950', '0.952', '0.955', '0.958', '0.961', '0.964', '0.967', '0.970', '0.973', '0.976',
    '0.979', '0.980',
)
EXPECTED_ETAS = tuple(Fr(s) for s in EXPECTED_ETA_STRS)
EXPECTED_ETA_SET = frozenset(EXPECTED_ETAS)
ENDPOINTS = (Fr('0.75'), Fr('0.98'))
TABLE_ETAS = (Fr('0.75'), Fr('0.80'), Fr('0.85'), Fr('0.90'), Fr('0.95'), Fr('0.98'))
# Continuum interval guarantees of Theorem (cover), as exact rationals.
BANDS = (
    (Fr('0.75'), Fr('0.80'), Fr('5.3') / 100),
    (Fr('0.80'), Fr('0.85'), Fr('8.5') / 100),
    (Fr('0.85'), Fr('0.90'), Fr('10.7') / 100),
    (Fr('0.90'), Fr('0.95'), Fr('12.2') / 100),
    (Fr('0.95'), Fr('0.98'), Fr('8.0') / 100),
)
# Grid-point gain range 8.8--18.0 percent; table efficiencies 9.2 percent and above.
# Do not globally replace 8.8 with 9.2.
GRID_GAIN_MIN = Fr('8.8') / 100
TABLE_GAIN_MIN = Fr('9.2') / 100

if len(EXPECTED_ETAS) != 112:
    raise RuntimeError('EXPECTED_ETAS does not have 112 entries')
if len(set(EXPECTED_ETAS)) != 112:
    raise RuntimeError('EXPECTED_ETAS contains duplicate rationals')
if not (EXPECTED_ETAS[0] == ENDPOINTS[0] and EXPECTED_ETAS[-1] == ENDPOINTS[1]):
    raise RuntimeError('EXPECTED_ETAS endpoints are not 0.75 and 0.98')
if not all(e in EXPECTED_ETA_SET for e in TABLE_ETAS):
    raise RuntimeError('TABLE_ETAS is not a subset of the grid')
if not all(EXPECTED_ETAS[i] < EXPECTED_ETAS[i + 1] for i in range(111)):
    raise RuntimeError('EXPECTED_ETAS is not strictly increasing')
