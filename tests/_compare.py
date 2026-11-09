"""Comparison helper: protpy no longer rounds descriptor values (v1.6.0), so expected
values written to 2-5 decimal places (from earlier, rounded outputs) are compared with a tolerance."""
import numpy as np


def frames_close(actual, expected, atol=5e-3):
    """Same columns, and values equal to within atol."""
    return (list(actual.columns) == list(expected.columns)
            and np.allclose(actual.values.astype(float), expected.values.astype(float), atol=atol))
