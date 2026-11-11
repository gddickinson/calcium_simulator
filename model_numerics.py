"""Numerical helpers for CalciumModel: the 9-point Laplacian and stable explicit sub-step counts.

The arithmetic here is written out term by term (no scipy convolve) so that the JavaScript port in
web/static/ca_engine.js reproduces it to the last bit.
"""
import math

import numpy as np

# Isotropic 9-point Laplacian (Patra & Karttunen 2006), multiply by 1/dx^2.
LAPLACIAN_KERNEL = np.array([[1.0, 4.0, 1.0],
                             [4.0, -20.0, 4.0],
                             [1.0, 4.0, 1.0]]) / 6.0


def laplacian9(u, dx):
    """9-point Laplacian of u (um^-2 times u's units) with zero-flux (mirror) boundaries,
    i.e. scipy.ndimage 'reflect' mode."""
    p = np.pad(u, 1, mode="symmetric")
    n, s = p[:-2, 1:-1], p[2:, 1:-1]
    w, e = p[1:-1, :-2], p[1:-1, 2:]
    nw, ne, sw, se = p[:-2, :-2], p[:-2, 2:], p[2:, :-2], p[2:, 2:]
    return (4.0 * (((n + s) + w) + e) + (((nw + ne) + sw) + se) - 20.0 * u) / (6.0 * dx * dx)


def diffusion_rate(D, dx):
    """Largest decay rate (s^-1) the explicit 9-point update puts on a grid point: 20 D / (6 dx^2)."""
    return (20.0 * D) / (6.0 * dx * dx)


def substeps(dt, rate):
    """Number of forward-Euler sub-steps so that h * rate <= 1 with h = dt / n. With `rate` the sum of every
    first-order loss rate on a grid point, the centre weight 1 - h * rate stays >= 0: the update is stable
    and keeps concentrations non-negative."""
    return max(1, math.ceil(dt * rate))
