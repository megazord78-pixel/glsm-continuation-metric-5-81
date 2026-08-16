"""Worker-local performance caches for the independent interval metric engine.

The formulas and Neumann series are unchanged.  Only midpoint inverses and
static parsed inputs that are reused with the identical matrix are cached.
"""

from __future__ import annotations

from functools import lru_cache


class CachedTextPath:
    def __init__(self, path):
        self.path = path
        self.text = path.read_text(encoding="utf-8")

    def read_text(self, encoding="utf-8"):
        if encoding.lower().replace("-", "") != "utf8":
            raise ValueError("cached parameter text is UTF-8 only")
        return self.text


def install_worker_caches():
    """Patch only the current child process and return instrumentation state."""
    from flint import arb, arb_mat
    import regeneration.models.string_5_81_glsm_transverse_tube_independent_formula_verifier as metric

    if getattr(metric, "_worker_cache_installed", False):
        return metric._worker_cache_stats

    metric.raw_inputs = lru_cache(maxsize=1)(metric.raw_inputs)
    metric.PARAMETERS = CachedTextPath(metric.PARAMETERS)
    state = {
        "solve_call_count": 0,
        "preconditioner_build_count": 0,
        "preconditioner_reuse_count": 0,
        "last_matrix": None,
        "last_component_signature": None,
        "last_factorization": None,
    }

    def cached_verified_solve(matrix, rhs, component_radii, arb_type, arb_mat_type, order=30):
        state["solve_call_count"] += 1
        signature = tuple(str(value) for value in component_radii)
        if matrix is state["last_matrix"] and signature == state["last_component_signature"]:
            preconditioner, scaled, rho = state["last_factorization"]
            state["preconditioner_reuse_count"] += 1
        else:
            n = matrix.nrows()
            midpoint = arb_mat_type([[matrix[i, j].mid() for j in range(n)] for i in range(n)])
            preconditioner = midpoint.inv()
            identity = arb_mat_type(n, n)
            for i in range(n):
                identity[i, i] = 1
            remainder = identity - preconditioner * matrix
            scaled = arb_mat_type([[
                remainder[i, j] * component_radii[j] / component_radii[i]
                for j in range(n)] for i in range(n)
            ])
            rho = max(sum((abs(scaled[i, j]) for j in range(n)), arb_type(0)) for i in range(n))
            if not rho < 1:
                raise ValueError(f"independent scaled solve norm is not below one: {rho}")
            state["last_matrix"] = matrix
            state["last_component_signature"] = signature
            state["last_factorization"] = (preconditioner, scaled, rho)
            state["preconditioner_build_count"] += 1
        initial_unscaled = preconditioner * rhs
        initial = arb_mat_type([[
            initial_unscaled[i, column] / component_radii[i]
            for column in range(rhs.ncols())
        ] for i in range(matrix.nrows())])
        initial_norm = max(
            sum((abs(initial[i, column]) for column in range(rhs.ncols())), arb_type(0))
            for i in range(matrix.nrows())
        )
        total = initial
        term = initial
        for _ in range(order):
            term = scaled * term
            total = total + term
        rho_upper = arb_type(rho.upper())
        tail = rho_upper ** (order + 1) / (1 - rho_upper) * initial_norm
        solved = arb_mat_type([[
            component_radii[i] * (total[i, column] + arb_type(0, tail))
            for column in range(rhs.ncols())
        ] for i in range(matrix.nrows())])
        return solved, rho, initial_norm, tail

    metric.verified_solve = cached_verified_solve
    metric._worker_cache_installed = True
    metric._worker_cache_stats = state
    return state


def snapshot_worker_cache_stats():
    import regeneration.models.string_5_81_glsm_transverse_tube_independent_formula_verifier as metric
    state = getattr(metric, "_worker_cache_stats", None)
    if state is None:
        return None
    return {key: value for key, value in state.items() if not key.startswith("last_")}
