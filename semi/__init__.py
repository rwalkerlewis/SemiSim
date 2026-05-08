"""
kronos-semi: JSON-driven FEM device simulator on FEniCSx.

Public API:
    from semi import schema, run, materials, constants, scaling

Example:
    from semi import schema, run
    cfg = schema.load('benchmarks/pn_1d/pn_junction.json')
    result = run.run(cfg)
    # result.psi_phys, result.n_phys, result.p_phys are numpy arrays

The FEM-heavy modules (mesh, physics, solver, run) require dolfinx.
The pure-Python modules (constants, materials, scaling, doping, schema)
do not and can be imported and tested independently.
"""

__version__ = "0.24.0"

# Pure-Python imports, always safe (doping requires numpy so is not imported here)
from . import constants, materials, scaling, schema

# dolfinx-dependent modules are NOT imported at package level. Import them
# explicitly when you have dolfinx available:
#     from semi import run, mesh, solver
#     from semi.physics import poisson

__all__ = ["constants", "materials", "scaling", "schema", "__version__"]
