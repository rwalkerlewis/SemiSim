"""
semisim_server: HTTP API for the SemiSim engine.

This package is deliberately separate from `semi/` to preserve the pure-engine
boundary documented in `docs/ARCHITECTURE.md`. The server process itself does
not import dolfinx, UFL, or PETSc; all FEM-heavy work happens inside worker
subprocesses.

Public API:
    from semisim_server.app import build_app, main
"""
