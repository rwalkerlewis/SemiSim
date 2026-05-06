"""
Small-signal AC sweep runner (M14, Slotboom-form).

Linearises the steady-state drift-diffusion system around a converged
DC operating point and solves the complex linear system

    (J + j*omega*M) * delta_u(omega) = -dF/dV * delta_V

at each frequency in a user-specified sweep.

Primary unknowns
----------------
The primary unknowns are the same Slotboom variables (psi, phi_n,
phi_p) used by ``run_bias_sweep`` and ``run_transient`` (after the
M13.1 Slotboomisation; see ADR 0014 for the transient case and the
"Errata #2" section of ADR 0011 for why the AC runner was switched
back to Slotboom in the audit-Cases-02/05 closure).

Why Slotboom (and not (n,p)-primary form):

- The DC operating point is solved in Slotboom by ``run_bias_sweep``.
  Computing the Jacobian J in (n, p)-primary form at a Slotboom-
  converged-then-Boltzmann-converted state gives an inconsistent
  linearisation: the (n, p)-form discrete residual is non-zero on
  the converted state by O(h^p) discretisation mismatch, and the
  ``b = J @ delta_u_DC`` identity ``J*du/dV = -dF/dV`` no longer
  holds. Audit Case 05 (forward bias) hit this with a 19% magnitude
  disagreement vs ``bias_sweep`` centred-difference dI/dV; see ADR
  0011 Errata #2 for the diagnosis.
- The terminal current evaluator used by ``bias_sweep`` /
  ``postprocess.evaluate_current_at_contact`` is in Slotboom form
  (``q*mu_n*n*grad(phi_n).n_outward``). Linearising that same
  expression by ``ufl.derivative`` for the AC small-signal current
  is exact and avoids any Slotboom <-> (n,p) sign / convention
  mismatch.

The mass matrix is no longer a simple lumped diagonal on the (n, p)
rows. In Slotboom form the time-derivative term ``d n / dt`` becomes
``(n / V_t) (d psi/dt - d phi_n/dt)`` and similarly for holes; the
chain-rule mass matrix is a sparse 3x3 block. We assemble it directly
from ``ufl.derivative`` of the lumped continuity time term -- the same
trick the M13.1 transient runner uses for its BDF time term (ADR 0014
section "Implementation").

PETSc complex-vs-real
---------------------
The dolfinx-real PETSc build (the only build available in CI today)
uses real scalars. M14 implements the standard real 2x2 block
reformulation:

    [   J     -omega*M ] [ Re delta_u ]   [ Re b ]
    [ omega*M    J     ] [ Im delta_u ] = [ Im b ]

The block size is 2 * (3N) = 6N where N is the per-block DOF count.
A direct LU (MUMPS) factorisation is used at each frequency. For the
benchmark mesh sizes shipped with M14 (~80-200 cells in 1D), this is
sub-second per frequency.

Sign convention
---------------
Y is reported with the "positive = current INTO device" convention
used by ``bias_sweep`` / ``postprocess.evaluate_current_at_contact``.
The Slotboom-form linearised terminal current produces this sign
directly (no global negation needed), so an ideal capacitor terminal
has Y = +j*omega*C and C = +Im(Y)/(2*pi*f) is positive for a positive
depletion capacitance.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def _resolve_frequencies(freq_spec: dict) -> np.ndarray:
    """Expand the schema's frequency-sweep object to a numpy array of Hz."""
    kind = freq_spec["type"]
    if kind == "list":
        return np.asarray(freq_spec["values"], dtype=float)
    if kind == "logspace":
        start = float(freq_spec["start"])
        stop = float(freq_spec["stop"])
        n = int(freq_spec["n_points"])
        return np.logspace(np.log10(start), np.log10(stop), n)
    if kind == "linspace":
        start = float(freq_spec["start"])
        stop = float(freq_spec["stop"])
        n = int(freq_spec["n_points"])
        return np.linspace(start, stop, n)
    raise ValueError(f"Unknown ac.frequencies.type {kind!r}")


def _build_dc_subcfg(cfg: dict, contact_name: str, voltage: float) -> dict:
    """Build a bias_sweep cfg pinned at a single DC bias point.

    Returns a deep copy of `cfg` with solver.type set to "bias_sweep" and
    a single-point voltage_sweep on `contact_name` at `voltage`. Any
    pre-existing voltage_sweep on other contacts is removed (the sweep
    runner picks the first ohmic contact with a voltage_sweep, so we
    have to be the only one).
    """
    import copy

    sub = copy.deepcopy(cfg)
    snes_overrides = (cfg.get("solver", {}).get("snes") or {})
    user_cont = (cfg.get("solver", {}).get("continuation") or {})
    target_abs = abs(float(voltage))
    nominal_step = max(min(target_abs, 0.1), 1.0e-3)
    cont = dict(user_cont)
    cont["max_step"] = max(float(cont.get("max_step", nominal_step)), nominal_step)
    cont["min_step"] = min(
        float(cont.get("min_step", 1.0e-4)),
        cont["max_step"] * 0.1,
        max(nominal_step * 0.01, 1.0e-7),
    )
    cont.setdefault("max_halvings", 10)
    cont.setdefault("easy_iter_threshold", 4)
    cont.setdefault("grow_factor", 1.5)
    sub["solver"] = {
        "type": "bias_sweep",
        "continuation": cont,
    }
    # Propagate the user's backend selection into the DC sub-solve so
    # the GPU path applies to the AC operating-point computation too.
    src_solver = cfg.get("solver", {}) or {}
    if "backend" in src_solver:
        sub["solver"]["backend"] = src_solver["backend"]
    if "compute" in src_solver:
        sub["solver"]["compute"] = copy.deepcopy(src_solver["compute"])
    if snes_overrides:
        sub["solver"]["snes"] = dict(snes_overrides)
    found = False
    for c in sub["contacts"]:
        if c["name"] == contact_name:
            c["voltage"] = float(voltage)
            c["voltage_sweep"] = {
                "start": float(voltage),
                "stop": float(voltage),
                "step": nominal_step,
            }
            found = True
        else:
            c.pop("voltage_sweep", None)
    if not found:
        raise ValueError(
            f"ac_sweep dc_bias.contact {contact_name!r} not found in cfg.contacts"
        )
    return sub


def run_ac_sweep(cfg: dict[str, Any], *, progress_callback=None):
    """
    Small-signal AC sweep around a DC operating point (Slotboom form).

    Parameters
    ----------
    cfg : dict
        Validated config dict. ``solver.type`` must be ``"ac_sweep"`` and
        ``solver.dc_bias`` and ``solver.ac`` must be provided per the
        schema (see ``schemas/input.v1.json``).
    progress_callback : callable or None
        If provided, called with ``{"type": "ac_step", "freq_index": k,
        "frequency": f}`` after each frequency point.

    Returns
    -------
    AcSweepResult
        See :class:`semi.results.AcSweepResult`.
    """
    import math

    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import assemble_matrix

    from ..bcs import build_dd_dirichlet_bcs, resolve_contacts
    from ..doping import build_profile
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import build_dd_block_residual, make_dd_block_spaces
    from ..physics.slotboom import n_from_slotboom, p_from_slotboom
    from ..postprocess import resolve_contact_facets
    from ..results import AcSweepResult
    from ..runners.bias_sweep import run_bias_sweep
    from ..scaling import make_scaling_from_config
    from ._common import reference_material

    solver_cfg = cfg.get("solver", {})
    if solver_cfg.get("type") != "ac_sweep":
        raise ValueError("run_ac_sweep requires solver.type == 'ac_sweep'")
    dc_cfg = solver_cfg["dc_bias"]
    ac_cfg = solver_cfg["ac"]
    contact_name = dc_cfg["contact"]
    V_DC = float(dc_cfg["voltage"])
    ac_amplitude = float(ac_cfg.get("amplitude", 1.0))
    frequencies = _resolve_frequencies(ac_cfg["frequencies"])

    # ------------------------------------------------------------------
    # 1. Solve the DC operating point at V_DC and at V_DC + eps_V in
    #    Slotboom variables. The pair lets us form the DC sensitivity
    #    delta_u_DC = du_slot/dV by finite difference, which gives
    #    J * delta_u_DC = -dF/dV at u_0 because
    #
    #        F_slot(u_slot(V), V) == 0    for both V points (modulo SNES)
    #
    #    so d/dV F_slot(u_slot(V), V) = J*du/dV + dF/dV = 0. The
    #    bias_sweep runner is reused for these solves so we do not
    #    duplicate the multi-region BC and SNES plumbing.
    # ------------------------------------------------------------------
    eps_V = 1.0e-3

    dc_sub = _build_dc_subcfg(cfg, contact_name, V_DC)
    res_dc = run_bias_sweep(dc_sub)
    dc_iter = int(res_dc.solver_info.get("iterations", 0))

    dc_sub_eps = _build_dc_subcfg(cfg, contact_name, V_DC + eps_V)
    res_dc_eps = run_bias_sweep(dc_sub_eps)
    dc_iter = max(dc_iter, int(res_dc_eps.solver_info.get("iterations", 0)))

    # ------------------------------------------------------------------
    # 2. Build a fresh (psi, phi_n, phi_p) Slotboom function-space
    #    stack on the same mesh and seed it with the bias_sweep DC
    #    operating point at V_DC. We keep our own mesh build so the
    #    runner is independent of the dc-subcfg's plotter / artifact
    #    side effects.
    # ------------------------------------------------------------------
    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])
    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi
    V_phi_n = spaces.V_phi_n
    V_phi_p = spaces.V_phi_p

    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    psi.x.array[:] = res_dc.psi.x.array
    phi_n.x.array[:] = res_dc.phi_n.x.array
    phi_p.x.array[:] = res_dc.phi_p.x.array
    for fn in (psi, phi_n, phi_p):
        fn.x.scatter_forward()

    # DC sensitivity delta_u_DC per unit dV in Slotboom variables
    # (concatenated [psi, phi_n, phi_p]).
    delta_psi_DC = (res_dc_eps.psi.x.array - res_dc.psi.x.array) / eps_V
    delta_phi_n_DC = (res_dc_eps.phi_n.x.array - res_dc.phi_n.x.array) / eps_V
    delta_phi_p_DC = (res_dc_eps.phi_p.x.array - res_dc.phi_p.x.array) / eps_V

    N_psi = psi.x.array.size
    N_n = phi_n.x.array.size

    # ------------------------------------------------------------------
    # 3. Material constants and physics-model coefficients (same as
    #    bias_sweep / transient).
    # ------------------------------------------------------------------
    phys = cfg.get("physics", {})
    mob = phys.get("mobility", {})
    mu_n_SI = float(mob.get("mu_n", 1400.0)) * 1.0e-4
    mu_p_SI = float(mob.get("mu_p", 450.0)) * 1.0e-4
    mu_n_hat = mu_n_SI / sc.mu0
    mu_p_hat = mu_p_SI / sc.mu0

    rec = dict(phys.get("recombination", {}))
    tau_n_s = float(rec.get("tau_n", 1.0e-7))
    tau_p_s = float(rec.get("tau_p", 1.0e-7))
    E_t_eV = float(rec.get("E_t", 0.0))
    tau_n_hat_v = tau_n_s / sc.t0
    tau_p_hat_v = tau_p_s / sc.t0
    E_t_over_Vt = E_t_eV / sc.V0
    # M16.6: merge physics.tunneling into recomb_cfg; see bias_sweep
    # for the rationale.
    tun = phys.get("tunneling", {})
    for key, value in tun.items():
        rec[key] = value

    # ------------------------------------------------------------------
    # 4. Steady-state Slotboom residual and its Jacobian. Identical to
    #    the residual built by `build_dd_block_residual` (the same form
    #    bias_sweep solves), so the AC linearisation is *exactly* the
    #    DC SNES Jacobian. UFL automatic differentiation handles dR/du
    #    (SRH recombination) and the n*grad(phi_n) drift coupling
    #    consistently.
    # ------------------------------------------------------------------
    stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}

    F_list = build_dd_block_residual(
        spaces, N_hat_fn, sc, ref_mat.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n_hat_v, tau_p_hat_v, E_t_over_Vt,
        recomb_cfg=rec,
        statistics_cfg=stat_cfg,
    )
    a_blocks = [
        [ufl.derivative(F_list[i], u_j) for u_j in (psi, phi_n, phi_p)]
        for i in range(3)
    ]

    # ------------------------------------------------------------------
    # 5. Build BCs at the DC operating point. The block Jacobian is
    #    assembled with these BCs so BC rows are zeroed with diag=1
    #    (standard PETSc lifting). All three Slotboom unknowns shift
    #    with the applied bias at the swept ohmic contact:
    #
    #        psi_bc   = psi_eq_hat + V/V_t
    #        phi_n_bc = V/V_t
    #        phi_p_bc = V/V_t
    #
    #    so delta_u_DC at the swept-contact BC dofs is +1/V_t for all
    #    three rows, and the AC linearised system inherits the same BC
    #    perturbation through `b = J @ delta_u_DC`.
    # ------------------------------------------------------------------
    static_voltages = {
        c["name"]: float(c.get("voltage", 0.0))
        for c in cfg["contacts"]
        if c["type"] == "ohmic"
    }
    static_voltages[contact_name] = V_DC

    contacts_dc = resolve_contacts(cfg, facet_tags=facet_tags, voltages=static_voltages)
    bcs_dc = build_dd_dirichlet_bcs(
        spaces, msh, facet_tags, contacts_dc, sc, ref_mat, N_raw_fn,
    )

    # Compiled forms for assembly.
    a_forms = [[fem.form(a_blocks[i][j]) for j in range(3)] for i in range(3)]

    A_J = assemble_matrix(a_forms, bcs=bcs_dc, diag=1.0)
    A_J.assemble()

    # ------------------------------------------------------------------
    # 6. Mass matrix in Slotboom form.
    #
    #    The transient continuity equations in Slotboom form (see
    #    semi.runners.transient._build_transient_residual and ADR 0014
    #    section "Implementation") add the lumped time term
    #
    #        F_mass_n = (alpha_0/dt_hat) * n_ufl * v_n * dx_lump
    #        F_mass_p = (alpha_0/dt_hat) * p_ufl * v_p * dx_lump
    #
    #    In frequency domain d/dt -> j*omega and (1/dt_phys) ->
    #    j*omega, so the bilinear coefficient on the time term becomes
    #    j*omega*t0 (with t0 the physical time scale used to non-
    #    dimensionalise dt_hat). The mass operator M is the Jacobian of
    #    (n_ufl * v_n + p_ufl * v_p) * dx_lump with respect to (psi,
    #    phi_n, phi_p); UFL produces the chain-rule entries
    #
    #        d(n_ufl)/d(psi)   = +n_ufl       d(n_ufl)/d(phi_n) = -n_ufl
    #        d(p_ufl)/d(phi_p) = +p_ufl       d(p_ufl)/d(psi)   = -p_ufl
    #
    #    automatically. Lumping (vertex quadrature) localises each entry
    #    to a single DOF, matching the M13 / M13.1 lumped-mass contract.
    #    Poisson row has no time term (M_psi_block == 0).
    # ------------------------------------------------------------------

    from petsc4py import PETSc as _PETSc

    ni_hat_c = fem.Constant(msh, _PETSc.ScalarType(sc.n_i / sc.C0))
    eta_offset_n_val = (
        sc.eta_offset_n
        if stat_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    eta_offset_p_val = (
        sc.eta_offset_p
        if stat_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    n_ufl = n_from_slotboom(
        psi, phi_n, ni_hat_c,
        statistics_cfg=stat_cfg, eta_offset_n=eta_offset_n_val,
    )
    p_ufl = p_from_slotboom(
        psi, phi_p, ni_hat_c,
        statistics_cfg=stat_cfg, eta_offset_p=eta_offset_p_val,
    )

    v_psi = ufl.TestFunction(V_psi)
    v_n = ufl.TestFunction(V_phi_n)
    v_p = ufl.TestFunction(V_phi_p)

    dx_lump = ufl.dx(metadata={
        "quadrature_rule": "vertex",
        "quadrature_degree": 1,
    })

    # Poisson row has no time term, but UFL refuses to compile a
    # bare 0*test*dx without a function-valued integrand. Multiply by
    # a zero Constant on the mesh so the form is well-formed and
    # contributes a zero block.
    zero_psi_c = fem.Constant(msh, _PETSc.ScalarType(0.0))
    F_mass = [
        zero_psi_c * v_psi * ufl.dx,
        n_ufl * v_n * dx_lump,
        p_ufl * v_p * dx_lump,
    ]
    m_blocks = [
        [ufl.derivative(F_mass[i], u_j) for u_j in (psi, phi_n, phi_p)]
        for i in range(3)
    ]
    m_forms = [[fem.form(m_blocks[i][j]) for j in range(3)] for i in range(3)]

    # diag=0.0 zeros BC rows (and cols) without setting an identity on
    # the diagonal. The full system at BC rows then reads
    #   (J + j*omega*M)_bc = (1 + 0) * delta_u_bc = delta_u_bc
    # i.e. the BC perturbation per-unit-V is imposed identically at
    # every frequency, which matches the lifted Dirichlet convention.
    A_M = assemble_matrix(m_forms, bcs=bcs_dc, diag=0.0)
    A_M.assemble()

    # ------------------------------------------------------------------
    # 7. RHS direction b = J @ delta_u_DC. By construction (see step 1)
    #    this equals -dF/dV at u_0; at BC rows it equals delta_u_DC at
    #    the BC dof, which is exactly the BC perturbation per-unit-V.
    # ------------------------------------------------------------------
    delta_u_DC_3N = np.concatenate([delta_psi_DC, delta_phi_n_DC, delta_phi_p_DC])
    b_RHS = np.zeros_like(delta_u_DC_3N)
    _matvec_into(A_J, delta_u_DC_3N, b_RHS)

    # ------------------------------------------------------------------
    # 8. Frequency loop. At each omega build the 6N x 6N real 2x2 block
    #    system, solve via direct LU, and compute terminal current
    #    perturbation Y(omega).
    # ------------------------------------------------------------------
    Y_list: list[complex] = []
    Z_list: list[complex] = []
    C_list: list[float] = []
    G_list: list[float] = []

    sweep_facet_info = resolve_contact_facets(cfg, msh, facet_tags, contact_name)

    psi_pert = fem.Function(V_psi, name="delta_psi")
    phi_n_pert = fem.Function(V_phi_n, name="delta_phi_n")
    phi_p_pert = fem.Function(V_phi_p, name="delta_phi_p")

    # Pre-build the linearised terminal-current form and a 1*ds area
    # form. UFL ufl.derivative on the Slotboom-form bias_sweep current
    # expression (q*mu_n*n*grad(phi_n).n_outward + q*mu_p*p*grad(phi_p)
    # .n_outward) handles the chain rule on n_ufl(psi, phi_n) and the
    # grad(phi_n) drift coupling cleanly.
    cond_form, area_form = _make_terminal_current_forms(
        psi, phi_n, phi_p,
        psi_pert, phi_n_pert, phi_p_pert,
        sc, ref_mat, mu_n_SI, mu_p_SI,
        sweep_facet_info, msh,
        statistics_cfg=stat_cfg,
    )

    for k, f in enumerate(frequencies):
        omega = 2.0 * math.pi * float(f)

        x_real, x_imag = _solve_2x2_real_block(A_J, A_M, omega * sc.t0, b_RHS)

        # Conduction current perturbation at the contact.
        I_cond_re = _eval_terminal_current(
            x_real, N_psi, N_n,
            psi_pert, phi_n_pert, phi_p_pert,
            cond_form, area_form, msh,
        )
        I_cond_im = _eval_terminal_current(
            x_imag, N_psi, N_n,
            psi_pert, phi_n_pert, phi_p_pert,
            cond_form, area_form, msh,
        )

        # Displacement current at the contact:
        #   I_disp = j*omega*Q_psi(delta_psi)  in the Slotboom (= INTO)
        # convention where Q_psi is +eps*grad(delta_psi).n_outward and
        # is the negative of the displacement-current density flowing
        # OUT. So I_disp_into = +j*omega*Q_psi.
        Q_psi_re = _eval_displacement_charge(
            x_real[:N_psi], psi_pert, sc, ref_mat.epsilon_r,
            sweep_facet_info, msh,
        )
        Q_psi_im = _eval_displacement_charge(
            x_imag[:N_psi], psi_pert, sc, ref_mat.epsilon_r,
            sweep_facet_info, msh,
        )

        # I_total = I_cond + j*omega*Q_psi
        I_total_re = I_cond_re - omega * Q_psi_im
        I_total_im = I_cond_im + omega * Q_psi_re

        Y = complex(I_total_re, I_total_im)
        Y_list.append(Y)

        if abs(Y) < 1.0e-300:
            Z_list.append(complex(0.0, 0.0))
        else:
            Z_list.append(1.0 / Y)

        if f > 0.0:
            C_list.append(Y.imag / (2.0 * math.pi * float(f)))
        else:
            C_list.append(0.0)
        G_list.append(Y.real)

        if progress_callback is not None:
            progress_callback({
                "type": "ac_step",
                "freq_index": k,
                "frequency": float(f),
            })

    return AcSweepResult(
        frequencies=[float(f) for f in frequencies],
        Y=Y_list,
        Z=Z_list,
        C=C_list,
        G=G_list,
        dc_bias={"contact": contact_name, "voltage": V_DC},
        meta={
            "n_freqs": len(frequencies),
            "perturbation_voltage": float(eps_V),
            "ac_amplitude": ac_amplitude,
            "dc_solver_iterations": int(dc_iter),
            "backend": "real-2x2-block",
            "displacement_current_included": True,
            "primary_unknowns": "slotboom",
        },
    )


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _matvec_into(A, x_arr: np.ndarray, out_arr: np.ndarray) -> None:
    """y = A @ x for monolithic block PETSc Mats, via local arrays."""
    x_vec = A.createVecRight()
    y_vec = A.createVecLeft()
    x_vec.array[:] = x_arr
    A.mult(x_vec, y_vec)
    out_arr[:] = y_vec.array
    x_vec.destroy()
    y_vec.destroy()


def _solve_2x2_real_block(A_J, A_M, omega_scaled: float, b_RHS: np.ndarray):
    """Solve the real 2x2 block AC system [J -wM; wM J][x;y]=[b;0].

    Parameters
    ----------
    A_J : PETSc.Mat
        Steady-state Jacobian (3N x 3N) with BC rows replaced by
        identity (diag=1.0).
    A_M : PETSc.Mat
        Mass matrix (3N x 3N) with BC rows zeroed (diag=0.0). Sparse;
        in Slotboom form not just diagonal -- the chain-rule on n_ufl,
        p_ufl couples (phi_n, psi) and (phi_p, psi).
    omega_scaled : float
        Already pre-multiplied by t0 (non-dimensional time scale) so
        that omega is in physical rad/s when forming the 2x2 real
        block and (J + j*omega*M) is dimensionally consistent with
        d/dt -> j*omega in the scaled continuity equation.
    b_RHS : np.ndarray
        Real-part RHS of length 3N; imaginary part is zero by
        construction (the BC perturbation has no imaginary component).

    Returns
    -------
    (x_real, x_imag) : tuple of np.ndarray
        Length-3N numpy arrays, the real and imaginary parts of the
        solution. Direct LU via MUMPS at each frequency.
    """
    from petsc4py import PETSc
    from scipy.sparse import bmat, csr_matrix

    N3 = b_RHS.size
    big_size = 2 * N3

    J_indptr, J_indices, J_data = A_J.getValuesCSR()
    M_indptr, M_indices, M_data = A_M.getValuesCSR()

    J_sp = csr_matrix((J_data, J_indices, J_indptr), shape=(N3, N3))
    M_sp = csr_matrix((M_data, M_indices, M_indptr), shape=(N3, N3))
    wM = (omega_scaled * M_sp).tocsr()

    big = bmat([[J_sp, -wM], [wM, J_sp]], format="csr")

    A_big = PETSc.Mat().createAIJ(
        size=(big_size, big_size),
        csr=(big.indptr.astype(PETSc.IntType),
             big.indices.astype(PETSc.IntType),
             big.data.astype(PETSc.RealType)),
        comm=A_J.getComm(),
    )
    A_big.assemble()

    b_big = A_big.createVecLeft()
    x_big = A_big.createVecRight()
    b_big.array[:N3] = b_RHS
    b_big.array[N3:] = 0.0

    ksp = PETSc.KSP().create(A_big.getComm())
    ksp.setOperators(A_big)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.solve(b_big, x_big)

    x_arr_full = x_big.array.copy()
    x_real = x_arr_full[:N3].copy()
    x_imag = x_arr_full[N3:].copy()

    ksp.destroy()
    A_big.destroy()
    b_big.destroy()
    x_big.destroy()

    return x_real, x_imag


def _make_terminal_current_forms(
    psi, phi_n, phi_p,
    psi_pert, phi_n_pert, phi_p_pert,
    sc, ref_mat, mu_n_SI: float, mu_p_SI: float,
    facet_info, msh,
    *, statistics_cfg: dict | None = None,
):
    """Build the linearised Slotboom-form terminal-current form.

    The unlinearised (steady-state) terminal current density at the
    contact is, in Slotboom form (matching
    ``semi.postprocess.evaluate_current_at_contact``):

        J_total = q*mu_n*n*grad(phi_n).n_outward
                + q*mu_p*p*grad(phi_p).n_outward

    where n = ni*exp(psi - phi_n), p = ni*exp(phi_p - psi) (in scaled
    Slotboom form). UFL ``ufl.derivative`` against (psi, phi_n, phi_p)
    with directions (psi_pert, phi_n_pert, phi_p_pert) produces the
    linearised expression with all chain-rule terms (in particular
    the n*(delta_psi - delta_phi_n) carrier-density coupling and the
    n_0 * grad(delta_phi_n) drift coupling) included consistently.

    Returns
    -------
    (cond_form, area_form) : tuple of compiled fem.form
        cond_form integrates the linearised total current density
        over the contact face. area_form integrates 1.0 over the same
        face (used to convert the integral to a per-unit-area value
        consistent with the existing 1D / 2D current convention).
    """
    import ufl
    from dolfinx import fem
    from dolfinx.mesh import meshtags
    from petsc4py import PETSc

    from ..constants import Q
    from ..physics.slotboom import n_from_slotboom, p_from_slotboom

    tag = int(facet_info["tag"])
    facets = facet_info["facets"]
    fdim = msh.topology.dim - 1
    if len(facets) == 0:
        return None, None

    values = np.full(len(facets), tag, dtype=np.int32)
    indices = np.asarray(facets, dtype=np.int32)
    sort = np.argsort(indices)
    facet_mt = meshtags(msh, fdim, indices[sort], values[sort])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)
    n_vec = ufl.FacetNormal(msh)

    ni_hat_c = fem.Constant(msh, PETSc.ScalarType(ref_mat.n_i / sc.C0))

    # Slotboom carrier densities in physical (per-m^3) units. Threads
    # the M16.4 statistics dispatch so the FD path picks up the
    # Blakemore prefactor in the terminal-current expression.
    eta_offset_n_t = (
        sc.eta_offset_n
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    eta_offset_p_t = (
        sc.eta_offset_p
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    n_phys = n_from_slotboom(
        psi, phi_n, ni_hat_c,
        statistics_cfg=statistics_cfg, eta_offset_n=eta_offset_n_t,
    ) * sc.C0
    p_phys = p_from_slotboom(
        psi, phi_p, ni_hat_c,
        statistics_cfg=statistics_cfg, eta_offset_p=eta_offset_p_t,
    ) * sc.C0
    grad_phi_n_phys = sc.V0 * ufl.grad(phi_n)
    grad_phi_p_phys = sc.V0 * ufl.grad(phi_p)

    # Steady-state Slotboom terminal current density (matches
    # postprocess.evaluate_current_at_contact: positive = current INTO
    # device when n_vec is outward).
    Jn = Q * mu_n_SI * n_phys * ufl.dot(grad_phi_n_phys, n_vec)
    Jp = Q * mu_p_SI * p_phys * ufl.dot(grad_phi_p_phys, n_vec)
    J_total_expr = (Jn + Jp) * ds

    # Linearise wrt each Slotboom unknown with the matching
    # perturbation Function as the direction. ufl.derivative handles
    # the n_ufl(psi, phi_n) and p_ufl(phi_p, psi) chain rules.
    dJ_dpsi = ufl.derivative(J_total_expr, psi, psi_pert)
    dJ_dphi_n = ufl.derivative(J_total_expr, phi_n, phi_n_pert)
    dJ_dphi_p = ufl.derivative(J_total_expr, phi_p, phi_p_pert)
    delta_J_form = dJ_dpsi + dJ_dphi_n + dJ_dphi_p

    cond_form = fem.form(delta_J_form)
    area_form = fem.form(1.0 * ds)
    return cond_form, area_form


def _eval_terminal_current(
    x_3N: np.ndarray, N_psi: int, N_n: int,
    psi_pert, phi_n_pert, phi_p_pert,
    cond_form, area_form, msh,
) -> float:
    """Evaluate the linearised terminal current for one part of x.

    `x_3N` is the concatenated [delta_psi | delta_phi_n | delta_phi_p]
    array (real or imaginary part of the AC solution). The pert
    Functions are filled in place and the pre-compiled cond_form is
    assembled to a scalar. Returns a per-unit-area current density in
    A/m^2 (1D conventions: divides by the 1.0*ds integral so the
    "area" scalar is 1).
    """
    from dolfinx import fem
    from mpi4py import MPI

    if cond_form is None or area_form is None:
        return 0.0

    psi_pert.x.array[:] = x_3N[:N_psi]
    phi_n_pert.x.array[:] = x_3N[N_psi:N_psi + N_n]
    phi_p_pert.x.array[:] = x_3N[N_psi + N_n:]
    for fn in (psi_pert, phi_n_pert, phi_p_pert):
        fn.x.scatter_forward()

    I_local = fem.assemble_scalar(cond_form)
    A_local = fem.assemble_scalar(area_form)
    I_val = msh.comm.allreduce(I_local, op=MPI.SUM)
    A_val = msh.comm.allreduce(A_local, op=MPI.SUM)
    if A_val == 0.0:
        A_val = 1.0
    return float(I_val / A_val)


def _eval_displacement_charge(
    delta_psi_arr, psi_pert_fn, sc, eps_r,
    facet_info, msh,
) -> float:
    """Evaluate Q_psi = eps * integral(grad(delta_psi) . n) dS / area.

    The displacement-current contribution at frequency omega is
    j*omega*Q_psi (in A/m^2 when divided by the contact area, as is
    done here for consistency with the conduction-current evaluator).
    The FE FacetNormal at boundary facets is outward, so Q_psi is the
    contact-outward integral of eps*grad(delta_psi).
    """
    import ufl
    from dolfinx import fem
    from dolfinx.mesh import meshtags
    from mpi4py import MPI

    from ..constants import EPS0

    psi_pert_fn.x.array[:] = delta_psi_arr
    psi_pert_fn.x.scatter_forward()

    tag = int(facet_info["tag"])
    facets = facet_info["facets"]
    fdim = msh.topology.dim - 1
    if len(facets) == 0:
        return 0.0
    values = np.full(len(facets), tag, dtype=np.int32)
    indices = np.asarray(facets, dtype=np.int32)
    sort = np.argsort(indices)
    facet_mt = meshtags(msh, fdim, indices[sort], values[sort])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)
    n_vec = ufl.FacetNormal(msh)

    if isinstance(eps_r, (int, float)):
        eps_r_v = float(eps_r)
    else:
        eps_r_v = eps_r

    grad_dpsi_phys = sc.V0 * ufl.grad(psi_pert_fn)
    D_n_form = EPS0 * eps_r_v * ufl.dot(grad_dpsi_phys, n_vec)
    form_Q = fem.form(D_n_form * ds)
    form_A = fem.form(1.0 * ds)
    Q_local = fem.assemble_scalar(form_Q)
    A_local = fem.assemble_scalar(form_A)
    Q_val = msh.comm.allreduce(Q_local, op=MPI.SUM)
    A_val = msh.comm.allreduce(A_local, op=MPI.SUM)
    if A_val == 0.0:
        A_val = 1.0
    return float(Q_val / A_val)
