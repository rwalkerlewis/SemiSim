"""
Boundary-condition construction for kronos-semi.

This module is the single place that turns a validated JSON config plus a
mesh into the dolfinx `dirichletbc` objects consumed by the equilibrium
Poisson and coupled drift-diffusion solvers. It joins the pure-Python
core tier: dolfinx is imported lazily inside the functions that need it
so the module itself loads cleanly on the lint / pure-python CI matrix
where dolfinx is absent.

Public API:

    ContactBC                 pure-data description of one resolved contact
    resolve_contacts(cfg)     validate + resolve facet refs, no dolfinx
    build_psi_dirichlet_bcs   psi-only BCs (equilibrium Poisson path)
    build_dd_dirichlet_bcs    (psi, phi_n, phi_p) BCs (coupled DD path)

Both build functions operate on the resolved `list[ContactBC]` and
delegate facet-tag resolution / kind validation to `resolve_contacts`.
Bias-sweep drivers call `resolve_contacts` once per step with a
`voltages` override dict to push the next ramp value through the same
code path as the seed solve.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

# Schema-allowed contact kinds that produce a Dirichlet BC. "insulating"
# is a natural (Neumann) boundary and is therefore skipped during
# resolution; it never produces a ContactBC.
_VALID_DIRICHLET_KINDS = frozenset({"ohmic", "gate", "schottky"})


@dataclass
class ContactBC:
    """
    Pure-data description of one contact's BC, resolved from JSON but
    independent of dolfinx. The FEM layer consumes this plus the mesh
    to build dolfinx.fem.dirichletbc objects.

    `barrier_height_eV` is the metal-semiconductor barrier phi_B in
    eV (Sze 3rd ed Table 5; Pt-on-n-Si is 0.85 eV). It is populated
    only on Schottky contacts (M16.5); ohmic, gate, and insulating
    contacts leave it `None`.
    """
    name: str
    kind: str
    facet_tag: int
    V_applied: float
    work_function: float | None = None
    barrier_height_eV: float | None = None


def resolve_contacts(
    cfg: dict[str, Any],
    facet_tags: Any | None = None,
    voltages: dict[str, float] | None = None,
) -> list[ContactBC]:
    """
    Walk the contacts section of the JSON config and resolve every
    Dirichlet contact into a `ContactBC`.

    Parameters
    ----------
    cfg
        Validated JSON config dict (the same object that `semi.run.run`
        receives).
    facet_tags
        Optional dolfinx `MeshTags` for the facet entities. When
        provided, the resolved tag is verified to map to at least one
        facet in the mesh. Pure-Python callers (tests, schema
        validators) can pass `None` to skip the mesh check.
    voltages
        Optional mapping `contact name -> applied voltage` that
        overrides the per-contact `voltage` field. Bias-sweep drivers
        use this to push the next ramp value through the same resolver
        without mutating the config.

    Returns
    -------
    list[ContactBC]
        One entry per Dirichlet contact, in the order they appear in
        the config. Insulating contacts are skipped (natural BC).

    Raises
    ------
    ValueError
        If a contact's `type` is not one of the schema-allowed kinds.
    RuntimeError
        If a string `facet` reference does not resolve to a tag listed
        under `mesh.facets_by_plane`, or if the resolved tag has no
        facets in `facet_tags` (when `facet_tags` is not None).
    """
    tag_by_name = {
        p["name"]: int(p["tag"])
        for p in cfg["mesh"].get("facets_by_plane", [])
    }

    out: list[ContactBC] = []
    for contact in cfg["contacts"]:
        kind = contact["type"]
        name = contact["name"]
        if kind == "insulating":
            continue
        if kind not in _VALID_DIRICHLET_KINDS:
            raise ValueError(
                f"Unknown contact kind {kind!r} for contact {name!r}; "
                f"expected one of {sorted(_VALID_DIRICHLET_KINDS)} "
                f"or 'insulating'."
            )

        facet_ref = contact["facet"]
        if isinstance(facet_ref, str):
            if facet_ref not in tag_by_name:
                raise RuntimeError(
                    f"Contact {name!r} references facet name "
                    f"{facet_ref!r}, which is not declared under "
                    f"mesh.facets_by_plane."
                )
            tag = tag_by_name[facet_ref]
        else:
            tag = int(facet_ref)

        if facet_tags is not None:
            facets = facet_tags.find(tag)
            if len(facets) == 0:
                raise RuntimeError(
                    f"Contact {name!r} maps to facet tag {tag} but no "
                    f"facets carry that tag in the mesh."
                )

        if voltages is not None and name in voltages:
            V_applied = float(voltages[name])
        else:
            V_applied = float(contact.get("voltage", 0.0))

        wf = contact.get("workfunction")
        phi_b = contact.get("barrier_height_eV")
        out.append(ContactBC(
            name=name,
            kind=kind,
            facet_tag=tag,
            V_applied=V_applied,
            work_function=(float(wf) if wf is not None else None),
            barrier_height_eV=(
                float(phi_b) if kind == "schottky" and phi_b is not None
                else None
            ),
        ))
    return out


def _schottky_psi_eq(contact: ContactBC, ref_mat, sc) -> float:
    """
    Equilibrium psi at a Schottky facet in scaled units (psi / V_t).

    The metal Fermi level sits ``barrier_height_eV`` below the
    semiconductor conduction band on the metal side. In the kronos-semi
    Slotboom convention psi is measured from the intrinsic Fermi level
    (so n = n_i exp(psi/V_t) under Boltzmann), which gives
    ``q psi = E_F - E_i``. At the Schottky interface E_F = E_F_metal and
    E_C = phi_B, while E_C - E_i = V_t ln(N_C / n_i) (Boltzmann limit;
    Sze 3rd ed eq 1.20), so

        psi_Schottky_eq = V_t * ln(N_C / n_i) - phi_B
        psi_eq_hat      = psi_Schottky_eq / V_t
                        = ln(N_C / n_i) - phi_B / V_t.

    The expression assumes the Boltzmann limit; under Fermi-Dirac
    statistics the ``V_t ln(N_C / n_i)`` term is replaced by the FD
    inversion of the bulk equilibrium ``n_eq`` value. M16.5 ships under
    Boltzmann statistics by default, so the FD generalization is a
    follow-up. See Sze 3rd ed Section 3.4 and ADR 0015 for the
    derivation and V&V scope.
    """
    if contact.barrier_height_eV is None:
        raise ValueError(
            f"Schottky contact {contact.name!r} has no barrier_height_eV"
        )
    N_C = ref_mat.Nc
    if N_C is None or N_C <= 0.0:
        raise ValueError(
            f"Schottky psi BC requires reference material {ref_mat.name!r}"
            f" to expose a positive conduction-band density of states"
            f" N_C; got {N_C}"
        )
    return float(math.log(N_C / ref_mat.n_i) - contact.barrier_height_eV / sc.V0)


def build_psi_dirichlet_bcs(
    V_psi,
    msh,
    facet_tags,
    contacts: list[ContactBC],
    sc,
    ref_mat,
    N_raw_fn,
) -> list:
    """
    Build the Dirichlet BC list for psi.

    Ohmic contacts: psi_hat = asinh(N_net_local / (2 n_i)) + V_applied / V_t,
    where N_net_local is the doping evaluated at the contact facet centroid.

    Gate contacts: psi_hat = (V_gate - phi_ms) / V_t, where phi_ms is the
    metal-semiconductor work function difference carried on the ContactBC
    as `work_function` (JSON key `workfunction`). Gate contacts sit on the
    oxide side of the Si/SiO2 interface and contribute no Slotboom BCs;
    `build_dd_dirichlet_bcs` therefore ignores them on purpose.

    Schottky contacts: psi_hat = ln(N_C / n_i) - phi_B / V_t + V_applied / V_t.
    The metal Fermi level is the Dirichlet anchor for psi (Sze 3rd ed
    Section 3.4); the continuity rows (phi_n / phi_p) at the Schottky
    facet are handled by the thermionic-emission Robin form added in
    Phase C, not by a Dirichlet here. ADR 0015 documents the V&V scope.
    """
    from dolfinx import fem
    from petsc4py import PETSc

    fdim = msh.topology.dim - 1
    bcs = []
    for c in contacts:
        facets = facet_tags.find(c.facet_tag)
        if len(facets) == 0:
            raise RuntimeError(
                f"No facets with tag {c.facet_tag} for contact {c.name!r}"
            )
        dofs = fem.locate_dofs_topological(V_psi, fdim, facets)
        if c.kind == "ohmic":
            N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
            psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
            psi_bc = psi_eq_hat + c.V_applied / sc.V0
            bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs, V_psi))
        elif c.kind == "gate":
            phi_ms = float(c.work_function) if c.work_function is not None else 0.0
            psi_bc = (c.V_applied - phi_ms) / sc.V0
            bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs, V_psi))
        elif c.kind == "schottky":
            psi_bc = _schottky_psi_eq(c, ref_mat, sc) + c.V_applied / sc.V0
            bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs, V_psi))
        else:
            raise ValueError(f"Unknown contact kind {c.kind!r} for {c.name!r}")
    return bcs


def build_dd_dirichlet_bcs(
    spaces,
    msh,
    facet_tags,
    contacts: list[ContactBC],
    sc,
    ref_mat,
    N_raw_fn,
) -> list:
    """
    Build the Dirichlet BC list for the (psi, phi_n, phi_p) block
    system on each ohmic contact (the coupled drift-diffusion path).

    Per ohmic contact the scaled values are

        psi_hat   = asinh(N_net_local / (2 n_i)) + V_applied / V_t
        phi_n_hat = phi_p_hat = V_applied / V_t

    which is the Shockley boundary: quasi-Fermi potentials pinned to
    the applied bias, electrostatic potential offset by the local
    majority-side value at equilibrium.

    Schottky contacts (M16.5) contribute a single psi Dirichlet
    anchored at the metal Fermi level
    (`psi_hat = _schottky_psi_eq(c, ref_mat, sc) + V_applied / V_t`).
    The continuity rows (phi_n / phi_p) are not Dirichlet-pinned at the
    Schottky facet; the thermionic-emission Robin form on those rows is
    added in `semi.physics.drift_diffusion.build_dd_block_residual`.
    See ADR 0015.
    """
    from dolfinx import fem
    from petsc4py import PETSc

    fdim = msh.topology.dim - 1
    bcs = []
    for c in contacts:
        if c.kind not in ("ohmic", "gate", "schottky"):
            continue
        facets = facet_tags.find(c.facet_tag)
        if len(facets) == 0:
            raise RuntimeError(
                f"No facets with tag {c.facet_tag} for contact {c.name!r}"
            )
        dofs_psi = fem.locate_dofs_topological(spaces.V_psi, fdim, facets)

        if c.kind == "gate":
            # Gate contacts sit on the oxide side of the Si / SiO2
            # interface; the Slotboom continuity blocks assemble only
            # on the semiconductor submesh, so gate facets carry no
            # (phi_n, phi_p) DOFs. The psi Dirichlet must still be
            # applied here because the block solver consumes a single
            # bcs list (the alternative `build_psi_dirichlet_bcs`
            # entry point is used by the equilibrium and mos_cv
            # runners which solve psi alone).
            phi_ms = float(c.work_function) if c.work_function is not None else 0.0
            psi_bc = (c.V_applied - phi_ms) / sc.V0
            bcs.append(fem.dirichletbc(
                PETSc.ScalarType(psi_bc), dofs_psi, spaces.V_psi))
            continue

        if c.kind == "schottky":
            psi_bc = _schottky_psi_eq(c, ref_mat, sc) + c.V_applied / sc.V0
            bcs.append(fem.dirichletbc(
                PETSc.ScalarType(psi_bc), dofs_psi, spaces.V_psi))
            continue

        # Ohmic: Shockley boundary on the (psi, phi_n, phi_p) block.
        N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
        psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
        V_hat = c.V_applied / sc.V0
        psi_bc = psi_eq_hat + V_hat
        phi_bc = V_hat

        dofs_n = fem.locate_dofs_topological(spaces.V_phi_n, fdim, facets)
        dofs_p = fem.locate_dofs_topological(spaces.V_phi_p, fdim, facets)

        bcs.append(fem.dirichletbc(
            PETSc.ScalarType(psi_bc), dofs_psi, spaces.V_psi))
        bcs.append(fem.dirichletbc(
            PETSc.ScalarType(phi_bc), dofs_n, spaces.V_phi_n))
        bcs.append(fem.dirichletbc(
            PETSc.ScalarType(phi_bc), dofs_p, spaces.V_phi_p))
    return bcs


def _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn) -> float:
    """
    Evaluate the (raw, dimensional) net doping at the centroid of the
    first facet in `facets`. Used by both BC builders to set the
    majority-side equilibrium psi value at an ohmic contact.
    """
    tdim = msh.topology.dim
    msh.topology.create_connectivity(fdim, 0)
    f2v = msh.topology.connectivity(fdim, 0)
    coords = msh.geometry.x
    verts = f2v.links(int(facets[0]))
    centroid = coords[verts, :tdim].mean(axis=0)
    pt = centroid.reshape(tdim, 1)
    return float(N_raw_fn(pt)[0])
