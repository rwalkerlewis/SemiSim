"""
Insert "Geometry & Mesh" visualization sections into the Colab notebooks
01-04. Each notebook becomes a self-disclosing FEM walkthrough: the device
geometry, region tags, contact facets, and FEM mesh nodes are rendered
*before* the solver runs.

Idempotent: detects whether the geometry section already exists (by a
marker comment in the inserted code cell) and skips re-insertion.

Notebook 05 already has its own 3-panel mesh visualization (Si/SiO2 region
shading, gate facet, axisymmetric metric), so it is not modified here.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

MARKER = "# SemiSim: geometry+mesh visualization"


# ----------------------------------------------------------------------
# Reusable code blocks
# ----------------------------------------------------------------------

CODE_1D = MARKER + r"""
# Build the FEM mesh from the JSON config and visualize it alongside the
# doping profile and ohmic contact locations. This runs *before* the solve
# so the geometry being simulated is unambiguous.
import numpy as np
import matplotlib.pyplot as plt

from semi.mesh import build_mesh
from semi.doping import build_profile

msh, cell_tags, facet_tags = build_mesh(cfg)

x_nodes = msh.geometry.x[:, 0]
order = np.argsort(x_nodes)
xs = x_nodes[order]

# Net doping N_D - N_A on each mesh node (m^-3 from the JSON, plotted in cm^-3)
N_net = build_profile(cfg["doping"])(xs[None, :])

# Locate ohmic contacts on the x-axis from the facets_by_plane block
plane_by_name = {p["name"]: p for p in cfg["mesh"]["facets_by_plane"]}

fig, axes = plt.subplots(
    2, 1, figsize=(8.5, 4.6), sharex=True,
    gridspec_kw={"height_ratios": [3, 1]},
)
ax = axes[0]
ax.plot(xs * 1e6, N_net * 1e-6, "k-", lw=1.6)
ax.axhline(0.0, color="gray", lw=0.5)
ax.set_yscale("symlog", linthresh=1e10)
ax.set_ylabel(r"net doping $N_D - N_A$ (cm$^{-3}$)")
ax.set_title(
    f"{cfg['name']}: 1D geometry, doping, and FEM mesh"
    f"  ({xs.size} nodes, {xs.size - 1} elements)"
)
ax.grid(alpha=0.3)

ymax = ax.get_ylim()[1]
for c in cfg["contacts"]:
    plane = plane_by_name[c["facet"]]
    xc_um = plane["value"] * 1e6
    ax.axvline(xc_um, color="C3", lw=1.2, ls="--", alpha=0.8)
    ax.text(xc_um, ymax, f"  {c['name']} ({c['type']})",
            color="C3", va="top", ha="left", fontsize=9)

# Mesh-node strip
ax2 = axes[1]
ax2.plot(xs * 1e6, np.zeros_like(xs), "|", color="C0", ms=10, mew=1.0)
ax2.set_yticks([])
ax2.set_xlabel(r"$x$ ($\mu$m)")
ax2.set_title(f"FEM mesh nodes (interval mesh, {xs.size} dofs)")
plt.tight_layout()
plt.show()
"""


CODE_2D_MOS = MARKER + r"""
# Build the FEM mesh from the JSON config and visualize the multi-region
# 2D layout (Si vs SiO2), the FEM triangulation, and the gate / body
# contact facets, *before* the solve.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Patch

from semi.mesh import build_mesh

msh, cell_tags, facet_tags = build_mesh(cfg)

# Mesh geometry
X = msh.geometry.x
tdim = msh.topology.dim
msh.topology.create_connectivity(tdim, 0)
c2v = msh.topology.connectivity(tdim, 0)
n_cells = msh.topology.index_map(tdim).size_local
cells_np = np.array([c2v.links(c) for c in range(n_cells)])
tri = Triangulation(X[:, 0] * 1e6, X[:, 1] * 1e9, triangles=cells_np)

# Region colours: Si below the interface, SiO2 above
region_color = {1: "#dfe9f5", 2: "#f7e3c0"}  # silicon, oxide
cell_color_idx = np.zeros(n_cells, dtype=int)
if cell_tags is not None:
    tag_arr = np.zeros(n_cells, dtype=int)
    tag_arr[cell_tags.indices] = cell_tags.values
    cell_color_idx = tag_arr

fig, ax = plt.subplots(figsize=(7.2, 5.4))
# Shade by region tag using tripcolor with face colors
tag_to_facecolor = np.array([
    region_color.get(int(t), "#e0e0e0") for t in cell_color_idx
])
# Use tripcolor with shading='flat' and a categorical colormap
unique_tags = sorted(set(int(t) for t in cell_color_idx))
tag_to_idx = {t: i for i, t in enumerate(unique_tags)}
facevals = np.array([tag_to_idx[int(t)] for t in cell_color_idx], dtype=float)
cmap = plt.matplotlib.colors.ListedColormap(
    [region_color.get(t, "#e0e0e0") for t in unique_tags]
)
ax.tripcolor(tri, facecolors=facevals, cmap=cmap, edgecolors="k",
             linewidth=0.15, shading="flat")

# Contact facets
plane_by_name = {p["name"]: p for p in cfg["mesh"]["facets_by_plane"]}
ext_um = cfg["mesh"]["extents"]
x_lo, x_hi = ext_um[0]
contact_styles = {"gate": ("C3", "-", 3.0), "ohmic": ("C0", "-", 3.0)}
for c in cfg["contacts"]:
    plane = plane_by_name[c["facet"]]
    color, ls, lw = contact_styles.get(c["type"], ("C2", "-", 2.0))
    if plane["axis"] == 1:
        y_nm = plane["value"] * 1e9
        ax.plot([x_lo * 1e6, x_hi * 1e6], [y_nm, y_nm],
                color=color, ls=ls, lw=lw, label=f"{c['name']} ({c['type']})")
    else:
        x_um = plane["value"] * 1e6
        y_lo, y_hi = ext_um[1]
        ax.plot([x_um, x_um], [y_lo * 1e9, y_hi * 1e9],
                color=color, ls=ls, lw=lw, label=f"{c['name']} ({c['type']})")

# Region legend patches
region_patches = [
    Patch(facecolor=region_color[t], edgecolor="k",
          label=f"region tag {t}: " + next(
              n for n, r in cfg["regions"].items() if r["tag"] == t))
    for t in unique_tags if t in region_color
]
contact_handles, contact_labels = ax.get_legend_handles_labels()
ax.legend(handles=region_patches + contact_handles,
          loc="center left", bbox_to_anchor=(1.02, 0.5),
          fontsize=9, frameon=False)

ax.set_xlabel(r"$x$ ($\mu$m)")
ax.set_ylabel(r"$y$ (nm)")
ax.set_title(
    f"{cfg['name']}: 2D geometry, regions, contacts, and FEM mesh"
    f"\n({X.shape[0]} nodes, {n_cells} triangles)"
)
ax.set_aspect("auto")
plt.tight_layout()
plt.show()
"""


def code_3d(cfg_var: str) -> str:
    """Return a self-contained 3D-mesh visualization cell for the given
    config variable name (e.g. 'cfg_builtin' or 'cfg_gmsh')."""
    return MARKER + rf"""
# Build the FEM mesh from the JSON config and render the 3D bounding
# surface coloured by facet tag (left contact, right contact, sidewalls)
# *before* the solve. This makes the V-I geometry unambiguous and lets
# you compare the builtin box mesh against the gmsh fixture below.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from dolfinx import mesh as dmesh

from semi.mesh import build_mesh

_cfg = {cfg_var}
_label = {cfg_var!r}

msh, cell_tags, facet_tags = build_mesh(_cfg)

# Boundary triangles
tdim = msh.topology.dim
fdim = tdim - 1
msh.topology.create_connectivity(fdim, tdim)
msh.topology.create_connectivity(fdim, 0)
boundary_facets = dmesh.exterior_facet_indices(msh.topology)
f2v = msh.topology.connectivity(fdim, 0)
X = msh.geometry.x

tris = np.stack([X[f2v.links(int(f))] for f in boundary_facets])

# Tag map for boundary facets (left=10, right=20, anything else = sidewall)
tag_lookup = {{}}
if facet_tags is not None:
    for idx, val in zip(facet_tags.indices, facet_tags.values):
        tag_lookup[int(idx)] = int(val)

color_left = "#1f77b4"
color_right = "#d62728"
color_side = "#cccccc"
colors = []
for f in boundary_facets:
    t = tag_lookup.get(int(f), 0)
    if t == 10:
        colors.append(color_left)
    elif t == 20:
        colors.append(color_right)
    else:
        colors.append(color_side)

fig = plt.figure(figsize=(9.5, 4.6))
ax = fig.add_subplot(111, projection="3d")
pc = Poly3DCollection(tris * 1e6, facecolors=colors, edgecolors="k",
                      linewidths=0.15, alpha=0.85)
ax.add_collection3d(pc)

xmin, xmax = X[:, 0].min() * 1e6, X[:, 0].max() * 1e6
ymin, ymax = X[:, 1].min() * 1e6, X[:, 1].max() * 1e6
zmin, zmax = X[:, 2].min() * 1e6, X[:, 2].max() * 1e6
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
ax.set_box_aspect((xmax - xmin, ymax - ymin, zmax - zmin))
# View angle chosen so both the left (tag 10) and right (tag 20) contact
# faces are visible at the ends of the long axis.
ax.view_init(elev=22, azim=-130)
ax.set_xlabel(r"$x$ ($\mu$m)")
ax.set_ylabel(r"$y$ ($\mu$m)")
ax.set_zlabel(r"$z$ ($\mu$m)")

n_cells = msh.topology.index_map(tdim).size_local
src = _cfg["mesh"]["source"]
ax.set_title(
    f"{{_label}}: 3D mesh ({{src}} source) - "
    f"{{X.shape[0]}} nodes, {{n_cells}} tetrahedra, "
    f"{{tris.shape[0]}} boundary triangles\n"
    f"left contact (tag 10) blue, right contact (tag 20) red, sidewalls grey"
)
plt.tight_layout()
plt.show()
"""


# ----------------------------------------------------------------------
# Insertion logic
# ----------------------------------------------------------------------

def _has_marker(nb) -> bool:
    return any(MARKER in c.get("source", "") for c in nb.cells)


def _insert_after(nb, idx: int, *cells) -> None:
    for offset, cell in enumerate(cells, start=1):
        nb.cells.insert(idx + offset, cell)


def patch_notebook_01(path: Path) -> None:
    nb = nbf.read(path, as_version=4)
    if _has_marker(nb):
        print(f"  [skip] {path.name}: already has geometry/mesh cells")
        return

    md = nbf.v4.new_markdown_cell(
        "## 3.5. Geometry, doping, and FEM mesh\n\n"
        "Before the solve runs, build the dolfinx mesh from the JSON config "
        "and inspect what is being simulated: the doping profile $N_D - N_A$, "
        "the ohmic contact locations from `facets_by_plane`, and every node "
        "in the 1D interval mesh. If the geometry on this plot does not match "
        "what you expect, the JSON is wrong and there is no point running the "
        "solver."
    )
    code = nbf.v4.new_code_cell(CODE_1D)
    _insert_after(nb, 7, md, code)
    nbf.write(nb, path)
    print(f"  [ok]   {path.name}: inserted geometry/mesh after cell 7")


def patch_notebook_02(path: Path) -> None:
    nb = nbf.read(path, as_version=4)
    if _has_marker(nb):
        print(f"  [skip] {path.name}: already has geometry/mesh cells")
        return

    md = nbf.v4.new_markdown_cell(
        "### Geometry, doping, and FEM mesh\n\n"
        "Same `build_mesh + build_profile` introspection as Notebook 01: the "
        "device under test is the same symmetric pn junction, so we render "
        "it once before the forward sweep and re-use it for the reverse "
        "sweep below."
    )
    # Notebook 02 uses cfg_fwd, not cfg → patch the source string
    code_src = CODE_1D.replace("(cfg)", "(cfg_fwd)") \
                      .replace("cfg['name']", "cfg_fwd['name']") \
                      .replace('cfg["doping"]', 'cfg_fwd["doping"]') \
                      .replace('cfg["mesh"]', 'cfg_fwd["mesh"]') \
                      .replace('cfg["contacts"]', 'cfg_fwd["contacts"]')
    code = nbf.v4.new_code_cell(code_src)
    _insert_after(nb, 7, md, code)
    nbf.write(nb, path)
    print(f"  [ok]   {path.name}: inserted geometry/mesh after cell 7")


def patch_notebook_03(path: Path) -> None:
    nb = nbf.read(path, as_version=4)
    if _has_marker(nb):
        print(f"  [skip] {path.name}: already has geometry/mesh cells")
        return

    md = nbf.v4.new_markdown_cell(
        "### Geometry, regions, contacts, and FEM mesh\n\n"
        "Before any solve, build the dolfinx mesh and visualize the 2D layout: "
        "the silicon body (region tag 1, 0–500 nm), the 5 nm SiO$_2$ "
        "layer on top (region tag 2), the body and gate contact facets, "
        "and every triangle in the FEM mesh. The aspect ratio is auto, not "
        "1:1, because the oxide is two orders of magnitude thinner than the "
        "silicon — the structured mesh is biased to put most of its 2020 "
        "triangles inside that 5 nm sliver."
    )
    code = nbf.v4.new_code_cell(CODE_2D_MOS)
    _insert_after(nb, 7, md, code)
    nbf.write(nb, path)
    print(f"  [ok]   {path.name}: inserted geometry/mesh after cell 7")


def patch_notebook_04(path: Path) -> None:
    nb = nbf.read(path, as_version=4)
    if _has_marker(nb):
        print(f"  [skip] {path.name}: already has geometry/mesh cells")
        return

    # Original layout: cell 7 = builtin load, cell 12 = gmsh load.
    # Insert after gmsh first (so builtin insertion does not shift its index).

    md_gmsh = nbf.v4.new_markdown_cell(
        "### Geometry and gmsh-fixture mesh\n\n"
        "Same 1 µm × 200 nm × 200 nm device, but the mesh is now an "
        "unstructured tet mesh from the gmsh fixture in "
        "`benchmarks/resistor_3d/`. Boundary triangles are coloured by "
        "facet tag: left contact (tag 10), right contact (tag 20), "
        "sidewalls (untagged, grey)."
    )
    code_gmsh = nbf.v4.new_code_cell(code_3d("cfg_gmsh"))
    _insert_after(nb, 12, md_gmsh, code_gmsh)

    md_builtin = nbf.v4.new_markdown_cell(
        "### Geometry and builtin box mesh\n\n"
        "Before the solve, build the dolfinx mesh from the JSON config and "
        "render the 3D bounding surface. The structured `box` source uses "
        "extents/resolution from the JSON (here 64 × 16 × 16), and the "
        "left/right contact facets at $x = 0$ and $x = L$ are tagged 10 "
        "and 20 respectively (`facets_by_plane`). Sidewalls are no-flux."
    )
    code_builtin = nbf.v4.new_code_cell(code_3d("cfg_builtin"))
    _insert_after(nb, 7, md_builtin, code_builtin)

    nbf.write(nb, path)
    print(f"  [ok]   {path.name}: inserted geometry/mesh for builtin and gmsh")


def main() -> None:
    here = Path(__file__).resolve().parent.parent
    nb_dir = here / "notebooks"
    print(f"Patching notebooks in {nb_dir}")
    patch_notebook_01(nb_dir / "01_pn_junction_1d.ipynb")
    patch_notebook_02(nb_dir / "02_pn_junction_bias.ipynb")
    patch_notebook_03(nb_dir / "03_mos_cv.ipynb")
    patch_notebook_04(nb_dir / "04_resistor_3d.ipynb")


if __name__ == "__main__":
    main()
