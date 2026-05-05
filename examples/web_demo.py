"""
Web interface demo for kronos-semi.

Starts the kronos-server in the background, opens a browser tab pointing at
a self-contained HTML dashboard with a benchmark selector dropdown, then:

  1. User selects a benchmark from the dropdown
  2. Submits the chosen benchmark via POST /solve
  3. Streams SNES progress over the WebSocket (/runs/{id}/stream)
  4. Polls GET /runs/{id} until the run finishes and renders the result

Run from the repo root:

    python3 examples/web_demo.py

The server is shut down automatically when the script exits (Ctrl-C or
natural completion).

Requirements: the kronos-semi:ci Docker image must be present locally.
"""
from __future__ import annotations

import http
import http.client
import http.server
import json
import os
import subprocess
import sys
import textwrap
import threading
import time
import urllib.parse
import webbrowser
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"

SERVER_HOST = "127.0.0.1"
SERVER_PORT = 8765
DASHBOARD_PORT = 8766
STARTUP_TIMEOUT = 30


# ---------------------------------------------------------------------------
# Benchmark discovery
# ---------------------------------------------------------------------------

def _discover_benchmarks() -> dict[str, Path]:
    """Return {display_name: path} for every *.json found under benchmarks/."""
    found: dict[str, Path] = {}
    if not BENCHMARKS_DIR.is_dir():
        return found
    for p in sorted(BENCHMARKS_DIR.glob("*/*.json")):
        # Key is "subdir/filename" so names are unambiguous when a subdir has
        # multiple JSON files.
        key = f"{p.parent.name}/{p.name}"
        found[key] = p
    return found

# ---------------------------------------------------------------------------
# HTML dashboard (served by a tiny stdlib HTTP server)
# ---------------------------------------------------------------------------

DASHBOARD_HTML = textwrap.dedent(f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>kronos-semi live demo</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.3/dist/chart.umd.min.js"></script>
<style>
  body {{font-family:monospace;background:#111;color:#eee;padding:1.5rem;max-width:1100px;}}
  h1 {{color:#7dd3fc;}}
  h2 {{color:#a5f3fc;margin-top:1.5rem;}}
  label {{display:block;margin-bottom:.3rem;color:#94a3b8;}}
  select {{
    background:#1e293b;color:#eee;border:1px solid #334155;
    padding:.4rem .7rem;border-radius:4px;font-size:1rem;min-width:320px;
  }}
  #status {{color:#fde68a;font-size:1.1rem;margin:.6rem 0;}}
  #log {{
    background:#1e293b;border:1px solid #334155;padding:1rem;
    height:280px;overflow-y:auto;white-space:pre-wrap;font-size:.85rem;
  }}
  table {{border-collapse:collapse;margin-top:.5rem;}}
  th,td {{border:1px solid #334155;padding:.35rem .7rem;text-align:right;}}
  th {{background:#1e293b;color:#7dd3fc;}}
  tr:nth-child(even) td {{background:#1e2a3a;}}
  .row {{display:flex;align-items:flex-end;gap:1rem;flex-wrap:wrap;margin-top:.5rem;}}
  #submit-btn {{
    background:#0ea5e9;color:#fff;border:none;padding:.5rem 1.2rem;
    border-radius:4px;cursor:pointer;font-size:1rem;
  }}
  #submit-btn:disabled {{background:#334155;cursor:default;}}
  #json-preview {{
    background:#1e293b;border:1px solid #334155;padding:.8rem;
    max-height:160px;overflow-y:auto;white-space:pre-wrap;font-size:.78rem;
    color:#94a3b8;margin-top:.5rem;display:none;
  }}
  #toggle-preview {{
    background:none;border:none;color:#7dd3fc;cursor:pointer;
    font-size:.85rem;text-decoration:underline;padding:0;margin-top:.4rem;
  }}
  .chart-grid {{display:grid;grid-template-columns:repeat(auto-fill,minmax(440px,1fr));gap:1.2rem;margin-top:.8rem;}}
  .chart-card {{background:#1e293b;border:1px solid #334155;border-radius:6px;padding:1rem;}}
  .chart-title {{color:#7dd3fc;font-size:.9rem;margin-bottom:.5rem;}}
  canvas {{max-height:280px;}}
  .mesh-card {{background:#1e293b;border:1px solid #334155;border-radius:6px;padding:1rem;margin-top:.8rem;}}
  .mesh-canvas {{width:100%;height:320px;display:block;cursor:crosshair;}}
</style>
</head>
<body>
<h1>kronos-semi &mdash; live web interface demo</h1>
<p>
  Server: <a href="http://{SERVER_HOST}:{SERVER_PORT}/docs" target="_blank">
    http://{SERVER_HOST}:{SERVER_PORT}/docs
  </a>
</p>

<h2>Select benchmark</h2>
<label for="bench-select">Benchmark JSON</label>
<div class="row">
  <select id="bench-select" onchange="onBenchmarkChange()">
    <option value="">-- loading... --</option>
  </select>
  <button id="submit-btn" onclick="submitJob()" disabled>Submit job</button>
</div>
<button id="toggle-preview" onclick="togglePreview()" style="display:none">
  show JSON preview
</button>
<div id="json-preview"></div>

<p id="status">&nbsp;</p>

<h2>SNES progress stream</h2>
<div id="log">(progress will appear here after submission)</div>

<h2>Mesh</h2>
<div id="mesh-container"></div>

<h2>Results</h2>
<div id="result-container"></div>

<script>
const DASH = "http://127.0.0.1:{DASHBOARD_PORT}";
const BASE = "http://{SERVER_HOST}:{SERVER_PORT}";
const WS_BASE = "ws://{SERVER_HOST}:{SERVER_PORT}";
let currentCfg = null;
let runId = null;
let pollTimer = null;
const _charts = {{}};

function setStatus(msg, color="#fde68a") {{
  const el = document.getElementById("status");
  el.textContent = msg;
  el.style.color = color;
}}

function appendLog(msg) {{
  const el = document.getElementById("log");
  el.textContent += msg + "\\n";
  el.scrollTop = el.scrollHeight;
}}

function resetSubmit() {{
  document.getElementById("submit-btn").disabled = false;
}}

async function loadBenchmarkList() {{
  try {{
    const resp = await fetch(DASH + "/benchmarks");
    if (!resp.ok) {{ setStatus("Could not load benchmark list (HTTP " + resp.status + ")", "#f87171"); return; }}
    const names = await resp.json();
    const sel = document.getElementById("bench-select");
    sel.innerHTML = '<option value="">-- select a benchmark --</option>';
    for (const name of names) {{
      const opt = document.createElement("option");
      opt.value = name;
      opt.textContent = name;
      sel.appendChild(opt);
    }}
    document.getElementById("submit-btn").disabled = false;
    setStatus("Select a benchmark and click Submit.", "#94a3b8");
  }} catch (err) {{
    setStatus("Could not reach dashboard server: " + err.message, "#f87171");
  }}
}}

async function onBenchmarkChange() {{
  const name = document.getElementById("bench-select").value;
  currentCfg = null;
  document.getElementById("json-preview").style.display = "none";
  document.getElementById("toggle-preview").style.display = "none";
  if (!name) return;
  try {{
    const resp = await fetch(DASH + "/benchmarks/" + encodeURIComponent(name));
    if (!resp.ok) {{ setStatus("Could not load " + name + " (HTTP " + resp.status + ")", "#f87171"); return; }}
    currentCfg = await resp.json();
    document.getElementById("json-preview").textContent = JSON.stringify(currentCfg, null, 2);
    document.getElementById("toggle-preview").style.display = "inline";
    setStatus("Ready to submit: " + name, "#94a3b8");
  }} catch (err) {{
    setStatus("Error loading benchmark: " + err.message, "#f87171");
  }}
}}

function togglePreview() {{
  const el = document.getElementById("json-preview");
  const btn = document.getElementById("toggle-preview");
  const visible = el.style.display !== "none";
  el.style.display = visible ? "none" : "block";
  btn.textContent = visible ? "show JSON preview" : "hide JSON preview";
}}

async function submitJob() {{
  if (!currentCfg) {{
    setStatus("Please select a benchmark first.", "#f87171");
    return;
  }}
  document.getElementById("submit-btn").disabled = true;
  document.getElementById("log").textContent = "";
  document.getElementById("result-container").innerHTML = "";
  document.getElementById("mesh-container").innerHTML = "";
  setStatus("Submitting job...");
  appendLog("[submit] POST " + BASE + "/solve");

  let resp;
  try {{
    resp = await fetch(BASE + "/solve", {{
      method: "POST",
      headers: {{"Content-Type": "application/json"}},
      body: JSON.stringify(currentCfg),
    }});
  }} catch (err) {{
    setStatus("Network error (CORS or connection refused): " + err.message, "#f87171");
    appendLog("[error] " + err.message);
    resetSubmit();
    return;
  }}

  if (!resp.ok) {{
    let detail = resp.statusText;
    try {{ const body = await resp.json(); detail = JSON.stringify(body); }} catch(_) {{}}
    setStatus("Submit failed (HTTP " + resp.status + "): " + detail, "#f87171");
    appendLog("[error] HTTP " + resp.status + " " + detail);
    resetSubmit();
    return;
  }}

  let data;
  try {{ data = await resp.json(); }} catch (err) {{
    setStatus("Bad response from server: " + err.message, "#f87171");
    resetSubmit();
    return;
  }}

  runId = data.run_id;
  setStatus("Job queued: " + runId);
  appendLog("[queued] run_id = " + runId);
  openStream(runId);
  startPoll(runId);
}}

function openStream(id) {{
  try {{
    const ws = new WebSocket(WS_BASE + "/runs/" + id + "/stream");
    ws.onopen = () => appendLog("[ws] connected");
    ws.onmessage = (e) => {{
      try {{
        const evt = JSON.parse(e.data);
        if (evt.type === "step_done") {{
          appendLog(
            "  step " + String(evt.bias_step).padStart(2) +
            "  V=" + evt.V_applied.toFixed(4) +
            "  iters=" + evt.iterations
          );
        }} else if (evt.type === "run_done") {{
          const ok = evt.status === "completed";
          setStatus("Run finished: " + evt.status, ok ? "#4ade80" : "#f87171");
          stopPoll();
          ws.close();
          resetSubmit();
          fetchResults(id);
        }} else {{
          appendLog("[" + (evt.type || "event") + "] " + e.data);
        }}
      }} catch(_) {{
        appendLog(e.data);
      }}
    }};
    ws.onerror = (e) => appendLog("[ws] error (will fall back to polling)");
    ws.onclose = () => appendLog("[ws] closed");
  }} catch (err) {{
    appendLog("[ws] could not open: " + err.message);
  }}
}}

function startPoll(id) {{
  stopPoll();
  pollTimer = setInterval(async () => {{
    try {{
      const resp = await fetch(BASE + "/runs/" + id);
      if (!resp.ok) return;
      const m = await resp.json();
      const st = m.status;
      if (st === "completed" || st === "failed" || st === "error") {{
        stopPoll();
        setStatus("Run finished: " + st, st === "completed" ? "#4ade80" : "#f87171");
        resetSubmit();
        renderMeshDiagram(id);
        renderManifest(id, m);
      }} else {{
        setStatus("Running... (status: " + st + ")");
      }}
    }} catch(_) {{}}
  }}, 2000);
}}

function stopPoll() {{
  if (pollTimer) {{ clearInterval(pollTimer); pollTimer = null; }}
}}

async function fetchResults(id) {{
  try {{
    const resp = await fetch(BASE + "/runs/" + id);
    if (!resp.ok) return;
    renderMeshDiagram(id);
    renderManifest(id, await resp.json());
  }} catch(_) {{}}
}}

async function renderMeshDiagram(id) {{
  const container = document.getElementById("mesh-container");
  container.innerHTML = "";
  try {{
    const resp = await fetch(BASE + "/runs/" + id + "/mesh?max_nodes=5000");
    if (!resp.ok) return;
    const data = await resp.json();
    const card = document.createElement("div");
    card.className = "mesh-card";
    const title = document.createElement("div");
    title.className = "chart-title";
    const nodeCount = data.nodes.length;
    const cellCount = data.cells.length;
    title.textContent = "Mesh — " + nodeCount + " nodes, " + cellCount + " elements"
      + (data.downsampled ? " (downsampled)" : "");
    card.appendChild(title);
    const cvs = document.createElement("canvas");
    cvs.className = "mesh-canvas";
    card.appendChild(cvs);
    container.appendChild(card);
    _drawMesh(cvs, data);
  }} catch(_) {{}}
}}

function _drawMesh(cvs, data) {{
  const nodes = data.nodes;   // [[x,y]] or [[x]] for 1-D
  const cells = data.cells;   // [[i,j,...]]
  const dim = data.dim;

  // Fit canvas to CSS size
  const rect = cvs.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  cvs.width  = Math.round(rect.width  * dpr) || 800;
  cvs.height = Math.round(rect.height * dpr) || 320;
  const ctx = cvs.getContext("2d");
  ctx.scale(dpr, dpr);
  const W = cvs.width / dpr;
  const H = cvs.height / dpr;

  // Determine coordinate range
  const pad = 28;
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (const n of nodes) {{
    if (n[0] < xMin) xMin = n[0];
    if (n[0] > xMax) xMax = n[0];
    if (dim >= 2) {{
      if (n[1] < yMin) yMin = n[1];
      if (n[1] > yMax) yMax = n[1];
    }}
  }}
  if (dim < 2) {{ yMin = -0.5; yMax = 0.5; }}
  const xSpan = xMax - xMin || 1;
  const ySpan = yMax - yMin || 1;
  // Aspect-ratio-preserving scale
  const scaleX = (W - 2*pad) / xSpan;
  const scaleY = (H - 2*pad) / ySpan;
  const scale = Math.min(scaleX, scaleY);
  const offX = pad + ((W - 2*pad) - xSpan * scale) / 2;
  const offY = pad + ((H - 2*pad) - ySpan * scale) / 2;
  const tx = x => offX + (x - xMin) * scale;
  const ty = y => H - offY - (y - yMin) * scale;

  // Background
  ctx.fillStyle = "#0f172a";
  ctx.fillRect(0, 0, W, H);

  // Draw edges/cells
  ctx.strokeStyle = "rgba(56,189,248,0.35)";
  ctx.lineWidth = 0.6;
  if (dim === 1) {{
    // 1D: draw line segments along x-axis
    ctx.beginPath();
    for (const c of cells) {{
      ctx.moveTo(tx(nodes[c[0]][0]), ty(0));
      ctx.lineTo(tx(nodes[c[1]][0]), ty(0));
    }}
    ctx.stroke();
  }} else {{
    // 2D/3D: draw cell outlines (triangles/quads/tets projected)
    for (const c of cells) {{
      ctx.beginPath();
      ctx.moveTo(tx(nodes[c[0]][0]), ty(nodes[c[0]][1]));
      for (let k = 1; k < c.length; k++) {{
        ctx.lineTo(tx(nodes[c[k]][0]), ty(nodes[c[k]][1]));
      }}
      ctx.closePath();
      ctx.stroke();
    }}
  }}

  // Draw nodes
  const nodeR = dim === 1 ? 2.5 : (nodes.length > 2000 ? 0 : 1.5);
  if (nodeR > 0) {{
    ctx.fillStyle = "#7dd3fc";
    for (const n of nodes) {{
      const nx = tx(n[0]);
      const ny = dim >= 2 ? ty(n[1]) : ty(0);
      ctx.beginPath();
      ctx.arc(nx, ny, nodeR, 0, 2 * Math.PI);
      ctx.fill();
    }}
  }}

  // Axis labels (physical units auto-scaled)
  const _fmt = v => {{
    const av = Math.abs(v);
    if (av === 0) return "0";
    if (av < 1e-6) return (v*1e9).toFixed(1)+"nm";
    if (av < 1e-3) return (v*1e6).toFixed(1)+"μm";
    if (av < 1)    return (v*1e3).toFixed(1)+"mm";
    return v.toFixed(3)+"m";
  }};
  ctx.fillStyle = "#94a3b8";
  ctx.font = "11px monospace";
  ctx.textAlign = "left";
  ctx.fillText(_fmt(xMin), pad, H - 6);
  ctx.textAlign = "right";
  ctx.fillText(_fmt(xMax), W - pad, H - 6);
  if (dim >= 2) {{
    ctx.save();
    ctx.translate(12, H - pad);
    ctx.rotate(-Math.PI/2);
    ctx.textAlign = "left";
    ctx.fillText(_fmt(yMin), 0, 0);
    ctx.restore();
    ctx.save();
    ctx.translate(12, pad);
    ctx.rotate(-Math.PI/2);
    ctx.textAlign = "right";
    ctx.fillText(_fmt(yMax), 0, 0);
    ctx.restore();
  }}
}}

function renderManifest(id, manifest) {{
  const container = document.getElementById("result-container");
  container.innerHTML = "";
  // Destroy old charts to prevent canvas reuse errors.
  for (const key of Object.keys(_charts)) {{
    _charts[key].destroy();
    delete _charts[key];
  }}

  // IV charts -- one per sweep contact listed in manifest.sweeps
  const sweeps = manifest.sweeps || [];
  const ivContacts = sweeps.map(s => s.contact).filter(Boolean);
  // Fallback guess if sweeps not listed.
  const contactsToTry = ivContacts.length
    ? [...new Set(ivContacts)]
    : ["anode", "cathode", "gate", "source", "drain"];

  const ivGrid = document.createElement("div");
  ivGrid.className = "chart-grid";
  container.appendChild(ivGrid);

  for (const contact of contactsToTry) {{
    fetchIVChart(id, contact, ivGrid);
  }}

  // Field profile charts (only for 1-D runs).
  const fields = (manifest.fields || []).filter(f => f.rank === 0);
  if (fields.length) {{
    const fieldGrid = document.createElement("div");
    fieldGrid.className = "chart-grid";
    container.appendChild(fieldGrid);
    for (const field of fields) {{
      fetchFieldChart(id, field.name, field.units || "", fieldGrid);
    }}
  }}
}}

async function fetchIVChart(id, contact, grid) {{
  try {{
    const resp = await fetch(BASE + "/runs/" + id + "/iv/" + contact);
    if (!resp.ok) return;
    const csv = await resp.text();
    const rows = parseCSV(csv);
    if (rows.length < 2) return;
    const headers = rows[0].map(h => h.toLowerCase().trim());

    // Detect AC sweep (f_Hz,C,G,...) vs IV (V,J_n,J_p,J_total) vs transient (t,V,...)
    const isAC = headers.includes("f_hz") || headers.includes("c");
    const isTrans = headers[0] === "t" && headers.includes("v");

    if (isAC) {{
      const fIdx = headers.findIndex(h => h === "f_hz");
      const cIdx = headers.findIndex(h => h === "c");
      const gIdx = headers.findIndex(h => h === "g");
      if (fIdx < 0 || cIdx < 0) return;
      const cData = rows.slice(1)
        .map(r => ({{ x: parseFloat(r[fIdx]), y: parseFloat(r[cIdx]) }}))
        .filter(p => isFinite(p.x) && isFinite(p.y));
      const gData = gIdx >= 0 ? rows.slice(1)
        .map(r => ({{ x: parseFloat(r[fIdx]), y: parseFloat(r[gIdx]) }}))
        .filter(p => isFinite(p.x) && isFinite(p.y)) : [];
      if (!cData.length) return;
      addChart(grid, contact + " C(f)", "f (Hz)", "C (F/m²)", [
        {{ label: "C", data: cData, borderColor: "#34d399", backgroundColor: "rgba(52,211,153,0.12)", pointRadius: 3, tension: 0.3 }},
        ...(gData.length ? [{{ label: "G", data: gData, borderColor: "#fb923c", backgroundColor: "rgba(251,146,60,0.12)", pointRadius: 3, tension: 0.3 }}] : []),
      ], true);
    }} else if (isTrans) {{
      const tIdx = 0;
      const jIdx = headers.findIndex(h => h.includes("total") || h === "j");
      if (jIdx < 0) return;
      const data = rows.slice(1)
        .map(r => ({{ x: parseFloat(r[tIdx]) * 1e9, y: Math.abs(parseFloat(r[jIdx])) }}))
        .filter(p => isFinite(p.x) && isFinite(p.y) && p.y > 0);
      if (!data.length) return;
      addChart(grid, contact + " I(t)", "t (ns)", "|J| (A/m²)", [{{
        label: contact, data, borderColor: "#a78bfa",
        backgroundColor: "rgba(167,139,250,0.12)", pointRadius: 2, tension: 0.3,
      }}], true);
    }} else {{
      const vIdx = headers.findIndex(h => h.startsWith("v"));
      const jIdx = headers.findIndex(h => h.includes("total") || h.startsWith("j"));
      if (vIdx < 0 || jIdx < 0) return;
      const data = rows.slice(1)
        .map(r => ({{ x: parseFloat(r[vIdx]), y: Math.abs(parseFloat(r[jIdx])) }}))
        .filter(p => isFinite(p.x) && isFinite(p.y) && p.y > 0);
      if (!data.length) return;
      addChart(grid, contact + " IV", "V (V)", "|J| (A/m²)", [{{
        label: contact, data, borderColor: "#38bdf8",
        backgroundColor: "rgba(56,189,248,0.12)", pointRadius: 3, tension: 0.3,
      }}], true);
    }}
  }} catch(_) {{}}
}}

async function fetchFieldChart(id, name, units, grid) {{
  try {{
    const resp = await fetch(BASE + "/runs/" + id + "/fields/" + name + "/profile");
    if (!resp.ok) return;
    const csv = await resp.text();
    const rows = parseCSV(csv);
    if (rows.length < 2) return;
    const data = rows.slice(1).map(r => ({{ x: parseFloat(r[0]) * 1e6, y: parseFloat(r[1]) }}));
    const colors = [
      "#34d399","#f472b6","#fb923c","#a78bfa","#facc15","#60a5fa","#f87171"
    ];
    const colorIdx = Object.keys(_charts).length % colors.length;
    addChart(grid, name + (units ? " (" + units + ")" : ""), "x (μm)", name + (units ? " [" + units + "]" : ""), [{{
      label: name,
      data,
      borderColor: colors[colorIdx],
      backgroundColor: colors[colorIdx].replace(")", ",0.1)").replace("rgb","rgba"),
      pointRadius: 0,
      tension: 0.2,
    }}], false);
  }} catch(_) {{}}
}}

function addChart(grid, title, xLabel, yLabel, datasets, logY) {{
  const card = document.createElement("div");
  card.className = "chart-card";
  const titleEl = document.createElement("div");
  titleEl.className = "chart-title";
  titleEl.textContent = title;
  card.appendChild(titleEl);
  const canvas = document.createElement("canvas");
  card.appendChild(canvas);
  grid.appendChild(card);
  const key = title;
  _charts[key] = new Chart(canvas, {{
    type: "line",
    data: {{ datasets }},
    options: {{
      animation: false,
      parsing: false,
      plugins: {{ legend: {{ labels: {{ color: "#cbd5e1" }} }} }},
      scales: {{
        x: {{
          type: "linear",
          title: {{ display: true, text: xLabel, color: "#94a3b8" }},
          ticks: {{ color: "#94a3b8" }},
          grid: {{ color: "#1e3a5f" }},
        }},
        y: {{
          type: logY ? "logarithmic" : "linear",
          title: {{ display: true, text: yLabel, color: "#94a3b8" }},
          ticks: {{ color: "#94a3b8" }},
          grid: {{ color: "#1e3a5f" }},
        }},
      }},
    }},
  }});
}}

function parseCSV(text) {{
  return text.trim().split("\\n").map(line => line.split(","));
}}

function renderIV(csv) {{
  // Legacy table fallback -- kept for non-chart contexts.
  const lines = csv.trim().split("\\n");
  const headers = lines[0].split(",");
  let html = "<table><tr>" + headers.map(h => "<th>" + h + "</th>").join("") + "</tr>";
  for (let i = 1; i < lines.length; i++) {{
    const cells = lines[i].split(",");
    html += "<tr>" + cells.map(c => "<td>" + parseFloat(c).toExponential(4) + "</td>").join("") + "</tr>";
  }}
  document.getElementById("result-container").innerHTML = html + "</table>";
}}

function renderIVFromJson(rows) {{
  if (!rows.length) return;
  const headers = Object.keys(rows[0]);
  let html = "<table><tr>" + headers.map(h => "<th>" + h + "</th>").join("") + "</tr>";
  for (const row of rows) {{
    html += "<tr>" + headers.map(h =>
      "<td>" + (typeof row[h] === "number" ? row[h].toExponential(4) : row[h]) + "</td>"
    ).join("") + "</tr>";
  }}
  document.getElementById("result-container").innerHTML = html + "</table>";
}}

window.addEventListener("DOMContentLoaded", loadBenchmarkList);
</script>
</body>
</html>
""")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _wait_for_server(host: str, port: int, timeout: float) -> bool:
    """Poll GET /health until 200 or timeout."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            conn = http.client.HTTPConnection(host, port, timeout=1)
            conn.request("GET", "/health")
            r = conn.getresponse()
            if r.status == 200:
                return True
        except Exception:
            pass
        time.sleep(0.4)
    return False


class _DashboardHandler(http.server.BaseHTTPRequestHandler):
    benchmarks: dict[str, Path] = {}
    html_bytes: bytes = b""

    def do_GET(self):  # noqa: N802
        path = urllib.parse.urlparse(self.path).path

        if path == "/":
            self._respond(200, "text/html; charset=utf-8", self.html_bytes)

        elif path == "/benchmarks":
            data = json.dumps(list(self.benchmarks.keys())).encode()
            self._respond(200, "application/json", data)

        elif path.startswith("/benchmarks/"):
            raw = path[len("/benchmarks/"):]
            name = urllib.parse.unquote(raw)
            if name not in self.benchmarks:
                self._respond(404, "application/json", b'{"detail":"not found"}')
                return
            try:
                bench_path = self.benchmarks[name]
                cfg = json.loads(bench_path.read_bytes())
                # Rewrite relative mesh paths to absolute Docker-side paths so
                # the server (running in Docker at /workspaces/kronos-semi) can
                # find them.  The repo root on the host maps to
                # /workspaces/kronos-semi inside the container.
                mesh = cfg.get("mesh", {})
                if mesh.get("source") == "file":
                    mp = mesh.get("path", "")
                    if mp and not mp.startswith("/"):
                        bench_dir = bench_path.parent
                        abs_host = (bench_dir / mp).resolve()
                        # Translate host absolute path → Docker absolute path.
                        try:
                            rel = abs_host.relative_to(REPO_ROOT)
                            cfg["mesh"]["path"] = f"/workspaces/kronos-semi/{rel}"
                        except ValueError:
                            pass  # path outside repo, leave as-is
                content = json.dumps(cfg).encode()
                self._respond(200, "application/json", content)
            except OSError:
                self._respond(500, "application/json", b'{"detail":"read error"}')

        else:
            self._respond(404, "text/plain", b"not found")

    def _respond(self, status: int, content_type: str, body: bytes) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *_):  # suppress access log noise
        pass


def _serve_dashboard(
    html_bytes: bytes, benchmarks: dict[str, Path], port: int
) -> http.server.HTTPServer:
    _DashboardHandler.html_bytes = html_bytes
    _DashboardHandler.benchmarks = benchmarks
    server = http.server.HTTPServer(("127.0.0.1", port), _DashboardHandler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("kronos-semi web interface demo")
    print("=" * 40)

    benchmarks = _discover_benchmarks()
    if not benchmarks:
        sys.exit(f"No benchmark JSON files found under {BENCHMARKS_DIR}")
    print(f"[demo] Found {len(benchmarks)} benchmark(s)")

    html_bytes = DASHBOARD_HTML.encode()

    _serve_dashboard(html_bytes, benchmarks, DASHBOARD_PORT)
    dashboard_url = f"http://127.0.0.1:{DASHBOARD_PORT}/"
    print(f"[demo] Dashboard served at {dashboard_url}")

    # Start kronos-server inside the Docker image so FEM + server deps are
    # available without any host-side installation.
    docker_cmd = [
        "docker", "run", "--rm",
        "-p", f"{SERVER_PORT}:{SERVER_PORT}",
        "-v", f"{REPO_ROOT}:/workspaces/kronos-semi",
        "-w", "/workspaces/kronos-semi",
        "-e", f"KRONOS_SERVER_HOST=0.0.0.0",
        "-e", f"KRONOS_SERVER_PORT={SERVER_PORT}",
        "-e", "KRONOS_SERVER_WORKERS=1",
        "-e", "KRONOS_SERVER_RUNS_DIR=/workspaces/kronos-semi/runs",        "-e", f"KRONOS_SERVER_CORS_ORIGINS=http://127.0.0.1:{DASHBOARD_PORT},http://localhost:{DASHBOARD_PORT}",        "-e", "MPLBACKEND=Agg",
        "kronos-semi:ci",
        "python", "-m", "uvicorn", "kronos_server.app:build_app",
        "--host", "0.0.0.0", "--port", str(SERVER_PORT),
        "--factory", "--log-level", "warning", "--reload",
    ]
    print(f"[demo] Starting server in Docker (port {SERVER_PORT}) ...")
    server_proc = subprocess.Popen(
        docker_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    threading.Thread(
        target=lambda: [_ for _ in server_proc.stdout], daemon=True
    ).start()

    print(f"[demo] Waiting for server on {SERVER_HOST}:{SERVER_PORT} ...", end="", flush=True)
    if not _wait_for_server(SERVER_HOST, SERVER_PORT, STARTUP_TIMEOUT):
        server_proc.terminate()
        sys.exit("\nServer did not start in time.")
    print(" ready")

    print(f"[demo] Opening browser: {dashboard_url}")
    webbrowser.open(dashboard_url)

    print()
    print("Select a benchmark from the dropdown and click Submit.")
    print("Press Ctrl-C to stop the server and exit.")
    print()

    try:
        server_proc.wait()
    except KeyboardInterrupt:
        print("\n[demo] Shutting down server...")
        server_proc.terminate()
        try:
            server_proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            server_proc.kill()
    finally:
        print("[demo] Done.")


if __name__ == "__main__":
    main()
