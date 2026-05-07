# kronos-semi examples catalogue

Each subdirectory is a self-contained practical device
configuration. Examples are illustrative, not V&V gates; for
analytical correctness gates, see [`benchmarks/`](../benchmarks).
The smoke verifier for each example asserts only that the run
completes, the output JSON is well-formed, and the recorded
I-V values are finite and non-NaN. Numerical correctness against
analytical references is benchmark territory.

## Examples shipped in v0.23.x

| Name | Device class | Runner | M16 features exercised |
|:-----|:-------------|:-------|:------------------------|
| [`nmos_idvgs`](nmos_idvgs/) | n-channel MOSFET | `bias_sweep` | M16.1 Caughey-Thomas, M16.2 Lombardi, M16.4 Fermi-Dirac |
| [`schottky_iv_temperature`](schottky_iv_temperature/) | Pt-on-n-Si Schottky diode | `bias_sweep` | M16.5 Schottky thermionic emission (under three temperatures) |
| [`power_diode_reverse_recovery`](power_diode_reverse_recovery/) | long-base pn rectifier | `transient` | M16.7 voltage_t (piecewise-linear), M16.3 Auger |

## How to run an example

Each example directory has a `README.md` with the exact
`docker compose run` invocation. The general form is:

```
docker compose run --rm benchmark <example_name>
```

The CLI in [`scripts/run_benchmark.py`](../scripts/run_benchmark.py)
looks under `benchmarks/` first and falls back to `examples/`,
so the same invocation works for both.

Output (plots and result JSON) lands in `results/<example_name>/`.

## How to add an example

1. Create a new subdirectory under `examples/` named
   `<device_class>_<test_type>/` in snake case.
2. Drop a JSON config in that directory using only schema
   features that already ship in v2.7.0; do not bump the
   schema for an example.
3. Write a `README.md` matching the load-bearing structure
   convention (see
   [`docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md`](../docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md)
   for the section ordering).
4. Register a smoke verifier `verify_<name>` in
   [`scripts/run_benchmark.py`](../scripts/run_benchmark.py).
   Smoke checks only: run completed, no NaN, qualitative
   ordering. No tight numerical gates.
5. Add a CI matrix entry in
   [`.github/workflows/ci.yml`](../.github/workflows/ci.yml)
   under `docker-fem-benchmarks`. Do not set
   `allow-failure: "true"`.
6. Add a one-line registration test in
   [`tests/test_examples_register.py`](../tests/test_examples_register.py)
   so the verifier registration is covered in the gated suite.

Each new example is its own PR on a `dev/examples-<name>`
branch and does not require a milestone number.

## Future examples (not in the v0.23.x catalogue)

These device classes are deliberately deferred to follow-up PRs:

- **NPN BJT Gummel plot.** Awaits 3-terminal infrastructure
  maturity (resistor_3d covers the ohmic-three-terminal case
  but the BJT base-injection geometry needs its own attention).
- **NMOS C-V at multiple body biases.** Companion to
  `nmos_idvgs` once the body-bias parametric sweep mechanism
  is comfortable in the bias_sweep runner.
- **Tunnel diode forward I-V.** Awaits independent validation
  of the M16.6 Kane formula in the forward-bias regime;
  reverse-bias breakdown (zener_1d) is the only validated
  BBT case today.
- **SiC Schottky diode.** Awaits SiC material parameters in
  `semi/materials.py`.
