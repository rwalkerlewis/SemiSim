"""CLI entry point: semi-run <input.json> [--out runs/] [--run-id <id>]"""
from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        prog="semi-run",
        description="Run a SemiSim simulation and write a result artifact.",
    )
    parser.add_argument("input_json", type=Path, help="Path to input JSON file.")
    parser.add_argument(
        "--out", type=Path, default=Path("runs"),
        help="Output directory for run artifacts (default: runs/).",
    )
    parser.add_argument(
        "--run-id", type=str, default=None, dest="run_id",
        help="Override the auto-generated run ID.",
    )
    args = parser.parse_args()

    try:
        from semi.io.artifact import write_artifact
        from semi.run import run
        from semi.schema import load as schema_load

        cfg = schema_load(str(args.input_json))
        result = run(cfg)
        run_dir = write_artifact(
            result,
            out_dir=args.out,
            run_id=args.run_id,
            input_json_path=args.input_json,
        )
        print(run_dir)
    except Exception as exc:
        print(f"semi-run: error: {exc}", file=sys.stderr)
        sys.exit(1)
