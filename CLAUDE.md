# Project notes

## GWF pipeline

This project uses [GWF](https://gwf.app/) to orchestrate analyses on Slurm and locally. The pipeline lives in `workflow.py`.

For a compact local reference to the GWF API, CLI, and common patterns used here, see [`claude-gwf-ref.md`](./claude-gwf-ref.md). Consult it before writing or modifying targets, templates, `gwf.map`/`collect` chains, or backend/executor configuration.

## iTRAILS

The pipeline runs [iTRAILS](https://itrails.readthedocs.io/en/docs-stable/) (TRAILS coalescent HMM). For a local reference to its CLI, config format, model parameters, and exact output file names, see [`claude-itrails-ref.md`](./claude-itrails-ref.md). Consult it before touching the templates in `itrails_workflow.py`. iTRAILS runs in the isolated pixi environment `itrails-tool` (see `pixi.toml`); specs invoke it as `pixi run -e itrails-tool itrails-...`.
