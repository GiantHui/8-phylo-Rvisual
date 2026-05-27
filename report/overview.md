# Pipeline Overview

Generated: 2026-05-27

## Purpose

This directory contains the latest reusable ML-tree visualization pipeline. It takes three inputs:

1. A tree file.
2. A sample metadata table.
3. A Markdown visual-design specification.

The current example is configured for the Genghis Khan Y-chromosome project and focuses on the C2a1a3 target clade.

## Directory Layout

```text
8-phylo-Rvisual/
├── conf/
│   ├── Config.yaml
│   ├── VisualDesign.md
│   └── VisualDesign.template.md
├── pipe/
│   └── run_chengjisihan.sh
├── script/
│   └── load_config.sh
├── python/
│   └── prepare_annotations.py
├── src/
│   └── visualize_tree.R
├── output/
│   └── chengjisihan/
└── report/
    ├── overview.md
    └── report.md
```

## Current Run

```bash
bash pipe/run_chengjisihan.sh
```

The latest run completed successfully and generated:

- annotation TSV files in `output/chengjisihan/annotation/`
- full-tree figures in `output/chengjisihan/figures/full_tree/`
- target-subtree figures in `output/chengjisihan/figures/subtree/`
- this updated report directory

## Reuse Workflow

For a new tree visualization project:

1. Edit `conf/Config.yaml` to point to the new tree, metadata table, and design Markdown.
2. Copy `conf/VisualDesign.template.md` and update the scientific questions and `visual_design_config` block.
3. Run `bash pipe/run_chengjisihan.sh`.
4. Check `output/<project>/annotation/design_config.tsv` to confirm the design settings used by the run.

## Notes

The pipeline no longer relies on old exploratory R scripts. The only shell helper retained under `script/` is `load_config.sh`, which is used by the main control script.
