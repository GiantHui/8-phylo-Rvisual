# Genghis Khan Y-Chromosome ML Tree Visualization Report

Generated: 2026-05-27

## Pipeline Status

The current workflow is a three-input, reusable tree-visualization pipeline:

1. Tree file: `inputs.tree_file` in `conf/Config.yaml`.
2. Sample metadata table: `inputs.meta_file` in `conf/Config.yaml`.
3. Visual design Markdown: `inputs.design_file` in `conf/Config.yaml`.

The design file is `conf/VisualDesign.md`. It stores the scientific questions, visual-design principles, and a machine-readable `visual_design_config` block. The parsed settings are copied to `output/chengjisihan/annotation/design_config.tsv`, and the source design document is copied to `output/chengjisihan/annotation/visual_design.md` for provenance.

Run command:

```bash
bash pipe/run_chengjisihan.sh
```

## Current Design Configuration

| Setting | Value |
|---|---|
| Display name | Genghis Khan Y-Chromosome ML Tree |
| Target clade | C2a1a3 |
| Subtree output prefix | c2a1a3 |
| Full-tree labels | B, C2, C2a1a3, C2a1a1, C2a1a2, C2a1b, C2b, O, Q, N, J, R, L, D, E, G, H, I |
| Subtree clade labels | C2a1a3, C2a1a3a1, C2a1a3a2, C2a1a3a6 |
| Fine subtree labels | C2a1a3a1a1, C2a1a3a1a2, C2a1a3a6b |
| Published prefix override | BDC |
| Tip label format | population_region_haplogroup |

## Outputs

Annotation outputs:

- `output/chengjisihan/annotation/tip_annotations.tsv`
- `output/chengjisihan/annotation/branch_colors.tsv`
- `output/chengjisihan/annotation/pop_group_colors.tsv`
- `output/chengjisihan/annotation/region_group_colors.tsv`
- `output/chengjisihan/annotation/language_colors.tsv`
- `output/chengjisihan/annotation/design_config.tsv`
- `output/chengjisihan/annotation/visual_design.md`

Figure outputs:

- `output/chengjisihan/figures/full_tree/full_tree_overview.png`
- `output/chengjisihan/figures/full_tree/full_tree_overview.pdf`
- `output/chengjisihan/figures/subtree/c2a1a3_subtree.png`
- `output/chengjisihan/figures/subtree/c2a1a3_subtree.pdf`

## Dataset Summary

The latest annotation table contains 337 tree tips, all matched to metadata. One hundred eleven tips are marked as unpublished/new according to `Status`/`data_published`; BDC-prefixed samples are now present in the main metadata table and are marked as published.

Major haplogroup groups in the full tree:

| Group | Tips |
|---|---:|
| C2a1a3a1 | 62 |
| C2a1a3_base | 52 |
| R | 36 |
| C2a1a3a2 | 25 |
| J | 25 |
| C2a1a3a6 | 19 |
| C2a1a1 | 19 |
| C2b | 15 |
| C2a1a2 | 13 |
| Q | 12 |
| L | 12 |
| N | 10 |
| O | 10 |
| C2a1b | 5 |
| G | 5 |
| D | 5 |
| I | 3 |
| B | 3 |
| E | 3 |
| H | 2 |
| C2a1a3a4 | 1 |

## Target Subtree Summary

The target clade `C2a1a3` contains 159 tips.

Subclade composition:

| Subclade group | Tips |
|---|---:|
| C2a1a3a1 | 62 |
| C2a1a3_base | 52 |
| C2a1a3a2 | 25 |
| C2a1a3a6 | 19 |
| C2a1a3a4 | 1 |

District grouping in the target subtree:

| District group | Tips |
|---|---:|
| Northern China | 103 |
| China | 21 |
| Central Asia | 18 |
| Southern China | 9 |
| South Asia | 5 |
| West Asia | 2 |
| Unknown Region | 1 |

Language composition in the target subtree:

| Language | Tips |
|---|---:|
| Sinitic | 77 |
| Tungusic | 25 |
| Mongolic | 25 |
| Turkic | 20 |
| Indo-European | 7 |
| NA | 5 |

## Figure Design

The full-tree overview uses branch colors to show haplogroup groups and marks macro-haplogroups specified in `VisualDesign.md` with arrow callouts pointing to MRCA nodes. Legends are embedded in blank figure space.

The target subtree separates three visual encodings:

- Branch color: clade identity.
- Tip label and tip point color: district/geographic group.
- Bold tip labels: focal downstream clades listed in `bold_tip_clades`; currently `C2a1a3a6` and downstream tips.
- Annotation strips: `Status`, `District`, and `Clade`, separated by blank spacer columns. The `Status` strip is closest to the tip labels.

The full tree uses the nearest annotation strip for `Status`, followed by a spacer and `Language`; strip legends are embedded in plot whitespace rather than exported as separate legend files. Clade labels in the subtree are placed in nearby blank space and connected to MRCA nodes using thin black arrows. Tip label bolding is controlled by haplogroup membership, not publication status. Fonts use the Arial family; PDF output embeds Arial TrueType fonts so text remains editable.

## Current Scientific Interpretation

C2a1a3 remains the focal candidate lineage in this configured analysis, representing 159 of 337 tips. Within it, C2a1a3a1 is the largest downstream group, followed by the C2a1a3 base group, C2a1a3a2, and C2a1a3a6. The subtree visualization is designed to compare clade structure with geographic distribution: branch colors preserve phylogenetic grouping, while ID colors reveal district-level geography.

BDC samples are no longer treated as private/new by default in this project. They are handled as published samples because `published_prefixes: BDC` is defined in `VisualDesign.md`.

## Reuse Notes

For a new project, copy `conf/VisualDesign.template.md` to a project-specific Markdown file, update `inputs.design_file` in `conf/Config.yaml`, then change:

- `target_clade`
- `output_prefix`
- `full_tree_clade_labels`
- `subtree_clade_labels`
- `subtree_fine_clade_labels`
- `bold_tip_clades`
- `haplogroup_group_rules`
- `published_prefixes`

The pipeline will record the parsed design in `design_config.tsv` and copy the source Markdown to the output annotation directory.
