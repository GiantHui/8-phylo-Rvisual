# Visual Design Specification

## Scientific Questions

1. What was the impact of Genghis Khan-associated patrilineal elite lineages on Eurasian populations?
2. How did C2a1a3 and its downstream Y-chromosome subclades shape paternal population history?
3. Is C2a1a3 a plausible representative lineage of Genghis Khan descendants, and what phylogenetic evidence supports this?

## Design Principles

- Integrate all annotations before visualization; do not separate the total tree from its legends.
- Use branch color for clade identity and tip-label color for district/geographic grouping.
- Keep tip symbols small enough to avoid overlap.
- Put clade labels in nearby blank space; if displaced from the MRCA, connect them using thin black arrows.
- Use Arial family fonts. Vector PDF text should remain editable and should not be converted into Type3 glyphs or paths.
- Keep the figure within an A4 portrait-oriented footprint.

## Key Visual Requirements

- The target subtree is C2a1a3 and its downstream lineages.
- The full tree should also label non-C macro-haplogroups such as O, Q, N, J, R, L, D, E, G, H, and I.
- BDC samples should be treated as published samples in this project.
- Tip labels should follow `population_region_fullHaplogroup`.
- In the subtree, ID color should represent district/geographic group, not clade identity.
- In the subtree, annotation strips should include separate `District` and `Clade` strips with clear spacing.

## Machine-Readable Settings

The pipeline reads only the following fenced block. Keep one `key: value` pair per line.

```visual_design_config
display_name: Genghis Khan Y-Chromosome ML Tree
target_clade: C2a1a3
output_prefix: c2a1a3
full_tree_clade_labels: B,C2,C2a1a3,C2a1a1,C2a1a2,C2a1b,C2b,O,Q,N,J,R,L,D,E,G,H,I
subtree_clade_labels: C2a1a3,C2a1a3a1,C2a1a3a2,C2a1a3a6
subtree_fine_clade_labels: C2a1a3a1a1,C2a1a3a1a2,C2a1a3a6b
haplogroup_group_rules: C2a1a3a1=C2a1a3a1;C2a1a3a2=C2a1a3a2;C2a1a3a4=C2a1a3a4;C2a1a3a6=C2a1a3a6;C2a1a3=C2a1a3_base;C2a1a1=C2a1a1;C2a1a2=C2a1a2;C2a1b=C2a1b;C2b=C2b
published_prefixes: BDC
tip_label_format: population_region_haplogroup
```
