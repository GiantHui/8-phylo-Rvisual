# Visual Design Specification

## Scientific Questions

1. Write the primary biological or historical question here.
2. Write secondary questions that should guide annotation and highlighting.

## Design Principles

- Integrate sample annotations before visualization.
- Use branch color for clade identity.
- Use tip-label color for geographic, population, phenotype, or other sample-level grouping.
- Keep symbols and labels small enough to avoid overlap.
- Place clade labels in blank space and connect displaced labels to MRCA nodes using thin arrows.
- Use Arial family fonts and keep vector PDF text editable.

## Key Visual Requirements

- Describe the target subtree.
- Describe which macro-clades should be labelled in the full tree.
- Describe which sample groups should be highlighted.
- Describe any published/private/new-sample rules.
- Describe the intended label format.

## Machine-Readable Settings

The pipeline reads only the following fenced block. Keep one `key: value` pair per line.

```visual_design_config
display_name: Example ML Tree
target_clade: C2a1a3
output_prefix: c2a1a3
full_tree_clade_labels: B,C2,C2a1a3,O,Q,N,J,R,L,D,E,G,H,I
subtree_clade_labels: C2a1a3,C2a1a3a1,C2a1a3a2
subtree_fine_clade_labels: C2a1a3a1a1,C2a1a3a1a2
haplogroup_group_rules: C2a1a3a1=C2a1a3a1;C2a1a3a2=C2a1a3a2;C2a1a3=C2a1a3_base
published_prefixes: BDC
tip_label_format: population_region_haplogroup
```
