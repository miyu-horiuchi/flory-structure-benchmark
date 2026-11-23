# Exploratory Lake grouped-scaffold sensitivity — 2026-11-23 UTC

Added after seeing the sealed primary results and before training any grouped Lake folds. This is a separate post-hoc analysis, not part of the sealed prespec. The 24 Lake compounds use the same no-TPX2 pKd endpoint as the sealed coupled run. It shares a global fit with cooperativity and is circular on either split.

Use four folds of six compounds, the smallest permitted k that achieves exact equality here. Use the same intact scaffold-group balancing rule as the primary exploratory plan: sort by decreasing size and then canonical typed scaffold SMILES, assign each whole group to the fold with fewest compounds, and break ties by fold ID. Five folds could not be exactly equal with 24 compounds. No labels or feature values determine the assignment.

Use the sealed Lake recipe apart from splitting: the same A/B arms and LoRA/frozen variants, five seeds, fold mean, cooperativity alone, ten row shuffles, ten Gaussian shams, fold-local standardisation, metrics and paired compound bootstrap. Keep results separate from the sealed headline.
