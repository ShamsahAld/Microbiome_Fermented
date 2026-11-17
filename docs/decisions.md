# Decisions

- AGP cohort: stool samples only (`sample_type == "Stool"`).
- Consumers = fermented_plant_frequency "Daily" or "Regularly (3-5 times/week)"; non-consumers = "Never". Middle categories dropped for a clear contrast.
- Samples with < 1,000 reads dropped.
- One sample per person (`host_subject_id`), keeping the one with the most reads, to avoid non-independent samples.
- Result: 2,368 consumers, 6,642 non-consumers.
- AGP only records fermented *plant* intake; dairy ferments (kefir, laban) are not captured here.
- ASVs kept if present in >= 1% of cohort samples: 1,940 of 301,708 ASVs, retaining 96.4% of reads.
- Taxonomy: QIIME 2 2024.10 (Singularity image quay.io/qiime2/amplicon:2024.10), classify-sklearn with the Greengenes2 2024.09 V4 (515F/806R) classifier (SHA256 6cafb03b…). 1,591 of 1,940 ASVs (82%) assigned to genus.
- ASVs without a genus labelled by their deepest known rank (e.g. Lachnospiraceae_unclassified) instead of being dropped.
- Genera kept if present in >= 10% of samples: 239 of 505 genera, retaining 96.1% of reads.
