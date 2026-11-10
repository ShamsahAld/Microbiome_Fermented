# Decisions

- AGP cohort: stool samples only (`sample_type == "Stool"`).
- Consumers = fermented_plant_frequency "Daily" or "Regularly (3-5 times/week)"; non-consumers = "Never". Middle categories dropped for a clear contrast.
- Samples with < 1,000 reads dropped.
- One sample per person (`host_subject_id`), keeping the one with the most reads, to avoid non-independent samples.
- Result: 2,368 consumers, 6,642 non-consumers.
- AGP only records fermented *plant* intake; dairy ferments (kefir, laban) are not captured here.
