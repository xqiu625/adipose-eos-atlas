# Revision figures (Reviewer #1 response)

Figure- and report-generators for the PNAS Nexus revision response. They read the
analysis outputs produced by `scripts/45`–`58b` (run on the HPCC) and emit the
Nature-style figures and the self-contained HTML response.

| Generator | Figure / report |
|-----------|-----------------|
| `make_benchmark_figures.py` | Integration-method benchmark (scib-metrics) |
| `make_robustness_figures.py` | Identity (Fig 1) + cluster robustness (Fig 2) |
| `make_integration_figure.py` | ComBat reproduces conclusions (Fig 3) |
| `make_broadref_figure.py` | Broad-reference co-clustering (Fig 4) |
| `make_shah_figure.py` | Shah 2023 bulk concordance (Fig 5) |
| `generate_integration_benchmark_report.py` | Integration-benchmark HTML report |
| `generate_robustness_report.py` | Cluster-robustness analysis memo (HTML) |
| `generate_response_letter_html.py` | Point-by-point Reviewer #1 response (HTML) |
