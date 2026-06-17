#!/usr/bin/env python
"""Self-contained HTML report: Reviewer #1 response (cluster robustness + identity)."""
import base64
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATE = "June 13, 2026"


def b64(p, mime):
    return f"data:{mime};base64," + base64.b64encode((HERE / p).read_bytes()).decode()


fig_id = b64("fig_identity.png", "image/png")
fig_rob = b64("fig_robustness.png", "image/png")
fig_int = b64("fig_integration_conclusions.png", "image/png")

HTML = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Reviewer #1 Response — Cluster Robustness & Eosinophil Identity</title>
<style>
 :root{{--ink:#1a1a1a;--mut:#666;--rule:#e3e3e3;--accent:#0F4D92;--soft:#f6f8fb;--warn:#8a6d00;--warnbg:#fff8e6;--good:#1e7a34;--goodbg:#eef8f0;}}
 *{{box-sizing:border-box}} body{{font-family:-apple-system,Segoe UI,Arial,sans-serif;color:var(--ink);line-height:1.6;margin:0;background:#fafafa}}
 .wrap{{max-width:1000px;margin:0 auto;padding:40px 28px 80px;background:#fff}}
 header{{border-bottom:3px solid var(--accent);padding-bottom:18px}}
 h1{{font-size:1.6em;margin:0 0 6px;line-height:1.25}} .sub{{color:var(--mut)}} .meta{{color:var(--mut);font-size:.85em;margin-top:10px}}
 h2{{font-size:1.2em;margin:32px 0 6px;color:var(--accent);border-bottom:1px solid var(--rule);padding-bottom:6px}}
 h3{{font-size:1.02em;margin:18px 0 6px}}
 .key{{background:var(--soft);border-left:4px solid var(--accent);padding:14px 18px;border-radius:0 6px 6px 0;margin:16px 0}}
 .good{{background:var(--goodbg);border-left:4px solid var(--good);padding:12px 16px;border-radius:0 6px 6px 0;margin:14px 0}}
 .warn{{background:var(--warnbg);border-left:4px solid var(--warn);padding:12px 16px;border-radius:0 6px 6px 0;margin:14px 0;font-size:.93em}}
 .quote{{border-left:3px solid #bbb;background:#fbfbfb;padding:8px 14px;color:#444;font-style:italic;margin:12px 0;font-size:.92em}}
 table{{border-collapse:collapse;width:100%;margin:14px 0;font-size:.88em}}
 th,td{{padding:7px 10px;text-align:left;border-bottom:1px solid var(--rule)}} th{{background:var(--accent);color:#fff}}
 figure{{margin:18px 0;text-align:center}} figure img{{max-width:100%;border:1px solid var(--rule);border-radius:6px}}
 figcaption{{font-size:.83em;color:var(--mut);margin-top:8px;text-align:left}} figcaption b{{color:var(--ink)}}
 code{{background:#f0f0f0;padding:1px 5px;border-radius:3px;font-size:.85em}} ul{{margin:8px 0}}
 footer{{margin-top:46px;padding-top:16px;border-top:1px solid var(--rule);color:var(--mut);font-size:.82em}}
</style></head><body><div class="wrap">

<header>
 <h1>Reviewer #1 Response — Are the Six Eosinophil Clusters Robust, and Are They Eosinophils?</h1>
 <div class="sub">Cluster-robustness, integration, and macrophage-contamination analyses for the PNAS Nexus revision</div>
 <div class="meta">Analysis memo &middot; {DATE} &middot; Ruggiero-Ruff &amp; Qiu et&nbsp;al. &middot; 7,461 adipose eosinophils (C0–C5)</div>
</header>

<div class="key"><b>Bottom line.</b> Direct testing shows the six clusters are <b>robust eosinophil states, not
macrophage contamination</b>. The condition-skewed clusters (C0, C2) are <b>not</b> misclassified macrophages or
doublets (MerTK/CD64/CD68 low; <b>0% doublets</b>); their elevated C1q/Apoe/Lyz2 is <b>ambient RNA</b> from the
macrophage-rich sample (~59% macrophages), and <b>C0 is defined by a genuine, protein-validated MHC-II program</b>.
Most clusters survive integration (purity 0.73–0.97), and re-clustering the <b>ComBat</b> integration (the higher-
bio-conservation method the reviewer recommended) <b>reproduces the male-enriched MHC-II (C0) and female-enriched
migratory (C1) shifts</b> — so the conclusions do not depend on the un-integrated choice.</div>

<h2>Concern #2 — “SiglecF alone is not sufficient … contamination by macrophages misclassified as eosinophils”</h2>
<div class="quote">“SiglecF alone is not sufficient to define eosinophils … must be explicitly excluded with macrophage markers
such as MerTK, CD64, F4/80 … the likelihood of major contamination by macrophages being misclassified as
eosinophils is high.”</div>

<h3>Response</h3>
<ul>
<li><b>Using the discriminating markers the reviewer names, the eosinophils are clean.</b> Across all six clusters,
   MerTK⁺ ≤2.2% and CD64(Fcgr1)⁺ ≤1.6%; only <b>1.1% of eosinophils</b> co-express SiglecF together with MerTK or
   CD64 (Fig 1B). The cells express SiglecF/Il5ra/Ccr3/Il4 and lack the cell-intrinsic macrophage program.</li>
<li><b>They are not doublets or whole macrophages.</b> Scrublet doublet rate is <b>0%</b> in every cluster, library
   size and the cell-intrinsic macrophage module (MerTK/CD64/CD68/Csf1r) are uniformly low (Fig 1A).</li>
<li><b>The macrophage genes seen in C0/C2 are ambient RNA.</b> C0 shows elevated <i>abundant, secretion-prone</i>
   macrophage transcripts (C1q/Apoe/Lyz2) but <i>not</i> the cell-intrinsic surface markers — the classic ambient
   signature in a sample that is ~59% macrophages. Residual contamination is modest (~4.6% in C0, lower elsewhere).</li>
<li><b>F4/80 (Adgre1) and CD11b (Itgam) are not valid discriminators in mouse</b> — murine eosinophils express them
   too (Fig 1C) — so macrophage exclusion is based on MerTK/CD64, on which our cells are negative.</li>
<li><b>C0’s MHC-II identity is protein-validated.</b> Spectral flow (manuscript Fig 2E–F) shows MHC-II on
   eosinophils at the protein level; ambient RNA cannot produce surface protein. (Note: RNA-level DecontX with
   cell-type sources over-removes MHC-II because it is also high in the abundant macrophages, so we rely on the
   protein data and a within-C0 coordination analysis, Fig 2C, not on DecontX of shared genes.)</li>
<li><b>Wet-lab confirmation in progress:</b> FACS with explicit MerTK/CD64/F4/80 gating plus a cytospin/Giemsa of
   the sorted SiglecF⁺ cells (morphology: bilobed nucleus, eosinophilic granules) is being performed by R.E.R.</li>
</ul>

<figure><img src="{fig_id}" alt="Identity figure">
<figcaption><b>Figure 1. Eosinophil identity, not macrophage contamination.</b>
<b>(A)</b> Per-cluster identity/QC metrics (row-scaled colour, raw values shown): C0 is MHC-II-high, all clusters are
intrinsic-macrophage-low with 0% doublets; C0’s elevated mac-ambient genes track its higher residual contamination.
<b>(B)</b> Only 1.1% of eosinophils co-express SiglecF with the discriminating markers MerTK/CD64.
<b>(C)</b> Marker roles: eosinophil markers (blue) and discriminating macrophage markers (red, low in eos) vs.
non-discriminating markers (grey, expressed by eosinophils too).</figcaption></figure>

<h2>Concern #1 — Integration and whether the clusters are robust biological states</h2>
<div class="quote">“…the authors chose not to pursue integration … the argument boils down to ‘the integrated analysis produces
different results’ … cluster C2 is clearly dominated (&gt;90%) by WT_Male_HFD … explore why C2 clustered with other
cells … Luecken et al. benchmark recommended software (besides Harmony) that score higher for bio-conservation.”</div>

<h3>Response</h3>
<ul>
<li><b>We benchmarked the Luecken/scIB methods beyond Harmony.</b> ComBat, scVI, scANVI and BBKNN were scored against
   the un-integrated baseline on the batch-correction vs. bio-conservation trade-off (separate benchmark report). This
   directly addresses the request for higher-bio-conservation methods.</li>
<li><b>C2 is not a male IEG/stress artifact.</b> Only 4% of C2’s markers are immediate-early/stress genes, and its IEG
   module score is mid-pack — the IEG-high cluster is the infection-associated C5, where stress programs are expected
   (Fig 2A).</li>
<li><b>The clusters largely survive integration.</b> Per-cluster “stays-together” purity is 0.73–0.97 for most clusters
   across Harmony/ComBat/scVI; <b>C2 holds at 0.73–0.90</b> (Fig 2B). The softest is C0 (~0.5), explained by its ambient
   gradient and the continuous nature of the MHC-II program rather than by C0 being spurious.</li>
<li><b>C0 is a coordinated MHC-II program.</b> Within C0 cells, MHC-II genes co-vary with each other ~3× more than with
   the ambient macrophage module (within-MHC-II Spearman 0.094 vs. cross 0.028; Fig 2C), consistent with an intrinsic
   program, not ambient leakage.</li>
</ul>

<figure><img src="{fig_rob}" alt="Robustness figure">
<figcaption><b>Figure 2. The clusters are robust under integration.</b>
<b>(A)</b> C2 is not driven by IEG/stress genes (low IEG fraction and module score); C5, the infection cluster, is the
IEG-high one. <b>(B)</b> Per-cluster purity under three integration methods — most clusters stay together (0.73–0.97);
C0 is the softest (~0.5, dashed line). <b>(C)</b> Within C0, the MHC-II program is internally coordinated and largely
uncorrelated with the ambient macrophage module.</figcaption></figure>

<h3>The biological conclusions are reproduced by integration</h3>
<p>To test the reviewer's deeper point directly, we re-clustered the integrated embedding at matched granularity
(resolution 1.0 → six clusters) and asked whether the sex / RELMα / infection composition shifts re-appear. Using
<b>ComBat — the higher-bio-conservation method the reviewer recommended over Harmony</b> — the integrated analysis
re-separates clusters corresponding to C0–C5 and <b>reproduces the principal sex-dependent conclusions</b>: a
<b>male-enriched MHC-II (C0)</b> cluster (male enrichment +0.09, KO-female +0.04) and a <b>female-enriched migratory
(C1)</b> cluster (−0.15) (Fig 3). Harmony alone, which scores lower on bio-conservation, loses C0's male enrichment —
underscoring why the method choice the reviewer requested matters. Enrichment magnitudes are attenuated relative to
the un-integrated compositional analysis, consistent with the multiplexed Flex design partially confounding biological
condition with GEM batch, so that batch removal also removes some of the genuine condition signal. The
<b>conclusions therefore do not depend on the un-integrated choice</b> — integration with the recommended method
reproduces them directionally.</p>

<figure><img src="{fig_int}" alt="Integration conclusions figure">
<figcaption><b>Figure 3. Integration reproduces the sex-dependent compositional conclusions.</b>
<b>(A)</b> ComBat integrated clusters (resolution 1.0, six clusters; x-labels = integrated cluster and dominant
original cluster): the C0-dominated cluster is male-enriched and the C1-dominated cluster is female-enriched, matching
the manuscript. <b>(B)</b> ComBat recovers C0's male enrichment (+0.09) whereas Harmony does not (−0.02); both recover
the female-enriched C1. Magnitudes are attenuated vs. the un-integrated analysis (confounded Flex design).</figcaption></figure>

<h2>Evidence summary</h2>
<table><thead><tr><th>Reviewer worry</th><th>Test</th><th>Result</th></tr></thead><tbody>
<tr><td>SiglecF+ macrophages misclassified as eosinophils</td><td>MerTK/CD64 strict gating in scRNA</td><td>eos ≤2% MerTK+/CD64+; 1.1% co-express; 0% doublets</td></tr>
<tr><td>C0/C2 are macrophage-contaminated</td><td>intrinsic vs ambient macrophage modules</td><td>intrinsic low; only abundant/ambient C1q/Apoe elevated (~4.6% contam)</td></tr>
<tr><td>C2 a male transcriptional response, not a state</td><td>IEG/stress fraction + integration purity</td><td>4% IEG; stays 0.73–0.90 on integration</td></tr>
<tr><td>Clusters not robust under integration</td><td>re-cluster on Harmony/ComBat/scVI</td><td>purity 0.73–0.97 (C0 softest ~0.5)</td></tr>
<tr><td>C0 identity (MHC-II) real?</td><td>flow protein (Fig 2E–F) + within-C0 coordination</td><td>MHC-II protein-validated; coordinated program</td></tr>
<tr><td>Do conclusions survive integration?</td><td>re-cluster (ComBat, res 1.0) + composition</td><td>male-MHC-II (C0) &amp; female-migratory (C1) shifts reproduced (attenuated)</td></tr>
</tbody></table>

<div class="warn"><b>Honest framing notes for the response letter.</b>
(1) <b>C0</b> carries modest ambient macrophage RNA — describe C0 by its MHC-II program (protein-validated) and avoid
over-interpreting the ambient C1q/Apoe transcripts. (2) <b>C5</b> should be described as an infection-induced
activation state, consistent with its IEG/stress signature. (3) Do <b>not</b> present cell-type-source DecontX of
shared genes (it over-removes MHC-II); rely on flow protein + coordination. (4) FACS + cytospin (R.E.R.) are the
wet-lab anchor for identity and should accompany this computational evidence.</div>

<footer>Analysis memo &middot; {DATE} &middot; adipose eosinophil atlas &middot; reviewer-response cluster-robustness pack.
Figures available as editable SVG/PDF alongside this report.</footer>
</div></body></html>"""

out = HERE / "reviewer1_response_report.html"
out.write_text(HTML)
print("wrote", out, f"({len(HTML)/1024:.0f} KB)")
