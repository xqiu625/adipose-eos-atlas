#!/usr/bin/env python
"""Self-contained HTML point-by-point response letter to Reviewer #1 (figures embedded)."""
import base64
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATE = "June 21, 2026"


def b64(p, mime="image/png"):
    return f"data:{mime};base64," + base64.b64encode((HERE / p).read_bytes()).decode()


F1 = b64("fig_identity.png")
F2 = b64("fig_robustness.png")
F3 = b64("fig_integration_conclusions.png")
F4 = b64("fig_broad_reference.png")
F5 = b64("fig_shah_concordance.png")

HTML = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Response to Reviewer #1</title>
<style>
 :root{{--ink:#1a1a1a;--mut:#666;--rule:#e3e3e3;--accent:#0F4D92;--soft:#f6f8fb;--warn:#8a6d00;--warnbg:#fff8e6;}}
 *{{box-sizing:border-box}} body{{font-family:Georgia,'Times New Roman',serif;color:var(--ink);line-height:1.65;margin:0;background:#fafafa}}
 .wrap{{max-width:880px;margin:0 auto;padding:44px 30px 80px;background:#fff}}
 header{{border-bottom:3px solid var(--accent);padding-bottom:16px;margin-bottom:6px}}
 h1{{font-size:1.5em;margin:0 0 6px}} .meta{{color:var(--mut);font-size:.85em;font-family:Arial,sans-serif}}
 h2{{font-size:1.16em;margin:30px 0 8px;color:var(--accent);font-family:Arial,sans-serif}}
 .quote{{border-left:3px solid #c4c4c4;background:#fbfbfb;padding:9px 16px;color:#555;font-style:italic;margin:12px 0;font-size:.95em}}
 .resp{{margin:10px 0}} .resp b,.resp strong{{color:#111}}
 .lead{{background:var(--soft);border-left:4px solid var(--accent);padding:14px 18px;border-radius:0 6px 6px 0;margin:16px 0;font-size:.97em}}
 .note{{background:var(--warnbg);border-left:4px solid var(--warn);padding:10px 15px;border-radius:0 6px 6px 0;margin:14px 0;font-size:.85em;font-family:Arial,sans-serif}}
 ul{{margin:8px 0}} li{{margin:5px 0}}
 figure{{margin:18px 0;text-align:center}} figure img{{max-width:100%;border:1px solid var(--rule);border-radius:6px}}
 figcaption{{font-size:.8em;color:var(--mut);margin-top:7px;text-align:left;font-family:Arial,sans-serif}} figcaption b{{color:var(--ink)}}
 .ph{{color:#b08900;font-family:Arial,sans-serif;font-size:.92em}} code{{background:#f0f0f0;padding:1px 5px;border-radius:3px;font-size:.85em}}
 footer{{margin-top:40px;padding-top:14px;border-top:1px solid var(--rule);color:var(--mut);font-size:.8em;font-family:Arial,sans-serif}}
</style></head><body><div class="wrap">

<header>
 <h1>Response to Reviewer #1</h1>
 <div class="meta">Ruggiero-Ruff &amp; Qiu et&nbsp;al. &middot; PNAS Nexus revision &middot; {DATE}</div>
</header>

<div class="note"><b>Draft.</b> Computational responses are final and backed by the figures below.
Items marked <span class="ph">[FACS]</span>/<span class="ph">[cytospin]</span> depend on R.E.R.'s
experiments; <span class="ph">Fig S__</span> are placeholders to be numbered in the revised manuscript.</div>

<div class="lead">We thank Reviewer #1 for the rigorous and constructive critique. The two central concerns —
whether the six eosinophil clusters are robust biological states, and whether the SiglecF-enriched
dataset is contaminated by macrophages — prompted substantial new analyses and experiments. In brief:
(i) using the discriminating macrophage markers the reviewer specifies — and, most directly, by
co-clustering our cells with published eosinophils against macrophages/mast/neutrophils — the clusters
are eosinophils, not misclassified macrophages; and (ii) we benchmarked the Luecken/scIB methods beyond
Harmony as recommended, and the principal conclusions are reproduced when the data are integrated with
the higher-bio-conservation method (ComBat). We respond point by point below.</div>

<h2>1. SiglecF enrichment and the risk of macrophage contamination</h2>
<div class="quote">“…their enrichment is based on a surface marker that is shared with macrophages, the likelihood of
major contamination by macrophages being misclassified as eosinophils is high.” … “SiglecF alone is not
sufficient to define eosinophils … must be explicitly excluded with macrophage markers such as MerTK,
CD64, F4/80…”</div>
<div class="resp"><p>We agree that SiglecF alone is insufficient, and we now exclude macrophages explicitly at
three independent levels:</p>
<ul>
<li><b>Protein / sorting <span class="ph">[FACS]</span> (new Fig S__).</b> FACS of the adipose SVF with
explicit macrophage exclusion — eosinophils gated as CD45⁺ SiglecF⁺ CCR3⁺ <b>MerTK⁻ CD64⁻ F4/80⁻</b> — and
quantification of the eosinophil vs. macrophage fraction before and after SiglecF enrichment.</li>
<li><b>Morphology <span class="ph">[cytospin]</span> (new Fig S__).</b> Cytospin/Giemsa of the sorted
SiglecF⁺ cells showing the defining eosinophil morphology (bilobed nucleus, eosinophilic granules).</li>
<li><b>Transcriptome (Fig 1 below).</b> In the sequenced data, the eosinophil clusters are negative for the
discriminating macrophage markers <b>MerTK (<i>Mertk</i>, ≤2.2%)</b> and <b>CD64 (<i>Fcgr1</i>, ≤1.6%)</b>;
only <b>1.1%</b> of eosinophils co-express SiglecF with either; the <b>Scrublet doublet rate is 0%</b> in
every cluster; and the macrophages in the same dataset form a separate, marker-positive compartment — a
misclassified macrophage would fall there, but none do.</li>
</ul>
<p>We note that <b>F4/80 (<i>Adgre1</i>) and CD11b (<i>Itgam</i>) are expressed by murine eosinophils</b> and
are not valid for macrophage exclusion; our exclusion is based on MerTK and CD64. The macrophage transcripts
detectable in two clusters (C0, C2) are highly expressed, secretion-prone genes (<i>C1qa/b</i>, <i>Apoe</i>,
<i>Lyz2</i>) — the expected <b>ambient-RNA</b> signature in a sample that is ~59% macrophages — not the
cell-intrinsic surface markers; residual contamination is modest (~4.6% in C0). C0's identity rests on
<b>MHC-II, which we validate at the protein level by spectral flow (Fig 2E–F)</b>; ambient RNA cannot
produce surface protein.</p></div>

<figure><img src="{F1}" alt="Identity figure">
<figcaption><b>Fig 1.</b> (A) Per-cluster identity/QC metrics: C0 is MHC-II-high, all clusters are
intrinsic-macrophage-low with 0% doublets. (B) Only 1.1% of eosinophils co-express SiglecF with MerTK/CD64.
(C) Marker roles: eosinophil markers (blue), discriminating macrophage markers (red, low in eos), and
non-discriminating markers (grey, expressed by eosinophils too).</figcaption></figure>

<div class="resp"><p><b>Most directly, we co-clustered our cells with published eosinophils against the
alternative lineages.</b> We built a broad reference — our full dataset (eosinophils C0–C5 plus our
macrophages, mast cells, neutrophils, dendritic cells, NK and B cells) co-embedded with external
eosinophils (Hu, Gurtner) — and asked where our cells fall (Fig 4). Our C0–C5 sit squarely among the
reference eosinophils: the fraction of each cell's neighbours that are reference eosinophils is
<b>0.23–0.31 for C0–C5 but only 0.02 for macrophages</b> and ≤0.07 for mast cells, neutrophils, and B
cells (Fig 4A–B). A misclassified macrophage would sit with macrophages (fraction ≈0.02); our
eosinophils sit 10–15× higher, with the reference eosinophils. This is the most direct refutation of
the contamination concern. The same analysis (Fig 4C) shows C0–C5 correlate strongly with Hu adipose
eosinophil subgroups (Spearman r = 0.75–0.81 with Eos3/Eos4/Eos5; C3 and C5 best match the homeostatic
Eos3) — i.e., they are bona fide adipose eosinophil states, not artefacts.</p></div>

<figure><img src="{F4}" alt="Broad-reference figure">
<figcaption><b>Fig 4. Broad-reference co-clustering.</b> (A) Joint embedding of our cells with published
eosinophils (Hu, Gurtner): our eosinophils (red) overlap the reference eosinophils (Hu orange, Gurtner
green), while our macrophages (blue) and other lineages (grey) occupy a separate region. (B) Fraction
of each cell's neighbours that are reference eosinophils, by lineage — C0–C5 high, macrophages ≈0.02.
(C) Pseudobulk Spearman correlation of C0–C5 with Hu Eos1–8 (best match boxed): strong correlation with
the adipose subgroups Eos3/Eos4/Eos5.</figcaption></figure>

<h2>2. Fixation / processing and selective cell recovery</h2>
<div class="quote">“If the intended eosinophil-enriched population is substantially altered during fixation and
processing, then the final single-cell dataset may reflect selective survival/recovery of certain cells…”</div>
<div class="resp"><p>We acknowledge this limitation and address it two ways. <b>(i) Orthogonal, single-cell-
independent check.</b> We compared the pseudobulk of our scRNA eosinophils against published bulk RNA-seq of
fresh-sorted eosinophils (Shah, Knights et&nbsp;al. 2023; GEO GSE217551). Our FLEX/fixed eosinophil
transcriptome correlates strongly with the independent bulk (<b>Spearman r = 0.88</b> over all genes; Fig 5A),
and the concordance is <b>eosinophil-specific</b> — on the variable genes (housekeeping removed) our
eosinophils correlate <b>r = 0.66</b> with bulk eosinophils versus only <b>0.24</b> for our macrophages and
0.44 for neutrophils (Fig 5B). This indicates fixation/FLEX processing did not grossly distort the eosinophil
transcriptome. <b>(ii)</b> <span class="ph">[FACS]</span> we compare composition in fresh vs. fixed/processed
material to bound any selective loss. We have also revised the Discussion to state that the recovered
proportions reflect post-fixation composition; our conclusions concern the transcriptional states present and
their reproducible shifts across biological replicates, not absolute in vivo frequencies. <i>(The bulk
comparison speaks to the average transcriptome; the relative proportions of the C0–C5 states are addressed by
the staged FACS.)</i></p></div>

<figure><img src="{F5}" alt="Shah concordance figure">
<figcaption><b>Fig 5. Our eosinophil transcriptome is representative.</b> (A) Pseudobulk of our scRNA
eosinophils vs. independent bulk-sorted gonadal-adipose eosinophils (Shah 2023, GSE217551); Spearman r = 0.88
over all common genes. (B) The concordance is eosinophil-specific: on variable genes our eosinophils correlate
with bulk eosinophils (0.66) far more than our macrophages (0.24) or neutrophils (0.44) do.</figcaption></figure>

<h2>3. Integration and whether the six clusters are robust states</h2>
<div class="quote">“…the argument boils down to ‘the integrated analysis produces different results’… These
differences are expected and are indeed why integration was recommended… Luecken et al. benchmark
recommended software (besides Harmony) that score higher for bio-conservation.”</div>
<div class="resp"><p>We agree the previous justification (ARI/Jaccard differences) was insufficient and have
replaced it with the following, which directly addresses the request:</p>
<ul>
<li><b>Benchmark of methods beyond Harmony (new Fig S__).</b> We benchmarked the Luecken/scIB panel —
ComBat, scVI, scANVI, BBKNN — against the un-integrated baseline on the batch-correction vs. bio-conservation
trade-off. <b>ComBat had the highest bio-conservation among the unsupervised methods.</b></li>
<li><b>The clusters are robust under integration (Fig 2).</b> Across Harmony/ComBat/scVI, per-cluster purity
is <b>0.73–0.97</b> for most clusters (C0 softest, ~0.5); C2 holds at 0.73–0.90.</li>
<li><b>The biological conclusions are reproduced by integration (Fig 3).</b> At matched granularity
(resolution 1.0 → six clusters), <b>ComBat re-separates clusters corresponding to C0–C5 and reproduces the
principal sex-dependent conclusions</b>: a male-enriched MHC-II (C0) cluster (+0.09; KO-female +0.04) and a
female-enriched migratory (C1) cluster (−0.15). Harmony, scoring lower on bio-conservation, loses C0's male
enrichment. Magnitudes are <b>attenuated</b> relative to the un-integrated analysis, consistent with the
multiplexed 10x Flex design partially confounding biological condition with GEM batch.</li>
</ul>
<p>We further document that there is <b>no residual technical batch effect</b> — within-condition replicates
processed in different GEM reactions mix in the embedding (revised Fig S7) — supporting that the un-integrated
structure reflects biology. We therefore retain the un-integrated clustering as the primary representation
(now properly justified) and present the ComBat integration as a validation that the conclusions do not depend
on that choice.</p></div>

<figure><img src="{F2}" alt="Robustness figure">
<figcaption><b>Fig 2.</b> (A) C2 is not an IEG/stress artifact (the IEG-high cluster is the infection cluster
C5). (B) Per-cluster purity under three integration methods — most clusters stay together (0.73–0.97), C0
softest. (C) Within C0, the MHC-II program is internally coordinated, not driven by the ambient macrophage
module.</figcaption></figure>

<figure><img src="{F3}" alt="Integration conclusions figure">
<figcaption><b>Fig 3.</b> (A) ComBat integration (res 1.0, six clusters) recovers the male-enriched MHC-II
(C0) and female-enriched migratory (C1) clusters. (B) ComBat recovers C0's male enrichment (+0.09); Harmony
does not (−0.02). Magnitudes are attenuated vs. the un-integrated analysis (confounded Flex design).</figcaption></figure>

<h2>4. Cluster C2 and per-cluster condition composition</h2>
<div class="quote">“…cluster C2 is clearly dominated (&gt;90%) by WT_Male_HFD … Per-cluster composition should
track cell types or states shared by all conditions … explore why C2 clustered with other cells … Do they
share marker genes … or upregulate a handful of genes (a transcriptional response)?”</div>
<div class="resp"><p>We examined C2 directly, as suggested. C2 is <b>not a male immediate-early/stress
artifact</b> — only 4% of its markers are IEG/stress genes and its IEG module score is mid-pack (the IEG-high
cluster is the infection-associated C5). C2 is <b>not macrophage contamination</b> (MerTK/CD64 low; 0%
doublets) and <b>stays 73–90% together under integration</b>. Its male bias is therefore a property of a
genuine eosinophil state, consistent with the sexually dimorphic eosinophil biology in this model. We have
revised the text to describe C0–C5 as eosinophil <b>states whose abundance is modulated by sex, RELMα, and
infection</b> rather than as condition-independent cell types; the ComBat integration recovers the same
sex-dependent axis (Section 3).</p></div>

<h2>5. Sample size and the weight placed on “unique” clusters</h2>
<div class="quote">“…given the limitations of this study, including small sample size, it seems premature to
place such biological weight on what is said to be unique clusters of eosinophils here.”</div>
<div class="resp"><p>We have tempered the claims accordingly, referring throughout to transcriptional
clusters/states whose abundance shifts with sex, RELMα, and infection, and we now list the modest sample size
explicitly as a limitation. The states are supported by biological-replicate reproducibility (sccomp),
orthogonal protein (spectral flow) and spatial (RNAscope) validation, and cross-dataset correspondence with
independent eosinophil datasets (Hu et al.; Gurtner et al.).</p></div>

<h2>Editorial requirements</h2>
<div class="resp"><ul>
<li><b>Page limit (≤12 pp):</b> the revised Research Report conforms.</li>
<li><b>Data availability:</b> sequencing data are in GEO (<b>GSE311941</b>); code is on GitHub
(github.com/xqiu625/adipose-eos-atlas); the statement is in the Methods.</li>
<li><b>Image integrity / alt-text:</b> no image features were enhanced, moved, or removed; adjustments were
applied to whole images; alt-text has been added under each figure legend.</li>
</ul></div>

<footer>Response to Reviewer #1 &middot; {DATE} &middot; figures embedded; editable SVG/PDF available alongside.</footer>
</div></body></html>"""

out = HERE / "response_to_reviewer1.html"
out.write_text(HTML)
print("wrote", out, f"({len(HTML)/1024:.0f} KB)")
