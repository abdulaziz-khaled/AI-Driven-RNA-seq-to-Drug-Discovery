# Decoding Hashimoto's Thyroiditis — AI-Driven RNA-seq to Drug Discovery

> **End-to-end transcriptomic analysis of Hashimoto's Thyroiditis — from raw expression data to drug target identification and structure-based virtual screening — conducted entirely by Biomni, the AI research agent on the Phylo platform.**

---

## Table of Contents

1. [About Phylo & Biomni](#about-phylo--biomni)
2. [Project Overview](#project-overview)
3. [Dataset](#dataset)
4. [Analysis Pipeline](#analysis-pipeline)
5. [Key Results](#key-results)
   - [Differential Expression](#differential-expression)
   - [Pathway Enrichment](#pathway-enrichment)
   - [Drug Target Prioritization](#drug-target-prioritization)
   - [Structure-Based Virtual Screening](#structure-based-virtual-screening)
   - [ADMET Analysis](#admet-analysis)
   - [Binding Energy Analysis](#binding-energy-analysis)
6. [Figures](#figures)
7. [Output Files](#output-files)
8. [Nextflow Pipeline](#nextflow-pipeline)
9. [Software & Versions](#software--versions)
10. [Limitations](#limitations)
11. [Citation](#citation)

---

## About Phylo & Biomni

**[Phylo](https://phylo.bio)** is an AI-native biomedical research platform designed to accelerate scientific discovery by deploying autonomous AI agents across the full drug discovery and genomics pipeline. Rather than providing static tools, Phylo orchestrates intelligent agents that reason, plan, execute multi-step analyses, and critically interpret results — functioning as a collaborative scientific partner rather than a passive software suite.

**Biomni** is Phylo's flagship AI research agent. Built on state-of-the-art large language models and equipped with a comprehensive suite of bioinformatics tools, databases, and computational chemistry capabilities, Biomni can:

- Design and execute complete RNA-seq pipelines (QC → alignment → DE → enrichment)
- Perform pathway enrichment analysis (GSEA, ORA) across MSigDB, KEGG, and GO databases
- Prioritize drug targets using multi-dimensional scoring against Open Targets, DepMap, and clinical databases
- Conduct structure-based virtual screening: pocket identification, molecular docking (AutoDock Vina), and contact analysis
- Predict ADMET properties using deep learning models (DeepPurpose MPNN)
- Compute binding energy decomposition (MM interaction energy, per-residue analysis)
- Generate publication-quality figures and reproducible Nextflow pipelines

All analyses in this repository were carried out autonomously by Biomni on the Phylo platform, with the researcher providing high-level direction and scientific oversight. The platform URL for this project session is: `https://biomni.phylo.bio`

---

## Project Overview

This project applies a full AI-driven drug discovery pipeline to **Hashimoto's Thyroiditis (HT)**, the most common autoimmune thyroid disease, using publicly available bulk RNA-seq data from thyroid tissue biopsies.

The analysis spans five major stages:

```
Bulk RNA-seq Data (HRA001684)
        │
        ▼
1. Differential Expression Analysis
   └─ 2,242 DEGs identified (padj < 0.05, |log2FC| ≥ 1)
        │
        ▼
2. Pathway Enrichment (GSEA + ORA)
   └─ 1,006 significant gene sets; immune activation dominant
        │
        ▼
3. Drug Target Prioritization
   └─ IL12B selected as primary target (ustekinumab repositioning)
        │
        ▼
4. Structure-Based Virtual Screening
   └─ 10 compounds docked; ZINC52448673 & ZINC55760081 as leads
        │
        ▼
5. ADMET + Binding Energy Analysis
   └─ ZINC55760081 recommended as development candidate
```

**Disease:** Hashimoto's Thyroiditis (autoimmune hypothyroidism)  
**Target:** IL-12B (p40 subunit of IL-12/IL-23)  
**Proposed drug:** Ustekinumab repositioning (anti-IL12/23 p40 mAb)  
**Analysis date:** April 2026

---

## Dataset

| Parameter | Value |
|-----------|-------|
| **Accession** | NGDC GSA-Human [HRA001684](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001684) |
| **Publication** | Zhang et al. (2022) *Nature Communications* 13:810 |
| **DOI** | [10.1038/s41467-022-28120-2](https://doi.org/10.1038/s41467-022-28120-2) |
| **Total samples** | 66 thyroid tissue biopsies |
| **HT (cases)** | 16 |
| **Non-HT (controls)** | 50 |
| **Library batches** | 3 (n=16, 18, 32) |
| **Genes measured** | 28,029 |
| **Genes with valid stats** | 18,983 |
| **Tissue** | Thyroid |
| **Species** | *Homo sapiens* (GRCh38) |
| **Expression data** | DESeq2 size-factor normalized counts (authors' MOESM7) |
| **DE method** | Authors' LME + DESeq2 (batch-corrected: `~ batch + condition`) |

> **Note:** Raw FASTQ files (~54 GB/sample; ~1.6 TB total) were not re-downloaded. Analysis was performed on the authors' processed supplementary data. Re-running DESeq2 on size-factor normalized values would be methodologically incorrect; the authors' pre-computed statistics were used directly.

---

## Analysis Pipeline

### Stage 1 — Differential Expression

- **Significance threshold:** padj < 0.05 AND |log2FC| ≥ 1 (Benjamini-Hochberg correction)
- **Batch correction:** Linear mixed-effects model (3 library preparation batches)
- **Visualizations:** PCA, Volcano plot, MA plot, Top-50 DEG heatmap (ComplexHeatmap, Z-scored log2(normalized+1))

### Stage 2 — Functional Enrichment

- **GSEA:** Genes ranked by `sign(log2FC) × −log10(padj)`; run via clusterProfiler v4.14 / fgsea
- **Databases:** MSigDB Hallmark (50 gene sets), KEGG Legacy (186 gene sets), GO:BP (7,583 gene sets)
- **ORA:** Hypergeometric test; upregulated DEGs (n=2,107) and downregulated DEGs (n=135) vs background (n=18,951)
- **FDR threshold:** < 0.05

### Stage 3 — Drug Target Prioritization

- **Candidate pool:** Top 46 pathway-central DEGs (ranked by GSEA leading-edge + ORA pathway co-appearance)
- **Scoring dimensions (5):**

| Dimension | Weight | Rationale |
|-----------|--------|-----------|
| Druggability (SM/AB tractability + clinical stage) | 30% | Feasibility of therapeutic intervention |
| Pathway centrality | 25% | Network importance in HT biology |
| Expression magnitude (\|log2FC\|) | 20% | Biological relevance |
| Disease relevance (Open Targets association) | 15% | Evidence linking target to thyroid autoimmunity |
| Autoimmune precedent (approved drugs) | 10% | Safety profile and translational readiness |

- **Database:** Open Targets Platform v26.03 (GraphQL API)

### Stage 4 — Structure-Based Virtual Screening

- **Target structure:** PDB [3HMX](https://www.rcsb.org/structure/3HMX) (IL-12B/p40 in complex with ustekinumab Fab; Chain A = target)
- **Cross-validation:** PDB [8YI7](https://www.rcsb.org/structure/8YI7) (IL-12B/p40 in complex with IL-12Rβ1)
- **Epitope pocket:** 31-residue unified pocket (Fab contacts ∪ IL-12Rβ1 contacts); 8 core hotspot residues: TRP15, PRO17, ASP18, ALA19, GLU45, LYS58, GLU59, PHE60
- **Docking:** AutoDock Vina v1.2.5; box center (5.99, −4.57, −1.37) Å, 30×30×30 Å; exhaustiveness=32, seed=42
- **Receptor preparation:** Fab chains H/L stripped; Chain A+B only used for docking
- **Compounds screened:** 10 platform-generated ligands (batch validation)
- **Contact analysis:** 4.5 Å cutoff for core hotspot hits; 12 Å shell for binding site definition

### Stage 5 — ADMET & Binding Energy

- **ADMET:** DeepPurpose MPNN model; RDKit physicochemical pre-screening (Lipinski Ro5, Veber rules)
- **SMILES extraction:** OpenBabel PDBQT → SMILES conversion
- **Binding energy:** Single-point MM interaction energy (ε=4r distance-dependent dielectric, soft-core LJ vdW); per-residue energy decomposition; MMFF94 ligand charges; AMBER-like backbone charges

---

## Key Results

### Differential Expression

| Category | Count |
|----------|-------|
| Total genes tested | 18,983 |
| Significant (padj < 0.05) | 11,098 |
| **Significant DEGs (padj < 0.05, \|log2FC\| ≥ 1)** | **2,242** |
| Upregulated in HT | 2,107 (94.0%) |
| Downregulated in HT | 135 (6.0%) |

The striking 15.6× up/down asymmetry is biologically expected: HT is a gain-of-function immune infiltration disease. The upregulated genes reflect lymphocyte invasion; the downregulated genes (MYOC, GLIS3-AS1, GLDC) are thyroid follicular cell markers displaced by infiltrating immune cells.

**Top upregulated genes (HT vs Non-HT):**

| Gene | log2FC | padj | Biological Role |
|------|--------|------|-----------------|
| TERT | 9.83 | 7.8×10⁻⁹ | Telomerase; lymphocyte proliferation |
| SERPINA9 | 9.18 | 2.7×10⁻⁶ | B-cell germinal center marker |
| FCAMR | 8.32 | 2.5×10⁻¹⁰ | IgA/IgM Fc receptor; B cells |
| SPIC | 7.90 | 4.6×10⁻¹⁰ | B-cell transcription factor |
| CR2 (CD21) | 7.59 | 2.4×10⁻⁹ | B-cell co-receptor; complement |
| FCRL4 | 7.35 | 1.4×10⁻¹⁴ | Fc receptor-like; memory B cells |
| TCL1A | 7.27 | 2.7×10⁻⁹ | B-cell survival/proliferation |
| PAX5 | 7.03 | 3.3×10⁻¹⁰ | Master B-cell transcription factor |
| IL12B | 3.07 | 3.4×10⁻⁹ | IL-12/23 p40 subunit; Th1 driver |

**Top downregulated genes:**

| Gene | log2FC | padj | Biological Role |
|------|--------|------|-----------------|
| MYOC | −2.19 | 3.4×10⁻⁴ | Myocilin; extracellular matrix |
| AMER3 | −2.19 | 5.0×10⁻³ | APC membrane recruitment |
| GLDC | −2.03 | 0.030 | Glycine decarboxylase; metabolism |
| GLIS3-AS1 | −1.96 | 7.2×10⁻³ | GLIS3 antisense; thyroid development |
| TNN | −1.68 | 1.5×10⁻⁵ | Tenascin-N; extracellular matrix |

---

### Pathway Enrichment

**GSEA summary (1,006 significant gene sets, FDR < 0.05):**

| Database | Activated | Suppressed |
|----------|-----------|------------|
| GO:BP | 884 | 51 |
| KEGG | 44 | 9 |
| Hallmark | 15 | 3 |

**Top activated Hallmark pathways:**

| Pathway | NES | FDR |
|---------|-----|-----|
| HALLMARK_ALLOGRAFT_REJECTION | 2.16 | 2.7×10⁻⁹ |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 2.04 | 2.7×10⁻⁹ |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | 2.00 | 2.7×10⁻⁹ |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | 1.98 | 2.7×10⁻⁹ |
| HALLMARK_INFLAMMATORY_RESPONSE | 1.93 | 2.7×10⁻⁹ |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 1.78 | 2.7×10⁻⁹ |

**Top suppressed Hallmark pathways:**

| Pathway | NES | FDR |
|---------|-----|-----|
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | −1.59 | 1.1×10⁻⁵ |
| HALLMARK_ADIPOGENESIS | −1.36 | 8.6×10⁻⁴ |
| HALLMARK_MYOGENESIS | −1.19 | 0.029 |

Hierarchical clustering of Hallmark pathways by leading-edge Jaccard similarity reveals 4 biological clusters: (1) immune activation, (2) cell cycle/lymphocyte proliferation, (3) metabolic suppression, (4) myogenesis.

---

### Drug Target Prioritization

**Top 5 ranked candidates:**

| Rank | Gene | Score | Approved Drug | HT Association |
|------|------|-------|---------------|----------------|
| 1 | IFNG | 0.827 | Emapalumab | Hypothyroidism (OT) — **overranked*** |
| **2** | **IL12B** | **0.814** | **Ustekinumab** | **Autoimmune thyroid disease (OT)** |
| 3 | TNF | 0.781 | Etanercept/Adalimumab | None — **contraindicated*** |
| 4 | BTK | 0.754 | Ibrutinib/Rilzabrutinib | None (data gap) |
| 5 | JAK3 | 0.697 | Ritlecitinib/Tofacitinib | None |

> *IFNG (rank 1) is overranked: its disease score is inflated by an Open Targets "hypothyroidism" association (downstream outcome, not the autoimmune process). TNF (rank 3) carries a documented safety contraindication — anti-TNF therapy can paradoxically induce destructive thyroiditis.

**Primary recommendation: IL12B / Ustekinumab repositioning**

IL12B (encoding the p40 subunit shared by IL-12 and IL-23) was selected based on convergent evidence:
- Strong differential expression: log2FC = +3.07, padj = 3.4×10⁻⁹ (~8× higher in HT)
- High pathway centrality: rank 6 of 2,242 DEGs (appearing in 161 GSEA leading edges + 168 ORA pathways)
- Direct Open Targets association: "autoimmune thyroid disease" (score = 0.204)
- Genetic evidence: IL12B SNPs associated with HT clinical presentation and severity
- Elevated serum IL-12p40 in HT patients confirmed by independent ELISA studies
- Approved drug (ustekinumab) with autoimmune precedent: psoriasis, PsA, Crohn's disease, UC
- No safety concern of worsening thyroid autoimmunity (unlike TNF inhibitors)
- Dual mechanism: blocks both IL-12 (Th1) and IL-23 (Th17) — relevant to HT's mixed Th1/Th17 phenotype

**IL12B co-expression within HT samples (n=16):**

| Gene | Within-HT r | p-value | Biological Link |
|------|-------------|---------|-----------------|
| CXCL9 | 0.78 | <0.001 | IFN-γ-induced chemokine, T-cell recruitment |
| CXCL10 | 0.76 | <0.001 | IFN-γ-induced chemokine |
| GZMB | 0.75 | 0.001 | Cytotoxic T-cell effector |
| IFNG | 0.72 | 0.002 | Th1 effector cytokine |
| CD8A | 0.68 | 0.004 | Cytotoxic T-cell marker |
| PRF1 | 0.66 | 0.006 | Perforin, cytotoxic killing |

---

### Structure-Based Virtual Screening

**Epitope pocket definition (3HMX Chain A numbering):**
- Fab contact residues (22): 15,17,18,19,20,21,23,40,41,42,43,45,47,54,55,56,58,59,60,61,62,66
- IL-12Rβ1 contact residues (17): 14,15,16,17,18,19,45,58,59,60,84,86,93,194,195,196,197
- **Unified pocket (31 residues):** union of both sets
- **Core hotspots (8):** TRP15, PRO17, ASP18, ALA19, GLU45, LYS58, GLU59, PHE60

**Batch docking results (10 compounds, clean receptor — Fab chains stripped):**

| Rank | ZINC ID | Vina Score | Core Hits | Rβ1 Hits | Epi Frac | Verdict |
|------|---------|------------|-----------|----------|----------|---------|
| 1 | **ZINC52448673** | −10.41 | **6/8** | 13/17 | 100% | STRONG HIT |
| 2 | **ZINC55760081** | −7.75 | **6/8** | 11/17 | 73% | STRONG HIT |
| 3 | ZINC53075763 | −8.26 | 2/8 | 9/17 | 64% | Moderate |
| 4 | ZINC52672394 | −7.49 | 2/8 | 8/17 | 90% | Moderate |
| 5 | ZINC52004842 | −8.25 | 2/8 | 2/17 | 93% | Moderate |
| 6–10 | Various | −7.49 to −10.18 | 1/8 | 5/17 | 44–78% | Weak |

> **Platform claim verification:** The docking platform claimed "all 8 core hotspot residues engaged." Independent re-docking confirmed 6/8 — ASP18 and GLU45 are NOT contacted at a strict 4.5 Å cutoff. The platform likely uses a looser cutoff (~5.5–6 Å) or reports pocket-defining residues rather than strict contacts.

> **False positive flag:** ZINC49219928 scored −10.18 kcal/mol (raw) but only 1/8 core hits — a deceptive high-scoring pose that misses the functional epitope.

---

### ADMET Analysis

SMILES extracted via OpenBabel from PDBQT files; ADMET predicted using DeepPurpose MPNN model.

| Property | ZINC52448673 | ZINC55760081 | Threshold |
|----------|-------------|-------------|-----------|
| MW (Da) | 555.7 ⚠ | 530.4 ⚠ | ≤500 |
| LogP | 5.53 ⚠ | **2.83 ✓** | ≤5 |
| TPSA (Å²) | 60.9 ✓ | 152.1 ⚠ | ≤140 |
| HIA | 98.5% ✓ | 95.9% ✓ | >80% |
| CYP3A4 inhibition | **85.1% ✗** | **18.1% ✓** | <50% |
| Pgp inhibitor | 85.0% ⚠ | 13.3% ✓ | <50% |
| Clinical toxicity | 31.7% ✓ | 27.9% ✓ | <50% |
| QED | 0.294 | 0.391 | >0.5 ideal |

**SMILES:**
- **ZINC52448673:** `c1c(cc2-c3c(Cc2c1)cccc3)C(=O)N1CC[C@]2(C1)CN(CC2)C(=O)[C@H]1CC(=O)N(C1)c1cccc2c1cccc2`
- **ZINC55760081:** `O=C1N(N=C(/C/1=C/c1oc(cc1)c1ccc2c(c1)C(=O)NC2=O)C(F)(F)F)c1ccc(cc1)S(=O)(=O)N`

**Verdict:** ZINC55760081 is the better drug development candidate despite a weaker docking score — its CYP3A4 profile (18% vs 85%) eliminates the major DDI risk that disqualifies ZINC52448673.

---

### Binding Energy Analysis

Single-point MM interaction energy (distance-dependent dielectric ε=4r, soft-core LJ vdW, MMFF94 ligand charges, AMBER-like backbone charges):

| Metric | ZINC52448673 | ZINC55760081 |
|--------|-------------|-------------|
| Vina docking score | −10.41 kcal/mol | −7.75 kcal/mol |
| ΔE_elec (ε=4r) | −29.20 kcal/mol | −17.46 kcal/mol |
| ΔE_vdW (soft-core) | −28.82 kcal/mol | −28.10 kcal/mol |
| **ΔE_MM total** | **−58.02 kcal/mol** | **−45.56 kcal/mol** |
| ΔG_SA (burial) | −6.22 kcal/mol | −6.34 kcal/mol |

**Key driving residues:**
- **ZINC52448673:** TRP15 (−17.0), GLU59 (−12.8), PHE60 (−4.2) — anchor residue binding mode
- **ZINC55760081:** ASP93 (−16.2), GLU59 (−6.6), LYS197 (−5.5) — more distributed interaction profile; comparable vdW to compound A

> **MM-GBSA caveat:** Full MM-GBSA (with GB solvation) requires energy minimization of the protein-ligand complex. Single-point GB calculations on docked poses produce unphysical desolvation penalties due to suboptimal geometry. The MM interaction energy (ΔE_MM) reported here is the reliable metric. For rigorous ΔG_bind, the recommended next step is a short energy minimization (500–1000 steps) followed by single-point MM-GBSA, or ensemble MM-GBSA from MD snapshots.

---

## Figures

| File | Description |
|------|-------------|
| `PCA_plot.png` | PCA of 66 samples (PC1=46.5%, PC2=6.9%); colored by condition, shaped by batch |
| `Volcano_plot.png` | Volcano plot: 2,107 up (red), 135 down (blue), 16,706 NS (grey) |
| `MA_plot.png` | log2 mean expression vs log2FC; significant DEGs highlighted |
| `Heatmap_top50_DEGs.png` | Top 50 DEGs, Z-scored, annotated by condition + batch |
| `DEG_landscape_analysis.png` | 4-panel: LFC distribution, gene families, pathway centrality scatter, B-cell heatmap |
| `GSEA_dotplot.png` | GSEA dotplot (top 20 activated + 20 suppressed pathways) |
| `GSEA_running_score.png` | GSEA running score plots for top 3 Hallmark pathways |
| `ORA_upregulated_barplot.png` | ORA barplot (top 25 enriched pathways in upregulated DEGs) |
| `pathway_jaccard_heatmap.png` | Hallmark pathway Jaccard similarity heatmap (4 clusters) |
| `pathway_NES_barplot.png` | Hallmark NES barplot by cluster |
| `pathway_cross_analysis.png` | Cross-database pathway comparison |
| `leading_edge_genes_bubble.png` | Top 25 leading-edge genes bubble chart |
| `drug_target_prioritization.png` | Drug target ranking (stacked bar + summary table) |
| `IL12B_deep_dive.png` | IL12B expression, co-expression network, and target profile |
| `batch_docking_validation.svg/.png` | 3-panel batch docking validation figure (10 compounds) |
| `admet_analysis.svg/.png` | 6-panel ADMET comparison (ZINC52448673 vs ZINC55760081) |
| `mm_binding_energy.svg/.png` | MM binding energy components + per-residue decomposition |

---

## Output Files

| File | Description |
|------|-------------|
| `DEG_table_HT_vs_NonHT.csv` | 2,242 significant DEGs (padj<0.05, \|log2FC\|≥1) with full statistics |
| `enrichment_gsea_results.csv` | All 1,006 significant GSEA gene sets (FDR<0.05) |
| `enrichment_ora_up_results.csv` | ORA results for upregulated DEGs (4,681 rows; 867 with FDR<0.05) |
| `enrichment_ora_down_results.csv` | ORA results for downregulated DEGs (0 significant at FDR<0.05) |
| `drug_target_ranking.csv` | Full ranked table of 46 drug target candidates with composite scores |
| `batch_docking_ranked.csv` | Ranked docking results for all 10 screened compounds |
| `nextflow_pipeline_HRA001684.zip` | Complete Nextflow DSL2 pipeline for HPC processing from raw FASTQs |
| `execution_trace.ipynb` | Full computational notebook with all code, outputs, and intermediate results |

---

## Nextflow Pipeline

A complete Nextflow DSL2 pipeline is provided in `nextflow_pipeline_HRA001684.zip` for processing raw FASTQ files from scratch.

**Pipeline steps:** FastQC → MultiQC → Trimmomatic → STAR (GRCh38) → featureCounts (Gencode v38) → DESeq2

**Quick start (SLURM HPC):**

```bash
unzip nextflow_pipeline_HRA001684.zip
cd nextflow_pipeline

# Download FASTQs from NGDC
bash bin/download_hra001684.sh --outdir data/fastq

# Build STAR index (one-time, ~30 GB, ~2 h)
STAR --runMode genomeGenerate \
     --genomeDir reference/star_index_GRCh38_gencode38 \
     --genomeFastaFiles reference/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile reference/gencode.v38.annotation.gtf \
     --runThreadN 16

# Run pipeline
nextflow run main.nf -profile slurm -params-file params.yaml
```

> **Important:** Verify sample-to-condition assignments in `params.yaml` against the NGDC HRA001684 metadata before running. The accession contains 74 runs but the paper reports 66 samples; 8 runs may be technical replicates or QC-excluded samples. The pipeline assigns batches by run order (groups of 12) as a proxy — replace with actual batch metadata from the NGDC submission.

**Nextflow pipeline containers:**
- FastQC: `quay.io/biocontainers/fastqc:0.12.1`
- Trimmomatic: `quay.io/biocontainers/trimmomatic:0.39`
- STAR: `quay.io/biocontainers/star:2.7.11a`
- featureCounts (Subread): `quay.io/biocontainers/subread:2.0.6`
- DESeq2: `bioconductor/bioconductor_docker:RELEASE_3_18`

---

## Software & Versions

| Software | Version | Purpose |
|----------|---------|---------|
| R | 4.4.x | Core analysis environment |
| clusterProfiler | 4.14.0 | GSEA and ORA |
| enrichplot | 1.26.1 | Enrichment visualization |
| msigdbr | 7.5.x | MSigDB gene sets |
| ComplexHeatmap | 2.22.0 | Heatmap generation |
| ggplot2 | 3.x | Visualization |
| ggrepel | 0.9.x | Label collision avoidance |
| MSigDB | v2023.2 | Hallmark, KEGG Legacy, GO:BP |
| AutoDock Vina | 1.2.5 | Molecular docking |
| RDKit | latest | Cheminformatics |
| OpenBabel | latest | PDBQT → SMILES conversion |
| DeepPurpose | latest | ADMET prediction (MPNN) |
| PDBFixer | latest | Protein structure preparation |
| Open Targets Platform | v26.03 | Drug target druggability data |

---

## Limitations

1. **No raw count re-analysis:** DESeq2 was not re-run from raw counts (unavailable). The authors' pre-computed statistics were used. Results cannot be independently validated at the count level.

2. **Heatmap normalization:** log2(size-factor normalized + 1) was used instead of VST (raw counts required for VST). Z-scoring mitigates this substantially.

3. **GSEA ranking metric:** `sign(log2FC) × −log10(padj)` used instead of the DESeq2 Wald statistic (not available in supplementary data). A 4.32% tie rate in the ranked list was flagged by fgsea.

4. **ORA for downregulated genes:** Underpowered (n=135 genes). No significant pathways detected; this is a power limitation, not evidence of no enrichment.

5. **Sample imbalance:** 16 HT vs 50 Non-HT. The imbalance is inherent to the study design and reduces power for detecting HT-specific signals.

6. **Bulk RNA-seq limitation:** Expression reflects the average of all cell types in the thyroid biopsy. Many "upregulated" genes reflect cell-type composition changes (more lymphocytes) rather than transcriptional changes within a given cell type. Single-cell RNA-seq would be needed to resolve cell-type-specific expression.

7. **Drug target scoring:** The composite scoring framework uses subjective weights. IFNG's overranking demonstrates that disease association scores from Open Targets can conflate related but distinct phenotypes. Scores should be interpreted as a prioritization aid, not a definitive ranking.

8. **Docking reproducibility:** AutoDock Vina uses stochastic sampling. Pose reproducibility was confirmed using exhaustiveness=32, seed=42 for the top 2 leads.

9. **MM-GBSA without minimization:** Single-point MM-GBSA on docked poses is unreliable for absolute ΔG values. The MM interaction energy (ΔE_MM) is reported as the reliable metric.

10. **No functional validation:** All findings are computational. Experimental validation (e.g., IL-12 blockade in a murine HT model, or analysis of ustekinumab-treated patients with comorbid HT) would be required before clinical translation.

---

## Citation

If using these results, the pipeline, or the data files, please cite the original dataset:

```
Zhang Y, et al. (2022). Single-cell RNA sequencing identifies key transcriptional
signatures characterizing the cellular landscape of Hashimoto's thyroiditis.
Nature Communications 13, 810. https://doi.org/10.1038/s41467-022-28120-2
```

And the tools used:

```
DESeq2:         Love MI, et al. (2014) Genome Biology 15:550
clusterProfiler: Yu G, et al. (2012) OMICS 16(5):284-7
GSEA:           Subramanian A, et al. (2005) PNAS 102(43):15545-50
MSigDB:         Liberzon A, et al. (2015) Cell Systems 1(6):417-25
ComplexHeatmap: Gu Z (2022) iMeta 1(3):e43
AutoDock Vina:  Eberhardt J, et al. (2021) J Chem Inf Model 61(8):3891-3898
Open Targets:   Ochoa D, et al. (2023) Nucleic Acids Res 51(D1):D1353-D1359
```

---

*This project was conducted entirely by **Biomni**, the AI research agent on the [Phylo](https://phylo.bio) platform. Analysis date: April 2026.*
