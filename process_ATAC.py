#!/usr/bin/env python
# coding: utf-8

"""
Script: process_combined.py
Author: Luisa Amaral
Description:
    This script processes raw scATAC-seq data using SnapATAC2, leveraging RNA-based annotations
    from multiome (female) samples. It includes:
      - Fragment & h5ad creation
      - Metadata annotation transfer from RNA to ATAC
      - Subclustering
      - Peak calling per cell type
      - Merging peaks
      - Differential accessibility analysis
      - UMAP & visualization
"""

# === SECTION 1: Imports and Setup ===
import snapatac2 as snap
import scanpy as sc
import pandas as pd
import numpy as np
import os
from pathlib import Path

# === SECTION 2: Create Fragment & h5ad Files ===
def preprocess_fragment(bam_file, sample_name, output_dir):
    frag_file = f"{output_dir}/{sample_name}.frag"
    h5ad_file = f"{output_dir}/{sample_name}.h5ad"
    tsse_file = f"{output_dir}/{sample_name}.tsse.pdf"
    frag_plot_file = f"{output_dir}/{sample_name}.frag.pdf"
    log_file_path = f"{output_dir}/{sample_name}_out.txt"

    if not os.path.exists(frag_file):
        snap.pp.make_fragment_file(bam_file, barcode_tag='BX', output_file=frag_file)

    if not os.path.exists(h5ad_file):
        with open(log_file_path, "w") as file:
            data = snap.pp.import_data(frag_file, chrom_sizes=snap.genome.mm10, file=h5ad_file,
                                       chunk_size=250000, sorted_by_barcode=False)
            snap.pl.frag_size_distr(data, interactive=False, out_file=frag_plot_file)
            snap.metrics.tsse(data, snap.genome.mm10)
            snap.pl.tsse(data, interactive=False, out_file=tsse_file)
            snap.pp.filter_cells(data, min_counts=1000, max_counts=100000, min_tsse=7)
            snap.pp.add_tile_matrix(data, bin_size=500)
            snap.pp.select_features(data)
            snap.pp.scrublet(data)
            snap.pp.filter_doublets(data)
            file.write(str(data))
            data.close()

# Preprocess male BAMs
track_male = pd.read_csv("male_files.csv")
for i, row in track_male.iterrows():
    preprocess_fragment(row['path'], row['Sample_id'], output_dir="male_h5ads")

# Preprocess female fragments
track_female = pd.read_csv("tracksheet_with_paths_oct1723.csv")
def preprocess_female_fragment(fragment_file, sample_name, nfrag, tsse, nfrag_max):
    h5ad_file = f"female_h5ads_filt_oct1723/{sample_name}.h5ad"
    if not os.path.exists(h5ad_file):
        data = snap.pp.import_data(fragment_file, chrom_sizes=snap.genome.mm10, file=h5ad_file,
                                   chunk_size=250000, sorted_by_barcode=False)
        snap.metrics.tsse(data, snap.genome.mm10)
        snap.pp.filter_cells(data, min_counts=nfrag, max_counts=nfrag_max, min_tsse=tsse)
        snap.pp.add_tile_matrix(data, bin_size=500)
        snap.pp.select_features(data)
        snap.pp.scrublet(data)
        snap.pp.filter_doublets(data)
        data.close()

for _, row in track_female.iterrows():
    frag_path = f"/{row['Path-luisa']}/outs/atac_fragments.tsv.gz"
    preprocess_female_fragment(frag_path, row['Sample_name'], row['nfrag_cutoff'], row['TSS_cutoff'], row['nfrag_max'])

# === SECTION 3: Combine All h5ads into One Dataset ===
adatas = []
for filename in os.listdir("male_h5ads"):
    if filename.endswith(".h5ad"):
        adatas.append(("Male:" + filename[:-5], Path("male_h5ads/" + filename)))
for filename in os.listdir("female_h5ads_filt_oct1723"):
    if filename.endswith(".h5ad"):
        adatas.append(("Female:" + filename[:-5], Path("female_h5ads_filt_oct1723/" + filename)))

combined = snap.AnnDataSet(adatas=adatas, filename="combined.h5ads")
combined.close()

# === SECTION 4: Integration and UMAP ===
data = snap.read_dataset("combined.h5ads")
data.obs['batch'] = [x.split(":")[0] for x in data.obs_names]
snap.pp.mnc_correct(data, batch="batch")
snap.tl.spectral(data)
snap.tl.umap(data, use_rep="X_spectral_mnn")
snap.pp.knn(data, use_rep="X_spectral_mnn")
snap.tl.leiden(data, resolution=1.0)
data.obs['leiden_1'] = data.obs['leiden']
data.close()

# === SECTION 5: RNA Annotation Transfer and Subclustering ===
data = snap.read_dataset("combined.h5ads")
meta = pd.read_csv("meta.csv")
data.obs['best_celltype'] = meta['best_celltype']
data.obs['region'] = meta['region']
data.obs['age'] = meta['age']
data.obs['celltype_final'] = meta['celltype_final']
data.close()

# === SECTION 6: Peak Calling with MACS3 ===
data = snap.read_dataset("combined.h5ads")
snap.tl.macs3(data, n_jobs=20, replicate='rep', groupby='celltype_final', key_added='celltype_final_rep',
              blacklist="/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed")
peaks = snap.tl.merge_peaks(data.uns['celltype_final_rep'], chrom_sizes=snap.genome.mm10)
peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks['Peaks'], file="pmat_celltype_final_rep.h5ad")
data.close()

# === SECTION 7: Differential Accessibility ===
data = snap.read_dataset("combined.h5ads")
peak_mat = snap.read("pmat_celltype_final_rep.h5ad")
data.obs['celltype_age'] = data.obs['celltype_final'] + ":" + data.obs['age']
uniq_cts = data.obs['celltype_final'].unique()

for ct in uniq_cts:
    group1 = ct + ":2mo"
    group2 = ct + ":18mo"
    g1 = data.obs['celltype_age'] == group1
    g2 = data.obs['celltype_age'] == group2
    if g1.sum() < 25 or g2.sum() < 25:
        continue
    peaks_selected = peaks[ct].to_numpy()
    diff_peaks = snap.tl.diff_test(
        peak_mat,
        cell_group1=g1,
        cell_group2=g2,
        min_log_fc=0.15, min_pct=0.025,
        features=peaks_selected
    )
    if diff_peaks.shape[0] == 0:
        continue
    pd.DataFrame(diff_peaks).to_csv(f"diff_peaks_{ct.replace(' ', '_')}_2vs18.csv")

data.close()
peak_mat.close()

# === SECTION 8: Save Final Metadata and Matrices ===
data = snap.read_dataset("combined.h5ads")
data.obs['umap_x'] = data.obsm['X_umap'][:, 0]
data.obs['umap_y'] = data.obsm['X_umap'][:, 1]
data.obs.to_csv("final_metadata.csv")
data.close()

print("✅ Processing complete.")
