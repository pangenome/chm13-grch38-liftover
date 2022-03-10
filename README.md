# wfmash based liftover of GRCh38 onto CHM13

This is an experimental liftover process designed as a counterpoint to methods based on minimap2.

We use `wfmash` to generate the alignment, and `paf2chain` to convert it to a chain file.

No filtering is done, and so the simulated regions of GRCh38 are likely to match CHM13 in unusual ways.

We used `wfmash` version `a36ab5f`, with the specific guix build `/gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash`.

## using 5kb segments and 95% identity (wfmash defaults)

GRCh38 onto CHM13 as a reference.
First autosomes and X, then Y.

```
i=grch38_onto_chm13_autosome_X_wfmash-0.7.0+a36ab5f_p95_s5k; sbatch -p workers -c 48 --wrap 'time /gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash -t 48 -p 95 -s 5k chm13v2.0_chr_auto_X.fasta GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_chr_auto_X.fasta >'$i.paf
i=grch38_onto_chm13_Y_wfmash-0.7.0+a36ab5f_p95_s5k; sbatch -p workers -c 48 --wrap 'time /gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash -t 48 -p 95 -s 5k chm13v2.0_chrY.fasta chrY_maskedcentromeres.fa >'$i.paf
```

Same but CHM13 as query and GRCh38 as reference.

```
i=chm13_onto_grch38_autosome_X_wfmash-0.7.0+a36ab5f_p95_s5k; sbatch -p workers -c 48 --wrap 'time /gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash -t 48 -p 95 -s 5k GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_chr_auto_X.fasta chm13v2.0_chr_auto_X.fasta >'$i.paf
i=chm13_onto_grch38_Y_wfmash-0.7.0+a36ab5f_p95_s5k; sbatch -p workers -c 48 --wrap 'time /gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash -t 48 -p 95 -s 5k chrY_maskedcentromeres.fa chm13v2.0_chrY.fasta >'$i.paf
```

We then merged the two PAFs for each.

```
cat grch38_onto_chm13_autosome_X_wfmash-0.7.0+a36ab5f_p95_s5k.paf grch38_onto_chm13_Y_wfmash-0.7.0+a36ab5f_p95_s5k.paf >grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k.paf
cat chm13_onto_grch38_autosome_X_wfmash-0.7.0+a36ab5f_p95_s5k.paf chm13_onto_grch38_Y_wfmash-0.7.0+a36ab5f_p95_s5k.paf >chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k.paf
```

## filtering with rustybam

We apply a filtering strategy from Mitchell Vollger based on [rustybam](https://mrvollger.github.io/rustybam/).
This forces same-chromosome mappings, and then breaks and trims the PAF to represent only 1:1 mappings.

```
awk '$1==$6' grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k.paf | rb break-paf --max-size 10000  | rb trim-paf -r | rb invert | rb trim-paf -r | rb invert > grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k.trim.paf
awk '$1==$6' chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k.paf | rb break-paf --max-size 10000  | rb trim-paf -r | rb invert | rb trim-paf -r | rb invert > chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k.trim.paf
```

These files are included here, gzipped.
