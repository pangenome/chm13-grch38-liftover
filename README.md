# wfmash based liftover of CHM13-T2T onto GRCh38

This is an experimental liftover process designed as a counterpoint to methods based on minimap2.

We use `wfmash` to generate the alignment, and `paf2chain` to convert it to a chain file.

No filtering is done, and so the simulated regions of GRCh38 are likely to match CHM13 in unusual ways.

We used `wfmash` version `26ca311`, with the specific guix build `/gnu/store/8i0kwx7xk7liiw57ihi36n7n3i6k6sp6-wfmash-0.7.0+26ca311-4/bin/wfmash`.

## using 5kb segments and 90% identity

```
i=19a; sbatch -p workers -c 48 --wrap 'wfmash -t 48 -p 90 -s 5k grch38.fa.gz chm13.fa.gz >'$i.paf
```

## using 1kb segments and 90% identity

```
i=23a; sbatch -p workers -c 48 --wrap 'wfmash -t 48 -p 90 -s 1k grch38.fa.gz chm13.fa.gz >'$i.paf
```

## downloads

Alignments in PAF format:

```
https://f004.backblazeb2.com/file/pangenome/T2T/liftover/chm13_vs_grch38_wfmash_p90s5k_19a.paf.gz
https://f004.backblazeb2.com/file/pangenome/T2T/liftover/chm13_vs_grch38_wfmash_p90s1k_23a.paf.gz
```

Chain files made by `paf2chain`:

```
https://f004.backblazeb2.com/file/pangenome/T2T/liftover/chm13_vs_grch38_wfmash_p90s5k_19a.chain.gz
https://f004.backblazeb2.com/file/pangenome/T2T/liftover/chm13_vs_grch38_wfmash_p90s1k_23a.chain.gz
```

## crude evaluation

We can examine approximately what fraction of GRCh38 is covered by these alignments using BED files in this repo:

We cover approximately 99.6% of GRCh38 with the 5kb segment length and >100% with the 1kb segment length (some regions of GRCh38 are mapped to by more than one region in CHM13).

```
-> % Rscript -e $(bedtools intersect -a <(zcat chm13_vs_grch38_wfmash_p90s5k_19a.paf.gz | awk '{ print $6, $8, $9 }' | tr ' ' '\t' | bedtools sort ) -b <(bedtools subtract -a grch38.fa.gz.fai.bed -b Modeled_regions_for_GRCh38.bed) | awk '{ sum += $3 - $2; } END { print sum }')'/'$(bedtools subtract -a grch38.fa.gz.fai.bed -b Modeled_regions_for_GRCh38.bed | bedtools subtract -a /dev/stdin -b grch38.fa.gz.Ns.bed | awk '{ sum += $3 - $2; } END { print sum }')
[1] 0.9957952
```

```
-> % Rscript -e $(bedtools intersect -a <(zcat chm13_vs_grch38_wfmash_p90s1k_23a.paf.gz | awk '{ print $6, $8, $9 }' | tr ' ' '\t' | bedtools sort ) -b <(bedtools subtract -a grch38.fa.gz.fai.bed -b Modeled_regions_for_GRCh38.bed) | awk '{ sum += $3 - $2; } END { print sum }')'/'$(bedtools subtract -a grch38.fa.gz.fai.bed -b Modeled_regions_for_GRCh38.bed | bedtools subtract -a /dev/stdin -b grch38.fa.gz.Ns.bed | awk '{ sum += $3 - $2; } END { print sum }')
[1] 1.016974
```
