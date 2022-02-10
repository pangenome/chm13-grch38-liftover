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
