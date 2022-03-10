# wfmash based liftover of GRCh38 onto CHM13

This is an experimental liftover process designed as a counterpoint to methods based on minimap2.

We use `wfmash` to generate the alignment, and `paf2chain` to convert it to a chain file.

We used `wfmash` version `a36ab5f`, with the specific guix build `/gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash`.

```shell
run_wfmash=/gnu/store/6h8zlg7kbiidsmin62bbg372in2l3wkb-wfmash-0.7.0+a36ab5f-24/bin/wfmash
run_rustybam=/home/guarracino/tools/rustybam/target/release/rustybam
run_paf2chain=/home/guarracino/tools/paf2chain/target/release/paf2chain-f68eecaade2f9a0c7adfb8baf822b5a5865594a0
```

## aligning with 5 kb segments and 95% identity and exploring different block lengths (wfmash defaults are 5 kb segment length, 95% identity, and 25 kb block length)

GRCh38 onto CHM13 as a reference. First autosomes and X, then Y.

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"
        i=grch38_onto_chm13_autosome_X_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        sbatch -p workers -c 48 --wrap '\time -v '$run_wfmash' -t 48 -p 95 -s '$s' -l '$l' chm13v2.0_chr_auto_X.fasta GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_chr_auto_X.fasta >'$i.paf
            
        j=grch38_onto_chm13_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        sbatch -p workers -c 48 --wrap '\time -v '$run_wfmash' -t 48 -p 95 -s '$s' -l '$l' chm13v2.0_chrY.fasta chrY_maskedcentromeres.fa >'$j.paf
    done
done
```

Same but CHM13 as query and GRCh38 as reference.

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"
    
        i=chm13_onto_grch38_autosome_X_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        sbatch -p workers -c 48 --wrap '\time -v '$run_wfmash' -t 48 -p 95 -s '$s' -l '$l' GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_chr_auto_X.fasta chm13v2.0_chr_auto_X.fasta >'$i.paf
    
        j=chm13_onto_grch38_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        sbatch -p workers -c 48 --wrap '\time -v '$run_wfmash' -t 48 -p 95 -s '$s' -l '$l' chrY_maskedcentromeres.fa chm13v2.0_chrY.fasta >'$j.paf
    done
done
```

We then merged the two PAFs for each.

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"

        i=grch38_onto_chm13_autosome_X_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        j=grch38_onto_chm13_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        cat $i.paf $j.paf >grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.paf
        
        i=chm13_onto_grch38_autosome_X_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        j=chm13_onto_grch38_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l};
        cat $i.paf $j.paf >chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.paf
    done
done
```

## filtering with rustybam

We apply a filtering strategy from Mitchell Vollger based on [rustybam](https://mrvollger.github.io/rustybam/).
This forces same-chromosome mappings, and then breaks and trims the PAF to represent only 1:1 mappings.

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"
    
        awk '$1==$6' grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.paf |\
          $run_rustybam break-paf --max-size 10000  |\
          $run_rustybam trim-paf -r |\
          $run_rustybam invert |\
          $run_rustybam trim-paf -r |\
          $run_rustybam invert > grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf
        
        awk '$1==$6' chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.paf |\
          $run_rustybam break-paf --max-size 10000  |\
          $run_rustybam trim-paf -r |\
          $run_rustybam invert |\
          $run_rustybam trim-paf -r |\
          $run_rustybam invert > chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf
    done
done
```

## obtaining the chains

We convert the PAF files to CHAIN files with `paf2chain`:

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"
      
        $run_paf2chain -i grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf > grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.chain
        $run_paf2chain -i chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf > chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.chain
    done
done
```


## quality inspection

With `SafFire`:

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"
    
        $run_rustybam stats grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf -p > grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.stats
        $run_rustybam stats chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf -p > chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.stats
    done
done
```

With a few statistics:

```shell
for s in 5k; do   
    for x in 1 3 5; do
        s_no_k=${s::-1}
        l_no_k=$(echo $s_no_k '*' $x | bc)
        l=${l_no_k}k
          
        echo "-s $s -l $l"
    
        f=grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf
        echo "$f"
        cat $f | sed s/gi:f:// | awk '{ sum += $13 * $10; tot += $10; block += $11; } END { print "matches",tot; print "block",block; print "gap.id", sum/tot; print "block.id",tot/block * 100; }';
        echo -n 'ins/del gaps '; cat $f | egrep -o '[[:digit:]]+I[[:digit:]]+D' | tr 'I' '\n' | tr -d 'D' | awk '{ s+= 1 } END { print s }'
        
        f=chm13_onto_grch38_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s${s}_l${l}.trim.paf
        echo "$f"
        cat $f | sed s/gi:f:// | awk '{ sum += $13 * $10; tot += $10; block += $11; } END { print "matches",tot; print "block",block; print "gap.id", sum/tot; print "block.id",tot/block * 100; }';
        echo -n 'ins/del gaps '; cat $f | egrep -o '[[:digit:]]+I[[:digit:]]+D' | tr 'I' '\n' | tr -d 'D' | awk '{ s+= 1 } END { print s }'
    done
done
```

We shared the results of `-s 5k -l 25k -p 95` on Globus, at `/team-liftover/v1_wfmash/`.
