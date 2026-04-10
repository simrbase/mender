#!/bin/bash
## ============================================================
## USEFUL QUERIES FOR SPLIT GENE OUTPUT FILES
## Run this script to get a summary of all four output files
## or copy individual sections as needed
## ============================================================
##
## Files:
##   split_genes_summary.txt   -- one row per gene pair, best hit only
##   split_genes_detail.txt    -- one row per tiling hit per pair
##   split_genes_merge.txt     -- one row per merge candidate (chained pairs)
##   validated_merge.txt       -- merge table + IsoSeq spanning evidence
##
## Column reference:
##   split_genes_summary.txt
##     1:g1  2:g2  3:g1_coord  4:g2_coord  5:genomic_dist
##     6:g1_desc  7:g2_desc  8:best_hit  9:hit_desc
##     10:combined_cov_pct  11:tiling_gap  12:g1_subj  13:g2_subj
##     14:g1_cov  15:g2_cov  16:g1_evalue  17:g2_evalue  18:num_tiling_hits
##
##   split_genes_detail.txt
##     1:g1  2:g2  3:genomic_dist  4:hit  5:hit_desc
##     6:combined_cov_pct  7:tiling_gap  8:g1_subj  9:g2_subj
##     10:g1_cov  11:g2_cov  12:g1_evalue  13:g2_evalue
##
##   split_genes_merge.txt
##     1:merge_id  2:num_genes  3:genes_in_order  4:coords_in_order  5:ref
##     6:gene_descs  7:best_hit  8:hit_desc  9:best_combined_cov_pct
##     10:best_gap  11:evalues  12:junction_tiling_hits
##     13:junction_genomic_dist  14:skipped_genes
##     15:max_tiling_hits  16:num_pairs_in_chain  17:flag
##
##   validated_merge.txt (merge table + 3 extra columns)
##     1:merge_id  2:num_genes  3:genes_in_order  4:coords_in_order  5:ref
##     6:gene_descs  7:best_hit  8:hit_desc  9:best_combined_cov_pct
##     10:best_gap  11:evalues  12:junction_tiling_hits
##     13:junction_genomic_dist  14:skipped_genes
##     15:max_tiling_hits  16:num_pairs_in_chain  17:flag
##     18:spanning_isoseq_count  19:spanning_isoseq_detail  20:isoseq_flag
##
## ============================================================

echo ""
echo "============================================================"
echo "SPLIT_GENES_SUMMARY.TXT"
echo "One row per gene pair. Best supporting hit only."
echo "============================================================"

echo ""
echo "--- Total gene pairs reported ---"
wc -l split_genes_summary.txt

echo ""
echo "--- Distribution of num_tiling_hits (confidence score) ---"
echo "How many different reference proteins independently confirm the same tiling split."
echo "Both genes must share a reference protein AND their alignments must tile end-to-end"
echo "on that protein. A high count means many reference proteins all show the same"
echo "split boundary — consistent with a domain boundary conserved across paralogs and"
echo "orthologs. A low count means few references happen to span both halves."
echo "  1 hit  = thin evidence — rely on IsoSeq for confidence"
echo "  2-4    = moderate, consistent tiling across a small set of references"
echo "  5+     = strong, conserved domain boundary confirmed across many proteins"
echo "  10+    = very strong"
echo "Note: low hit counts can reflect phylogenetic distance to the reference or"
echo "a small gene family — not necessarily a false positive. Always cross-check"
echo "with IsoSeq for 1-hit cases. If two genes preferentially hit DIFFERENT"
echo "reference proteins (paralogs with diverged domains), their pair is invisible"
echo "to this method entirely — a known blind spot."
awk -F'\t' 'NR>1 {print $18}' split_genes_summary.txt | sort -n | uniq -c

echo ""
echo "--- Pairs with strong tiling support (>= 5 independent hits) ---"
awk -F'\t' 'NR>1 && $18 >= 5' split_genes_summary.txt | wc -l

echo ""
echo "--- Only directly adjacent pairs (genomic_dist=1, no gene between them) ---"
echo "The remainder have an annotated gene sitting between them (check skipped_genes"
echo "in the merge table — it may be another fragment of the same split gene)"
awk -F'\t' 'NR>1 && $5==1' split_genes_summary.txt | wc -l

echo ""
echo "--- Pairs with genomic_dist > 1 (a gene sits between them) ---"
awk -F'\t' 'NR>1 && $5>1' split_genes_summary.txt | wc -l

echo ""
echo "--- Top 10 pairs by num_tiling_hits ---"
echo "g1 | g2 | hit_desc | combined_cov_pct | tiling_gap | num_tiling_hits"
sort -t$'\t' -k18 -rn split_genes_summary.txt | \
    awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$9"\t"$10"\t"$11"\t"$18}' | head -10

echo ""
echo "--- Top 10 pairs by combined_cov_pct ---"
echo "g1 | g2 | hit_desc | combined_cov_pct | num_tiling_hits"
sort -t$'\t' -k10 -rn split_genes_summary.txt | \
    awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$9"\t"$10"\t"$18}' | head -10

echo ""
echo "============================================================"
echo "SPLIT_GENES_MERGE.TXT"
echo "One row per merge candidate. Pairs sharing a gene are chained."
echo "Trimming and splitting applied based on junction asymmetry."
echo "============================================================"

echo ""
echo "--- Total merge candidates ---"
wc -l split_genes_merge.txt

echo ""
echo "--- Size distribution (num_genes per merge candidate) ---"
echo "Most should be 2-gene pairs; multi-gene chains are most interesting"
awk -F'\t' 'NR>1 {print $2}' split_genes_merge.txt | sort -n | uniq -c

echo ""
echo "--- Flag distribution ---"
echo "noflag        = directly adjacent, 2+ tiling hits, coverage >= 60%"
echo "STRONG        = all junctions >= 3 hits (chains only)"
echo "SINGLE_HIT    = all junctions have exactly 1 tiling hit"
echo "WEAK_END      = low-evidence terminal junction survived trimming"
echo "WEAK_INTERNAL = low-evidence internal junction survived splitting"
echo "LOW_COV       = combined coverage < 60%"
echo "SKIPPED_GENE  = a non-adjacent gene sits inside this merge locus"
echo ""
echo "Note: '?' in junction_tiling_hits means two consecutive genes in the chain"
echo "have NO direct tiling pair — they joined transitively through a shared neighbor."
echo "Possible causes (all same-strand since strand filter is applied at pair build time):"
echo "  - the two genes hit different reference proteins (empty intersection)"
echo "  - their alignments heavily overlap on the same reference (gap >> wiggle)"
echo "Both indicate these two consecutive genes are not direct tiling fragments."
awk -F'\t' 'NR>1 {print $17}' split_genes_merge.txt | sort | uniq -c | sort -rn

echo ""
echo "--- STRONG candidates (highest confidence protein evidence, chains only) ---"
awk -F'\t' 'NR>1 && $17~/STRONG/' split_genes_merge.txt | wc -l

echo ""
echo "--- Top 20 STRONG candidates by max_tiling_hits ---"
echo "merge_id | num_genes | hit_desc | combined_cov | junction_hits | max_tiling | flag"
awk -F'\t' 'NR>1 && $17~/STRONG/' split_genes_merge.txt | \
    sort -t$'\t' -k15 -rn | \
    awk -F'\t' '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$15"\t"$17}' | head -20

echo ""
echo "--- Multi-gene chains (3+ genes) ---"
echo "These are the most biologically interesting — multiple consecutive split genes"
awk -F'\t' 'NR>1 && $2 > 2' split_genes_merge.txt | wc -l

echo ""
echo "--- Multi-gene chains that are STRONG ---"
echo "merge_id | num_genes | hit_desc | combined_cov | junction_hits | junction_dists | max_tiling | flag"
awk -F'\t' 'NR>1 && $2 > 2 && $17~/STRONG/' split_genes_merge.txt | \
    cut -f1,2,8,9,12,13,15,17

echo ""
echo "--- SKIPPED_GENE candidates (gene sits inside merge locus) ---"
echo "The skipped gene may be another fragment that failed filters"
echo "Check skipped_genes column in GFF and blast output"
awk -F'\t' 'NR>1 && $17~/SKIPPED_GENE/' split_genes_merge.txt | wc -l

echo ""
echo "--- SKIPPED_GENE candidates with STRONG protein evidence ---"
echo "merge_id | num_genes | hit_desc | combined_cov | skipped_genes | max_tiling | flag"
awk -F'\t' 'NR>1 && $17~/STRONG/ && $17~/SKIPPED_GENE/' split_genes_merge.txt | \
    cut -f1,2,8,9,14,15,17

echo ""
echo "--- WEAK_END candidates (terminal gene may not belong) ---"
echo "Review junction_tiling_hits and junction_genomic_dist to see which end is weak"
awk -F'\t' 'NR>1 && $17~/WEAK_END/' split_genes_merge.txt | wc -l

echo ""
echo "--- Junction tiling hits for all multi-gene chains ---"
echo "merge_id | num_genes | hit_desc | junction_hits | junction_dists | max_tiling | flag"
awk -F'\t' 'NR>1 && $2 > 2' split_genes_merge.txt | \
    cut -f1,2,8,12,13,15,17 | head -30

echo ""
echo "--- Chains with internal 1-hit junctions (potential mis-chaining) ---"
echo "Format: merge_id | num_genes | hit_desc | junction_hits | junction_dists | flag"
awk -F'\t' 'NR>1 && $12!="NA"' split_genes_merge.txt | \
    awk -F'\t' '{
        n=split($12,a,"|");
        for(i=2;i<=n-1;i++) {
            if(a[i]+0==1) {
                print $1"\t"$2"\t"$8"\t"$12"\t"$13"\t"$17;
                break
            }
        }
    }'

echo ""
echo "--- Chains with ? junctions (transitive joins with no direct tiling evidence) ---"
echo "Two consecutive genes joined by no direct tiling pair."
echo "They connected via a shared neighbor gene, not a direct protein hit."
echo "Causes: different reference proteins (empty hit intersection), or"
echo "heavy alignment overlap on the same reference (gap >> wiggle)."
echo "Format: merge_id | num_genes | hit_desc | junction_hits | max_tiling | flag"
awk -F'\t' 'NR>1 && $12~/\?/' split_genes_merge.txt | \
    cut -f1,2,8,12,15,17

echo ""
echo "--- Chains where a terminal gene has no hit to the chain's best reference ---"
echo "? in the evalue column means that gene contributed NO direct tiling evidence"
echo "to the chain's best reference protein. It joined via a different protein entirely."
echo "This is the strongest indicator that a terminal gene may not belong in the merge."
echo "Format: merge_id | num_genes | hit_desc | evalues | junction_hits | flag"
awk -F'\t' 'NR>1' split_genes_merge.txt | awk -F'\t' '{
    n=split($11,ev,",");
    if (ev[1]=="?" || ev[n]=="?")
        print $1"\t"$2"\t"$8"\t"$11"\t"$12"\t"$17
}'

echo ""
echo "============================================================"
echo "SPLIT_GENES_DETAIL.TXT"
echo "One row per supporting tiling hit per pair."
echo "Use this to see all evidence for a specific pair."
echo "============================================================"

echo ""
echo "--- Total tiling hit records ---"
wc -l split_genes_detail.txt

echo ""
echo "--- Unique gene pairs with detail records ---"
awk -F'\t' 'NR>1 {print $1"\t"$2}' split_genes_detail.txt | sort -u | wc -l

echo ""
echo "--- Top 10 hits by combined_cov_pct across all pairs ---"
echo "g1 | g2 | hit_desc | combined_cov_pct | tiling_gap"
sort -t$'\t' -k6 -rn split_genes_detail.txt | \
    awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | head -10

echo ""
echo "============================================================"
echo "VALIDATED_MERGE.TXT"
echo "Merge table with IsoSeq spanning read evidence added."
echo "Note: absence of spanning reads is NOT a red flag —"
echo "IsoSeq came from specific developmental stages."
echo "PARTIAL_SPAN is the meaningful structural signal."
echo "============================================================"

echo ""
echo "--- IsoSeq flag distribution ---"
echo "FULL_SPAN    = at least one read spans all genes — strong confirmation"
echo "PARTIAL_SPAN = reads exist but dont reach a terminal gene — structural signal"
echo "               meaningful even with stage-specific IsoSeq — a read long"
echo "               enough to span genes 1+2 should also cross gene 3 if they"
echo "               are truly one transcriptional unit"
echo "none         = no spanning reads found (may just be expression stage)"
awk -F'\t' 'NR>1 {print $20}' validated_merge.txt | sort | uniq -c

echo ""
echo "--- Spanning read count distribution ---"
echo "How many candidates have how many spanning reads"
awk -F'\t' 'NR>1 {print $18}' validated_merge.txt | sort -n | uniq -c | head -20

echo ""
echo "--- Flag distribution in validated merge ---"
awk -F'\t' 'NR>1 {print $17}' validated_merge.txt | sort | uniq -c | sort -rn

echo ""
echo "--- Candidates by flag + isoseq_flag combination ---"
echo "Shows full picture of protein evidence vs IsoSeq evidence"
awk -F'\t' 'NR>1 {print $17"\t"$20}' validated_merge.txt | \
    sort | uniq -c | sort -rn | head -20

echo ""
echo "--- Gold standard candidates ---"
echo "STRONG protein + FULL_SPAN IsoSeq + >=3 spanning reads + >=80% coverage"
awk -F'\t' 'NR>1 && $17~/STRONG/ && $9+0>=80 && $20=="FULL_SPAN" && $18+0>=3' \
    validated_merge.txt | wc -l

echo ""
echo "--- Gold standard list sorted by spanning read count ---"
echo "merge_id | num_genes | hit_desc | combined_cov | max_tiling | flag | spanning | isoseq_flag"
awk -F'\t' 'NR>1 && $17~/STRONG/ && $9+0>=80 && $20=="FULL_SPAN" && $18+0>=3' \
    validated_merge.txt | sort -t$'\t' -k18 -rn | \
    cut -f1,2,8,9,15,17,18,20 | head -30

echo ""
echo "--- Recommended merge set (what merge_split_genes.pl would process) ---"
echo "FULL_SPAN or none isoseq, skip SKIPPED_GENE and LOW_COV"
awk -F'\t' 'NR>1 && $17!~/SKIPPED_GENE/ && $17!~/LOW_COV/' \
    validated_merge.txt | wc -l

echo ""
echo "--- Breakdown of recommended set by isoseq_flag ---"
awk -F'\t' 'NR>1 && $17!~/SKIPPED_GENE/ && $17!~/LOW_COV/' \
    validated_merge.txt | awk -F'\t' '{print $20}' | sort | uniq -c

echo ""
echo "--- PARTIAL_SPAN cases — terminal gene may not belong in merge ---"
echo "Especially important when protein flag is also STRONG"
echo "merge_id | num_genes | genes | hit_desc | max_tiling | flag | spanning_count | isoseq_flag"
awk -F'\t' 'NR>1 && $20=="PARTIAL_SPAN"' validated_merge.txt | \
    cut -f1,2,3,8,15,17,18,20

echo ""
echo "--- WEAK_END chains where the terminal gene has no hit to chain's best ref ---"
echo "These are the highest priority WEAK_END cases: protein evidence AND"
echo "evalue both say the terminal gene does not belong. Manual review recommended."
echo "Format: merge_id | num_genes | hit_desc | evalues | junction_hits | flag | isoseq_flag"
awk -F'\t' 'NR>1 && $17~/WEAK_END/' validated_merge.txt | awk -F'\t' '{
    n=split($11,ev,",");
    if (ev[1]=="?" || ev[n]=="?")
        print $1"\t"$2"\t"$8"\t"$11"\t"$12"\t"$17"\t"$20
}'

echo ""
echo "--- STRONG protein but PARTIAL_SPAN IsoSeq ---"
echo "These are your highest priority manual review cases"
echo "Protein says merge all genes, IsoSeq says a terminal gene may not belong"
awk -F'\t' 'NR>1 && $17~/STRONG/ && $20=="PARTIAL_SPAN"' validated_merge.txt | \
    cut -f1,2,3,8,9,15,17,18,20

echo ""
echo "--- Multi-gene chains with IsoSeq support ---"
echo "merge_id | num_genes | hit_desc | junction_hits | max_tiling | flag | spanning | isoseq_flag"
awk -F'\t' 'NR>1 && $2>2 && $18+0>0' validated_merge.txt | \
    sort -t$'\t' -k18 -rn | \
    cut -f1,2,8,12,15,17,18,20 | head -20

echo ""
echo "--- Candidates without IsoSeq support ---"
echo "May still be real — check if gene is expressed in your IsoSeq stage"
awk -F'\t' 'NR>1 && $18==0' validated_merge.txt | wc -l

echo ""
echo "--- SKIPPED_GENE in validated merge ---"
echo "merge_id | num_genes | hit_desc | skipped_genes | max_tiling | flag | spanning | isoseq_flag"
awk -F'\t' 'NR>1 && $17~/SKIPPED_GENE/' validated_merge.txt | \
    cut -f1,2,8,14,15,17,18,20 | head -20

echo ""
echo "============================================================"
echo "GENE LOOKUP EXAMPLES"
echo "Replace CCA3gXXXXXX with your gene of interest"
echo "Replace GENENAME with gene symbol e.g. JMJD1C"
echo "Replace merge_NNN with merge ID e.g. merge_502"
echo "============================================================"

echo ""
echo "--- Find a gene in all files ---"
echo "  grep CCA3gXXXXXX split_genes_summary.txt"
echo "  grep CCA3gXXXXXX split_genes_detail.txt"
echo "  grep CCA3gXXXXXX split_genes_merge.txt"
echo "  grep CCA3gXXXXXX validated_merge.txt | cut -f1,2,3,8,15,17,18,20"

echo ""
echo "--- Find a gene by name ---"
echo "  grep GENENAME split_genes_merge.txt | cut -f1,2,3,8,9,12,13,14,15,17"
echo "  grep GENENAME validated_merge.txt   | cut -f1,2,3,8,15,17,18,20"

echo ""
echo "--- Trace a merge through all files ---"
echo "  grep ^merge_NNN  validated_merge.txt | cut -f1,2,3,8,15,17,18,20"
echo "  grep merge_id=merge_NNN new_merges.gff | cut -f1,3,4,5,9"
echo "  grep ^merge_NNN  new_merges.log"

echo ""
echo "--- Show full detail for a specific pair ---"
echo "  awk -F'\t' '\$1==\"CCA3gXXX\" && \$2==\"CCA3gYYY\"' split_genes_detail.txt"

echo ""
echo "--- Show all tiling hits for one gene ---"
echo "  grep CCA3gXXXXXX split_genes_detail.txt | cut -f1,2,4,5,6,7"

echo ""
echo "--- Check skipped genes in GFF ---"
echo "  grep SKIPPED_GENE_ID genes.gff | grep -P '\tgene\t'"

echo ""
echo "--- Check skipped gene in blast output ---"
echo "  grep SKIPPED_TRANSCRIPT_ID diamond.out | head -5"

echo ""
echo "--- Sanity check merged GFF gene counts ---"
echo "  grep -c \$'\\tgene\\t' original.gff     # before"
echo "  grep -c \$'\\tgene\\t' new_merges.gff    # after"
echo "  grep -c \$'\\tgene\\t' removed_genes.gff # removed"
echo "  # after = before - removed + merged_count"
