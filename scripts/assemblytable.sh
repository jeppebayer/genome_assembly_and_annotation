#!/bin/bash

repeatmaskfulltable="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/repeatmasking/RepeatMasker/EntNic.fasta.fullmask.tbl"
buscotable="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/HiC_scaffolding/juicer/3ddna_in10000_r6/finalize/busco_purged.nodebris.fasta/short_summary.specific.arthropoda_odb10.busco_purged.nodebris.fasta.txt"
annotationgtf="/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/annotation/braker2/braker.gtf"

buscotablevalues=$(
	awk \
	'BEGIN{
	FS = "\t"
	OFS = "|"
	}
	{
	if (NR == 2)
		{split($0, line2v1, ":")
		split(line2v1[2], line2v2, " ")
		buscoset = line2v2[1]}
	if ($2 ~ /^C:.*\[S:.*,D:.*\],F:.*,M:.*,n:[0-9]*/)
		{busco = $2}
	if ($3 == "Number of scaffolds")
		{nchrom = $2}
	if ($3 == "Number of contigs")
		{ncontigs = $2}
	if ($3 == "Total length")
		{genome_size = $2
		if (length(genome_size) > 3 && length(genome_size) <= 6)
			{unit_gs = "KB"
			genome_size_rounded = int(genome_size / 1000 * 10 + 0.5) / 10}
		if (length(genome_size) > 6 && length(genome_size) <= 9)
			{unit_gs = "MB"
			genome_size_rounded = int(genome_size / 1000000 * 10 + 0.5) / 10}
		if (length(genome_size) > 9)
			{unit_gs = "GB"
			genome_size_rounded = int(genome_size / 1000000000 * 10 + 0.5) / 10}
		}
	if ($3 == "Percent gaps")
		{gaps = $2}
	if ($3 == "Contigs N50")
		{split($2, n50info, " ")
		n50 = n50info[1]
		unit_n50 = n50info[2]}
	if ($2 ~ /^busco:/)
		{split($2, buscoversioninfo, " ")
		buscoversion = "v"buscoversioninfo[2]}
	}
	END{
	print unit_gs, genome_size_rounded, gaps, ncontigs, unit_n50, n50, nchrom, busco, buscoset, buscoversion
	}' \
	"$buscotable"
)

readarray -d "|" -t buscoarray <<< "$buscotablevalues"
unit_gs=${buscoarray[0]}
genome_size=${buscoarray[1]}
gaps=${buscoarray[2]}
ncontigs=${buscoarray[3]}
unit_n50=${buscoarray[4]}
n50=${buscoarray[5]}
nchrom=${buscoarray[6]}
busco=${buscoarray[7]}
buscoset=${buscoarray[8]}
buscoversion=${buscoarray[9]//$'\n'/}

annotationgtfvalues=$(
	awk \
	'function unit(size) {
		len = length(size)
		if (len > 0 && len <= 3)
			{return "BP"}
		if (len > 3 && len <= 6)
			{return "KB"}
		if (len > 6 && len <= 9)
			{return "MB"}
		if (len > 9)
			{return "GB"}
		}
	function roundnodecimal(num) {
		if (num - int(num) >= 0.5)
			{return int(num) + 1}
		else
			{return int(num)}
		}
	function round1decimalunit(num, unit) {
		if (unit == "KB")
			{return int(num / 1000 * 10 + 0.5) / 10}
		if (unit == "MB")
			{return int(num / 1000000 * 10 + 0.5) / 10}
		if (unit == "GB")
			{return int(num / 1000000000 * 10 + 0.5) / 10}
		}
	function round1decimal(num) {
		return int(num * 10 + 0.5) / 10
		}
	BEGIN{
	FS = "\t"
	OFS = "|"
	}
	{
	if ($3 == "gene")
		{gene_sum += 1}
	if ($9 ~ /.t1/)
		{if ($3 == "exon")
			{exon_length += $5 - $4 + 1
			exon_sum += 1}
		if ($3 == "intron")
			{intron_length += $5 - $4 + 1
			intron_sum += 1}
		}
	}
	END{
	print gene_sum, unit(roundnodecimal(exon_length / exon_sum)), roundnodecimal(exon_length / exon_sum), unit(exon_length), round1decimalunit(exon_length, unit(exon_length)), exon_sum, unit(roundnodecimal(intron_length / intron_sum)), roundnodecimal(intron_length / intron_sum), unit(intron_length), round1decimalunit(intron_length, unit(intron_length)), intron_sum, round1decimal(exon_sum / gene_sum)
	}' \
	"$annotationgtf"
)

readarray -d "|" -t annotationgtfarray <<< "$annotationgtfvalues"
ngenes=${annotationgtfarray[0]}
unit_mel=${annotationgtfarray[1]}
meanexon=${annotationgtfarray[2]}
unit_tel=${annotationgtfarray[3]}
totalexon=${annotationgtfarray[4]}
nexon=${annotationgtfarray[5]}
unit_mit=${annotationgtfarray[6]}
meanintron=${annotationgtfarray[7]}
unit_til=${annotationgtfarray[8]}
totalintron=${annotationgtfarray[9]}
nintron=${annotationgtfarray[10]}
exonpergene=${annotationgtfarray[11]//$'\n'/}

repeatmaskfulltablevalues=$(
	awk \
	'BEGIN{
	FS = "\t"
	OFS = "|"
	}
	{
	if ($1 == "bases masked:")
		{}
	}'
)

markdownformat() {
cat << EOF

| Genome Assembly |  |
| --- | :---: |
| Specimen | *$species_name* |
| Isolate | $isolate |
| Sequence coverage | $seq_coverage |
| Genome size ($unit_gs) | $genome_size |
| Gaps | $gaps |
| Number of contigs | $ncontigs |
| Contig N50 ($unit_n50) | $n50 |
| Number of chromosomes | $nchrom |
| Number of protein-coding genes | $ngenes |
| Mean exon length ($unit_mel) | $meanexon |
| Total exon length ($unit_tel) | $totalexon |
| Number of exons |  $nexon |
| Mean intron length ($unit_mit) | $meanintron |
| Total intron length ($unit_til) | $totalintron |
| Number of introns | $nintron |
| Number of exons per gene | $exonpergene |
| GC content | $gc |
| Repeat content | $repeattotal |
| \$nbsp; \$nbsp; \$nbsp; \$nbsp; Unclassified | $repeatunclass |
| BUSCO* genome score | $busco |

*BUSCO scores based on the $buscoset BUSCO set using $buscoversion C=complete  [S=single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison. 

EOF
}

tsvformat() {
cat << EOF

Genome assembly
Specimen	$species_name
Isolate	$isolate
Sequence coverage	$seq_coverage
Genome size ($unit_gs)	$genome_size
Gaps	$gaps
Number of contigs	$ncontigs
Contig N50 ($unit_n50)	$n50
Number of chromosomes	$nchrom
Number of protein-coding genes	$ngenes
Mean exon length (BP)	$meanexon
Total exon length ($unit_tel)	$totalexon
Number of exons	$nexon
Mean intron length (BP)	$meanintron
Total intron length ($unit_til)	$totalintron
Number of introns	$nintron
Number of exons per gene	$exonpergene
GC content	$gc
Repeat content	$repeattotal
Unclassified	$repeatunclass
BUSCO* genome score	$busco

*BUSCO scores based on the $buscoset BUSCO set using v5.5.0 C=complete  [S=single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison. 

EOF
}

markdownformat

tsvformat