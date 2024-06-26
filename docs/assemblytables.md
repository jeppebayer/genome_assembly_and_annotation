# Genome Assembly Tables

## Population genetics

### *Entomobrya nicoleti*

Assembly file: [EntNic.softmasked.fna](../../genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/repeatmasking/EntNic.softmasked.fna)

| Genome Assembly |  |
| --- | :---: |
| Specimen | *Entomobrya nicoleti* |
| Isolate | Pooled 20 individuals (males + females) |
| Sequence coverage | 69 |
| Genome size (Mb)* | 254.7 |
| Gaps | 0.604% |
| Number of contigs | 3,196 |
| Contig N50 (Kb) | 150 |
| Number of chromosomes | 7 |
| Number of protein-coding genes | 31,573 |
| Mean exon length (Bp) | 367 |
| Total exon length (Mp) | 40.2 (15.79%) |
| Number of exons | 109,588 |
| Mean intron length (Bp) | 409 |
| Total intron length (Mb) | 31.9 (12.53%) |
| Number of introns | 78,015 |
| Number of exons per gene | 3.5 |
| GC content | 34.27% |
| Repeat content | 20.77% |
| &nbsp; &nbsp; &nbsp; &nbsp; Unclassified | 14.62% |
| BUSCO** genome score*** | C:50.9%[S:48.9%,D:2.0%],F:2.4%,M:46.7%,n:1013 |

**Assembled* genome size. Total genome size including contigs which could not be Hi-C scaffolded is 367.7 Mb.  
**BUSCO scores based on the arthropoda_odb10 BUSCO set using v5.5.0 C=complete  [S=single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison.  
***BUSCO score for the assembled genome. The BUSCO score for the total genome including contigs which could not be Hi-C scaffolded: C:93.9%[S:90.5%,D:3.4%],F:1.8%,M:4.3%,n:1013.

Files:

- [BUSCO table no debris](../../genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/HiC_scaffolding/juicer/3ddna_in10000_r6/finalize/busco_purged.nodebris.fasta/short_summary.specific.arthropoda_odb10.busco_purged.nodebris.fasta.txt)
- [BUSCO table total genome](../../genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/HiC_scaffolding/juicer/3ddna_in10000_r6/finalize/busco_purged_HiC.fasta/short_summary.specific.arthropoda_odb10.busco_purged_HiC.fasta.txt)
- [RepeatMasker](../../genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/repeatmasking/RepeatMasker/EntNic.fasta.fullmask.tbl)
- [Annotation](../../genome_assembly_and_annotation/steps/collembola/Entomobrya_nicoleti/annotation/braker2/braker.gtf)
  - Number of genes
  
    ```awk
    awk 'BEGIN{FS = OFS = "\t"} {if ($3 == "gene") {sum += 1}} END{print sum}' annotation.gtf
    ```

  - Total exon/intron size

    ```awk
    awk 'BEGIN{FS = OFS = "\t"} {if ($9 ~ /.t1/) {if ($3 == "exon") {sum += $5 - $4 + 1}}} END{print sum}' annotation.gtf
    ```

  - Total number of exons/introns

    ```awk
    awk 'BEGIN{FS = OFS = "\t"} {if ($9 ~ /.t1/) {if ($3 == "exon") {sum += 1}}} END{print sum}' annotation.gtf
    ```

- [PacBio HiFi data](../../old_setup/data/EntNic_HiFi/mixed99/EntNic_from_mixed_sample_99.fq)

  - Sequence coverage

    ```awk
    seqtk size HiFidata.fq | awk '{FS = OFS = "\t"} {print $2 / 367698499}'
    ```

### *Isotoma viridis*

### *Orchesella villosa*

### *Pogonognathellus flavescens*

### *Apion virens*

### *Bathyphantes gracilis*

## Museomics

### *Gonepteryx rhamni*

### *Lycaena_virgaureae*

### *Thymelicus lineola*