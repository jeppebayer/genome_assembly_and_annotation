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

| Genome Assembly |  |
| --- | :---: |
| Specimen | *Orchesella villosa* |
| Isolate |  |
| Sequence coverage |  |
| Genome size (MB) | 295 |
| Gaps | 0.172% |
| Number of contigs | 1078 |
| Contig N50 (KB) | 649 |
| Number of chromosomes | 6 |
| Number of protein-coding genes | 22713 |
| Mean exon length (BP) | 304 |
| Total exon length (MB) | 30.6 (10.38%) |
| Number of exons |  100827 |
| Mean intron length (BP) | 438 |
| Total intron length (MB) | 34.2 (11.61%) |
| Number of introns | 78114 |
| Number of exons per gene | 4.4 |
| GC content |  |
| Repeat content | 20.89% |
| &nbsp; &nbsp; &nbsp; &nbsp; Unclassified | 15.38% |
| BUSCO* genome score | C:94.5%[S:89.5%,D:5.0%],F:2.2%,M:3.3%,n:1013 |

*BUSCO scores based on the arthropoda_odb10 BUSCO set using v5.5.0 C=complete  [S=single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison.

RepeatMasker results: /faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Orchesella_villosa/repeatmasking/RepeatMasker/OrcVil.fasta.fullmask.tbl

BUSCO result: /faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Orchesella_villosa/HiC_scaffolding/juicer/3ddna_in10000_r6/finalize/busco_OrcVil.purged.nodebris.fasta/short_summary.specific.arthropoda_odb10.busco_OrcVil.purged.nodebris.fasta.txt

Annotation result: /faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Orchesella_villosa/annotation/braker3/braker.gtf

### *Pogonognathellus flavescens*

| Genome Assembly |  |
| --- | :---: |
| Specimen | *Pogonognathellus flavescens* |
| Isolate |  |
| Sequence coverage |  |
| Genome size (MB) | 327.4 |
| Gaps | 0.191% |
| Number of contigs | 1311 |
| Contig N50 (KB) | 544 |
| Number of chromosomes | 8 |
| Number of protein-coding genes | 13989 |
| Mean exon length (BP) | 230 |
| Total exon length (MB) | 17.3 (5.28%) |
| Number of exons |  75122 |
| Mean intron length (BP) | 373 |
| Total intron length (MB) | 22.8 (6.97%) |
| Number of introns | 61133 |
| Number of exons per gene | 5.4 |
| GC content |  |
| Repeat content | 21.23% |
| &nbsp; &nbsp; &nbsp; &nbsp; Unclassified | 16.60% |
| BUSCO* genome score | C:94.5%[S:91.4%,D:3.1%],F:1.5%,M:4.0%,n:1013 |

*BUSCO scores based on the arthropoda_odb10 BUSCO set using v5.5.0 C=complete  [S=single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison.

RepeatMasker results: /faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Pogonognathellus_flavescens/repeatmasking/RepeatMasker/PogFla.fasta.fullmask.tbl

BUSCO result: /faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Pogonognathellus_flavescens/HiC_scaffolding/juicer/3ddna_in10000_r6/finalize/busco_PogFla.purged.nodebris.fasta/short_summary.specific.arthropoda_odb10.busco_PogFla.purged.nodebris.fasta.txt

Annotation result: /faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/collembola/Pogonognathellus_flavescens/annotation/braker3/braker.gtf

### *Apion virens*

### *Bathyphantes gracilis*

## Museomics

### *Gonepteryx rhamni*

### *Lycaena_virgaureae*

### *Thymelicus lineola*
