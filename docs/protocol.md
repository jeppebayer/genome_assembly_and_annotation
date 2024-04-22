# PROTOCOL

- [PROTOCOL](#protocol)
	- [Draft assembly](#draft-assembly)
	- [Hi-C scaffolding](#hi-c-scaffolding)
	- [Annotation and masking of repetitive content](#annotation-and-masking-of-repetitive-content)
	- [Gene annotation](#gene-annotation)
	- [Software](#software)
		- [Draft assembly](#draft-assembly-1)
	- [References](#references)
	- [Acknowledgements](#acknowledgements)

## Draft assembly

The workflow template file is located [here](../scripts/draft_assembly/workflow_source/workflow_templates.py)

The workflow source file is located [here](../scripts/draft_assembly/workflow_source/workflow_source.py)

Prior to draft assembly of PacBio HiFi sequence data, `HiFiAdapterRemoval` is used to remove any potetial remaining adapter sequences. For the draft assembly `hifiasm` is used for which there are several options:

1. Create a haplotype phased assembly. This is feasable if the HiFi data is from a single individual. Two separate files will be created, one containing each haplotype.
2. Create a primary assembly. When the HiFi data is pooled from multiple individuals this is the proper choice. A primary and a secondary sequence file will be created. The primary sequence is the one to continue work on.
3. Create either a haplotype phased assembly or a primary assembly supported with Hi-C sequence data. Theoretically, if Hi-C data is available, this should create the better assembly, testing on pooled data for a collembola species showed no improvement.

## Hi-C scaffolding

## Annotation and masking of repetitive content

## Gene annotation

## Software

### Draft assembly

- HiFiAdapterFilt v2.0.1
- hifiasm v0.19.8
- BUSCO v5.5.0
- purge_dups v1.2.6
- minimap2 v2.26

## References

## Acknowledgements