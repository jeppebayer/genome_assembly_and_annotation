#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../../../workflow_source/'))
from workflow_source import *

gwf = braker3_annotation_workflow()

# cd /home/anneaa/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/scripts/annotation/configurations/Tethyida/Tethya_wilhelma
# conda activate annotation
# gwf status -f summary
