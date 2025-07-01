#!/bin/env python3
import sys, os, argparse, pandas as pd

def parse_args(args = None):
	description = 'Combine files to create input file for blobtoolkit window-stats.'

	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--freq', help='Frequence file produced by fasta_windows', required=True)
	parser.add_argument('--mononuc', help='Mononucleotide file produced by fasta_windows', required=True)
	parser.add_argument('--depth', help='Depth file produced by blotk depth', required=True)
	parser.add_argument('--countbuscos', help='Count of BUSCO genes within regions produced by blobtoolkit count-busco-genes', required=True)
	parser.add_argument('--prefix', help='Prefix for output TSV', required=True)
	return parser.parse_args(args)

def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)

def merge_files(freqFile: str, mononucFile: str, depthFile: str, countBuscoFile: str):
	freqFileIn = pd.read_csv(freqFile, sep='\t').rename(columns={"ID": "sequence", "GC_prop": "gc", "Prop_Ns": "n", "N": "ncount"})
	mononucFileIn = pd.read_csv(mononucFile, sep='\t').rename(columns={"ID": "sequence", "GC_prop": "gc", "Prop_Ns": "n", "N": "ncount"})
	depthFileIn = pd.read_csv(depthFile, sep='\t', names=['sequence', 'start', 'end', f'{os.path.basename(depthFile).replace("regionsCoverage.bed", "")}_cov'])
	countBuscoFileIn = pd.read_csv(countBuscoFile, sep='\t').rename(columns={'ID': 'sequence'})
	combined = freqFileIn.merge(mononucFileIn).merge(countBuscoFileIn).merge(depthFileIn)
	return combined

def main(args = None):
	args = parse_args(args)
	make_dir(os.path.dirname(args.prefix))
	merge_files(args.freq, args.mononuc, args.depth, args.countbuscos).to_csv(f'{args.prefix}.tsv', sep='\t', index=False)

if __name__ == '__main__':
	sys.exit(main())