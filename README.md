All the files assume a certain folder structure which will need to be updated if used in any other settings.

These following were scripts used to obtain TSS and TTS features from PlantCaduceus and a2z, upon which the models were trained. The DNA sequences and TPM data are sourced from the PlantCaduceus training data.
Unlike PlantCaduceus we only used about a 5000bp TSS/TTS sequence length.

prepare.sequences.py: This script reads PlantCaduceus training data file "data.csv", and outputs sequences to be used for PlantCaduceus and a2z feature extraction.

sequence2embedding.caduceus.py: This script extracts features from PlantCaduceus for the sequences. It requires a CUDA capable machine.

sequence2embedding.a2z.py: This script extracts features from a2z for the sequences. It runs on CPU.

The following scripts were used for RNAseq processing.

RNAseq.pipeline.sh: This shell script is used to process raw RNAseq reads
RNAseq.R: This shell script is used to calculate TPMs for each gene in each sample using the output from "RNAseq.pipeline.sh"
RNAseq.QC.py: This script parses QC from the output of "RNAseq.pipeline.sh", to determine which lines fail the minimum 70% aligned reads. 
