These are the scripts used for obtaining the TSS and TTS features from PlantCaduceus and a2z, upon which the models were trained. The DNA sequences and TPM data are sourced from the PlantCaduceus training data.
Unlike PlantCaduceus we only used about a 5000bp TSS/TTS sequence length.

prepare.sequences.py: This script reads PlantCaduceus training data file "data.csv", and outputs sequences to be used for PlantCaduceus and a2z feature extraction.

sequence2embedding.caduceus.py: This script extracts features from PlantCaduceus for the sequences. It requires a CUDA capable machine.

sequence2embedding.a2z.py: This script extracts features from a2z for the sequences. It runs on CPU.
