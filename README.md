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


The following scripts were used for PhytroExpr processing.

PhytoExpr.scores.R: This script is used to process predictions in the file bdi.newexppred.sequences.csv and outputs scores.tsv, which contains the scores for each variant, calculated as the score of the alt variant subtracted by the score of the ref variant. Additionally the file scores.ref.alt.tsv is output which contains the ref and alt variant scores by themselves.

The following scripts were used for PhytroExpr processing.

PhytoExpr.scores.R: This script is used to process predictions in the file bdi.newexppred.sequences.csv and outputs scores.tsv, which contains the scores for each variant, calculated as the score of the alt variant subtracted by the score of the ref variant. Additionally the file scores.ref.alt.tsv is output which contains the ref and alt variant scores by themselves.

The following script was used for generating REF and ALT a2z scores for BD21.3 variants.

a2z.ocr.py: Uses provided model (model-accessibility-full.h5) to generate a2z scores for REF and ALT alleles of all variants in provided .bim file (snps.combined.bim). Sequences are obtained from provided fasta file (BdistachyonBd21_3_537_v1.0.fa) to obtain sequences. It uses functions from a2z source, which is expected to be located in the folder (a2z). A zip containing the files is present in this github (a2z.zip) .

The following scripts were used for generating caduceus and a2z embeddings for BD21.3 variants.

FASTA.py, bedbug.py, mRNA.py: Tools needed by make.bd.data.py

make.bd.data.py: This script generates unique TTS and TSS sequences containing variants for each gene, for each plant. It needs a fasta file defining sequences (BdistachyonBd21_3_537_v1.0.fa), gff defining genes (annotation/BdistachyonBd21_3_537_v1.2.gene.gff3), bed/bim/fam files defining all variants (snps.combined.*), a file (gene.id.translation.tsv) providing gene id translations between BD21 and BD21.3. The script generates sequences in the file (data.bd.csv) which generally imitates the format of the file (data.csv) which contains the published Phytoexpr training data, and stats into the file (stats.tsv) counting number of unique sequences per transcript.
make.bd.sequences.py: This scripts takes the output from make.bd.data.py (data.bd.tsv) and prepares the sequences for a2z embedding extraction in the file (bd.sequences.a2z.tsv), and for caduceus embedding extraction in the file (bd.sequences.caduceus.tsv). 
make.bd.embedding.a2z.py: This script takes the prepared sequences in the file (bd.sequences.a2z.tsv), extracts a2z embeddings which are output into (embeddings.bd.a2z.h5).
make.bd.embedding.caduceus.py: This script takes the prepared sequences in the file (bd.sequences.cadeuceus.tsv), extracts a2z embeddings which are output into (embeddings.bd.cadeuceus.h5).
