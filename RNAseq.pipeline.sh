#!/bin/sh
#List of lines to process
LINES='100 130 156 181 26 53 802 827 852 102 131 157 182 27 54 803 828 853 103 132 158 183 28 55 804 829 854 104 133 159 185 2 56 805 82 855 105 134 15 186 30 57 806 830 856 106 135 160 187 31 58 807 831 857 107 136 162 188 32 59 808 832 85 108 137 163 189 33 5 809 833 86 109 138 164 18 34 60 80 834 87 10 139 165 190 35 61 810 835 88 110 13 166 191 37 62 811 836 89 111 140 167 192 38 63 812 837 8 112 141 168 193 39 64 813 839 90 113 142 169 194 3 65 814 83 91 114 143 16 196 40 67 815 840 92 117 144 170 197 41 69 816 841 94 118 145 171 198 42 6 817 842 95 119 146 172 199 44 70 818 843 96 11 147 173 19 45 71 819 844 97 121 148 174 1 46 72 81 845 98 122 149 175 200 47 73 820 846 99 123 150 176 20 48 74 821 847 9 124 151 177 21 49 76 822 848 126 152 178 22 4 77 823 849 128 153 179 23 50 79 824 84 129 154 17 24 51 7 825 850 12 155 180 25 52 801 826 851'
#The folder containing the tool binaries used for processing
TFOLDER=/home/camous/sieve/Camous/TOOLS
#The input folder that contains 00_fastq, where the RNAseqs are read from
IFOLDER=/home/camous/sieve/90-1080311147/00_fastq
#The putput folder that contains 01_output into which files are written
OFOLDER=/home/camous/sieve/Camous/90-1080311147/01_output

$TFOLDER/cufflinks/gffread /Volumes/N1/INPUT/GENOME/BD/REF/annotation/BdistachyonBd21_3_537_v1.2.gene.gff3 -T -o $OFOLDER/BD.gtf
python3 $TFOLDER/hisat2/hisat2_extract_splice_sites.py $OFOLDER/BD.gtf > $OFOLDER/BD.ss
python3 $TFOLDER/hisat2/hisat2_extract_exons.py $OFOLDER/BD.gtf > $OFOLDER/BD.exon
$TFOLDER/hisat2/hisat2-build -p 8 --exon $OFOLDER/BD.exon --ss $OFOLDER/BD.ss /Volumes/N1/INPUT/GENOME/BD/REF/assembly/BdistachyonBd21_3_537_v1.0.fa $OFOLDER/BD
for i in $LINES
do
	fastp -w 8 -i $IFOLDER'/'$i'A1_R1_001.fastq.gz' -I $IFOLDER'/'$i'A1_R2_001.fastq.gz' -o $OFOLDER'/'$i'_R1.fq.gz' -O $OFOLDER'/'$i'_R2.fq.gz'  -j $OFOLDER'/'$i'_fastp.json' -h $OFOLDER'/'$i'_fastp.html'
    $TFOLDER/hisat2/hisat2 -p 8 -x $OFOLDER/BD -1 $OFOLDER'/'$i'_R1.fq.gz' -2 $OFOLDER'/'$i'_R2.fq.gz' | samtools view -bSh > $OFOLDER'/'$i'_aligned.bam'
	$TFOLDER/subread/featureCounts -T 8 -p -a $OFOLDER'/BD.gtf' -o $OFOLDER'/'$i'_aligned.FC' $OFOLDER'/'$i'_aligned.bam'
done
