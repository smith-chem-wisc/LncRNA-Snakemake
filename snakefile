# Samantha Shrum 2021.8.5
# (Several code contributions by A. Cesnik)

# Snakemake Pipeline for lncRNA analysis


GENOME_VERSION = "GRCh38"
ENSEMBL_VERSION = "96"
GENEMODEL_VERSION = GENOME_VERSION + "." + ENSEMBL_VERSION
GENOME_FA = f"lib/Genome/{GENOME_VERSION}.fa"
GENOME_GFF = f"lib/Genome/{GENEMODEL_VERSION}.gff3"
GENOME_GTF = f"lib/Genome/{GENEMODEL_VERSION}.gtf"

SLNCKY_ANNOTATIONS = "hg38"

REPS, THROWAWAY = glob_wildcards("{reps}//{rep}_R1.fastq")
SAMPLES = ["R1", "R2"]

rule all:
    input:
        "output/donefiles/folders_made.done",
        "output/IsoMat.txt"

rule make_folders:
    output:
        touch("output/donefiles/folders_made.done")

    shell:
        "(mkdir -p output/log && mkdir -p output/donefiles)"

rule trim_fastq:
    input:
        "output/donefiles/folders_made.done",
        i1 = "{rep}/{rep}_R1.fastq",
        i2 = "{rep}/{rep}_R2.fastq"

    output:
        o1 = "{rep}/{rep}_R1.trimmed.fastq",
        o2 = "{rep}/{rep}_R2.trimmed.fastq"

    log:        
        "output/log/1_trim_fastq_{rep}.log"

    shell:
        "fastp -q 20 -I {input.i1} -i {input.i2} -O {output.o1} -o {output.o2} -w 16 &> {log}"

rule genome_construction_1:
    input:
        "output/donefiles/folders_made.done"

    output:
        touch("output/donefiles/genome_construction_1.done")

    log:
        "output/log/2_genome_construction_1.log"

    shell:
        "STAR --runThreadN 23"
        " --runMode genomeGenerate"
        " --genomeFastaFiles {GENOME_FA}"
        " --sjdbGTFfile {GENOME_GTF} &> {log}"

rule rep_prep:
    input:
        "output/donefiles/genome_construction_1.done"
    output:
        touch("output/donefiles/{rep}_rep_prep.done")
    shell:
        "cp -r GenomeDir {wildcards.rep}/GenomeDir"

rule genome_alignment_1:
    input:
        "output/donefiles/{rep}_rep_prep.done",
        i1 = "{rep}/{rep}_R1.trimmed.fastq",
        i2 = "{rep}/{rep}_R2.trimmed.fastq"
        
    output:
        o = "{rep}/first.Aligned.sortedByCoord.out.bam"

    log:
        "output/log/3_{rep}_genome_alignment_1.log"

    shell:
        "STAR --runMode alignReads"
        " --runThreadN 10"
        " --readFilesIn {input.i1} {input.i2}"
        " --genomeDir ./{wildcards.rep}/GenomeDir"
        " --sjdbGTFfile {GENOME_GTF}"
        " --outSAMtype BAM SortedByCoordinate"
        " --outFilterIntronMotifs RemoveNoncanonicalUnannotated"
        " --outTmpDir {wildcards.rep}/tmp"
        " --outFileNamePrefix {wildcards.rep}/first. &> {log}"    


rule genome_construction_2:
    input:
        "{rep}/first.Aligned.sortedByCoord.out.bam"

    output:
        touch("{rep}/genome_construction_2.done")

    log:
        "output/log/4_{rep}_genome_construction_2.log"

    shell:
        "STAR --runMode genomeGenerate"
        " --genomeDir ./{wildcards.rep}/GenomeDir"
        " --genomeFastaFiles {GENOME_FA}"
        " --sjdbFileChrStartEnd ./GenomeDir/sjdbList.out.tab"
        " --sjdbGTFtagExonParentTranscript Parent"
        " --sjdbOverhang 100"
        " --runThreadN 10"
        " --limitSjdbInsertNsj 1200000"
        " --outTmpDir {wildcards.rep}/tmp &> {log}"

rule genome_alignment_2:
    input:
        "{rep}/genome_construction_2.done",
        i1 = "{rep}/{rep}_R1.trimmed.fastq",
        i2 = "{rep}/{rep}_R2.trimmed.fastq"
        
    output:
        "{rep}/second.Aligned.sortedByCoord.out.bam"

    log:
        "output/log/5_{rep}_genome_alignment_2.log"

    shell:
        "STAR --runMode alignReads"
        " --genomeDir ./{wildcards.rep}/GenomeDir"
        " --outSJfilterReads Unique"
        " --runThreadN 10"
        " --outSAMstrandField intronMotif"
        " --outFilterIntronMotifs RemoveNoncanonical"
        " --outSAMtype BAM SortedByCoordinate"
        " --outBAMcompression 10"
        " --outSAMattrRGline ID:1 PU:platform  PL:illumina SM:sample LB:library"
        " --outSAMmapqUnique 60"
        " --readFilesIn {input.i1} {input.i2}"
        " --outTmpDir {wildcards.rep}/tmp"
        " --outFileNamePrefix {wildcards.rep}/second. &> {log}"

rule merge_prep:
    input:
        "{rep}/second.Aligned.sortedByCoord.out.bam"
    output:
        "{rep}/second.Aligned.sortedByCoord.out.gtf"
    log:
        "output/log/6_{rep}_merge_prep.log"
    shell:
        "stringtie {input} -p 10 -G {GENOME_GFF} -o {output} -c 2.5 -m 300 -f .01 -A ./{wildcards.rep}/gene_abund.tab &> {log}"

rule merge:
    input:
        expand("{rep}/second.Aligned.sortedByCoord.out.gtf", rep = REPS)

    output:
        "output/merged.gtf"

    log:
        "output/log/7_merge.log"

    shell:
        "stringtie --merge -o {output} -c 2.5 -m 300 -f .01 -p 10 -i {input}"

# A. Cesnik
# https://github.com/acesnik/LncRNADiscovery/blob/master/rules/isoforms.smk
rule build_gtf_sharp:
    output:
        "lib/GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll"

    log:
        "output/log/8_build_gtf_sharp.log"

    shell:
        "(cd lib/GtfSharp && "
        "dotnet restore && "
        "dotnet build -c Release GtfSharp.sln) &> {log}"

# A. Cesnik
# https://github.com/acesnik/LncRNADiscovery/blob/master/rules/isoforms.smk
rule filter_transcripts_add_cds:
    input:
        gtfsharp="lib/GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll",
        gtf= "output/merged.gtf",

    output:
        "output/merged.filtered.withcds.gtf"

    log:
        "output/log/9_filter_transcripts_add_cds.log"

    shell:
        "dotnet {input.gtfsharp} -f {GENOME_FA} -g {input.gtf} -r {GENOME_GFF} &> {log}"

rule to_bed12:
    input:
        "output/merged.filtered.withcds.gtf"

    output:
        temp("temp.genePred"),
        temp("temp.bed"),
        "output/merged_trimmed.bed12"

    log:
        "output/log/10_to_bed12.log"

    shell:
        "(gtfToGenePred {input} temp.genePred && "
        "genePredToBed temp.genePred temp.bed && "
        "python3.8 lib/convert_ensembl2ucsc.py temp.bed output/merged.bed12 && "
        "python3.8 lib/modify_bed12.py) &> {log}"

# A. Cesnik
# https://github.com/acesnik/LncRNADiscovery/blob/master/rules/lncRNAs.smk
rule annotate_lncrnas:
    input:
        "lib/annotations/{SLNCKY_ANNOTATIONS}.fa",
        slncky = "lib/slncky/slncky.v1.0",
        bed12 = "output/merged_trimmed.bed12"

    output:
        "output/annotated.canonical_to_lncs.txt",
        "output/annotated.cluster_info.txt",
        "output/annotated.filtered_info.txt",
        "output/annotated.lncs.bed",
        "output/annotated.lncs.info.txt",
        "output/annotated.orfs.txt",
        "output/annotated.orthologs.top.txt",
        "output/annotated.orthologs.txt"

    params: 
        ref = SLNCKY_ANNOTATIONS

    threads:
        8

    conda:
        "lib/env.yml"

    log:
        "output/log/11_annotate_lncrnas.log"

    shell:
        "python {input.slncky} --threads {threads} {input.bed12} {params.ref} output/annotated &> {log}"


rule pre_RSEM_gtf_clean:
    input:
        "output/merged.gtf"

    output:
        "output/merged_strand_clean.gtf"

    log:
        "output/log/12_clean_strand.log"

    shell:
        "python3.8 lib/gtf_clean_strand.py &> {log}"

rule RSEM_create_reference:
    input:
        "output/merged_strand_clean.gtf"

    output:
        "output/rsemRef.transcripts.fa"

    log:
        "output/log/13_RSEM_create_reference.log"

    shell:
        "rsem-prepare-reference"
        " --num-threads 8"
        " --star"
        " --gtf {input} {GENOME_FA} {output} &> {log}"

rule RSEM_run:
    input:
        i1 = "{rep}/{rep}_R1.trimmed.fastq",
        i2 = "{rep}/{rep}_R2.trimmed.fastq",
        i3 = "output/rsemRef.transcripts.fa"

    output:
        "{rep}/rsemRefOut.isoforms.results"

    log:
        "output/log/14_{rep}_RSEM_run.log"

    shell:
        "rsem-calculate-expression"
        " --no-bam-output"
        " --time"
        " --star"
        " --num-threads 8"
        " --paired-end {input.i1} {input.i2} output/rsemRef {output} &> {log}"

rule RSEM_DEA:
    input:
        expand("{rep}/rsemRefOut.isoforms.results", rep = REPS)

    output:
        "output/IsoMat.txt",

    log:
        "output/log/15_RSEM_DEA.log"

    shell:
        "rsem-generate-data-matrix {input} > {output}"