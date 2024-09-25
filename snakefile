######
# CrossM pipeline
# Antton Alberdi
# 2024/09/25
######

# If genomes are not split: seqkit split -i --id-regexp "^([^@]+)" --by-id-prefix "" --out-dir . drep.0.95.fa
# module purge && module load snakemake/7.20.0 mamba/1.3.1
# export XDG_CACHE_HOME=/maps/projects/mjolnir1/people/jpl786/.cache
# snakemake -j 20 --cluster 'sbatch -o log/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600


# List genome and target wildcards
genomes, = glob_wildcards("resources/genomes/{genome}.fa")

# Expand target files
rule all:
    input:
        expand("results/03_cross/{genome}.fq", genome=genomes),
        expand("results/03_cross/{genome}.bed", genome=genomes)

rule concatenate_fasta:
    input:
        expand("resources/genomes/{genome}.fa", genome=genomes)
    output:
        "results/01_genomes/allgenomes.fna"
    params:
        jobname="allgenomes.con"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        cat {input} > {output}
        """

rule kmer_db:
    input:
         "results/01_genomes/allgenomes.fna"
    output:
        "results/02_kmers/allgenomes.tsv"
    params:
        jobname="allgenomes.db",
        kmersize=21,
        minkmer=5
    threads:
        8
    resources:
        mem_gb=16,
        time=60
    conda:
        "workflow/envs/env.yml"
    shell:
        """
        jellyfish count -m {params.kmersize} -s 3300M -o {output} --out-counter-len 1 -L {params.minkmer} -t {threads} --text {input}
        """

rule browse_kmers:
    input:
        fasta="resources/genomes/{genome}.fa",
        kmerdb="results/02_kmers/allgenomes.tsv"
    output:
        fastq="results/03_cross/{genome}.fq",
        bed="results/03_cross/{genome}.bed"
    params:
        jobname="{genome}.kmer",
        kmersize=21
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    conda:
        "workflow/envs/env.yml"
    shell:
        """
        python workflow/scripts/kmer_browse.py \
                --fasta_file {input.fasta} \
                --kmer_count_file {input.kmerdb} \
                --kmer_size {params.kmersize} \
                --output_fastq {output.fastq} \
                --output_bed {output.bed}
        """
