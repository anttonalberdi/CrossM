######
# CrossM pipeline
# Antton Alberdi
# 2024/09/25
######

# If genomes are not split: seqkit split -i --id-regexp "^([^@]+)" --by-id-prefix "" --out-dir . drep.0.95.fa
# module purge && module load snakemake/7.20.0 mamba/1.3.1
# snakemake -j 20 --cluster 'sbatch -o logs/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600


# List genome and target wildcards
genomes, = glob_wildcards("resources/genomes/{genome}.fa")

# Expand target files
rule all:
    input:
        "results/02_kmers/allgenomes.tsv"
        #expand("results/probes/{target}.tsv", target=targets)

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
        python workflow/scripts/concatenate_fasta.py {output} {input}
        """

# Jellyfish count renamed fasta file for downstream kmer count
rule build_jellyfish:
    input:
         "results/01_genomes/allgenomes.fna"
    output:
        "results/02_kmers/allgenomes.tsv"
    params:
        jobname="allgenomes.jf",
        kmersize=21,
        minkmer=5
    resources:
        mem_gb=16,
        time=60
    shell:
        """
        module load jellyfish/2.2.10
        jellyfish count -m {params.kmersize} -s 3300M -o {output} --out-counter-len 1 -L {params.minkmer} --text {input}
        """
