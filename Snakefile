# Snakemake pipeline definition file

# metadata
__author__ = "Thom Griffioen"
__copyright__ = "Copyright 2019 Thom Griffioen"
__email__ = "t.griffioen@nioo.knaw.nl"
__license__ = "MIT"

from snakemake.utils import min_version

min_version("5.2.0")

# configuration
configfile: "config.yml"  # global configuration
configfile: "samples_f8.yml"  # input sample files
localrules: all, multiqc, checkm_mkroot


# project constants
PROJECT = config["project"]
DIRECTIONS = config["directions"]
SAMPLES = config["data"]

# output files
OUTFILES = list()
if config["parts-enabled"]["qc"]:
    OUTFILES.append("{project}/multiqc/qc_report.html")
if config["parts-enabled"]["krona"]:
    OUTFILES.append("{project}/krona/{sample}.html")
if config["parts-enabled"]["diamond"]:
    OUTFILES.append("{project}/diamond/{sample}.daa")
if config["parts-enabled"]["quast"]:
    OUTFILES.append("{project}/quast/")
if config["parts-enabled"]["bwa"]:
    OUTFILES.append("{project}/samtools/{sample}/alignment.bam")
if config["parts-enabled"]["checkm"]:
    OUTFILES.append("{project}/checkm/{sample}/checkm_lineagewf.tsv")
if config["parts-enabled"]["prodigal"]:
    OUTFILES.append("{project}/prodigal/{sample}.gff")
if config["parts-enabled"]["mmseqs"]:
    OUTFILES.append("{project}/mmseqs/{sample}/results.tsv")
OUTFILES.append("{project}/mg_pheno/{sample}/alignment.bam")

# rule to target all output files
rule all:
    input:
        expand(OUTFILES, project=PROJECT, sample=SAMPLES, direction=DIRECTIONS)


#
# QUALITY CONTROL
#
rule fastqc:
    input:
        lambda wildcards: config["data"][wildcards.sample][wildcards.direction]
    output:
# apparently *has* to end with '_fastqc.zip' for multiqc to detect it ...
# https://github.com/ewels/MultiQC/blob/v1.7/multiqc/modules/fastqc/fastqc.py#L49
        "{project}/fastqc/{sample}_{direction}_fastqc.zip"
    conda:
        "py3_env.yml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 300
    params:
        out_dir = "{project}/fastqc/"
    shell:
        'TMPDIR="$(mktemp -d)" '
        '&& fastqc --noextract -q -f fastq -t "{threads}" -o "${{TMPDIR}}" "{input}" '
        '&& mv -f "${{TMPDIR}}"/*.zip "{output}" '
        '&& rm -rf "${{TMPDIR}}"'


rule multiqc:
    input:
        expand("{{project}}/fastqc/{sample}_{direction}_fastqc.zip", sample=SAMPLES, direction=DIRECTIONS)
    output:
        html = "{project}/multiqc/qc_report.html",
        zip = "{project}/multiqc/qc_report_data.zip"
    conda:
        "py3_env.yml"
    threads: 1
    resources:
        mem_mb = 1000
    params:
        fastqc_results = "{project}/fastqc/",
        outdir = "{project}/multiqc/",
        html = "qc_report.html"
    shell:
        'multiqc -n "{params.html}" -o "{params.outdir}" -z "{params.fastqc_results}"'


#
# TAXONOMIC CLASSIFICATION
#
rule kaiju:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]["R1"],
        reverse = lambda wildcards: config["data"][wildcards.sample]["R2"]
    output:
        "{project}/kaiju/{sample}.tsv"
    conda:
        "py3_env.yml"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 80000 + attempt * 20000
    params:
        kaiju_files = '-t "{0}/nodes.dmp" -f "{0}/kaiju_db_nr_euk.fmi"'.format(config["static-files"]["kaiju-db"]),
        kaiju_args = '-a greedy -e 3'
    shell:
        'kaiju -z {threads} {params.kaiju_args} {params.kaiju_files} -i "{input.forward}" -j "{input.reverse}" -v -o "{output}"'


rule krona:
    input:
        rules.kaiju.output
    output:
        html = "{project}/krona/{sample}.html",
        krona = "{project}/krona/{sample}.krona"
    conda:
        "py3_env.yml"
    threads: 1
    resources:
        mem_mb = 5000
    params:
        kaiju_files = '-t "{0}/nodes.dmp" -n "{0}/names.dmp"'.format(config["static-files"]["kaiju-db"])
    shell:
        'kaiju2krona {params.kaiju_files} -i "{input}" -o "{output.krona}" '
        '&& ktImportText -o "{output.html}" "{output.krona}"'


#
# READ ALIGNMENT
#
rule diamond:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]["R1"],
        reverse = lambda wildcards: config["data"][wildcards.sample]["R2"]
    output:
        "{project}/diamond/{sample}.daa"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 20000 + attempt * 20000
    conda:
        "py3_env.yml"
    params:
        diamond_args = "--evalue 0.001",
        db = config["static-files"]["diamond-db"]
    shell:
        'zcat "{input.forward}" "{input.reverse}" '
        '| diamond blastx {params.diamond_args} --threads {threads} --block-size 5.0 --db "{params.db}" --outfmt 100 --out "{output}"'


#
# ASSEMBLY
#
rule spades:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]["R1"],
        reverse = lambda wildcards: config["data"][wildcards.sample]["R2"]
    output:
        contigs = "{project}/spades/{sample}/contigs.fasta",
        scaffolds = "{project}/spades/{sample}/scaffolds.fasta"
    threads: 40
    resources:
        mem_mb = 800000
    conda:
        "py3_env.yml"
    params:
        out_dir = '{project}/spades/{sample}/'
    shell:
        'metaspades.py -k 21,33,55 -1 "{input.forward}" -2 "{input.reverse}" --threads {threads} --memory 800 --only-assembler -o "{params.out_dir}"'


rule quast:
    input:
        expand("{{project}}/spades/{sample}/contigs.fasta", sample=SAMPLES)
    output:
        directory("{project}/quast/")
    threads: 8
    resources:
        mem_mb = 5000
    conda:
        "py3_env.yml"
    params:
        sample_names = ','.join(SAMPLES)
    shell:
        'metaquast.py --min-contig 0 --max-ref-number 0 --no-icarus --threads {threads} --labels "{params.sample_names}" --output-dir "{output}" {input}'


#
# BINNING
#
rule depth:
    input:
        assembly = rules.spades.output.contigs,
        forward = lambda wildcards: config["data"][wildcards.sample]["R1"],
        reverse = lambda wildcards: config["data"][wildcards.sample]["R2"]
    output:
        fasta = "{project}/samtools/{sample}/filtered.fa.gz",
        bam = "{project}/samtools/{sample}/alignment.bam",
        depth = "{project}/samtools/{sample}/depth.txt"
    threads: 16
    resources:
        mem_mb = lambda wildcards, threads: max(threads * 2000, 5000)
    conda:
        "py3_env.yml"
    params:
        min_length = config["min-contig-length"],
        index = "{project}/samtools/{sample}/index",
        min_mq = 20
    shell:
        'bioawk -c fastx \'length($seq)<{params.min_length} {{exit 0;}} {{print ">" $name,$comment "\\n" $seq;}}\' "{input.assembly}" '
        '| pigz --processes {threads} --stdout > "{output.fasta}" '
        '&& bwa index -p "{params.index}" "{output.fasta}" '
        '&& bwa mem -t {threads} "{params.index}" "{input.forward}" "{input.reverse}" '
        '| samtools view --threads 1 -b - '
        '| samtools sort --threads {threads} --output-fmt bam -o "{output.bam}" - '
        '&& samtools index -b -@ {threads} "{output.bam}" '
        '&& jgi_summarize_bam_contig_depths --minMapQual {params.min_mq} --referenceFasta "{output.fasta}" --outputDepth "{output.depth}" "{output.bam}"'


rule metabat:
    input:
        assembly = rules.spades.output.contigs,
        depth = rules.depth.output.depth
    output:
        directory("{project}/metabat/{sample}/")
    threads: 8
    resources:
        mem_mb = 5000
    conda:
        "py3_env.yml"
    params:
        out_dir = "{project}/metabat/{sample}/{sample}",
        min_length = rules.depth.params.min_length
    shell:
        'metabat2 --minContig {params.min_length} --inFile "{input.assembly}" --abdFile "{input.depth}" --outFile "{params.out_dir}" --numThreads {threads}'


rule checkm_mkroot:
    output:
        directory("{project}/checkm_root/")
    threads: 1
    resources:
        mem_mb = 1000
    conda:
        "py2_env.yml"
    params:
        dat_file = '{project}/checkm_root/checkm_db.tgz'
    shell:
        'wget --quiet --output-document "{params.dat_file}" https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz '
        '&& tar --directory "{output}" --extract --file "{params.dat_file}" '
        '&& checkm data setRoot "{output}"'


rule checkm:
    input:
        metabat = rules.metabat.output,
        checkm_root = rules.checkm_mkroot.output
    output:
        bin_stats = "{project}/checkm/{sample}/checkm_lineagewf.tsv",
        bin_dir = directory("{project}/checkm/{sample}/bins/"),
        lineage_markers = "{project}/checkm/{sample}/lineage.ms"
    threads: 8
    resources:
        mem_mb = 50000
    conda:
        "py2_env.yml"
    params:
        out_dir = "{project}/checkm/{sample}/"
    shell:
        'checkm lineage_wf --force_overwrite --tab_table --file "{output.bin_stats}" --threads {threads} --extension "fa" "{input.metabat}" "{params.out_dir}"'


#
# CONTIG ALIGNMENT
#
rule prodigal:
    input:
        rules.spades.output.contigs
    output:
        nucl = "{project}/prodigal/{sample}.fna",
        prot = "{project}/prodigal/{sample}.faa",
        gff = "{project}/prodigal/{sample}.gff"
    threads: 1
    resources:
        mem_mb = 500
    conda:
        "py3_env.yml"
    params:
        min_length = config["min-contig-length"]
    shell:
        'bioawk -c fastx \'length($seq)<{params.min_length} {{exit 0;}} {{print ">" $name,$comment "\\n" $seq;}}\' "{input}" '
        '| prodigal -q -f gff -p meta -d "{output.nucl}" -a "{output.prot}" -o "{output.gff}"'


rule make_ril_phenotype:
    input:
        fwd = lambda wildcards: config["data"][wildcards.sample]["R1"],
        rev = lambda wildcards: config["data"][wildcards.sample]["R2"]
    output:
        bam = "{project}/mg_pheno/{sample}/alignment.bam",
        bed = "{project}/mg_pheno/{sample}/depth.bed",
        gene_cov = "{project}/mg_pheno/{sample}/coverage.csv",
    threads: 16
    resources:
        mem_mb = lambda wildcards, threads: max(threads * 2000, 5000)
    conda:
        "py3_env.yml"
    params:
        index = "results_p/prodigal/db/all_bwa"
    shell:
        'bwa mem -t {threads} "{params.index}" "{input.fwd}" "{input.rev}" '
        '| samtools view --threads 1 -b - '
        '| samtools sort --threads {threads} --output-fmt bam -o "{output.bam}" - '
        '&& samtools index -b -@ {threads} "{output.bam}" '
        '&& bedtools genomecov -bga -ibam "{output.bam}" > "{output.bed}" '
        '&& python3 bed_stats.py "{output.bed}" > "{output.gene_cov}"'

rule mmseqs:
    input:
        rules.prodigal.output.prot
    output:
        "{project}/mmseqs/{sample}/results.tsv"
    threads: 32
    resources:
        mem_mb = 100000
    conda:
        "py3_env.yml"
    params:
        db = config["static-files"]["mmseqs-db"],
        tmp_dir = "/tmp",
        out_fmt = "query,target,pident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits,qcov,tcov"
        # Query name
        # Target name
        # Percentage identity
        # Alignment length
        # Mismatches
        # Query start position
        # Query end position
        # Target start position
        # Target end position
        # Expect value
        # Bit score
        # Query coverage
        # Target coverage
    shell:
        'mmseqs easy-search --threads {threads} --search-type 1 --format-mode 0 --format-output "{params.out_fmt}" --sort-results 1 -v 2 "{input}" "{params.db}" "{output}" "{params.tmp_dir}"'
        # 'mmseqs createdb --dbtype 1 -v 2 --compressed 1 "{input}" "{output.query_db}" '
#        '&& mmseqs createindex --threads {threads} --search-type 1 --compressed 1 "{output.query_db}" "{params.tmp_dir}" '
#         '&& mmseqs search --threads {threads} --search-type 1 --sort-results 1 --compressed 1 -v 2 "{output.query_db}" "{params.db}" "{output.results_db}" "{params.tmp_dir}" '
#         '&& mmseqs convertalis --threads {threads} --format-mode 0 --format-output "{params.out_fmt}" -v 2 "{output.query_db}" "{params.db}" "{output.results_db}" "{output.alignment}"'
