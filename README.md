# Metagenome-QTL

Create the conda environment needed for running `snakemake` and [`gen_conf.py`](gen_conf.py).
```sh
$ conda create --file environment.yml --name YOUR_ENV_NAME
$ conda activate YOUR_ENV_NAME
```

Generate a YAML file containing metadata concerning your shotgun sequence files.
```sh
$ python3 gen_conf.py PATH_TO_RAW_FASTQ/ > samples.yml
```

Edit the options in [`config.yml`](config.yml) to fit your analysis.
The `project` field is the path where your results will be written to.
Options under the `parts-enabled` field control which parts of the pipeline to run.

Run the pipeline using the [`run_pipeline.sh`](run_pipeline.sh) wrapper script.
```sh
$ bash run_pipeline.sh [SNAKEMAKE_ARGS]
```
This script will automatically run the pipeline using SLURM when executed on the `nioo0004.nioo.int` master node.
You can change the contents of [`slurm.yml`](slurm.yml) to modify the behaviour of the scheduler if needed.
