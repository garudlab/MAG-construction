
The pipeline depends on singularity (I used 3.8.6).
Install the python requirements into a conda environment along with singularity

One of the singularity envs I built myself, I've attached containers/samtools_bwa.def to the repo. You can remote build the container `singularity build --remote containers/samtools_bwa.sif containers/samtools_bwa.def`, but you need to sign up for an account. You can also build it on your desktop. I can also transfer you the .sif (~250MB).

Create the folder `logs_cluster/`

I ran the pipeline as such
```snakemake -j 20 --use-singularity --singularity-args "--bind /u:/u" --keep-going --keep-incomplete --rerun-triggers mtime --singularity-prefix containers/```

Disclaimer that I actually used a conda env for DAS_tool, following bhattlab. During the code cleanup I changed it over to singularity for consistency and I tested that it worked, but if its giving you problems you can go back to conda. Uncomment the `#conda: "../envs/das_tool.yaml"` line, comment out the singularity container and add `--use-conda` to your snakemake command.

Set your file paths in `config.yaml`. Right now it needs a metadata table and it assumes co-assembly, although we can change that if need be.

This depends on 3 datafiles, a host genome, checkm, kraken. I have copies downloaded and can transfer them if need be.
