
The pipeline depends on singularity (I used 3.8.6).
Install the python requirements into a conda environment along with singularity

One of the singularity envs I built myself, I've attached containers/samtools_bwa.def to the repo. You can remote build the container `singularity build --remote containers/samtools_bwa.sif containers/samtools_bwa.def`, but you need to sign up for an account. You can also build it on your desktop. I can also transfer you the .sif (~250MB).

Create the folder `logs_cluster/`

I ran the pipeline as such
```snakemake -j 20 --use-singularity --singularity-args "--bind /u:/u" --keep-going --keep-incomplete --rerun-triggers mtime --singularity-prefix containers/```

Disclaimer that I actually used a conda env for DAS_tool, following bhattlab. During the code cleanup I changed it over to singularity for consistency and I tested that it worked, but if its giving you problems you can go back to conda. Uncomment the `#conda: "../envs/das_tool.yaml"` line, comment out the singularity container and add `--use-conda` to your snakemake command.

Set your file paths in `config.yaml`. Right now it needs a metadata table and it assumes co-assembly, although we can change that if need be.

This depends on 3 datafiles, a host genome, checkm, kraken. I have copies downloaded and can transfer them if need be.


### Issues when running on a large set of genomes
I encountered several technical issues when trying to assemble MAGs from a set of >1800 metagenomic samples. Snakemake runs into serious problems when the DAG construction is too large, and its scheduler can throw errors, meaning that you need to babysit the pipeline and repeatedly restart it, which is undesireable. The following helped:

* Run snakemake on a compute node. While snakemake claims to be meant to run on head node, frequently there are thread assignment python errors which crash the pipeline. It's much better to login to a compute node. To keep the window open I was using tmux. There was an occasional issue where I would be mysteriously logged out of the compute node on my tmux window, and I was unable to resolve this. I should contact IT if the problem occurs again in the future
* Update snakemake to the latest version (9.23). This fixed some of the issues with crashing snakemake
* Use `--scheduler greedy` to speed up the DAG construction
* Use `--allowed-rules` to focus the pipeline only on the relevant section and prevent the whole DAG from being constructed
* Remove all `checkpoint` from the pipeline
* Combine rules that occur sequentially into a single rule (TODO: concut rules)

In addition I found that several of the final steps were computationally impossible to run on a set of ~80,000 assembled genomes. Whenever there is going to be a series of pairwise computations on a set of this size we are going to run into problems. To solve this I upgraded checkm to checkm2 and I upgraded dRep to v4 which can use skani instead of ANI to greatly speed up the computation.
