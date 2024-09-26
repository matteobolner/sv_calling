# sv_calling
Repository containing a Snakemake pipeline for structural variant calling


@nickp60 Have you bound the location of these external files? For example, I run snakemake from my home directory, which is automatically bound by snakemake + singularity (it binds $HOME), but I wanted it to also have access to my scratch user directory where my data and results are stored.

To do this, I ran: snakemake --use-singularity --use-conda --singularity-args "--bind /scratch/<myusername>/"
