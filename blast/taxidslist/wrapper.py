import subprocess
from more_itertools import always_iterable

# Getting taxa for query
params = snakemake.params

# Getting taxid lists
for k, v in params.items():
    for i in always_iterable(v):
        p = subprocess.run(
            ["get_species_taxids.sh", "-t", str(i)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
        if "esearch error" in p.stdout:
            raise Exception(p.stdout)
        else:
            with open(snakemake.output[k], "wt") as f:
                f.write(p.stdout)
