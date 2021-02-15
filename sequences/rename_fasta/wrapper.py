from Bio import SeqIO
from hashlib import blake2s
from datetime import datetime


host = snakemake.params.get("host", "")
location = snakemake.params.get("location", "")
sample = snakemake.params.get("sample", "")
collection_date = snakemake.params.get("collection_date", "")
hexdigest = snakemake.params.get("hexdigest", False)

assert len(sample) > 0, "Sample name is missing"
assert len(location) > 0, "Location is missing"
assert len(collection_date) > 0, "Collection date is missing"
assert len(host) > 0, "Host is missing"

h = blake2s()
sample = h.hexdigest()[:10] if hexdigest else sample
isolate = f"SARS-CoV-2/{host}/{location}/{sample}/{collection_date}"

with open(snakemake.input[0], "r") as input_handle, open(
    snakemake.output[0], "w"
) as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        h.update(sample.encode("utf-8"))
        record.id = isolate
        record.description = ""
        SeqIO.write(record, output_handle, "fasta")
