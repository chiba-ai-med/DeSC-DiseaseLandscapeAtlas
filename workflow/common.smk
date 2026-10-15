import pandas as pd

configfile: "config/config.yaml"

DISEASES_TSV = pd.read_csv(config["data"]["disease_families_tsv"], sep="\t")["family_id"].tolist()
DISEASES_EXTRA = config.get("disease_families_extra", [])
DISEASES = sorted(set(DISEASES_TSV + DISEASES_EXTRA))

PREPROC = config["preproc"]
DIM = str(config["pca"]["dim"])

CONTAINER_ANALYSIS = config["containers"]["analysis"]
CONTAINER_LANDSCAPER = config["containers"]["landscaper"]
CONTAINER_ONLINEPCA = config["containers"]["onlinepca"]
CONTAINER_MMJL = config["containers"]["mmjl"]

TEKIYO = config["data"]["tekiyo_path"]
HOSP = config["data"]["hospitalization_path"]
