import csv

# === Configuration ===
noncode_file = "../mapping/noncode.tsv"
ensembl_file = "../mapping/ensembl.tsv"

# Prefix rules for human and mouse
HUMAN_NONCODE_PREFIX = "NONHSAG"
HUMAN_ENSEMBL_PREFIX = "ENSG"
MOUSE_NONCODE_PREFIX = "NONMMUG"
MOUSE_ENSEMBL_PREFIX = "ENSMUSG"

def strip_version(id_str: str) -> str:
    """
    Remove version suffix after '.' if present.
    Example: ENSG000001.10 -> ENSG000001
             NONHSAG0001.2 -> NONHSAG0001
    """
    return id_str.split('.')[0]

def read_bridge(path, id_col_index=5):
    """
    Read RNAcentral bridge file.
    Returns a dict: {RNAcentral_ID: set of database IDs}.
    - Column 0 = RNAcentral ID
    - Column id_col_index = database ID (NONCODE or Ensembl)
    """
    mapping = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) <= id_col_index:
                continue
            rna_id = row[0].strip()
            db_id = strip_version(row[id_col_index].strip())
            if not rna_id or not db_id:
                continue
            mapping.setdefault(rna_id, set()).add(db_id)
    return mapping

def join_mappings(rna2noncode, rna2ensembl):
    """
    Inner-join on RNAcentral ID.
    Returns list of tuples: (noncode_id, ensembl_id, rna_id).
    - If one RNAcentral ID maps to multiple IDs, all combinations are included.
    - Duplicates are removed.
    """
    pairs = []
    for rna_id in set(rna2noncode) & set(rna2ensembl):
        for n in rna2noncode[rna_id]:
            for e in rna2ensembl[rna_id]:
                # Ensure IDs are version-stripped
                pairs.append((strip_version(n), strip_version(e), rna_id))
    return list(set(pairs))  # remove duplicates

def is_human(noncode_id, ensembl_id):
    """Return True if the pair is human (NONHSAG ↔ ENSG)."""
    return noncode_id.startswith(HUMAN_NONCODE_PREFIX) and ensembl_id.startswith(HUMAN_ENSEMBL_PREFIX)

def is_mouse(noncode_id, ensembl_id):
    """Return True if the pair is mouse (NONMMUG ↔ ENSMUSG)."""
    return noncode_id.startswith(MOUSE_NONCODE_PREFIX) and ensembl_id.startswith(MOUSE_ENSEMBL_PREFIX)

# === Main workflow ===
rna2noncode = read_bridge(noncode_file, id_col_index=5)
rna2ensembl = read_bridge(ensembl_file, id_col_index=5)

pairs = join_mappings(rna2noncode, rna2ensembl)

# Filter by species
human_pairs = [p for p in pairs if is_human(p[0], p[1])]
mouse_pairs = [p for p in pairs if is_mouse(p[0], p[1])]
all_pairs = human_pairs + mouse_pairs

# Save results
with open("../mapping/mapping_human.csv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["noncode_id", "ensembl_id", "rnacentral_id"])
    w.writerows(human_pairs)

with open("../mapping/mapping_mouse.csv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["noncode_id", "ensembl_id", "rnacentral_id"])
    w.writerows(mouse_pairs)

print("Total mappings:", len(all_pairs))
print("Human mappings:", len(human_pairs))
print("Mouse mappings:", len(mouse_pairs))
