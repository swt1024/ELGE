import os

def fix_gtf_wrapped_lines_inplace(gtf_path):
    """
    Fix wrapped/broken lines in a GTF file in-place (overwrite the original file).

    Problem:
    - A valid GTF line should have at least 9 tab-separated columns.
    - Sometimes the 9th column (attributes) is accidentally split into a new line.
      In that case, the next line usually has < 9 columns and should be appended
      to the previous line's attributes column.
    """
    tmp_path = gtf_path + ".tmp"

    fixed_lines = []
    buffer = None  # Store the current complete (or pending) GTF record (9 columns)

    with open(gtf_path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")

            # Keep empty lines and comment lines, but flush buffer first
            if not line.strip() or line.startswith("#"):
                if buffer is not None:
                    fixed_lines.append("\t".join(buffer))
                    buffer = None
                fixed_lines.append(line)
                continue

            parts = line.split("\t")

            if len(parts) >= 9:
                # New valid GTF record: flush the previous buffer first
                if buffer is not None:
                    fixed_lines.append("\t".join(buffer))
                buffer = parts[:9]  # Keep only the first 9 columns
            else:
                # This is likely a continuation of the previous attributes column
                if buffer is None:
                    # Edge case: continuation line appears without a previous record
                    fixed_lines.append(line)
                else:
                    cont = line.strip()
                    # Append with a space to avoid merging words together
                    buffer[8] = (buffer[8].rstrip() + " " + cont).strip()

        # Flush the last buffer at EOF
        if buffer is not None:
            fixed_lines.append("\t".join(buffer))

    # Write to a temporary file first
    with open(tmp_path, "w", encoding="utf-8") as out:
        out.write("\n".join(fixed_lines) + "\n")

    # Atomically replace the original file
    os.replace(tmp_path, gtf_path)


# Example usage (this will overwrite the original file)
fix_gtf_wrapped_lines_inplace("../human/gtf/NONCODEv6_human_hg38_lncRNA.gtf")
print("Done: original GTF file has been overwritten.")
