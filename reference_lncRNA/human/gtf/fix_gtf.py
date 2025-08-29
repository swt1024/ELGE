#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fix abnormal line breaks in the attribute column (9th column) of a GTF file, such as:
gene_id "XXX
"; transcript_id "YYY";

Features:
1) Merge "logical lines" that are incorrectly split into multiple physical lines.
2) Specifically fix cases where a newline/extra whitespace appears between the quote and semicolon (e.g., `"value\n";` → `";`).
3) Preserve other fields and whitespace structure as much as possible.
4) Output a repair log containing start/end line numbers, original text, and fixed text.

Usage:
python fix_gtf_with_log.py input.gtf output.gtf --log fixed.log
"""

import re
import argparse

def is_complete_logical_line(buf: str) -> bool:
    """
    Check whether the current buffer has formed a complete "logical" GTF line:
    - Must have at least 9 columns (8 tab separators),
    - The attribute column (9th) must have an even number of quotes (avoid unmatched quotes).
    Note: We do not strictly check if each attribute ends with a semicolon, because this will be fixed later.
    """
    if buf.count("\t") < 8:
        return False
    # Allow TAB > 8 (some tools may have tabs inside the attribute column)
    # Only require that the 9th column exists and quotes are balanced
    attr = buf.split("\t", 8)[-1]
    return (attr.count('"') % 2 == 0)

def minimal_attr_fix(attr: str) -> tuple[str, bool]:
    """
    Perform minimal fixes on the attribute column:
    1) Fix newline/whitespace between the quote and semicolon:
       "value \n   ;" → "\"value\"; "
    2) Remove extra whitespace before semicolon.
    3) Ensure at least one space after each semicolon.
    4) **New**: Trim leading/trailing spaces inside all quoted values
       (fix cases like: gene_id "xxx ").
    """
    original = attr

    # 1) Fix newline/whitespace between quote and semicolon
    attr = re.sub(r'"\s*;\s*', r'"; ', attr)

    # 2) Remove extra whitespace before semicolon
    attr = re.sub(r'\s+;', r';', attr)

    # 4) Trim spaces inside quotes (without affecting valid spaces in the middle)
    attr = re.sub(r'"([^"]*)"', lambda m: f'"{m.group(1).strip()}"', attr)

    # 3) Add one space after semicolon if followed by a non-space character
    attr = re.sub(r';(?=\S)', r'; ', attr)

    changed = (attr != original)
    return attr, changed


def process_gtf(input_path: str, output_path: str, log_path: str):
    total_lines = 0
    fixed_count = 0
    logical_line_no = 0  # Count of logical lines
    log_entries = []

    with open(input_path, "r", encoding="utf-8", newline="") as fin, \
         open(output_path, "w", encoding="utf-8", newline="\n") as fout:

        buf = ""                # Buffer for the current logical line (may span multiple physical lines)
        buf_start_line = None   # Start line number of the current buffer
        pending_original = ""   # Original concatenated text (for logging)

        for lineno, raw in enumerate(fin, start=1):
            total_lines += 1
            line = raw.rstrip("\n")

            # If comment line and not currently buffering a logical line, write directly
            if line.startswith("#") and not buf:
                fout.write(line + "\n")
                continue

            # Start or continue building a logical line
            if not buf:
                buf = line
                buf_start_line = lineno
                pending_original = line
            else:
                # For continuation lines: join with a space for logging
                buf += " " + line
                pending_original += "\n" + line

            # Check if we have a complete logical line
            if not is_complete_logical_line(buf):
                continue

            # At this point: buf is a complete logical line
            logical_line_no += 1

            # Split into first 8 columns + attribute column
            if buf.count("\t") >= 8:
                head, attr = buf.split("\t", 8)[0:8], buf.split("\t", 8)[8]
                fixed_attr, changed = minimal_attr_fix(attr)

                fixed_line = "\t".join(head + [fixed_attr])

                # Log if original differs from fixed (including multi-line merge cases)
                if pending_original != fixed_line or changed:
                    fixed_count += 1
                    log_entries.append({
                        "start": buf_start_line,
                        "end": lineno,
                        "orig": pending_original,
                        "fixed": fixed_line
                    })

                # Output fixed line
                fout.write(fixed_line + "\n")

            else:
                # Fallback: output as is for rare malformed cases
                fout.write(buf + "\n")

            # Reset buffer
            buf = ""
            buf_start_line = None
            pending_original = ""

        # If leftover buffer at end of file, try to output and log it
        if buf:
            logical_line_no += 1
            if buf.count("\t") >= 8:
                head, attr = buf.split("\t", 8)[0:8], buf.split("\t", 8)[8]
                fixed_attr, changed = minimal_attr_fix(attr)
                fixed_line = "\t".join(head + [fixed_attr])

                if pending_original != fixed_line or changed:
                    fixed_count += 1
                    log_entries.append({
                        "start": buf_start_line,
                        "end": total_lines,
                        "orig": pending_original,
                        "fixed": fixed_line
                    })

                fout.write(fixed_line + "\n")
            else:
                fout.write(buf + "\n")


def main():
    parser = argparse.ArgumentParser(description="Fix GTF files with attribute column broken across lines, and generate a repair log.")
    parser.add_argument("input", help="Input GTF file")
    parser.add_argument("output", help="Fixed GTF output file")
    args = parser.parse_args()

    process_gtf(args.input, args.output, args.log)


if __name__ == "__main__":
    main()
