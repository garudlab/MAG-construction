import argparse
import csv
from tqdm import tqdm
from pathlib import Path

STATS_TO_KEEP = [
    "# contigs",
    "Total length",
    "Largest contig",
    "N50",
    "N90",
    "L50",
    "GC (%)",
]

def parse_quast_report(report_path):
    stats = {}
    with open(report_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            key, value = parts[0], parts[1]
            if key in STATS_TO_KEEP:
                stats[key] = value
    return stats

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dastool_dir", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("quast_summary.csv"))
    args = parser.parse_args()


    rows = []
    for report_path in tqdm(args.dastool_dir.glob("*/quast/*/report.tsv")):
        batch_or_sample = report_path.parents[2].name
        bin_name = report_path.parent.name
        stats = parse_quast_report(report_path)
        rows.append({
            "batch_or_sample": batch_or_sample,
            "bin": bin_name,
            **stats,
        })

    fieldnames = ["batch_or_sample", "bin"] + STATS_TO_KEEP
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")

if __name__ == "__main__":
    main()
