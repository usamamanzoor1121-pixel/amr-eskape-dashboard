"""
Script 1 (final): Fetch ESKAPE AMR data from NCBI Pathogen Detection FTP
Correctly handles dynamic versioned filenames
"""

import pandas as pd
import requests
import re
import os
from tqdm import tqdm

OUTPUT_DIR = "data/raw"
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results"

# Correct NCBI folder names confirmed from FTP
ESKAPE_FOLDERS = {
    "Staphylococcus aureus":   "Staphylococcus_aureus",
    "Klebsiella pneumoniae":   "Klebsiella",
    "Acinetobacter baumannii": "Acinetobacter",
    "Pseudomonas aeruginosa":  "Pseudomonas_aeruginosa",
    "Enterococcus faecium":    "Enterococcus_faecium",
    "Enterobacter cloacae":    "Enterobacter_cloacae",
}

def get_amr_filename(folder_name):
    """Discover the versioned AMR filename dynamically from FTP listing."""
    url = f"{BASE_URL}/{folder_name}/latest_snps/AMR/"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        # Find the .amr.metadata.tsv file
        match = re.search(r'href="(PDG\d+\.\d+\.amr\.metadata\.tsv)"', response.text)
        if match:
            return match.group(1)
        else:
            print(f"  Could not find AMR TSV file in: {url}")
            return None
    except Exception as e:
        print(f"  Failed to list AMR folder for {folder_name}: {e}")
        return None


def download_file(url, output_path, chunk_size=65536):
    """Download file with progress bar."""
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()
    total = int(response.headers.get("content-length", 0))
    with open(output_path, "wb") as f, tqdm(
        desc=os.path.basename(output_path)[:40],
        total=total,
        unit="B",
        unit_scale=True,
    ) as bar:
        for chunk in response.iter_content(chunk_size=chunk_size):
            f.write(chunk)
            bar.update(len(chunk))


def fetch_eskape_amr():
    all_dfs = []

    for common_name, folder_name in ESKAPE_FOLDERS.items():
        print(f"\n{'='*50}")
        print(f"Pathogen: {common_name}")

        # Step 1: discover filename
        amr_filename = get_amr_filename(folder_name)
        if not amr_filename:
            continue

        # Step 2: download if not cached
        out_path = os.path.join(OUTPUT_DIR, f"{folder_name}_amr.tsv")
        if not os.path.exists(out_path):
            url = f"{BASE_URL}/{folder_name}/latest_snps/AMR/{amr_filename}"
            print(f"  Downloading: {amr_filename}")
            try:
                download_file(url, out_path)
            except Exception as e:
                print(f"  Download failed: {e}")
                continue
        else:
            print(f"  Already cached: {out_path}")

        # Step 3: load and tag
        try:
            df = pd.read_csv(out_path, sep="\t", low_memory=False)
            df["pathogen"] = common_name
            all_dfs.append(df)
            print(f"  Records loaded : {len(df):,}")
            print(f"  Columns        : {len(df.columns)}")
        except Exception as e:
            print(f"  Failed to load: {e}")

    return all_dfs


def main():
    print("ESKAPE AMR Data Fetch — NCBI Pathogen Detection")
    print("=" * 50)

    all_dfs = fetch_eskape_amr()

    if not all_dfs:
        print("\nNo data fetched.")
        return

    # Combine all pathogens
    combined = pd.concat(all_dfs, ignore_index=True)

    # Summary
    print(f"\n{'='*50}")
    print("COMBINED DATASET SUMMARY")
    print(f"{'='*50}")
    print(f"Total records    : {len(combined):,}")
    print(f"Pathogens        : {combined['pathogen'].nunique()}")
    print(f"Columns          : {len(combined.columns)}")
    print(f"\nRecords per pathogen:")
    print(combined["pathogen"].value_counts().to_string())
    print(f"\nAll columns:")
    for col in combined.columns:
        print(f"  {col}")

    # Save
    out_path = os.path.join(OUTPUT_DIR, "eskape_amr_raw.csv")
    combined.to_csv(out_path, index=False)
    print(f"\nSaved to: {out_path}")


if __name__ == "__main__":
    main()
