"""
Script 2 (fixed): Clean and parse ESKAPE AMR dataset
Handles abbreviated phenotype codes: r/s/i/ns/nd/dd
"""

import pandas as pd
import numpy as np
import re
import os

INPUT_PATH = "data/raw/eskape_amr_raw.csv"
OUTPUT_DIR = "data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load ───────────────────────────────────────────────────────────────────────

print("Loading raw data...")
df = pd.read_csv(INPUT_PATH, low_memory=False)
print(f"Loaded: {len(df):,} records")

# ── Keep useful columns ────────────────────────────────────────────────────────

KEEP_COLS = [
    "pathogen", "scientific_name", "geo_loc_name", "collection_date",
    "isolation_source", "host", "host_disease", "epi_type",
    "number_drugs_resistant", "number_drugs_intermediate",
    "number_drugs_susceptible", "number_drugs_tested", "number_amr_genes",
    "AST_phenotypes", "AMR_genotypes", "AMR_genotypes_core", "biosample_acc",
]

df = df[[c for c in KEEP_COLS if c in df.columns]].copy()

# ── Extract year ───────────────────────────────────────────────────────────────

def extract_year(val):
    if pd.isna(val):
        return np.nan
    match = re.search(r'(19[5-9]\d|20[0-2]\d)', str(val))
    return int(match.group()) if match else np.nan

print("Extracting year...")
df["year"] = df["collection_date"].apply(extract_year)

# Filter realistic years only
df = df[df["year"].isna() | df["year"].between(1950, 2026)]
print(f"After year filter: {len(df):,} records")

# ── Clean country ──────────────────────────────────────────────────────────────

def extract_country(val):
    if pd.isna(val):
        return "Unknown"
    country = str(val).split(":")[0].split(",")[0].strip()
    return country if country else "Unknown"

print("Cleaning geographic data...")
df["country"] = df["geo_loc_name"].apply(extract_country)

# ── Phenotype code mapping ─────────────────────────────────────────────────────

PHENOTYPE_MAP = {
    "r"  : "resistant",
    "s"  : "susceptible",
    "i"  : "intermediate",
    "ns" : "nonsusceptible",
    "nd" : "not_defined",
    "dd" : "dose_dependent",
    # full words as fallback
    "resistant"      : "resistant",
    "susceptible"    : "susceptible",
    "intermediate"   : "intermediate",
    "nonsusceptible" : "nonsusceptible",
}

RESISTANT_CODES = {"resistant", "nonsusceptible"}

# ── Parse AST phenotypes ───────────────────────────────────────────────────────

print("Parsing AST phenotypes...")

def parse_ast(val):
    if pd.isna(val) or str(val).strip() in ("", "nan"):
        return {}
    result = {}
    for item in str(val).split(","):
        item = item.strip()
        if "=" in item:
            parts     = item.split("=", 1)
            antibiotic = parts[0].strip().lower()
            raw_code   = parts[1].strip().lower()
            phenotype  = PHENOTYPE_MAP.get(raw_code, raw_code)
            result[antibiotic] = phenotype
    return result

ast_records = []
for _, row in df.iterrows():
    ast_dict = parse_ast(row.get("AST_phenotypes"))
    if not ast_dict:
        continue
    for antibiotic, phenotype in ast_dict.items():
        ast_records.append({
            "biosample_acc"   : row.get("biosample_acc"),
            "pathogen"        : row["pathogen"],
            "country"         : row["country"],
            "year"            : row["year"],
            "isolation_source": row.get("isolation_source"),
            "host"            : row.get("host"),
            "antibiotic"      : antibiotic,
            "phenotype"       : phenotype,
        })

ast_long = pd.DataFrame(ast_records)
ast_long["is_resistant"] = ast_long["phenotype"].isin(RESISTANT_CODES)

print(f"AST long records   : {len(ast_long):,}")
print(f"Unique antibiotics : {ast_long['antibiotic'].nunique()}")
print(f"Phenotype counts:")
print(ast_long["phenotype"].value_counts().to_string())

# ── Resistance rates ───────────────────────────────────────────────────────────

print("\nComputing resistance rates...")

resistance_rates = (
    ast_long
    .groupby(["pathogen", "antibiotic"])
    .agg(
        total_tested    = ("is_resistant", "count"),
        total_resistant = ("is_resistant", "sum"),
    )
    .reset_index()
)
resistance_rates["resistance_rate"] = (
    resistance_rates["total_resistant"] / resistance_rates["total_tested"] * 100
).round(2)
resistance_rates = resistance_rates[resistance_rates["total_tested"] >= 100]

# ── Resistance trends over time ────────────────────────────────────────────────

print("Computing trends...")

trend_df = (
    ast_long[ast_long["year"].between(2000, 2026)]
    .groupby(["pathogen", "antibiotic", "year"])
    .agg(
        total_tested    = ("is_resistant", "count"),
        total_resistant = ("is_resistant", "sum"),
    )
    .reset_index()
)
trend_df["resistance_rate"] = (
    trend_df["total_resistant"] / trend_df["total_tested"] * 100
).round(2)
trend_df = trend_df[trend_df["total_tested"] >= 10]

# ── Country-level resistance ───────────────────────────────────────────────────

print("Computing country resistance...")

country_df = (
    ast_long[ast_long["country"] != "Unknown"]
    .groupby(["pathogen", "country", "antibiotic"])
    .agg(
        total_tested    = ("is_resistant", "count"),
        total_resistant = ("is_resistant", "sum"),
    )
    .reset_index()
)
country_df["resistance_rate"] = (
    country_df["total_resistant"] / country_df["total_tested"] * 100
).round(2)
country_df = country_df[country_df["total_tested"] >= 20]

# ── Pathogen summary ───────────────────────────────────────────────────────────

summary = df.groupby("pathogen").agg(
    total_isolates      = ("biosample_acc", "count"),
    isolates_with_AST   = ("AST_phenotypes", lambda x: x.notna().sum()),
    mean_drugs_resistant= ("number_drugs_resistant", "mean"),
    mean_drugs_tested   = ("number_drugs_tested", "mean"),
    mean_amr_genes      = ("number_amr_genes", "mean"),
    year_min            = ("year", "min"),
    year_max            = ("year", "max"),
).reset_index()

for col in ["mean_drugs_resistant", "mean_drugs_tested", "mean_amr_genes"]:
    summary[col] = summary[col].round(2)

# ── Save ───────────────────────────────────────────────────────────────────────

print("\nSaving files...")
df.to_csv(f"{OUTPUT_DIR}/eskape_clean.csv", index=False)
ast_long.to_csv(f"{OUTPUT_DIR}/ast_long.csv", index=False)
resistance_rates.to_csv(f"{OUTPUT_DIR}/resistance_rates.csv", index=False)
trend_df.to_csv(f"{OUTPUT_DIR}/resistance_trends.csv", index=False)
country_df.to_csv(f"{OUTPUT_DIR}/country_resistance.csv", index=False)
summary.to_csv(f"{OUTPUT_DIR}/pathogen_summary.csv", index=False)

print("\n" + "="*50)
print("PROCESSING COMPLETE")
print("="*50)
print(f"Clean isolates     : {len(df):,}")
print(f"AST long records   : {len(ast_long):,}")
print(f"Resistance pairs   : {len(resistance_rates)}")
print(f"Trend records      : {len(trend_df):,}")
print(f"Country records    : {len(country_df):,}")

print("\nPathogen summary:")
print(summary.to_string(index=False))

print("\nTop 15 highest resistance rates (min 100 tested):")
print(
    resistance_rates
    .sort_values("resistance_rate", ascending=False)
    .head(15)
    [["pathogen", "antibiotic", "total_tested", "resistance_rate"]]
    .to_string(index=False)
)
