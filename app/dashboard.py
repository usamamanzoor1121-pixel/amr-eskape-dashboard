"""
AMR ESKAPE Surveillance Dashboard
Interactive Streamlit app for exploring AMR trends across ESKAPE pathogens
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

# ── Page config ────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="ESKAPE AMR Dashboard",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Load data ──────────────────────────────────────────────────────────────────

DATA_DIR = "data/processed"

@st.cache_data
def load_data():
    resistance_rates = pd.read_csv(f"{DATA_DIR}/resistance_rates.csv")
    trend_df         = pd.read_csv(f"{DATA_DIR}/resistance_trends.csv")
    country_df       = pd.read_csv(f"{DATA_DIR}/country_resistance.csv")
    summary          = pd.read_csv(f"{DATA_DIR}/pathogen_summary.csv")
    ast_long         = pd.read_csv(f"{DATA_DIR}/ast_long.csv")
    return resistance_rates, trend_df, country_df, summary, ast_long

resistance_rates, trend_df, country_df, summary, ast_long = load_data()

PATHOGENS = sorted(resistance_rates["pathogen"].unique().tolist())

# ── Colour palette ─────────────────────────────────────────────────────────────

PATHOGEN_COLORS = {
    "Staphylococcus aureus"   : "#e63946",
    "Klebsiella pneumoniae"   : "#2a9d8f",
    "Acinetobacter baumannii" : "#e9c46a",
    "Pseudomonas aeruginosa"  : "#457b9d",
    "Enterococcus faecium"    : "#f4a261",
    "Enterobacter cloacae"    : "#a8dadc",
}

# ── Header ─────────────────────────────────────────────────────────────────────

st.title("🦠 ESKAPE Pathogen AMR Surveillance Dashboard")
st.markdown(
    "Antimicrobial resistance trends across **ESKAPE pathogens** "
    "using NCBI Pathogen Detection data (504,096 isolates · 6 pathogens · 110 antibiotics)"
)
st.divider()

# ── Sidebar ────────────────────────────────────────────────────────────────────

st.sidebar.header("Filters")

selected_pathogens = st.sidebar.multiselect(
    "Select Pathogens",
    options=PATHOGENS,
    default=PATHOGENS,
)

min_tested = st.sidebar.slider(
    "Minimum isolates tested (per antibiotic)",
    min_value=50,
    max_value=500,
    value=100,
    step=50,
)

st.sidebar.divider()
st.sidebar.markdown("**Data Source**")
st.sidebar.markdown("[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/)")
st.sidebar.markdown("Updated: June 2026")
st.sidebar.divider()
st.sidebar.markdown("Built by **Usama Manzoor**")
st.sidebar.markdown("[GitHub](https://github.com/usamamanzoor1121-pixel)")

# ── Tab layout ─────────────────────────────────────────────────────────────────

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Overview",
    "🔬 Resistance Rates",
    "📈 Trends Over Time",
    "🌍 Geographic Distribution",
    "🧬 AMR Gene Burden",
])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — OVERVIEW
# ══════════════════════════════════════════════════════════════════════════════

with tab1:
    st.subheader("Dataset Overview")

    # KPI cards
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Isolates",   f"{summary['total_isolates'].sum():,}")
    col2.metric("Pathogens",        f"{len(summary)}")
    col3.metric("Antibiotics Tested", f"{resistance_rates['antibiotic'].nunique()}")
    col4.metric("Countries",        f"{country_df['country'].nunique()}")

    st.divider()

    # Isolates per pathogen bar chart
    fig_iso = px.bar(
        summary.sort_values("total_isolates", ascending=True),
        x="total_isolates",
        y="pathogen",
        orientation="h",
        color="pathogen",
        color_discrete_map=PATHOGEN_COLORS,
        title="Total Isolates per Pathogen",
        labels={"total_isolates": "Number of Isolates", "pathogen": ""},
        text="total_isolates",
    )
    fig_iso.update_traces(texttemplate="%{text:,}", textposition="outside")
    fig_iso.update_layout(showlegend=False, height=350)
    st.plotly_chart(fig_iso, use_container_width=True)

    st.divider()

    # Summary table
    st.subheader("Pathogen Summary Statistics")
    display_summary = summary.copy()
    display_summary.columns = [
        "Pathogen", "Total Isolates", "Isolates with AST",
        "Mean Drugs Resistant", "Mean Drugs Tested",
        "Mean AMR Genes", "Year Min", "Year Max"
    ]
    st.dataframe(display_summary, use_container_width=True, hide_index=True)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — RESISTANCE RATES
# ══════════════════════════════════════════════════════════════════════════════

with tab2:
    st.subheader("Resistance Rates by Pathogen and Antibiotic")

    filtered_rates = resistance_rates[
        (resistance_rates["pathogen"].isin(selected_pathogens)) &
        (resistance_rates["total_tested"] >= min_tested)
    ].copy()

    if filtered_rates.empty:
        st.warning("No data for current filters.")
    else:
        # Top resistant combinations
        top_resistant = filtered_rates.sort_values(
            "resistance_rate", ascending=False
        ).head(20)

        fig_top = px.bar(
            top_resistant,
            x="resistance_rate",
            y="antibiotic",
            color="pathogen",
            orientation="h",
            color_discrete_map=PATHOGEN_COLORS,
            title="Top 20 Pathogen-Antibiotic Resistance Rates",
            labels={
                "resistance_rate": "Resistance Rate (%)",
                "antibiotic": "Antibiotic",
                "pathogen": "Pathogen",
            },
            hover_data=["total_tested", "total_resistant"],
        )
        fig_top.update_layout(height=550, yaxis={"categoryorder": "total ascending"})
        st.plotly_chart(fig_top, use_container_width=True)

        st.divider()

        # Heatmap — pathogen vs antibiotic
        st.subheader("Resistance Heatmap")

        pivot = filtered_rates.pivot_table(
            index="antibiotic",
            columns="pathogen",
            values="resistance_rate",
            aggfunc="mean",
        ).fillna(0)

        fig_heat = px.imshow(
            pivot,
            color_continuous_scale="RdYlGn_r",
            title="Resistance Rate Heatmap (%) — Antibiotic × Pathogen",
            labels={"color": "Resistance Rate (%)"},
            aspect="auto",
            zmin=0,
            zmax=100,
        )
        fig_heat.update_layout(height=700)
        st.plotly_chart(fig_heat, use_container_width=True)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — TRENDS OVER TIME
# ══════════════════════════════════════════════════════════════════════════════

with tab3:
    st.subheader("Resistance Trends Over Time")

    col_left, col_right = st.columns(2)

    with col_left:
        selected_pathogen_trend = st.selectbox(
            "Select Pathogen",
            options=selected_pathogens,
            key="trend_pathogen",
        )

    # Get top antibiotics for this pathogen by number of records
    top_abx = (
        trend_df[trend_df["pathogen"] == selected_pathogen_trend]
        .groupby("antibiotic")["total_tested"]
        .sum()
        .sort_values(ascending=False)
        .head(10)
        .index.tolist()
    )

    with col_right:
        selected_abx_trend = st.multiselect(
            "Select Antibiotics",
            options=top_abx,
            default=top_abx[:5],
            key="trend_abx",
        )

    trend_filtered = trend_df[
        (trend_df["pathogen"] == selected_pathogen_trend) &
        (trend_df["antibiotic"].isin(selected_abx_trend))
    ]

    if trend_filtered.empty:
        st.warning("No trend data for this selection.")
    else:
        fig_trend = px.line(
            trend_filtered,
            x="year",
            y="resistance_rate",
            color="antibiotic",
            markers=True,
            title=f"Resistance Trends — {selected_pathogen_trend}",
            labels={
                "resistance_rate": "Resistance Rate (%)",
                "year": "Year",
                "antibiotic": "Antibiotic",
            },
        )
        fig_trend.update_layout(height=450)
        st.plotly_chart(fig_trend, use_container_width=True)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 — GEOGRAPHIC DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════

with tab4:
    st.subheader("Geographic Distribution of Resistance")

    col_l, col_r = st.columns(2)
    with col_l:
        geo_pathogen = st.selectbox(
            "Select Pathogen",
            options=selected_pathogens,
            key="geo_pathogen",
        )
    
    geo_abx_options = (
        country_df[country_df["pathogen"] == geo_pathogen]["antibiotic"]
        .value_counts()
        .head(20)
        .index.tolist()
    )
    with col_r:
        geo_antibiotic = st.selectbox(
            "Select Antibiotic",
            options=geo_abx_options,
            key="geo_abx",
        )

    geo_filtered = country_df[
        (country_df["pathogen"] == geo_pathogen) &
        (country_df["antibiotic"] == geo_antibiotic)
    ].sort_values("resistance_rate", ascending=False)

    if geo_filtered.empty:
        st.warning("No geographic data for this selection.")
    else:
        fig_geo = px.bar(
            geo_filtered.head(30),
            x="country",
            y="resistance_rate",
            color="resistance_rate",
            color_continuous_scale="RdYlGn_r",
            title=f"{geo_pathogen} — {geo_antibiotic} resistance by country",
            labels={
                "resistance_rate": "Resistance Rate (%)",
                "country": "Country",
            },
            hover_data=["total_tested"],
        )
        fig_geo.update_layout(height=450, xaxis_tickangle=-45)
        st.plotly_chart(fig_geo, use_container_width=True)

        st.dataframe(
            geo_filtered[["country", "total_tested", "total_resistant", "resistance_rate"]],
            use_container_width=True,
            hide_index=True,
        )

# ══════════════════════════════════════════════════════════════════════════════
# TAB 5 — AMR GENE BURDEN
# ══════════════════════════════════════════════════════════════════════════════

with tab5:
    st.subheader("AMR Gene Burden Across Pathogens")

    clean_df = pd.read_csv(f"{DATA_DIR}/eskape_clean.csv")
    clean_filtered = clean_df[
        (clean_df["pathogen"].isin(selected_pathogens)) &
        (clean_df["number_amr_genes"].notna())
    ]

    fig_box = px.box(
        clean_filtered,
        x="pathogen",
        y="number_amr_genes",
        color="pathogen",
        color_discrete_map=PATHOGEN_COLORS,
        title="Distribution of AMR Gene Count per Isolate",
        labels={
            "number_amr_genes": "Number of AMR Genes",
            "pathogen": "",
        },
        points=False,
    )
    fig_box.update_layout(
        height=450,
        showlegend=False,
        xaxis_tickangle=-20,
    )
    st.plotly_chart(fig_box, use_container_width=True)

    st.divider()

    # Top AMR genes across all pathogens
    st.subheader("Most Common AMR Genes")

    ast_genes = pd.read_csv(f"{DATA_DIR}/ast_long.csv")
    amr_gene_df = clean_df[
        clean_df["pathogen"].isin(selected_pathogens) &
        clean_df["AMR_genotypes"].notna()
    ].copy()

    gene_records = []
    for _, row in amr_gene_df.iterrows():
        genes = str(row["AMR_genotypes"]).split(",")
        for g in genes:
            g = g.strip()
            if g and g != "nan":
                gene_records.append({"pathogen": row["pathogen"], "gene": g})

    if gene_records:
        gene_df = pd.DataFrame(gene_records)
        top_genes = (
            gene_df.groupby(["pathogen", "gene"])
            .size()
            .reset_index(name="count")
            .sort_values("count", ascending=False)
        )
        top_genes_display = top_genes.groupby("pathogen").head(5)

        fig_genes = px.bar(
            top_genes_display,
            x="count",
            y="gene",
            color="pathogen",
            orientation="h",
            color_discrete_map=PATHOGEN_COLORS,
            title="Top 5 AMR Genes per Pathogen",
            labels={"count": "Isolate Count", "gene": "AMR Gene"},
            facet_col="pathogen",
            facet_col_wrap=3,
        )
        fig_genes.update_layout(height=600, showlegend=False)
        fig_genes.update_yaxes(matches=None)
        st.plotly_chart(fig_genes, use_container_width=True)
    else:
        st.info("No AMR gene data available for selected pathogens.")

