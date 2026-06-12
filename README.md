# ESKAPE Pathogen AMR Surveillance Dashboard

Interactive antimicrobial resistance surveillance dashboard for ESKAPE pathogens built with Python, Plotly, and Streamlit, deployed on AWS EC2.

## Live Demo

http://3.81.108.239:8501

## Overview

- 504,096 isolates across 6 ESKAPE pathogens
- 110 antibiotics with resistance phenotype data
- 16 countries with geographic resistance patterns
- Data source: NCBI Pathogen Detection (updated June 2026)

## Features

- Resistance rate heatmap across pathogen-antibiotic combinations
- Temporal resistance trends from 2000 to 2026
- Geographic distribution of resistance by country
- AMR gene burden analysis per pathogen

## ESKAPE Pathogens Covered

- Enterococcus faecium
- Staphylococcus aureus
- Klebsiella pneumoniae
- Acinetobacter baumannii
- Pseudomonas aeruginosa
- Enterobacter cloacae

## Tech Stack

- Data: NCBI Pathogen Detection FTP
- Analysis: Python, pandas, numpy
- Visualization: Plotly, Streamlit
- Cloud: AWS EC2 t3.micro with Elastic IP
- Version Control: GitHub

## Local Setup

git clone https://github.com/usamamanzoor1121-pixel/amr-eskape-dashboard.git
cd amr-eskape-dashboard
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
streamlit run app/dashboard.py

## Project Structure

amr-eskape-dashboard/
├── data/
│   ├── raw/               Raw NCBI downloads, not tracked
│   └── processed/         Cleaned analysis files, not tracked
├── scripts/
│   ├── fetch_amr_data.py  NCBI data acquisition
│   └── clean_amr_data.py  Data cleaning and processing
├── app/
│   └── dashboard.py       Streamlit dashboard
├── deploy.sh              EC2 deployment script
└── requirements.txt

## Author

Usama Manzoor
JSMU Diagnostic Laboratory and Blood Bank, Karachi, Pakistan
GitHub: https://github.com/usamamanzoor1121-pixel
LinkedIn: https://www.linkedin.com/in/usama-manzoor-042595182/
