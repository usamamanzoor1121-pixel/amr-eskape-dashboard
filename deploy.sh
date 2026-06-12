#!/bin/bash
# Run this ON the EC2 server after SSH-ing in
# Sets up Python, installs dependencies, runs the dashboard

set -e

echo "=== Updating system ==="
sudo apt-get update -y
sudo apt-get install -y python3-pip python3-venv git screen

echo "=== Cloning repo ==="
git clone https://github.com/usamamanzoor1121-pixel/amr-eskape-dashboard.git
cd amr-eskape-dashboard

echo "=== Setting up virtual environment ==="
python3 -m venv venv
source venv/bin/activate

echo "=== Installing dependencies ==="
pip install --upgrade pip
pip install pandas numpy requests plotly streamlit boto3 tqdm

echo "=== Starting dashboard ==="
screen -dmS dashboard streamlit run app/dashboard.py \
    --server.port 8501 \
    --server.address 0.0.0.0 \
    --server.headless true

echo "=== Done ==="
echo "Dashboard running on port 8501"
