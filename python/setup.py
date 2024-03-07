# Advice from: https://solutions.posit.co/write-code/minimum-viable-python/installing-packages/
# Create virtual environment
python3 -m venv .venv

# Activate virtual environment
source .venv/bin/activate

# Update pip
python3 -m pip install -U pip setuptools wheel

# Install sympy
python3 -m pip install sympy
