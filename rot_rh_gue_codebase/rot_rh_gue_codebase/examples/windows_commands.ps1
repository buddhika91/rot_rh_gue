# Install
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e .

# Smoke test
python scripts/run_frozen_audit.py `
  --N-list 20000,24000 `
  --jacobi-dims 64,96 `
  --controls 1 `
  --out-prefix smoke_frozen_audit

# Main endpoint-to-96000 frozen audit
python scripts/run_frozen_audit.py `
  --N-list 20000,24000,28000,32000,36000,40000,44000,48000,56000,64000,80000,96000 `
  --jacobi-dims 64,96,128,160,192,224 `
  --controls 32 `
  --seed 2027 `
  --out-prefix reflect_frozen_endpoint_96000_32_seed2027

# Control family separation
python scripts/control_family_scoreboard.py `
  --prefix reflect_frozen_endpoint_96000_32_seed2027
