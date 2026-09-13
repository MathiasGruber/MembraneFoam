#!/bin/bash
# Run from any directory after sourcing OpenCFD OpenFOAM v2606.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
mkdir -p "$root/runs"
work=$(mktemp -d "$root/runs/test-XXXXXXXX")
echo "Validation directory: $work"
python3 examples/run.py channel "$work/fo" --nx 24 --nz 24 --end 5000
python3 examples/run.py channel "$work/ro" --nx 24 --nz 48 --mode RO --end 5000
python3 examples/run.py channel "$work/transient" --nx 24 --nz 24 \
    --transient 0.0058 --time-step 0.0002
python3 examples/run.py channel "$work/mpi" --nx 24 --nz 24 --end 5000 --ranks 4
python3 - "$work" <<'PY'
import csv,json,re,sys
from pathlib import Path
root=Path(sys.argv[1]); results={c:json.loads((root/c/'result.json').read_text()) for c in ('fo','ro','mpi')}
for case,r in results.items():
    assert r['converged'], (case,'not converged')
    assert r['relative_mass_imbalance'] < 1e-6, (case,r)
    assert r['relative_salt_imbalance'] is not None and r['relative_salt_imbalance'] < 0.01, (case,'salt balance exceeds 1% of membrane transfer',r)
    assert abs(r['membrane_salt_imbalance_kg_s']) < 1e-15, (case,'membrane creates salt',r)
for case in ('fo','mpi'):
    with (root/case/'result-surface.csv').open() as source:
        rows=list(csv.DictReader(source))
    for sign in (-1,1):
        area=sum(float(row['area_m2']) for row in rows if sign*float(row['nz'])>0.5)
        assert abs(area-0.03*0.015)<1e-12,(case,sign,area)
for side in ('positive_z','negative_z'):
    a=results['fo']['surface_integrals'][side]['mean_mass_fraction']
    b=results['mpi']['surface_integrals'][side]['mean_mass_fraction']
    assert abs(a-b)/a<1e-5,(side,a,b)
a=results['fo']['water_mass_flux_kg_m2_h'];b=results['mpi']['water_mass_flux_kg_m2_h']
assert abs(a-b)/a<1e-5,(a,b)
transient=json.loads((root/'transient/result.json').read_text())
assert transient['run_type']=='transient' and not transient['converged'],transient
log=(root/'transient/log.solver').read_text()
initial=float(re.search(r'Initial salt mass \[kg\] = ([\d.eE+-]+)',log)[1])
# 0.5 mol/L NaCl × 58.44 g/mol × 0.03 × 0.015 × 0.001 m³.
assert abs(initial-1.3149e-5)<1e-10,(initial,'incorrect initial salt inventory')
audits=[dict((key,float(value)) for key,value in re.findall(r'(\w+)=([\d.eE+-]+)',line))
        for line in log.splitlines() if line.startswith('MEMBRANE_TRANSIENT ')]
# Repeated 0.0002 steps finish just below 0.0058 in floating-point arithmetic.
assert len(audits)==15,('missing transient balances',len(audits))
assert (root/'transient/0.0058/m_A').is_file(), 'Final off-schedule checkpoint missing'
assert abs(audits[-1]['time']-0.0058)<1e-12, audits[-1]
assert max(r['relative_salt_balance'] for r in audits)<0.01,audits
assert max(r['relative_mass_balance'] for r in audits)<1e-6,audits
print('Serial FO/RO, transient inventory/conservation, and four-rank FO checks passed')
PY
