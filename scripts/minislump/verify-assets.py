#!/usr/bin/env python3
"""Reject mixed mesh/transform or modified reproduction assets before debugging."""
import hashlib
import json
from pathlib import Path

assets = Path(__file__).resolve().parents[2] / 'src/demo/DFCMiniSlump'
manifest = json.loads((assets / 'sha256.json').read_text())
for name, expected in manifest.items():
    actual = hashlib.sha256((assets / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f'Asset changed: {name}; update the reproduction and its manifest deliberately')
transform = json.loads((assets / 'cone_chrono_transform.json').read_text())
if manifest['minislump-Cut4.obj'] != transform['mesh_sha256']:
    raise SystemExit('Cone transform does not match mesh')
print('MiniSlump asset hashes verified')
