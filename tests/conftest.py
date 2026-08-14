from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "examples" / "cardiac_benchmark"
sys.path.insert(0, str(BENCHMARK))
