#!/usr/bin/env bash
# Install the app and the exact public CoupFE core commit declared in pyproject.toml.
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
python -m pip install -e "${APP_DIR}[dev]"
python "${APP_DIR}/.github/scripts/check_runtime_core.py"
python -m pytest -q "${APP_DIR}/tests"
