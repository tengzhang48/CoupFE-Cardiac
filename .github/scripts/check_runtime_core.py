"""Verify that the installed CoupFE came from the audited public Git revision."""
from __future__ import annotations

import importlib.metadata
import json
from pathlib import Path

import coupfe


PUBLIC_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"
PUBLIC_CORE_REF = "e2f42ed5772850a0a23a2ce434f430c287eae5c8"


def main():
    distribution = importlib.metadata.distribution("coupfe")
    module_file = getattr(coupfe, "__file__", None)
    if module_file is None:
        raise SystemExit("the imported coupfe module has no source file")
    try:
        imported_path = Path(module_file).resolve(strict=True)
        distributed_path = Path(
            distribution.locate_file("coupfe/__init__.py")
        ).resolve(strict=True)
    except (AttributeError, OSError, RuntimeError) as error:
        raise SystemExit(
            "could not bind the imported coupfe module to its installed distribution"
        ) from error
    if imported_path != distributed_path:
        raise SystemExit(
            "the imported coupfe module is not supplied by the distribution whose "
            "PEP 610 record is being checked"
        )
    raw = distribution.read_text("direct_url.json")
    if raw is None:
        raise SystemExit(
            "installed CoupFE has no direct_url.json; install the exact VCS dependency "
            "declared by CoupFE-Cardiac"
        )
    try:
        record = json.loads(raw)
    except json.JSONDecodeError as error:
        raise SystemExit("installed CoupFE has invalid direct_url.json") from error
    vcs = record.get("vcs_info", {})
    if record.get("url") != PUBLIC_CORE_URL:
        raise SystemExit(
            f"installed CoupFE URL is {record.get('url')!r}; expected {PUBLIC_CORE_URL!r}"
        )
    if vcs.get("vcs") != "git" or vcs.get("commit_id") != PUBLIC_CORE_REF:
        raise SystemExit(
            "installed CoupFE does not resolve to the audited public commit "
            f"{PUBLIC_CORE_REF}"
        )

    print(f"CoupFE runtime verified: {PUBLIC_CORE_URL}@{PUBLIC_CORE_REF}")


if __name__ == "__main__":
    main()
