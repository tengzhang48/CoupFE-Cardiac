"""Fail-closed result writing for cardiac time histories."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np


def save_completed(path, *, completed_steps, expected_steps, **payload):
    """Write an archive only after every requested step has completed."""
    if completed_steps != expected_steps:
        raise RuntimeError(
            f"incomplete solve: completed {completed_steps}/{expected_steps} steps"
        )
    destination = Path(path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{destination.name}.",
            suffix=".tmp.npz",
            dir=destination.parent,
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
        np.savez(
            temporary,
            completed_steps=int(completed_steps),
            expected_steps=int(expected_steps),
            converged=True,
            **payload,
        )
        os.replace(temporary, destination)
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
    return destination
