"""Distributed companion for Cardiac Benchmark 1 Step 0 and Step 2 Case B.

The retained ``run.py`` remains the serial example.  This companion now
supports the like-for-like closed Step 0 configuration for Cases A and B:
five-block closed geometry, standard Hex8 material with pointwise ``kappa`` or
an explicitly selected condensed local-pressure law, assembled consistent
mass, fixed ``eta=100 Pa s``, and GP-direct fibers from analytic or validated
Laplace transmural coordinates.  The closed MPI contract requires the
validated Laplace field; the analytic rule remains a serial diagnostic option.
The historical truncated polar-ring Q1/P0, lumped-mass, CG1 configuration
remains available through explicit command-line choices.

Closed Step 0 Cases A/B and Step 2 Case B have an explicit source-matched
generalized-alpha path. It uses the fixed Simula parameters, stages consistent
inertia at ``alpha_m`` and all force terms and loads at ``alpha_f``, and
selects a dedicated generated material whose viscous rate is reconstructed
from stage velocity. Historical and closed backward-Euler identities remain
separate.

Only rank zero prints progress and writes a completed result archive.  Every
rank compiles into an isolated build directory and owns only its local compiled
element state.  The global mesh is intentionally replicated for this benchmark
companion and for the replicated Robin/follower-pressure surface evaluation.
"""
from __future__ import annotations

import argparse
import atexit
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import geometry as geom
import run as serial_driver
from benchmark_parameters import (
    benchmark_configuration,
    benchmark_metadata,
    validate_load_histories,
    validate_runtime_material_properties,
)
from boundary_audit import audit_geometry, audit_pressure, audit_robin
from distributed_local_pressure import FusedQ1P0Batch
from distributed_mass import assemble_owned_consistent_mass
from distributed_material import RankLocalStandardMaterialBatch
from distributed_solver import (
    DIRECT_SUPERLU_DIST_PROFILE,
    FGMRES_ASM_ILU1_PROFILE,
    FGMRES_ASM_LU_PROFILE,
    FGMRES_GAMG_RIGID_PROFILE,
    FGMRES_GAMG_RIGID_REBUILD_PROFILE,
    DistributedPetscSnesSolver,
    DistributedSnesSolveError,
    _load_parallel_runtime,
    linear_solver_requires_node_aligned_ownership,
)
from generalized_alpha import SOURCE_MATCHED_GENERALIZED_ALPHA
from local_pressure import (
    LOG_MEAN_DILATATION_LAW,
    PAPER_MEAN_DILATATION_LAW,
)
from material import (
    CardiacHex8,
    CardiacHex8GeneralizedAlpha,
    HO_PARAMS,
    MATERIAL_MODEL_ID,
    build_kernel,
    build_generalized_alpha_kernel,
    struct_tensors,
)
from pressure import FollowerPressureOperator
from result_io import save_completed
from robin import RobinOperator
from sampling import interpolate_displacement, locate_hex8_point

from coupfe.assembly.distributed import element_partition, partition_elements
from coupfe.mesh import KernelMeshView
from coupfe.runtime.compiled_element import CompiledElement


LOCAL_PRESSURE_FORMULATIONS = frozenset(
    {"local-pressure", "local-pressure-paper"}
)
LOCAL_PRESSURE_ARCHIVE_LAWS = {
    "local-pressure": "linear-reference-volume-mean-log-j-v1",
    "local-pressure-paper": (
        "paper-j2-of-reference-volume-weighted-geometric-mean-j-v1"
    ),
}


def _local_pressure_law(formulation):
    return {
        "local-pressure": LOG_MEAN_DILATATION_LAW,
        "local-pressure-paper": PAPER_MEAN_DILATATION_LAW,
    }[formulation]


def _sanitize_failure_json(value):
    """Return a strict-JSON copy without hiding an original solver failure."""
    if isinstance(value, dict):
        return {
            str(key): _sanitize_failure_json(item)
            for key, item in value.items()
        }
    if isinstance(value, (tuple, list)):
        return [_sanitize_failure_json(item) for item in value]
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return number if np.isfinite(number) else None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def _init_fiber_svars_subset(
    element,
    mesh,
    schema,
    per_gp,
    global_element_ids,
    *,
    sampling="cg1",
    tbar=None,
    asign=-1.0,
):
    """Initialize rank-local structural tensors exactly like the serial path."""
    if sampling not in {"cg1", "gp-direct"}:
        raise ValueError(f"unsupported fiber sampling {sampling!r}")
    element_ids = np.asarray(global_element_ids, dtype=np.int64)
    if element_ids.ndim != 1:
        raise ValueError("global_element_ids must be one-dimensional")
    if len(element_ids) != len(element.svars):
        raise ValueError(
            "global_element_ids must match the rank-local compiled element count"
        )
    if len(element_ids) and (
        np.min(element_ids) < 0 or np.max(element_ids) >= mesh.n_elem
    ):
        raise ValueError("global_element_ids contains an element outside the mesh")
    n_gp = element.svars.shape[1] // per_gp
    offsets = {name: entry["offset"] for name, entry in schema.items()}
    shape = serial_driver._hex8_gp_shape()
    if n_gp != len(shape):
        raise ValueError(f"expected 8 Hex8 Gauss points, found {n_gp}")
    if tbar is None:
        tbar = mesh.param[:, 0]
    tbar = np.asarray(tbar, dtype=float)
    if tbar.shape != (mesh.n_node,) or not np.all(np.isfinite(tbar)):
        raise ValueError("fiber transmural coordinates must be finite per-node data")
    for local_index, global_index in enumerate(element_ids):
        element_nodes = mesh.elems[global_index]
        if sampling == "cg1":
            fiber_nodes = mesh.fiber_node[element_nodes]
            sheet_nodes = mesh.sheet_node[element_nodes]
        for gauss_point in range(n_gp):
            if sampling == "cg1":
                fiber = shape[gauss_point] @ fiber_nodes
                sheet = shape[gauss_point] @ sheet_nodes
            else:
                fiber, sheet, _normal = geom.sample_structural_frame(
                    mesh,
                    element_nodes,
                    shape[gauss_point],
                    tbar,
                    asign=asign,
                )
            ff, ss, fssym = struct_tensors(fiber, sheet)
            base = gauss_point * per_gp
            for name, values in (("ff", ff), ("ss", ss), ("fssym", fssym)):
                offset = base + offsets[name]
                element.svars[
                    local_index, offset : offset + 9
                ] = values.ravel()
    element.svars_trial = element.svars.copy()


def _rank_build_directory(requested, rank):
    """Return an isolated build directory and its optional temporary owner."""
    if requested is None:
        owner = tempfile.TemporaryDirectory(
            prefix=f"coupfe-cardiac-step0-mpi-rank-{rank:05d}-"
        )
        atexit.register(owner.cleanup)
        return Path(owner.name), owner
    root = Path(requested).expanduser()
    rank_directory = root / f"rank-{rank:05d}"
    rank_directory.mkdir(parents=True, exist_ok=True)
    return rank_directory, None


def _build_rank_local_batch(
    mesh,
    element_ids,
    coordinates,
    dt,
    build_directory,
    rank,
    *,
    formulation,
    element_evaluation,
    fiber_sampling,
    fiber_tbar,
    fiber_asign,
    integrator="be",
    material_parameters=None,
):
    """Compile one rank's selected standard-material or Q1/P0 batch."""
    if formulation not in LOCAL_PRESSURE_FORMULATIONS | {"std-kappa"}:
        raise ValueError(f"unsupported MPI formulation {formulation!r}")
    if integrator == "be":
        problem = CardiacHex8()
        kernel_builder = build_kernel
        pressure_stage_alpha_f = None
    elif integrator == "generalized-alpha":
        problem = CardiacHex8GeneralizedAlpha()
        kernel_builder = build_generalized_alpha_kernel
        pressure_stage_alpha_f = SOURCE_MATCHED_GENERALIZED_ALPHA.alpha_f
    else:
        raise ValueError(f"unsupported MPI integrator {integrator!r}")
    properties = np.asarray(problem._mat.props_array, dtype=float).copy()
    property_names = list(problem._mat.props.keys())
    if material_parameters is not None:
        selected = dict(material_parameters)
        if set(selected) != set(property_names):
            raise RuntimeError(
                "selected material parameters do not match the generated kernel"
            )
        properties[:] = [float(selected[name]) for name in property_names]
    kappa_index = property_names.index("kappa")
    eta_index = property_names.index("eta")
    tension_index = property_names.index("Ta")
    bulk_modulus = float(properties[kappa_index])
    if properties[eta_index] != 100.0:
        raise RuntimeError("the benchmark material must retain eta=100 Pa s")
    if formulation in LOCAL_PRESSURE_FORMULATIONS:
        properties[kappa_index] = 0.0
    schema = problem._mat.state_schema
    per_gp = sum(entry["size"] for entry in schema.values())
    _source, module = kernel_builder(
        tmpdir=str(build_directory),
        formulation="standard",
        module_name=(
            "cardiac_hex8_"
            f"{formulation.replace('-', '_')}_"
            f"{integrator.replace('-', '_')}_mpi_rank_{rank:05d}"
        ),
    )
    element = CompiledElement(
        module,
        props=properties,
        dof_per_node=3,
        n_svars=per_gp * 8,
        mcrd=3,
        n_elem=len(element_ids),
        dt=dt,
        backend="native",
        state_schema=schema,
    )
    _init_fiber_svars_subset(
        element,
        mesh,
        schema,
        per_gp,
        element_ids,
        sampling=fiber_sampling,
        tbar=fiber_tbar,
        asign=fiber_asign,
    )
    if formulation in LOCAL_PRESSURE_FORMULATIONS:
        batch = FusedQ1P0Batch(
            element,
            coordinates,
            bulk_modulus=bulk_modulus,
            global_element_ids=element_ids,
            evaluation_mode=element_evaluation,
            pressure_stage_alpha_f=pressure_stage_alpha_f,
            time_integrator=(
                "generalized-alpha-source-matched"
                if integrator == "generalized-alpha"
                else "be"
            ),
            pressure_law=_local_pressure_law(formulation),
        )
    else:
        batch = RankLocalStandardMaterialBatch(
            element,
            coordinates,
            global_element_ids=element_ids,
            evaluation_mode=element_evaluation,
            time_integrator=(
                "generalized-alpha-source-matched"
                if integrator == "generalized-alpha"
                else "be"
            ),
        )
    return batch, element, tension_index, bulk_modulus


def _petsc_owned_row_range(
    PETSc,
    ndof,
    communicator,
    *,
    dof_per_node=3,
    node_aligned=False,
):
    """Query the ownership layout that the persistent solver will reuse."""
    block_size = dof_per_node if node_aligned else 1
    probe = PETSc.Vec().createMPI(
        ndof, bsize=block_size, comm=communicator
    )
    try:
        return tuple(int(value) for value in probe.getOwnershipRange())
    finally:
        probe.destroy()


def _collective_setup_error(mpi_comm, rank, error, context):
    record = None if error is None else (
        rank,
        f"{type(error).__name__}: {error}",
    )
    records = mpi_comm.allgather(record)
    failures = [entry for entry in records if entry is not None]
    if failures:
        failed_rank, message = min(failures, key=lambda entry: entry[0])
        raise RuntimeError(
            f"{context} failed on MPI rank {failed_rank}: {message}"
        )


def _gather_global_element_field(
    mpi_comm, rank, n_element, element_ids, local_values, trailing_shape
):
    """Gather a rank-local element field and restore global element order on root."""
    local_values = np.asarray(local_values, dtype=float)
    expected = (len(element_ids), *trailing_shape)
    if local_values.shape != expected or not np.all(np.isfinite(local_values)):
        raise RuntimeError(
            f"rank-local element field has shape {local_values.shape}; expected {expected}"
        )
    pieces = mpi_comm.gather(
        (np.asarray(element_ids, dtype=np.int64), local_values), root=0
    )
    if rank != 0:
        return None
    result = np.empty((n_element, *trailing_shape), dtype=float)
    populated = np.zeros(n_element, dtype=bool)
    for identifiers, values in pieces:
        if np.any(populated[identifiers]):
            raise RuntimeError("MPI element partition contains duplicate global cells")
        result[identifiers] = values
        populated[identifiers] = True
    if not np.all(populated):
        missing = np.flatnonzero(~populated)
        raise RuntimeError(
            "MPI element partition did not cover every global cell; first missing "
            f"element={int(missing[0])}"
        )
    return result


def _validate_mass_partition_records(records, ndof):
    """Validate and summarize rank-ordered owned-row mass metadata."""
    if not isinstance(records, (list, tuple)) or not records:
        raise ValueError("mass partition records must be a nonempty sequence")
    expected_start = 0
    representation = records[0].get("representation")
    partition = records[0].get("partition")
    local_nnz = []
    touching_counts = []
    row_ranges = []
    for rank, record in enumerate(records):
        if not isinstance(record, dict):
            raise ValueError(f"mass partition record {rank} is not a mapping")
        if record.get("representation") != representation or record.get(
            "partition"
        ) != partition:
            raise ValueError("MPI ranks disagree on mass representation")
        start = record.get("row_start")
        stop = record.get("row_end")
        columns = record.get("global_columns")
        nnz = record.get("local_nnz")
        touching = record.get("touching_elements")
        numeric = (start, stop, columns, nnz, touching)
        if any(
            isinstance(value, bool) or not isinstance(value, (int, np.integer))
            for value in numeric
        ):
            raise ValueError("mass partition metadata must contain integers")
        start, stop, columns, nnz, touching = map(int, numeric)
        if start != expected_start or stop < start or columns != ndof:
            raise ValueError(
                "owned mass row ranges must form one contiguous global partition"
            )
        if nnz < 0 or touching < 0:
            raise ValueError("mass partition counts must be nonnegative")
        expected_start = stop
        row_ranges.append((start, stop))
        local_nnz.append(nnz)
        touching_counts.append(touching)
    if expected_start != ndof:
        raise ValueError("owned mass row ranges do not cover the global system")
    return {
        "representation": representation,
        "partition": partition,
        "row_ranges": row_ranges,
        "local_nnz": local_nnz,
        "touching_elements": touching_counts,
    }


def _mpi_load_histories(
    case,
    dt,
    t_end,
    load_horizon=None,
    *,
    integrator="be",
    benchmark_step=0,
):
    """Return endpoint labels, applied loads, and their evaluation times.

    Backward Euler delegates unchanged to the retained serial schedule.  The
    source-matched generalized-alpha path follows the local Simula driver: the
    first dynamic endpoint is ``dt``, but activation/pressure are sampled on a
    grid shifted back by ``alpha_f*dt``.  The ODE initial value remains at
    ``t=0`` and its integration span remains the fixed load horizon, exactly as
    in ``benchmark1.py``.  This matters for the adaptive Radau trajectory even
    when the physical load is initially zero.
    """
    if integrator == "be":
        times, tau, pressure, horizon = (
            serial_driver._benchmark_load_histories(
                case,
                dt,
                t_end,
                load_horizon,
                benchmark_step=benchmark_step,
            )
        )
        return times, tau, pressure, horizon, times.copy()
    if integrator != "generalized-alpha":
        raise ValueError(f"unsupported MPI integrator {integrator!r}")
    configuration = benchmark_configuration(benchmark_step, case)

    times, schedule_times, horizon = serial_driver._load_schedule_grid(
        dt, t_end, load_horizon
    )
    alpha_f = SOURCE_MATCHED_GENERALIZED_ALPHA.alpha_f
    shifted_schedule = schedule_times[1:] - alpha_f * float(dt)
    count = len(times)
    tau = np.zeros(count, dtype=float)
    pressure = np.zeros(count, dtype=float)
    evaluation_times = np.zeros(count, dtype=float)
    evaluation_times[1:] = shifted_schedule[: count - 1]
    if configuration.active_stress_enabled:
        tau[1:] = serial_driver.tau_of_t(
            shifted_schedule,
            p=configuration.activation_parameters,
            t_span=(0.0, horizon),
        )[: count - 1]
    if configuration.pressure_enabled:
        pressure[1:] = serial_driver.p_of_t(
            shifted_schedule,
            p=configuration.pressure_parameters,
            t_span=(0.0, horizon),
        )[: count - 1]
    validate_load_histories(configuration, times, tau, pressure)
    return times, tau, pressure, horizon, evaluation_times


def _parser():
    parser = argparse.ArgumentParser(
        description=(
            "Distributed Cardiac Benchmark 1 companion; launch with mpiexec. "
            "Use closed-multiblock with consistent mass, GP-direct fibers, "
            "and an explicit std-kappa or local-pressure formulation for "
            "Step 0 Case A/B or Step 2 Case B."
        )
    )
    parser.add_argument(
        "--case",
        choices=["A", "B"],
        default="A",
        help=(
            "case identity within --benchmark-step; Step 0 A is active-only, "
            "Step 0 B is pressure-only, and Step 2 B combines both"
        ),
    )
    parser.add_argument(
        "--benchmark-step",
        type=int,
        choices=[0, 2],
        default=0,
        help=(
            "Benchmark 1 protocol step; Step 2 currently supports Case B and "
            "applies active stress and ventricular pressure together"
        ),
    )
    parser.add_argument("--nt", type=int, default=2)
    parser.add_argument("--nmu", type=int, default=12)
    parser.add_argument("--ntheta", type=int, default=16)
    parser.add_argument(
        "--mesh-topology",
        choices=["polar-ring", "closed-multiblock"],
        default="closed-multiblock",
        help=(
            "polar-ring preserves the historical truncated companion; "
            "closed-multiblock selects the noncollapsed five-block domain "
            "(default)"
        ),
    )
    parser.add_argument("--ncore", type=int, default=20)
    parser.add_argument("--nradial", type=int, default=17)
    parser.add_argument("--core-half-width", type=float, default=0.36)
    parser.add_argument("--tip-refine", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=2.0e-3)
    parser.add_argument("--tend", type=float, default=1.0)
    parser.add_argument(
        "--integrator",
        choices=["be", "generalized-alpha"],
        default="be",
        help=(
            "be preserves the historical MPI path; generalized-alpha selects "
            "the fixed source-matched Step 0 Cases A/B / Step 2 Case B "
            "alpha_m=0.2, alpha_f=0.4 path"
        ),
    )
    parser.add_argument(
        "--load-horizon",
        type=float,
        default=None,
        metavar="SECONDS",
        help=(
            "integrate the benchmark load through this fixed horizon, then use "
            "the exact prefix through --tend; use 1.0 s for shortened rank gates"
        ),
    )
    parser.add_argument(
        "--raw-helix",
        action="store_true",
        help="use asign=+1 instead of the pinned-formula-matched default",
    )
    parser.add_argument(
        "--perturb",
        type=float,
        default=0.0,
        help="break mesh symmetry with N(0,sigma) node noise [m] (diagnostic)",
    )
    parser.add_argument(
        "--apex-offset",
        type=float,
        default=0.2,
        dest="apex_offset",
        help="non-degenerate open-apex offset in radians (must be positive)",
    )
    parser.add_argument(
        "--formulation",
        choices=["local-pressure", "local-pressure-paper", "std-kappa"],
        default="local-pressure",
        help=(
            "local-pressure preserves the historical log condensed Q1/P0 path; "
            "local-pressure-paper applies the benchmark paper volume law to "
            "the element mean dilatation; "
            "std-kappa retains the paper pointwise kappa material term"
        ),
    )
    parser.add_argument(
        "--mass",
        choices=["lumped", "consistent"],
        default="lumped",
        help=(
            "lumped preserves the historical MPI companion; consistent uses "
            "complete owned-row Q1-Hex8 CSR assembly"
        ),
    )
    parser.add_argument(
        "--material-eta",
        type=float,
        default=100.0,
        metavar="PA_S",
        help="paper viscosity; only eta=100 Pa s is accepted",
    )
    parser.add_argument(
        "--tbar-laplace",
        default=None,
        metavar="PATH.npy",
        help="validated per-node Laplace transmural field and sibling metadata",
    )
    parser.add_argument(
        "--fiber-sampling",
        choices=["cg1", "gp-direct"],
        default="cg1",
        help=(
            "cg1 preserves interpolated nodal frames; gp-direct evaluates the "
            "fiber rule at each Hex8 Gauss point"
        ),
    )
    parser.add_argument(
        "--out",
        default=None,
        help="rank-zero completed archive (distinct from the serial default)",
    )
    parser.add_argument(
        "--build-dir",
        help="persistent build root; a unique rank-NNNNN child is always used",
    )
    parser.add_argument(
        "--element-evaluation",
        choices=["joint", "split"],
        default="joint",
        help=(
            "joint cached R/K evaluation (default), or explicit split "
            "residual/tangent work"
        ),
    )
    parser.add_argument(
        "--linear-solver-profile",
        choices=[
            DIRECT_SUPERLU_DIST_PROFILE,
            FGMRES_GAMG_RIGID_PROFILE,
            FGMRES_GAMG_RIGID_REBUILD_PROFILE,
            FGMRES_ASM_LU_PROFILE,
            FGMRES_ASM_ILU1_PROFILE,
        ],
        default=DIRECT_SUPERLU_DIST_PROFILE,
        help=(
            "explicit PETSc linear-solver/preconditioner profile; iterative "
            "profiles preserve the fixed SNES and independent acceptance gate"
        ),
    )
    parser.add_argument(
        "--epicardial-normal-mode",
        choices=["facet", "nodal-smoothed"],
        default="facet",
        dest="epicardial_normal_mode",
        help=(
            "epicardial Robin normal field: 'facet' (default, historical per-Gauss "
            "facet normal, bit-identical) or 'nodal-smoothed' (area-weighted "
            "Q1-interpolated nodal normal; a mechanism diagnostic, not a "
            "benchmark-faithful operator)"
        ),
    )
    return parser


def _validate_arguments(parser, args):
    """Fail before loading MPI for invalid or parameter-changing inputs."""
    try:
        args.benchmark_configuration = benchmark_configuration(
            args.benchmark_step, args.case
        )
    except ValueError as error:
        parser.error(str(error))
    if args.nt < 1:
        parser.error("--nt must be positive")
    if args.nmu < 1 or args.ntheta < 3:
        parser.error("--nmu must be positive and --ntheta must be at least 3")
    if args.ncore < 4 or args.ncore % 2:
        parser.error("--ncore must be an even integer >= 4")
    if args.nradial < 1:
        parser.error("--nradial must be positive")
    if not 0.1 < args.core_half_width < 1.0 / np.sqrt(2.0):
        parser.error("--core-half-width must lie in (0.1, 1/sqrt(2))")
    if not np.isfinite(args.tip_refine) or not 1.0 <= args.tip_refine <= 8.0:
        parser.error("--tip-refine must lie in [1, 8]")
    if not np.isfinite(args.apex_offset) or args.apex_offset < 0.0:
        parser.error("--apex-offset must be finite and nonnegative")
    if (
        args.mesh_topology == "polar-ring"
        and args.formulation in LOCAL_PRESSURE_FORMULATIONS
        and args.apex_offset <= 0.0
    ):
        parser.error(
            "a local-pressure formulation on the polar-ring mesh requires "
            "--apex-offset > 0"
        )
    if not np.isfinite(args.perturb) or args.perturb < 0.0:
        parser.error("--perturb must be finite and nonnegative")
    if not np.isfinite(args.material_eta) or args.material_eta != 100.0:
        parser.error(
            "the benchmark physical parameters are fixed: --material-eta "
            "must equal the paper value 100 Pa s"
        )
    if args.integrator == "generalized-alpha":
        if args.load_horizon is None:
            args.load_horizon = 1.0
        elif args.load_horizon != 1.0:
            parser.error(
                "source-matched generalized-alpha requires the canonical "
                "1 s load-integration horizon"
            )
    # Validate both exact grids before importing PETSc or compiling kernels.
    try:
        serial_driver._load_schedule_grid(
            args.dt, args.tend, args.load_horizon
        )
    except ValueError as error:
        parser.error(str(error))

    historical_contract = (
        args.integrator == "be"
        and args.mesh_topology == "polar-ring"
        and args.formulation == "local-pressure"
        and args.mass == "lumped"
        and args.fiber_sampling == "cg1"
        and args.tbar_laplace is None
    )
    closed_contract = (
        args.case in {"A", "B"}
        and args.mesh_topology == "closed-multiblock"
        and args.formulation in LOCAL_PRESSURE_FORMULATIONS | {"std-kappa"}
        and args.mass == "consistent"
        and args.fiber_sampling == "gp-direct"
        and args.tbar_laplace is not None
        and not args.raw_helix
        and args.perturb == 0.0
    )
    if historical_contract:
        args.mpi_implementation = (
            DistributedPetscSnesSolver.HISTORICAL_Q1P0_IMPLEMENTATION
        )
    elif closed_contract:
        generalized_alpha_case_supported = (
            (args.benchmark_step == 0 and args.case in {"A", "B"})
            or (args.benchmark_step == 2 and args.case == "B")
        )
        if (
            args.integrator == "generalized-alpha"
            and not generalized_alpha_case_supported
        ):
            parser.error(
                "the source-matched generalized-alpha MPI path is currently "
                "validated for Step 0 Cases A/B and Step 2 Case B"
            )
        # These implementation identifiers retain the Step 0 provenance of
        # the numerical discretization. They do not label the physical load
        # case: benchmark_step/configuration_id in the archive are authoritative.
        args.mpi_implementation = {
            ("be", "std-kappa"): (
                DistributedPetscSnesSolver.CLOSED_STD_KAPPA_IMPLEMENTATION
            ),
            ("be", "local-pressure"): (
                DistributedPetscSnesSolver.CLOSED_LOCAL_PRESSURE_IMPLEMENTATION
            ),
            ("be", "local-pressure-paper"): (
                DistributedPetscSnesSolver.
                CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION
            ),
            ("generalized-alpha", "std-kappa"): (
                DistributedPetscSnesSolver.
                CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION
            ),
            ("generalized-alpha", "local-pressure"): (
                DistributedPetscSnesSolver.
                CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION
            ),
            ("generalized-alpha", "local-pressure-paper"): (
                DistributedPetscSnesSolver.
                CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION
            ),
        }[(args.integrator, args.formulation)]
    else:
        parser.error(
            "the MPI driver accepts either the historical polar-ring/"
            "local-pressure/lumped/CG1 contract or the closed Step 0/"
            "{std-kappa,local-pressure,local-pressure-paper}/consistent/"
            "GP-direct/Laplace-tbar "
            "contract; mixed discretizations are not labeled as an "
            "implementation"
        )


def _benchmark_reproduction_profile(args):
    if args.benchmark_step != 2:
        return "not-applicable"
    canonical = (
        args.case == "B"
        and args.integrator == "generalized-alpha"
        and args.formulation == "std-kappa"
        and args.dt == 0.001
        and args.load_horizon == 1.0
        and args.mesh_topology == "closed-multiblock"
        and args.mass == "consistent"
        and args.fiber_sampling == "gp-direct"
        and args.tbar_laplace is not None
        and not args.raw_helix
        and args.perturb == 0.0
    )
    if not canonical:
        return "diagnostic-noncanonical"
    return (
        "paper-source-matched-full-cycle"
        if args.tend == 1.0
        else "paper-source-matched-prefix"
    )


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    _validate_arguments(parser, args)
    configuration = args.benchmark_configuration
    reproduction_profile = _benchmark_reproduction_profile(args)
    if args.out is None:
        formulation_tag = args.formulation.replace("-", "_")
        args.out = (
            f"benchmark1_step{args.benchmark_step}_case{args.case}_mpi_"
            f"{formulation_tag}_full.npz"
            if args.benchmark_step != 0
            else f"case{args.case}_mpi_{formulation_tag}_full.npz"
        )
    times, tau, pressure_history, load_horizon, load_evaluation_times = (
        _mpi_load_histories(
            args.case,
            args.dt,
            args.tend,
            args.load_horizon,
            integrator=args.integrator,
            benchmark_step=args.benchmark_step,
        )
    )

    PETSc, MPI = _load_parallel_runtime()
    petsc_comm = PETSc.COMM_WORLD
    mpi_comm = MPI.COMM_WORLD
    rank = int(petsc_comm.getRank())
    size = int(petsc_comm.getSize())
    node_aligned_ownership = linear_solver_requires_node_aligned_ownership(
        args.linear_solver_profile
    )
    if mpi_comm.Get_rank() != rank or mpi_comm.Get_size() != size:
        raise RuntimeError("petsc4py and mpi4py COMM_WORLD communicators disagree")

    if args.mesh_topology == "closed-multiblock":
        mesh = geom.build_closed_mesh(
            args.nt,
            args.ncore,
            args.nradial,
            core_half_width=args.core_half_width,
            tip_refine=args.tip_refine,
            flip_helix=not args.raw_helix,
        )
    else:
        mesh = geom.build_mesh(
            args.nt,
            args.nmu,
            args.ntheta,
            flip_helix=not args.raw_helix,
            apex_offset=args.apex_offset,
        )

    fiber_tbar = np.asarray(mesh.param[:, 0], dtype=float).copy()
    tbar_definition = "analytic_parametric"
    tbar_source_filename = ""
    tbar_source_sha256 = ""
    tbar_metadata_filename = ""
    tbar_metadata_sha256 = ""
    tbar_metadata_schema = ""
    if args.tbar_laplace:
        tbar_path = Path(args.tbar_laplace).expanduser().resolve()
        if not tbar_path.is_file():
            parser.error(
                f"--tbar-laplace does not name a file: {args.tbar_laplace}"
            )
        try:
            tbar_provenance = serial_driver._load_laplace_tbar_metadata(
                mesh, tbar_path
            )
            fiber_tbar = serial_driver._apply_laplace_tbar(
                mesh,
                tbar_path,
                asign=-1.0 if not args.raw_helix else 1.0,
            )
        except (OSError, ValueError) as error:
            parser.error(str(error))
        tbar_definition = "laplace_presolved"
        tbar_source_filename = tbar_path.name
        tbar_source_sha256 = tbar_provenance["field_sha256"]
        tbar_metadata_filename = tbar_provenance["metadata_filename"]
        tbar_metadata_sha256 = tbar_provenance["metadata_sha256"]
        tbar_metadata_schema = tbar_provenance["metadata_schema"]

    if args.perturb > 0.0:
        rng = np.random.default_rng(0)
        mesh.nodes[:] = mesh.nodes + rng.normal(
            0.0, args.perturb, mesh.nodes.shape
        )
        if mesh.topology == "closed_multiblock_disk":
            geom.update_structural_frames(
                mesh,
                fiber_tbar,
                asign=-1.0 if not args.raw_helix else 1.0,
            )

    geometry_audit = None
    geometry_error = None
    if rank == 0:
        try:
            geometry_audit = audit_geometry(
                mesh,
                require_closed=args.mesh_topology == "closed-multiblock",
            )
        except Exception as error:
            geometry_error = error
    _collective_setup_error(
        mpi_comm, rank, geometry_error, "pre-solve geometry audit"
    )
    geometry_audit = mpi_comm.bcast(geometry_audit, root=0)
    pre_solve_audit = {"geometry": geometry_audit}

    view = KernelMeshView(mesh.nodes, mesh.elems, dof_per_node=3)
    partition = partition_elements(view, size)
    element_ids = np.flatnonzero(partition == rank).astype(np.int64)
    my_gm, my_coordinates, ndof = element_partition(view, rank, size)
    local_counts = np.asarray(mpi_comm.allgather(len(element_ids)), dtype=np.int64)
    if int(local_counts.sum()) != mesh.n_elem:
        raise RuntimeError("MPI element partition does not cover the global mesh")

    build_directory, _build_owner = _rank_build_directory(args.build_dir, rank)
    setup_error = None
    local_setup = None
    try:
        local_setup = _build_rank_local_batch(
            mesh,
            element_ids,
            my_coordinates,
            args.dt,
            build_directory,
            rank,
            formulation=args.formulation,
            element_evaluation=args.element_evaluation,
            fiber_sampling=args.fiber_sampling,
            fiber_tbar=fiber_tbar,
            fiber_asign=-1.0 if not args.raw_helix else 1.0,
            integrator=args.integrator,
            material_parameters=configuration.material_parameters,
        )
    except Exception as error:
        setup_error = error
    _collective_setup_error(
        mpi_comm, rank, setup_error, "rank-local cardiac kernel setup"
    )
    batch, element, tension_index, bulk_modulus = local_setup
    material_property_names = tuple(
        (
            CardiacHex8GeneralizedAlpha()
            if args.integrator == "generalized-alpha"
            else CardiacHex8()
        )._mat.props.keys()
    )

    mass_partition_summary = None
    if args.mass == "consistent":
        row_start, row_end = _petsc_owned_row_range(
            PETSc,
            ndof,
            petsc_comm,
            dof_per_node=3,
            node_aligned=node_aligned_ownership,
        )
        mass_error = None
        mass = None
        try:
            mass = assemble_owned_consistent_mass(
                mesh.nodes,
                mesh.elems,
                serial_driver.DENSITY,
                ndof,
                row_start,
                row_end,
                dof_per_node=3,
            )
        except Exception as error:
            mass_error = error
        _collective_setup_error(
            mpi_comm, rank, mass_error, "owned-row consistent-mass assembly"
        )
        mass_records = mpi_comm.allgather(mass.metadata())
        try:
            mass_partition_summary = _validate_mass_partition_records(
                mass_records, ndof
            )
        except Exception as error:
            mass_error = error
        _collective_setup_error(
            mpi_comm, rank, mass_error, "consistent-mass partition audit"
        )
        mass_representation = "consistent_q1_hex8"
    else:
        mass = serial_driver._hex_lumped_mass(
            mesh.nodes,
            mesh.elems,
            serial_driver.DENSITY,
            ndof,
            dof_per_node=3,
        )
        mass_representation = "lumped_row_sum"
        row_start, row_end = _petsc_owned_row_range(
            PETSc,
            ndof,
            petsc_comm,
            dof_per_node=3,
            node_aligned=node_aligned_ownership,
        )
        lumped_ranges = mpi_comm.allgather((row_start, row_end))
        expected_start = 0
        for start, stop in lumped_ranges:
            if start != expected_start or stop < start:
                raise RuntimeError(
                    "PETSc lumped-mass ownership is not a contiguous partition"
                )
            expected_start = stop
        if expected_start != ndof:
            raise RuntimeError(
                "PETSc lumped-mass ownership does not cover the global system"
            )
        mass_partition_summary = {
            "representation": mass_representation,
            "partition": "replicated-diagonal-owned-entry-insertion",
            "row_ranges": lumped_ranges,
            "local_nnz": mpi_comm.allgather(
                int(np.count_nonzero(mass[row_start:row_end]))
            ),
            "touching_elements": [0] * size,
        }

    _epi_mode = {
        "facet": "normal",
        "nodal-smoothed": "normal-smoothed",
    }[args.epicardial_normal_mode]
    robin = RobinOperator(
        mesh.nodes,
        ndof,
        [
            (
                mesh.facets_base,
                serial_driver.A_TOP,
                serial_driver.B_TOP,
                "full",
            ),
            (
                mesh.facets_epi,
                serial_driver.A_EPI,
                serial_driver.B_EPI,
                _epi_mode,
            ),
        ],
        dof_per_node=3,
    )
    boundary_error = None
    robin_audit = None
    if rank == 0:
        try:
            robin_audit = audit_robin(robin)
        except Exception as error:
            boundary_error = error
    _collective_setup_error(
        mpi_comm, rank, boundary_error, "pre-solve Robin audit"
    )
    pre_solve_audit["robin"] = mpi_comm.bcast(robin_audit, root=0)

    pressure = None
    if configuration.pressure_enabled:
        interior = mesh.nodes[mesh.elems[mesh.facets_endo_elem]].mean(axis=1)
        pressure = FollowerPressureOperator(
            mesh.nodes,
            ndof,
            mesh.facets_endo,
            dof_per_node=3,
            interior=interior,
        )
        pressure_audit = None
        boundary_error = None
        if rank == 0:
            try:
                pressure_audit = audit_pressure(mesh, dof_per_node=3)
            except Exception as error:
                boundary_error = error
        _collective_setup_error(
            mpi_comm, rank, boundary_error, "pre-solve pressure audit"
        )
        pre_solve_audit["pressure"] = mpi_comm.bcast(
            pressure_audit, root=0
        )

    if args.mesh_topology == "polar-ring" and args.apex_offset > 0.0:
        print(
            "WARNING: polar-ring with a positive apex offset selects the "
            "historical truncated open-tip domain, which is not the Benchmark 1 "
            "geometry; retained results from it are archive records, not "
            "current benchmark evidence.",
            flush=True,
        )

    if rank == 0:
        print(
            f"MPI mesh: {mesh.n_node} nodes, {mesh.n_elem} hexes, ndof={ndof}, "
            f"ranks={size}",
            flush=True,
        )
        print(
            f"BENCHMARK 1 STEP {args.benchmark_step} CASE {args.case} "
            f"({configuration.load_contract}): "
            f"tau_peak={tau.max():.1f}Pa  "
            f"p_peak={pressure_history.max():.1f}Pa",
            flush=True,
        )
        if args.mesh_topology == "closed-multiblock":
            print(
                "geometry: closed square-to-disk topology; no apex cut surface",
                flush=True,
            )
        print(
            "pre-solve boundary audit: "
            + json.dumps(pre_solve_audit, sort_keys=True),
            flush=True,
        )
        if args.formulation in LOCAL_PRESSURE_FORMULATIONS:
            print(
                "formulation: rank-local standard Hex8 material (kappa=0) + "
                "exactly condensed Q1/P0 pressure; law="
                f"{_local_pressure_law(args.formulation)}",
                flush=True,
            )
        else:
            print(
                "formulation: rank-local standard Hex8 material with pointwise "
                f"kappa={HO_PARAMS['kappa']:.3e} Pa",
                flush=True,
            )
        if args.mass == "consistent":
            print(
                "mass: consistent Q1-Hex8, complete PETSc-owned CSR rows",
                flush=True,
            )
        else:
            print("mass: lumped (row-summed historical companion)", flush=True)
        print(
            "fiber sampling: "
            + {
                "cg1": "cg1_gram_schmidt",
                "gp-direct": "gp_direct_rule",
            }[args.fiber_sampling],
            flush=True,
        )
        print(
            f"material viscosity: eta={args.material_eta:.9g} Pa s "
            "(paper parameter fixed)",
            flush=True,
        )
        print(
            "element evaluation: "
            f"{args.element_evaluation} "
            f"(compiled material residual-only available="
            f"{batch.material_residual_only_available})",
            flush=True,
        )
        print(
            "nonlinear solver: persistent PETSc COMM_WORLD SNES newtonls + bt, "
            f"linear profile={args.linear_solver_profile}",
            flush=True,
        )
        if args.integrator == "generalized-alpha":
            print(
                "time integration: source-matched generalized-alpha "
                "alpha_m=0.2 alpha_f=0.4 gamma=0.7 beta=0.36; "
                "material/Robin/active load at 1-alpha_f, inertia at "
                "1-alpha_m",
                flush=True,
            )
        else:
            print("time integration: backward Euler", flush=True)

    p0_location = p1_location = None
    history_u0 = history_u1 = None
    if rank == 0:
        p0_location = locate_hex8_point(mesh.nodes, mesh.elems, geom.P0)
        p1_location = locate_hex8_point(mesh.nodes, mesh.elems, geom.P1)
        history_u0 = np.zeros((len(times), 3), dtype=float)
        history_u1 = np.zeros((len(times), 3), dtype=float)
        print(
            "Hex8 output sampling: "
            f"p0 elem={p0_location.element_index} "
            f"|recon|={p0_location.reconstruction_error * 1e3:.3e}mm; "
            f"p1 elem={p1_location.element_index} "
            f"|recon|={p1_location.reconstruction_error * 1e3:.3e}mm",
            flush=True,
        )

    load_history = validate_load_histories(
        configuration, times, tau, pressure_history
    )
    peak_index = int(np.argmax(np.abs(load_history)))
    peak_displacement = np.zeros(ndof, dtype=float)
    peak_pressure_displacement = np.zeros(ndof, dtype=float)
    accepted_displacement = np.zeros(ndof, dtype=float)
    completed_steps = 0
    step_diagnostics = []
    solver = None
    solver_configuration = None
    started = time.time()
    try:
        solver = DistributedPetscSnesSolver(
            ndof,
            my_gm,
            my_coordinates,
            batch,
            mass,
            dof_per_node=3,
            robin=robin,
            pressure=pressure,
            linear_solver_profile=args.linear_solver_profile,
            node_coordinates=mesh.nodes,
            implementation=args.mpi_implementation,
            integrator=args.integrator,
        )
        for step in range(1, len(times)):
            previous_displacement = accepted_displacement
            element.props[tension_index] = float(tau[step])
            if pressure is not None:
                pressure.p = float(pressure_history[step])
            try:
                displacement, diagnostics = solver.solve_step(
                    t=float(times[step]), dt=float(args.dt)
                )
            except DistributedSnesSolveError as error:
                failure_record = _sanitize_failure_json(
                    {
                        "schema": "coupfe-cardiac-mpi-failure-v1",
                        "step": step,
                        "expected_steps": len(times) - 1,
                        "time_s": times[step],
                        "dt_s": args.dt,
                        "load_pa": load_history[step],
                        "active_tension_pa": tau[step],
                        "pressure_pa": pressure_history[step],
                        "completed_steps": completed_steps,
                        "linear_solver_profile": args.linear_solver_profile,
                        "error": str(error),
                        "solver_diagnostics": (
                            None
                            if error.diagnostics is None
                            else error.diagnostics.as_dict()
                        ),
                    }
                )
                if rank == 0:
                    print(
                        "FAILED before state commit: "
                        + json.dumps(
                            failure_record,
                            sort_keys=True,
                            separators=(",", ":"),
                            allow_nan=False,
                        ),
                        file=sys.stderr,
                        flush=True,
                    )
                raise
            diagnostic_record = diagnostics.as_dict()
            diagnostic_record.update(
                {
                    "active_tension_pa": float(tau[step]),
                    "pressure_pa": float(pressure_history[step]),
                }
            )
            # Run on every rank so a serialization defect cannot become a
            # rank-zero-only failure followed by a collective deadlock.
            json.dumps(diagnostic_record, allow_nan=False)
            step_diagnostics.append(diagnostic_record)
            completed_steps = step
            accepted_displacement = displacement.copy()
            if rank == 0:
                history_u0[step] = interpolate_displacement(
                    displacement, p0_location, dof_per_node=3
                )
                history_u1[step] = interpolate_displacement(
                    displacement, p1_location, dof_per_node=3
                )
            if step == peak_index:
                peak_displacement = displacement.copy()
                if args.integrator == "generalized-alpha":
                    peak_pressure_displacement = (
                        SOURCE_MATCHED_GENERALIZED_ALPHA.force_stage(
                            previous_displacement, displacement
                        )
                    )
                else:
                    peak_pressure_displacement = displacement.copy()
            if rank == 0 and (
                step % max(1, len(times) // 25) == 0 or step < 3
            ):
                maximum = float(
                    np.max(np.abs(displacement.reshape(-1, 3)))
                )
                if configuration.load_contract == "active-stress-plus-pressure":
                    load_label = (
                        f"tau={tau[step]:9.1f}Pa "
                        f"p={pressure_history[step]:8.1f}Pa"
                    )
                else:
                    load_label = f"load={load_history[step]:9.1f}Pa"
                print(
                    f"  step {step:4d} t={times[step]:.3f}s  "
                    f"{load_label}  "
                    f"nit={diagnostics.nonlinear_iterations:2d}  "
                    f"|u|max={maximum * 1e3:7.3f}mm  "
                    f"u0=({history_u0[step, 0] * 1e3:+.2f},"
                    f"{history_u0[step, 1] * 1e3:+.2f},"
                    f"{history_u0[step, 2] * 1e3:+.2f})mm",
                    flush=True,
                )
        solver_configuration = solver.configuration()
    finally:
        if solver is not None:
            solver.close()

    if completed_steps != len(times) - 1:
        # Normally unreachable because a failed step propagates above.  Keep a
        # second explicit guard immediately before all result-field gathers.
        raise RuntimeError(
            f"incomplete MPI solve: completed {completed_steps}/{len(times) - 1}"
        )

    local_peak = peak_displacement[my_gm]
    local_pressure_peak = peak_pressure_displacement[my_gm]
    field_error = None
    local_det_f = local_pressure = None
    try:
        local_det_f = batch.deformation_jacobians(local_peak)
        if args.formulation in LOCAL_PRESSURE_FORMULATIONS:
            local_pressure = batch.element_pressure(local_pressure_peak)
    except Exception as error:
        field_error = error
    _collective_setup_error(
        mpi_comm, rank, field_error, "accepted peak-field reconstruction"
    )
    det_f_peak = _gather_global_element_field(
        mpi_comm,
        rank,
        mesh.n_elem,
        element_ids,
        local_det_f,
        (8,),
    )
    if args.formulation in LOCAL_PRESSURE_FORMULATIONS:
        element_pressure_peak = _gather_global_element_field(
            mpi_comm,
            rank,
            mesh.n_elem,
            element_ids,
            local_pressure,
            (),
        )
    else:
        element_pressure_peak = (
            np.empty(0, dtype=float) if rank == 0 else None
        )

    if rank != 0:
        return None

    print(f"elapsed {time.time() - started:.1f}s", flush=True)
    solver_configuration_json = json.dumps(
        solver_configuration,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    step_diagnostics_json = json.dumps(
        step_diagnostics,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    runtime_metadata = serial_driver._runtime_metadata()
    runtime_metadata.update(
        {
            "driver": "examples/cardiac_benchmark/run_mpi.py",
            "mpi_enabled": True,
            "mpi_ranks": size,
            "mpi_implementation": args.mpi_implementation,
        }
    )
    formulation_label = {
        "local-pressure": "hex8_local_pressure_p0_condensed_logj",
        "local-pressure-paper": (
            "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
        ),
        "std-kappa": "hex8_standard_pointwise_kappa",
    }[args.formulation]
    fiber_sampling_label = {
        "cg1": "cg1_gram_schmidt",
        "gp-direct": "gp_direct_rule",
    }[args.fiber_sampling]
    validate_runtime_material_properties(
        configuration,
        material_property_names,
        element.props,
        condensed_local_pressure=args.formulation in LOCAL_PRESSURE_FORMULATIONS,
        active_tension_pa=tau[-1],
    )
    benchmark_archive_metadata = benchmark_metadata(
        configuration,
        material_parameters=configuration.material_parameters,
        activation_parameters=configuration.activation_parameters,
        pressure_parameters=configuration.pressure_parameters,
    )
    output = save_completed(
        args.out,
        completed_steps=completed_steps,
        expected_steps=len(times) - 1,
        **runtime_metadata,
        **benchmark_archive_metadata,
        benchmark_reproduction_profile=reproduction_profile,
        times=times,
        tau=tau,
        pres=pressure_history,
        u0=history_u0,
        u1=history_u1,
        p0=geom.P0,
        p1=geom.P1,
        U_peak=peak_displacement,
        n_peak=peak_index,
        nodes=mesh.nodes,
        elems=mesh.elems,
        fiber=mesh.fiber,
        facets_endo=mesh.facets_endo,
        case=args.case,
        dt=args.dt,
        integrator=args.integrator,
        generalized_alpha_alpha_m=(
            SOURCE_MATCHED_GENERALIZED_ALPHA.alpha_m
            if args.integrator == "generalized-alpha"
            else 0.0
        ),
        generalized_alpha_alpha_f=(
            SOURCE_MATCHED_GENERALIZED_ALPHA.alpha_f
            if args.integrator == "generalized-alpha"
            else 0.0
        ),
        generalized_alpha_gamma=(
            SOURCE_MATCHED_GENERALIZED_ALPHA.gamma
            if args.integrator == "generalized-alpha"
            else 0.0
        ),
        generalized_alpha_beta=(
            SOURCE_MATCHED_GENERALIZED_ALPHA.beta
            if args.integrator == "generalized-alpha"
            else 0.0
        ),
        generalized_alpha_stage_contract=(
            "simula-source-matched-v1"
            if args.integrator == "generalized-alpha"
            else "not-applicable"
        ),
        load_evaluation_times_s=load_evaluation_times,
        formulation=formulation_label,
        density=serial_driver.DENSITY,
        material_kernel_formulation="standard",
        material_model_id=MATERIAL_MODEL_ID,
        material_kappa_pa=(
            0.0
            if args.formulation in LOCAL_PRESSURE_FORMULATIONS
            else bulk_modulus
        ),
        local_pressure_bulk_modulus_pa=(
            bulk_modulus
            if args.formulation in LOCAL_PRESSURE_FORMULATIONS
            else 0.0
        ),
        local_pressure_volume_law=(
            LOCAL_PRESSURE_ARCHIVE_LAWS[args.formulation]
            if args.formulation in LOCAL_PRESSURE_FORMULATIONS
            else "not-applicable"
        ),
        mass_representation=mass_representation,
        det_f_gauss_peak=det_f_peak,
        element_pressure_peak_pa=element_pressure_peak,
        element_pressure_peak_stage=(
            "alpha_f_force_stage"
            if args.integrator == "generalized-alpha"
            else "endpoint"
        ),
        nonlinear_solver="petsc-snes-mpi",
        solver_configuration_json=solver_configuration_json,
        nonlinear_step_diagnostics_json=step_diagnostics_json,
        mpi_world_size=size,
        mpi_local_element_counts=local_counts,
        mpi_partition="coupfe.partition_elements",
        mpi_build_layout="isolated-rank-directories",
        mpi_factor_solver_type=(
            "superlu_dist"
            if args.linear_solver_profile == DIRECT_SUPERLU_DIST_PROFILE
            else "none"
        ),
        mpi_linear_solver_profile=args.linear_solver_profile,
        mpi_mass_partition=mass_partition_summary["partition"],
        mpi_mass_owned_row_ranges=np.asarray(
            mass_partition_summary["row_ranges"], dtype=np.int64
        ),
        mpi_mass_local_nnz=np.asarray(
            mass_partition_summary["local_nnz"], dtype=np.int64
        ),
        mpi_mass_touching_element_counts=np.asarray(
            mass_partition_summary["touching_elements"], dtype=np.int64
        ),
        element_evaluation_mode=args.element_evaluation,
        compiled_material_residual_only_available=(
            batch.material_residual_only_available
        ),
        epicardial_normal_mode=args.epicardial_normal_mode,
        pre_solve_audit_json=json.dumps(
            pre_solve_audit,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ),
        fiber_sampling=fiber_sampling_label,
        fiber_sampling_option=args.fiber_sampling,
        fiber_direction_reconstruction=geom.structural_direction_reconstruction(
            mesh.topology
        ),
        tbar_definition=tbar_definition,
        tbar_source_filename=tbar_source_filename,
        tbar_source_sha256=tbar_source_sha256,
        tbar_metadata_filename=tbar_metadata_filename,
        tbar_metadata_sha256=tbar_metadata_sha256,
        tbar_metadata_schema=tbar_metadata_schema,
        isotropic=False,
        material_eta_pa_s=args.material_eta,
        viscous_term_active=bool(args.material_eta > 0.0),
        parameter_variant="benchmark_eta",
        point_sampling="hex8_reference_isoparametric",
        p0_sampling_element=p0_location.element_index,
        p0_sampling_natural=np.asarray(p0_location.natural_coordinates),
        p0_sampling_weights=np.asarray(p0_location.weights),
        p0_sampling_reconstruction_error_m=p0_location.reconstruction_error,
        p1_sampling_element=p1_location.element_index,
        p1_sampling_natural=np.asarray(p1_location.natural_coordinates),
        p1_sampling_weights=np.asarray(p1_location.weights),
        p1_sampling_reconstruction_error_m=p1_location.reconstruction_error,
        viscous_rate=(
            "velocity_consistent_green_lagrange_at_alpha_f_stage"
            if args.integrator == "generalized-alpha"
            else "backward_difference"
        ),
        mesh_topology=mesh.topology,
        n_side=mesh.n_side,
        n_core=mesh.n_core,
        n_radial=mesh.n_radial,
        core_half_width=mesh.core_half_width,
        tip_refine=mesh.tip_refine,
        apex_offset=(
            args.apex_offset if args.mesh_topology == "polar-ring" else 0.0
        ),
        flip_helix=not args.raw_helix,
        perturb=args.perturb,
        t_end=args.tend,
        load_horizon=load_horizon,
        a_top=serial_driver.A_TOP,
        b_top=serial_driver.B_TOP,
        a_epi=serial_driver.A_EPI,
        b_epi=serial_driver.B_EPI,
        n_t=mesh.n_t,
        n_mu=mesh.n_mu,
        n_theta=mesh.n_theta,
    )
    print(f"saved -> {output}", flush=True)
    return output


if __name__ == "__main__":
    main()
