from __future__ import annotations

import json
import importlib

import numpy as np
import pytest

import activation
import benchmark_parameters as benchmark
import post
import run as serial_driver
import run_mpi
from coupfe.operators.base import Residual, Tangent
from distributed_solver import DistributedPetscSnesSolver
from generalized_alpha import SOURCE_MATCHED_GENERALIZED_ALPHA


def test_step2_case_b_parameters_are_exact_and_immutable():
    configuration = benchmark.benchmark_configuration(2, "B")

    assert configuration.identity == (
        "benchmark-1-step-2-case-B-active-stress-plus-pressure"
    )
    assert configuration.load_contract == "active-stress-plus-pressure"
    assert configuration.active_stress_enabled is True
    assert configuration.pressure_enabled is True
    assert dict(configuration.material_parameters) == {
        "a": 295.0,
        "b": 8.023,
        "a_f": 92360.0,
        "b_f": 16.026,
        "a_s": 12405.0,
        "b_s": 11.12,
        "a_fs": 1080.0,
        "b_fs": 11.436,
        "kappa": 1.0e6,
        "k_sw": 100.0,
        "eta": 100.0,
        "Ta": 0.0,
    }
    assert configuration.activation_parameters["sigma_0"] == 100000.0
    assert dict(configuration.pressure_parameters) == activation.PRES_PARAMS
    with pytest.raises(TypeError):
        configuration.material_parameters["a"] = 1.0


def test_cardiac_modules_remain_importable_as_namespace_packages():
    packaged_activation = importlib.import_module(
        "examples.cardiac_benchmark.activation"
    )
    packaged_material = importlib.import_module(
        "examples.cardiac_benchmark.material"
    )
    assert packaged_activation.ACT_PARAMS["sigma_0"] == 150000.0
    assert packaged_material.HO_PARAMS["a"] == 59.0


def test_step0_parameter_and_load_identities_are_unchanged():
    step0_a = benchmark.benchmark_configuration(0, "A")
    step0_b = benchmark.benchmark_configuration(0, "B")

    assert dict(step0_a.material_parameters) == dict(step0_b.material_parameters)
    assert dict(step0_a.activation_parameters) == activation.ACT_PARAMS
    assert dict(step0_a.pressure_parameters) == activation.PRES_PARAMS
    assert step0_a.load_contract == "active-stress-only"
    assert step0_b.load_contract == "pressure-only"

    times_a, tau_a, pressure_a, _ = serial_driver._benchmark_load_histories(
        "A", 0.002, 0.32, 1.0
    )
    times_b, tau_b, pressure_b, _ = serial_driver._benchmark_load_histories(
        "B", 0.002, 0.32, 1.0
    )
    schedule = serial_driver._time_grid(0.002, 1.0)
    np.testing.assert_array_equal(times_a, times_b)
    np.testing.assert_array_equal(
        tau_a, activation.tau_of_t(schedule)[: len(times_a)]
    )
    np.testing.assert_array_equal(pressure_a, np.zeros_like(pressure_a))
    np.testing.assert_array_equal(tau_b, np.zeros_like(tau_b))
    np.testing.assert_array_equal(
        pressure_b, activation.p_of_t(schedule)[: len(times_b)]
    )


def test_step2_case_b_load_history_combines_active_stress_and_pressure():
    dt = 0.002
    times, tau, pressure, horizon = serial_driver._benchmark_load_histories(
        "B", dt, 0.32, 1.0, benchmark_step=2
    )
    configuration = benchmark.benchmark_configuration(2, "B")
    schedule = serial_driver._time_grid(dt, 1.0)

    expected_tau = activation.tau_of_t(
        schedule, p=configuration.activation_parameters
    )[: len(times)]
    expected_pressure = activation.p_of_t(
        schedule, p=configuration.pressure_parameters
    )[: len(times)]
    np.testing.assert_array_equal(tau, expected_tau)
    np.testing.assert_array_equal(pressure, expected_pressure)
    assert horizon == 1.0
    assert np.max(tau) > 0.0
    assert np.max(pressure) > 0.0
    np.testing.assert_array_equal(
        benchmark.validate_load_histories(configuration, times, tau, pressure),
        tau + pressure,
    )

    step0_tau = serial_driver._benchmark_load_histories(
        "A", dt, 0.32, 1.0
    )[1]
    # The ODE is linear in sigma_0; adaptive Radau tolerances make separately
    # integrated trajectories agree at solver tolerance rather than bitwise.
    np.testing.assert_allclose(
        tau, (2.0 / 3.0) * step0_tau, rtol=1.0e-4, atol=1.0e-10
    )
    np.testing.assert_array_equal(
        pressure,
        serial_driver._benchmark_load_histories("B", dt, 0.32, 1.0)[2],
    )


def test_step2_case_b_short_run_is_exact_prefix_of_fixed_load_horizon():
    short = serial_driver._benchmark_load_histories(
        "B", 0.002, 0.32, 1.0, benchmark_step=2
    )
    complete = serial_driver._benchmark_load_histories(
        "B", 0.002, 1.0, 1.0, benchmark_step=2
    )
    count = len(short[0])
    for short_history, complete_history in zip(short[:3], complete[:3]):
        np.testing.assert_array_equal(short_history, complete_history[:count])

    mpi = run_mpi._mpi_load_histories(
        "B", 0.002, 0.32, 1.0, benchmark_step=2, integrator="be"
    )
    for observed, expected in zip(mpi[:4], short):
        np.testing.assert_array_equal(observed, expected)
    np.testing.assert_array_equal(mpi[4], mpi[0])


def test_step2_case_b_generalized_alpha_uses_shifted_joint_loads():
    dt = 0.001
    times, tau, pressure, horizon, evaluation_times = run_mpi._mpi_load_histories(
        "B",
        dt,
        0.32,
        1.0,
        benchmark_step=2,
        integrator="generalized-alpha",
    )
    configuration = benchmark.benchmark_configuration(2, "B")
    np.testing.assert_array_equal(
        evaluation_times[1:], times[1:] - 0.4 * dt
    )
    expected_tau = np.zeros_like(times)
    expected_pressure = np.zeros_like(times)
    expected_tau[1:] = activation.tau_of_t(
        evaluation_times[1:],
        p=configuration.activation_parameters,
        t_span=(0.0, horizon),
    )
    expected_pressure[1:] = activation.p_of_t(
        evaluation_times[1:],
        p=configuration.pressure_parameters,
        t_span=(0.0, horizon),
    )
    np.testing.assert_array_equal(tau, expected_tau)
    np.testing.assert_array_equal(pressure, expected_pressure)
    assert np.max(tau) > 0.0
    assert np.max(pressure) > 0.0


def _closed_generalized_alpha_arguments(*, benchmark_step, case):
    return [
        "--benchmark-step",
        str(benchmark_step),
        "--case",
        case,
        "--integrator",
        "generalized-alpha",
        "--mesh-topology",
        "closed-multiblock",
        "--formulation",
        "std-kappa",
        "--mass",
        "consistent",
        "--fiber-sampling",
        "gp-direct",
        "--tbar-laplace",
        "closed-tbar.npy",
    ]


def test_generalized_alpha_cli_allows_step0_cases_and_step2_case_b():
    parser = run_mpi._parser()
    step2 = parser.parse_args(
        _closed_generalized_alpha_arguments(benchmark_step=2, case="B")
    )
    run_mpi._validate_arguments(parser, step2)
    assert step2.load_horizon == 1.0
    assert step2.mpi_implementation == (
        DistributedPetscSnesSolver.CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION
    )
    assert run_mpi._benchmark_reproduction_profile(step2) == (
        "diagnostic-noncanonical"
    )

    canonical = parser.parse_args(
        _closed_generalized_alpha_arguments(benchmark_step=2, case="B")
        + ["--dt", "0.001", "--tend", "0.32"]
    )
    run_mpi._validate_arguments(parser, canonical)
    assert run_mpi._benchmark_reproduction_profile(canonical) == (
        "paper-source-matched-prefix"
    )

    step0 = parser.parse_args(
        _closed_generalized_alpha_arguments(benchmark_step=0, case="B")
    )
    run_mpi._validate_arguments(parser, step0)
    assert step0.load_horizon == 1.0
    assert step0.mpi_implementation == (
        DistributedPetscSnesSolver.CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION
    )
    assert run_mpi._benchmark_reproduction_profile(step0) == "not-applicable"


@pytest.mark.parametrize("formulation", ["local-pressure", "local-pressure-paper"])
def test_step2_condensed_pressure_variants_are_not_labeled_paper_source_matched(
    formulation,
):
    parser = run_mpi._parser()
    arguments = _closed_generalized_alpha_arguments(benchmark_step=2, case="B")
    arguments[arguments.index("std-kappa")] = formulation
    args = parser.parse_args(arguments + ["--dt", "0.001", "--tend", "1.0"])

    run_mpi._validate_arguments(parser, args)

    assert run_mpi._benchmark_reproduction_profile(args) == (
        "diagnostic-noncanonical"
    )


def test_step0_case_b_generalized_alpha_uses_shifted_pressure_only_load():
    dt = 0.001
    times, tau, pressure, horizon, evaluation_times = (
        run_mpi._mpi_load_histories(
            "B",
            dt,
            0.32,
            1.0,
            benchmark_step=0,
            integrator="generalized-alpha",
        )
    )
    expected_pressure = np.zeros_like(times)
    expected_pressure[1:] = activation.p_of_t(
        evaluation_times[1:],
        p=benchmark.benchmark_configuration(0, "B").pressure_parameters,
        t_span=(0.0, horizon),
    )
    np.testing.assert_array_equal(
        evaluation_times[1:], times[1:] - 0.4 * dt
    )
    np.testing.assert_array_equal(tau, np.zeros_like(tau))
    np.testing.assert_array_equal(pressure, expected_pressure)


def test_generalized_alpha_follower_pressure_uses_force_stage_and_chain_rule():
    class RecordingPressure:
        def __init__(self):
            self.residual_call = None
            self.tangent_call = None

        def residual(self, displacement, _state, time, dt):
            self.residual_call = (np.asarray(displacement).copy(), time, dt)
            return Residual(np.array([0]), np.array([7.0]))

        def tangent(self, displacement, _state, time, dt):
            self.tangent_call = (np.asarray(displacement).copy(), time, dt)
            return Tangent(np.array([0]), np.array([0]), np.array([5.0]))

    solver = object.__new__(DistributedPetscSnesSolver)
    solver.generalized_alpha = SOURCE_MATCHED_GENERALIZED_ALPHA
    solver._active = {"t": 0.1, "stage_time": 0.096, "dt": 0.01}
    solver._u_prev = np.array([2.0])
    solver._v_prev = np.array([0.0])
    solver._a_prev = np.array([0.0])
    solver.robin = None
    solver.pressure = RecordingPressure()
    endpoint = np.array([10.0])

    _robin, residual = solver._boundary_residuals(endpoint)
    _robin_tangent, tangent = solver._boundary_tangents(endpoint)
    expected_stage = SOURCE_MATCHED_GENERALIZED_ALPHA.force_stage(
        solver._u_prev, endpoint
    )
    np.testing.assert_array_equal(solver.pressure.residual_call[0], expected_stage)
    np.testing.assert_array_equal(solver.pressure.tangent_call[0], expected_stage)
    assert solver.pressure.residual_call[1:] == (0.096, 0.01)
    assert solver.pressure.tangent_call[1:] == (0.096, 0.01)
    np.testing.assert_array_equal(residual.values, [7.0])
    np.testing.assert_allclose(
        tangent.values,
        5.0 * SOURCE_MATCHED_GENERALIZED_ALPHA.force_displacement_tangent(),
    )


def test_step2_case_b_cli_identity_fails_closed_for_other_cases(capsys):
    with pytest.raises(ValueError, match="only for Case B"):
        benchmark.benchmark_configuration(2, "A")

    with pytest.raises(SystemExit) as serial_error:
        serial_driver.main(["--benchmark-step", "2", "--case", "A"])
    assert serial_error.value.code == 2
    assert "only for Case B" in capsys.readouterr().err

    parser = run_mpi._parser()
    invalid = parser.parse_args(["--benchmark-step", "2", "--case", "A"])
    with pytest.raises(SystemExit) as mpi_error:
        run_mpi._validate_arguments(parser, invalid)
    assert mpi_error.value.code == 2

    valid = parser.parse_args(
        ["--benchmark-step", "2", "--case", "B",
         "--mesh-topology", "polar-ring"]
    )
    run_mpi._validate_arguments(parser, valid)
    assert valid.benchmark_configuration == benchmark.benchmark_configuration(2, "B")


def test_step2_case_b_metadata_is_exact_and_parameter_drift_fails_closed():
    configuration = benchmark.benchmark_configuration(2, "B")
    metadata = benchmark.benchmark_metadata(
        configuration,
        material_parameters=configuration.material_parameters,
        activation_parameters=configuration.activation_parameters,
        pressure_parameters=configuration.pressure_parameters,
    )

    assert metadata["benchmark_step"] == 2
    assert metadata["benchmark_load_contract"] == "active-stress-plus-pressure"
    assert metadata["benchmark_active_stress_enabled"] is True
    assert metadata["benchmark_pressure_enabled"] is True
    material = json.loads(metadata["benchmark_material_parameters_json"])
    activation_parameters = json.loads(
        metadata["benchmark_activation_parameters_json"]
    )
    pressure_parameters = json.loads(
        metadata["benchmark_pressure_parameters_json"]
    )
    assert material["a_f"] == 92360.0
    assert material["a_s"] == 12405.0
    assert material["a_fs"] == 1080.0
    assert activation_parameters["sigma_0"] == 100000.0
    assert pressure_parameters == dict(configuration.pressure_parameters)
    source_manifest = json.loads(
        metadata["benchmark_runtime_source_manifest_json"]
    )
    assert tuple(sorted(source_manifest)) == tuple(
        sorted(benchmark.RUNTIME_SOURCE_FILES)
    )
    assert len(metadata["benchmark_runtime_source_sha256"]) == 64

    wrong_material = dict(configuration.material_parameters)
    wrong_material["a"] = 59.0
    with pytest.raises(RuntimeError, match="material parameter 'a'"):
        benchmark.benchmark_metadata(
            configuration,
            material_parameters=wrong_material,
            activation_parameters=configuration.activation_parameters,
            pressure_parameters=configuration.pressure_parameters,
        )
    wrong_activation = dict(configuration.activation_parameters)
    wrong_activation["sigma_0"] = 150000.0
    with pytest.raises(RuntimeError, match="activation parameter 'sigma_0'"):
        benchmark.benchmark_metadata(
            configuration,
            material_parameters=configuration.material_parameters,
            activation_parameters=wrong_activation,
            pressure_parameters=configuration.pressure_parameters,
        )


def test_report_metadata_validates_step2_identity_and_exact_parameters(tmp_path):
    configuration = benchmark.benchmark_configuration(2, "B")
    metadata = benchmark.benchmark_metadata(
        configuration,
        material_parameters=configuration.material_parameters,
        activation_parameters=configuration.activation_parameters,
        pressure_parameters=configuration.pressure_parameters,
    )
    path = tmp_path / "step2-metadata.npz"
    np.savez(path, case="B", **metadata)
    with np.load(path, allow_pickle=False) as archive:
        validated = post._validate_benchmark_configuration_metadata(
            archive, path, "B"
        )
    assert validated["step"] == 2
    assert validated["load_contract"] == "active-stress-plus-pressure"
    assert validated["material_parameters"]["a"] == 295.0
    assert validated["activation_parameters"]["sigma_0"] == 100000.0
    assert validated["runtime_source_sha256"] == metadata[
        "benchmark_runtime_source_sha256"
    ]

    wrong = dict(metadata)
    activation_parameters = json.loads(
        wrong["benchmark_activation_parameters_json"]
    )
    activation_parameters["sigma_0"] = 150000.0
    wrong["benchmark_activation_parameters_json"] = json.dumps(
        activation_parameters
    )
    np.savez(path, case="B", **wrong)
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="inconsistent with its step/case"):
            post._validate_benchmark_configuration_metadata(archive, path, "B")

    wrong = dict(metadata)
    wrong["benchmark_runtime_source_sha256"] = "0" * 64
    np.savez(path, case="B", **wrong)
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="runtime_source_sha256"):
            post._validate_benchmark_configuration_metadata(archive, path, "B")


def test_legacy_reporter_routes_step2_to_the_blinded_comparator():
    with pytest.raises(ValueError, match="compare_step2b_case_b.py"):
        post.load_reference("step_2B", reference_dir="unused")

@pytest.mark.parametrize("condensed", [False, True])
def test_step2_case_b_compiled_property_identity_is_checked(condensed):
    configuration = benchmark.benchmark_configuration(2, "B")
    names = tuple(configuration.material_parameters)
    properties = np.asarray(
        [configuration.material_parameters[name] for name in names], dtype=float
    )
    if condensed:
        properties[names.index("kappa")] = 0.0
    properties[names.index("Ta")] = 1234.5
    benchmark.validate_runtime_material_properties(
        configuration,
        names,
        properties,
        condensed_local_pressure=condensed,
        active_tension_pa=1234.5,
    )

    properties[names.index("a_s")] = 2481.0
    with pytest.raises(RuntimeError, match="property 'a_s'"):
        benchmark.validate_runtime_material_properties(
            configuration,
            names,
            properties,
            condensed_local_pressure=condensed,
            active_tension_pa=1234.5,
        )


def test_load_contract_rejects_missing_half_of_step2_case_b():
    configuration = benchmark.benchmark_configuration(2, "B")
    times = np.array([0.0, 0.2, 0.3])
    tau = np.array([0.0, 1.0, 2.0])
    pressure = np.array([0.0, 3.0, 4.0])

    # A Step 2 run may be truncated before a load rises, so enabled histories
    # are allowed to be zero. The immutable flags and exact ODE provenance in
    # the metadata are what distinguish such a run; disabled Step 0 loads must
    # remain identically zero.
    benchmark.validate_load_histories(configuration, times, tau, pressure)
    with pytest.raises(RuntimeError, match="pressure is nonzero"):
        benchmark.validate_load_histories(
            benchmark.benchmark_configuration(0, "A"), times, tau, pressure
        )
    with pytest.raises(RuntimeError, match="active stress is nonzero"):
        benchmark.validate_load_histories(
            benchmark.benchmark_configuration(0, "B"), times, tau, pressure
        )
