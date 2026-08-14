"""Validate the releasable cardiac source tree, wheel, and source archive."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import stat
import subprocess
import tarfile
import zipfile
import xml.etree.ElementTree as ElementTree
from collections import Counter
from fnmatch import fnmatchcase
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Optional


# Current runtime/dependency Core release, anonymously reachable from public main.
APPROVED_PUBLIC_CORE_REF: Optional[str] = (
    "e2f42ed5772850a0a23a2ce434f430c287eae5c8"
)
# Immutable Core provenance recorded by the retained result and figure inputs.
HISTORICAL_RETAINED_CORE_REF = (
    "454f73ce2de284262b214a2b37bd676c6aca3c0a"
)
LICENSE_EXPRESSION = "Apache-2.0 AND CC-BY-4.0 AND MIT"
PUBLIC_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"
RETAINED_RESULT_APP_REF = "44cbfed9e09d4150203faae3087f2e4617d1fc47"
TRUNCATED_POLAR_ARCHIVE_DIRECTORY = (
    "examples/cardiac_benchmark/results/archive/truncated_polar"
)
TRUNCATED_POLAR_ARCHIVE_README = (
    f"{TRUNCATED_POLAR_ARCHIVE_DIRECTORY}/README.md"
)
RETAINED_RESULT_JSON = (
    f"{TRUNCATED_POLAR_ARCHIVE_DIRECTORY}/case_a/case_a_reduced.json"
)
RETAINED_RESULT_JSON_SHA256 = (
    "3b97d6cede4c98a2d79e373f676029d8a4286f9a96e986942c9991f48f2c58bf"
)
RETAINED_RESULT_STDOUT = (
    f"{TRUNCATED_POLAR_ARCHIVE_DIRECTORY}/case_a/case_a_reduced_stdout.txt"
)
RETAINED_RESULT_STDOUT_SHA256 = (
    "d2ff64aa29191114a7d41a0780b933f4d762981c85764c105ceeb5f6f07b26d4"
)
RETAINED_RESULT_NPZ_SHA256 = (
    "d475f8cc4ddde8655cf93ba1fc1e5402757d32cc00283f35f2611f21b83fc642"
)
HISTORICAL_RETAINED_RESULT_SOURCE_SHA256 = {
    "examples/cardiac_benchmark/activation.py": (
        "c1b39cd93992272e29b55907379a58cae7a08b59ddb481918b236cf2a32146f4"
    ),
    "examples/cardiac_benchmark/diagnose.py": (
        "b001e9fa333900362ee00b017ac8d145112a22c261a5ff7f5c6c348f573d4f01"
    ),
    "examples/cardiac_benchmark/geometry.py": (
        "c9052074ecdc52ac96779ff630d41febfdf458712a7895e49728b19d1271c0cb"
    ),
    "examples/cardiac_benchmark/material.py": (
        "f4a38289eda6d9178ab5ed8bb1f429b1f17100868a25bf820ffb9ea6f9ab2e08"
    ),
    "examples/cardiac_benchmark/newmark.py": (
        "4a195a6400b41c8a9f90a0070a14100a99ca007b2d4ea6a5b28d1ab96fc0287a"
    ),
    "examples/cardiac_benchmark/pressure.py": (
        "b3b3425fcd5331adb69ebd6d528a0e8b3b84271fec12fcad1883b326ee9990c9"
    ),
    "examples/cardiac_benchmark/result_io.py": (
        "4448a009f3801fa442877a1c858649b33f78c7211ff38ae19fb1430dccc78b7c"
    ),
    "examples/cardiac_benchmark/robin.py": (
        "0d33a8ce57c321ded2bb10a5fdc61d02a30aa3582839393e7110dbdab1b8d998"
    ),
    "examples/cardiac_benchmark/run.py": (
        "2e7c8564b003f484c93c227c9e2f63dc58f30dc4b45cfd5eb532e69ec138d5a2"
    ),
    "examples/cardiac_benchmark/solver.py": (
        "2445d7de9b3d5005d5b114d442ab9b6489fd3709d8cce37e5b8f3948fb48679d"
    ),
}
RETAINED_RESULT_COMMAND = (
    "python",
    "examples/cardiac_benchmark/run.py",
    "--case",
    "A",
    "--integrator",
    "be",
    "--nt",
    "1",
    "--nmu",
    "2",
    "--ntheta",
    "4",
    "--apex-offset",
    "0.2",
    "--dt",
    "0.002",
    "--tend",
    "1.0",
    "--build-dir",
    "build",
    "--out",
    "case_a_reduced.npz",
)
RETAINED_RESULT_CLAIM_BOUNDARY = (
    "The retained output checks a reduced executable pipeline.",
    "It is not full Case A or paper-curve validation.",
    "It uses a coarse open-apex mesh and global Delaunay-tetra point sampling.",
    "No external Zenodo reference curves are redistributed or compared.",
    "The broad centroid det(F) range requires mesh and resolution studies before quantitative use.",
)

# These archived reports are immutable, source-identified records of the
# historical truncated ``polar_ring`` geometry. They remain release/regression
# evidence for that configuration, but are not current Benchmark 1 validation
# evidence. Add one entry for every archived reviewed report/log pair; the
# public inventories and semantic validation come from this single list.
CURRENT_RESULT_APP_REF = "62ad760d2a1731bb9668897863ac026d3768194e"
DOMAIN_RECOVERY_APP_REF = "e07993bcf1166bd20eb87370c0b458552753e7ee"
CORRECTED_SWITCH_APP_REF = "6839c13b5bc80ec06c897684c51f503e80bd4b19"
COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
)
CURRENT_REPORT_SCHEMA = "coupfe-cardiac-reference-comparison-v2"
CURRENT_RESULT_SCHEMA = "coupfe-cardiac-result-v1"
TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS = (
    {
        "report": "case_a_fbar_1x2x4_dt0p002.report.json",
        "report_sha256": (
            "4278cfc2f282fc7a15a4b67d58e0b65ab425ab81e746067dbb786c7ea080f7c1"
        ),
        "predecessor_report_sha256": (
            "918353314c391063f79a6c101b000fb8b3ff21ca5b2d53e5b1e4f9161478c323"
        ),
        "log": "case_a_fbar_1x2x4_dt0p002.raw.stdout.txt",
        "log_sha256": (
            "db1b16b520781a18730f42f60378988d9adbe3d5be95adfb6dc76e4a129510a6"
        ),
        "result": "case_a_fbar_1x2x4_dt0p002.npz",
        "result_sha256": (
            "c8085d2e3ab92c572e617fab75b38ecb27cb83f7966f3b89c97261570dfceed4"
        ),
        "result_size_bytes": 162700,
        "case": "A",
        "reference_case": "step_0A",
        "app_ref": CURRENT_RESULT_APP_REF,
        "formulation": "hex8_fbar",
        "solver": "core-newton",
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 240,
        "mesh": {
            "n_t": 1,
            "n_mu": 2,
            "n_theta": 4,
            "nodes": 24,
            "elements": 8,
            "degrees_of_freedom": 72,
        },
    },
    {
        "report": (
            "case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json"
        ),
        "report_sha256": (
            "103b67171bd15326983295940e267e0d0c7884481f05fbed2bbd1d3042ce2ddf"
        ),
        "log": (
            "case_a_fbar_1x2x4_dt0p002_corrected_switch.raw.stdout.txt"
        ),
        "log_sha256": (
            "76aefb1b9500420ece439f28a635c41cb868a39208969f67234b3a9b1935e1a7"
        ),
        "result": "case_a_fbar_1x2x4_dt0p002_corrected_switch.npz",
        "result_sha256": (
            "39013056c9ae9b90e24525894af2a9b52790b5f6b1604db5afe0e91624453dfc"
        ),
        "result_size_bytes": 163645,
        "case": "A",
        "reference_case": "step_0A",
        "app_ref": CORRECTED_SWITCH_APP_REF,
        "core_ref": APPROVED_PUBLIC_CORE_REF,
        "formulation": "hex8_fbar",
        "solver": "core-newton",
        "method_metadata": True,
        "current_claim_boundary": True,
        "material_model_id": COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID,
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 240,
        "mesh": {
            "n_t": 1,
            "n_mu": 2,
            "n_theta": 4,
            "nodes": 24,
            "elements": 8,
            "degrees_of_freedom": 72,
            "topology": "polar_ring",
        },
    },
    {
        "report": "case_b_local_pressure_2x12x16_dt0p002.report.json",
        "report_sha256": (
            "7e06f6f88c9db3970e93a689da762e6ab8d9c38cd75d668747cc4c286a2020ae"
        ),
        "predecessor_report_sha256": (
            "945252b6fce3b3fb796cc554037e24eecde942a2454e6cc27bfd3f0bb69ce6c0"
        ),
        "log": "case_b_local_pressure_2x12x16_dt0p002.raw.stdout.txt",
        "log_sha256": (
            "a8d2c09adf2e291c2b3b1966290e4aa3f72b3f151729e5a884293c95f80bfec6"
        ),
        "result": "case_b_local_pressure_2x12x16_dt0p002.npz",
        "result_sha256": (
            "04a99ccdc1c84d8d6d7d425fd869b44ee4983431fa1c112a18e58414651b1549"
        ),
        "result_size_bytes": 1077108,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": CURRENT_RESULT_APP_REF,
        "formulation": "hex8_local_pressure_p0_condensed_logj",
        "solver": "petsc-snes",
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 241,
        "ring_rotation_relative_degrees": (
            17.030798282727744,
            7.641732032376751,
            -10.581144840595954,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 12,
            "n_theta": 16,
            "nodes": 624,
            "elements": 384,
            "degrees_of_freedom": 1872,
        },
    },
    {
        "report": "case_b_local_pressure_2x12x16_dt0p004.report.json",
        "report_sha256": (
            "abd3c28c18e8f44cc9c5529f7661d69158d516c7070a893656b3cb437556de79"
        ),
        "predecessor_report_sha256": (
            "4c09013f127fe213c34e2d1354d7c723b0275fb2f0002108a024f48c848cb7ad"
        ),
        "log": "case_b_local_pressure_2x12x16_dt0p004.raw.stdout.txt",
        "log_sha256": (
            "00073af16e0e93ccf9a78ed028847874f3950a0c188d82497d52464955a6538a"
        ),
        "result": "case_b_local_pressure_2x12x16_dt0p004.npz",
        "result_sha256": (
            "c5278d167f9fe2d860161525761f78f30ea7bc9d7d94cce57478368a41bec84d"
        ),
        "result_size_bytes": 601664,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": CURRENT_RESULT_APP_REF,
        "formulation": "hex8_local_pressure_p0_condensed_logj",
        "solver": "petsc-snes",
        "dt_s": 0.004,
        "steps": 250,
        "peak_index": 120,
        "ring_rotation_relative_degrees": (
            17.03286621981411,
            7.64216954308651,
            -10.579987979929172,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 12,
            "n_theta": 16,
            "nodes": 624,
            "elements": 384,
            "degrees_of_freedom": 1872,
        },
    },
    {
        "report": "case_b_local_pressure_2x24x32_dt0p002.report.json",
        "report_sha256": (
            "f707b258da14fe7f3a2a990e22f7b9c50286e444bc392b8a4d290b7b33f1209c"
        ),
        "predecessor_report_sha256": (
            "2a2b900eb5f396b300f2e3fb383cabd1f1253b900c536263e4ec48ffebea7190"
        ),
        "log": "case_b_local_pressure_2x24x32_dt0p002.raw.stdout.txt",
        "log_sha256": (
            "efb7196363818b1a1744cc38e5463f309fcfeb25ce358cb8f0a53cb62909984e"
        ),
        "result": "case_b_local_pressure_2x24x32_dt0p002.npz",
        "result_sha256": (
            "26088c11a9b9a9ca7af625c85bcc10569bb324efb4a14f02040082e912ee8356"
        ),
        "result_size_bytes": 1376376,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": CURRENT_RESULT_APP_REF,
        "formulation": "hex8_local_pressure_p0_condensed_logj",
        "solver": "petsc-snes",
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 241,
        "ring_rotation_relative_degrees": (
            45.585039667432596,
            56.825136160398955,
            63.07751853519841,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 24,
            "n_theta": 32,
            "nodes": 2400,
            "elements": 1536,
            "degrees_of_freedom": 7200,
        },
    },
    {
        "report": "case_b_local_pressure_2x24x32_dt0p004.report.json",
        "report_sha256": (
            "334ba3032b9232295e8b56fa4f21971da97cb52aefa9b41cce72c9a0067ef38d"
        ),
        "predecessor_report_sha256": (
            "cea4c7ab7b0ada3c8c5a67b4ce001114cc04490531352e83589271aa06a64047"
        ),
        "log": "case_b_local_pressure_2x24x32_dt0p004.raw.stdout.txt",
        "log_sha256": (
            "954b037c5ff3512099adc3f0dccaa7a7f250150d1c9dbe839cf4cc5846b5715b"
        ),
        "result": "case_b_local_pressure_2x24x32_dt0p004.npz",
        "result_sha256": (
            "55edf5691632b2470faeb6bcc462d6398a610e3537d772281c655b0f000cd579"
        ),
        "result_size_bytes": 965144,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": DOMAIN_RECOVERY_APP_REF,
        "formulation": "hex8_local_pressure_p0_condensed_logj",
        "solver": "petsc-snes",
        "function_domain_diagnostics": True,
        "dt_s": 0.004,
        "steps": 250,
        "peak_index": 120,
        "ring_rotation_relative_degrees": (
            45.358965079195656,
            56.582440012318926,
            62.94400687196972,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 24,
            "n_theta": 32,
            "nodes": 2400,
            "elements": 1536,
            "degrees_of_freedom": 7200,
        },
    },
    {
        "report": "case_b_local_pressure_2x36x48_dt0p002.report.json",
        "report_sha256": (
            "d409ecaac0c0abf418fce3ab2f0549979d38e02b9d73381d56982f7fc4e3bf14"
        ),
        "predecessor_report_sha256": (
            "83d513e933b1144c67b6066c03cb27fba0d212d8477c521a548acdd6e7917956"
        ),
        "log": "case_b_local_pressure_2x36x48_dt0p002.raw.stdout.txt",
        "log_sha256": (
            "a08cabdca2ae6246b232b1c78202af0d6f78aa1834451be673b380b4e22a4077"
        ),
        "result": "case_b_local_pressure_2x36x48_dt0p002.npz",
        "result_sha256": (
            "6bcb5e0b98a044335d452be021cb765223a7ed3c21e4a2fc89496ba8d66e911a"
        ),
        "result_size_bytes": 1997944,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": DOMAIN_RECOVERY_APP_REF,
        "formulation": "hex8_local_pressure_p0_condensed_logj",
        "solver": "petsc-snes",
        "function_domain_diagnostics": True,
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 241,
        "ring_rotation_relative_degrees": (
            -46.410399086878016,
            -46.32253057989533,
            -48.67835724664141,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 36,
            "n_theta": 48,
            "nodes": 5328,
            "elements": 3456,
            "degrees_of_freedom": 15984,
        },
    },
    {
        "report": "case_b_fbar_2x24x32_dt0p002.report.json",
        "report_sha256": (
            "77ec1d5735d6353dad386ce84ea1052a4d5b47e2046cfd140e405f3061d7c2f9"
        ),
        "predecessor_report_sha256": (
            "5595ce3b69a65e41aab5bdc52f85a3d7de7f213bf23b5dfbbed673691288b49f"
        ),
        "log": "case_b_fbar_2x24x32_dt0p002.raw.stdout.txt",
        "log_sha256": (
            "184451eb047f24e0bf6c982a0a5a5857e3bbe4e9a87c8f28b775066dac3a9384"
        ),
        "result": "case_b_fbar_2x24x32_dt0p002.npz",
        "result_sha256": (
            "87636b9626e88b7e6b1d0d408659028f412c25b1accb62ca842cfe41541e84d2"
        ),
        "result_size_bytes": 1369028,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": CURRENT_RESULT_APP_REF,
        "formulation": "hex8_fbar",
        "solver": "petsc-snes",
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 241,
        "ring_rotation_relative_degrees": (
            48.43437798171243,
            59.460521088807724,
            63.27314210826961,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 24,
            "n_theta": 32,
            "nodes": 2400,
            "elements": 1536,
            "degrees_of_freedom": 7200,
        },
    },
    {
        "report": "case_b_fbar_2x36x48_dt0p002.report.json",
        "report_sha256": (
            "3415100632097b9f9fbf452b346a18fc1a585263d40ab501d43638af1b5450e4"
        ),
        "predecessor_report_sha256": (
            "f889ec8374869d9f23cedd7be83520ffc730682e7629598a53165d6a60039f2c"
        ),
        "log": "case_b_fbar_2x36x48_dt0p002.raw.stdout.txt",
        "log_sha256": (
            "9dd36c7e9e326bf1f3e76af381c37e53ea0e9cf3cb916abdd37d9e4f4e04a98c"
        ),
        "result": "case_b_fbar_2x36x48_dt0p002.npz",
        "result_sha256": (
            "d1c25bb6a4cd62cf662bed00f621643098ca00e8654ddf469b491d56af13c38c"
        ),
        "result_size_bytes": 1835588,
        "case": "B",
        "reference_case": "step_0B",
        "app_ref": CURRENT_RESULT_APP_REF,
        "formulation": "hex8_fbar",
        "solver": "petsc-snes",
        "dt_s": 0.002,
        "steps": 500,
        "peak_index": 241,
        "ring_rotation_relative_degrees": (
            -43.72087401574177,
            -46.40910968910197,
            -44.69799642640277,
        ),
        "mesh": {
            "n_t": 2,
            "n_mu": 36,
            "n_theta": 48,
            "nodes": 5328,
            "elements": 3456,
            "degrees_of_freedom": 15984,
        },
    },
)

CURRENT_RESULT_SOURCE_SHA256 = {
    "examples/cardiac_benchmark/activation.py": (
        "c1b39cd93992272e29b55907379a58cae7a08b59ddb481918b236cf2a32146f4"
    ),
    "examples/cardiac_benchmark/geometry.py": (
        "9bfba31ebfed7895416533097d784bf01bc67bc14813de2d42adec5aa8cfbaad"
    ),
    "examples/cardiac_benchmark/local_pressure.py": (
        "1f6738689f4f05890439aef1c804b54099b73d3e2575eb171f923b7a13850285"
    ),
    "examples/cardiac_benchmark/material.py": (
        "3f45b1d59a4008c90b72405becaf910a6d09a10747af5a0817e1695463e0520c"
    ),
    "examples/cardiac_benchmark/newmark.py": (
        "4a195a6400b41c8a9f90a0070a14100a99ca007b2d4ea6a5b28d1ab96fc0287a"
    ),
    "examples/cardiac_benchmark/pressure.py": (
        "b3b3425fcd5331adb69ebd6d528a0e8b3b84271fec12fcad1883b326ee9990c9"
    ),
    "examples/cardiac_benchmark/result_io.py": (
        "05b19f3f501430e3d20f963cb86a2a040c18f8eda4b8ef92ea9ee3c205f18267"
    ),
    "examples/cardiac_benchmark/robin.py": (
        "0d33a8ce57c321ded2bb10a5fdc61d02a30aa3582839393e7110dbdab1b8d998"
    ),
    "examples/cardiac_benchmark/run.py": (
        "627e7e0f68c1f7e82caf804b174c6f4ac7a9125e8037294f20ac8703e39a1d8b"
    ),
    "examples/cardiac_benchmark/sampling.py": (
        "b6a9828f86c3a635df7c7967589cec8734ccfece629a1721590852e6ea50f32b"
    ),
    "examples/cardiac_benchmark/solver.py": (
        "0f1823ec2fe8148341d575750262513f0090c30011d688424934b95d58c30c8d"
    ),
}

DOMAIN_RECOVERY_SOURCE_SHA256 = {
    **CURRENT_RESULT_SOURCE_SHA256,
    "examples/cardiac_benchmark/local_pressure.py": (
        "2930e211df9f3bbe481504a1f6a866193a4e2f145ae1b0004c6f72f23f46c1b2"
    ),
    "examples/cardiac_benchmark/solver.py": (
        "7cec77dba72d646d115a13d338bca973ee95ebe27301521a40838a9726818598"
    ),
}

CORRECTED_SWITCH_SOURCE_SHA256 = {
    "examples/cardiac_benchmark/activation.py": (
        "c1b39cd93992272e29b55907379a58cae7a08b59ddb481918b236cf2a32146f4"
    ),
    "examples/cardiac_benchmark/geometry.py": (
        "9bfba31ebfed7895416533097d784bf01bc67bc14813de2d42adec5aa8cfbaad"
    ),
    "examples/cardiac_benchmark/local_pressure.py": (
        "2930e211df9f3bbe481504a1f6a866193a4e2f145ae1b0004c6f72f23f46c1b2"
    ),
    "examples/cardiac_benchmark/material.py": (
        "253242fd97f225e7cbfd1ba202bd146973a481314160f493fefcbc3e92e54eea"
    ),
    "examples/cardiac_benchmark/newmark.py": (
        "4a195a6400b41c8a9f90a0070a14100a99ca007b2d4ea6a5b28d1ab96fc0287a"
    ),
    "examples/cardiac_benchmark/pressure.py": (
        "b3b3425fcd5331adb69ebd6d528a0e8b3b84271fec12fcad1883b326ee9990c9"
    ),
    "examples/cardiac_benchmark/result_io.py": (
        "05b19f3f501430e3d20f963cb86a2a040c18f8eda4b8ef92ea9ee3c205f18267"
    ),
    "examples/cardiac_benchmark/robin.py": (
        "0d33a8ce57c321ded2bb10a5fdc61d02a30aa3582839393e7110dbdab1b8d998"
    ),
    "examples/cardiac_benchmark/run.py": (
        "c92f8348f76dd67f341ff13ad253aed07692627c9b378abbdf17c162658a54fd"
    ),
    "examples/cardiac_benchmark/sampling.py": (
        "b6a9828f86c3a635df7c7967589cec8734ccfece629a1721590852e6ea50f32b"
    ),
    "examples/cardiac_benchmark/solver.py": (
        "7cec77dba72d646d115a13d338bca973ee95ebe27301521a40838a9726818598"
    ),
}

# Exact executable source shipped by this release. This is deliberately
# separate from the immutable source-checkpoint maps above: changing the
# current driver must not relabel a historical result as if it were rerun.
CURRENT_RELEASE_SOURCE_SHA256 = {
    **DOMAIN_RECOVERY_SOURCE_SHA256,
    "examples/cardiac_benchmark/local_pressure.py": (
        "6447e63b4810363d8cc26fa6deaabe6d729e2542dbab003328a79f365c4b3d26"
    ),
    "examples/cardiac_benchmark/activation.py": (
        "61251c4991c6282a903eb768102c14fb09671059ee93aada0ddc3e8bdf86c2ae"
    ),
    "examples/cardiac_benchmark/benchmark_parameters.py": (
        "842ea66f915dc2153df653ebf9dfcf7a8bf0dcd34beb32ea7ddcc3498e8f41d3"
    ),
    "examples/cardiac_benchmark/boundary_audit.py": (
        "b6a2db0d85e6efb62c5fa9ad098b725928a2ee12475f1df8979b837f214cfb77"
    ),
    "examples/cardiac_benchmark/robin.py": (
        "712b9a25937136a5ec57cfd65a680814f58c092b9911eaf3258a8f4c59f2aed2"
    ),
    "examples/cardiac_benchmark/consistent_mass.py": (
        "add9e0e54d3fc7ebb1ce84ffa1ecc774d0ba8a366bb2751d141d3d5fefd64629"
    ),
    "examples/cardiac_benchmark/distributed_local_pressure.py": (
        "7021901afba7fc57a4d0dbc28da2bdab7984686fab2e17d93ed313e507291cff"
    ),
    "examples/cardiac_benchmark/distributed_mass.py": (
        "ca77a000a2be6518c22e94dfa436ebe4cef0d026d9bec112b530df037f17aeb7"
    ),
    "examples/cardiac_benchmark/distributed_material.py": (
        "475b7209d34ebd5cc7054102dfcafc756ca18b869acdf689a93fa5aba8c76e9b"
    ),
    "examples/cardiac_benchmark/distributed_solver.py": (
        "fb0582788591db37e76dd0916b4c617a3293f57efc7840413649498261bf848c"
    ),
    "examples/cardiac_benchmark/generalized_alpha.py": (
        "5197b3f22df6852043e0f987ec7a433eb3a97fd5abe963baa34f41aa05a22518"
    ),
    "examples/cardiac_benchmark/geometry.py": (
        "4ed0ec4c4e6421a1a9495eb591c2e316b1d32d806de5ba1101ddb6daf090030f"
    ),
    "examples/cardiac_benchmark/material.py": (
        "ab71054e9a035527e70f6c359422cc9dc276fa214cbd086111b9771aa055d0b3"
    ),
    "examples/cardiac_benchmark/run.py": (
        "a44632e1a90257fa0ba6331d393bffbde36cac89aa04551866829979455c14af"
    ),
    "examples/cardiac_benchmark/run_mpi.py": (
        "150b50a4bc452580b82f40699d428bcbb4d215cecc456c9a9ee8aa897b16a57d"
    ),
    "examples/cardiac_benchmark/solver.py": (
        "d9dd5c9fe5e50870c4e137f8991be27177a61eff24c8b281a1d2686af1070bfc"
    ),
    "examples/cardiac_benchmark/structural_directions.py": (
        "b09f678d67d271d087961c7059f2b8570326f8c5227b43a307a9a3c075cb926d"
    ),
    "examples/cardiac_benchmark/tbar_laplace.py": (
        "1e9850ea43acac269059b33787c7d888c2b24b4365bfb0dae1b92d52364995e2"
    ),
    "examples/cardiac_benchmark/viscous_evidence.py": (
        "63350a12a81ad256bfc898d8b8c3c811f315cbbd7192b16ac57c5390db49f095"
    ),
}

# A report names its source checkpoint through ``spec["app_ref"]``.  Keeping
# hashes by checkpoint lets an immutable historical report coexist with later
# reruns while still requiring every packaged executable source set to match a
# reviewed result-producing checkpoint.  Add a new entry after the rerun
# source commit exists; do not replace the historical entry.
CURRENT_RESULT_SOURCE_CHECKPOINTS = {
    CURRENT_RESULT_APP_REF: CURRENT_RESULT_SOURCE_SHA256,
    DOMAIN_RECOVERY_APP_REF: DOMAIN_RECOVERY_SOURCE_SHA256,
    CORRECTED_SWITCH_APP_REF: CORRECTED_SWITCH_SOURCE_SHA256,
}
# Reporting utilities do not produce simulation state. Result-producing
# checkpoint hashes intentionally exclude them, while this separate map binds
# the current selection and comparison code.
CURRENT_REPORTING_SOURCE_SHA256 = {
    "examples/cardiac_benchmark/post.py": (
        "75448a510d8e5abd4a9eb20f8b3bcfba3b4030e9a6bee3039dcea63b98d9272c"
    ),
    "examples/cardiac_benchmark/compare_fenics_case_b.py": (
        "2a9ab7e9415756be4908d293f985c6cc89838a4608c090143db94bc7b33f6768"
    ),
    "examples/cardiac_benchmark/compare_mpi_rank_gate.py": (
        "dfa13431014bc7487063cfddf5692a3a7d409ebb1622fbd255a587dc7374409b"
    ),
    "examples/cardiac_benchmark/compare_step2b_case_b.py": (
        "16a5982150674fa572f56009a476907b529182f5bff4b5ca7c023d3ec0cf39bf"
    ),
    "examples/cardiac_benchmark/plot_step2b_case_b.py": (
        "6d605acab25a3f300accb3fde07baf0479a0d74da9e896e7e75c1625c931be7d"
    ),
    "examples/cardiac_benchmark/step2b_case_b_reference_hashes.json": (
        "8392e28a1e6971b0c5c26e79aa2d0e86a48935cb7174c2454920b08a486631ac"
    ),
    "examples/cardiac_benchmark/step2b_case_b_runtime_source_hashes.json": (
        "d39ccbd6e67d7517b31a536f0c34472afb770bfd76916a3572e20d52acf39a41"
    ),
    "examples/cardiac_benchmark/compare_tip_refine_full_cycle.py": (
        "d091c302650fedf5714526788a0bed0f61aa2e86690df11deb76fb83f2649268"
    ),
    "examples/cardiac_benchmark/plot_step2b_current_rerun.py": (
        "f6834780bce112b9fc0d1c42bd5161b9f581001c8faa5c5ee0642350a80b92e7"
    ),
}
PETSC_FUNCTION_DOMAIN_REJECTION_API = "nonfinite residual for PETSc BT"

REFERENCE_ARCHIVE_IDENTITY = {
    "filename": "benchmark_article_data.zip",
    "size_bytes": 23180741494,
    "md5": "75602be4777c4ca2262c2bcfd2134b15",
    "sha256": "134951af5e38d147b0223f0a83666eb3fe1b75acb5bfa9f1b9aa30f255f8f1f5",
}
REFERENCE_FIGURES_PY_IDENTITY = {
    "filename": "results_time_curves/figures.py",
    "size_bytes": 24076,
    "sha256": "f8f519b357349341207faea4b57bfafc1a311aadccefe4968ecbdb37339c8a5b",
}
REFERENCE_SELECTION_POLICY = (
    "Select the exact 10 Case A/Case B team files explicitly loaded by upstream "
    "results_time_curves/figures.py. Reject missing selected files and unexpected "
    "matching files. Accept the unselected base-name SimVascular file only when "
    "byte-identical to selected SimVascular P2, then exclude it as a duplicate alias."
)
REFERENCE_TEAM_ORDER = (
    "carpentry",
    "ambit",
    "4C",
    "simula",
    "chimera",
    "cheart",
    "lifex",
    "simvascular_p1p1",
    "simvascular_p2",
    "comsol",
)
REFERENCE_MANIFEST_VARIABLES = {
    "step_0A": "TEAMS_DATASETS_0A",
    "step_0B": "TEAMS_DATASET_0B",
}
REFERENCE_TEAM_SHA256 = {
    "step_0A": {
        "4C": "69e2793b50aaf3de6a749cee0776817a2c1561a670776265873f603463e16a6f",
        "ambit": "5ff4c59ddcdf4da4257c919237081cead7460212c03f569e36cda027c66d3b64",
        "carpentry": "945f3b5efb8a5fff284c80de751909faa2fdb7da6e324bb2e4f079c92f538dbb",
        "cheart": "44da014817ca9fa847954cade3509e6a8bd24c13bc8a43280a430b8e7d0bd0e1",
        "chimera": "91478419a1a1226b1db86233aca009f2db3d9fcd8d78b718689f1184cd4fc603",
        "comsol": "86f314345095e14936791ddd58f980807713c2a7e11528ba319d9ac0d05a310f",
        "lifex": "bbf79b72dc73580cb4ce38f8f253c8f71c224ff43ecd05381aec01f9c23ab38f",
        "simula": "e5804d9f7bb9f99690512a55fbab30c45cffcf3fa812e014113fb4779f386575",
        "simvascular_p1p1": "eb2f63bd90484c27d1fbbfdcdf7b148b8a1a3df04151c7e63af4c5234abd578a",
        "simvascular_p2": "7b183de36449f8dd3874abfa1f1eccdcf3f6e4ed16d2e8535b7c4f3daa4f81bd",
    },
    "step_0B": {
        "4C": "d066fdc316d92c5e76c2da08aeb8612b9c3f7d5fc2a6f364a970adaeed505a99",
        "ambit": "e9d287ffad46fd4a96f7b0f9f6187067d1ff7e392dbf188368ed9ce857e0895c",
        "carpentry": "5a223507bdd54daf0d775a8c14b07ab26933923063b171cd16daf61625b2a7cf",
        "cheart": "08460c6708673b57f07ad1475db88381da9608b4ef5f34248c2558a3bec38da6",
        "chimera": "9062c2204d21b6bc3146711f061cd82af5be19543f4fcce2ff2cfde351b1cbcd",
        "comsol": "80121d17056e46cf7a2419fe25f70c14afac7f1599e2de08db33124582191c88",
        "lifex": "f7e72fae62055f2b202196e18f5e8527fc4476383d17b376f6c105f13e7aa4f4",
        "simula": "4dd3aace9580dccd803e0e70bb74524bc061015b8dfc27274a00f747369f16c4",
        "simvascular_p1p1": "d82772ff96ce900cf45cc8cca50985f5bb9042476742637adc77765b363db9e3",
        "simvascular_p2": "1ca6f64f7eb63f1d7aef2023211c9f14d2d7fdfc4c9d232fc914d7a3c3ed6cf5",
    },
}
REFERENCE_TEAM_SIZE_BYTES = {
    "4C": 8651,
    "ambit": 23037,
    "carpentry": 8651,
    "cheart": 8651,
    "chimera": 15172,
    "comsol": 7835,
    "lifex": 7801,
    "simula": 8651,
    "simvascular_p1p1": 8651,
    "simvascular_p2": 8651,
}
REFERENCE_EXCLUDED_ALIAS_SHA256 = {
    "step_0A": "7b183de36449f8dd3874abfa1f1eccdcf3f6e4ed16d2e8535b7c4f3daa4f81bd",
    "step_0B": "1ca6f64f7eb63f1d7aef2023211c9f14d2d7fdfc4c9d232fc914d7a3c3ed6cf5",
}
CORRECTION_REASON = (
    "The predecessor selected 11 files with a wildcard and double-counted the "
    "byte-identical SimVascular base-name alias and SimVascular P2 curve. This "
    "successor selects the 10 files explicitly used by upstream "
    "results_time_curves/figures.py."
)
CORRECTION_PREDECESSOR_REVISION = "7f8e726cd2a79ae2ad13ebac4d9c39bca5cec8b2"
CURRENT_REPORT_BOUNDED_CLAIM = (
    "This report compares one completed, source-identified CoupFE-Cardiac run "
    "with the separately distributed benchmark time curves on a common grid. "
    "It is example-level evidence, not clinical, device, or broad solver "
    "validation. It does not by itself establish mesh/time convergence, a "
    "closed nondegenerate ventricular geometry, or a unique twist direction."
)
CURRENT_REPORT_DISCRETIZATION_BOUNDED_CLAIM = (
    "This report compares one completed, source-identified CoupFE-Cardiac run "
    "with the separately distributed benchmark time curves on a common grid. "
    "It is example-level evidence, not clinical, device, or broad solver "
    "validation. It does not by itself establish mesh/time convergence, "
    "equivalence to the reference P2 discretization, or a unique twist "
    "direction. A closed-topology statement, when present, is bounded to the "
    "retained pre-solve geometry and boundary audits."
)
CURRENT_RUNTIME_VERSIONS = {
    "coupfe_version": "0.0.1",
    "numpy_version": "1.26.4",
    "python_version": "3.10.8",
    "scipy_version": "1.15.2",
}
CORE_NEWTON_CONFIGURATION = {
    "independent_acceptance_atol": 1.0e-14,
    "max_it": 40,
    "name": "core-newton",
    "rtol": 1.0e-8,
}
PETSC_SNES_CONFIGURATION = {
    "atol": 1.0e-10,
    "dirichlet_support": "none",
    "factor_solver_type": "petsc",
    "ksp_type": "preonly",
    "line_search_configuration_api": "namespaced PETSc option",
    "line_search_type": "bt",
    "matrix_scope": "one solver instance per application run",
    "max_it": 60,
    "name": "petsc-snes",
    "pc_type": "lu",
    "petsc4py_version": "3.18.4",
    "petsc_version": "3.18.4",
    "rtol": 1.0e-9,
    "settings_source": "recovered 2026-06-27 Case B development adapter",
    "snes_type": "newtonls",
    "stol": 1.0e-12,
}

LEGAL_MARKERS = {
    "LICENSE": (
        "Apache License",
        "Version 2.0, January 2004",
    ),
    "NOTICE": (
        "Apache-2.0 AND CC-BY-4.0 AND MIT",
        "examples/cardiac_benchmark/activation.py",
        "examples/cardiac_benchmark/fiber_crosscheck.py",
        "Henrik Finsberg, Joakim Sundnes, and Jonas van den Brink",
        "identified as modifications of the cited dataset",
    ),
    "THIRD_PARTY_NOTICES.md": (
        "Incorporated source: Finsberg, Sundnes, and van den Brink cardiac benchmark",
        "`examples/cardiac_benchmark/activation.py` — CC BY 4.0",
        "Henrik Finsberg, Joakim Sundnes",
        "Jonas van den Brink",
        "be92da5dbc1fd26d424bf88ef7db13b4",
        "df6e5f03e644cb055ba3649f901030ba1d18840ff2ea94f6c71fd12bded28185",
        "uses that conservative license basis",
        "Modification indication",
        "325d17d850c2e2032abb85a4191a5795d3008ab7",
        "Incorporated source: cardiac benchmark toolkit",
        "e8d47553cfc83eb274eba3e177de33148e7f441c",
        "do redistribute transformed CC-BY-4.0 material",
        "10.5281/zenodo.14260459",
    ),
    "docs/LICENSE.md": (
        "Creative Commons Attribution 4.0 International",
        "LICENSES/CC-BY-4.0.txt",
        "Apache-2.0 AND CC-BY-4.0 AND MIT",
    ),
    "LICENSES/CC-BY-4.0.txt": (
        "Attribution 4.0 International",
        "Creative Commons Attribution 4.0 International Public License",
        "Section 3 -- License Conditions.",
    ),
    "LICENSES/Reidmen-cardiac_benchmark_toolkit-MIT.txt": (
        "MIT License",
        "Copyright (c) 2023 Reidmen Arostica",
        "The above copyright notice and this permission notice shall be included",
    ),
}

CORE_VERIFIER_MARKERS = {
    ".github/scripts/check_runtime_core.py": (
        f'PUBLIC_CORE_URL = "{PUBLIC_CORE_URL}"',
        f'PUBLIC_CORE_REF = "{APPROVED_PUBLIC_CORE_REF}"',
    ),
    "examples/cardiac_benchmark/post.py": (
        f'PUBLIC_CORE_URL = "{PUBLIC_CORE_URL}"',
        f'PUBLIC_CORE_REF = "{APPROVED_PUBLIC_CORE_REF}"',
    ),
}

RETAINED_FIGURE_RENDERER = (
    "examples/cardiac_benchmark/plot_retained_comparisons.py"
)
RETAINED_FIGURE_RENDERER_SHA256 = (
    "a8ec416f2333eecc14e4c0bcf7cdffde9b2c6fb026fa6cc0743c0d998518da7d"
)
STEP0A_RETAINED_REPORT_SCHEMA = (
    "coupfe-cardiac-step0a-retained-comparison-v1"
)
STEP0A_RETAINED_REPORT = (
    "examples/cardiac_benchmark/results/"
    "case_a_local_pressure_4x36x32_dt0p001.report.json"
)
STEP0A_RETAINED_REPORT_SHA256 = (
    "bbd26f3b30819ff2b67ffb48c9ad52cc9825c7fa0486e3984673c9e349bf82b1"
)
STEP0A_RETAINED_REPORT_SIZE_BYTES = 62674
STEP0A_RETAINED_APP_REF = "016a4f9eec6f2a4c74d10c734ddff3e24cf343de"
STEP0B_PREFIX_DIAGNOSTIC_REPORT = (
    "examples/cardiac_benchmark/results/"
    "step0b_case_b_clean_frame_0p32.report.json"
)
STEP0B_PREFIX_DIAGNOSTIC_REPORT_SHA256 = (
    "faef788ede2b42a175ed422292f64d916a0a30088faabff1e9252408388c0a3f"
)
STEP0B_PREFIX_DIAGNOSTIC_REPORT_SIZE_BYTES = 12385
STEP0B_PREFIX_DIAGNOSTIC_SCHEMA = (
    "coupfe-cardiac-step0b-prefix-diagnostic-v1"
)
STEP0B_PREFIX_DIAGNOSTIC_APP_REF = (
    "056c02df7c2a56bbc36e41973b7c8a8d8c917e2a"
)
STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_SHA256 = (
    "bec13f9ab1bbd9e50116e05b7342501e6ff8992c5710f62b83bcd02a089b6cf1"
)
STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST = {
    "examples/cardiac_benchmark/activation.py": (
        "0d3075b4b7866577e3a7676679efa152954c443ef6bcd0c57116b99baeef30a5"
    ),
    "examples/cardiac_benchmark/benchmark_parameters.py": (
        "842ea66f915dc2153df653ebf9dfcf7a8bf0dcd34beb32ea7ddcc3498e8f41d3"
    ),
    "examples/cardiac_benchmark/boundary_audit.py": (
        "b6a2db0d85e6efb62c5fa9ad098b725928a2ee12475f1df8979b837f214cfb77"
    ),
    "examples/cardiac_benchmark/consistent_mass.py": (
        "add9e0e54d3fc7ebb1ce84ffa1ecc774d0ba8a366bb2751d141d3d5fefd64629"
    ),
    "examples/cardiac_benchmark/distributed_local_pressure.py": (
        "7021901afba7fc57a4d0dbc28da2bdab7984686fab2e17d93ed313e507291cff"
    ),
    "examples/cardiac_benchmark/distributed_mass.py": (
        "ca77a000a2be6518c22e94dfa436ebe4cef0d026d9bec112b530df037f17aeb7"
    ),
    "examples/cardiac_benchmark/distributed_material.py": (
        "475b7209d34ebd5cc7054102dfcafc756ca18b869acdf689a93fa5aba8c76e9b"
    ),
    "examples/cardiac_benchmark/distributed_solver.py": (
        "fb0582788591db37e76dd0916b4c617a3293f57efc7840413649498261bf848c"
    ),
    "examples/cardiac_benchmark/generalized_alpha.py": (
        "5197b3f22df6852043e0f987ec7a433eb3a97fd5abe963baa34f41aa05a22518"
    ),
    "examples/cardiac_benchmark/geometry.py": (
        "413e36086b5c7ac26d8a15e3d2ba15299a9d4c10e9ca6ea88e21f0329eededdb"
    ),
    "examples/cardiac_benchmark/local_pressure.py": (
        "6447e63b4810363d8cc26fa6deaabe6d729e2542dbab003328a79f365c4b3d26"
    ),
    "examples/cardiac_benchmark/material.py": (
        "ab71054e9a035527e70f6c359422cc9dc276fa214cbd086111b9771aa055d0b3"
    ),
    "examples/cardiac_benchmark/pressure.py": (
        "b3b3425fcd5331adb69ebd6d528a0e8b3b84271fec12fcad1883b326ee9990c9"
    ),
    "examples/cardiac_benchmark/result_io.py": (
        "05b19f3f501430e3d20f963cb86a2a040c18f8eda4b8ef92ea9ee3c205f18267"
    ),
    "examples/cardiac_benchmark/robin.py": (
        "402e2a9cf6d48b1e9de5e9d40e16dc6caaa3dbd451a55d57285534c87fba7e3d"
    ),
    "examples/cardiac_benchmark/run.py": (
        "5fe3faa02589ac99de5e0f232b61c7f61bdd7f7b6da5d4d9f3f3d98235afa537"
    ),
    "examples/cardiac_benchmark/run_mpi.py": (
        "e01385b65464ab47c8d17990f9c2d8ee868b62c481e4206317b1303c6d3c8690"
    ),
    "examples/cardiac_benchmark/sampling.py": (
        "b6a9828f86c3a635df7c7967589cec8734ccfece629a1721590852e6ea50f32b"
    ),
    "examples/cardiac_benchmark/structural_directions.py": (
        "b09f678d67d271d087961c7059f2b8570326f8c5227b43a307a9a3c075cb926d"
    ),
    "examples/cardiac_benchmark/tbar_laplace.py": (
        "84f30f2e943d6611b7177bb4bbe160b8ef8c27d728e5384d697f1ba0e21c3725"
    ),
}
STEP0B_PREFIX_DIAGNOSTIC_MESH_SPLIT_APP_REF = (
    "a2006b78104109c625ea3c502753b5cff15452d4"
)
TIP_REFINE_FULL_CYCLE_REPORT = (
    "examples/cardiac_benchmark/results/"
    "step0b_tip6p0_full_cycle_comparison.report.json"
)
TIP_REFINE_FIGURE_RENDERER = (
    "examples/cardiac_benchmark/compare_tip_refine_full_cycle.py"
)
TIP_REFINE_FIGURE_RENDERER_SHA256 = (
    "d091c302650fedf5714526788a0bed0f61aa2e86690df11deb76fb83f2649268"
)
STEP2B_RERUN_FIGURE_RENDERER = (
    "examples/cardiac_benchmark/plot_step2b_current_rerun.py"
)
STEP2B_RERUN_FIGURE_RENDERER_SHA256 = (
    "f6834780bce112b9fc0d1c42bd5161b9f581001c8faa5c5ee0642350a80b92e7"
)
STEP2B_RERUN_FULL_CYCLE_REPORT = (
    "examples/cardiac_benchmark/results/"
    "step2b_current_rerun_comparison.report.json"
)
STEP0B_PREFIX_DIAGNOSTIC_BENCHMARK_IDENTITY = {
    "active_stress_enabled": False,
    "benchmark": "Benchmark 1",
    "case": "step_0B",
    "configuration_id": "benchmark-1-step-0-case-B-pressure-only",
    "identity_scope": "explicit archive metadata",
    "pressure_enabled": True,
}
STEP0B_PREFIX_DIAGNOSTIC_CONFIGURATION = {
    "density_kg_m3": 1000.0,
    "dt_s": 0.001,
    "element_evaluation": "joint",
    "fiber_direction_reconstruction": "toolkit-physical-coordinate-u-v-v1",
    "fiber_sampling": "gp_direct_rule",
    "formulation": "hex8_local_pressure_p0_condensed_logj",
    "generalized_alpha": {
        "alpha_f": 0.4,
        "alpha_m": 0.2,
        "beta": 0.36,
        "gamma": 0.7,
        "stage_contract": "simula-source-matched-v1",
    },
    "linear_solver_profile": "fgmres-gamg-rigid-rebuild",
    "load_horizon_s": 1.0,
    "local_pressure_bulk_modulus_pa": 1_000_000.0,
    "local_pressure_volume_law": "linear-reference-volume-mean-log-j-v1",
    "mass": "consistent_q1_hex8",
    "material_eta_pa_s": 100.0,
    "mesh_topology": "closed_multiblock_disk",
    "point_sampling": "hex8_reference_isoparametric",
    "t_end_s": 0.32,
    "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
}
STEP0B_PREFIX_DIAGNOSTIC_COMPLETION_RECORD = {
    "completed_steps": 320,
    "converged": True,
    "expected_steps": 320,
    "function_domain_rejections_total": 0,
    "mpi_ranks": 8,
    "pre_solve_audits_passed": True,
    "status": "complete",
}
STEP0B_PREFIX_DIAGNOSTIC_DECISION = {
    "full_1s_status": "paused",
    "reason": (
        "The 0.20--0.32 s snap window already isolates a strong "
        "surface-resolution association, quantifies a consistent "
        "Robin-faceting mechanism, and leaves a local axial difference "
        "unresolved. Extending the same configuration cannot identify the "
        "remaining cause."
    ),
    "restart_condition": (
        "Define a new controlled hypothesis that changes the unresolved "
        "discrete surface comparison without changing or tuning the prescribed "
        "Robin law."
    ),
}
STEP0A_RETAINED_RESULT_IDENTITY = {
    "filename": "caseA_ga_local_pressure_rank8_t100.npz",
    "sha256": (
        "ba9b31ec533398be1f39fc9a898e72f77d9587c90f9b7d9e00ce91e4d2ae6a6c"
    ),
    "size_bytes": 16232720,
}
STEP0A_RETAINED_MANIFEST_IDENTITY = {
    "filename": "manifest.json",
    "sha256": (
        "db769328c9ba13079311cecaa33f1bbea0c1b1d9b33c5681118cf116b829b938"
    ),
    "size_bytes": 4080,
}
STEP0A_RETAINED_STDOUT_IDENTITY = {
    "filename": "caseA_ga_local_pressure_rank8_t100.stdout.txt",
    "sha256": (
        "beaef006e462bdfad4ce2e827620488aed633471d5907fd4debb3fab23187331"
    ),
    "size_bytes": 19467,
    "retained_in_repository": False,
    "reason": (
        "The raw transcript contains machine-local host and absolute-path text; "
        "only its external identity is retained here."
    ),
}
STEP0A_RETAINED_SEMANTIC_SHA256 = {
    "benchmark_identity": (
        "092852f3a3ff3fa02d97493374f63e141ebb39622ca78b6e864b924a9495c76a"
    ),
    "configuration": (
        "ca9115eb0f9e8640ebeb8a0c62161cf7df9c79051a50a4ac0db6c42b54cdfc3c"
    ),
    "completion": (
        "28fdae670bde507250c6a18936671df262d94a5beae5cc2df22508e705de7a0b"
    ),
    "ours_p0": (
        "31adddf571314d06e21861f13ed380a22246b8fffe64e70cadb92465b82fc2c4"
    ),
    "ours_p1": (
        "9abb1f154f9a5329c017f6dff340c00908e6199f748b29b568d97fa17b92d778"
    ),
    "mean_p0": (
        "b2dd9a1f3ea69f9b61f944a9f0c9678394ba50abe1a69c6138a9f067c3d5d973"
    ),
    "mean_p1": (
        "15a781e98e8b15c20db9b830ebc95442ab81e7ddbc141eea41b820a49f0cf21b"
    ),
}
STEP0A_RETAINED_TEAM_RED = {
    "p0": {
        "4C": 0.19760364288386148,
        "ambit": 0.1733154059485442,
        "carpentry": 0.2135461138894636,
        "cheart": 0.21265528065453948,
        "chimera": 0.1916964702132682,
        "comsol": 0.25163103019350447,
        "lifex": 0.18889640696579738,
        "simula": 1.4741510753649665,
        "simvascular_p1p1": 0.2981284488222884,
        "simvascular_p2": 0.2844364207240889,
    },
    "p1": {
        "4C": 0.2056416606201286,
        "ambit": 0.23357195772587572,
        "carpentry": 0.2735626233814094,
        "cheart": 0.21997041761934633,
        "chimera": 0.1992894246212024,
        "comsol": 0.24289587335477628,
        "lifex": 0.23888566289508656,
        "simula": 1.5450455180763714,
        "simvascular_p1p1": 0.3780582253566899,
        "simvascular_p2": 0.33243950884108947,
    },
}
STEP2B_FIGURE_RENDERER = (
    "examples/cardiac_benchmark/plot_step2b_case_b.py"
)
STEP2B_FIGURE_RENDERER_SHA256 = (
    "6d605acab25a3f300accb3fde07baf0479a0d74da9e896e7e75c1625c931be7d"
)
STEP2B_RAW_STDOUT = (
    "examples/cardiac_benchmark/results/"
    "step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt"
)
STEP2B_RAW_STDOUT_SHA256 = (
    "2bc7ddf633e1b905a1e8b42551ccf22eda4702f122cccd4370dd4c561a2a381c"
)
RETAINED_FIGURE_SPECS = {
    "docs/figures/case_a_comparison.svg": {
        "sha256": (
            "aff6a4c15b42967b6cd0f015de6a2362b41ba15e3e2bcf6d639564a08cb1da66"
        ),
        "size_bytes": 78136,
        "report": STEP0A_RETAINED_REPORT,
        "report_sha256": STEP0A_RETAINED_REPORT_SHA256,
        "title": "Benchmark 1, Case A — displacement comparison",
        "description": (
            "Six line charts compare CoupFE-Cardiac with the benchmark all-team "
            "mean from zero to one second. Rows show p0 and p1; columns show x, y, "
            "and z displacement in millimetres. Reported relative discrepancies "
            "are 0.334 at p0 and 0.502 at p1. Source: "
            "case_a_local_pressure_4x36x32_dt0p001.report.json, application "
            "016a4f9, Core e2f42ed, benchmark DOI 10.5281/zenodo.14260459. The "
            "Step 0A identity is explicitly labelled legacy-inferred. The chart "
            "is not a validation or pass claim."
        ),
        "visible_markers": (
            "closed t/core/radial 4×36×32 Hex8 Q1/P0 local pressure",
            "consistent-mass generalized-alpha",
            "Report: case_a_local_pressure_4x36x32_dt0p001.report.json",
            "app 016a4f9",
            "Core e2f42ed",
            "Step 0A legacy-inferred",
            "benchmark DOI 10.5281/zenodo.14260459",
            "CC-BY-4.0",
        ),
    },
    "docs/figures/archive/truncated_polar/case_a_comparison.svg": {
        "sha256": (
            "08b803fece500e0ebd6dc2b1e0b82736700d4fa465875f2ca204392515b8f431"
        ),
        "size_bytes": 77331,
        "report": (
            "examples/cardiac_benchmark/results/archive/truncated_polar/case_a/"
            "case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json"
        ),
        "report_sha256": (
            "103b67171bd15326983295940e267e0d0c7884481f05fbed2bbd1d3042ce2ddf"
        ),
        "title": "Benchmark 1, Case A — displacement comparison",
        "description": (
            "Six line charts compare CoupFE-Cardiac with the benchmark all-team "
            "mean from zero to one second. Rows show p0 and p1; columns show x, y, "
            "and z displacement in millimetres. Reported relative discrepancies "
            "are 0.506 at p0 and 0.679 at p1. Source: "
            "case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json, application "
            "6839c13, Core e2f42ed, benchmark DOI 10.5281/zenodo.14260459. The "
            "chart is not a validation or pass claim."
        ),
        "visible_markers": (
            "Report: case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json",
            "app 6839c13",
            "Core e2f42ed",
            "benchmark DOI 10.5281/zenodo.14260459",
            "CC-BY-4.0",
        ),
    },
    "docs/figures/archive/truncated_polar/case_b_comparison.svg": {
        "sha256": (
            "1814d84f78140606d3dfa748ffb7a601943765163faaa16875ac2eb2741d5a1b"
        ),
        "size_bytes": 78712,
        "report": (
            "examples/cardiac_benchmark/results/archive/truncated_polar/case_b/"
            "case_b_local_pressure_2x36x48_dt0p002.report.json"
        ),
        "report_sha256": (
            "d409ecaac0c0abf418fce3ab2f0549979d38e02b9d73381d56982f7fc4e3bf14"
        ),
        "title": "Benchmark 1, Case B — displacement comparison",
        "description": (
            "Six line charts compare CoupFE-Cardiac with the benchmark all-team "
            "mean from zero to one second. Rows show p0 and p1; columns show x, y, "
            "and z displacement in millimetres. Reported relative discrepancies "
            "are 0.551 at p0 and 0.654 at p1. Source: "
            "case_b_local_pressure_2x36x48_dt0p002.report.json, application "
            "e07993b, Core 454f73c, benchmark DOI 10.5281/zenodo.14260459. The "
            "chart is not a validation or pass claim."
        ),
        "visible_markers": (
            "Report: case_b_local_pressure_2x36x48_dt0p002.report.json",
            "app e07993b",
            "Core 454f73c",
            "benchmark DOI 10.5281/zenodo.14260459",
            "CC-BY-4.0",
        ),
    },
    "docs/figures/step2_case_b_comparison.svg": {
        "sha256": (
            "6a4c1b16fb2af6098d923380401fa1768ef85bd6062d670eee0114fdee694935"
        ),
        "size_bytes": 128994,
        "renderer": STEP2B_FIGURE_RENDERER,
        "report": (
            "examples/cardiac_benchmark/results/"
            "step2_case_b_std_kappa_2x20x17_dt0p001.report.json"
        ),
        "report_sha256": (
            "098e316daaea369a2a595cf43829d28597e53d2ff5a38cf32388e01c8dfa74aa"
        ),
        "title": "Benchmark 1, Step 2 Case B — displacement comparison",
        "description": (
            "Six line charts compare the full-cycle CoupFE-Cardiac result with "
            "the official ten-team range, all-team mean, and named Simula curve "
            "from zero to one second. Rows show p0 and p1; columns show x, y, and "
            "z displacement in millimetres. Full-history relative L2 error against "
            "the all-team mean is 9.8% overall, 9.1% at p0, and 10.9% at p1; "
            "benchmark-paper Eq. 21 RED is 28.3% at p0 and 35.3% at p1. The p1 z "
            "plateau has the opposite sign from all ten published curves. Source: "
            "step2_case_b_std_kappa_2x20x17_dt0p001.report.json, application "
            "e9b7d90 with an exact dirty-tree runtime manifest, Core e2f42ed, "
            "benchmark DOI 10.5281/zenodo.14260459. The chart is not a validation "
            "or pass claim."
        ),
        "visible_markers": (
            "paper Eq. 21 RED=28.3%/35.3%",
            "Report: step2_case_b_std_kappa_2x20x17_dt0p001.report.json",
            "app e9b7d90 dirty tree",
            "runtime manifest 6b96395761dd",
            "Core e2f42ed",
            "benchmark DOI 10.5281/zenodo.14260459",
            "CC-BY-4.0",
        ),
    },
    "docs/figures/step0b_tip_refine_full_cycle.svg": {
        "sha256": (
            "3cd45f62657d6c13574c5be6e64926e3ba0fb58dcdb119864576dc7ca0037a9b"
        ),
        "size_bytes": 426150,
        "renderer": TIP_REFINE_FIGURE_RENDERER,
        "report": TIP_REFINE_FULL_CYCLE_REPORT,
        "report_sha256": (
            "604fd83ca4c2773d2fca6dcb2488834925d81c9ba2328639aea6ba2f3c91d808"
        ),
        "title": (
            "Benchmark 1, Step 0 Case B - full-cycle displacement comparison"
        ),
        "description": (
            "Six line charts compare CoupFE-Cardiac 2x20x17 and 4x20x17 "
            "tip_refine=6.0 full-cycle trajectories with the retained local "
            "FEniCS Step 0B curves and the published ten-team envelope from "
            "zero to one second. Rows show p0 and p1; columns show x, y, and "
            "z displacement in millimetres. CoupFE-FEniCS comparisons use "
            "matched timestamps (FEniCS spans 0.001-0.999 s). The chart is "
            "not a validation or pass claim."
        ),
        "visible_markers": (
            "closed t/core/radial 2x20x17 and 4x20x17",
            "Hex8 Q1/P0 local pressure, tip_refine=6.0",
            "consistent-mass generalized-alpha",
            "Report: step0b_tip6p0_full_cycle_comparison.report.json",
            "app 2L ae2c2eb, 4L 2458e7c",
            "Core 454f73c",
            "benchmark DOI 10.5281/zenodo.14260459",
            "CC-BY-4.0",
        ),
    },
    "docs/figures/step2b_current_rerun_comparison.svg": {
        "sha256": (
            "7408e1324613f8f76ba389bd22a8fbb9d1b2e458af9257f46272074c2c0a2243"
        ),
        "size_bytes": 393630,
        "renderer": STEP2B_RERUN_FIGURE_RENDERER,
        "report": STEP2B_RERUN_FULL_CYCLE_REPORT,
        "report_sha256": (
            "3a4efd96db303a7dc0bdcf7d7e6c27ec6b91ed579d429b37c812138d1bac20cc"
        ),
        "title": (
            "Benchmark 1, Step 2 Case B - retained corrected-setup diagnostic"
        ),
        "description": (
            "Six line charts compare a retained CoupFE-Cardiac Step 2 Case B "
            "corrected-setup diagnostic (straight-wall geometry, physical "
            "frame, generalized-alpha Q1/P0 local pressure, tip_refine=6.0) "
            "with the official ten-team Step 2 Case B envelope from zero to "
            "one second. Rows show p0 and p1; columns show x, y, and z "
            "displacement in millimetres. The archive records application "
            "97d4474 and Core 454f73c and is not a current-release result. "
            "The chart is a provenance-incomplete diagnostic, not a "
            "reproduction, validation, or pass claim."
        ),
        "visible_markers": (
            "closed t/core/radial 2x20x17 Hex8 Q1/P0 local pressure, "
            "tip_refine=6.0",
            "consistent-mass generalized-alpha",
            "retained corrected-setup diagnostic",
            "Report: step2b_current_rerun_comparison.report.json",
            "app 97d4474",
            "Core 454f73c",
            "benchmark DOI 10.5281/zenodo.14260459",
            "CC-BY-4.0",
        ),
    },
}

HISTORICAL_BENCHMARK_RESULT_SOURCES = {
    "examples/cardiac_benchmark/results/README.md",
    TRUNCATED_POLAR_ARCHIVE_README,
    RETAINED_RESULT_JSON,
    RETAINED_RESULT_STDOUT,
}
STEP2B_BENCHMARK_RESULT_SOURCES = {
    "examples/cardiac_benchmark/results/"
    "step2_case_b_std_kappa_2x20x17_dt0p001.report.json",
    STEP2B_RAW_STDOUT,
}
STEP0A_RETAINED_BENCHMARK_RESULT_SOURCES = {
    STEP0A_RETAINED_REPORT,
}
STEP0B_PREFIX_DIAGNOSTIC_RESULT_SOURCES = {
    STEP0B_PREFIX_DIAGNOSTIC_REPORT,
}
BENCHMARK_RESULT_SOURCES = (
    HISTORICAL_BENCHMARK_RESULT_SOURCES
    | STEP0A_RETAINED_BENCHMARK_RESULT_SOURCES
    | STEP0B_PREFIX_DIAGNOSTIC_RESULT_SOURCES
    | STEP2B_BENCHMARK_RESULT_SOURCES
    | {TIP_REFINE_FULL_CYCLE_REPORT}
    | {STEP2B_RERUN_FULL_CYCLE_REPORT}
    | {
        (
            f"{TRUNCATED_POLAR_ARCHIVE_DIRECTORY}/"
            f"case_{spec['case'].lower()}/{spec[field]}"
        )
        for spec in TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS
        for field in ("report", "log")
    }
)
BENCHMARK_SOURCES = {
    f"examples/cardiac_benchmark/{name}"
    for name in {
        "README.md",
        ".gitignore",
        "activation.py",
        "benchmark_parameters.py",
        "boundary_audit.py",
        "compare_fenics_case_b.py",
        "compare_mpi_rank_gate.py",
        "compare_step2b_case_b.py",
        "compare_tip_refine_full_cycle.py",
        "consistent_mass.py",
        "diagnose.py",
        "fiber_crosscheck.py",
        "generalized_alpha.py",
        "geometry.py",
        "local_pressure.py",
        "material.py",
        "mesh_quality.py",
        "newmark.py",
        "post.py",
        "plot_retained_comparisons.py",
        "plot_step2b_case_b.py",
        "plot_step2b_current_rerun.py",
        "pressure.py",
        "publish_step0_comparison.py",
        "result_io.py",
        "robin.py",
        "run.py",
        "run_mpi.py",
        "sampling.py",
        "solver.py",
        "structural_directions.py",
        "tbar_laplace.py",
        "viscous_evidence.py",
        "distributed_local_pressure.py",
        "distributed_mass.py",
        "distributed_material.py",
        "distributed_solver.py",
        "fenics_case_b_reference_hashes.example.json",
        "step2b_case_b_reference_hashes.json",
        "step2b_case_b_runtime_source_hashes.json",
    }
} | BENCHMARK_RESULT_SOURCES
MPI_SOURCES = {
    f"examples/mpi_smoke/{name}"
    for name in {
        "distributed_cardiac_dynamics.py",
        "distributed_cardiac_generalized_alpha.py",
        "distributed_cardiac_passive.py",
        "distributed_cardiac_pressure.py",
        "distributed_cardiac_scaling.py",
        "distributed_cardiac_viscous.py",
        "README.md",
    }
}
TEST_SOURCES = {
    f"tests/{name}"
    for name in {
        "conftest.py",
        "test_cardiac_closed_geometry.py",
        "test_cardiac_consistent_mass.py",
        "test_cardiac_fast.py",
        "test_cardiac_distributed_local_pressure.py",
        "test_cardiac_local_pressure.py",
        "test_cardiac_mpi.py",
        "test_cardiac_mpi_closed_configuration.py",
        "test_cardiac_mpi_companion.py",
        "test_cardiac_mpi_rank_gate.py",
        "test_cardiac_fenics_comparison.py",
        "test_cardiac_generalized_alpha.py",
        "test_cardiac_reporting.py",
        "test_cardiac_sampling.py",
        "test_cardiac_slow.py",
        "test_cardiac_step0_publisher.py",
        "test_cardiac_step2_case_b.py",
        "test_cardiac_step2b_comparison.py",
        "test_cardiac_structural_directions.py",
        "test_cardiac_tbar_laplace.py",
        "test_cardiac_tip_refine_full_cycle_comparison.py",
        "test_cardiac_viscous_evidence.py",
    }
}
DOC_SOURCES = {
    f"docs/{name}"
    for name in {
        "API.md",
        "BENCHMARK_COMPARISON.md",
        "BENCHMARK_REPRODUCTION_STATUS.md",
        "BENCHMARK_TEST_MATRIX.md",
        "CASE_A_STATUS.md",
        "CASE_SPECIFICATIONS.md",
        "CASE_B_DEBUGGING_POSTMORTEM.md",
        "CASE_B_FENICS_COMPARISON.md",
        "CASE_B_MESH_ERROR_LAYERS.md",
        "CASE_B_MPI_RANK_GATE.md",
        "CASE_B_ROBIN_NODAL_SMOOTHED_NORMAL_DIAGNOSTIC.md",
        "CASE_B_STATUS.md",
        "CONTROLLED_BENCHMARK_RUNS.md",
        "LICENSE.md",
        "MESH_REFINEMENT_GUIDE.md",
        "PHYSICAL_FRAME_RECONSTRUCTION.md",
        "RELEASE_CHECKLIST.md",
        "STEP2_CASE_B_REPRODUCTION_LOG.md",
        "lessons_learned.md",
        "figures/README.md",
        "figures/archive/truncated_polar/README.md",
        "figures/archive/truncated_polar/case_a_comparison.svg",
        "figures/archive/truncated_polar/case_b_comparison.svg",
        "figures/case_a_comparison.svg",
        "figures/step0b_tip_refine_full_cycle.svg",
        "figures/step2_case_b_comparison.svg",
        "figures/step2b_current_rerun_comparison.svg",
    }
}
SKILL_SOURCES = {
    "skills/SKILL.md",
    "skills/cardiac.md",
}
LICENSE_SOURCES = {
    "LICENSES/CC-BY-4.0.txt",
    "LICENSES/Reidmen-cardiac_benchmark_toolkit-MIT.txt",
}
GITHUB_SOURCES = {
    ".github/scripts/check_release_artifacts.py",
    ".github/scripts/check_runtime_core.py",
    ".github/workflows/ci.yml",
    ".github/workflows/extended.yml",
}
ROOT_SOURCES = {
    ".gitattributes",
    ".gitignore",
    "CONTRIBUTING.md",
    "LICENSE",
    "MANIFEST.in",
    "NOTICE",
    "README.md",
    "THIRD_PARTY_NOTICES.md",
    "pyproject.toml",
    "setup.sh",
}
PUBLIC_SOURCE_INVENTORIES = {
    ".github": GITHUB_SOURCES,
    "LICENSES": LICENSE_SOURCES,
    "docs": DOC_SOURCES,
    "examples": BENCHMARK_SOURCES
    | MPI_SOURCES
    | {"examples/README.md", "examples/REFERENCES.md"},
    "skills": SKILL_SOURCES,
    "tests": TEST_SOURCES,
}
PUBLIC_RELEASE_FILES = ROOT_SOURCES | set().union(
    *PUBLIC_SOURCE_INVENTORIES.values()
)
SDIST_METADATA_FILES = {
    f"coupfe_cardiac.egg-info/{name}"
    for name in {
        "PKG-INFO",
        "SOURCES.txt",
        "dependency_links.txt",
        "requires.txt",
        "top_level.txt",
    }
}
EXPECTED_SDIST_FILES = PUBLIC_RELEASE_FILES | SDIST_METADATA_FILES | {
    "PKG-INFO",
    "setup.cfg",
}
REQUIRED_SDIST_FILES = EXPECTED_SDIST_FILES

FORBIDDEN_PARTS = {
    ".claude",
    ".codex",
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "notes",
}
FORBIDDEN_NAMES = {
    "AGENTS.md",
    "CLAUDE.md",
}
FORBIDDEN_NAME_PATTERNS = {
    "CURRENT_STATUS*.md",
    "HANDOFF*.md",
    "NOTE_FROM_*.md",
    "NOTE_TO_*.md",
    "PORT_PROMPT*.md",
    "RELEASE_READINESS*.md",
    "REVIEW_NOTES*.md",
    "TRUST_NET_FINDINGS*.md",
    "*_for_claude*.md",
    "*_for_gpt*.md",
    "*_for_kimi*.md",
    "*.tar.bz2",
    "*.tar.gz",
    "*.tar.xz",
    "*.tgz",
}
FORBIDDEN_SUFFIXES = {
    ".7z",
    ".a",
    ".bin",
    ".bz2",
    ".cab",
    ".class",
    ".dll",
    ".dmg",
    ".doc",
    ".docx",
    ".dylib",
    ".exe",
    ".gif",
    ".gz",
    ".h5",
    ".hdf5",
    ".iso",
    ".jar",
    ".jpeg",
    ".jpg",
    ".lib",
    ".lz",
    ".lz4",
    ".mod",
    ".nbc",
    ".nbi",
    ".npz",
    ".npy",
    ".o",
    ".obj",
    ".pdf",
    ".pickle",
    ".pkl",
    ".png",
    ".ppt",
    ".pptx",
    ".pyc",
    ".pyd",
    ".pyo",
    ".rar",
    ".so",
    ".svg",
    ".tar",
    ".tgz",
    ".vtk",
    ".vtu",
    ".webp",
    ".whl",
    ".xls",
    ".xlsx",
    ".xz",
    ".zip",
    ".zst",
    ".zstd",
}
APPROVED_SVG_PATHS = set(RETAINED_FIGURE_SPECS)
TEXT_SUFFIXES = {
    "",
    ".cfg",
    ".csv",
    ".f90",
    ".for",
    ".in",
    ".ini",
    ".json",
    ".md",
    ".py",
    ".rst",
    ".sh",
    ".svg",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}


def _validate_names(names: list[str], artifact: Path) -> None:
    duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicates:
        raise SystemExit(f"{artifact.name} contains duplicate entries: {duplicates}")

    unsafe = []
    for name in names:
        path = PurePosixPath(name)
        normalized_name = name[:-1] if name.endswith("/") else name
        raw_parts = normalized_name.split("/")
        if (
            not name
            or "\\" in name
            or re.match(r"^[A-Za-z]:", name)
            or path.is_absolute()
            or ".." in path.parts
            or any(part in {"", ".", ".."} for part in raw_parts)
        ):
            unsafe.append(name)
    if unsafe:
        raise SystemExit(f"{artifact.name} contains unsafe paths: {sorted(unsafe)}")


def _retired_terms() -> tuple[str, ...]:
    # Concatenation prevents this guard from rejecting its own source.
    return (
        "p" + "q" + "p",
        "material" + "_" + "mixed",
        "--" + "mixed",
    )


def _is_forbidden_path(name: str) -> bool:
    path = PurePosixPath(name)
    parts = tuple(part.casefold() for part in path.parts)
    basename = path.name.casefold()
    suffix = path.suffix.casefold()
    return bool(
        set(parts).intersection(part.casefold() for part in FORBIDDEN_PARTS)
        or basename in {item.casefold() for item in FORBIDDEN_NAMES}
        or any(
            fnmatchcase(basename, pattern.casefold())
            for pattern in FORBIDDEN_NAME_PATTERNS
        )
        or any(term in part for term in _retired_terms() for part in parts)
        or (suffix in FORBIDDEN_SUFFIXES and name not in APPROVED_SVG_PATHS)
    )


def _reject_forbidden_files(names: set[str], artifact: Path) -> None:
    rejected = sorted(name for name in names if _is_forbidden_path(name))
    if rejected:
        raise SystemExit(f"{artifact.name} contains forbidden entries: {rejected}")


def _require_files(available: set[str], required: set[str], artifact: Path) -> None:
    missing = sorted(required - available)
    if missing:
        raise SystemExit(f"{artifact.name} is missing required files: {missing}")


def _validate_exact_subtree(
    names: set[str], artifact: Path, prefix: str, expected: set[str]
) -> None:
    """Require a reviewed public subtree to match its static inventory."""
    present = {
        name
        for name in names
        if PurePosixPath(name).parts
        and PurePosixPath(name).parts[0] == prefix
    }
    missing = sorted(expected - present)
    unexpected = sorted(present - expected)
    if missing or unexpected:
        raise SystemExit(
            f"{artifact.name} public {prefix} inventory mismatch: "
            f"missing={missing}, unexpected={unexpected}"
        )


def _validate_exact_root_files(
    names: set[str], artifact: Path, expected: set[str]
) -> None:
    """Require reviewed top-level files and reject unreviewed additions."""
    present = {
        name for name in names if len(PurePosixPath(name).parts) == 1
    }
    missing = sorted(expected - present)
    unexpected = sorted(present - expected)
    if missing or unexpected:
        raise SystemExit(
            f"{artifact.name} public root inventory mismatch: "
            f"missing={missing}, unexpected={unexpected}"
        )


def _validate_public_source_inventories(
    names: set[str],
    artifact: Path,
    *,
    root_files: set[str],
    allowed_generated_subtrees: Optional[set[str]] = None,
) -> None:
    """Validate every reviewed subtree and the top-level boundary."""
    _validate_exact_root_files(names, artifact, root_files)
    for prefix, expected in PUBLIC_SOURCE_INVENTORIES.items():
        _validate_exact_subtree(names, artifact, prefix, expected)

    allowed_prefixes = set(PUBLIC_SOURCE_INVENTORIES)
    allowed_prefixes.update(allowed_generated_subtrees or set())
    unexpected_prefixes = sorted(
        {
            path.parts[0]
            for name in names
            if len((path := PurePosixPath(name)).parts) > 1
            and path.parts[0] not in allowed_prefixes
        }
    )
    if unexpected_prefixes:
        raise SystemExit(
            f"{artifact.name} contains unreviewed top-level subtrees: "
            f"{unexpected_prefixes}"
        )


def _require_markers(
    payloads: dict[str, bytes],
    required: dict[str, tuple[str, ...]],
    artifact: Path,
) -> None:
    """Require release-critical license boundaries and provenance records."""
    for name, markers in required.items():
        if name not in payloads:
            raise SystemExit(f"{artifact.name} is missing legal record {name}")
        text = payloads[name].decode("utf-8")
        missing = [marker for marker in markers if marker not in text]
        if missing:
            raise SystemExit(
                f"{artifact.name}:{name} is missing required legal markers: "
                f"{missing}"
            )


def _load_strict_json(payload: bytes, name: str, artifact: Path):
    """Decode JSON while rejecting duplicate keys and non-finite constants."""

    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate key {key!r}")
            result[key] = value
        return result

    def finite_constant(value):
        raise ValueError(f"non-finite numeric constant {value!r}")

    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=finite_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise SystemExit(
            f"{artifact.name}:{name} is not strict finite JSON: {exc}"
        ) from exc


def _finite_number(value) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _validate_retained_result(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Validate the immutable historical reduced result and its audit trail."""
    required = {
        RETAINED_RESULT_JSON,
        TRUNCATED_POLAR_ARCHIVE_README,
        RETAINED_RESULT_STDOUT,
    }
    missing = sorted(required - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing retained-result inputs: {missing}"
        )

    record = _load_strict_json(
        payloads[RETAINED_RESULT_JSON], RETAINED_RESULT_JSON, artifact
    )
    if not isinstance(record, dict):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} must contain a JSON object"
        )
    if hashlib.sha256(payloads[RETAINED_RESULT_JSON]).hexdigest() != (
        RETAINED_RESULT_JSON_SHA256
    ):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} does not match its "
            "reviewed historical SHA-256"
        )

    exact_top_level = {
        "schema": "coupfe-cardiac-retained-result-v1",
        "evidence_label": "checked reduced execution demonstration",
        "case": "A",
        "license": "CC-BY-4.0",
        "command": list(RETAINED_RESULT_COMMAND),
        "claim_boundary": list(RETAINED_RESULT_CLAIM_BOUNDARY),
    }
    mismatched = {
        key: record.get(key)
        for key, expected in exact_top_level.items()
        if record.get(key) != expected
    }
    if mismatched:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} has altered retained-result "
            f"scope or command: {sorted(mismatched)}"
        )

    source = record.get("source")
    if not isinstance(source, dict):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} is missing source provenance"
        )
    expected_source_identity = {
        "app_revision": RETAINED_RESULT_APP_REF,
        "app_tree_state": "clean",
        "core_revision": HISTORICAL_RETAINED_CORE_REF,
        "core_tree_state": "clean",
        "core_url": PUBLIC_CORE_URL,
    }
    identity_mismatch = {
        key: source.get(key)
        for key, expected in expected_source_identity.items()
        if source.get(key) != expected
    }
    if identity_mismatch:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} has altered source "
            f"provenance: {sorted(identity_mismatch)}"
        )

    recorded_hashes = source.get("source_files_sha256")
    if recorded_hashes != HISTORICAL_RETAINED_RESULT_SOURCE_SHA256:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} must record the exact "
            "historical load-bearing source hashes"
        )

    artifacts = record.get("artifacts")
    expected_artifacts = {
        "generated_npz_distributed": False,
        "generated_npz_sha256": RETAINED_RESULT_NPZ_SHA256,
        "normalized_stdout": PurePosixPath(RETAINED_RESULT_STDOUT).name,
        "normalized_stdout_sha256": RETAINED_RESULT_STDOUT_SHA256,
        "repeat_run_byte_identical": True,
    }
    if artifacts != expected_artifacts:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} has altered artifact "
            "provenance"
        )
    stdout_hash = hashlib.sha256(payloads[RETAINED_RESULT_STDOUT]).hexdigest()
    if stdout_hash != RETAINED_RESULT_STDOUT_SHA256:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_STDOUT} does not match its "
            "reviewed SHA-256"
        )
    stdout = payloads[RETAINED_RESULT_STDOUT].decode("utf-8")
    expected_command_line = "$ " + " ".join(RETAINED_RESULT_COMMAND) + "\n"
    if not stdout.startswith(expected_command_line):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_STDOUT} does not begin with the "
            "recorded command"
        )
    for marker in (
        "CASE A: tau_peak=118042.8Pa  p_peak=0.0Pa",
        "step  500 t=1.000s",
        "saved -> case_a_reduced.npz",
        "finished 500/500 steps and wrote a converged result archive",
    ):
        if marker not in stdout:
            raise SystemExit(
                f"{artifact.name}:{RETAINED_RESULT_STDOUT} is missing "
                f"completion marker {marker!r}"
            )

    configuration = record.get("configuration")
    if not isinstance(configuration, dict):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} is missing configuration"
        )
    expected_time = {
        "completed_steps": 500,
        "dt_s": 0.002,
        "expected_steps": 500,
        "integrator": "be",
        "t_end_s": 1.0,
    }
    expected_mesh = {
        "apex_offset": 0.2,
        "degrees_of_freedom": 72,
        "elements": 8,
        "n_mu": 2,
        "n_t": 1,
        "n_theta": 4,
        "nodes": 24,
    }
    expected_labels = {
        "fiber_sampling": "cg1_gram_schmidt",
        "formulation": "hex8_fbar",
        "point_sampling": "global_delaunay_tetra",
        "viscous_rate": "backward_difference",
    }
    label_mismatch = {
        key: configuration.get(key)
        for key, expected in expected_labels.items()
        if configuration.get(key) != expected
    }
    if (
        configuration.get("time") != expected_time
        or configuration.get("mesh") != expected_mesh
        or label_mismatch
    ):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} does not describe the "
            "reviewed reduced configuration"
        )

    history = record.get("retained_history")
    if not isinstance(history, dict):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} is missing retained history"
        )
    if history.get("sample_interval_s") != 0.02 or history.get(
        "source_step_stride"
    ) != 10:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} has altered history sampling"
        )
    samples = history.get("samples")
    if not isinstance(samples, list) or len(samples) != 51:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} must retain exactly 51 "
            "history samples"
        )

    expected_sample_keys = {
        "active_tension_pa",
        "pressure_pa",
        "time_s",
        "u0_m",
        "u1_m",
    }
    for index, sample in enumerate(samples):
        if not isinstance(sample, dict) or set(sample) != expected_sample_keys:
            raise SystemExit(
                f"{artifact.name}:{RETAINED_RESULT_JSON} history sample {index} "
                "has an unexpected schema"
            )
        expected_time_value = index * 0.02
        if not _finite_number(sample["time_s"]) or not math.isclose(
            float(sample["time_s"]),
            expected_time_value,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise SystemExit(
                f"{artifact.name}:{RETAINED_RESULT_JSON} history sample {index} "
                "has an inconsistent time"
            )
        for field in ("active_tension_pa", "pressure_pa"):
            if not _finite_number(sample[field]):
                raise SystemExit(
                    f"{artifact.name}:{RETAINED_RESULT_JSON} history sample "
                    f"{index} has non-finite {field}"
                )
        if sample["active_tension_pa"] < 0.0 or sample["pressure_pa"] != 0.0:
            raise SystemExit(
                f"{artifact.name}:{RETAINED_RESULT_JSON} history sample {index} "
                "is inconsistent with retained Case A loading"
            )
        for field in ("u0_m", "u1_m"):
            vector = sample[field]
            if (
                not isinstance(vector, list)
                or len(vector) != 3
                or not all(_finite_number(value) for value in vector)
            ):
                raise SystemExit(
                    f"{artifact.name}:{RETAINED_RESULT_JSON} history sample "
                    f"{index} has invalid {field}"
                )

    result = record.get("result")
    if not isinstance(result, dict) or result.get("converged") is not True:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} is not a converged result"
        )
    peak_index = max(
        range(len(samples)), key=lambda index: samples[index]["active_tension_pa"]
    )
    peak = samples[peak_index]
    summary_pairs = {
        "peak_active_tension_pa": peak["active_tension_pa"],
        "peak_pressure_pa": peak["pressure_pa"],
        "peak_time_s": peak["time_s"],
        "u0_at_peak_m": peak["u0_m"],
        "u1_at_peak_m": peak["u1_m"],
        "u0_final_m": samples[-1]["u0_m"],
        "u1_final_m": samples[-1]["u1_m"],
    }
    inconsistent_summary = sorted(
        key for key, expected in summary_pairs.items() if result.get(key) != expected
    )
    if result.get("peak_step") != 240 or peak_index != 24 or inconsistent_summary:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} summary is inconsistent "
            f"with its retained history: {inconsistent_summary}"
        )
    for field in (
        "max_abs_nodal_component_at_peak_m",
        "max_nodal_norm_at_peak_m",
    ):
        if not _finite_number(result.get(field)) or result[field] <= 0.0:
            raise SystemExit(
                f"{artifact.name}:{RETAINED_RESULT_JSON} has invalid summary "
                f"field {field!r}"
            )
    diagnostics = result.get("diagnostics_at_peak")
    if (
        not isinstance(diagnostics, dict)
        or not diagnostics
        or not all(_finite_number(value) for value in diagnostics.values())
        or not (
            diagnostics.get("centroid_det_f_min", math.inf)
            < 1.0
            < diagnostics.get("centroid_det_f_max", -math.inf)
        )
    ):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_RESULT_JSON} must retain the finite, "
            "broad deformation diagnostic behind its coarse-mesh warning"
        )

    _require_markers(
        payloads,
        {
            TRUNCATED_POLAR_ARCHIVE_README: (
                "Truncated-polar result archive",
                "non-benchmark geometry",
                "not current Benchmark 1 validation evidence",
                "Historical reduced Case A",
                "global Delaunay-tetra policy",
                "The generated CoupFE NPZ archives are not committed",
                "CC BY 4.0",
            ),
        },
        artifact,
    )


def _require_exact_keys(value, expected: set[str], description: str, artifact: Path):
    if not isinstance(value, dict) or set(value) != expected:
        actual = sorted(value) if isinstance(value, dict) else type(value).__name__
        raise SystemExit(
            f"{artifact.name}:{description} has an unexpected schema: {actual}"
        )


def _truncated_polar_archive_path(filename: str) -> str:
    if filename.startswith("case_a_"):
        case_directory = "case_a"
    elif filename.startswith("case_b_"):
        case_directory = "case_b"
    else:
        raise ValueError(f"not a Case A/B truncated-polar artifact: {filename}")
    return f"{TRUNCATED_POLAR_ARCHIVE_DIRECTORY}/{case_directory}/{filename}"


def _validate_no_absolute_paths(value, description: str, artifact: Path) -> None:
    if isinstance(value, str):
        if Path(value).is_absolute() or PureWindowsPath(value).is_absolute():
            raise SystemExit(
                f"{artifact.name}:{description} contains an absolute filesystem path"
            )
        return
    if isinstance(value, list):
        for item in value:
            _validate_no_absolute_paths(item, description, artifact)
    elif isinstance(value, dict):
        for item in value.values():
            _validate_no_absolute_paths(item, description, artifact)


def _integer_number(value, *, minimum: int = 0) -> bool:
    return (
        _finite_number(value)
        and float(value).is_integer()
        and int(value) >= minimum
    )


def _finite_vector(value, length: int) -> bool:
    return (
        isinstance(value, list)
        and len(value) == length
        and all(_finite_number(component) for component in value)
    )


def _expected_reference_team_files(reference_case: str) -> list[dict]:
    expected = []
    for team in REFERENCE_TEAM_ORDER:
        digest = REFERENCE_TEAM_SHA256[reference_case][team]
        start = 0.001 if team in {"ambit", "simula"} else 0.0
        end = 0.999 if team == "simula" else 1.0
        if reference_case == "step_0B" and team == "4C":
            end = 0.999999999999906
        expected.append(
            {
                "team": team,
                "filename": (
                    f"monoventricular_nonblinded_{reference_case}_group_"
                    f"{team}.pickle"
                ),
                "sha256": digest,
                "size_bytes": REFERENCE_TEAM_SIZE_BYTES[team],
                "source_sample_count": 101,
                "source_time_start": start,
                "source_time_end": end,
            }
        )
    return expected


def _expected_reference_selection(reference_case: str) -> dict:
    alias_filename = (
        f"monoventricular_nonblinded_{reference_case}_group_simvascular.pickle"
    )
    target_filename = (
        f"monoventricular_nonblinded_{reference_case}_group_"
        "simvascular_p2.pickle"
    )
    return {
        "policy": REFERENCE_SELECTION_POLICY,
        "upstream_figures_py_identity": REFERENCE_FIGURES_PY_IDENTITY,
        "upstream_manifest_variable": REFERENCE_MANIFEST_VARIABLES[reference_case],
        "selected_count": len(REFERENCE_TEAM_ORDER),
        "selected_files": [
            f"monoventricular_nonblinded_{reference_case}_group_{team}.pickle"
            for team in REFERENCE_TEAM_ORDER
        ],
        "excluded_aliases": [
            {
                "filename": alias_filename,
                "sha256": REFERENCE_EXCLUDED_ALIAS_SHA256[reference_case],
                "size_bytes": 8651,
                "identical_to_selected_filename": target_filename,
                "reason": (
                    "Byte-identical base-name alias of selected SimVascular P2; "
                    "excluded because upstream figures.py does not select it."
                ),
            }
        ],
    }


def _validate_curve(value, description: str, artifact: Path) -> list[list[float]]:
    if (
        not isinstance(value, list)
        or len(value) != 101
        or not all(_finite_vector(row, 3) for row in value)
    ):
        raise SystemExit(
            f"{artifact.name}:{description} must be a finite 101-by-3 curve"
        )
    return value


def _relative_discrepancy(curve, mean) -> float:
    terms = []
    for ours, reference in zip(curve, mean):
        numerator = math.sqrt(
            sum((float(left) - float(right)) ** 2 for left, right in zip(ours, reference))
        )
        denominator = math.sqrt(sum(float(value) ** 2 for value in reference))
        terms.append(numerator / (denominator + 1.0e-30))
    return sum(terms) / len(terms)


def _validate_current_source_hashes(
    payloads: dict[str, bytes], artifact: Path, app_ref: str = CURRENT_RESULT_APP_REF
) -> None:
    expected_hashes = CURRENT_RESULT_SOURCE_CHECKPOINTS.get(app_ref)
    if expected_hashes is None:
        raise SystemExit(
            f"{artifact.name} has no reviewed source hashes for checkpoint {app_ref}"
        )
    missing = sorted(set(expected_hashes) - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing result source inputs for checkpoint "
            f"{app_ref}: {missing}"
        )
    changed = sorted(
        name
        for name, expected in expected_hashes.items()
        if hashlib.sha256(payloads[name]).hexdigest() != expected
    )
    if changed:
        raise SystemExit(
            f"{artifact.name} result sources differ from checkpoint "
            f"{app_ref}: {changed}"
        )


def _validate_result_source_checkpoints(
    payloads: dict[str, bytes], artifact: Path
) -> str:
    """Require packaged executable sources to match one reviewed checkpoint."""
    if not CURRENT_RESULT_SOURCE_CHECKPOINTS:
        raise SystemExit("release guard has no reviewed result source checkpoints")

    paths = None
    for app_ref, hashes in CURRENT_RESULT_SOURCE_CHECKPOINTS.items():
        if (
            not isinstance(app_ref, str)
            or re.fullmatch(r"[0-9a-f]{40}", app_ref) is None
        ):
            raise SystemExit(
                f"release guard has invalid result source checkpoint {app_ref!r}"
            )
        if not isinstance(hashes, dict) or not hashes:
            raise SystemExit(
                f"release guard has no source hashes for checkpoint {app_ref}"
            )
        if paths is None:
            paths = set(hashes)
        elif set(hashes) != paths:
            raise SystemExit(
                "release guard result source checkpoints cover different files"
            )
        malformed = sorted(
            name
            for name, digest in hashes.items()
            if not isinstance(name, str)
            or not isinstance(digest, str)
            or re.fullmatch(r"[0-9a-f]{64}", digest) is None
        )
        if malformed:
            raise SystemExit(
                f"release guard has malformed source hashes for checkpoint "
                f"{app_ref}: {malformed}"
            )

    assert paths is not None
    missing = sorted(paths - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing reviewed result source inputs: {missing}"
        )
    matches = [
        app_ref
        for app_ref, hashes in CURRENT_RESULT_SOURCE_CHECKPOINTS.items()
        if all(
            hashlib.sha256(payloads[name]).hexdigest() == digest
            for name, digest in hashes.items()
        )
    ]
    if len(matches) != 1:
        raise SystemExit(
            f"{artifact.name} result sources must match exactly one reviewed "
            f"checkpoint; matches={matches}"
        )
    return matches[0]


def _validate_current_release_source_hashes(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Bind current executable code without rewriting historical provenance."""
    malformed = sorted(
        name
        for name, digest in CURRENT_RELEASE_SOURCE_SHA256.items()
        if not isinstance(name, str)
        or not isinstance(digest, str)
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    )
    if malformed:
        raise SystemExit(
            f"release guard has malformed current source hashes: {malformed}"
        )
    missing = sorted(set(CURRENT_RELEASE_SOURCE_SHA256) - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing current executable source inputs: {missing}"
        )
    changed = sorted(
        name
        for name, digest in CURRENT_RELEASE_SOURCE_SHA256.items()
        if hashlib.sha256(payloads[name]).hexdigest() != digest
    )
    if changed:
        raise SystemExit(
            f"{artifact.name} current executable sources differ from reviewed "
            f"release code: {changed}"
        )


def _validate_current_reporting_source_hashes(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Bind reporting-only code separately from simulation checkpoints."""
    missing = sorted(set(CURRENT_REPORTING_SOURCE_SHA256) - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing current reporting source inputs: {missing}"
        )
    changed = sorted(
        name
        for name, digest in CURRENT_REPORTING_SOURCE_SHA256.items()
        if hashlib.sha256(payloads[name]).hexdigest() != digest
    )
    if changed:
        raise SystemExit(
            f"{artifact.name} current reporting sources differ from reviewed "
            f"selection code: {changed}"
        )


def _validate_current_report_configuration(
    configuration: dict, spec: dict, report_name: str, artifact: Path
) -> None:
    expected_configuration_keys = {
        "apex_offset_rad",
        "dt_s",
        "fiber_sampling",
        "flip_helix",
        "formulation",
        "integrator",
        "mesh",
        "model_parameters",
        "nonlinear_solver",
        "point_sampling",
        "sampling_points",
        "t_end_s",
        "viscous_rate",
    }
    if spec.get("method_metadata", False):
        expected_configuration_keys |= {
            "fiber_sampling_option",
            "isotropic",
            "mass_representation",
            "method_metadata_origin",
            "material_eta_pa_s",
            "parameter_variant",
            "tbar",
            "viscous_term_active",
        }
    _require_exact_keys(
        configuration,
        expected_configuration_keys,
        f"{report_name}:result.configuration",
        artifact,
    )
    expected_labels = {
        "apex_offset_rad": 0.2,
        "dt_s": spec["dt_s"],
        "fiber_sampling": "cg1_gram_schmidt",
        "flip_helix": True,
        "formulation": spec["formulation"],
        "integrator": "be",
        "nonlinear_solver": spec["solver"],
        "point_sampling": "hex8_reference_isoparametric",
        "t_end_s": 1.0,
        "viscous_rate": "backward_difference",
    }
    if spec.get("method_metadata", False):
        expected_labels.update(
            {
                "fiber_sampling_option": "cg1",
                "isotropic": False,
                "mass_representation": "lumped_row_sum",
                "method_metadata_origin": (
                    "reviewed-predecessor-source-checkpoint"
                ),
                "material_eta_pa_s": 100.0,
                "parameter_variant": "benchmark_eta",
                "tbar": {
                    "definition": "analytic_parametric",
                    "metadata_filename": "",
                    "metadata_schema": "",
                    "metadata_sha256": "",
                    "source_filename": "",
                    "source_sha256": "",
                },
                "viscous_term_active": True,
            }
        )
    altered = sorted(
        key for key, expected in expected_labels.items()
        if configuration.get(key) != expected
    )
    if altered or configuration.get("mesh") != spec["mesh"]:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered run configuration: "
            f"{altered or ['mesh']}"
        )

    common_model = {
        "base_robin_damping_pa_s_m": 5000.0,
        "base_robin_stiffness_pa_m": 100000.0,
        "density_kg_m3": 1000.0,
        "epicardial_robin_damping_pa_s_m": 5000.0,
        "epicardial_robin_stiffness_pa_m": 100000000.0,
        "mesh_perturbation_std_m": 0.0,
    }
    if spec["formulation"] == "hex8_local_pressure_p0_condensed_logj":
        expected_model = {
            **common_model,
            "local_pressure_bulk_modulus_pa": 1000000.0,
            "material_kappa_pa": 0.0,
            "material_kernel_formulation": "standard",
        }
    else:
        expected_model = {
            **common_model,
            "local_pressure_bulk_modulus_pa": 0.0,
            "material_kappa_pa": 1000000.0,
            "material_kernel_formulation": "fbar_mechanics",
        }
    if "material_model_id" in spec:
        expected_model["material_model_id"] = spec["material_model_id"]
    if configuration.get("model_parameters") != expected_model:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered model parameters"
        )

    sampling = configuration.get("sampling_points")
    _require_exact_keys(
        sampling,
        {"p0", "p1"},
        f"{report_name}:result.configuration.sampling_points",
        artifact,
    )
    for point in ("p0", "p1"):
        record = sampling[point]
        _require_exact_keys(
            record,
            {
                "element",
                "natural_coordinates",
                "reconstruction_error_m",
                "weights",
            },
            f"{report_name}:sampling_points.{point}",
            artifact,
        )
        natural = record["natural_coordinates"]
        weights = record["weights"]
        if (
            not _integer_number(record["element"])
            or record["element"] >= spec["mesh"]["elements"]
            or not _finite_vector(natural, 3)
            or any(abs(float(value)) > 1.0 + 1.0e-8 for value in natural)
            or not _finite_vector(weights, 8)
            or any(float(value) < -1.0e-10 for value in weights)
            or not math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-10)
            or not _finite_number(record["reconstruction_error_m"])
            or record["reconstruction_error_m"] < 0.0
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has invalid {point} Hex8 sampling metadata"
            )


def _validate_current_report_solver(
    result: dict, spec: dict, report_name: str, artifact: Path
) -> None:
    solver = spec["solver"]
    configuration = result["solver_configuration"]
    records_domain_rejections = bool(
        spec.get("function_domain_diagnostics", False)
    )
    if records_domain_rejections and solver != "petsc-snes":
        raise SystemExit(
            f"release guard marks non-PETSc report {report_name} with "
            "function-domain diagnostics"
        )
    if solver == "core-newton":
        expected_configuration = dict(CORE_NEWTON_CONFIGURATION)
    else:
        expected_configuration = dict(PETSC_SNES_CONFIGURATION)
        if records_domain_rejections:
            expected_configuration["function_domain_rejection_api"] = (
                PETSC_FUNCTION_DOMAIN_REJECTION_API
            )
    if spec.get("method_metadata", False):
        expected_configuration.update(
            {
                "compiled_material_residual_only_available": True,
                "element_evaluation_mode": "joint",
            }
        )
    if configuration != expected_configuration:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered {solver} configuration"
        )
    diagnostics = result["nonlinear_step_diagnostics"]
    if not isinstance(diagnostics, list) or len(diagnostics) != spec["steps"]:
        raise SystemExit(
            f"{artifact.name}:{report_name} does not retain one nonlinear "
            "diagnostic per completed step"
        )
    core_keys = {"dt", "nonlinear_iterations", "time"}
    petsc_keys = core_keys | {
        "assembly_seconds",
        "final_residual_norm",
        "initial_residual_norm",
        "ksp_converged_reason",
        "linear_iterations",
        "petsc_function_norm",
        "residual_acceptance_threshold",
        "residual_history",
        "snes_converged_reason",
        "solve_seconds",
    }
    if records_domain_rejections:
        petsc_keys |= {
            "function_domain_rejections",
            "last_function_domain_error",
        }
    for index, record in enumerate(diagnostics, start=1):
        _require_exact_keys(
            record,
            core_keys if solver == "core-newton" else petsc_keys,
            f"{report_name}:nonlinear diagnostic {index}",
            artifact,
        )
        if (
            not _finite_number(record["time"])
            or not math.isclose(
                record["time"], index * spec["dt_s"],
                rel_tol=1.0e-12, abs_tol=1.0e-14,
            )
            or record["dt"] != spec["dt_s"]
            or not _integer_number(record["nonlinear_iterations"])
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has invalid diagnostic at step {index}"
            )
        if solver == "core-newton":
            continue

        if records_domain_rejections:
            rejection_count = record["function_domain_rejections"]
            last_error = record["last_function_domain_error"]
            if (
                not _integer_number(rejection_count)
                or (rejection_count == 0 and last_error is not None)
                or (
                    rejection_count > 0
                    and (not isinstance(last_error, str) or not last_error.strip())
                )
            ):
                raise SystemExit(
                    f"{artifact.name}:{report_name} has invalid function-domain "
                    f"rejection evidence at step {index}"
                )

        nonnegative = (
            "initial_residual_norm",
            "final_residual_norm",
            "petsc_function_norm",
            "assembly_seconds",
            "solve_seconds",
        )
        if any(
            not _finite_number(record[key]) or record[key] < 0.0
            for key in nonnegative
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has non-finite PETSc evidence "
                f"at step {index}"
            )
        threshold = max(
            expected_configuration["atol"],
            expected_configuration["rtol"] * record["initial_residual_norm"],
        )
        if (
            not math.isclose(
                record["residual_acceptance_threshold"], threshold,
                rel_tol=1.0e-15, abs_tol=0.0,
            )
            or record["final_residual_norm"] > threshold
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} fails the independently "
                f"recomputed PETSc residual rule at step {index}"
            )
        if (
            not _integer_number(record["snes_converged_reason"], minimum=1)
            or not _integer_number(record["ksp_converged_reason"])
            or not _integer_number(record["linear_iterations"])
            or (
                record["linear_iterations"] > 0
                and record["ksp_converged_reason"] == 0
            )
            or not isinstance(record["residual_history"], list)
            or not record["residual_history"]
            or any(
                not _finite_number(value) or value < 0.0
                for value in record["residual_history"]
            )
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has invalid PETSc convergence "
                f"diagnostics at step {index}"
            )


def _validate_current_report_histories(
    result: dict, spec: dict, report_name: str, artifact: Path
) -> None:
    histories = result["retained_histories"]
    expected_keys = {
        "active_tension_pa",
        "cavity_pressure_pa",
        "times_s",
        "u0_m",
        "u1_m",
    }
    _require_exact_keys(
        histories, expected_keys, f"{report_name}:result.retained_histories", artifact
    )
    sample_count = spec["steps"] + 1
    times = histories["times_s"]
    if (
        not isinstance(times, list)
        or len(times) != sample_count
        or any(
            not _finite_number(value)
            or not math.isclose(
                value, index * spec["dt_s"], rel_tol=1.0e-12, abs_tol=1.0e-14
            )
            for index, value in enumerate(times)
        )
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has an invalid retained time grid"
        )
    for field in ("active_tension_pa", "cavity_pressure_pa"):
        values = histories[field]
        if (
            not isinstance(values, list)
            or len(values) != sample_count
            or not all(_finite_number(value) for value in values)
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has invalid {field} history"
            )
    for field in ("u0_m", "u1_m"):
        values = histories[field]
        if (
            not isinstance(values, list)
            or len(values) != sample_count
            or not all(_finite_vector(value, 3) for value in values)
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has invalid {field} history"
            )

    peak_field = "active_tension_pa" if spec["case"] == "A" else "cavity_pressure_pa"
    inactive_field = (
        "cavity_pressure_pa" if spec["case"] == "A" else "active_tension_pa"
    )
    if any(value != 0.0 for value in histories[inactive_field]):
        raise SystemExit(
            f"{artifact.name}:{report_name} has loading inconsistent with Case "
            f"{spec['case']}"
        )
    peak_index = max(
        range(sample_count), key=lambda index: histories[peak_field][index]
    )
    if peak_index != spec["peak_index"]:
        raise SystemExit(
            f"{artifact.name}:{report_name} peak index is inconsistent with its history"
        )
    peak = result["peak"]
    expected_peak_keys = {
        "active_tension_pa",
        "available",
        "cavity_pressure_pa",
        "index",
        "time_s",
        "u0_m",
        "u1_m",
    }
    _require_exact_keys(peak, expected_peak_keys, f"{report_name}:result.peak", artifact)
    expected_peak = {
        "available": True,
        "index": peak_index,
        "time_s": times[peak_index],
        "active_tension_pa": histories["active_tension_pa"][peak_index],
        "cavity_pressure_pa": histories["cavity_pressure_pa"][peak_index],
        "u0_m": histories["u0_m"][peak_index],
        "u1_m": histories["u1_m"][peak_index],
    }
    if peak != expected_peak:
        raise SystemExit(
            f"{artifact.name}:{report_name} peak summary disagrees with its history"
        )


def _validate_current_report_summaries(
    result: dict, spec: dict, report_name: str, artifact: Path
) -> None:
    det_f = result["det_f_gauss_peak_summary"]
    summary_keys = {"available", "count", "maximum", "mean", "minimum", "shape"}
    _require_exact_keys(det_f, summary_keys, f"{report_name}:det_f summary", artifact)
    elements = spec["mesh"]["elements"]
    if (
        det_f["available"] is not True
        or det_f["shape"] != [elements, 8]
        or det_f["count"] != elements * 8
        or not all(_finite_number(det_f[key]) for key in ("minimum", "mean", "maximum"))
        or not 0.0 < det_f["minimum"] <= det_f["mean"] <= det_f["maximum"]
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} lacks positive finite 8-GP det(F) evidence"
        )

    pressure = result["element_pressure_peak_pa_summary"]
    local = spec["formulation"] == "hex8_local_pressure_p0_condensed_logj"
    if not local:
        if pressure != {"available": False}:
            raise SystemExit(
                f"{artifact.name}:{report_name} unexpectedly reports local pressure"
            )
        return
    _require_exact_keys(
        pressure, summary_keys, f"{report_name}:element pressure summary", artifact
    )
    if (
        pressure["available"] is not True
        or pressure["shape"] != [elements]
        or pressure["count"] != elements
        or not all(_finite_number(pressure[key]) for key in ("minimum", "mean", "maximum"))
        or not pressure["minimum"] <= pressure["mean"] <= pressure["maximum"]
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has invalid local-pressure evidence"
        )


def _validate_peak_ring_rotation(
    result: dict, spec: dict, report_name: str, artifact: Path
) -> None:
    profile = result["peak_circumferential_ring_rotation"]
    if spec["case"] != "B":
        if profile != {
            "available": False,
            "reason": "Reported only for Case B results.",
        }:
            raise SystemExit(
                f"{artifact.name}:{report_name} has unexpected Case A ring rotation"
            )
        return

    expected_keys = {
        "angle_definition",
        "apex_longitudinal_ring_index",
        "available",
        "base_longitudinal_ring_index",
        "centering",
        "circumferential_nodes_per_ring",
        "coordinate_system",
        "layers",
        "method",
        "peak_time_s",
        "relative_angle_definition",
        "sign_convention",
    }
    _require_exact_keys(
        profile, expected_keys, f"{report_name}:peak ring rotation", artifact
    )
    expected_labels = {
        "available": True,
        "method": "least-squares rotation of centered circumferential rings",
        "coordinate_system": (
            "global Cartesian; long axis +x; circumferential plane yz"
        ),
        "centering": (
            "subtract separate reference and deformed yz centroids per ring"
        ),
        "angle_definition": (
            "atan2(sum(Y_y*Yd_z - Y_z*Yd_y), "
            "sum(Y_y*Yd_y + Y_z*Yd_z))"
        ),
        "sign_convention": (
            "positive is right-handed about global +x (+y toward +z)"
        ),
        "relative_angle_definition": (
            "apex minus base, wrapped to [-180, 180) degrees"
        ),
        "peak_time_s": result["peak"]["time_s"],
        "base_longitudinal_ring_index": 0,
        "apex_longitudinal_ring_index": spec["mesh"]["n_mu"],
        "circumferential_nodes_per_ring": spec["mesh"]["n_theta"],
    }
    altered = sorted(
        key for key, expected in expected_labels.items()
        if profile.get(key) != expected
    )
    if altered:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered ring-rotation definition: "
            f"{altered}"
        )

    layers = profile["layers"]
    expected_relative = spec["ring_rotation_relative_degrees"]
    n_t = spec["mesh"]["n_t"]
    if not isinstance(layers, list) or len(layers) != n_t + 1:
        raise SystemExit(
            f"{artifact.name}:{report_name} has incomplete transmural ring rotations"
        )
    layer_keys = {
        "apex_rotation_degrees",
        "base_rotation_degrees",
        "relative_apex_minus_base_degrees",
        "transmural_coordinate",
        "transmural_node_layer_index",
    }
    for index, (layer, expected_angle) in enumerate(
        zip(layers, expected_relative)
    ):
        _require_exact_keys(
            layer,
            layer_keys,
            f"{report_name}:ring rotation layer {index}",
            artifact,
        )
        base = layer["base_rotation_degrees"]
        apex = layer["apex_rotation_degrees"]
        relative = layer["relative_apex_minus_base_degrees"]
        recomputed = (float(apex) - float(base) + 180.0) % 360.0 - 180.0
        if (
            layer["transmural_node_layer_index"] != index
            or layer["transmural_coordinate"] != index / n_t
            or any(
                not _finite_number(value) or not -180.0 <= value < 180.0
                for value in (base, apex, relative)
            )
            or not math.isclose(
                relative, recomputed, rel_tol=1.0e-12, abs_tol=1.0e-12
            )
            or not math.isclose(
                relative, expected_angle, rel_tol=1.0e-12, abs_tol=1.0e-12
            )
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has invalid ring rotation at "
                f"transmural node layer {index}"
            )


def _validate_current_report_reference_and_comparison(
    report: dict, spec: dict, report_name: str, artifact: Path
) -> None:
    reference = report["reference"]
    _require_exact_keys(
        reference,
        {
            "canonical_grid_s",
            "case",
            "doi",
            "endpoint_offset_tolerance_s",
            "license",
            "mean_curves_m",
            "published_archive_identity",
            "selection",
            "team_files",
        },
        f"{report_name}:reference",
        artifact,
    )
    expected_labels = {
        "case": spec["reference_case"],
        "doi": "10.5281/zenodo.14260459",
        "endpoint_offset_tolerance_s": 0.001000001,
        "license": "CC-BY-4.0",
        "published_archive_identity": REFERENCE_ARCHIVE_IDENTITY,
        "selection": _expected_reference_selection(spec["reference_case"]),
        "team_files": _expected_reference_team_files(spec["reference_case"]),
    }
    altered = sorted(
        key for key, expected in expected_labels.items() if reference.get(key) != expected
    )
    if altered:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered external reference "
            f"identity: {altered}"
        )
    grid = reference["canonical_grid_s"]
    if (
        not isinstance(grid, list)
        or len(grid) != 101
        or any(
            not _finite_number(value)
            or not math.isclose(value, index * 0.01, rel_tol=0.0, abs_tol=1.0e-15)
            for index, value in enumerate(grid)
        )
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} does not use the canonical 101-point grid"
        )

    means = reference["mean_curves_m"]
    comparison = report["comparison"]
    _require_exact_keys(means, {"p0", "p1"}, f"{report_name}:mean curves", artifact)
    _require_exact_keys(
        comparison,
        {"metric", "ours_on_canonical_grid_m", "red"},
        f"{report_name}:comparison",
        artifact,
    )
    if comparison["metric"] != "relative discrepancy (benchmark paper Eq. 21)":
        raise SystemExit(f"{artifact.name}:{report_name} has altered comparison metric")
    ours = comparison["ours_on_canonical_grid_m"]
    red = comparison["red"]
    _require_exact_keys(ours, {"p0", "p1"}, f"{report_name}:ours curves", artifact)
    _require_exact_keys(red, {"p0", "p1"}, f"{report_name}:RED", artifact)
    team_names = set(REFERENCE_TEAM_SHA256[spec["reference_case"]])
    for point in ("p0", "p1"):
        mean_curve = _validate_curve(
            means[point], f"{report_name}:reference mean {point}", artifact
        )
        ours_curve = _validate_curve(
            ours[point], f"{report_name}:CoupFE {point}", artifact
        )
        record = red[point]
        _require_exact_keys(record, {"ours", "teams"}, f"{report_name}:RED {point}", artifact)
        recomputed = _relative_discrepancy(ours_curve, mean_curve)
        if (
            not _finite_number(record["ours"])
            or record["ours"] < 0.0
            or not math.isclose(record["ours"], recomputed, rel_tol=1.0e-12, abs_tol=1.0e-14)
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has an unreproducible {point} RED"
            )
        teams = record["teams"]
        if (
            not isinstance(teams, dict)
            or set(teams) != team_names
            or any(not _finite_number(value) or value < 0.0 for value in teams.values())
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has incomplete reference-team RED evidence"
            )


def _validate_archived_truncated_polar_report(
    payloads: dict[str, bytes], artifact: Path, spec: dict
) -> None:
    report_name = _truncated_polar_archive_path(spec["report"])
    log_name = _truncated_polar_archive_path(spec["log"])
    missing = sorted({report_name, log_name} - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing archived truncated-polar evidence: "
            f"{missing}"
        )
    if hashlib.sha256(payloads[report_name]).hexdigest() != spec["report_sha256"]:
        raise SystemExit(
            f"{artifact.name}:{report_name} does not match its reviewed SHA-256"
        )
    if hashlib.sha256(payloads[log_name]).hexdigest() != spec["log_sha256"]:
        raise SystemExit(
            f"{artifact.name}:{log_name} does not match its reviewed SHA-256"
        )

    report = _load_strict_json(payloads[report_name], report_name, artifact)
    _validate_no_absolute_paths(report, report_name, artifact)
    expected_report_keys = {
        "bounded_claim",
        "comparison",
        "reference",
        "result",
        "schema",
    }
    if "predecessor_report_sha256" in spec:
        expected_report_keys.add("correction")
    _require_exact_keys(report, expected_report_keys, report_name, artifact)
    expected_bounded_claim = (
        CURRENT_REPORT_DISCRETIZATION_BOUNDED_CLAIM
        if spec.get("current_claim_boundary", False)
        else CURRENT_REPORT_BOUNDED_CLAIM
    )
    if (
        report["schema"] != CURRENT_REPORT_SCHEMA
        or report["bounded_claim"] != expected_bounded_claim
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered schema or claim boundary"
        )
    expected_correction = None
    if "predecessor_report_sha256" in spec:
        expected_correction = {
            "predecessor_repository_revision": CORRECTION_PREDECESSOR_REVISION,
            "reason": CORRECTION_REASON,
            "supersedes_report_sha256": spec["predecessor_report_sha256"],
        }
    if report.get("correction") != expected_correction:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered correction lineage"
        )
    result = report["result"]
    expected_result_keys = {
        "case",
        "configuration",
        "det_f_gauss_peak_summary",
        "element_pressure_peak_pa_summary",
        "filename",
        "nonlinear_step_diagnostics",
        "normalized_run_log",
        "peak",
        "peak_circumferential_ring_rotation",
        "reference_case",
        "result_schema",
        "retained_histories",
        "runtime_versions",
        "sha256",
        "size_bytes",
        "solver_configuration",
        "source_identity",
    }
    if spec.get("method_metadata", False):
        expected_result_keys.add("pre_solve_audit")
    _require_exact_keys(
        result,
        expected_result_keys,
        f"{report_name}:result",
        artifact,
    )
    if spec.get("method_metadata", False) and result["pre_solve_audit"] is not None:
        raise SystemExit(
            f"{artifact.name}:{report_name} unexpectedly assigns a closed-mesh "
            "pre-solve audit to the polar-ring record"
        )
    expected_result_identity = {
        "case": spec["case"],
        "filename": spec["result"],
        "reference_case": spec["reference_case"],
        "result_schema": CURRENT_RESULT_SCHEMA,
        "runtime_versions": CURRENT_RUNTIME_VERSIONS,
        "sha256": spec["result_sha256"],
        "size_bytes": spec["result_size_bytes"],
        "source_identity": {
            "app": {
                "revision": spec["app_ref"],
                "source_kind": "git-checkout",
                "tree_state": "clean",
            },
            "core": {
                "revision": spec.get(
                    "core_ref", HISTORICAL_RETAINED_CORE_REF
                ),
                "source_kind": "git-checkout",
                "source_url": PUBLIC_CORE_URL,
                "tree_state": "clean",
            },
        },
    }
    altered = sorted(
        key for key, expected in expected_result_identity.items()
        if result.get(key) != expected
    )
    if altered:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered result identity: {altered}"
        )
    _validate_current_report_configuration(result["configuration"], spec, report_name, artifact)
    _validate_current_report_solver(result, spec, report_name, artifact)
    _validate_current_report_histories(result, spec, report_name, artifact)
    _validate_current_report_summaries(result, spec, report_name, artifact)
    _validate_peak_ring_rotation(result, spec, report_name, artifact)
    _validate_current_report_reference_and_comparison(report, spec, report_name, artifact)

    log_metadata = result["normalized_run_log"]
    expected_log_metadata = {
        "filename": spec["log"],
        "normalization": "UTF-8; CRLF/CR converted to LF; final LF added when nonempty",
        "sha256": spec["log_sha256"],
        "size_bytes": len(payloads[log_name]),
    }
    if log_metadata != expected_log_metadata:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered normalized-log identity"
        )
    try:
        log = payloads[log_name].decode("utf-8")
    except UnicodeDecodeError as exc:
        raise SystemExit(f"{artifact.name}:{log_name} is not UTF-8") from exc
    if "\x00" in log or "\r" in log or (log and not log.endswith("\n")):
        raise SystemExit(f"{artifact.name}:{log_name} is not normalized text")
    markers = (
        f"mesh: {spec['mesh']['nodes']} nodes, {spec['mesh']['elements']} hexes, "
        f"ndof={spec['mesh']['degrees_of_freedom']}",
        f"CASE {spec['case']}:",
        "nonlinear solver:",
        "Hex8 output sampling:",
        f"step {spec['steps']:4d} t=1.000s",
        f"saved -> {spec['result']}",
    )
    missing_markers = [marker for marker in markers if marker not in log]
    if missing_markers:
        raise SystemExit(
            f"{artifact.name}:{log_name} is missing completion markers: "
            f"{missing_markers}"
        )


def _validate_truncated_polar_archive(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Bind archived non-benchmark reports to exact source, logs, and semantics."""
    for spec in TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS:
        app_ref = spec.get("app_ref")
        if app_ref not in CURRENT_RESULT_SOURCE_CHECKPOINTS:
            raise SystemExit(
                f"release guard report {spec.get('report')!r} cites unreviewed "
                f"source checkpoint {app_ref!r}"
            )
    _validate_current_release_source_hashes(payloads, artifact)
    _validate_current_reporting_source_hashes(payloads, artifact)
    for spec in TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS:
        _validate_archived_truncated_polar_report(payloads, artifact, spec)


def _sensitive_fragments() -> tuple[str, ...]:
    # Keep literal release secrets and machine-specific paths out of this file too.
    return (
        "/" + "home/",
        "/" + "media/",
        "/" + "mnt/",
        "/" + "root/",
        "/" + "tmp/",
        "/" + "Users/",
        "/" + "var/folders/",
        "C:" + "\\Users\\",
        "git" + "@jetstream",
        "gh" + "p_",
        "github" + "_pat_",
        "gl" + "pat-",
        "sk" + "-proj-",
        "xo" + "xb-",
        "xo" + "xp-",
        "BEGIN " + "PRIVATE KEY",
        "BEGIN RSA" + " PRIVATE KEY",
        "BEGIN OPENSSH" + " PRIVATE KEY",
        "AK" + "IA",
    )


def _credential_assignment_patterns() -> tuple[re.Pattern[str], ...]:
    words = (
        "pass" + "word",
        "api" + r"[_-]?" + "key",
        "access" + r"[_-]?" + "token",
        "client" + r"[_-]?" + "secret",
    )
    return tuple(
        re.compile(
            rf"(?i)\b{word}\b\s*[:=]\s*[\"']?[A-Za-z0-9_./+=-]{{8,}}"
        )
        for word in words
    )


def _validate_text(name: str, payload: bytes, artifact: Path) -> None:
    path = PurePosixPath(name)
    if path.suffix.casefold() not in TEXT_SUFFIXES:
        raise SystemExit(
            f"{artifact.name}:{name} has an unapproved file type "
            f"{path.suffix or '<none>'}"
        )
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise SystemExit(f"{artifact.name}:{name} is not valid UTF-8 text") from exc

    hits = [fragment for fragment in _sensitive_fragments() if fragment in text]
    if re.search(r"[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}\b", text):
        hits.append("personal-email-address")
    if re.search(r"(?<![.\w])" + "git" + r"@[A-Za-z0-9.-]+:", text):
        hits.append("ssh-git-url")
    if re.search(r"[A-Za-z]:" + r"\\Users\\", text):
        hits.append("windows-user-path")
    if any(pattern.search(text) for pattern in _credential_assignment_patterns()):
        hits.append("credential-assignment")
    folded = text.casefold()
    if any(term in folded for term in _retired_terms()):
        hits.append("retired-private-or-mixed-path")
    if hits:
        raise SystemExit(
            f"{artifact.name}:{name} contains private, credential, or retired-path "
            f"material: {sorted(set(hits))}"
        )


def _xml_local_name(name: str) -> str:
    return name.rsplit("}", 1)[-1]


def _validate_retained_svg(
    payload: bytes, name: str, spec: dict, artifact: Path
) -> None:
    """Validate one inert, accessible, source-labelled SVG figure."""
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise SystemExit(f"{artifact.name}:{name} is not UTF-8 SVG text") from exc

    without_declaration = re.sub(
        r"^\s*<\?xml[^?]*\?>", "", text, count=1, flags=re.IGNORECASE
    )
    forbidden_markup = (
        "<!DOCTYPE",
        "<!ENTITY",
        "<?xml-stylesheet",
        "<script",
        "javascript:",
        "file:",
        "data:",
    )
    folded = without_declaration.casefold()
    hits = [token for token in forbidden_markup if token.casefold() in folded]
    if "<?" in without_declaration:
        hits.append("processing-instruction")
    if hits:
        raise SystemExit(
            f"{artifact.name}:{name} contains active or external SVG markup: "
            f"{sorted(set(hits))}"
        )

    try:
        root = ElementTree.fromstring(payload)
    except ElementTree.ParseError as exc:
        raise SystemExit(f"{artifact.name}:{name} is not well-formed XML") from exc
    svg_namespace = "http://www.w3.org/2000/svg"
    if root.tag != f"{{{svg_namespace}}}svg":
        raise SystemExit(f"{artifact.name}:{name} does not have an SVG root")
    if root.get("role") != "img" or root.get("aria-labelledby") != (
        "figure-title figure-description"
    ):
        raise SystemExit(
            f"{artifact.name}:{name} lacks the reviewed accessible image role"
        )

    title_elements = root.findall(f"{{{svg_namespace}}}title")
    description_elements = root.findall(f"{{{svg_namespace}}}desc")
    if (
        len(title_elements) != 1
        or title_elements[0].get("id") != "figure-title"
        or "".join(title_elements[0].itertext()) != spec["title"]
    ):
        raise SystemExit(
            f"{artifact.name}:{name} lacks its exact accessible SVG title"
        )
    if (
        len(description_elements) != 1
        or description_elements[0].get("id") != "figure-description"
        or "".join(description_elements[0].itertext()) != spec["description"]
    ):
        raise SystemExit(
            f"{artifact.name}:{name} lacks its exact accessible SVG description"
        )

    elements = list(root.iter())
    forbidden_elements = {"foreignObject", "iframe", "image", "script"}
    present_forbidden = sorted(
        {
            _xml_local_name(element.tag)
            for element in elements
            if _xml_local_name(element.tag) in forbidden_elements
        }
    )
    foreign_elements = sorted(
        {
            element.tag
            for element in elements
            if not element.tag.startswith(f"{{{svg_namespace}}}")
        }
    )
    if present_forbidden or foreign_elements:
        raise SystemExit(
            f"{artifact.name}:{name} contains forbidden or foreign SVG elements: "
            f"{present_forbidden + foreign_elements}"
        )

    ids = [element.get("id") for element in elements if element.get("id")]
    if len(ids) != len(set(ids)):
        raise SystemExit(f"{artifact.name}:{name} contains duplicate SVG IDs")
    available_ids = set(ids)
    referenced_ids = set()
    for element in elements:
        for attribute, value in element.attrib.items():
            local_attribute = _xml_local_name(attribute)
            folded_value = value.strip().casefold()
            if local_attribute.casefold().startswith("on"):
                raise SystemExit(
                    f"{artifact.name}:{name} contains an SVG event handler"
                )
            if local_attribute in {"href", "src"}:
                if not value.startswith("#") or len(value) == 1:
                    raise SystemExit(
                        f"{artifact.name}:{name} contains an external SVG link"
                    )
                referenced_ids.add(value[1:])
            elif local_attribute in {"about", "base", "resource"} and value:
                raise SystemExit(
                    f"{artifact.name}:{name} contains external resource metadata"
                )
            if (
                "javascript:" in folded_value
                or "file:" in folded_value
                or "data:" in folded_value
                or re.search(r"https?://", folded_value)
                or Path(value).is_absolute()
                or PureWindowsPath(value).is_absolute()
            ):
                raise SystemExit(
                    f"{artifact.name}:{name} contains an external or absolute SVG value"
                )
            css_urls = re.findall(r"url\(([^)]*)\)", value, flags=re.IGNORECASE)
            if folded_value.count("url(") != len(css_urls) or any(
                re.fullmatch(r"\s*#[^)\s]+\s*", target) is None
                for target in css_urls
            ):
                raise SystemExit(
                    f"{artifact.name}:{name} contains a non-local CSS URL"
                )
            referenced_ids.update(target.strip()[1:] for target in css_urls)

    missing_ids = sorted(referenced_ids - available_ids)
    if missing_ids:
        raise SystemExit(
            f"{artifact.name}:{name} references missing SVG IDs: {missing_ids}"
        )
    visible_text = " ".join(" ".join(root.itertext()).split())
    missing_markers = [
        marker for marker in spec["visible_markers"] if marker not in visible_text
    ]
    if missing_markers:
        raise SystemExit(
            f"{artifact.name}:{name} lacks reviewed provenance labels: "
            f"{missing_markers}"
        )


def _canonical_semantic_sha256(value) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_step0a_retained_report_semantics(
    report: dict, report_name: str, artifact: Path
) -> None:
    """Validate the compact fine-mesh Case A comparison independently."""
    _validate_no_absolute_paths(report, report_name, artifact)
    _require_exact_keys(
        report,
        {
            "benchmark_identity",
            "bounded_claim",
            "comparison",
            "reference",
            "result",
            "schema",
        },
        report_name,
        artifact,
    )
    if (
        report["schema"] != STEP0A_RETAINED_REPORT_SCHEMA
        or report["bounded_claim"] != CURRENT_REPORT_DISCRETIZATION_BOUNDED_CLAIM
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered Step 0A schema or "
            "claim boundary"
        )

    benchmark_identity = report["benchmark_identity"]
    _require_exact_keys(
        benchmark_identity,
        {
            "benchmark",
            "case",
            "claim_boundary",
            "explicit_archive_identity_fields",
            "inference_basis",
            "load_history_audit",
            "missing_archive_fields",
            "status",
        },
        f"{report_name}:benchmark_identity",
        artifact,
    )
    if (
        benchmark_identity.get("benchmark") != "Benchmark 1"
        or benchmark_identity.get("case") != "step_0A"
        or benchmark_identity.get("status") != "legacy-inferred"
        or benchmark_identity.get("explicit_archive_identity_fields") is not False
        or _canonical_semantic_sha256(benchmark_identity)
        != STEP0A_RETAINED_SEMANTIC_SHA256["benchmark_identity"]
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} alters the reviewed legacy-inferred "
            "Step 0A identity"
        )

    result = report["result"]
    _require_exact_keys(
        result,
        {
            "campaign_records",
            "case",
            "completion",
            "configuration",
            "filename",
            "reference_case",
            "result_schema",
            "sha256",
            "size_bytes",
            "source_identity",
        },
        f"{report_name}:result",
        artifact,
    )
    expected_result_identity = {
        **STEP0A_RETAINED_RESULT_IDENTITY,
        "case": "A",
        "reference_case": "step_0A",
        "result_schema": CURRENT_RESULT_SCHEMA,
        "source_identity": {
            "app": {
                "revision": STEP0A_RETAINED_APP_REF,
                "source_kind": "git-checkout",
                "tree_state": "clean",
            },
            "core": {
                "revision": APPROVED_PUBLIC_CORE_REF,
                "source_kind": "git-checkout",
                "source_url": PUBLIC_CORE_URL,
                "tree_state": "clean",
            },
        },
        "campaign_records": {
            "manifest": STEP0A_RETAINED_MANIFEST_IDENTITY,
            "stdout": STEP0A_RETAINED_STDOUT_IDENTITY,
        },
    }
    altered_identity = sorted(
        key
        for key, expected in expected_result_identity.items()
        if result.get(key) != expected
    )
    if altered_identity:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered fine Case A result "
            f"identity: {altered_identity}"
        )

    configuration = result["configuration"]
    if (
        _canonical_semantic_sha256(configuration)
        != STEP0A_RETAINED_SEMANTIC_SHA256["configuration"]
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered retained configuration"
        )
    mesh = configuration.get("mesh")
    expected_mesh = {
        "core_half_width": 0.36,
        "degrees_of_freedom": 89655,
        "elements": 23616,
        "n_core": 36,
        "n_mu": 0,
        "n_radial": 32,
        "n_side": 0,
        "n_t": 4,
        "n_theta": 0,
        "nodes": 29885,
        "topology": "closed_multiblock_disk",
    }
    expected_method = {
        "formulation": "hex8_local_pressure_p0_condensed_logj",
        "dt_s": 0.001,
        "t_end_s": 1.0,
        "load_horizon_s": 1.0,
        "integrator": "generalized-alpha",
        "generalized_alpha": {
            "alpha_f": 0.4,
            "alpha_m": 0.2,
            "beta": 0.36,
            "gamma": 0.7,
            "stage_contract": "simula-source-matched-v1",
        },
        "mass_representation": "consistent_q1_hex8",
        "material_eta_pa_s": 100.0,
        "viscous_term_active": True,
        "fiber_sampling": "gp_direct_rule",
        "point_sampling": "hex8_reference_isoparametric",
    }
    altered_method = sorted(
        key
        for key, expected in expected_method.items()
        if configuration.get(key) != expected
    )
    if mesh != expected_mesh or altered_method:
        raise SystemExit(
            f"{artifact.name}:{report_name} does not retain the reviewed closed "
            f"4x36x32 local-pressure/generalized-alpha method"
        )
    model = configuration.get("model_parameters")
    if (
        not isinstance(model, dict)
        or model.get("local_pressure_bulk_modulus_pa") != 1000000.0
        or model.get("local_pressure_volume_law")
        != "linear-reference-volume-mean-log-j-v1"
        or model.get("material_kappa_pa") != 0.0
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered local-pressure volume model"
        )
    mpi = configuration.get("mpi")
    if (
        not isinstance(mpi, dict)
        or mpi.get("ranks") != 8
        or mpi.get("threads_per_rank") != 1
        or mpi.get("linear_solver_profile")
        != "fgmres-gamg-rigid-rebuild"
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered parallel solver contract"
        )

    completion = result["completion"]
    if (
        _canonical_semantic_sha256(completion)
        != STEP0A_RETAINED_SEMANTIC_SHA256["completion"]
        or completion.get("converged") is not True
        or completion.get("completed_steps") != 1000
        or completion.get("expected_steps") != 1000
        or completion.get("function_domain_rejections") != 0
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered completion evidence"
        )
    det_f = completion.get("det_f_gauss_peak_summary")
    if (
        not isinstance(det_f, dict)
        or det_f.get("shape") != [23616, 8]
        or det_f.get("count") != 23616 * 8
        or not all(
            _finite_number(det_f.get(key))
            for key in ("minimum", "mean", "maximum")
        )
        or not 0.0 < det_f["minimum"] <= det_f["mean"] <= det_f["maximum"]
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} lacks positive finite fine-mesh "
            "det(F) evidence"
        )

    reference = report["reference"]
    _require_exact_keys(
        reference,
        {
            "canonical_grid_s",
            "case",
            "doi",
            "license",
            "mean_curves_m",
            "published_archive_identity",
            "selection",
            "team_files",
        },
        f"{report_name}:reference",
        artifact,
    )
    expected_reference_identity = {
        "case": "step_0A",
        "doi": "10.5281/zenodo.14260459",
        "license": "CC-BY-4.0",
        "published_archive_identity": REFERENCE_ARCHIVE_IDENTITY,
        "selection": _expected_reference_selection("step_0A"),
        "team_files": _expected_reference_team_files("step_0A"),
    }
    altered_reference = sorted(
        key
        for key, expected in expected_reference_identity.items()
        if reference.get(key) != expected
    )
    if altered_reference:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered external reference "
            f"identity: {altered_reference}"
        )
    grid = reference["canonical_grid_s"]
    if (
        not isinstance(grid, list)
        or len(grid) != 101
        or any(
            not _finite_number(value)
            or not math.isclose(value, index * 0.01, rel_tol=0.0, abs_tol=1.0e-15)
            for index, value in enumerate(grid)
        )
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} does not retain the canonical "
            "101-point grid"
        )

    means = reference["mean_curves_m"]
    comparison = report["comparison"]
    _require_exact_keys(means, {"p0", "p1"}, f"{report_name}:means", artifact)
    _require_exact_keys(
        comparison,
        {"metric", "ours_on_canonical_grid_m", "red"},
        f"{report_name}:comparison",
        artifact,
    )
    if comparison.get("metric") != "relative discrepancy (benchmark paper Eq. 21)":
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered comparison metric"
        )
    ours = comparison["ours_on_canonical_grid_m"]
    red = comparison["red"]
    _require_exact_keys(ours, {"p0", "p1"}, f"{report_name}:ours", artifact)
    _require_exact_keys(red, {"p0", "p1"}, f"{report_name}:RED", artifact)
    for point in ("p0", "p1"):
        mean_curve = _validate_curve(
            means[point], f"{report_name}:reference mean {point}", artifact
        )
        ours_curve = _validate_curve(
            ours[point], f"{report_name}:CoupFE {point}", artifact
        )
        if (
            _canonical_semantic_sha256(mean_curve)
            != STEP0A_RETAINED_SEMANTIC_SHA256[f"mean_{point}"]
            or _canonical_semantic_sha256(ours_curve)
            != STEP0A_RETAINED_SEMANTIC_SHA256[f"ours_{point}"]
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has altered retained {point} curves"
            )
        record = red[point]
        _require_exact_keys(
            record, {"ours", "teams"}, f"{report_name}:RED {point}", artifact
        )
        recomputed = _relative_discrepancy(ours_curve, mean_curve)
        if (
            not _finite_number(record["ours"])
            or not math.isclose(
                record["ours"], recomputed, rel_tol=1.0e-12, abs_tol=1.0e-14
            )
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has an unreproducible {point} RED"
            )
        if record["teams"] != STEP0A_RETAINED_TEAM_RED[point]:
            raise SystemExit(
                f"{artifact.name}:{report_name} has altered reference-team "
                f"{point} RED evidence"
            )


def _validate_step0a_retained_comparison(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    if STEP0A_RETAINED_REPORT not in payloads:
        raise SystemExit(
            f"{artifact.name} is missing retained fine Case A comparison report"
        )
    payload = payloads[STEP0A_RETAINED_REPORT]
    if (
        len(payload) != STEP0A_RETAINED_REPORT_SIZE_BYTES
        or hashlib.sha256(payload).hexdigest() != STEP0A_RETAINED_REPORT_SHA256
    ):
        raise SystemExit(
            f"{artifact.name}:{STEP0A_RETAINED_REPORT} differs from the reviewed "
            "compact report"
        )
    report = _load_strict_json(payload, STEP0A_RETAINED_REPORT, artifact)
    if not isinstance(report, dict):
        raise SystemExit(
            f"{artifact.name}:{STEP0A_RETAINED_REPORT} must contain an object"
        )
    _validate_step0a_retained_report_semantics(
        report, STEP0A_RETAINED_REPORT, artifact
    )


def _validate_step0b_prefix_runtime_sources(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Validate the prefix's immutable historical runtime-source identity.

    The manifest describes the source bytes that produced the retained prefix,
    not the current checkout.  Current executable bytes are compared
    independently with ``CURRENT_RELEASE_SOURCE_SHA256``; comparing them with
    this historical map would freeze the package or relabel old evidence.
    """
    manifest = STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST
    if not isinstance(manifest, dict):
        raise SystemExit(
            "release guard has malformed Step 0B runtime-source manifest: "
            f"{type(manifest).__name__}"
        )
    malformed = sorted(
        name
        for name, digest in manifest.items()
        if not isinstance(name, str)
        or not isinstance(digest, str)
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    )
    if len(manifest) != 20 or malformed:
        raise SystemExit(
            "release guard has malformed Step 0B runtime-source manifest: "
            f"count={len(manifest)}, entries={malformed}"
        )
    if (
        not isinstance(STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_SHA256, str)
        or re.fullmatch(
            r"[0-9a-f]{64}", STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_SHA256
        )
        is None
    ):
        raise SystemExit(
            "release guard has malformed Step 0B runtime-source aggregate"
        )
    encoded = json.dumps(
        manifest,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    aggregate = hashlib.sha256(encoded).hexdigest()
    if aggregate != STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_SHA256:
        raise SystemExit(
            "release guard Step 0B runtime-source manifest digest differs "
            "from the reviewed aggregate"
        )
    unguarded = sorted(set(manifest) - set(CURRENT_RELEASE_SOURCE_SHA256))
    if unguarded:
        raise SystemExit(
            "release guard is missing current hashes for historical Step 0B "
            f"runtime paths: {unguarded}"
        )
    missing = sorted(set(manifest) - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing current Step 0B runtime paths: {missing}"
        )
    changed = sorted(
        name
        for name in manifest
        if hashlib.sha256(payloads[name]).hexdigest()
        != CURRENT_RELEASE_SOURCE_SHA256[name]
    )
    if changed:
        raise SystemExit(
            f"{artifact.name} current Step 0B runtime paths differ from the "
            f"current release hashes: {changed}"
        )


def _validate_step0b_prefix_diagnostic_semantics(
    report, report_name: str, artifact: Path
) -> None:
    """Validate the completed short-gate identity independently of file bytes."""
    if not isinstance(report, dict):
        raise SystemExit(
            f"{artifact.name}:{report_name} must contain a JSON object"
        )
    _validate_no_absolute_paths(report, report_name, artifact)
    _require_exact_keys(
        report,
        {
            "benchmark_identity",
            "bounded_claim",
            "comparison",
            "completion",
            "configuration",
            "decision",
            "excluded_attempts",
            "external_artifacts",
            "fenics_input_identities",
            "mesh_inputs",
            "mesh_split_diagnosis",
            "physical_frame_control",
            "reference_data",
            "robin_surface_mechanism",
            "schema",
            "source_identity",
        },
        report_name,
        artifact,
    )
    if report["schema"] != STEP0B_PREFIX_DIAGNOSTIC_SCHEMA:
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered Step 0B schema"
        )

    benchmark_identity = report["benchmark_identity"]
    _require_exact_keys(
        benchmark_identity,
        set(STEP0B_PREFIX_DIAGNOSTIC_BENCHMARK_IDENTITY),
        f"{report_name}:benchmark_identity",
        artifact,
    )
    if _canonical_semantic_sha256(benchmark_identity) != (
        _canonical_semantic_sha256(STEP0B_PREFIX_DIAGNOSTIC_BENCHMARK_IDENTITY)
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered Step 0B benchmark identity"
        )

    configuration = report["configuration"]
    _require_exact_keys(
        configuration,
        set(STEP0B_PREFIX_DIAGNOSTIC_CONFIGURATION),
        f"{report_name}:configuration",
        artifact,
    )
    if _canonical_semantic_sha256(configuration) != _canonical_semantic_sha256(
        STEP0B_PREFIX_DIAGNOSTIC_CONFIGURATION
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered Step 0B configuration"
        )

    expected_source_identity = {
        "application_revision": STEP0B_PREFIX_DIAGNOSTIC_APP_REF,
        "application_tree_state": "clean",
        "benchmark_runtime_source_sha256": (
            STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_SHA256
        ),
        "core_revision": APPROVED_PUBLIC_CORE_REF,
        "core_tree_state": "clean",
        "mesh_split_application_revision": (
            STEP0B_PREFIX_DIAGNOSTIC_MESH_SPLIT_APP_REF
        ),
    }
    source_identity = report["source_identity"]
    _require_exact_keys(
        source_identity,
        set(expected_source_identity),
        f"{report_name}:source_identity",
        artifact,
    )
    if _canonical_semantic_sha256(source_identity) != _canonical_semantic_sha256(
        expected_source_identity
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered Step 0B source identity"
        )

    completion = report["completion"]
    _require_exact_keys(
        completion,
        {"coarse_2x20x17", "wall_only_4x20x17"},
        f"{report_name}:completion",
        artifact,
    )
    for mesh in ("coarse_2x20x17", "wall_only_4x20x17"):
        record = completion[mesh]
        _require_exact_keys(
            record,
            set(STEP0B_PREFIX_DIAGNOSTIC_COMPLETION_RECORD),
            f"{report_name}:completion.{mesh}",
            artifact,
        )
        if _canonical_semantic_sha256(record) != _canonical_semantic_sha256(
            STEP0B_PREFIX_DIAGNOSTIC_COMPLETION_RECORD
        ):
            raise SystemExit(
                f"{artifact.name}:{report_name} has altered Step 0B {mesh} "
                "completion"
            )

    decision = report["decision"]
    _require_exact_keys(
        decision,
        set(STEP0B_PREFIX_DIAGNOSTIC_DECISION),
        f"{report_name}:decision",
        artifact,
    )
    if _canonical_semantic_sha256(decision) != _canonical_semantic_sha256(
        STEP0B_PREFIX_DIAGNOSTIC_DECISION
    ):
        raise SystemExit(
            f"{artifact.name}:{report_name} has altered Step 0B paused decision"
        )


def _validate_step0b_prefix_diagnostic(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Keep the reviewed short-gate decision separate from full-cycle evidence."""
    if STEP0B_PREFIX_DIAGNOSTIC_REPORT not in payloads:
        raise SystemExit(
            f"{artifact.name} is missing the Step 0B prefix diagnostic report"
        )
    payload = payloads[STEP0B_PREFIX_DIAGNOSTIC_REPORT]
    if (
        len(payload) != STEP0B_PREFIX_DIAGNOSTIC_REPORT_SIZE_BYTES
        or hashlib.sha256(payload).hexdigest()
        != STEP0B_PREFIX_DIAGNOSTIC_REPORT_SHA256
    ):
        raise SystemExit(
            f"{artifact.name}:{STEP0B_PREFIX_DIAGNOSTIC_REPORT} differs from "
            "the reviewed compact report"
        )
    _validate_step0b_prefix_runtime_sources(payloads, artifact)
    report = _load_strict_json(
        payload, STEP0B_PREFIX_DIAGNOSTIC_REPORT, artifact
    )
    _validate_step0b_prefix_diagnostic_semantics(
        report, STEP0B_PREFIX_DIAGNOSTIC_REPORT, artifact
    )


def _validate_retained_figures(
    payloads: dict[str, bytes], artifact: Path
) -> None:
    """Bind the reviewed renderers and SVGs to exact inert payloads."""
    required = {
        RETAINED_FIGURE_RENDERER,
        STEP2B_FIGURE_RENDERER,
        STEP2B_RAW_STDOUT,
        TIP_REFINE_FIGURE_RENDERER,
        STEP2B_RERUN_FIGURE_RENDERER,
        *RETAINED_FIGURE_SPECS,
        *(spec["report"] for spec in RETAINED_FIGURE_SPECS.values()),
    }
    missing = sorted(required - set(payloads))
    if missing:
        raise SystemExit(
            f"{artifact.name} is missing retained comparison figure inputs: {missing}"
        )
    renderer = payloads[RETAINED_FIGURE_RENDERER]
    if hashlib.sha256(renderer).hexdigest() != RETAINED_FIGURE_RENDERER_SHA256:
        raise SystemExit(
            f"{artifact.name}:{RETAINED_FIGURE_RENDERER} differs from the "
            "reviewed renderer"
        )
    renderer_text = renderer.decode("utf-8")
    renderer_markers = (
        "case_a_local_pressure_4x36x32_dt0p001.report.json",
        "case_b_local_pressure_2x36x48_dt0p002.report.json",
        STEP0A_RETAINED_REPORT_SHA256,
        "d409ecaac0c0abf418fce3ab2f0549979d38e02b9d73381d56982f7fc4e3bf14",
        STEP0A_RETAINED_APP_REF,
        "e07993bcf1166bd20eb87370c0b458552753e7ee",
        APPROVED_PUBLIC_CORE_REF,
        HISTORICAL_RETAINED_CORE_REF,
        "closed_multiblock_disk",
        "legacy-inferred",
        "10.5281/zenodo.14260459",
    )
    if any(marker not in renderer_text for marker in renderer_markers):
        raise SystemExit(
            f"{artifact.name}:{RETAINED_FIGURE_RENDERER} lacks exact source labels"
        )
    _validate_step0a_retained_comparison(payloads, artifact)
    step2b_renderer = payloads[STEP2B_FIGURE_RENDERER]
    if (
        hashlib.sha256(step2b_renderer).hexdigest()
        != STEP2B_FIGURE_RENDERER_SHA256
    ):
        raise SystemExit(
            f"{artifact.name}:{STEP2B_FIGURE_RENDERER} differs from the "
            "reviewed renderer"
        )
    step2b_renderer_text = step2b_renderer.decode("utf-8")
    step2b_renderer_markers = (
        "step2_case_b_std_kappa_2x20x17_dt0p001.report.json",
        "098e316daaea369a2a595cf43829d28597e53d2ff5a38cf32388e01c8dfa74aa",
        "23312a5e0147544eb9a4e6de004a166ada2722b70d3d39742f93aacd8a0fa0e6",
        "6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1",
        "e9b7d9084b24f7098170a221061eb159d0b090c1",
        APPROVED_PUBLIC_CORE_REF,
        "10.5281/zenodo.14260459",
        "paper_relative_discrepancy",
    )
    if any(marker not in step2b_renderer_text for marker in step2b_renderer_markers):
        raise SystemExit(
            f"{artifact.name}:{STEP2B_FIGURE_RENDERER} lacks exact source labels"
        )
    step2b_stdout = payloads[STEP2B_RAW_STDOUT]
    if hashlib.sha256(step2b_stdout).hexdigest() != STEP2B_RAW_STDOUT_SHA256:
        raise SystemExit(
            f"{artifact.name}:{STEP2B_RAW_STDOUT} differs from the reviewed "
            "normalized stdout"
        )
    tbar_prefix = b"tbar: injected Laplace field "
    if any(
        line.startswith(tbar_prefix + b"/")
        for line in step2b_stdout.splitlines()
    ):
        raise SystemExit(
            f"{artifact.name}:{STEP2B_RAW_STDOUT} contains a machine-local path"
        )
    for marker in (
        b"step 1000 t=1.000s",
        b"elapsed 1144.7s",
        b"saved -> runs/step2b/coarse_nt2_core20_rad17/mpi4_ga_full.npz",
    ):
        if marker not in step2b_stdout:
            raise SystemExit(
                f"{artifact.name}:{STEP2B_RAW_STDOUT} lacks completion marker "
                f"{marker!r}"
            )
    for name, spec in RETAINED_FIGURE_SPECS.items():
        renderer_name = spec.get("renderer", RETAINED_FIGURE_RENDERER)
        if renderer_name not in {
            RETAINED_FIGURE_RENDERER,
            STEP2B_FIGURE_RENDERER,
            TIP_REFINE_FIGURE_RENDERER,
            STEP2B_RERUN_FIGURE_RENDERER,
        }:
            raise SystemExit(
                f"{artifact.name}:{name} names an unreviewed figure renderer"
            )
        if renderer_name == STEP2B_RERUN_FIGURE_RENDERER:
            rerun_renderer = payloads[STEP2B_RERUN_FIGURE_RENDERER]
            if (
                hashlib.sha256(rerun_renderer).hexdigest()
                != STEP2B_RERUN_FIGURE_RENDERER_SHA256
            ):
                raise SystemExit(
                    f"{artifact.name}:{STEP2B_RERUN_FIGURE_RENDERER} differs "
                    "from the reviewed renderer"
                )
            rerun_renderer_text = rerun_renderer.decode("utf-8")
            rerun_renderer_markers = (
                "step2b_current_rerun_comparison.report.json",
                "coupfe-cardiac-step2b-current-rerun-comparison-v1",
                "load_publisher_reference",
                "team_times_s",
                "load_reviewed_corrected_run",
                "load_reviewed_legacy_report",
                "63a8de59b7b8b9ab309896ff69989d6ff89f6dfe2532151605486ad67967dd41",
                "svg.hashsalt",
                "10.5281/zenodo.14260459",
            )
            if any(
                marker not in rerun_renderer_text
                for marker in rerun_renderer_markers
            ):
                raise SystemExit(
                    f"{artifact.name}:{STEP2B_RERUN_FIGURE_RENDERER} lacks "
                    "exact source labels"
                )
        if renderer_name == TIP_REFINE_FIGURE_RENDERER:
            tip_renderer = payloads[TIP_REFINE_FIGURE_RENDERER]
            if (
                hashlib.sha256(tip_renderer).hexdigest()
                != TIP_REFINE_FIGURE_RENDERER_SHA256
            ):
                raise SystemExit(
                    f"{artifact.name}:{TIP_REFINE_FIGURE_RENDERER} differs "
                    "from the reviewed renderer"
                )
            tip_renderer_text = tip_renderer.decode("utf-8")
            tip_renderer_markers = (
                "step0b_tip6p0_full_cycle_comparison.report.json",
                "coupfe-cardiac-step0b-tip6p0-full-cycle-comparison-v2",
                "tip_refine",
                "10.5281/zenodo.14260459",
            )
            if any(
                marker not in tip_renderer_text
                for marker in tip_renderer_markers
            ):
                raise SystemExit(
                    f"{artifact.name}:{TIP_REFINE_FIGURE_RENDERER} lacks "
                    "exact source labels"
                )
        report_digest = hashlib.sha256(payloads[spec["report"]]).hexdigest()
        if report_digest != spec["report_sha256"]:
            raise SystemExit(
                f"{artifact.name}:{spec['report']} differs from the exact report "
                f"used by {name}"
            )
        payload = payloads[name]
        if (
            len(payload) != spec["size_bytes"]
            or hashlib.sha256(payload).hexdigest() != spec["sha256"]
        ):
            raise SystemExit(
                f"{artifact.name}:{name} differs from the reviewed SVG payload"
            )
        _validate_retained_svg(payload, name, spec, artifact)


def _core_dependency_entry(text: str, artifact: Path, member: str) -> str:
    if member == "pyproject.toml":
        entries = [
            line
            for line in text.splitlines()
            if re.match(
                r"^\s*[\"']coupfe(?:\[[^\]\r\n]+\])?(?:\s|@|[<>=!~])",
                line,
                flags=re.IGNORECASE,
            )
        ]
    else:
        entries = [
            line
            for line in text.splitlines()
            if re.match(
                r"^Requires-Dist:\s*coupfe(?:\[[^\]\r\n]+\])?(?:\s|@|[<>=!~])",
                line,
                flags=re.IGNORECASE,
            )
        ]
    if len(entries) != 1:
        raise SystemExit(
            f"{artifact.name}:{member} must declare exactly one CoupFE dependency"
        )
    return entries[0]


def _extract_core_pin(text: str, artifact: Path, member: str) -> str:
    dependency = _core_dependency_entry(text, artifact, member)
    pattern = re.compile(
        r"(?i)\bcoupfe(?:\[[^\]]+\])?\s*@\s*git\+https://github\.com/"
        r"tengzhang48/CoupFE\.git@([0-9a-f]{40})(?=[\s\"',;)]|$)"
    )
    pins = pattern.findall(dependency)
    if len(pins) != 1:
        raise SystemExit(
            f"{artifact.name}:{member} must contain exactly one CoupFE dependency "
            "pinned to the full 40-hex revision on the public HTTPS repository"
        )
    return pins[0].casefold()


def _validate_core_pin(
    text: str,
    artifact: Path,
    member: str,
    *,
    allow_unapproved_core_ref: bool,
) -> str:
    pin = _extract_core_pin(text, artifact, member)
    approved = APPROVED_PUBLIC_CORE_REF
    if approved is not None and not re.fullmatch(r"[0-9a-f]{40}", approved):
        raise SystemExit(
            "release guard configuration error: APPROVED_PUBLIC_CORE_REF must "
            "be None or one lowercase full 40-hex revision"
        )
    if approved is None:
        message = (
            f"{artifact.name}:{member} pins CoupFE revision {pin}, but the "
            "release guard has no approved public Core revision yet; set "
            "APPROVED_PUBLIC_CORE_REF to the final publicly reachable clean-root "
            "commit before release"
        )
    elif pin != approved:
        message = (
            f"{artifact.name}:{member} pins unapproved CoupFE revision {pin}; "
            f"the only approved public Core revision is {approved}"
        )
    else:
        return pin

    if not allow_unapproved_core_ref:
        raise SystemExit(message)
    print(f"WARNING: {message} (audit override enabled)")
    return pin


def _validate_source_tree(
    source_root: Path,
    *,
    allow_unapproved_core_ref: bool,
    allow_untracked_required: bool,
    allow_dirty_source: bool,
) -> tuple[int, str]:
    """Check tracked and non-ignored untracked release inputs.

    A forbidden file already deleted in the working tree is tolerated so the
    guard can inspect an uncommitted cleanup. It remains rejected if present.
    """

    status_result = subprocess.run(
        [
            "git",
            "-C",
            str(source_root),
            "status",
            "--porcelain=v1",
            "-z",
            "--untracked-files=all",
        ],
        check=True,
        capture_output=True,
    )
    if status_result.stdout:
        message = f"{source_root} is not a clean Git worktree"
        if not allow_dirty_source:
            raise SystemExit(message)
        print(f"WARNING: {message} (audit override enabled)")

    tracked_result = subprocess.run(
        ["git", "-C", str(source_root), "ls-files", "--cached", "-z"],
        check=True,
        capture_output=True,
    )
    untracked_result = subprocess.run(
        [
            "git",
            "-C",
            str(source_root),
            "ls-files",
            "--others",
            "--exclude-standard",
            "-z",
        ],
        check=True,
        capture_output=True,
    )
    tracked = {
        item.decode("utf-8", errors="strict")
        for item in tracked_result.stdout.split(b"\0")
        if item
    }
    untracked = {
        item.decode("utf-8", errors="strict")
        for item in untracked_result.stdout.split(b"\0")
        if item
    }
    names = sorted(tracked | untracked)
    _validate_names(names, source_root)

    missing = sorted(
        name
        for name in names
        if not (source_root / PurePosixPath(name)).exists()
        and not (source_root / PurePosixPath(name)).is_symlink()
        and not _is_forbidden_path(name)
    )
    if missing:
        raise SystemExit(
            f"{source_root} has tracked releasable paths missing from the "
            f"working tree: {missing}"
        )

    available = {
        name
        for name in names
        if (source_root / PurePosixPath(name)).exists()
        or (source_root / PurePosixPath(name)).is_symlink()
    }
    symlinks = sorted(
        name
        for name in available
        if (source_root / PurePosixPath(name)).is_symlink()
    )
    if symlinks:
        raise SystemExit(f"{source_root} contains symbolic links: {symlinks}")
    non_files = sorted(
        name
        for name in available
        if not (source_root / PurePosixPath(name)).is_file()
    )
    if non_files:
        raise SystemExit(
            f"{source_root} contains tracked submodules or special paths: {non_files}"
        )

    _validate_public_source_inventories(
        available,
        source_root,
        root_files=ROOT_SOURCES,
    )
    _require_files(available, PUBLIC_RELEASE_FILES, source_root)
    untracked_required = sorted(
        PUBLIC_RELEASE_FILES - tracked
    )
    if untracked_required:
        message = (
            f"{source_root} has required release inputs that are not tracked: "
            f"{untracked_required}"
        )
        if not allow_untracked_required:
            raise SystemExit(message)
        print(f"WARNING: {message} (audit override enabled)")

    _reject_forbidden_files(available, source_root)
    payloads = {}
    for name in sorted(available):
        payload = (source_root / PurePosixPath(name)).read_bytes()
        payloads[name] = payload
        _validate_text(
            name,
            payload,
            source_root,
        )
    _require_markers(payloads, LEGAL_MARKERS, source_root)
    _require_markers(payloads, CORE_VERIFIER_MARKERS, source_root)
    _validate_retained_result(payloads, source_root)
    _validate_truncated_polar_archive(payloads, source_root)
    _validate_retained_figures(payloads, source_root)
    _validate_step0b_prefix_diagnostic(payloads, source_root)
    _require_markers(
        payloads,
        {
            "pyproject.toml": (
                f'license = "{LICENSE_EXPRESSION}"',
                '"LICENSES/*.txt"',
            ),
            "examples/cardiac_benchmark/activation.py": (
                "SPDX-License-Identifier: CC-BY-4.0",
                "Changes made for CoupFE-Cardiac",
                "10.5281/zenodo.10875818",
                "be92da5dbc1fd26d424bf88ef7db13b4",
                "325d17d850c2e2032abb85a4191a5795d3008ab7",
            ),
            "examples/cardiac_benchmark/fiber_crosscheck.py": (
                "SPDX-License-Identifier: MIT",
                "two incorporated third-party source",
            ),
        },
        source_root,
    )
    core_ref = _validate_core_pin(
        (source_root / "pyproject.toml").read_text(encoding="utf-8"),
        source_root,
        "pyproject.toml",
        allow_unapproved_core_ref=allow_unapproved_core_ref,
    )
    return len(available), core_ref


def _validate_wheel(
    wheel: Path,
    *,
    allow_unapproved_core_ref: bool,
) -> tuple[int, str]:
    with zipfile.ZipFile(wheel) as archive:
        members = archive.infolist()
        names = [member.filename for member in members]
        _validate_names(names, wheel)

        symlinks = sorted(
            member.filename
            for member in members
            if stat.S_ISLNK((member.external_attr >> 16) & 0xFFFF)
        )
        if symlinks:
            raise SystemExit(f"{wheel.name} contains symbolic links: {symlinks}")

        files = {member.filename for member in members if not member.is_dir()}
        metadata_roots = {
            PurePosixPath(name).parts[0]
            for name in files
            if len(PurePosixPath(name).parts) == 2
            and PurePosixPath(name).parts[0].endswith(".dist-info")
            and PurePosixPath(name).name == "METADATA"
        }
        if len(metadata_roots) != 1:
            raise SystemExit(
                f"{wheel.name} must contain exactly one .dist-info/METADATA file"
            )
        dist_info = metadata_roots.pop()
        non_metadata = sorted(
            name
            for name in files
            if not PurePosixPath(name).parts
            or PurePosixPath(name).parts[0] != dist_info
        )
        if non_metadata:
            raise SystemExit(
                f"{wheel.name} violates the metadata-only wheel policy: "
                f"{non_metadata}"
            )

        required = {
            f"{dist_info}/METADATA",
            f"{dist_info}/RECORD",
            f"{dist_info}/WHEEL",
            f"{dist_info}/licenses/LICENSE",
            f"{dist_info}/licenses/LICENSES/CC-BY-4.0.txt",
            f"{dist_info}/licenses/LICENSES/Reidmen-cardiac_benchmark_toolkit-MIT.txt",
            f"{dist_info}/licenses/NOTICE",
            f"{dist_info}/licenses/THIRD_PARTY_NOTICES.md",
            f"{dist_info}/licenses/docs/LICENSE.md",
            f"{dist_info}/top_level.txt",
        }
        _require_files(files, required, wheel)
        unexpected = sorted(files - required)
        if unexpected:
            raise SystemExit(
                f"{wheel.name} contains unexpected metadata-only wheel files: "
                f"{unexpected}"
            )
        _reject_forbidden_files(files, wheel)
        payloads = {}
        for member in members:
            if not member.is_dir():
                payload = archive.read(member)
                payloads[member.filename] = payload
                _validate_text(member.filename, payload, wheel)

        wheel_legal_markers = {
            f"{dist_info}/licenses/{name}": markers
            for name, markers in LEGAL_MARKERS.items()
        }
        _require_markers(payloads, wheel_legal_markers, wheel)

        metadata_name = f"{dist_info}/METADATA"
        metadata = archive.read(metadata_name).decode("utf-8")
        if f"License-Expression: {LICENSE_EXPRESSION}" not in metadata:
            raise SystemExit(
                f"{wheel.name}:{metadata_name} is missing the expected license "
                "expression"
            )
        core_ref = _validate_core_pin(
            metadata,
            wheel,
            metadata_name,
            allow_unapproved_core_ref=allow_unapproved_core_ref,
        )
    return len(files), core_ref


def _validate_sdist(
    sdist: Path,
    *,
    allow_unapproved_core_ref: bool,
) -> tuple[int, str]:
    with tarfile.open(sdist, mode="r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        _validate_names(names, sdist)

        unsupported = sorted(
            member.name
            for member in members
            if not (member.isfile() or member.isdir())
        )
        if unsupported:
            raise SystemExit(
                f"{sdist.name} contains links or special-file entries: {unsupported}"
            )

        roots = {
            PurePosixPath(member.name).parts[0]
            for member in members
            if PurePosixPath(member.name).parts
        }
        if len(roots) != 1:
            raise SystemExit(f"{sdist.name} must have one top-level directory")
        root = roots.pop()
        file_members = {
            PurePosixPath(*PurePosixPath(member.name).parts[1:]).as_posix(): member
            for member in members
            if member.isfile()
            and len(PurePosixPath(member.name).parts) > 1
            and PurePosixPath(member.name).parts[0] == root
        }
        files = set(file_members)
        missing = sorted(EXPECTED_SDIST_FILES - files)
        unexpected = sorted(files - EXPECTED_SDIST_FILES)
        if missing or unexpected:
            raise SystemExit(
                f"{sdist.name} source-archive inventory mismatch: "
                f"missing={missing}, unexpected={unexpected}"
            )
        _reject_forbidden_files(files, sdist)

        payloads: dict[str, bytes] = {}
        for name, member in file_members.items():
            stream = archive.extractfile(member)
            if stream is None:
                raise SystemExit(f"{sdist.name}:{name} could not be read")
            payload = stream.read()
            payloads[name] = payload
            _validate_text(name, payload, sdist)
        _require_markers(payloads, LEGAL_MARKERS, sdist)
        _require_markers(payloads, CORE_VERIFIER_MARKERS, sdist)
        _validate_retained_result(payloads, sdist)
        _validate_truncated_polar_archive(payloads, sdist)
        _validate_retained_figures(payloads, sdist)
        _validate_step0b_prefix_diagnostic(payloads, sdist)
        _require_markers(
            payloads,
            {
                "pyproject.toml": (
                    f'license = "{LICENSE_EXPRESSION}"',
                    '"LICENSES/*.txt"',
                ),
                "PKG-INFO": (
                    f"License-Expression: {LICENSE_EXPRESSION}",
                    "License-File: LICENSES/CC-BY-4.0.txt",
                    "License-File: LICENSES/Reidmen-cardiac_benchmark_toolkit-MIT.txt",
                ),
                "examples/cardiac_benchmark/activation.py": (
                    "SPDX-License-Identifier: CC-BY-4.0",
                    "Changes made for CoupFE-Cardiac",
                    "10.5281/zenodo.10875818",
                    "be92da5dbc1fd26d424bf88ef7db13b4",
                    "325d17d850c2e2032abb85a4191a5795d3008ab7",
                ),
                "examples/cardiac_benchmark/fiber_crosscheck.py": (
                    "SPDX-License-Identifier: MIT",
                    "two incorporated third-party source",
                ),
            },
            sdist,
        )

        pyproject_pin = _validate_core_pin(
            payloads["pyproject.toml"].decode("utf-8"),
            sdist,
            "pyproject.toml",
            allow_unapproved_core_ref=allow_unapproved_core_ref,
        )
        pkg_info_pin = _validate_core_pin(
            payloads["PKG-INFO"].decode("utf-8"),
            sdist,
            "PKG-INFO",
            allow_unapproved_core_ref=allow_unapproved_core_ref,
        )
        if pyproject_pin != pkg_info_pin:
            raise SystemExit(
                f"{sdist.name} records different Core revisions in pyproject.toml "
                "and PKG-INFO"
            )
    return len(files), pyproject_pin


def validate(
    dist_dir: Path,
    source_root: Optional[Path] = None,
    *,
    allow_unapproved_core_ref: bool = False,
    allow_untracked_required: bool = False,
    allow_dirty_source: bool = False,
) -> None:
    source_result = (
        _validate_source_tree(
            source_root.resolve(),
            allow_unapproved_core_ref=allow_unapproved_core_ref,
            allow_untracked_required=allow_untracked_required,
            allow_dirty_source=allow_dirty_source,
        )
        if source_root is not None
        else None
    )
    wheels = sorted(dist_dir.glob("*.whl"))
    sdists = sorted(dist_dir.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise SystemExit(
            f"expected one wheel and one sdist in {dist_dir}, found "
            f"{len(wheels)} wheel(s) and {len(sdists)} sdist(s)"
        )

    wheel_count, wheel_ref = _validate_wheel(
        wheels[0], allow_unapproved_core_ref=allow_unapproved_core_ref
    )
    sdist_count, sdist_ref = _validate_sdist(
        sdists[0], allow_unapproved_core_ref=allow_unapproved_core_ref
    )
    refs = {
        "wheel": wheel_ref,
        "sdist": sdist_ref,
    }
    if source_result is not None:
        refs["source"] = source_result[1]
    if len(set(refs.values())) != 1:
        rendered = ", ".join(
            f"{kind}={ref}" for kind, ref in sorted(refs.items())
        )
        raise SystemExit(
            f"release inputs and artifacts record different Core revisions: "
            f"{rendered}"
        )

    source_count = source_result[0] if source_result is not None else None
    source_summary = (
        f"source tree ({source_count} files), " if source_count is not None else ""
    )
    print(
        f"validated {source_summary}metadata-only {wheels[0].name} "
        f"({wheel_count} files) and {sdists[0].name} ({sdist_count} files)"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dist_dir",
        nargs="?",
        type=Path,
        default=Path("dist"),
        help="directory containing exactly one wheel and one .tar.gz sdist",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path("."),
        help="git worktree checked before the artifacts (default: current directory)",
    )
    parser.add_argument(
        "--artifacts-only",
        action="store_true",
        help="skip the source-tree check",
    )
    parser.add_argument(
        "--allow-unapproved-core-ref-for-audit",
        "--allow-private-core-pin-for-audit",
        dest="allow_unapproved_core_ref_for_audit",
        action="store_true",
        help=(
            "inspect source/artifacts before the final public Core ref is "
            "approved (the older --allow-private-core-pin-for-audit spelling "
            "is retained as an alias); the result is not publishable"
        ),
    )
    parser.add_argument(
        "--allow-untracked-required-for-audit",
        action="store_true",
        help=(
            "inspect a review worktree whose required new release files are not "
            "yet tracked; the source result is not publishable"
        ),
    )
    parser.add_argument(
        "--allow-dirty-source-for-audit",
        action="store_true",
        help=(
            "inspect an intentionally modified review worktree; the source "
            "result is not publishable"
        ),
    )
    args = parser.parse_args()
    validate(
        args.dist_dir,
        source_root=None if args.artifacts_only else args.source_root,
        allow_unapproved_core_ref=args.allow_unapproved_core_ref_for_audit,
        allow_untracked_required=args.allow_untracked_required_for_audit,
        allow_dirty_source=args.allow_dirty_source_for_audit,
    )


if __name__ == "__main__":
    main()
