# Mesh refinement: uniform refinement versus tip grading

**Status:** reference document. Read this before interpreting any mesh label,
figure, or refinement number in this repository.

This project varies the mesh along **three independent axes**. Two of them are
genuine refinement (more elements). The third is **node relocation at a fixed
element count** and is *not* refinement, despite being named `tip_refine`.
Results from the three axes point in different directions, so attributing a
number to the wrong axis inverts its meaning. This document fixes the
vocabulary, tabulates every retained result by axis, and lists the specific
places where the existing naming is misleading.

---

## 1. The three axes

A closed-mesh label is `A x B x C`, optionally followed by `tip F`.

| Axis | Parameter | Label position | What it changes | Is it refinement? |
|---|---|---|---|---|
| **Through-wall** | `n_t` (`--nt`) | `A` | Hex8 element layers across the wall thickness | **Yes** — more elements |
| **Surface** | `n_core` (`--ncore`), `n_radial` (`--nradial`) | `B`, `C` | In-plane resolution of the endo/epi surfaces | **Yes** — more elements |
| **Tip grading** | `tip_refine` (`--tip-refine`) | `tip F` suffix | Moves existing nodes along the meridian toward the apex | **No** — element count is unchanged |

`B` is elements per side of the central square containing the apex; `C` is
radial element layers in each of the four outer square-to-disk blocks. Element
count is `n_t * (n_core^2 + 4*n_core*n_radial)`.

Absence of a `tip` token means `tip_refine = 1.0`, which reproduces the uniform
mesh bit-for-bit (gated by
`tests/test_cardiac_closed_geometry.py::test_tip_refine_one_is_the_uniform_mesh`).
In filenames the decimal point becomes `p`: `tip6p0`.

### Mesh sizes

| Mesh | Elements | Mean epicardial edge (mm) | Wall layer (mm) |
|---|---:|---:|---:|
| `2x20x17` | 3,520 | 4.00 | ~4.5 |
| `2x36x32` | 11,808 | 2.17 | ~4.5 |
| `4x32x48` | 28,672 | 1.91 | ~2.3 |
| `6x64x72` | 135,168 | 1.09 | ~1.5 |
| local FEniCS reference (P1 tets, P2 field) | 17,648 tets | 3.88 | -- |

The paper's nominal target edge length is approximately 5 mm, so `2x20x17` is
the mesh closest to the benchmark's own scale and `2x36x32` is already about
2.3x finer than nominal.

---

## 2. What `tip_refine` actually does

`tip_refine` is **global meridional node grading (r-adaptation)**, not local
h-refinement. It relocates nodes; it creates none.

From `examples/cardiac_benchmark/geometry.py:485-496`, the meridian coordinate
`rho` (0 at the apex, 1 at the base rim) is remapped before the wall point is
computed:

```python
def _cluster_meridian_to_apex(rho, strength):
    a = 1.0 / strength
    return rho * (a + (1.0 - a) * rho)
```

Properties, all verifiable from that expression:

- `g(0) = 0` and `g(1) = 1`: the apex and the base rim are fixed points.
- `g'(rho) = a + 2(1-a)*rho`, so meridian spacing next to the apex is scaled by
  `1/F`, while spacing at the base rim **grows** by up to `2 - 1/F`.
- Element count, connectivity, boundary labels, wall ruling and base rim are
  all unchanged.

Valid range is `[1.0, 8.0]`; the default is `1.0`.

**The name is misleading and the code says so in two places.** At `tip 6.0` on
`2x20x17`, **5,160 of 5,403 nodes move**, with a median displacement of 21.95 mm
and a maximum of 29.07 mm. That is essentially the whole mesh, not a tip. The
`--tip-refine` CLI help string ("local refinement near the apex") and the
`geometry.py:517` docstring ("changes only the local element size near the
apex") both overstate locality; `docs/CASE_B_MESH_ERROR_LAYERS.md:550-553`
corrects them: "`tip_refine` moves nodes throughout the meridian at fixed
element count. It is global r-adaptation, not local h-refinement."

Read `tip_refine` as **"cluster the meridian toward the apex, globally"**.

---

## 3. Results by axis

All values are 0.32 s prefix endpoints at landmark `p0` unless stated.
Reference points: local FEniCS `p0 x = -16.06 mm`; the ten-team envelope at
0.32 s is `p0 x` in `[-17.15, -15.37]` mm.

### 3a. Surface refinement — the dominant lever, moves AWAY from the reference

| Change | `p0` endpoint effect | `p1` endpoint effect |
|---|---:|---:|
| `2x20x17 -> 2x36x32` | **9.994 mm** | **8.772 mm** |

`p0 x` goes `-22.15 -> -32.10`, i.e. **further from** FEniCS `-16.06`. This is
the largest single mesh effect measured, roughly 15x the through-wall effect.

The mechanism is understood: surface refinement simultaneously changes the Q1
resolution **and** the discrete epicardial Robin normal. The spurious long-axis
rotational restraint from facet normals falls from `1.285103` to `0.396153`
N m/rad (epicardial part) under exactly this change, and the two are not
separable by refining alone. See `CASE_B_MESH_ERROR_LAYERS.md`.

### 3b. Through-wall refinement — small at the prefix, phase-mixed at full cycle

| Change | `p0` endpoint effect | `p1` endpoint effect |
|---|---:|---:|
| `2x20x17 -> 4x20x17` (coarse surface) | 0.669 mm | 0.635 mm |
| `2x36x32 -> 4x36x32` (fine surface) | 0.881 mm | 1.217 mm |

At full cycle, with `tip_refine = 6.0` held fixed, `n_t` 2 -> 4:

| Metric | `p0` | `p1` |
|---|---:|---:|
| full-cycle RMSE between the two meshes | 0.4076 mm | 0.3618 mm |
| max transient difference | 1.633 mm @0.250 s | 1.406 mm @0.247 s |
| RMSE vs FEniCS | 1.0927 -> **1.0169** mm (-6.9%) | 1.1685 -> **1.0730** mm (-8.2%) |
| max gap vs FEniCS | 5.552 -> **4.188** mm (-24.6%) | 5.024 -> **3.753** mm (-25.3%) |
| relaxation-phase RMSE | **worsens 16.2%** | **worsens 14.9%** |

This is phase-mixed: better through the snap, worse during relaxation. It is
numerical sensitivity, not convergence.

### 3c. Tip grading — moves TOWARD the reference at the coarse surface

At `2x20x17`, 0.32 s, facet normals, current frame:

| Run | `p0` gap vs FEniCS | `p0 x` | RMSE 0-0.32 s | max history error |
|---|---:|---:|---:|---:|
| uniform (`F = 1.0`) | 6.154 mm | -22.15 | 2.593 mm | 4.19 mm |
| `tip 2.5` | 3.608 mm | -19.59 | 1.535 mm | 4.93 mm |
| `tip 4.0` | 2.587 mm | -18.54 | 1.421 mm | 5.55 mm |
| `tip 6.0` | **1.998 mm** | **-17.92** | 1.424 mm | 5.55 mm |

Two cautions that must travel with this table:

- **The endpoint improves monotonically while the maximum history error grows**
  (4.19 -> 5.55 mm). Grading improves where the trajectory ends, not
  everywhere along it. Reporting the endpoint alone overstates the benefit.
- **The effect nearly vanishes on a finer surface.** At `4x32x48`, `tip 2.5`
  changes the endpoint by only **0.23 mm** (`p0`) and **0.095 mm** (`p1`), and
  the FEniCS prefix RMSE *worsens* (5.638 -> 5.866 mm). The ~4 mm of endpoint
  movement seen at `2x20x17` is a coarse-surface phenomenon.

---

## 4. The comparison, stated plainly

**Uniform surface refinement and tip grading move the solution in opposite
directions.** From `CASE_B_MESH_ERROR_LAYERS.md:173-177`:

> **Direction check**: uniform surface refinement moves the endpoint *away*
> from FEniCS (-22.2 -> -32.1 x); tip refinement moves it *toward*
> (-22.2 -> -19.6 -> -18.5). The opposite trends are consistent with multiple
> coupled numerical sensitivities. These controls do not provide a causal
> decomposition of the historical non-monotone behavior.

| | Surface (uniform) | Through-wall (uniform) | Tip grading |
|---|---|---|---|
| Adds elements | yes | yes | **no** |
| Endpoint effect, coarse surface | ~10 mm | ~0.7 mm | ~4 mm over `F=1..6` |
| Direction vs FEniCS endpoint | **away** | toward (full cycle) | **toward** |
| Effect on snap timing | earlier, overshoots past FEniCS | slightly earlier | slightly earlier |
| Behaviour on a fine surface | -- | ~1 mm | **nearly inert (0.23 mm)** |
| Run to full cycle? | **never** | yes (one pair) | never varied at full cycle |

**They agree on one thing:** every refinement of any kind moves the snap
earlier. FEniCS crosses `u_z = -5 mm` at 0.2423 s; the uniform coarse mesh at
0.2487 s; tip grading brings this to 0.2450-0.2461 s and through-wall
refinement to 0.2450 s, while uniform *surface* refinement overshoots to
0.2368 s.

**The single most important caveat:** because grading moves the endpoint toward
the reference while adding no resolution, a smaller endpoint gap at higher `F`
is **not** evidence of a more accurate mesh. An endpoint trend is not
convergence. No mesh-independence evidence exists for this problem; the
declared convergence ladder (surface `{2x20x17, 2x36x32}` x wall `{2, 4}`) has
not been run, and it deliberately excludes grading.

---

## 5. Known naming traps

These are real inconsistencies in the current tree. Until they are fixed, a
reader should treat the following as unreliable and consult this document.

**Trap 1 — the headline full-cycle artifact is named for tip refinement but
varies wall layers.** The figure `docs/figures/step0b_tip_refine_full_cycle.svg`,
the report `results/step0b_tip6p0_full_cycle_comparison.report.json`, the
generator `compare_tip_refine_full_cycle.py`, and the README caption ("Current
Step 0 Case B tip-refined full-cycle displacement comparison") are all named
after tip refinement. **Both runs in that comparison hold `tip_refine = 6.0`
fixed and vary `n_t` from 2 to 4.** Verified directly in the report JSON. The
-6.9%/-8.2% RMSE and -24.6%/-25.3% max-gap improvements are therefore
**through-wall** results, not tip results. The SVG's own subtitle ("Two and
four wall layers") is correct; the filename and caption are not.

**Trap 2 — `2x20x17` and `4x20x17` name two different mesh pairs.** In
`CASE_B_STATUS.md` the same two labels denote the uniform clean gate
(`tip_refine = 1.0`, endpoints -22.15 / -22.28 mm) in one section and the
graded full-cycle pair (`tip_refine = 6.0`, endpoints -17.92 / -18.09 mm)
roughly 85 lines earlier. The endpoints differ by about 4 mm. Always check for
a `tip` token.

**Trap 3 — the reproduction verdict never states that its trajectory is
graded.** `BENCHMARK_REPRODUCTION_STATUS.md` reports the headline Step 0B
numbers without the words "tip", "graded", or "uniform" appearing in its prose.
Its trajectory is `tip_refine = 6.0`.

**Trap 4 — two containment statistics look contradictory but are different
meshes.** `CASE_B_MESH_ERROR_LAYERS.md` reports 21.0% / 22.2% all-component
envelope containment (the **two-layer** run);
`BENCHMARK_REPRODUCTION_STATUS.md` reports 15.3% / 26.6% (the **four-layer**
run). Both are correct; neither names its mesh.

**Trap 5 — "tip" is overloaded.** It means meridional grading in `tip_refine`,
and it also means the *truncated open apex* of the superseded polar-ring
geometry ("open-tip", `--apex-offset 0.2`). These are unrelated. Grepping
"tip" returns both.

**Trap 6 — "fine" names at least three meshes:** `4x36x32` (the four-way
split's "fine"), `4x32x48`, and `6x64x72` ("global fine"). `2x36x32` appears as
"surface-only", "fine surface", and bare `2x36x32`.

**Trap 7 — two three-number conventions share the same notation.** Closed-mesh
labels are `n_t x n_core x n_radial`. Superseded polar labels are
`n_t x n_mu x n_theta` (2x12x16, 2x24x32, 2x36x48, ...). `2x36x48` (polar,
open apex) and `2x36x32` (closed, surface-only) differ by one character and
denote different domains.

**Trap 8 — grading was never included in the mechanism analysis.** The
six-rigid-mode table and the snap-association table cover only `2x20x17`,
`2x36x32`, `4x32x48` and FEniCS. **No graded mesh appears in either.** Both
tables sit next to tip narrative, but the twist-restraint mechanism was never
measured on a graded mesh.

---

## 6. What is not established

Quoting the repository's own claim boundaries, which this document does not
soften:

- "These controls do not prove that the two contributors are exhaustive, do not
  uniquely allocate the remaining error, and do not establish distinct solution
  branches."
- "**Monotone endpoint sensitivity, not a convergence sequence.** ... the data
  do not establish an error order, apex ownership, or convergence."
- "One pair does not establish a vanishing grading effect, an apex-converged
  mesh, or a distinct solution branch."
- "Keep facet normals as the benchmark-comparable primary operator ... **This
  mesh-independence evidence does not yet exist for CoupFE.**"
- "Do not treat any retained trajectory (CoupFE coarse/fine or local FEniCS) as
  a converged reference; **none exists yet**."
- "**An endpoint trend is not convergence.** ... Report endpoint, history,
  phase, and transient extrema together."

Every completed full-cycle run in this repository has `tip_refine = 6.0`. The
surface axis, which carries the largest effect, has never been run to a full
cycle.

---

## 7. Reporting rules

When citing any mesh result in this project:

1. **Always write the full label including the `tip` token**, or explicitly
   write "uniform". A bare `2x20x17` is ambiguous.
2. **Write the tip strength with its meaning**, not as a bare number. A label
   like `tip_refine = 6.0` does not say what 6.0 measures; "tip refinement
   strength 6.0 (apex spacing 1/6 of uniform)" does. The value is a clustering
   factor, never an element count or a refinement level.
3. **Name the axis** that changed: surface, through-wall, or tip. Never write
   "refined" unqualified.
4. **Report endpoint, history RMSE, and maximum transient together.** Tip
   refinement improves the first and worsens the third.
5. **Do not call it local**, and do not imply it adds elements — it relocates
   nodes at a fixed element count.
6. **Do not present an endpoint trend as convergence.**
