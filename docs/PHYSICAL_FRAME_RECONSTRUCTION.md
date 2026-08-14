# Physical-coordinate structural-frame reconstruction

## Finding

The closed multiblock mesh uses the benchmark toolkit's straight Cartesian
wall ruling.  An interior point of that ruling is generally not the point on
the interpolated ellipsoid identified by the mesh generator's stored
`(t, mu, theta)` values.  The old `gp-direct` path nevertheless interpolated
those stored parameters and passed them to `geometry.fiber_frame`.  The same
assumption was used to rebuild nodal and element frames after injecting a
Laplace transmural field, and the formula was duplicated in the serial and MPI
drivers.

The pinned toolkit instead evaluates the Laplace coordinate at the physical
point, reconstructs

```text
u = atan2(sqrt(y*y + z*z) / r_short(tbar), x / r_long(tbar))
v = 0 if u < 1e-7 else pi - atan2(z, -y)
```

and constructs fiber, sheet normal, and sheet from that reconstructed local
basis.  The source is
`Reidmen/cardiac_benchmark_toolkit@e8d47553cfc83eb274eba3e177de33148e7f441c`,
`ellipsoid_fiber_generation.py`.

## Implemented contract

- `structural_directions.py` is the single MIT-licensed NumPy adaptation of
  the pinned physical-coordinate formula, including its apex branch and its
  complete `(fiber, sheet, sheet_normal)` convention.
- On a closed Hex8, direct Gauss-point sampling evaluates
  `tbar_gp = N @ tbar_nodes` and `X_gp = N @ X_nodes`, then reconstructs the
  frame from `(tbar_gp, X_gp)`.
- Initial and Laplace-injected closed-mesh nodal frames use each physical node;
  stored element frames use the physical Q1 center.
- Serial and MPI state initialization call the same geometry evaluator.
- The historical `polar_ring` topology retains its prior parametric rule and
  apex fallbacks.  Archives identify the choice as
  `historical-parametric-mu-theta-v1`; closed archives identify
  `toolkit-physical-coordinate-u-v-v1`.

## Direction-change audit

This is a field-only audit; it did not run a mechanical trajectory.  It used
the corrected straight-wall meshes and their matching Q1 Laplace fields:

| Mesh | Nodes / Hex8 | Laplace field SHA-256 | max `|tbar-layer|` |
|---|---:|---|---:|
| `2x20x17` | 5,403 / 3,520 | `de24749a85b458c039a16a1e4b24422cf35d54cf853a40daa763e3137cb930a4` | 0.161807 |
| `4x36x32` | 29,885 / 23,616 | `7fbdbcf9a6b6c5135ef87fe998ec23a3dbf44957ccd3184584e8f6d60768c2b6` | 0.157794 |

Angles below compare the old stored-parametric GP rule with the new physical
rule over all 8 Gauss points.  Sheet and normal use sign-insensitive axis
angles because the old convention was almost exactly the negative of the
toolkit convention.

| Mesh / direction | median | p90 | p99 | maximum |
|---|---:|---:|---:|---:|
| `2x20x17` fiber | 0.0197° | 0.0906° | 0.3860° | 2.8229° |
| `2x20x17` sheet axis | 0.0403° | 0.2559° | 0.8579° | 2.8874° |
| `2x20x17` normal axis | 0.0515° | 0.2830° | 0.8995° | 1.5173° |
| `4x36x32` fiber | 0.0205° | 0.1142° | 0.3217° | 2.6916° |
| `4x36x32` sheet axis | 0.0465° | 0.2731° | 0.7767° | 2.7228° |
| `4x36x32` normal axis | 0.0574° | 0.3004° | 0.7924° | 0.9383° |

The directed mean changes for sheet and normal are about 179.9° on both
meshes: the prior `asign=-1` mapping aligned the fiber with the toolkit but
left sheet and sheet normal oppositely oriented.  The new code reproduces the
whole toolkit frame.  In the current Holzapfel--Ogden law, a global
`sheet -> -sheet` leaves the stress unchanged: `ss` is unchanged, while both
`fssym` and `I8fs` change sign in their product.  The small axis rotations are
the material change attributable to physical-coordinate reconstruction.

The largest nodal fiber change is 8.84° at the apex, where the old code used an
arbitrary `mu=-pi+0.05, theta=0` fallback.  It is highly localized: the nodal
fiber p99 is 0.172° on the coarse mesh and 0.257° on the fine mesh.  Direct
Gauss-point maxima occur in apex-adjacent cells and decrease slightly with
refinement.

## Lessons learned

1. A mesh generator's construction parameters are not an inverse coordinate
   map after a Cartesian wall interpolation; reconstruct from the physical
   point when that is the reference algorithm.
2. Matching only fiber dot products can hide an opposite sheet/normal
   convention.  Verify the complete frame and distinguish directed-vector
   angles from sign-insensitive material-axis angles.
3. Nodal, element-center, serial-GP, and MPI-GP paths must use one declared
   reconstruction rule.  Duplicated scientific formulas are a provenance and
   correctness risk.
4. Preserve the historical topology explicitly rather than silently changing
   old diagnostic examples.
5. Record the reconstruction identifier in every new archive.  A mesh and a
   Laplace-field hash alone do not identify the structural-direction field.

The measured fiber-axis change is small over almost the entire domain, so this
fix alone is not evidence that it explains a large Case B curve discrepancy.
It is nevertheless required for source fidelity.  The appropriate next gate
is a short, snap-through-focused rerun before considering a full 1 s run.
