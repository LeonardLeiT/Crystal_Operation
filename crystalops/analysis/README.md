# Structure analysis

`crystalops.analysis` performs numerical analysis using the shared data models
in `crystalops.core`. Results can be displayed by `crystalops.visualization`;
the numerical functions do not import a plotting backend.

The existing `aniso_alpha` and `aniso_lattice` functions remain available in
`crystalops.analysis.anisotropy`, with their original behavior retained.

## Grain-boundary curvature

The analysis APIs are:

- `identify_grain_boundary_atoms`: find unlike-grain contacts, grain pairs,
  and junction atoms using periodic neighbors.
- `reconstruct_grain_boundary_surfaces`: build shared interfaces from smooth
  Gaussian grain-occupancy fields and sample their geometry.
- `analyze_grain_boundary_curvature`: run identification, reconstruction,
  curvature evaluation, atom mapping, and area-weighted statistics.
- `mean_curvature_from_derivatives`: evaluate the mean curvature of an
  implicit surface from its gradient and Hessian.
- `segment_grains_with_ptm`: optionally recover current grain labels from
  OVITO PTM orientations and grain segmentation.

The grain-boundary implementation currently supports **orthogonal cells with
periodicity in all three directions**. It rejects other cells explicitly.
Grain labels must represent the current spatial grain assignment and cover
all atoms. Element types are never automatically treated as grain labels.

For neighboring grains `i < j`, define the interface by `q_j - q_i = 0` in
the region where those two grain fields dominate. The normal points from `i`
to `j`. Curvature uses `H = div(n)/2`; an outward-oriented sphere has
`H = +1/R`. The Gaussian length `smoothing_length` is a physical analysis
scale, not a target curvature. No curvature-flow relaxation is performed.

Analytic Gaussian derivatives supply curvature. An auxiliary sampling grid
locates the interface and approximates its area; `sampling_spacing` defaults
to `smoothing_length/3`. The Gaussian support is truncated at
`kernel_cutoff * smoothing_length`, default `5 * smoothing_length`. That
support must be smaller than half the shortest box edge for this implementation.

Junction competition and weak/undefined gradients invalidate local curvature.
The output records total and valid analyzed area, the excluded fraction,
local validity masks, and smoothing/sampling settings. Pairwise surface meshes
are intended for curvature quadrature and display; they are not a watertight
multi-region reconstruction suitable for genus or other topology measurements.

Atomic curvature is mapped from the nearest sampled point of the same
grain-pair interface within the boundary band. Its stored distance is to that
sample, rather than an exact closest-point distance to a continuous surface.
Unmapped or invalid atoms retain NaN and an explicit validity mask.

## File input and reproducible examples

LAMMPS `atomic` data does not store a grain-ID column in this reader. Supply a
CSV containing `atom_id,grain_id`, or explicitly request PTM grain recovery.
CSV rows are aligned by atom ID; duplicates, missing IDs, and extra IDs are
rejected. The original structure files are not rewritten.

From the repository root, analyze the supplied Cu structures with:

```bash
python -m example.crystalops.analysis.analyze_grain_boundary_curvature PSC.data --segment-grains --grain-algorithm minimum_spanning_tree --merging-threshold 3 --sigma 2.2 --spacing 0.75 --neighbor-cutoff 3.2 --output-dir analysis_results/PSC_sigma2p2
python -m example.crystalops.analysis.analyze_grain_boundary_curvature DSC.data --segment-grains --grain-algorithm minimum_spanning_tree --merging-threshold 3 --sigma 2.2 --spacing 0.75 --neighbor-cutoff 3.2 --output-dir analysis_results/DSC_sigma2p2
```

These commands illustrate parameter selection; convergence must be assessed
by holding the physical smoothing length fixed while refining spacing, then
checking sensitivity to the physical smoothing length separately. The supplied
Cu coordinates use angstroms, so `H` is in inverse angstroms and area is in
square angstroms. For other inputs, `--length-unit` labels the coordinate unit
without converting coordinates.

PTM recovery retains the pre-adoption core labels and the labels actually used
after noncrystalline/orphan atoms are assigned to neighboring grains. Its
parameters and diagnostics are saved. It does not use an ideal P/D surface or
the input filename to define the interface.

The supplied DSC structure contains coincident atoms, and its zero-misorientation
central regions can be merged or split by different segmentation methods. The
sample runs therefore retain label-recovery and tie-assignment sensitivity
checks. These are conditional geometric measurements, not unique intrinsic
atomic curvatures. The current nearest-core adoption is deterministic, resolves
equal distances by canonical grain ID, and retains the pre-adoption cores.

Prepared label files can be used to run the sample comparison independently of
segmentation:

```bash
python -m example.crystalops.analysis.compare_grain_boundary_curvature DSC.data PSC.data --labels-dir analysis_results/label_diagnostics --output-dir analysis_results/curvature --dsc-six-grain-labels analysis_results/label_diagnostics/DSC_graph_auto_6grain_labels.csv
```

## Outputs

`write_grain_boundary_results` in `crystalops.io.grain_boundary_results`
writes:

- `atoms.csv`: original atom IDs and coordinates, grain labels, boundary and
  junction flags, mapped curvature, validity, and mapping distance.
- `surface_faces.csv`: grain pair, sample location, area, local mean/Gaussian
  curvature, junction flag, and validity for each triangle.
- `surfaces.vtk`: an ASCII VTK PolyData mesh with the same cell measurements,
  readable by ParaView and compatible VTK readers.
- `summary.json`: parameters, label provenance, units, per-pair and network
  area-weighted statistics. Undefined JSON values are `null`.

The example additionally saves `grain_labels.csv`, including pre-adoption
grain cores and PTM diagnostics when segmentation was requested. CSV/VTK retain
NaN for undefined curvature together with explicit validity flags.

The main shape measures are area-weighted mean absolute curvature and RMS
curvature. A small signed average may result from cancellation. These
measurements alone do not establish that a grain-boundary network is a TPMS.

## Reconstruction calibration and affine self-audit

`mean_curvature` is the local geometric quantity `(k1+k2)/2`.
`mean_signed_H` is a network-wide signed area average; its value also depends
on grain-ID normal conventions. Positive and negative local curvature can
cancel even on a strongly nonminimal surface. The statistics therefore also
report absolute/RMS curvature and the cancellation fraction. No minimal-surface
classification is performed.

A small change on refining the sampling grid only verifies quadrature stability.
It does **not** verify the accuracy of the reconstructed physical dividing
surface or its ability to distinguish strain. Atomic sampling, thermal roughness,
hard assignment of boundary atoms, the smoothing length, and combining several
different interface types can dominate second derivatives. A matched geometric
control and a smoothing-scale study are required before interpreting minimality.

The exact catenoid and its affine image provide an independent mathematical
check; the labeled point-cloud control below additionally measures reconstruction
error. Its undeformed reference has exactly zero mean curvature, so any nonzero
reconstructed RMS is an estimator error floor for that specified control, not
evidence that the exact catenoid is nonminimal:

```bash
python -m example.crystalops.analysis.audit_curvature_reconstruction --output-dir analysis_results/self_audit/reproduction
```

This control fixes atom identities and labels before stretching; it does not
resegment the stretched cloud. Its error floor is not a universal error bar
for unrelated structures or different sigma/spacing/density settings.

For viewing exported values, use `valid=1`. Undefined atomic values are NaN,
not zero. VTK curvature is a triangle `CELL_DATA` property, not a vertex
`POINT_DATA` property.

## DXA defect mesh

`extract_dxa_defect_mesh` extracts the defect-region surface produced by OVITO
DXA directly from a complete atomistic structure. It returns the raw periodic
vertices/faces, a clipped display mesh without box-face caps, cell information,
parameters, and diagnostics. The default surface-smoothing iteration count is
zero. `crystal_structure` selects the reference atomic lattice: `fcc`, `hcp`,
`bcc`, `cubic_diamond`, or `hexagonal_diamond`; the diamond choices describe
atomic lattices.

Install extraction support with `pip install -e ".[defect-mesh]"`. From the
repository root:

```bash
python -m example.crystalops.analysis.extract_dxa_defect_mesh 1-1.data --crystal-structure fcc --length-unit angstrom --output-dir analysis_results/1-1_defect_mesh
```

The example writes an NPZ archive, a VTK mesh, and JSON metadata through
`crystalops.io.write_defect_mesh`. The NPZ can be loaded again with
`crystalops.io.read_defect_mesh`. Extraction uses the complete structure and
does not require grain labels or run the grain-boundary curvature pipeline.
`--length-unit` labels the original coordinate unit without conversion and
defaults to `input_length_unit`; the Cu example above explicitly uses angstroms.

For a static PNG and an offline HTML view that can be rotated and zoomed,
install `pip install -e ".[defect-mesh,mesh-viewer]"` and add `--plot`:

```bash
python -m example.crystalops.analysis.extract_dxa_defect_mesh 1-1.data --crystal-structure fcc --length-unit angstrom --output-dir analysis_results/1-1_defect_mesh --plot
```

`--show-atoms` adds the original atom positions to the views and implies
`--plot`. `--no-interactive` omits the HTML view. Display uses the clipped
defect mesh without periodic-boundary caps.

To extract and view the four supplied configurations together:

```bash
python -m example.crystalops.analysis.extract_dxa_defect_mesh_batch 1-1.6.data 1-1.data 1.2-1.data 1.4-1.data --crystal-structure fcc --surface-smoothing 0 --length-unit angstrom --output-dir analysis_results/dxa_defect_mesh
```

The batch example checks inputs and duplicate file stems before extraction,
saves each case under its filename stem, and creates `index.html` with a 2×2
comparison. Every card embeds the local offline HTML view and links its
standalone view and PNG; counts, cell-vector lengths, and smoothing settings
are shown. `--show-atoms` adds atoms to all views. The helper
`build_index(output_dir, records)` can build the page from existing artifacts
without repeating DXA.

## Curvature of a saved DXA defect mesh

`compute_defect_mesh_curvature` estimates `H = (k1 + k2) / 2` on the original
periodic DXA triangles using the cotangent Laplacian and the source mesh
normal. It uses the input lengths, with inverse length as the curvature unit.
Area statistics use triangles with three valid vertex estimates. The summary
includes the valid area fraction and counts for geometry or connectivity
issues. Invalid vertex values are NaN. The values describe the DXA defect
envelope at the smoothing setting recorded in the input archive.

```bash
python -m example.crystalops.analysis.analyze_defect_mesh_curvature analysis_results/dxa_defect_mesh/smoothing_2/1-1/defect_mesh.npz --output-dir analysis_results/dxa_defect_mesh/curvature/smoothing_2/1-1
```

`defect_mesh_curvature.npz` retains the native periodic geometry, per-vertex
curvature, face areas, and masks. `defect_mesh_curvature.json` records the
area-weighted signed mean, mean absolute curvature, and RMS curvature. For
visualization, color `defect_mesh_curvature.vtk` by its `mean_curvature`
cell field. The display mesh is cut at periodic boundaries with no cap faces.

Use `--fit-radius` to estimate curvature from a weighted quadratic fit over
each vertex's connected neighborhood. The radius is in input length units
and is measured along native mesh edges with periodic distances. Use a
common physical radius to compare configurations, and report fit radius and
supported area alongside the curvature. The estimate describes the fitted
geometry at that radius. Sharp folds and sizable fit residuals can limit its
interpretation as pointwise curvature. For the supplied Cu meshes:

```bash
python -m example.crystalops.analysis.analyze_defect_mesh_curvature analysis_results/dxa_defect_mesh/smoothing_2/1-1/defect_mesh.npz --fit-radius 8 --output-dir analysis_results/dxa_defect_mesh/curvature/local_fit_8A/smoothing_2/1-1
```

## Checking an exact Diamond template against boundary atoms

The example below selects two-grain contacts from atom coordinates and grain
IDs. It tests one exact Schwarz D period aligned with the actual rectangular
box. Registration fits translations on a subset of spatial blocks. Distances
on the held-out atoms and area-weighted distances in the reverse direction
assess agreement. With `--central-grain-id`, the example also tests contacts
between the remaining grains as a separate fit.

```bash
python -m example.crystalops.analysis.fit_diamond_to_grain_boundary_atoms 1-1.data --grain-labels analysis_results/new_data_tests/labels/1-1_grain_labels.csv --central-grain-id 1 --output-dir analysis_results/diamond_fit/1-1
```

The NPZ files retain atom IDs, the selected grain pair, coordinates, nearest
template points, and holdout masks. The JSON includes fit diagnostics and
the curvature of the exact template for those cell lengths. Template curvature
is a geometric prediction under this single-period model and its applicability
depends on the distance and coverage checks.
