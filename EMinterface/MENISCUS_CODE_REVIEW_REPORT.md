# Deep Code Review Report — Execution Trace of `tile_overlap.cpp`

## Goal of this report
This document gives a **strict execution-order walkthrough** of the simulation in `tile_overlap.cpp`, with emphasis on:
- where execution starts,
- which function is called next,
- which file owns that function,
- what physics/computation is performed there,
- and when control returns to `tile_overlap.cpp`.

It is written as a **function-by-function trace**, not only as a conceptual summary.

---

## 0) Entry point and architecture

### 0.1 Program entry
Execution starts in:
- `int main()` in `EMinterface/tile_overlap.cpp`.

At compile time, this file includes:
- `em_solver.cpp` (with `EMSOLVER_NO_MAIN`),
- geometry/math headers,
- meniscus/model/projection/frame utilities,
- tile-level EM wrapper,
- grid postprocessing.

So runtime control starts in `tile_overlap.cpp::main`, and from there jumps into functions defined in headers and in `em_solver.cpp`.

---

## 1) Step-by-step runtime execution order (exact flow)

Below is the practical sequence as the executable runs.

---

### STEP 1 — Read and set experiment/model parameters in `main`

**File:** `tile_overlap.cpp`  
**Function:** `main()`

The code defines all scenario inputs:
- cuvette geometry: `Rin`, `V`, `Dm`,
- beam geometry: `Rf`, misalignment `x_af`, tilts `theta`, `phi`,
- mesh resolution: `n_pr`, `n_rays`,
- target plane: `Cpp`, output grid `Nx`, `Ny`,
- EM media: `media.eps1_r`, `media.mu1_r`, `media.sigma1`, `media.eps2_r`, `media.mu2_r`, `media.sigma2`, frequency `omega`,
- optical excitation: `beam.A_inc`, `beam.isTE`, and beam direction sign `beam_sign`.

**Why first?** Because every downstream object (meniscus, transforms, solver, mapping) depends on these initial physical parameters.

---

### STEP 2 — Construct the meniscus model object

**Call from:** `tile_overlap.cpp::main`  
**To file/function:** `meniscus_model.hpp` → `ParaboloidMeniscus::ParaboloidMeniscus(...)`

Code path:
1. `em::ParaboloidMeniscus men(em::MeniscusParams{V, Dm, Rin});`
2. Constructor validates: `V>0`, `Rin>0`.
3. Computes and stores:
   - `apm_ = V/(pi Rin^2) - Dm/2`,
   - `k_ = Dm/(Rin^2)`.

Then `main` calls:
- `men.apm()` for checks,
- `men.f2(x_af,0.0)` to compute `z0`.

**Why this step here?** `z0` is needed to place the local frame origin consistently with the meniscus shape at beam offset `x_af`.

Control returns to `main` after each getter/evaluation.

---

### STEP 3 — Build coordinate transform system (`Frames`)

**Call from:** `tile_overlap.cpp::main`  
**To file/function:** `frames.hpp` → `Frames::Frames(const FrameParams&)`

Code path:
1. `em::Frames fr(em::FrameParams{x_af, theta, phi, z0});`
2. Inside constructor, two rotation matrices are precomputed:
   - forward: `R_x2_to_x3_ = Ry(theta) * Rz(phi)`,
   - inverse: `R_x3_to_x2_ = Rz(-phi) * Ry(-theta)`.

Then `main` calls:
- `fr.z3_axis_in_x1()`.

Inside `z3_axis_in_x1()`:
1. create `ez3 = (0,0,1)`,
2. transform it with inverse rotation into global-aligned frame,
3. normalize result.

Back in `main`, incident direction is set as:
- `s_in = normalize(z3_dir_x1 * beam_sign)`,
- `beam.s_in_global = s_in`.

**Why this step?** It translates your experimental beam-axis parameters (`x_af`, `theta`, `phi`) into the global direction used by every triangle-local solve.

---

### STEP 4 — Generate beam footprint mesh in (x3,y3)

**Call from:** `tile_overlap.cpp::main`  
**To file/function:** `beam_mesh.hpp` → `generate_disk_mesh(Rf,n_pr,n_rays)`

Inside `generate_disk_mesh`:
1. validates input,
2. inserts center node,
3. creates concentric rings of nodes,
4. creates triangle connectivity:
   - fan around center,
   - two triangles per annular quadrilateral.

Returns `Mesh2D mesh` with:
- `mesh.nodes` (2D points),
- `mesh.tris` (index triplets).

**Why this step?** It discretizes the beam aperture so each triangle can be treated as a local planar optical interaction patch.

---

### STEP 5 — Project every mesh node from x3-plane to real meniscus

**Loop in:** `tile_overlap.cpp::main` over `mesh.nodes`  
**Call to file/function:** `meniscus_projection.hpp` → `project_node_to_meniscus(fr, men, x3, y3)`

Inside `project_node_to_meniscus` for each node:
1. reads frame params (`x_af`, `theta`, `phi`, `z0`) via `fr.params()`,
2. writes parametric ray in global coordinates as function of `t = z3`:
   - `x1(t)`, `y1(t)`, `z1(t)`,
3. imposes intersection with meniscus equation `z1(t)=k(x1^2+y1^2)+apm`,
4. derives quadratic `A t^2 + B t + C = 0`,
5. calls `smallest_nonneg_root(A,B,C,t)` (same file) to pick first physical hit,
6. computes hit point `p1=(x1,y1,z1)`,
7. validates inside cuvette via `men.inside_cuvette(x1,y1,eps)`.

Returns `ProjectionResult{ok,t_z3,p1}`.

Back in `main`:
- if `ok`, store `nodes_on_meniscus[i]=p1`, mark `ok[i]=1`.

**Why this step?** This is the true geometry transfer from abstract beam disk to curved air-saline interface.

---

### STEP 6 — Initialize target irradiance grid

**File/function:** `postprocess_grid.hpp` → `Grid2D::Grid2D(...)`

`main` creates:
- `Grid2D grid(-Rin, Rin, -Rin, Rin, Nx, Ny)`.

Constructor validates bounds and allocates `val[Nx*Ny]` initialized to zero.

**Why now?** We are about to accumulate per-triangle intensity contributions on this map.

---

### STEP 7 — Triangle processing loop (core optical execution)

This is the central runtime block in `tile_overlap.cpp`: `for (const auto& tri_idx : mesh.tris)`.

For each triangle, execution chain is:

---

#### STEP 7.1 — Build 3D triangle on meniscus

In `main`:
1. fetch indices `i0,i1,i2`,
2. skip if any projected node invalid,
3. build `geom::Triangle tri` from `nodes_on_meniscus[...]`.

No cross-file call yet beyond struct construction.

---

#### STEP 7.2 — Solve local EM interface problem on that triangle

**Call from:** `tile_overlap.cpp::main`  
**To file/function:** `tile_wrapper.hpp` → `solve_on_triangle(tri, media, beam, true)`

This function is the orchestrator between geometry and EM solver.

Inside `solve_on_triangle`:

##### a) Build local frame per facet
- Calls `make_local_frame(tri, beam.s_in_global, ...)` (same file).
- Inside:
  - `geom::normal_unit(tri)` from `triangle.hpp`,
  - optional flip so normal matches incidence convention,
  - computes `theta_i = acos(dot(s_in,ez))`,
  - builds tangent basis `ex,ey`.

##### b) Build local incident field components
- Calls `make_incident_E_local(theta_i, A_inc, isTE)`.
- TE: field along local y.
- TM: field split on local x/z using `cos(theta_i), sin(theta_i)`.

##### c) Run EM boundary solve
- Instantiates `EMInterfaceSolver(...)` from `em_solver.cpp`.
- Constructor immediately calls `initialize_wave_vectors()`.

Inside `initialize_wave_vectors()` (`em_solver.cpp`):
1. from `theta_i` compute incident/reflected unit components,
2. compute medium scalar wavenumbers `K1`, `K2`,
3. enforce tangential conservation `Kp2x=K1*sin(theta_i)`,
4. compute `Kp2z = sqrt(K2^2-Kp2x^2)`, choose physical branch (`Im(Kp2z)<=0`),
5. assemble `kp1`, `kr1`, `kp2`.

Then `solve_on_triangle` calls `solver.solve()`.

Inside `EMInterfaceSolver::solve()`:
- dispatches to `solve_TE()` or `solve_TM()`.

Inside `solve_TE` / `solve_TM`:
- solve boundary-condition algebra for reflected and transmitted field amplitudes (`Ar1`, `Ap2`).

Then `calculate_power(...)` computes:
- Poynting flux components via `calc_Sz`,
- `R`, `T`, `A`, and diagnostic transmission angle.

Control returns to `solve_on_triangle` with local solution.

##### d) Rotate local outputs back to global
- `rotate_local_to_global(sol_local.Ap2, ex,ey,ez)`
- `rotate_local_to_global(kp2_local, ex,ey,ez)`

##### e) Build propagation/transport directions
- `d_phase ~ normalize(Re(kp2_global))`,
- `d_atten ~ normalize(-Im(kp2_global))`,
- computes average Poynting vector with `poynting_avg_real(...)`,
- `d_poynting ~ normalize(S_avg)`.

Returns `TileSolveResult` to `tile_overlap.cpp`.

**Why this big step?** Each triangle has different local normal → different local incidence angle → different local transmission/reflection/attenuation behavior.

---

#### STEP 7.3 — Compute power entering saline through triangle

Back in `tile_overlap.cpp`:
1. calls helper `tri_normal_pointing_with_incident(tri,s_in)` (same file), internally using `geom::normal_unit`.
2. computes face area `A_face = geom::area(tri)`.
3. computes incoming power estimate:
   \[
   P_{in}=\max(0,\,\mathbf S_{avg}\cdot\hat n_{int})\,A_{face}
   \]
4. skips if non-positive.

**Why?** Only positive forward power contributes to deposited energy.

---

#### STEP 7.4 — Project triangle footprint to target plane

In `main`, with direction `dir = res.d_poynting`:
1. if `|dir.z|` tiny, skip (nearly parallel to target plane),
2. local lambda `proj_to_plane` computes per-vertex line-plane hit to `z=Cpp`:
   \[
   t=(Cpp-p_z)/dir_z, \quad p'=p+t\,dir
   \]
3. get projected vertices `p0,p1,p2`, build 2D triangle `tri2d` in x-y.

Computes footprint area `A_fp` using shoelace expression; skip if too small.

**Why?** Converts transmitted 3D transport into where power lands on target plane.

---

#### STEP 7.5 — Compute in-medium attenuation to target

Still in `main`:
1. compute face centroid `c_face`, project centroid to plane (`pc`),
2. displacement `disp = pc.second - c_face`,
3. attenuation vector:
   - `alpha = -Im(res.kp2_global)`,
4. attenuation factor:
   \[
   att = \exp\big(-2\,\alpha\cdot disp\big)
   \]
5. clamp `[0,1]`.

Then:
- `P_on_plane = Pin * att`.

**Why?** Applies Beer–Lambert-like decay using local complex propagation constant and local propagation path.

---

#### STEP 7.6 — Deposit irradiance contribution on raster grid

In `main`:
1. uniform triangle intensity:
   \[
   I = P_{on\_plane}/A_{fp}
   \]
2. call `grid.add_triangle_uniform(tri2d, I)`.

Inside `Grid2D::add_triangle_uniform` (`postprocess_grid.hpp`):
1. compute triangle bounding box,
2. convert to pixel index bounds,
3. loop pixels in box,
4. compute pixel center,
5. test inclusion with `point_in_tri(...)`,
6. if inside: `grid.at(ix,iy)+=I`.

**Why?** This is the numerical integration/accumulation step producing final irradiance map.

---

### STEP 8 — Save and print diagnostics

Back in `tile_overlap.cpp::main` after triangle loop:
1. `grid.save_csv("field_map_cpp.csv")` (in `postprocess_grid.hpp`),
2. print summary: geometry, mesh stats, power totals.

Program exits.

---

## 2) Condensed call-order list (file-by-file jump map)

If we flatten execution into a call map:

1. `tile_overlap.cpp::main`
2. `meniscus_model.hpp::ParaboloidMeniscus::ParaboloidMeniscus`
3. `meniscus_model.hpp::ParaboloidMeniscus::apm` / `f2`
4. `frames.hpp::Frames::Frames`
5. `frames.hpp::Frames::z3_axis_in_x1`
6. `beam_mesh.hpp::generate_disk_mesh`
7. loop nodes:
   - `meniscus_projection.hpp::project_node_to_meniscus`
   - `meniscus_projection.hpp::smallest_nonneg_root`
   - `meniscus_model.hpp::ParaboloidMeniscus::inside_cuvette`
8. `postprocess_grid.hpp::Grid2D::Grid2D`
9. loop triangles:
   - `tile_wrapper.hpp::solve_on_triangle`
     - `tile_wrapper.hpp::make_local_frame`
       - `triangle.hpp::normal_unit`
     - `tile_wrapper.hpp::make_incident_E_local`
     - `em_solver.cpp::EMInterfaceSolver::EMInterfaceSolver`
       - `em_solver.cpp::initialize_wave_vectors`
     - `em_solver.cpp::EMInterfaceSolver::solve`
       - `em_solver.cpp::solve_TE` or `solve_TM`
       - `em_solver.cpp::calculate_power`
     - `tile_wrapper.hpp::rotate_local_to_global`
     - `tile_wrapper.hpp::poynting_avg_real`
   - `triangle.hpp::area` / `normal_unit`
   - `postprocess_grid.hpp::Grid2D::add_triangle_uniform`
     - `postprocess_grid.hpp::Grid2D::point_in_tri`
10. `postprocess_grid.hpp::Grid2D::save_csv`
11. return from `main`

---

## 3) Why this execution order is physically consistent

- **Geometry first** (meniscus + frames + projection) is mandatory before optics because local incidence depends on actual curved facet orientation.
- **Local EM solve per facet** is needed because each triangle has unique `theta_i` and local basis.
- **Transport to target + attenuation** happens after interface transmission, not before.
- **Raster accumulation last** converts many local contributions into the final map.

This ordering mirrors the physical pipeline:
\[
\text{Experiment parameters} \rightarrow \text{Interface geometry} \rightarrow \text{Local refraction/reflection} \rightarrow \text{Propagation in saline} \rightarrow \text{Power map on cells}
\]

---

## 4) Important implementation notes for your specific question

1. **Yes, after parameters the first major cross-file jump is to `meniscus_model.hpp`** when constructing `ParaboloidMeniscus`.  
2. **Then it jumps to `frames.hpp`**, where rotation matrices and beam-axis direction in global coordinates are prepared.  
3. **Then to `beam_mesh.hpp`** for tessellation.  
4. **Then repeatedly to `meniscus_projection.hpp`** to map mesh nodes onto the real curved interface.  
5. **Then repeatedly to `tile_wrapper.hpp` + `em_solver.cpp`** for per-triangle EM interface solve.  
6. **Then to `postprocess_grid.hpp`** for map accumulation and CSV export.

So your intuition about execution progression is correct; this report refines it into exact function-level transitions.
