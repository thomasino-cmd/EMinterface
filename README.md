# Code Review Report — Meniscus Optics Pipeline

## Scope
This review explains how the code models and computes optical propagation from an incident plane wave through an air–saline meniscus and onto a target plane (cuvette bottom proxy).

Primary files :
- `EMinterface/tile_overlap.cpp`
- `EMinterface/include/em/meniscus_model.hpp`
- `EMinterface/include/em/beam_mesh.hpp`
- `EMinterface/include/em/meniscus_projection.hpp`
- `EMinterface/include/em/frames.hpp`
- `EMinterface/include/em/tile_wrapper.hpp`
- `EMinterface/em_solver.cpp`
- `EMinterface/include/em/postprocess_grid.hpp`
- `EMinterface/include/geom/triangle.hpp`

---

## 1) Mesh Generation (Meniscus Discretization)

### 1.1 Geometric model of the meniscus
The meniscus is modeled as a paraboloid in global coordinates \((x_1,y_1,z_1)\):

\[
z_1 = f_2(x_1,y_1) = k(x_1^2+y_1^2)+a_{pm}
\]

with
\[
a_{pm} = \frac{V}{\pi R_{in}^2} - \frac{D_m}{2}, \qquad
k = \frac{D_m}{R_{in}^2}.
\]

This is implemented in `ParaboloidMeniscus`, via `apm_`, `k_`, and `f2(...)`.

### 1.2 Tessellation/discretization strategy
The surface is not meshed directly in curved 3D coordinates. Instead:
1. A **2D disk mesh** is generated in beam coordinates \((x_3,y_3,z_3=0)\).
2. Each mesh node is **projected along +z3** to intersect the paraboloid meniscus.
3. Triangles connect projected nodes, yielding a piecewise-planar approximation of the curved interface.

This workflow is done in `tile_overlap.cpp`:
- Mesh generation: `em::generate_disk_mesh(Rf, n_pr, n_rays)`
- Node projection: `em::project_node_to_meniscus(fr, men, x3, y3)`
- Triangle construction from projected nodes in the loop over `mesh.tris`.

### 1.3 Storage of vertices/faces
- Vertices in 2D generation stage: `Mesh2D::nodes` (`std::vector<geom::Vec2>`).
- Triangle connectivity: `Mesh2D::tris` (`std::vector<std::array<int,3>>`).
- Projected 3D interface vertices: `nodes_on_meniscus` (`std::vector<geom::Vec3>`), with validity mask `ok`.

---

## 2) Local Angle Calculation (Incident Wave vs. Each Triangle)

### 2.1 Triangle normal
For each interface triangle, the normal is computed by:
\[
\hat n = \frac{(\mathbf v_1-\mathbf v_0)\times(\mathbf v_2-\mathbf v_0)}{\left\|(\mathbf v_1-\mathbf v_0)\times(\mathbf v_2-\mathbf v_0)\right\|}
\]
implemented as `geom::normal_unit(...)` in `triangle.hpp`.

In `tile_wrapper.hpp`, `make_local_frame(...)` computes `ez = geom::normal_unit(tri, eps)`, then optionally flips it to align with incidence convention if `dot(s_in, ez) < 0`.

### 2.2 Local incidence angle
Given the incident unit direction `s_in_global_unit` and local normal `ez`, the code uses:
\[
\cos\theta_i = \hat s_{in}\cdot\hat n,
\qquad
\theta_i = \arccos(\cos\theta_i)
\]
implemented in `make_local_frame(...)` with `theta = std::acos(cos_th)`.

The local tangent basis (`ex`,`ey`) is then built from the tangential component of the incident direction.

---

## 3) Wave Vectors and Snell’s Law per Element

### 3.1 Snell in vector/wave-number form
In `EMInterfaceSolver::initialize_wave_vectors()` (`em_solver.cpp`), the transmitted wave vector is computed from tangential continuity:

- Keep tangential component:
\[
K_{p2x} = K_1\sin\theta_i
\]
- Enforce dispersion in medium 2:
\[
K_{p2z} = \sqrt{K_2^2 - K_{p2x}^2}
\]

with branch selection forcing physical decay (`Im(K_{p2z}) <= 0`), then setting `kp2 = (K_{p2x},0,K_{p2z})`.

Incident and reflected vectors are similarly built from sign-flipped z-components (`kp1`, `kr1`).

### 3.2 Reflected/transmitted field vectors
For each tile, `solve_on_triangle(...)` runs:
1. local frame construction,
2. local incident field vector (`make_incident_E_local(...)`),
3. interface solve via `EMInterfaceSolver` (TE/TM),
4. rotation back to global (`rotate_local_to_global(...)`).

TE and TM solutions are implemented in:
- `EMInterfaceSolver::solve_TE()`
- `EMInterfaceSolver::solve_TM()`

with boundary-condition algebra producing `Ar1` and `Ap2` vectors.

### 3.3 Fresnel coefficients?
There is **no explicit function named Fresnel coefficients** (e.g., no direct `r_s`, `t_s`, `r_p`, `t_p` function).
However, the solver computes equivalent physics by solving boundary-condition equations for field amplitudes (`Ar1`, `Ap2`) and then derives power coefficients:
\[
R = \frac{|S_{z,ref}|}{S_{z,inc}},\qquad
T = \frac{S_{z,tra}}{S_{z,inc}},\qquad
A = 1-R-T
\]
inside `calculate_power(...)`.

So Fresnel-like behavior is present implicitly through boundary-condition solve and Poynting-based flux ratios.

---

## 4) Beam Misalignment Handling

### 4.1 Parameters defining source axis/position
In `tile_overlap.cpp`, misalignment and pointing are parameterized by:
- `x_af`: lateral beam-axis offset,
- `theta`: tilt angle,
- `phi`: azimuthal rotation.

### 4.2 Coordinate transforms/offset usage
`Frames` (`frames.hpp`) defines transforms among coordinate systems and applies:
- translation by `x_af` and `z0` (`x2_from_x1`, `x1_from_x2`),
- rotations `Rz(phi)` and `Ry(theta)` (`R_x2_to_x3_`, `R_x3_to_x2_`).

In projection (`project_node_to_meniscus(...)`), the formulas for `x1(t), y1(t), z1(t)` explicitly include `x_af`, `theta`, and `phi`, so off-axis and tilted beams naturally shift where rays hit the meniscus.

---

## 5) Attenuation in Aqueous Solution

### 5.1 Propagation distance/path to target plane
For each triangle, vertices are projected from interface to the post-process plane `z1=Cpp` along energy transport direction `dir = res.d_poynting`.

For a point `p`:
\[
t = \frac{C_{pp}-p_z}{dir_z}, \qquad p' = p + t\,dir
\]
implemented in local lambda `proj_to_plane` in `tile_overlap.cpp`.

For attenuation, the code uses the centroid displacement:
\[
\Delta \mathbf r = \mathbf p_c' - \mathbf p_c
\]
where `c_face` is triangle centroid and `pc` is projected centroid.

### 5.2 Attenuation model
The code computes an attenuation vector from complex transmitted wave vector:
\[
\boldsymbol\alpha = -\operatorname{Im}(\mathbf k_{p2})
\]
and applies:
\[
\text{att} = \exp\big(-2\,\boldsymbol\alpha\cdot\Delta\mathbf r\big)
\]
then clamps to \([0,1]\).

This is a vectorized Beer–Lambert-like exponential decay using local complex propagation constant from each tile solution.

---

## 6) Power/Irradiance Mapping on Bottom Plane

### 6.1 Per-triangle power entering medium 2
For each projected meniscus triangle:
\[
P_{in} \approx \max\left(0,\langle\mathbf S\rangle\cdot\hat n\right)\,A_{face}
\]
where `res.S_avg` comes from `poynting_avg_real(...)`, `n_int` is interface normal oriented with incidence, and `A_face` from `geom::area(...)`.

### 6.2 Footprint mapping
Each triangle is projected to `z1=Cpp`, producing a 2D triangle (`tri2d`) in \((x_1,y_1)\). Its footprint area is `A_fp`.

After attenuation:
\[
P_{plane}=P_{in}\cdot \text{att},
\qquad
I = \frac{P_{plane}}{A_{fp}}
\]

Then `Grid2D::add_triangle_uniform(...)` adds this uniform intensity contribution to all cells whose centers lie inside the triangle.

### 6.3 Integration/summation of contributions
The final irradiance map is built by accumulation over all valid triangles. Power diagnostics are tracked by:
- `P_total_in += Pin`
- `P_total_on_plane += P_on_plane`

Finally, `grid.save_csv("field_map_cpp.csv")` outputs `(x,y,value)` where value is approximate intensity in `W/m^2`.

---

## 7) Practical Notes / Model Assumptions

- The meniscus surface is approximated with planar facets (triangle mesh resolution set by `n_pr`, `n_rays`).
- Energy footprint is treated as uniform over each projected triangle (no intra-triangle gradient).
- Interference between neighboring tiles is not explicitly superposed at field level; the map is power-additive.
- Attenuation is computed from local `Im(k)` and centroid path, which is a robust first-order approximation.

