# Final Report — Per-triangle outputs, medium-2 field quantities, and CSV interpretation

## Scope
This report focuses exclusively on **outputs** of the simulation pipeline:
- what is computed for each triangle,
- how the transmitted electric field in medium 2 is obtained,
- how amplitude/power per triangle is used,
- how contributions are accumulated,
- how to interpret the final CSV (`field_map_cpp.csv`).

Primary files involved:
- `EMinterface/tile_overlap.cpp`
- `EMinterface/include/em/tile_wrapper.hpp`
- `EMinterface/em_solver.cpp`
- `EMinterface/include/em/postprocess_grid.hpp`
- `EMinterface/include/geom/triangle.hpp`

---

## 1) What is the per-triangle computational output?

For each valid meniscus triangle, `tile_overlap.cpp` calls:

```cpp
em::TileSolveResult res = em::solve_on_triangle(tri, media, beam, true);
```

So the first “output block” per triangle is `TileSolveResult`.

From `tile_wrapper.hpp`, `TileSolveResult` contains:

1. **Complex transmitted electric field vector in global coordinates**
   - `Ap2_global` (phasor vector, medium 2).
2. **Complex transmitted wavevector in global coordinates**
   - `kp2_global`.
3. **Derived directions**
   - `d_phase` (from `Re(k)`),
   - `d_atten` (from `-Im(k)`),
   - `d_poynting` (from real time-average Poynting vector).
4. **Real time-average Poynting vector**
   - `S_avg` in `W/m^2`.
5. **Interface coefficients/diagnostics**
   - `theta_i`, `R`, `T`, `A`.

These are the local optical outputs before geometric transport/deposition.

---

## 2) How transmitted electric field in medium 2 is calculated

## Step A — Local frame + local incident field

In `solve_on_triangle(...)`:
1. build local frame (`ex,ey,ez`, `theta_i`) from triangle normal + incident direction,
2. build local incident electric field `Ap1_local`:
   - TE: `(0, A_inc, 0)`,
   - TM: `(A_inc cos(theta_i), 0, -A_inc sin(theta_i))`.

So, the incident wave amplitude is injected as complex phasor `A_inc` in the local polarization basis.

---

## Step B — Interface solve in `EMInterfaceSolver`

`EMInterfaceSolver` computes medium parameters and wavevectors, then solves boundary conditions.

### Medium quantities (`Medium` constructor)
For each medium:

\[
\varepsilon_{eff} = \varepsilon_r\varepsilon_0 - j\frac{\sigma}{\omega}
\]
\[
k = \omega\sqrt{\mu\varepsilon_{eff}}, \qquad \eta = \sqrt{\mu/\varepsilon_{eff}}
\]

Thus medium losses enter both through complex permittivity and complex wave number.

### Wavevector transmitted component
In `initialize_wave_vectors()`:

\[
K_{p2x}=K_1\sin\theta_i,
\qquad
K_{p2z}=\sqrt{K_2^2-K_{p2x}^2}
\]

with branch selection ensuring physical decay convention (`Im(Kp2z)<=0`).

### TE/TM boundary solve
- `solve_TE()` computes reflected/transmitted field components (`Ar1`, `Ap2`) from TE equations.
- `solve_TM()` computes analogous TM components.

So the **local transmitted electric field phasor** is `sol_local.Ap2`.

---

## Step C — Rotate local transmitted field to global

Back in `solve_on_triangle`:

```cpp
out.Ap2_global = rotate_local_to_global(sol_local.Ap2, fr.ex, fr.ey, fr.ez);
```

This is the per-triangle medium-2 electric field output in global coordinates.

### Interpretation
`Ap2_global` is the complex field amplitude vector at the interface facet (local phasor amplitude after transmission), not yet attenuated to the bottom plane.

---

## 3) What is done with that medium-2 amplitude?

Directly, `tile_overlap.cpp` does not sum complex fields from different triangles.
Instead, it uses `Ap2_global` and `kp2_global` to compute power-flow quantities:

- `S_avg` via `poynting_avg_real(E,k,...)`:
  \[
  \mathbf H = \frac{1}{\mu\omega}(\mathbf k\times \mathbf E),
  \qquad
  \langle\mathbf S\rangle = \frac{1}{2}\Re\left(\mathbf E\times\mathbf H^*\right)
  \]

So the transmitted amplitude influences output through `S_avg` (and directions derived from `k`).

No coherent interference term between triangles is accumulated at electric-field level.
The final map is a **power accumulation** map.

---

## 4) Per-triangle power chain to bottom plane

For each triangle, `tile_overlap.cpp` executes this chain:

## 4.1 Interface-crossing power

\[
P_{in} = \max\big(0,\langle\mathbf S\rangle\cdot\hat n_{int}\big)\,A_{face}
\]

Where:
- `A_face = geom::area(tri)`,
- `n_int` is interface normal oriented with incident direction.

Interpretation: watts entering medium 2 through that facet.

---

## 4.2 Geometric transport to target plane

Direction used: `dir = res.d_poynting`.

Project each triangle vertex (and centroid) to `z1=Cpp`:

\[
t = \frac{Cpp-p_z}{dir_z}, \qquad p' = p + t\,dir
\]

Footprint area on target plane: `A_fp`.

---

## 4.3 Lossy-medium attenuation to bottom plane

From complex transmitted wavevector:

\[
\boldsymbol\alpha = -\Im(\mathbf k_{p2,global})
\]

Centroid path:
\[
\Delta\mathbf r = \mathbf p_c' - \mathbf p_c
\]

Attenuation:

\[
att = \exp\left(-2\,\boldsymbol\alpha\cdot\Delta\mathbf r\right)
\]

Then:

\[
P_{on\_plane}=P_{in}\cdot att
\]

This is the actual per-triangle delivered power at the target plane in this model.

---

## 4.4 Convert delivered power to irradiance contribution

Each projected triangle contributes uniform irradiance:

\[
I = \frac{P_{on\_plane}}{A_{fp}}
\]

This `I` is rasterized by `grid.add_triangle_uniform(...)` and added to all pixels whose centers lie inside the footprint triangle.

---

## 5) Quantities that are accumulated globally

Inside loop:
- `P_total_in += Pin`
- `P_total_on_plane += P_on_plane`
- `grid` values receive `+I` on covered cells.

So there are two scalar global power diagnostics + one spatial map.

Interpretation:
- `P_total_in`: estimated total power entering medium 2.
- `P_total_on_plane`: estimated total power reaching plane `z1=Cpp` after attenuation.
- `grid(x,y)`: estimated local irradiance map at plane.

---

## 6) What exactly is written in the CSV

`grid.save_csv("field_map_cpp.csv")` outputs rows:

```text
x,y,value
```

for each grid cell center `(x,y)`.

- `x`: x-coordinate of cell center (meters).
- `y`: y-coordinate of cell center (meters).
- `value`: accumulated irradiance-like quantity from all triangles covering that cell.

In current implementation, `value` is in practice interpreted as **intensity / irradiance in W/m²** (as also printed in console message).

So each CSV row is one map pixel sample of bottom-plane irradiance.

---

## 7) How to read the CSV physically

Given one row `(x_i, y_j, value_ij)`:
- `value_ij` is the predicted irradiance at that bottom-plane location.
- high value = more delivered optical power density.
- zero/low value = little or no footprint overlap and/or strong attenuation.

To estimate total power from CSV approximately:

\[
P_{csv} \approx \sum_{cells} value_{ij}\,\Delta x\,\Delta y
\]

with `Δx = grid.dx()`, `Δy = grid.dy()`.

This should be broadly consistent with `P_total_on_plane` (differences can exist due to raster discretization and uniform-per-triangle deposition approximation).

---

## 8) Trace of output quantities from birth to file

End-to-end per triangle:

1. `A_inc` + media + local geometry → `Ap2_local`, `kp2_local`.
2. rotate → `Ap2_global`, `kp2_global`.
3. `Ap2_global`,`kp2_global` → `S_avg` + directions.
4. `S_avg`, normal, area → `Pin`.
5. projection geometry + `kp2_global` imaginary part → `att`.
6. `Pin * att` → `P_on_plane`.
7. `P_on_plane / A_fp` → `I`.
8. `I` accumulated in grid cells.
9. grid dumped to CSV as `x,y,value`.

This is exactly the output construction path used by the code.

---

## 9) Important output-model assumptions

1. **Power-domain accumulation** (not coherent E-field summation between triangles).
2. **Uniform irradiance per projected triangle footprint**.
3. **Single centroid attenuation path per triangle**.
4. **Transport direction based on Poynting direction (`d_poynting`)**.

These assumptions define how to interpret CSV values: they are engineering irradiance estimates from facet-wise energy transport with lossy attenuation.

---

## 10) Bottom-line interpretation

- The per-triangle transmitted field amplitude in medium 2 is computed as complex vector `Ap2_global`.
- That amplitude is converted to power flow through `S_avg`, then to facet power `Pin`.
- Loss in medium 2 is applied exponentially using `-Im(kp2_global)` over propagation path.
- Final CSV `value` is the accumulated irradiance (`W/m²`) on the bottom/postprocess plane.
