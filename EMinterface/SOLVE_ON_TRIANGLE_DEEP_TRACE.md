# Report tecnico ultra-dettagliato: esecuzione di `solve_on_triangle(...)`

## Obiettivo
Questo report descrive **in modo molto dettagliato** l'esecuzione della funzione:

- `em::solve_on_triangle(...)` in `EMinterface/include/em/tile_wrapper.hpp`

con:
- ordine esatto delle istruzioni,
- chiamate annidate (e file di destinazione),
- ritorni di controllo,
- motivazione fisica/numerica di ogni passaggio.

---

## 1) Firma della funzione e significato dei parametri

Funzione analizzata:

```cpp
inline TileSolveResult solve_on_triangle(const geom::Triangle& tri,
                                         const MediaParams& media,
                                         const BeamParams& beam,
                                         bool flip_normal_to_match_incident = true,
                                         double eps = 1e-12)
```

### Parametri in ingresso
- `tri`: triangolo locale sulla superficie menisco (3 vertici globali 3D).
- `media`: parametri elettromagnetici dei due mezzi (aria e soluzione), inclusa `omega`.
- `beam`: direzione incidente globale `s_in_global`, ampiezza `A_inc`, polarizzazione `isTE`.
- `flip_normal_to_match_incident`: forza la normale del triangolo verso il verso d'incidenza.
- `eps`: tolleranza numerica per normalizzazione e controlli di robustezza.

### Output
- `TileSolveResult`: campi trasmessi globali, vettore d'onda trasmesso globale, direzioni utili (`d_phase`, `d_atten`, `d_poynting`), Poynting medio, e diagnostiche (`theta_i`, `R`, `T`, `A`).

---

## 2) Esecuzione istruzione-per-istruzione dentro `solve_on_triangle(...)`

Di seguito la funzione viene “srotolata” nell'ordine reale di esecuzione.

---

### STEP A — Validazione frequenza

```cpp
if (media.omega <= 0.0) throw std::runtime_error("solve_on_triangle: omega must be > 0");
```

**Cosa fa:** blocca casi non fisici/non definiti.

**Perché:** tutte le formule EM usate dopo dipendono da `omega` (ad esempio `H = (1/(mu*omega)) (k x E)`). Se `omega<=0`, si creano divisioni improprie o modello non fisico.

---

### STEP B — Costruzione frame locale del triangolo

```cpp
LocalFrame fr = make_local_frame(tri, beam.s_in_global, flip_normal_to_match_incident, eps);
```

Questa è la prima chiamata importante. Entra in:
- file: `tile_wrapper.hpp`
- funzione: `make_local_frame(...)`

#### B.1 Normalizzazione direzione incidente
```cpp
s_in_global_unit = geom::normalize(s_in_global_unit, eps);
```
**Perché:** il prodotto scalare con la normale deve restituire direttamente `cos(theta_i)`.

#### B.2 Normale geometrica triangolo
```cpp
geom::Vec3 ez = geom::normal_unit(tri, eps);
```
Chiamata a:
- file: `triangle.hpp`
- funzione: `normal_unit(...)`
  - internamente usa `normal_unnormalized = cross(v1-v0, v2-v0)`
  - poi normalizza.

**Perché:** il triangolo discretizzato approssima una patch locale planare; serve la normale locale per angolo d'incidenza e piano d'incidenza.

#### B.3 Eventuale flip della normale
```cpp
if (flip_normal_to_match_incident && geom::dot(s_in_global_unit, ez) < 0.0) {
    ez = ez * -1.0;
}
```
**Perché:** convenzione coerente: normale orientata nel semispazio “visto” dal raggio incidente, evitando angoli ottusi non desiderati e segni incoerenti.

#### B.4 Angolo locale di incidenza
```cpp
const double cos_th = clamp01(geom::dot(s_in_global_unit, ez));
const double theta = std::acos(cos_th);
```
**Perché:** `theta_i` è la variabile fisica che alimenta il solver di interfaccia locale.

#### B.5 Asse tangenziale `ex` (piano d'incidenza)
```cpp
geom::Vec3 ex = s_in_global_unit - ez * cos_th;
```
- rimuove la componente normale dalla direzione incidente.

Se quasi incidenza normale:
```cpp
if (geom::is_zero(ex, eps)) {
    geom::Vec3 a = (std::abs(ez.z) < 0.9) ? geom::Vec3{0.0,0.0,1.0} : geom::Vec3{0.0,1.0,0.0};
    ex = geom::cross(a, ez);
}
```
**Perché:** a incidenza normale il piano d'incidenza non è univoco; si sceglie una tangente stabile numericamente.

Poi:
```cpp
ex = geom::normalize(ex, eps);
geom::Vec3 ey = geom::cross(ez, ex);
ey = geom::normalize(ey, eps);
```
**Perché:** costruisce terna ortonormale locale destra (`ex`,`ey`,`ez`).

#### B.6 Return
`make_local_frame` ritorna `LocalFrame{ex, ey, ez, theta}` a `solve_on_triangle`.

---

### STEP C — Costruzione campo incidente locale

```cpp
Vector3D Ap1_local = make_incident_E_local(fr.theta_i, beam.A_inc, beam.isTE);
```
Chiamata a:
- file: `tile_wrapper.hpp`
- funzione: `make_incident_E_local(...)`

#### C.1 Caso TE
```cpp
if (isTE) return Vector3D(0.0, A_inc, 0.0);
```
**Perché:** nel frame locale TE, il campo elettrico è ortogonale al piano di incidenza (asse `y`).

#### C.2 Caso TM
```cpp
const double c = std::cos(theta_i);
const double s = std::sin(theta_i);
return Vector3D(A_inc * c, 0.0, A_inc * (-s));
```
**Perché:** in TM il campo è nel piano di incidenza (`x-z`) con decomposizione coerente alla geometria locale.

Return a `solve_on_triangle` con `Ap1_local`.

---

### STEP D — Istanza solver elettromagnetico locale

```cpp
EMInterfaceSolver solver(..., media.omega, fr.theta_i, Ap1_local, beam.isTE);
```
Questa chiamata entra in:
- file: `em_solver.cpp`
- costruttore: `EMInterfaceSolver::EMInterfaceSolver(...)`

#### D.1 Nel costruttore: inizializzazione mezzi
Oggetti `Medium m1,m2` calcolano:
- `eps_eff = eps_r*eps0 - j sigma/omega`
- `k = omega*sqrt(mu*eps_eff)`
- `eta = sqrt(mu/eps_eff)`

**Perché:** servono numeri d'onda complessi e impedenze intrinseche per mezzi dissipativi.

#### D.2 Nel costruttore: chiamata immediata a `initialize_wave_vectors()`
Sempre in `em_solver.cpp`.

Passi interni:
1. `sin(theta_i)`, `cos(theta_i)`.
2. versori onda incidente/riflessa (`kp1v*`, `kr1v*`).
3. tangenziale trasmessa: `Kp2x = K1 * sin(theta_i)`.
4. componente normale trasmessa: `Kp2z = sqrt(K2^2 - Kp2x^2)`.
5. scelta ramo fisico: se `Im(Kp2z)>0` inverte segno (`Im<=0` per decadimento nel mezzo 2 con convenzione adottata).
6. compone vettori completi `kp1`, `kr1`, `kp2`.

**Perché fisico:** è la forma vettoriale della legge di Snell + condizione di radiazione/decadimento nel mezzo lossy.

Return al costruttore, poi ritorno a `solve_on_triangle` con oggetto `solver` pronto.

---

### STEP E — Risoluzione interfaccia (TE/TM)

```cpp
Solution sol_local = solver.solve();
```
Entra in:
- file: `em_solver.cpp`
- funzione: `EMInterfaceSolver::solve()`

#### E.1 Dispatch polarizzazione
```cpp
if (isTE) return solve_TE();
else return solve_TM();
```

---

#### E.2 Ramo TE (`solve_TE()`)
1. prende `Ap1y`.
2. usa `K1`, `mu1`, `mu2`, `kp1vz`, `kr1vz`, `kp2.z`.
3. calcola ampiezza riflessa `Ar1y` da condizioni al contorno.
4. calcola trasmessa `Ap2y = Ap1y + Ar1y`.
5. costruisce `Solution` con vettori campo riflesso/trasmesso.
6. chiama `calculate_power(sol)`.

---

#### E.3 Ramo TM (`solve_TM()`)
1. legge `Ap1x`, `Ap1z`.
2. prepara variabili (`K1`,`mu1`,`mu2`,`kp2x`,`kp2z`), quadrati.
3. costruisce denominatore/numeratore (`Denom`,`Num_Ar1x`).
4. calcola `Ar1x = Num_Ar1x / Denom`.
5. ricava `Ar1z`, `Ap2x`, `Ap2z`.
6. check singolarità `abs(Denom)<1e-15`.
7. costruisce `Solution`.
8. chiama `calculate_power(sol)`.

**Perché:** implementa condizioni di continuità tangenziale `E/H` per polarizzazione TM.

---

#### E.4 Calcolo potenze e coefficienti (`calculate_power`)
Internamente:
1. funzione locale `calc_Sz(E,k,m)`:
   - `H = (1/(mu*omega)) (k x E)`
   - `S = 0.5 E x conj(H)`
   - ritorna `Re(S.z)`.
2. calcola `Sz_inc`, `Sz_ref`, `Sz_tra`.
3. calcola:
   - `R = |Sz_ref|/Sz_inc`
   - `T = Sz_tra/Sz_inc`
   - `A = 1 - R - T`.
4. angolo trasmesso diagnostico da `atan2(Re(kp2.x), Re(kp2.z))`.

Return a `solve_TE`/`solve_TM`, poi ritorno da `solve()` a `solve_on_triangle`.

---

### STEP F — Lettura vettore d'onda trasmesso locale

```cpp
const Vector3D kp2_local = solver.get_kp2();
```

**Perché:** `Solution` non contiene sempre tutto ciò che serve per direzioni fase/attenuazione; qui serve esplicitamente `k` trasmesso complesso.

---

### STEP G — Rotazione risultati dal frame locale al globale

```cpp
out.Ap2_global = rotate_local_to_global(sol_local.Ap2, fr.ex, fr.ey, fr.ez);
out.kp2_global = rotate_local_to_global(kp2_local, fr.ex, fr.ey, fr.ez);
```

Chiamata a:
- file: `tile_wrapper.hpp`
- funzione: `rotate_local_to_global(...)`

Meccanica:
- combina componenti locali complesse `vx,vy,vz` con basi reali globali `ex,ey,ez`.

**Perché:** il solver EM lavora in frame locale di patch; la pipeline globale (`tile_overlap.cpp`) richiede vettori in coordinate globali per proiezione e mappa finale.

---

### STEP H — Copia diagnostiche locali

```cpp
out.theta_i = fr.theta_i;
out.R = sol_local.R;
out.T = sol_local.T;
out.A = sol_local.A;
```

**Perché:** queste metriche servono per debug/analisi fisica per triangolo.

---

### STEP I — Direzioni di fase e attenuazione

```cpp
const geom::Vec3 k_re = vec_re(out.kp2_global);
const geom::Vec3 k_im = vec_im(out.kp2_global);

out.d_phase = geom::safe_normalize(k_re, eps);
out.d_atten = geom::safe_normalize(k_im * (-1.0), eps);
```

Chiamate di supporto in `tile_wrapper.hpp`:
- `vec_re(...)` estrae parte reale di `k`,
- `vec_im(...)` estrae parte immaginaria,
- `safe_normalize(...)` evita instabilità su norme piccole.

**Perché fisico:**
- direzione fase ~ `Re(k)`;
- direzione decadimento ampiezza ~ `-Im(k)`.

---

### STEP J — Poynting medio reale e direzione trasporto energia

```cpp
const geom::Vec3 S = poynting_avg_real(out.Ap2_global, out.kp2_global, media.mu2_r, media.omega);
out.d_poynting = geom::safe_normalize(S, eps);
out.S_avg = S;
```

Chiamata a:
- file: `tile_wrapper.hpp`
- funzione: `poynting_avg_real(...)`

Dentro `poynting_avg_real`:
1. check `omega>0`.
2. `mu_abs = mu_r * MU0`.
3. `H = (1/(mu_abs*omega)) * (k x E)`.
4. `S = 0.5 * (E x conj(H))`.
5. ritorna parte reale.

**Perché:** `tile_overlap.cpp` deposita energia sul piano seguendo `d_poynting`, quindi questa è la direzione fisicamente più appropriata del trasporto di potenza media.

---

### STEP K — Return finale

```cpp
return out;
```

`solve_on_triangle` ritorna `TileSolveResult` al chiamante (`tile_overlap.cpp`), che prosegue con:
- stima `Pin` sulla faccia,
- proiezione footprint su `z=Cpp`,
- attenuazione,
- deposito su griglia.

---

## 3) Catena completa delle chiamate (call graph operativo)

Per una singola chiamata di `solve_on_triangle(...)`, la catena tipica è:

1. `tile_wrapper.hpp::solve_on_triangle`
2. `tile_wrapper.hpp::make_local_frame`
   - `triangle.hpp::normal_unit`
     - `triangle.hpp::normal_unnormalized`
     - `vec3.hpp::cross`, `vec3.hpp::normalize`
   - `vec3.hpp::dot`, `vec3.hpp::cross`, `vec3.hpp::normalize`
3. `tile_wrapper.hpp::make_incident_E_local`
4. `em_solver.cpp::EMInterfaceSolver::EMInterfaceSolver`
   - `em_solver.cpp::Medium::Medium` (m1,m2)
   - `em_solver.cpp::initialize_wave_vectors`
5. `em_solver.cpp::EMInterfaceSolver::solve`
   - `em_solver.cpp::solve_TE` **oppure** `em_solver.cpp::solve_TM`
   - `em_solver.cpp::calculate_power`
6. `em_solver.cpp::EMInterfaceSolver::get_kp2`
7. `tile_wrapper.hpp::rotate_local_to_global` (due volte)
8. `tile_wrapper.hpp::vec_re` / `vec_im`
9. `tile_wrapper.hpp::poynting_avg_real`
10. return `TileSolveResult`

---

## 4) Perché `solve_on_triangle(...)` è strutturata così

La struttura è intenzionale e riflette una separazione netta dei compiti:

1. **Geometria locale** (`make_local_frame`)  
   porta il triangolo curvo-discretizzato a un problema di piano locale.

2. **Fisica d'interfaccia locale** (`EMInterfaceSolver`)  
   applica condizioni al contorno TE/TM e mezzi lossy usando `theta_i` locale.

3. **Ritorno al mondo globale** (`rotate_local_to_global`)  
   rende i risultati utilizzabili dal resto della pipeline (proiezione su piano, mappa potenza).

4. **Direzioni fisiche operative** (`d_phase`,`d_atten`,`d_poynting`)  
   fornisce tre direzioni diverse utili per analisi e postprocessing.

In breve: `solve_on_triangle` è il “ponte” tra discretizzazione geometrica del menisco e trasporto energetico globale verso il piano cellule.

---

## 5) Nota pratica di verifica (debug consigliato)

Se vuoi verificare numericamente un triangolo specifico, i valori più utili da stampare nel chiamante sono:
- `fr.theta_i`,
- `out.R, out.T, out.A`,
- `out.kp2_global` (parti reali/immaginarie),
- `out.d_poynting` e `out.S_avg`.

Questa cinquina ti conferma:
- correttezza geometrica locale,
- correttezza energetica interfaccia,
- coerenza direzione di trasporto.
