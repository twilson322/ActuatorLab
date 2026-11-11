# Code Review: actuator-lab

Reviewed as a staff engineer on an actuator or joint-control team would review a pull request that adds a compliant-joint simulator, an FOC drive and seven position controllers. Scope: `include/actlab`, `src`, the two apps, the tests, the identification and analysis scripts, and the build.

## 1. Overall assessment

A strong, honest controls project. The plant, the drive and the controllers are separated cleanly, the controllers are scored with identical metrics and a shared tuning objective, robustness is measured with a 200-plant Monte Carlo and measured loop margins, and the README explains every failure.

The review found one real estimation defect, and it changed a reported result. The state-variable filter that feeds both adaptive controllers computed its acceleration output half a sample ahead of its other outputs. On an exact rigid-joint model this biased the inertia estimate by +6.5 % and invented a 3.1 N m s/rad viscous term. In the payload experiment the adaptive LQG's inertia estimate ended 6 % high, and its tracking was slightly better than it should have been. The fix puts the estimate within about 1 % of the truth. The README now reports the corrected numbers and says what changed.

The review also removed a duplicated copy of the whole adaptive front end (filter class and regression code existed twice), tied the Kalman filter's measurement noise to the configured encoders instead of hard-coded resolutions, and named the remaining unexplained constants. The deduplication was checked byte-for-byte against the committed results before the filter fix was applied.

## 2. Scorecard

| Area | Score | Basis |
|---|---|---|
| Correctness | 8 | FOC transforms, SVPWM, RK4 energy conservation, DARE residuals and the LQG's discretization are tested. One estimator defect fixed in this pass, with a test that pins it. |
| Architecture | 8 | Plant, drive, controllers and experiments are separate. Controllers share one interface and a factory. The adaptive front end is now one class. |
| Readability | 7 | Clear names and comments in the controllers. `experiments.cpp` is long (460 lines) with a string-keyed override chain. |
| Self-documentation | 8 | Units in parameter structs, reasons for clamps and thresholds, comments that record measured failures. |
| Maintainability | 7 | Tuned gains are copied from a grid-search CSV into `named_config` by hand. |
| Testing | 8 | 23 tests including physics invariants, RLS behavior, concurrency structures and closed-loop tracking. |
| Performance | 8 | 1 kHz host loop and 20 kHz current loop well within budget in simulation. Timing comes from a VM. |
| Reliability | 8 | Watchdog, anti-windup, speed clamp, covariance capping, plausibility bounds on adaptive estimates. |
| API design | 7 | `JointController` interface with `reset`, `update`, `commit` and `diag`. `diag(int)` is an untyped side channel. |
| Domain-specific correctness | 8 | Collocation and notch placement reasoned correctly, closed-loop identification uses an instrumental variable, ZOH discretization throughout, saturation modelled at every level. |
| Safety and robustness | 7 | Monte Carlo over plant uncertainty, measured margins, watchdog. No thermal or fault model. |

## 3. Problems

### Fixed in this review

| Severity | Location | What | Why | Fix and evidence |
|---|---|---|---|---|
| High | `SVF2::step` (now `include/actlab/output_rls.hpp`) | The acceleration output was $\dot x_2$ evaluated at the end of the sample, which includes the feedthrough $w_f^2 u_k$ of the newest input and leads $x_1$ and $x_2$ by about half a sample. | With an exact rigid-joint model the regression gave $J = 1.279$ for a true 1.2 and $b = 3.1$ N m s/rad for a true 0. In the payload experiment the adaptive LQG ended at 1.59 vs 1.5 kg m², and its slightly better tracking came from that bias. | Acceleration output is the mean over the sample, $(x_2[k] - x_2[k-1])/T$. Same test now gives $J$ within 1 % and $b$ within 0.1. Payload experiment rerun: 1.48 vs 1.5 kg m², tracking 0.30 / 0.48 mrad instead of 0.24 / 0.44. README updated with an explanation. |
| Medium | `src/controllers.cpp`, `src/lqg.cpp` | The filter class and the whole regression front end (four filters, gating, RLS setup with five literals) were duplicated. | The fix above would have had to be made twice, and the two copies could drift. | `OutputSideRls` with an `OutputRlsConfig` of named parameters. Byte-identical outputs on every track and adapt CSV before the filter fix. |
| Medium | `LqgController::synthesize` | Measurement noise used hard-coded 2^14 and 2^17 encoder resolutions. | A different encoder in `DriveParams` would have left the Kalman filter tuned for the old one. | `LqgConfig::encoder_motor_bits` and `encoder_link_bits`, defaulting to `DriveParams`. `Lqg.KalmanFilterFollowsTheEncoderResolution`. |
| Low | `named_config` | A missing identification file silently fell back to 13 and 29 Hz notches. | On another plant every notch would be in the wrong place with no message. | One warning on stderr (thread-safe, since the Monte Carlo workers call it). |
| Low | several | The 420 rad/s speed clamp, Coulomb smoothing widths, the LQG integrator limit, the estimator's plausibility bounds and the 5 % re-synthesis and notch-recentring thresholds were literals. | Tuning knobs without names. | Named constants or config fields with units and reasons. |

### Open

| Severity | Location | What | Why | Recommended fix |
|---|---|---|---|---|
| Medium | `exp_tune` | The grid search optimizes a time-domain objective with no margin constraint. | The README already shows the result: the tuned cascade sits at 22° phase margin and the Monte Carlo exposes it. | A phase-margin or $M_s$ constraint inside the tuner. |
| Low | `named_config` | Tuned gains are copied by hand from `results/tuning.csv`. | They drift from the tuner silently. | Generate the config from the tuning output, or check them against it in a test. |
| Low | `JointController::diag(int)` | Untyped diagnostic channel indexed by integers. | Callers must know which index means what for which controller. | A small struct of named optional diagnostics. |
| Low | adaptive cascade | Still unstable after the filter fix (diverges after the payload change). | A design limit, documented in the README: the rigid-body regressor is wrong while the 13 Hz mode is undamped. | Keep it as the documented negative result. |
| Low | plant | No thermal model, cogging, backlash or load-dependent gearbox stiffness. | Listed in the README's limitations. | As the README says. |

## 4. Readability and self-documentation

```cpp
// before (in two files)
svf_th_(10.0, 0.8, f_pos), svf_tau_(10.0, 0.8, f_pos), svf_fric_(10.0, 0.8, f_pos), svf_sin_(10.0, 0.8, f_pos),
rls_(RLS<4>::Vec(...), 5.0, 0.998, 1e4)
...
svf_fric_.step(std::tanh(svf_th_.dy() / 0.02));
if (std::abs(svf_th_.dy()) < 0.05) return;
...
const double lsb_m = 2.0 * kPi / 16384.0 / m_.N, lsb_l = 2.0 * kPi / 131072.0;

// after
rls_(OutputSideRls::Vec(...), f_pos)   // OutputRlsConfig names every one of those numbers
...
if (!rls_.update(th_l, m_.N * tau_applied_)) return;
...
const double lsb_m = 2.0 * kPi / std::pow(2.0, c_.encoder_motor_bits) / m_.N;
const double lsb_l = 2.0 * kPi / std::pow(2.0, c_.encoder_link_bits);
```

The comments on the notch placement, the friction-regressor collinearity and the speed clamp are good examples of explaining why.

## 5. Architecture and design

The split between plant (`Plant`, physics only), drive (`Drive`, FOC and current loop, watchdog) and controllers (`JointController`, position and velocity layers) mirrors real hardware: the controllers never see plant state, only quantized sensors. The `commit(tau_unsat, tau_sat)` call that feeds the saturated torque back to the controller is what makes back-calculation anti-windup and the RLS torque regressor correct. The RT harness reuses the same controllers on wall-clock threads with lock-free buffers.

## 6. Testing and verification

| Claim | Verified by |
|---|---|
| Plant conserves energy when lossless, holds static equilibrium, has the right locked-motor frequency | Math and tests |
| Clarke, Park and SVPWM are correct and saturate properly | Tests |
| Filters have the designed DC gain, notch depth and roll-off | Tests |
| DARE solutions satisfy the Riccati equation, ZOH matches the analytic double integrator | Tests |
| RLS converges, tracks a jump, and keeps its covariance bounded | Tests |
| The shared estimator identifies a rigid joint without bias | Test added in this review, fails on the old filter |
| Concurrency structures are race-free | Tests under ThreadSanitizer |
| Memory safety | All 23 tests under ASan and UBSan |
| Deduplication changed nothing | Byte-identical CSVs before the filter fix |
| Any of it on a real joint | Not verified |

## 7. Prioritized action plan

Must fix before commit: none outstanding.

Before sharing: margin constraint in the tuner, and generate or check the tuned gains automatically.

Recommended: typed diagnostics instead of `diag(int)`, and split `experiments.cpp` by experiment.

Optional: thermal model and load-dependent gearbox stiffness.

## 8. Final verdict

| Context | Verdict |
|---|---|
| Personal or portfolio project | STRONG |
| Shared codebase | ACCEPTABLE, after the tuner and gain-generation items |
| Simulation | STRONG |
| Hardware | REQUEST CHANGES: the controllers are structured for an STM32 port, but nothing has run on a real joint, and the tuned cascade's 22° margin needs fixing first. |
