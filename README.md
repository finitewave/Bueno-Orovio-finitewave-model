## Bueno-Orovio finitewave model
Two-dimensional implementation of the Bueno-Orovio–Cherry–Fenton (BOCF) model for simulating human ventricular tissue electrophysiology.

The BOCF model is a minimal phenomenological model developed to capture 
key ionic mechanisms and reproduce realistic human ventricular action potential 
dynamics, including restitution, conduction block, and spiral wave behavior. 
It consists of four variables: transmembrane potential (u), two gating variables (v, w), 
and one additional slow variable (s), representing calcium-related dynamics.

This model implementation can be used separately from the Finitewave, allowing for standalone simulations and testing of the model dynamics without the need for the entire framework.

### Reference
Bueno-Orovio, A., Cherry, E. M., & Fenton, F. H. (2008). 
Minimal model for human ventricular action potentials in tissue. J Theor Biol., 253(3), 544-60.

DOI: https://doi.org/10.1016/j.jtbi.2008.03.029

### How to use (quickstart)
```bash
python -m examples.bueno_orovio_example
```

### How to test
```bash
python -m pytest -q
```

### Repository structure
```text
.
├── bueno_orovio/                  # equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                     # model equations (pure functions)
├── implementation/                # 0D model implementation
│   ├── __init__.py
│   └── bueno_orovio_0d.py
├── example/
│   └── bueno_orovio_example.py    # minimal script to run a short trace
├── tests/
│   └── test.py                    # smoke test; reproducibility checks
├── .gitignore
├── LICENSE                        # MIT
├── pyproject.toml                 # configuration file
└── README.md                      # this file
```

### Variables
- `u = 0.0` — Membrane potential
- `v = 1.0` - Fast gating variable representing sodium channel inactivation.
- `w = 1.0` - Slow recovery variable representing calcium and potassium gating.
- `s = 0.0` - Slow variable related to calcium inactivation.

### Parameters
- `u_o       = 0.0`    - Resting membrane potential.
- `u_u       = 1.55`   - Peak potential (upper bound).
- `theta_v   = 0.3`    - Activation threshold for v.
- `theta_w   = 0.13`   - Activation threshold for w.
- `theta_v_m = 0.006`  - Threshold for switching time constants for v.
- `theta_o   = 0.006`  - Threshold for switching time constants for w.
- `tau_v1_m  = 60`     - Time constant for v below threshold.
- `tau_v2_m  = 1150`   - Time constant for v above threshold.
- `tau_v_p   = 1.4506` - Decay constant for v.
- `tau_w1_m  = 60`     - Base time constant for w.
- `tau_w2_m  = 15`    - Transition time constant for w.
- `k_w_m     = 65`     - Parameter controlling shape of τw curve.
- `u_w_m     = 0.03`   - Parameter controlling shape of τw curve.
- `tau_w_p   = 200`    - Decay constant for w above threshold.
- `tau_fi    = 0.11`   - Time constant for fast inward current (J_fi).
- `tau_o1    = 400`    - Time constant for outward current below threshold.
- `tau_o2    = 6`      - Time constant for outward current above threshold.
- `tau_so1   = 30.0181` - Time constant for repolarizing tail current below threshold.
- `tau_so2   = 0.9957` - Time constant for repolarizing tail current above threshold.
- `k_so      = 2.0458` - Parameter controlling nonlinearity in tau_so.
- `u_so      = 0.65`   - Parameter controlling nonlinearity in tau_so.
- `tau_s1    = 2.7342` - Time constant for s below threshold.
- `tau_s2    = 16`     - Time constant for s above threshold.
- `k_s       = 2.0994` - Parameter for tanh activation of s variable.
- `u_s       = 0.9087` - Parameter for tanh activation of s variable.
- `tau_si    = 1.8875` - Time constant for slow inward current (J_si).
- `tau_w_inf = 0.07`   - Slope of w∞ below threshold.
- `w_inf_    = 0.94`   - Asymptotic value of w∞ above threshold.


