"""
ops.py — mathematical core of the model.

This module provides functions to compute the model equations,
as well as functions to retrieve default parameters and initial
values for the state variables.

The BOCF (Bueno-Orovio-Cherry-Fenton) model is a minimal phenomenological model developed 
to capture key ionic mechanisms and reproduce realistic human ventricular action potential 
dynamics, including restitution, conduction block, and spiral wave behavior. 
It consists of four variables: transmembrane potential (u), two gating variables (v, w), 
and one additional slow variable (s), representing calcium-related dynamics.

This implementation corresponds to the EPI (epicardial) parameter set described in the paper.

References:
- Bueno-Orovio, A., Cherry, E. M., & Fenton, F. H. (2008). 
Minimal model for human ventricular action potentials in tissue. J Theor Biol., 253(3), 544-60.

DOI: https://doi.org/10.1016/j.jtbi.2008.03.029
"""

__all__ = (
    "get_variables",
    "get_parameters",
    "calc_rhs",  
    "calc_v",
    "calc_w",
    "calc_s",
    "calc_Jfi",
    "calc_Jso",
    "calc_Jsi",
    "calc_tau_v_m",
    "calc_tau_w_m",
    "calc_tau_so",
    "calc_tau_s",
    "calc_tau_o",
    "calc_v_inf",
    "calc_w_inf"
)

import math


def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.

    State Variables:
    u = 0.0 - Membrane potential
    v = 1.0 - Fast gating variable representing sodium channel inactivation.
    w = 1.0 - Slow recovery variable representing calcium and potassium gating.
    s = 0.0 - Slow variable related to calcium inactivation.
    """
    return {
        "u": 0.0, 
        "v": 1.0, 
        "w": 1.0, 
        "s": 0.0
    }


def get_parameters() -> dict[str, float]:
    """
    Returns default parameter values for the model.

    Parameters:
    u_o       = 0.0    - Resting membrane potential.
    u_u       = 1.55   - Peak potential (upper bound).
    theta_v   = 0.3    - Activation threshold for v.
    theta_w   = 0.13   - Activation threshold for w.
    theta_v_m = 0.006  - Threshold for switching time constants for v.
    theta_o   = 0.006  - Threshold for switching time constants for w.
    tau_v1_m  = 60     - Time constant for v below threshold.
    tau_v2_m  = 1150   - Time constant for v above threshold.
    tau_v_p   = 1.4506 - Decay constant for v.
    tau_w1_m  = 60     - Base time constant for w.
    tau_w2_m  = 15     - Transition time constant for w.
    k_w_m     = 65     - Parameter controlling shape of τw curve.
    u_w_m     = 0.03   - Parameter controlling shape of τw curve.
    tau_w_p   = 200    - Decay constant for w above threshold.
    tau_fi    = 0.11   - Time constant for fast inward current (J_fi).
    tau_o1    = 400    - Time constant for outward current below threshold.
    tau_o2    = 6      - Time constant for outward current above threshold.
    tau_so1   = 30.0181- Time constant for repolarizing tail current below threshold.
    tau_so2   = 0.9957 - Time constant for repolarizing tail current above threshold.
    k_so      = 2.0458 - Parameter controlling nonlinearity in tau_so.
    u_so      = 0.65   - Parameter controlling nonlinearity in tau_so.
    tau_s1    = 2.7342 - Time constant for s below threshold.
    tau_s2    = 16     - Time constant for s above threshold.
    k_s       = 2.0994 - Parameter for tanh activation of s variable.
    u_s       = 0.9087 - Parameter for tanh activation of s variable.
    tau_si    = 1.8875 - Time constant for slow inward current (J_si).
    tau_w_inf = 0.07   - Slope of w∞ below threshold.
    w_inf_    = 0.94   - Asymptotic value of w∞ above threshold.
    """

    return {
        "u_o": 0.0,
        "u_u": 1.55,
        "theta_v": 0.3,
        "theta_w": 0.13,
        "theta_v_m": 0.006,
        "theta_o": 0.006,
        "tau_v1_m": 60.0,
        "tau_v2_m": 1150.0,
        "tau_v_p": 1.4506,
        "tau_w1_m": 60.0,
        "tau_w2_m": 15.0,
        "k_w_m": 65.0,
        "u_w_m": 0.03,
        "tau_w_p": 200.0,
        "tau_fi": 0.11,
        "tau_o1": 400.0,
        "tau_o2": 6.0,
        "tau_so1": 30.0181,
        "tau_so2": 0.9957,
        "k_so": 2.0458,
        "u_so": 0.65,
        "tau_s1": 2.7342,
        "tau_s2": 16.0,
        "k_s": 2.0994,
        "u_s": 0.9087,
        "tau_si": 1.8875,
        "tau_w_inf": 0.07,
        "w_inf_": 0.94
    }


def calc_rhs(J_fi, J_so, J_si) -> float:
    """
    Computes the right-hand side of the model.
    """
    return -J_fi - J_so - J_si


def calc_v(v, u, theta_v, v_inf, tau_v_m, tau_v_p):
    """
    Calculates the fast inactivation gate variable `v`.

    The variable `v` models the fast sodium channel inactivation.
    It follows a piecewise ODE with different dynamics depending
    on whether the membrane potential `u` is above or below `theta_v`.

    Parameters
    ----------
    v : float
        Current value of the v gate.
    u : float
        Current membrane potential.
    theta_v : float
        Threshold for switching recovery behavior.
    v_inf : float
        Steady-state value of v.
    tau_v_m : float
        Time constant for activation (u < threshold).
    tau_v_p : float
        Time constant for decay (u >= threshold).

    Returns
    -------
    float
        Updated value of the v gate.
    """
    return (v_inf - v)/tau_v_m if (u - theta_v) < 0 else -v/tau_v_p


def calc_w(w, u, theta_w, w_inf, tau_w_m, tau_w_p):
    """
    Calculates the slow gating variable `w`.

    The variable `w` represents calcium/potassium channel gating.
    It has different recovery dynamics below and above the threshold `theta_w`.

    Parameters
    ----------
    w : float
        Current value of the w gate.
    u : float
        Membrane potential.
    theta_w : float
        Threshold for switching between time constants.
    w_inf : float
        Steady-state value of w.
    tau_w_m : float
        Time constant for approach to w_inf (u < threshold).
    tau_w_p : float
        Time constant for decay (u >= threshold).

    Returns
    -------
    float
        Updated value of the w gate.
    """
    return (w_inf - w)/tau_w_m if (u - theta_w) < 0 else -w/tau_w_p


def calc_s(s, u, tau_s, k_s, u_s, tanh=math.tanh):
    """
    Calculates the slow variable `s`, related to calcium dynamics.

    The variable `s` evolves toward a tanh-based steady-state function of `u`.

    Parameters
    ----------
    s : float
        Current value of the s variable.
    u : float
        Membrane potential.
    tau_s : float
        Time constant.
    k_s : float
        Slope of the tanh function.
    u_s : float
        Midpoint of the tanh function.
    tanh : function, optional
        Hyperbolic tangent function (default is math.tanh).

    Returns
    -------
    float
        Updated value of the s variable.
    """
    return ((1 + tanh(k_s*(u - u_s)))/2 - s)/tau_s


def calc_Jfi(u, v, theta_v, u_u, tau_fi):
    """
    Computes the fast inward sodium current (J_fi).

    Active when membrane potential exceeds `theta_v`.
    Models rapid depolarization due to sodium influx.

    Parameters
    ----------
    u : float
        Membrane potential.
    v : float
        Fast gating variable.
    theta_v : float
        Activation threshold.
    u_u : float
        Upper limit for depolarization.
    tau_fi : float
        Time constant of the fast inward current.

    Returns
    -------
    float
        Current value of J_fi.
    """
    H = 1.0 if (u - theta_v) >= 0 else 0.0
    return -v*H*(u - theta_v)*(u_u - u)/tau_fi


def calc_Jso(u, u_o, theta_w, tau_o, tau_so):
    """
    Computes the slow outward current (J_so).

    Consists of a linear repolarization component below `theta_w`
    and a constant component above.

    Parameters
    ----------
    u : float
        Membrane potential.
    u_o : float
        Resting potential (offset).
    theta_w : float
        Threshold for switching between components.
    tau_o : float
        Time constant below threshold.
    tau_so : float
        Time constant above threshold.

    Returns
    -------
    float
        Current value of J_so.
    """
    H = 1.0 if (u - theta_w) >= 0 else 0.0
    return (u - u_o)*(1 - H)/tau_o + H/tau_so


def calc_Jsi(u, w, s, theta_w, tau_si):
    """
    Computes the slow inward current (J_si), active during plateau phase.

    Active only when `u > theta_w` and controlled by gating variables `w` and `s`.

    Parameters
    ----------
    u : float
        Membrane potential.
    w : float
        Slow gating variable.
    s : float
        Calcium-related variable.
    theta_w : float
        Threshold for activation.
    tau_si : float
        Time constant of slow inward current.

    Returns
    -------
    float
        Current value of J_si.
    """
    H = 1.0 if (u - theta_w) >= 0 else 0.0
    return -H*w*s/tau_si


def calc_tau_v_m(u, theta_v_m, tau_v1_m, tau_v2_m):
    """
    Selects time constant for v gate depending on membrane potential.

    Returns `tau_v1_m` below `theta_v_m`, and `tau_v2_m` above.

    Returns
    -------
    float
        Time constant for v gate.
    """
    return tau_v1_m if (u - theta_v_m) < 0 else tau_v2_m


def calc_tau_w_m(u, tau_w1_m, tau_w2_m, k_w_m, u_w_m, tanh=math.tanh):
    """
    Computes smooth transition time constant for w gate using tanh.

    Parameters
    ----------
    u : float
        Membrane potential.
    tau_w1_m : float
        Time constant below transition.
    tau_w2_m : float
        Time constant above transition.
    k_w_m : float
        Slope of the transition.
    u_w_m : float
        Midpoint of the transition.
    tanh : function, optional
        Hyperbolic tangent function (default is math.tanh).

    Returns
    -------
    float
        Blended time constant for w gate.
    """
    return tau_w1_m + (tau_w2_m - tau_w1_m)*(1 + tanh(k_w_m*(u - u_w_m)))/2


def calc_tau_so(u, tau_so1, tau_so2, k_so, u_so, tanh=math.tanh):
    """
    Computes tau_so using a sigmoidal transition between two values.

    Parameters
    ----------
    u : float
        Membrane potential.
    tau_so1 : float
        Time constant below transition.
    tau_so2 : float
        Time constant above transition.
    k_so : float
        Slope of the transition.
    u_so : float
        Midpoint of the transition.
    tanh : function, optional
        Hyperbolic tangent function (default is math.tanh).

    Returns
    -------
    float
        Blended time constant tau_so.
    """
    return tau_so1 + (tau_so2 - tau_so1)*(1 + tanh(k_so*(u - u_so)))/2


def calc_tau_s(u, tau_s1, tau_s2, theta_w):
    """
    Selects tau_s based on threshold.

    Parameters
    ----------
    u : float
        Membrane potential.
    tau_s1 : float
        Time constant below threshold.
    tau_s2 : float
        Time constant above threshold.
    theta_w : float
        Threshold for switching.

    Returns
    -------
    float
        Selected time constant tau_s.
    """
    return tau_s1 if (u - theta_w) < 0 else tau_s2


def calc_tau_o(u, tau_o1, tau_o2, theta_o):
    """
    Selects tau_o based on threshold.

    Parameters
    ----------
    u : float
        Membrane potential.
    tau_o1 : float
        Time constant below threshold.
    tau_o2 : float
        Time constant above threshold.
    theta_o : float
        Threshold for switching.

    Returns
    -------
    float
        Selected time constant tau_o.
    """
    return tau_o1 if (u - theta_o) < 0 else tau_o2
    

def calc_v_inf(u, theta_v_m):
    """
    Computes the value of v based on membrane potential.

    Parameters
    ----------
    u : float
        Membrane potential.
    theta_v_m : float
        Threshold for v_inf.

    Returns
    -------
    float
        Steady-state value of v.
    """
    return 1.0 if u < theta_v_m else 0.0


def calc_w_inf(u, theta_o, tau_w_inf, w_inf_):
    """
    Computes the value of w based on membrane potential.

    Parameters
    ----------
    u : float
        Membrane potential.
    theta_o : float
        Threshold for w_inf.
    tau_w_inf : float
        Time constant for w_inf.
    w_inf_ : float
        Steady-state value of w above threshold.

    Returns
    -------
    float
        Steady-state value of w.
    """
    return 1 - u/tau_w_inf if (u - theta_o) < 0 else w_inf_