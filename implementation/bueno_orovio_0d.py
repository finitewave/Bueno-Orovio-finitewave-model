"""
This module provides a simple interface to run the model in a 0D setting,
i.e., without spatial dimensions. It includes class for defining stimulation protocols
and a class for the 0D model itself.

"""

from bueno_orovio import ops


class Stimulation:
    """
    Stimulus protocol for the 0D model.

    Parameters
    ----------
    t_start : float
        Start time (ms) of the first stimulus window.
    duration : float
        Duration (ms) of a single pulse.
    amplitude : float
        Pulse amplitude in the same units as du/dt contribution (typically "units/ms").

    Method
    ------
    stim(t: float) -> float
        Returns the instantaneous stimulus value at time t.

    """

    def __init__(self, t_start: float, duration: float, amplitude: float):
        self.t_start = t_start
        self.duration = duration
        self.amplitude = amplitude

    def stim(self, t: float) -> float:
        return self.amplitude if self.t_start <= t < self.t_start + self.duration else 0.0


class BuenoOrovio0D:
    """
    Bueno-Orovio OD implementation.

    Parameters
    ----------

    dt : float
        Time step size (ms).
    stimulations : list[Stimulation]
        List of stimulation protocols to apply during the simulation.

    Attributes
    ----------
    variables : dict[str, float]
        Current state variables of the model.
    parameters : dict[str, float]
        Model parameters.
    history : dict[str, list[float]]
        Time history of state variables for post-processing.
    
    Methods
    -------
    step(i: int)
        Perform a single time step update.
    run(t_max: float)
        Run the simulation up to time t_max.
    """
    def __init__(self, dt: float, stimulations: list[Stimulation]):
        self.dt = dt
        self.stimulations = stimulations
        self.variables = ops.get_variables()
        self.parameters = ops.get_parameters()
        self.history = {s: [] for s in self.variables}

    def step(self, i: int):
        """
        Perform a single time step update.

        Parameters
        ----------
        i : int
            Current time step index.
        """
        u_old = self.variables["u"]
        v_old = self.variables["v"]
        w_old = self.variables["w"]
        s_old = self.variables["s"]

        v_inf = ops.calc_v_inf(u_old, self.parameters["theta_v_m"])
        tau_v_m = ops.calc_tau_v_m(
            u_old,
            self.parameters["theta_v_m"],
            self.parameters["tau_v1_m"],
            self.parameters["tau_v2_m"],
        )
        dv = ops.calc_v(
            v_old,
            u_old,
            self.parameters["theta_v"],
            v_inf,
            tau_v_m,
            self.parameters["tau_v_p"],
        )

        w_inf = ops.calc_w_inf(
            u_old,
            self.parameters["theta_o"],
            self.parameters["tau_w_inf"],
            self.parameters["w_inf_"],
        )
        tau_w_m = ops.calc_tau_w_m(
            u_old,
            self.parameters["tau_w1_m"],
            self.parameters["tau_w2_m"],
            self.parameters["k_w_m"],
            self.parameters["u_w_m"],
        )
        dw = ops.calc_w(
            w_old,
            u_old,
            self.parameters["theta_w"],
            w_inf,
            tau_w_m,
            self.parameters["tau_w_p"],
        )

        tau_s = ops.calc_tau_s(
            u_old,
            self.parameters["tau_s1"],
            self.parameters["tau_s2"],
            self.parameters["theta_w"],
        )
        ds = ops.calc_s(
            s_old,
            u_old,
            tau_s,
            self.parameters["k_s"],
            self.parameters["u_s"],
        )

        J_fi = ops.calc_Jfi(
            u_old,
            v_old,
            self.parameters["theta_v"],
            self.parameters["u_u"],
            self.parameters["tau_fi"],
        )

        tau_o = ops.calc_tau_o(
            u_old,
            self.parameters["tau_o1"],
            self.parameters["tau_o2"],
            self.parameters["theta_o"],
        )
        tau_so = ops.calc_tau_so(
            u_old,
            self.parameters["tau_so1"],
            self.parameters["tau_so2"],
            self.parameters["k_so"],
            self.parameters["u_so"],
        )
        J_so = ops.calc_Jso(
            u_old,
            self.parameters["u_o"],
            self.parameters["theta_w"],
            tau_o,
            tau_so,
        )

        J_si = ops.calc_Jsi(
            u_old,
            w_old,
            s_old,
            self.parameters["theta_w"],
            self.parameters["tau_si"],
        )

        stim_current = sum(stim.stim(t=self.dt * i) for stim in self.stimulations)

        du = ops.calc_rhs(J_fi, J_so, J_si) + stim_current

        self.variables["u"] = u_old + self.dt * du
        self.variables["v"] = v_old + self.dt * dv
        self.variables["w"] = w_old + self.dt * dw
        self.variables["s"] = s_old + self.dt * ds

    def run(self, t_max: float):
        """
        Run the simulation up to time t_max.
        
        Parameters
        ----------
        t_max : float
            Maximum simulation time.
        """
        n_steps = int(round(t_max/self.dt))
        for i in range(n_steps):
            self.step(i)
            for s in self.variables:
                self.history[s].append(self.variables[s])