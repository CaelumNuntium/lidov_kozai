import numpy as np
from matplotlib import pyplot as plt


def to_interval(x, interval):
    l = interval[1] - interval[0]
    periods = (x - interval[0]) / l
    return x - np.floor(periods) * l


def plot_elements(name, t_arr, a_set, n_set, a_color, n_color, t_unit, a_unit, filename):
    plt.rcParams["figure.constrained_layout.use"] = True
    plt.xlabel(f"$t$, {t_unit}")
    if name == "a":
        plt.ylabel(f"$a$, {a_unit}")
    if name == "x" or name == "y":
        plt.ylabel(f"${name}$")
    if name == "g":
        plt.ylabel("$g$, degrees")
    if name == "Omega":
        plt.ylabel("$\\Omega$, degrees")
    if name == "p_long":
        plt.ylabel("$\\varpi$, degrees")
    plt.plot(t_arr, n_set[name], color=n_color, linewidth=0.75)
    plt.plot(t_arr, a_set[name], color=a_color, linewidth=0.75)
    f = plt.gcf()
    f.savefig(filename, dpi=450)
    plt.cla()
    plt.clf()


def read_config(config_name):
    parameters = {}
    with open(config_name, "r") as config:
        for st in config:
            if '#' in st:
                s1 = st.split('#')
            else:
                s1 = st
            if '=' in s1:
                parameter = s1.split('=')[0].strip()
                value = s1.split('=')[1].strip().strip('"').strip("'")
                parameters[parameter] = value
    return parameters


params = read_config("config.ini")
n = int(params["N_steps"]) + 1
stat_point_type = int(params["stat_point_type"])
a_u = params["a_unit"]
t_u = params["t_unit"]
a_elements = {"a": np.zeros(n), "x": np.zeros(n), "y": np.zeros(n), "Omega": np.zeros(n), "g": np.zeros(n),
              "p_long": np.zeros(n)}
n_elements = {"a": np.zeros(n), "x": np.zeros(n), "y": np.zeros(n), "Omega": np.zeros(n), "g": np.zeros(n),
              "p_long": np.zeros(n)}
times = np.zeros(n)
i = 0
with open("result_analytical.dat") as inp:
    for s in inp:
        if '#' not in s:
            s_arr = s.split()
            times[i] = float(s_arr[0])
            a_elements["a"][i] = float(s_arr[1])
            a_elements["x"][i] = float(s_arr[2])
            a_elements["y"][i] = float(s_arr[3])
            a_elements["g"][i] = float(s_arr[4])
            a_elements["Omega"][i] = to_interval(float(s_arr[5]), (-180, 180))
            a_elements["p_long"][i] = to_interval(float(s_arr[6]), (-180, 180))
            i = i + 1
i = 0
with open("result_numerical.dat") as inp:
    for s in inp:
        if '#' not in s:
            s_arr = s.split()
            n_elements["a"][i] = float(s_arr[1])
            n_elements["x"][i] = float(s_arr[2])
            n_elements["y"][i] = float(s_arr[3])
            n_elements["g"][i] = float(s_arr[4])
            n_elements["Omega"][i] = to_interval(float(s_arr[5]), (-180, 180))
            n_elements["p_long"][i] = to_interval(float(s_arr[6]), (-180, 180))
            i = i + 1
plot_elements("a", times, a_elements, n_elements, "red", "green", t_u, a_u, "semimajor_axis.png")
plot_elements("x", times, a_elements, n_elements, "red", "green", t_u, a_u, "x.png")
plot_elements("y", times, a_elements, n_elements, "red", "green", t_u, a_u, "y.png")
if stat_point_type == 1 or stat_point_type == 3:
    plot_elements("Omega", times, a_elements, n_elements, "red", "green", t_u, a_u, "asc_node.png")
if stat_point_type == 2:
    plot_elements("p_long", times, a_elements, n_elements, "red", "green", t_u, a_u, "pericenter_longitude.png")
