# Stationary points of the Lidov - Kozai problem

The program compares the analytical solution of the Lidov - Kozai problem for stationary points with the numerical solution of the three-body problem.

Types of stationary points:
1. $x = 0$, $y = 1 - \gamma$ $\implies$ $\dot{x} = 0$, $\dot{y} = 0$, $\dot{\Omega} \neq 0$, circular orbit;
2. $y = 0$, $x = 1 - \gamma$ $\implies$ $\dot{x} = 0$, $\dot{y} = 0$, $\dot{\varpi} \neq 0$;
3. $\cos(g) = 0$, $x = 1 - \sqrt{5 \gamma / 3}$, $y = 1 - \sqrt{3 \gamma / 5}$ $\implies$ $\dot{x} = 0$, $\dot{y} = 0$, $\dot{g} = 0$, $\dot{\Omega} \neq 0$.
## Installation

1. Clone the repository:
```bash
git clone https://github.com/CaelumNuntium/lidov_kozai.git
cd lidov_kozai
```

2. Build program:
```bash
make
```
## Usage

1. First, set the parameters (stationary point type, $\gamma$; $Gm_0$, $Gm_1$, $Gm_2$, $a'$, $e'$, $i'$, $g'$, $\Omega'$, $M'(0)$, $a(0)$, $M(0)$ ) in *config.ini*

2. Run the program:

    * On Linux:
    ```bash
    ./lidov_kozai
    ```
    
    * On Windows:
    ```cmd
    lidov_kozai.exe
    ```

3. Visualization:
```bash
python plot.py
```
Files *semimajor_axis.png*, *x.png*, *y.png*, *asc_node.png*, *pericenter_longitude.png* will contain graphs of the osculating elements $a$, $x$, $y$, $\Omega$, $\varpi$ of the satellite's orbit in the numerical and analytical models of motion.
