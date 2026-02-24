"""
APMA 4302 - HW2 Q3: Plot BVP solution and convergence.

ADDITIONS FOR HW2:
  (b) read_hdf5_vec: read PETSc HDF5 (petsc4py) for u, uexact, f.
  (b) plot_bvp_solution: plot numerical vs exact and pointwise error (twin axis).
  (c) run_bvp_and_get_error: run ./bvp for given m,k, parse "Relative error (2-norm):" from stdout.
  (c) run_convergence_study: gamma=0, k=1,5,10, m=40,80,...,1280; log-log error vs h; polyfit for order.

Usage:
  python plot_bvp.py                    # plot single solution from bvp_solution.h5
  python plot_bvp.py convergence [./bvp] # run convergence study and plot order
"""
import numpy as np
import subprocess
import re
import os

try:
    from petsc4py import PETSc
except (ImportError, AttributeError):
    PETSc = None


def read_hdf5_vec(filename, vec_name):
    """[3(b)] Read vector from PETSc HDF5 (bvp_solution.h5). Requires petsc4py."""
    if PETSc is None:
        raise RuntimeError("petsc4py required to read HDF5; source env_plot_bvp or set PYTHONPATH to PETSc lib")
    viewer = PETSc.Viewer().createHDF5(filename, 'r')
    vec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    vec.setName(vec_name)
    vec.load(viewer)
    arr = vec.getArray().copy()
    vec.destroy()
    viewer.destroy()
    return arr


def plot_bvp_solution(x, u_numeric, u_exact, outfile=None):
    """[3(b)] Plot numerical vs exact solution and pointwise error (right axis)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    ax1.plot(x, u_numeric, 'b-', label='Numerical Solution', linewidth=2)
    ax1.plot(x, u_exact, 'r--', label='Exact Solution', linewidth=2)
    ax1.set_xlabel('X', fontsize=14)
    ax1.set_ylabel('u(X)', fontsize=14)
    ax1.set_title('BVP Numerical vs Exact Solution\n' + r'$u(x)=\sin(k\pi x)+c/8$, $f=(k^2\pi^2+\gamma)\sin(k\pi x)+\gamma c/8$', fontsize=14)
    ax1.legend(fontsize=12)
    ax1.grid(True)
    err = u_numeric - u_exact
    ax2.plot(x, err, 'g--', label='Error', linewidth=1)
    ax2.set_ylabel('Error', fontsize=14)
    ax2.legend(loc='lower right', fontsize=12)
    emax = np.max(np.abs(err))
    if emax > 0:
        ax2.set_ylim(-emax * 1.5, emax * 1.5)
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"Saved {outfile}")
    else:
        plt.show()
    plt.close()


def run_bvp_and_get_error(bvp_exe, m, k, gamma=0.0, c=1.0, nproc=1):
    """[3(c)] Run BVP, parse 'Relative error (2-norm):' from stdout. Returns (h, rel_err) or (None, None)."""
    h = 1.0 / m
    cmd = [
        'mpiexec', '-np', str(nproc), bvp_exe,
        '-bvp_m', str(m), '-bvp_gamma', str(gamma), '-bvp_k', str(k), '-bvp_c', str(c),
        '-ksp_rtol', '1e-10', '-ksp_atol', '1e-12',
    ]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=120, cwd=os.path.dirname(os.path.abspath(bvp_exe)) or '.')
        match = re.search(r'Relative error \(2-norm\):\s*([\d.eE+-]+)', out.stdout or out.stderr or '')
        if match:
            return (h, float(match.group(1)))
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        pass
    return (h, None)


def run_convergence_study(bvp_exe='./bvp', m_list=None, k_list=None, gamma=0.0, nproc=1, outfile=None):
    """[3(c)] Error vs h (log-log), gamma=0, k=1,5,10, m=40,...,1280; report order via log-log fit."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    m_list = m_list or [40, 80, 160, 320, 640, 1280]  # m = 40, 80, 160, ..., 1280
    k_list = k_list or [1, 5, 10]
    results = {}
    for k in k_list:
        results[k] = []
        for m in m_list:
            h, err = run_bvp_and_get_error(bvp_exe, m, k, gamma=gamma, nproc=nproc)
            if err is not None:
                results[k].append((h, err))
        results[k] = sorted(results[k], key=lambda x: x[0])

    # Plot: error vs h (log-log) for each k
    fig, ax = plt.subplots(figsize=(8, 6))
    print('Q3(c) Order of convergence (gamma=0):')
    has_series = False
    for k in k_list:
        if not results[k]:
            continue
        has_series = True
        h_arr = np.array([x[0] for x in results[k]])
        e_arr = np.array([x[1] for x in results[k]])
        ax.loglog(h_arr, e_arr, 'o-', label=f'k={k}', base=2)
        if len(h_arr) >= 2:
            logh = np.log(h_arr)
            loge = np.log(e_arr + 1e-30)
            order = np.polyfit(logh, loge, 1)[0]
            print(f'  k={k}: order ≈ {order:.3f}')
    if not has_series:
        print('  (No BVP runs succeeded; run on cluster or with working MPI to get convergence data.)')
    ax.set_xlabel('h', fontsize=14)
    ax.set_ylabel('Relative error (2-norm)', fontsize=14)
    ax.set_title('Q3(c): Convergence of error vs h (gamma=0)', fontsize=14)
    if has_series:
        ax.legend(fontsize=12)
    ax.grid(True, which='major')
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"Saved {outfile}")
    else:
        plt.show()
    plt.close()
    return results


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'convergence':
        # Run convergence study: python plot_bvp.py convergence [path/to/bvp]
        bvp_exe = sys.argv[2] if len(sys.argv) > 2 else './bvp'
        out = os.environ.get('BVP_CONVERGENCE_OUT', 'bvp_convergence.pdf')
        nproc = int(os.environ.get('BVP_NPROC', '1'))
        run_convergence_study(bvp_exe=bvp_exe, outfile=out, nproc=nproc)
    else:
        # Default: read HDF5 and plot single solution (after one BVP run)
        h5_filename = 'bvp_solution.h5'
        if not os.path.isfile(h5_filename):
            print(f'No {h5_filename} found. Run: mpiexec -np 1 ./bvp -options_file options_file')
            print('Or run convergence study: python plot_bvp.py convergence [./bvp]')
            sys.exit(1)
        u = read_hdf5_vec(h5_filename, 'u')
        u_exact = read_hdf5_vec(h5_filename, 'uexact')
        x = np.linspace(0, 1, len(u))
        out = os.environ.get('BVP_PLOT_OUT', 'bvp_solution.pdf')
        plot_bvp_solution(x, u, u_exact, outfile=out)
