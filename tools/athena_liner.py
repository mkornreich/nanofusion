"""2-D (r-z) and 3-D MHD simulations of the thin, ice-lined liner with Athena++: does it hold together to ~170x?

The 0-D liner model (nanofusion/improved/zpinch.py) has no liner instability, so the simulator's best liner
(docs/liner_experiment.md) converges ~170x. Here Athena++ (Princeton's open MHD code) implodes the same target
with the magneto-Rayleigh-Taylor instability free to grow from a ripple on the liner's surface: in r-z, and in a
3-D quick look (half the r-z resolution, the full circle) that also lets helical modes grow.

Setup: tools/athena/liner.cpp (the problem generator) and tools/athena/athinput.liner (the input). Athena++
itself is not part of this repository: clone github.com/PrincetonUniversity/athena next to it (../athena),
copy liner.cpp into its src/pgen, and build with

    python3 configure.py --prob=liner --coord=cylindrical -b --nscalars=3 --flux=hlld -omp && make -j

Then

    python3 tools/athena_liner.py run [case ...]   # default: the three rippled 2-D cases (~1.5 h on 12 cores)
    python3 tools/athena_liner.py report           # docs/results/liner_athena/summary.md and results.json

Limits: ideal MHD, one ideal gas (gamma 5/3) for beryllium, ice and gas (no strength, no real equation of
state), the axial field frozen in (no Nernst loss), no radiation, conduction or burn, a reflecting wall at
r = 10 um instead of the axis, and 0.5 mm of liner with periodic ends. It tests whether the liner stays intact
and how far it squeezes the fuel, not the yield.

The yield: the "physics" cases add what a yield needs (tools/athena/liner.cpp, compiled with the physics
build): D-T burn, alpha heating, radiation, resistivity, electron heat conduction, separate ion and electron
temperatures, the Nernst effect, the current from a capacitor bank, and an equation of state with cold curves
and ionisation (tools/athena/liner_eos.cpp, liner_physics.hpp). They need a patched Athena++ in its own tree
(../athena-physics, so the ideal-MHD binary is left alone):

    git -C ../athena worktree add ../athena-physics          # or a second clone
    python3 tools/athena_liner.py build-physics               # patch, copy, configure, make
    python3 tools/athena_liner.py run "physics flat"          # ~1-2 h on 2 threads
    python3 tools/athena_liner.py report

The physics runs write restart files every 5 ns; `run` resumes a case from its newest one, and
`job "<case>"` prints a stand-alone script that does the same on another machine. The report then adds the yield, the burn-weighted temperatures and field, and the electric gain counted as the
0-D simulator counts it (docs/simulator.md), for each physics case. Each "physics flat, no X" case switches
one piece off, to see what it is worth.
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ATHENA = Path(os.environ.get("ATHENA_DIR", HERE.parents[1] / "athena"))
RUNS = ATHENA / "runs"
INPUT = HERE / "athena" / "athinput.liner"
OUT = HERE.parent / "docs" / "results" / "liner_athena"
R_GAS0 = 4.495 * (1 - 1 / 12) - 0.5833  # mm, the gas column's starting radius
T_UNIT_KEV = 12.95  # (p / rho) in code units -> keV, for D-T with one temperature for ions and electrons
# half the r-z resolution (6 um cells at the axis, 15.6 um along z), axial modes 1-8 and azimuthal modes 1-8
HALF = ["mesh/nx1=128", "mesh/x1rat=1.026169", "mesh/nx3=32", "meshblock/nx1=32", "meshblock/nx3=16",
        "problem/modes=8", "problem/mmax=8"]
PHYS = Path(os.environ.get("ATHENA_PHYSICS_DIR", HERE.parents[1] / "athena-physics"))
PHYS_INPUT = HERE / "athena" / "athinput.liner_physics"
PHYS_CONFIG = ["--prob=liner", "--coord=cylindrical", "-b", "--nscalars=5", "--flux=hlld",
               "--eos=general/liner_eos", "-omp"]
FLAT = ["mesh/nx3=8", "meshblock/nx3=8", "mesh/x3max=0.0625"]   # 1/8 of the length: no ripple to resolve
CASES = {  # name: rms ripple of the liner's outer surface (mm), extra Athena++ arguments, OpenMP threads
    "flat": {"ripple": 0.0, "extra": ["mesh/nx3=8", "meshblock/nx3=8", "mesh/x3max=0.0625"]},
    "ripple 0.1 um": {"ripple": 1e-4, "extra": []},
    "ripple 1 um": {"ripple": 1e-3, "extra": []},
    "ripple 10 um": {"ripple": 1e-2, "extra": []},
    # 3-D: 128 x 64 x 32 cells over the full circle, 144 seed modes (8 axial, 136 azimuthal or helical).
    # Run to 180 ns it took 11.7 h on 11 threads, 85% of it after stagnation (~140 ns): once beryllium reaches
    # the axis, the thin cells there force steps ~40x smaller. Add "time/tlim=150" to stop it after stagnation.
    "3-D quick look": {"ripple": 1e-3, "threads": 11,
                       "extra": [*HALF, "mesh/nx2=64", "meshblock/nx2=16", "output3/dt=20", "output4/dt=2.5"]},
    # its 2-D twin: the same grid in r-z and the same 8 axial modes, without the 136 others
    "2-D twin": {"ripple": 1e-3, "threads": 1, "extra": [*HALF, "problem/mseed=0"]},
    # with the physics a yield needs (athinput.liner_physics: everything on); then one piece off at a time
    "physics flat": {"ripple": 0.0, "physics": True, "threads": 2, "extra": FLAT},
    "physics flat, no alpha heating": {"ripple": 0.0, "physics": True, "threads": 2,
                                       "extra": [*FLAT, "problem/alpha=-1"]},
    "physics flat, alphas stop where made": {"ripple": 0.0, "physics": True, "threads": 2,
                                             "extra": [*FLAT, "problem/alpha=0"]},
    "physics flat, no radiation": {"ripple": 0.0, "physics": True, "threads": 2,
                                   "extra": [*FLAT, "problem/radiation=false"]},
    "physics flat, no conduction": {"ripple": 0.0, "physics": True, "threads": 2,
                                    "extra": [*FLAT, "problem/conduction=false"]},
    "physics flat, one temperature": {"ripple": 0.0, "physics": True, "threads": 2,
                                      "extra": [*FLAT, "problem/two_temp=false"]},
    "physics flat, no Nernst": {"ripple": 0.0, "physics": True, "threads": 2,
                                "extra": [*FLAT, "problem/nernst=false"]},
    "physics flat, no resistivity": {"ripple": 0.0, "physics": True, "threads": 2,
                                     "extra": [*FLAT, "problem/resistivity=false"]},
    "physics flat, prescribed current": {"ripple": 0.0, "physics": True, "threads": 2,
                                         "extra": [*FLAT, "problem/circuit=false"]},
    # the full r-z grid with a 1 um ripple (about 4x the flat case's cells)
    "physics ripple 1 um": {"ripple": 1e-3, "physics": True, "threads": 4, "extra": []},
}
DEFAULT_RUN = ["ripple 0.1 um", "ripple 1 um", "ripple 10 um"]


def folder(name: str) -> Path:
    runs = PHYS / "runs" if CASES[name].get("physics") else RUNS
    return runs / name.replace(" ", "_").replace(",", "")


def z_length_mm(name: str) -> float:
    """The simulated length of liner (mm): mesh/x3max, 0.5 unless the case changes it."""
    for arg in CASES[name]["extra"]:
        if arg.startswith("mesh/x3max="):
            return float(arg.split("=")[1])
    return 0.5


def build_physics() -> int:
    """Patch the Athena++ tree at PHYS (once), copy the liner files in, and build it from clean."""
    if not (PHYS / "configure.py").exists():
        print(f"no Athena++ tree at {PHYS}: git -C {ATHENA} worktree add {PHYS}")
        return 1
    patch = HERE / "athena" / "athena-physics.patch"
    applied = subprocess.run(["git", "apply", "--check", "-R", str(patch)], cwd=PHYS, capture_output=True)
    if applied.returncode != 0:
        subprocess.run(["git", "apply", str(patch)], cwd=PHYS, check=True)
    for src, dst in [("liner.cpp", "src/pgen/liner.cpp"), ("liner_eos.cpp", "src/eos/general/liner_eos.cpp"),
                     ("liner_physics.hpp", "src/eos/general/liner_physics.hpp")]:
        (PHYS / dst).write_text((HERE / "athena" / src).read_text())
    subprocess.run([sys.executable, "configure.py", *PHYS_CONFIG], cwd=PHYS, check=True)
    # Athena++'s Makefile does not track headers: always build from clean
    subprocess.run(["make", "clean"], cwd=PHYS, check=True, capture_output=True)
    subprocess.run(["make", f"-j{min(8, os.cpu_count() or 1)}"], cwd=PHYS, check=True)
    return 0


TLIM = 180  # ns


def command(name: str, threads: int | None = None, binary: str | None = None, inp: str | None = None,
            restart: str | None = None) -> list[str]:
    """The Athena++ command line for a case (binary and input default to this machine's; threads, the
    case's own count); with `restart`, resume from that restart file (its parameters, the mesh included)."""
    case = CASES[name]
    tree, default_inp = (PHYS, PHYS_INPUT) if case.get("physics") else (ATHENA, INPUT)
    binary = binary or str(tree / "bin" / "athena")
    # Athena++ takes its OpenMP thread count from mesh/num_threads (default 1), not OMP_NUM_THREADS
    threads = threads if threads is not None else case.get("threads", 4)
    common = [f"time/tlim={TLIM}", "time/ncycle_out=5000", f"mesh/num_threads={threads}"]
    if restart:
        return [binary, "-r", restart, *common]
    return [binary, "-i", inp or str(default_inp), f"problem/ripple={case['ripple']}", *common, *case["extra"]]


def newest_restart(where: Path) -> Path | None:
    files = sorted(where.glob("liner.*.rst"))
    return files[-1] if files else None


def shell_job(name: str, threads: int = 1, binary: str = "../../bin/athena",
              inp: str = "../../athinput.liner_physics") -> str:
    """A self-contained sh script that runs (or resumes) a case in its own directory, for another
    machine: start it from that directory; it resumes from the newest restart file if there is one."""
    fresh = shlex.join(command(name, threads, binary, inp))
    again = shlex.join(command(name, threads, binary, restart="RESTART")).replace("RESTART", '"$r"')
    return ('cd "$(dirname "$0")" || exit 1\n'
            'r=$(ls -t liner.*.rst 2>/dev/null | head -1)\n'
            f'if [ -n "$r" ]; then exec nice -n 19 {again} >> run.log 2>&1\n'
            f'else exec nice -n 19 {fresh} > run.log 2>&1; fi\n')


def run(names: list[str], threads: int = 4) -> None:
    """Start the cases in parallel (each in runs/<name>), and wait for them."""
    procs = []
    for name in names:
        where = folder(name)
        where.mkdir(parents=True, exist_ok=True)
        rst = newest_restart(where)   # resume a run that stopped part-way
        cmd = command(name, CASES[name].get("threads", threads), restart=str(rst) if rst else None)
        with open(where / "run.log", "a" if rst else "w") as log:  # the child keeps its own copy of the file
            procs.append(subprocess.Popen(cmd, cwd=where, stdout=log, stderr=subprocess.STDOUT))
    for p in procs:
        p.wait()


def history(where: Path) -> dict:
    lines = (where / "liner.hst").read_text().splitlines()
    header = next(line for line in lines if line.startswith("# [1]"))
    names = [part.split("=", 1)[1].strip() for part in header[1:].split("[")[1:] for part in [part.split("]", 1)[1]]]
    data = np.atleast_2d(np.loadtxt([line for line in lines if not line.startswith("#")]))
    # a resumed run writes again the rows after its restart file: keep the last of each time
    keep = np.ones(len(data), bool)
    t_min = np.inf
    for k in range(len(data) - 1, -1, -1):
        if data[k, 0] >= t_min:
            keep[k] = False
        t_min = min(t_min, data[k, 0])
    data = data[keep]
    return {name: data[:, k] for k, name in enumerate(names)}


def read_vtk(path: Path):
    sys.path.insert(0, str(ATHENA / "vis" / "python"))
    import athena_read  # Athena++'s own reader, from the Athena++ tree
    return athena_read.vtk(str(path))


def vtk_time(path: Path) -> float:
    with open(path, "rb") as f:
        f.readline()
        line = f.readline().decode(errors="replace")  # "# Athena++ data at time=... cycle=..."
    return float(line.split("time=")[1].split()[0])


def assemble(files: list[Path], names: tuple[str, ...]) -> tuple:
    """One output's VTK blocks -> faces (r, phi, z) and arrays on the whole (z, phi, r) grid."""
    blocks = [read_vtk(f) for f in files]
    # float64 and rounded, so that faces shared by neighbouring blocks match
    faces = [np.unique(np.round(np.concatenate([b[a] for b in blocks]).astype(float), 7)) for a in range(3)]
    cells = [max(len(f) - 1, 1) for f in faces]  # a slice output has one face in phi: one cell
    data = {n: np.zeros((cells[2], cells[1], cells[0])) for n in names}
    for b in blocks:
        start = [int(np.searchsorted(faces[a], round(float(b[a][0]), 7))) for a in range(3)]
        size = [max(len(b[a]) - 1, 1) for a in range(3)]
        sl = tuple(slice(start[a], start[a] + size[a]) for a in (2, 1, 0))
        for n in names:
            data[n][sl] = b[3][n].reshape(size[2], size[1], size[0])
    return faces, data


def snapshots(where: Path) -> list[dict]:
    """Per r-z slice (output 2): the liner's inner surface along z and the liner mass mixed into the fuel."""
    out = []
    stamps = sorted({p.name.split(".")[-2] for p in where.glob("liner.block0.out2.*.vtk")})
    for stamp in stamps:
        files = sorted(where.glob(f"liner.block*.out2.{stamp}.vtk"))
        (r_f, _, _), d = assemble(files, ("rho", "r0", "r1", "r2"))
        rho, liner, gas = d["rho"][:, 0], d["r0"][:, 0], d["r2"][:, 0]  # density; liner and gas shares
        rc = 0.5 * (r_f[1:] + r_f[:-1])
        area = np.pi * (r_f[1:] ** 2 - r_f[:-1] ** 2)
        # the liner's inner surface: the smallest radius where liner material dominates, along each z
        inner = np.array([rc[np.argmax(row > 0.5)] if (row > 0.5).any() else np.nan for row in liner])
        # beryllium that has got through the ice into the hot gas: liner mass in cells that are mostly gas,
        # over the gas mass (the ice lies between them, so numerical smearing at one boundary does not count)
        mostly_gas = gas > 0.5
        gas_mass = (rho * gas * area * mostly_gas).sum()
        mixed = (rho * liner * area * mostly_gas).sum() / gas_mass if gas_mass > 0 else None  # None: none left
        out.append({"time_ns": vtk_time(files[0]), "inner_mean_mm": float(np.nanmean(inner)),
                    "inner_spread_mm": float(np.nanmax(inner) - np.nanmin(inner)),
                    "mixed_share": None if mixed is None else float(mixed)})
    return out


def azimuthal(where: Path) -> list[dict]:
    """Per full 3-D snapshot (output 4): how far the liner's inner surface and the fuel column leave axisymmetry."""
    out = []
    stamps = sorted({p.name.split(".")[-2] for p in where.glob("liner.block0.out4.*.vtk")})
    for stamp in stamps:
        files = sorted(where.glob(f"liner.block*.out4.{stamp}.vtk"))
        (r_f, p_f, _), d = assemble(files, ("rho", "r0", "r2"))
        if len(p_f) < 5:
            continue  # not a 3-D run
        rc, pc = 0.5 * (r_f[1:] + r_f[:-1]), 0.5 * (p_f[1:] + p_f[:-1])
        is_liner = d["r0"] > 0.5  # (z, phi, r)
        inner = np.where(is_liner.any(-1), rc[np.argmax(is_liner, -1)], np.nan)  # (z, phi)
        ring = np.nanmean(inner, axis=1)  # the axisymmetric part, per z
        dev = np.nan_to_num(inner - np.nanmean(inner))
        spec = np.abs(np.fft.fft2(dev))  # over (z, phi)
        spec[:, 0] = 0.0  # leave out m = 0
        kz, km = np.unravel_index(np.argmax(spec), spec.shape)
        nz, nphi = spec.shape
        # the fuel column (gas): its mass-weighted centre and rms radius in each z slice
        w = d["rho"] * d["r2"] * (rc * np.diff(r_f))[None, None, :] * np.diff(p_f)[None, :, None]
        mass = w.sum((1, 2))
        x = (w * rc[None, None, :] * np.cos(pc)[None, :, None]).sum((1, 2)) / mass
        y = (w * rc[None, None, :] * np.sin(pc)[None, :, None]).sum((1, 2)) / mass
        r_rms = np.sqrt((w * rc[None, None, :] ** 2).sum((1, 2)) / mass)
        out.append({"time_ns": vtk_time(files[0]), "inner_mean_um": float(np.nanmean(inner) * 1e3),
                    # lines of sight (z, phi) with beryllium dominating the cell next to the axis: the liner has
                    # collapsed through the fuel there, and its "inner surface" is just the first cell
                    "liner_at_axis_share": float(is_liner[:, :, 0].mean()),
                    "inner_axisym_rms_um": float(np.sqrt(np.nanmean((ring - np.nanmean(ring)) ** 2)) * 1e3),
                    "inner_nonaxisym_rms_um": float(np.sqrt(np.nanmean((inner - ring[:, None]) ** 2)) * 1e3),
                    "top_mode_m_n": [int(km if km <= nphi // 2 else km - nphi), int(kz if kz <= nz // 2 else kz - nz)],
                    "fuel_rms_um": float(np.sqrt((w * rc[None, None, :] ** 2).sum() / w.sum()) * 1e3),
                    "fuel_offset_max_um": float(np.hypot(x, y).max() * 1e3),
                    "fuel_rms_spread_um": float((r_rms.max() - r_rms.min()) * 1e3)})
    return out


BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf",
        "#1c5cab", "#184f95", "#104281", "#0d366b"]  # one-hue sequential ramp, light = low density
SURFACE, INK, MUTED, ORANGE = "#fcfcfb", "#0b0b0b", "#6b6a66", "#eb6834"


def rz_slice(where: Path, t: float) -> tuple:
    """The r-z slice (output 2) nearest time t: its time, faces, density and gas share."""
    times = {f.name.split(".")[-2]: vtk_time(f) for f in where.glob("liner.block0.out2.*.vtk")}
    stamp = min(times, key=lambda k: abs(times[k] - t))
    (r_f, _, z_f), d = assemble(sorted(where.glob(f"liner.block*.out2.{stamp}.vtk")), ("rho", "r2"))
    return times[stamp], r_f, z_f, d["rho"][:, 0], d["r2"][:, 0]


def density_map(ax, where: Path, t: float, r_max: float, title: str):
    tt, r_f, z_f, rho, gas = rz_slice(where, t)
    if z_f[-1] - z_f[0] < 0.4:  # the flat run covers 1/8 of the length and is uniform along it: repeat it
        reps = round(0.5 / (z_f[-1] - z_f[0]))
        rho, gas = np.tile(rho, (reps, 1)), np.tile(gas, (reps, 1))
        z_f = np.linspace(0.0, 0.5, rho.shape[0] + 1)
    from matplotlib.colors import LinearSegmentedColormap
    mesh = ax.pcolormesh(r_f, z_f, np.log10(rho), cmap=LinearSegmentedColormap.from_list("blue", BLUE),
                         vmin=-4, vmax=1.5, shading="flat", rasterized=True)
    rc, zc = 0.5 * (r_f[1:] + r_f[:-1]), 0.5 * (z_f[1:] + z_f[:-1])
    if (gas > 0.5).any():
        ax.contour(rc, zc, gas, levels=[0.5], colors=ORANGE, linewidths=1.0)
    ax.set_xlim(0.0, r_max)
    ax.set_ylim(0.0, 0.5)
    style(ax, f"{title}, {tt:.0f} ns")
    return mesh


def style(ax, title: str) -> None:
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.tick_params(colors=MUTED, labelsize=7, length=2)
    for side in ax.spines.values():
        side.set_color(MUTED)
        side.set_linewidth(0.5)


def cross_section(ax, where: Path, t: float, r_max: float, title: str):
    """A cut across a 3-D run at mid-length (output 4), nearest time t: density, and the gas's edge."""
    from matplotlib.colors import LinearSegmentedColormap
    times = {f.name.split(".")[-2]: vtk_time(f) for f in where.glob("liner.block0.out4.*.vtk")}
    stamp = min(times, key=lambda k: abs(times[k] - t))
    (r_f, p_f, z_f), d = assemble(sorted(where.glob(f"liner.block*.out4.{stamp}.vtk")), ("rho", "r2"))
    k = (len(z_f) - 1) // 2
    rr, pp = np.meshgrid(r_f, p_f)
    mesh = ax.pcolormesh(rr * np.cos(pp), rr * np.sin(pp), np.log10(d["rho"][k]), vmin=-4, vmax=1.5,
                         cmap=LinearSegmentedColormap.from_list("blue", BLUE), shading="flat", rasterized=True)
    rc, pc = 0.5 * (r_f[1:] + r_f[:-1]), 0.5 * (p_f[1:] + p_f[:-1])
    pc = np.append(pc, pc[0] + 2 * np.pi)  # close the circle for the contour
    gas = np.vstack([d["r2"][k], d["r2"][k][:1]])
    rr, pp = np.meshgrid(rc, pc)
    if (gas > 0.5).any():
        ax.contour(rr * np.cos(pp), rr * np.sin(pp), gas, levels=[0.5], colors=ORANGE, linewidths=1.0)
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    style(ax, f"{title}, {times[stamp]:.1f} ns")
    return mesh


def colour_bar(fig, mesh, axes, pad: float = 0.07) -> None:
    bar = fig.colorbar(mesh, ax=axes, orientation="horizontal", fraction=0.025, pad=pad, aspect=50)
    bar.set_ticks([-4, -2, 0, 1])
    bar.set_ticklabels(["0.0001", "0.01", "1", "10"])
    bar.set_label("density (g/cm³)", fontsize=8, color=INK)
    bar.ax.tick_params(labelsize=7, colors=MUTED)


def figure_3d() -> list[Path]:
    """The 3-D quick look beside its 2-D twin (r-z slices), and cuts across the 3-D run."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = json.loads((OUT / "results.json").read_text())
    paths = []
    fig, axes = plt.subplots(2, 2, figsize=(11, 4.4), facecolor=SURFACE,
                             gridspec_kw={"width_ratios": [2.5, 0.6], "wspace": 0.12, "hspace": 0.6})
    for row, (name, label) in enumerate([("2-D twin", "2-D twin"), ("3-D quick look", "3-D, slice through the axis")]):
        density_map(axes[row, 0], folder(name), 120.0, 2.5, label)
        mesh = density_map(axes[row, 1], folder(name), res[name]["stagnation"]["time_ns"], 0.6, "stagnation (zoom)")
        axes[row, 0].set_ylabel("z (mm)", fontsize=8, color=INK)
    for ax in axes[-1]:
        ax.set_xlabel("radius (mm)", fontsize=8, color=INK)
    colour_bar(fig, mesh, axes, pad=0.16)
    fig.suptitle("The liner at half resolution: 2-D against 3-D (1 µm ripple; the edge of the D-T gas in orange)",
                 x=0.06, ha="left", fontsize=10, color=INK)
    paths.append(OUT / "density_3d_rz.png")
    fig.savefig(paths[-1], dpi=130, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)

    t_stag = res["3-D quick look"]["stagnation"]["time_ns"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.3), facecolor=SURFACE, gridspec_kw={"wspace": 0.25})
    for ax, (t, r_max) in zip(axes, [(120.0, 2.5), (t_stag - 5.0, 1.0), (t_stag, 0.4)]):
        mesh = cross_section(ax, folder("3-D quick look"), t, r_max, "3-D, across the liner")
        ax.set_xlabel("x (mm)", fontsize=8, color=INK)
    axes[0].set_ylabel("y (mm)", fontsize=8, color=INK)
    colour_bar(fig, mesh, axes, pad=0.14)
    fig.suptitle("The 3-D quick look, cut across at mid-length (the edge of the D-T gas in orange)",
                 x=0.06, ha="left", fontsize=10, color=INK)
    paths.append(OUT / "density_3d_cut.png")
    fig.savefig(paths[-1], dpi=130, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return paths


def figure_2d() -> Path:
    """Density at 120 ns (the same moment for every case) and at each case's stagnation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = json.loads((OUT / "results.json").read_text())
    cases = [("flat", "Smooth liner"), ("ripple 0.1 um", "0.1 µm ripple"), ("ripple 1 um", "1 µm ripple"),
             ("ripple 10 um", "10 µm ripple")]
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 7.6), facecolor=SURFACE,
                             gridspec_kw={"width_ratios": [2.5, 0.6], "wspace": 0.12, "hspace": 0.55})
    for row, (name, label) in enumerate(cases):
        density_map(axes[row, 0], folder(name), 120.0, 2.5, label)
        mesh = density_map(axes[row, 1], folder(name), res[name]["stagnation"]["time_ns"], 0.6,
                           "stagnation (zoom)")
        axes[row, 0].set_ylabel("z (mm)", fontsize=8, color=INK)
    for ax in axes[-1]:
        ax.set_xlabel("radius (mm)", fontsize=8, color=INK)
    colour_bar(fig, mesh, axes)
    fig.suptitle("The thin, ice-lined liner in 2-D (Athena++): density, with the edge of the D-T gas in orange",
                 x=0.06, ha="left", fontsize=10, color=INK)
    path = OUT / "density_2d.png"
    fig.savefig(path, dpi=130, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return path


def summarise(name: str) -> dict:
    where = folder(name)
    h = history(where)
    r_gas = np.sqrt(h["gas_r2"] / h["gas_mass"])
    r_ice = np.sqrt(h["ice_r2"] / h["ice_mass"])
    temp = h["gas_pv"] / h["gas_mass"] * T_UNIT_KEV
    k = int(np.argmin(r_gas))
    stag = {"time_ns": float(h["time"][k]), "gas_rms_um": float(r_gas[k] * 1e3),
            "convergence": float(r_gas[0] / r_gas[k]), "gas_T_keV": float(temp[k]),
            "gas_rho_max_g_cc": float(h["gas_rho_max"][k]), "ice_rms_um": float(r_ice[k] * 1e3),
            "ended_ns": float(h["time"][-1])}
    out = {"name": name, "ripple_um": CASES[name]["ripple"] * 1e3, "stagnation": stag,
           "snapshots": snapshots(where), "azimuthal": azimuthal(where),
           "trace": {"time_ns": h["time"][::20].tolist(), "gas_rms_um": (r_gas[::20] * 1e3).tolist(),
                     "gas_T_keV": temp[::20].tolist()}}
    if "fus_W" in h:
        out["burn"] = burn_summary(name, h)
    return out


LENGTH_MM = 10.0  # the real liner (athinput.liner_physics length_mm): the physics history is for all of it
DT_J_PER_KG = 17.59e6 * 1.602176634e-19 / (5.03 * 1.66053907e-27)  # D-T burnt completely


def electric_gain(yield_j: float, bank_j: float) -> float:
    """Electricity out / in as the 0-D simulator counts it (docs/simulator.md): a lithium blanket, heat to
    electricity, and the bank charged from the grid at the pulsed-power efficiency (the preheat laser left out)."""
    sys.path.insert(0, str(HERE.parent))
    from nanofusion.improved import simulator as S
    st = S.Settings()
    blanket = S.ALT.blanket_energy_per_fusion_mev(0.0, uranium=False) / S.Q_DT_MEV
    return yield_j * blanket * st.thermal_efficiency / (bank_j / st.efficiency("pulsed power"))


def burn_summary(name: str, h: dict) -> dict:
    """The fusion yield and the conditions it came from, from a physics run's history."""
    t = h["time"]

    def integral(y):  # over time, s
        return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t))) * 1e-9

    f = h["fus_W"]
    energy = max(integral(f), 1e-300)
    k = int(np.argmax(f))
    above = t[f >= 0.5 * f[k]] if f[k] > 0 else np.array([t[k]])
    fuel_kg = (h["ice_mass"][0] + h["gas_mass"][0]) * 1e-6 * LENGTH_MM / z_length_mm(name)  # mg -> kg
    yield_j = float(h["yield_J"][-1])
    bank = float(h["bank_J"][0])
    b_z, r_burn = integral(h["fus_Bz"]) / energy, integral(h["fus_r"]) / energy
    return {"yield_J": yield_j, "electric_gain": electric_gain(yield_j, bank), "bank_J": bank,
            "peak_fusion_W": float(f[k]), "bang_ns": float(t[k]), "burn_width_ns": float(above[-1] - above[0]),
            "Ti_keV": integral(h["fus_Ti"]) / energy, "Te_keV": integral(h["fus_Te"]) / energy,
            "Bz_T": b_z, "r_burn_um": r_burn * 1e6, "BR_Tm": b_z * 1.5 * r_burn,  # uniform column: R = 1.5 <r>
            "mixed_share": integral(h["fus_mixed_W"]) / energy, "burnup": yield_j / (fuel_kg * DT_J_PER_KG),
            "peak_current_MA": float(h["current_MA"].max()), "delivered_J": float(h["delivered_J"][-1]),
            "radiated_J": float(h["radiated_J"][-1]), "alpha_made_J": float(h["alpha_made_J"][-1]),
            "alpha_dep_J": float(h["alpha_dep_J"][-1]), "ended_ns": float(t[-1])}


def report() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    done = [n for n in CASES if (folder(n) / "liner.hst").exists()]
    results = {n: summarise(n) for n in done}
    (OUT / "results.json").write_text(json.dumps(results, indent=1))
    lines = ["# The thin, ice-lined liner in 2-D and 3-D (Athena++)", "",
             ("Produced by `tools/athena_liner.py` from Athena++ runs of `tools/athena/liner.cpp`. The design and "
              "the 0-D model's prediction (~177× convergence) are in [docs/liner_experiment.md]"
              "(../../liner_experiment.md). Convergence is the gas's starting rms radius over its smallest."), "",
             ("| Case | Liner ripple (rms) | Fuel convergence at stagnation | When | Gas temperature then | "
              "Liner's inner surface spread (r-z slice nearest stagnation) | Beryllium in the hot gas, then |"),
             "|---|---|---|---|---|---|---|"]
    for n, r in results.items():
        s = r["stagnation"]
        near = min(r["snapshots"], key=lambda x: abs(x["time_ns"] - s["time_ns"])) if r["snapshots"] else None
        lines.append(f"| {n} | {r['ripple_um']:g} µm | {s['convergence']:.0f}× ({s['gas_rms_um']:.0f} µm rms) | "
                     f"{s['time_ns']:.1f} ns | {s['gas_T_keV']:.2g} keV | "
                     + (f"{near['inner_spread_mm'] * 1e3:.0f} µm at {near['time_ns']:.0f} ns | "
                        + ("no mostly-gas cell left in the slice |" if near["mixed_share"] is None
                           else f"{near['mixed_share']:.0%} of the gas mass |") if near else "– | – |"))
    burns = [r for r in results.values() if "burn" in r]
    if burns:
        lines += ["", "## With the physics a yield needs", "",
                  ("Burn, alpha heating, radiation, resistivity, conduction, two temperatures, the Nernst effect, "
                   "the bank circuit and a real equation of state (tools/athena/liner.cpp, physics build); each "
                   "\"no X\" case switches one off. Energies are for the whole 1 cm liner. Electric gain is counted "
                   "as the 0-D simulator counts it (lithium blanket, 40% heat to electricity, the bank charged at 15%): "
                   "the 0-D model's figure for this design is 22."), "",
                  "| Case | Yield | Electric gain | Fuel burnt | Peak fusion power | Bang time | Burn width | "
                  "Convergence | Ended |", "|---|---|---|---|---|---|---|---|---|"]
        for r in burns:
            b, s = r["burn"], r["stagnation"]
            lines.append(f"| {r['name']} | {b['yield_J']:.3g} J | {b['electric_gain']:.3g} | {b['burnup']:.2%} | "
                         f"{b['peak_fusion_W']:.3g} W | {b['bang_ns']:.1f} ns | {b['burn_width_ns']:.2f} ns | "
                         f"{s['convergence']:.0f}× | {b['ended_ns']:.0f} ns |")
        lines += ["", "| Case | Burn-weighted T_i / T_e | B_z, B R | Burn radius | Burn in Be-mixed cells | "
                  "Peak current | Into the load | Radiated | Alpha energy deposited / made |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for r in burns:
            b = r["burn"]
            share = b["alpha_dep_J"] / b["alpha_made_J"] if b["alpha_made_J"] > 0 else float("nan")
            lines.append(f"| {r['name']} | {b['Ti_keV']:.2f} / {b['Te_keV']:.2f} keV | {b['Bz_T']:.3g} T, "
                         f"{b['BR_Tm']:.2g} T m | {b['r_burn_um']:.0f} µm | {b['mixed_share']:.0%} | "
                         f"{b['peak_current_MA']:.1f} MA | {b['delivered_J']:.3g} J | {b['radiated_J']:.3g} J | "
                         f"{share:.0%} |")
    three = [r for r in results.values() if r["azimuthal"]]
    for r in three:
        t_stag = r["stagnation"]["time_ns"]
        lines += ["", f"## {r['name']}: departures from axisymmetry", "",
                  ("| Time | Liner inner surface (mean) | Its rms ripple: axisymmetric | non-axisymmetric | "
                   "Strongest non-axisymmetric mode (m, n) | Fuel rms radius | Fuel column off axis (max) |"),
                  "|---|---|---|---|---|---|---|"]
        for a in r["azimuthal"]:
            if abs(a["time_ns"] - t_stag) > 6 and round(a["time_ns"]) % 20:
                continue  # every 20 ns, and every snapshot near stagnation
            if a.get("liner_at_axis_share", 0.0) > 0.5:
                surface = f"beryllium on the axis ({a['liner_at_axis_share']:.0%} of the length and circle) | – | – | –"
            else:
                surface = (f"{a['inner_mean_um']:.0f} µm | {a['inner_axisym_rms_um']:.1f} µm | "
                           f"{a['inner_nonaxisym_rms_um']:.1f} µm | {tuple(a['top_mode_m_n'])}")
            lines.append(f"| {a['time_ns']:.1f} ns | {surface} | {a['fuel_rms_um']:.0f} µm | "
                         f"{a['fuel_offset_max_um']:.1f} µm |")
    figures = []
    if all(n in results for n in ("flat", "ripple 0.1 um", "ripple 1 um", "ripple 10 um")):
        figures.append(figure_2d())
    if all(n in results for n in ("3-D quick look", "2-D twin")):
        figures += figure_3d()
    lines += [""] + [f"![{f.stem.replace('_', ' ')}]({f.name})" for f in figures]
    lines.append("")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run(sys.argv[2:] or DEFAULT_RUN)
        return 0
    if len(sys.argv) > 1 and sys.argv[1] == "build-physics":
        return build_physics()
    if len(sys.argv) > 2 and sys.argv[1] == "job":
        print(shell_job(sys.argv[2]), end="")
        return 0
    return report()


if __name__ == "__main__":
    sys.exit(main())
