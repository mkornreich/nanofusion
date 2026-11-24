"""Checks on the improved physics model (small crystals, so they run fast)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from nanofusion.improved import ImprovedProperties, build_crystal, native, run_improved
from nanofusion.improved import physics as P


def test_crystal_is_rock_salt_lid():
    pos, is_li = build_crystal(ImprovedProperties(shape="cube", size=4))
    assert len(pos) == 64 and is_li.sum() == 32
    d = np.linalg.norm(pos[:, None] - pos[None], axis=2)
    np.fill_diagonal(d, np.inf)
    nearest = d.argmin(axis=1)
    assert np.allclose(d.min(axis=1), P.SPACING)
    assert np.all(is_li != is_li[nearest])  # every nearest neighbour is the other species


def test_coulomb_explosion_conserves_energy_and_matches_theory():
    _, rep = run_improved(ImprovedProperties(size=0.7, intensity_w_cm2=0, duration_fs=60, frames=5))
    d = rep.dynamics
    total_ke = rep.deuterons * d["mean_D_keV"] + (rep.particles - rep.deuterons) * d["mean_Li_keV"]
    assert abs(d["energy_check_ion_KE_minus_potential_drop_keV"]) < 1e-5 * total_ke
    assert 0.6 < d["max_D_keV"] / d["coulomb_explosion_theory_max_D_keV"] < 1.1
    assert rep.fusion["total"] < 1e-15


def test_pulse_has_no_net_push():
    pulse = P.Pulse(e0=1e12, wavelength=3.6e-6, fwhm=6e-15, phase=math.radians(37))
    t = np.linspace(-60e-15, 60e-15, 200001)
    area = np.trapezoid([pulse.field(x) for x in t], t)
    assert abs(area) < 1e-6 * pulse.e0 * pulse.fwhm


def test_laser_energy_is_accounted_for():
    # Stripping at a field zero (phase 90 degrees) leaves the ions a net drift,
    # so the laser does real work; ion energy gain must equal it.
    _, rep = run_improved(ImprovedProperties(size=0.6, intensity_w_cm2=1e18, phase_deg=90,
                                             duration_fs=30, frames=5))
    d = rep.dynamics
    assert d["laser_work_on_ions_keV"] > 1.0
    assert d["energy_check_ion_KE_minus_potential_drop_keV"] == pytest.approx(d["laser_work_on_ions_keV"],
                                                                              rel=1e-3)


def test_cross_sections_against_known_values():
    # Bosch-Hale: total D+D near 70 mb at 100 keV CM; D+6Li Gamow constant ~115 keV^0.5
    assert P.sigma_dd_mb(100.0) == pytest.approx(70, rel=0.05)
    assert P.gamow_constant(1, 1, P.reduced_mass_u(P.MASS_U["D"], P.MASS_U["D"])) == pytest.approx(31.397, rel=1e-3)
    assert P.gamow_constant(1, 3, P.reduced_mass_u(P.MASS_U["D"], P.MASS_U["Li6"])) == pytest.approx(115.3, rel=1e-2)
    # Thick-target D+D yield of a 100 keV deuteron in LiD: of order 1e-7
    _, dd = P.thick_target_yield([100.0], "Li6")
    assert 1e-8 < dd[0] < 1e-6


def test_numpy_fallback_matches_kernel():
    if not native.native_available():
        pytest.skip("no C compiler")
    pos, is_li = build_crystal(ImprovedProperties(size=0.6))
    rng = np.random.default_rng(1)
    vel = rng.normal(scale=3e5, size=pos.shape)
    kind = is_li.astype(np.int32)
    mass = np.where(is_li, 7.016, 2.014) * P.AMU
    q = np.where(is_li, 3.0, 1.0) * P.E_CHARGE
    args = (pos, vel, q, mass, kind, 1e-24, (4.08e-10) ** 2, 1e27, P.S_DLI_KEV_MB["Li6"], 115.3)
    f1, s1, hd1, hl1 = native.coulomb_step(*args)
    f2, s2, hd2, hl2 = native._coulomb_step_numpy(*args)
    np.testing.assert_allclose(f1, f2, rtol=1e-9, atol=1e-20)
    np.testing.assert_allclose(s1, s2, rtol=1e-6)
    np.testing.assert_allclose(hd1, hd2, rtol=1e-9, atol=1e-6)
    np.testing.assert_allclose(hl1, hl2, rtol=1e-9, atol=1e-6)
    assert hl1.sum() > 0


def test_expanded_profile_and_pair_geometry():
    base = ImprovedProperties(size=1.0)
    pos_u, li_u = build_crystal(base)
    pos_g, li_g = build_crystal(ImprovedProperties(size=1.0, profile="expanded", core_density=0.5, expansion=3))
    assert len(pos_g) == len(pos_u) and li_g.sum() == li_u.sum()
    d = np.linalg.norm(pos_g[:, None] - pos_g[None], axis=2)
    np.fill_diagonal(d, np.inf)
    assert d.min() >= P.SPACING * 0.999  # never compressed below the solid spacing
    assert np.linalg.norm(pos_g, axis=1).max() > 1.5 * np.linalg.norm(pos_u, axis=1).max()
    pos_p, li_p = build_crystal(ImprovedProperties(size=1.0, geometry="pair", gap_nm=2.0))
    assert len(pos_p) == 2 * len(pos_u)
    assert (pos_p[:, 0] < 0).sum() == len(pos_u)


def test_dt_cross_section_and_catchers():
    assert P.sigma_dt_mb(64.0) == pytest.approx(5000, rel=0.05)  # D-T peak ~5 b near 64 keV
    y = P.beam_target([100.0], "D", "LiT", "Li6")
    assert 3e-6 < y["D+T"][0] < 1e-4  # neutron-generator scale: ~1e-5 per 100 keV deuteron
    assert P.beam_target([1000.0], "Li", "LiD", "Li6")["Li+D"][0] > 0


def test_pair_exposure_reproduces_in_crystal_rate():
    # The exposure histogram times sigma at each bin should reproduce the in-crystal D+D yield
    # computed directly (they differ only by binning and by the small Coulomb energy term).
    _, rep = run_improved(ImprovedProperties(size=0.8, lithium="Li7", intensity_w_cm2=1e21, phase_deg=90,
                                             duration_fs=20, frames=3))
    edges = np.array(rep.pair_exposure["edges_keV"])
    centres = np.sqrt(edges[1:] * edges[:-1])
    mu = P.reduced_mass_u(P.MASS_U["D"], P.MASS_U["Li7"])
    bg = P.gamow_constant(1, 3, mu)
    from_hist = np.sum(P.sigma_dli_mb(centres, P.S_DLI_KEV_MB["Li7"], bg) * P.MB_TO_M2
                       * np.array(rep.pair_exposure["D_Li"]))
    direct = rep.fusion["in_crystal_D_Li"]
    assert direct > 0 and 1e-3 < from_hist / direct < 1e3


def test_report_has_size_table():
    _, rep = run_improved(ImprovedProperties(size=0.6, intensity_w_cm2=1e19, duration_fs=10, frames=3))
    radii = [row["radius_nm"] for row in rep.size_table]
    assert 15.0 in radii and radii == sorted(radii)
    big = rep.size_table[-1]
    assert big["LiD"]["total"] > rep.size_table[0]["LiD"]["total"]  # bigger crystals fuse more


def test_energy_model_against_textbook_values():
    from nanofusion.improved import energy as E

    # Bosch-Hale D-T reactivity, and the Q = 1 (scientific breakeven) triple product at 15 keV
    assert E.reactivity("D", "T", 10.0) == pytest.approx(1.136e-22, rel=0.03)
    assert E.reactivity("D", "T", 20.0) == pytest.approx(4.33e-22, rel=0.03)
    triple = 12 * 15.0**2 / (E.reactivity("D", "T", 15.0) * 17590)
    assert 4e20 < triple < 8e20
    # Spitzer/Stix slowing time of 2 u ions at T_e = 10 keV, n_e = 1e20 m^-3: ~0.74 s
    p = E.Plasma({"D": 5e19, "T": 5e19}, t_e=10, t_i=10)
    assert p.slowing_time("D") == pytest.approx(0.74, rel=0.05)
    # Beam into hot D-T plasma multiplies energy about 1-3x (neutral-beam-driven fusion); cold matter far less
    assert 0.5 < E.beam_in_plasma(200, "D", p)["gain"] < 3
    cold = E.beam_in_plasma(200, "D", E.Plasma({"D": 5e19, "T": 5e19}, t_e=0.3, t_i=0.3))["gain"]
    assert cold < 0.2


def test_crystal_ball_burn_scales_like_inertial_fusion():
    from nanofusion.improved import energy as E

    small = E.array_burn("DT", 25, 1.0, 1e-4)
    big = E.array_burn("DT", 25, 1.0, 3e-2)
    assert small["gain_after_radiation"] < 0.5 < 1 < big["gain_after_radiation"]
    # Lithium-7 deuteride radiates more than it fuses at any temperature reached here
    assert E.array_burn("7LiD", 50, 1.0, 1e-2)["fusion_to_brem"] < 1


def test_proton_boron_cross_section_and_avalanche():
    from nanofusion.improved import energy as E

    # Resonances: ~0.1 b at 148 keV and ~1.1-1.2 b near 612 keV (CM); the fit sections join smoothly
    assert P.sigma_pb11_mb(148.0) == pytest.approx(100, rel=0.2)
    assert 900 < P.sigma_pb11_mb(612.0) < 1300
    for edge in (400.0, 700.0):
        lo, hi = P.s_factor_pb11_mev_b(edge - 1e-6), P.s_factor_pb11_mev_b(edge + 1e-6)
        assert lo == pytest.approx(hi, rel=0.01)
    # Reactivity of order 1e-22..1e-21 m^3/s at a few hundred keV
    assert 1e-22 < E.reactivity("p", "B11", 300.0) < 1e-21
    # No chain reaction: knock-on protons trigger far less than one fusion per fusion
    assert E.pb11_avalanche("decaborane", nuclear_boost=5)["k"] < 1e-4
    hot = E.Plasma({"p": 6.4e28, "B11": 4.6e28}, t_e=100, t_i=100)
    assert E.pb11_avalanche(hot_plasma=hot, nuclear_boost=5)["k"] < 0.1


def test_screening_raises_low_energy_yields_only():
    from nanofusion.improved import energy as E

    n = P.lid_molecule_density()
    bare = E.thick_target([3.0, 300.0], "D", "LiD", {"T": n})["reactions"]
    screened = E.thick_target([3.0, 300.0], "D", "LiD", {"T": n}, screening_ev=1000)["reactions"]
    assert screened[0] > 10 * bare[0]  # big boost for slow ions
    assert screened[1] == pytest.approx(bare[1], rel=0.05)  # negligible where the fusion actually happens


def test_reproduces_texas_petawatt_cluster_record():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location("energy_experiments",
                                                  Path(__file__).parents[1] / "tools" / "energy_experiments.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    lines, data = [], {}
    mod.reproduce_cluster_experiment(lines, data)
    predicted = [row["neutrons"] for row in data["texas_petawatt"]["rows"]]
    assert min(predicted) < 1.9e7 < max(predicted)  # measured record lies inside the published input ranges
    assert 1e-7 < data["texas_petawatt"]["energy_ratio"] < 1e-6


def test_hollow_shell_implosion_focuses_and_obeys_space_charge():
    from nanofusion.improved.implosion import (
        ImplosionProperties,
        build_shell,
        run_implosion,
    )

    props = ImplosionProperties(cavity_nm=3, electron_ratio=8, fuel="DT", reaction_radius_nm=0.05)
    shell = build_shell(props)
    charge = len(shell["pos"])  # D and T both carry one charge
    assert abs(shell["electrons"] / charge - 8) < 0.1
    assert np.all(np.linalg.norm(shell["pos"], axis=1) >= 3e-9)
    rep = run_implosion(props)
    assert abs(rep["dynamics"]["energy_drift_relative"]) < 1e-4
    th = rep["theory"]
    foci = {f["species"]: f for f in rep["focus"] if f["pass"] == 1}
    # Every ion of a species reaches the centre together, near pi / (2 omega); tritium later than deuterium.
    for s in ("D", "T"):
        assert abs(foci[s]["time_fs"] / th["focus_time_fs"][s] - 1) < 0.25  # space charge delays the bounce a little
    ratio = math.sqrt(P.MASS_U["T"] / P.MASS_U["D"])
    assert th["focus_time_fs"]["T"] / th["focus_time_fs"]["D"] == pytest.approx(ratio)
    assert foci["T"]["time_fs"] > 1.15 * foci["D"]["time_fs"]  # so D and T miss each other at the focus
    col = rep["fusion"]["column_per_ion_m2"]
    assert col["D sees T"] < 0.2 * col["D sees D"]
    # The ions' own repulsion stops the focus near where a smooth shell of the same charge would bounce.
    assert 0.5 < foci["D"]["median_radius_over_R"] / th["thin_shell_bounce_radius_over_R"] < 3
    # The matter each ion crosses stays below the space-charge bound 4 eps0 E / (e^2 r_c).
    e = rep["dynamics"]["mean_ion_KE_at_max_keV"] * P.J_PER_KEV
    bound = 4 * P.EPS0 * e / (P.E_CHARGE**2 * foci["D"]["median_radius_nm"] * 1e-9)
    assert 0.1 < rep["fusion"]["column_per_ion_m2"]["D sees D"] / bound < 2
    assert rep["energy"]["gain_vs_field_energy"] < 1e-6


def test_disk_explodes_sideways_and_dot_target_beams():
    # A single-layer LiD disk: its mid-plane is a symmetry plane, so every ion stays in it (a flat fan).
    _, rep = run_improved(ImprovedProperties(shape="disk", size=2.0, thickness_nm=0.2, intensity_w_cm2=0,
                                             duration_fs=40, frames=3))
    v = rep.final_state["velocities"]
    assert np.all(np.abs(v[:, 2]) < 1e-6 * np.linalg.norm(v, axis=1).max())
    # Double layer: lithium-only disk with a deuterium dot on its +z face. The dot is pushed forward as a beam.
    props = ImprovedProperties(shape="disk", size=3.0, thickness_nm=0.6, dot_nm=1.2, lithium="Li7",
                               intensity_w_cm2=0, duration_fs=60, frames=3)
    pos, is_li = build_crystal(props)
    assert pos[~is_li, 2].min() > pos[is_li, 2].max()  # the dot sits on top of the disk
    _, rep = run_improved(props)
    v = rep.final_state["velocities"][~rep.final_state["is_lithium"]]
    angle = np.degrees(np.arccos(v[:, 2] / np.linalg.norm(v, axis=1)))
    assert np.mean(angle < 30) > 0.8
    assert abs(rep.dynamics["energy_check_ion_KE_minus_potential_drop_keV"]) < 1e-4 * rep.particles


def test_muon_catalysis_and_lattice_confinement_models():
    from nanofusion.improved import alternatives as A

    # The cycle formula gives the Los Alamos record (150 per muon) at 2.1e8 /s with 0.45% sticking,
    # and needs 284 fusions to repay a 5 GeV muon (the published breakeven figure).
    assert A.fusions_per_muon(2.1e8, 0.0045) == pytest.approx(150, rel=0.01)
    assert A.fusions_needed(5.0) == pytest.approx(284, rel=0.01)
    assert A.fusions_per_muon(1e15, 0.0045) == pytest.approx(1 / 0.0045, rel=1e-3)  # sticking caps it
    # Photoneutrons from a 2.9 MeV bremsstrahlung endpoint: mean 145 keV, knocked-on deuterons 64 keV (NASA's values)
    e = A.photoneutron_energies()
    assert e.mean() == pytest.approx(145, rel=0.03)
    res = A.lattice_fusions_per_neutron(e, "ErD3", 8e28, 0.0, collisions=1.0)
    assert res["mean_deuteron_keV"] == pytest.approx(64, rel=0.05)
    screened = A.lattice_fusions_per_neutron(e, "ErD3", 8e28, 3000.0, collisions=1.0)
    assert res["fusions_per_neutron"] < screened["fusions_per_neutron"] < 1.3 * res["fusions_per_neutron"]
    assert screened["fusions_per_neutron"] < 1e-6


def test_cold_squeeze_muon_factory_and_polarisation():
    from nanofusion.improved import alternatives as A
    from nanofusion.improved import energy as E

    # Degenerate electrons: the relativistic formula reduces to 3/5 E_F at low density.
    n = 1e30
    e_f = (A.HBAR_C_MEV_FM * 1e-15) ** 2 * (3 * math.pi**2 * n) ** (2 / 3) / (2 * A.M_E_MEV) * 1e3
    assert A.degenerate_electron_energy_kev(n) == pytest.approx(0.6 * e_f, rel=1e-3)
    # Pycnonuclear burning switches on steeply between NIF-class and white-dwarf densities.
    rates = [A.pycnonuclear_rate_per_nucleus(rho, "DT") for rho in (1e3, 1e5, 1e7)]
    assert rates[0] < 1e-20 < rates[1] < 1 < 1e9 < rates[2]
    assert A.cold_compression(1e7)["ignition_radius_m"] < 1e-5 < A.cold_compression(1e6)["ignition_radius_m"]
    # A crystal muon factory can't beat the accelerator cost it is built on.
    m = A.crystal_muon_factory(4.0, 4.7, 1.0)
    assert m["laser_GeV_per_muon"] >= 4.7
    assert A.crystal_muon_factory(0.2, 4.7, 1.0)["laser_GeV_per_muon"] == math.inf  # below the pion threshold
    # Spin polarisation scales only the D-T cross section, and is off by default.
    base = float(E.sigma_m2("D", "T", 64.0))
    try:
        E.DT_SPIN_FACTOR = 1.5
        assert float(E.sigma_m2("D", "T", 64.0)) == pytest.approx(1.5 * base)
        assert float(E.sigma_m2("D", "D", 64.0)) == pytest.approx(float(P.sigma_dd_mb(64.0)) * P.MB_TO_M2)
    finally:
        E.DT_SPIN_FACTOR = 1.0


def test_dense_matter_rates_cover_all_regimes():
    from nanofusion.improved import alternatives as A
    from nanofusion.improved import energy as E

    # Hot, dilute limit: the ordinary Maxwellian D-T rate (per nucleus n <sigma v> / 4).
    n = 1e3 / (2.515 * 1.66053906660e-27)
    ratio = A.dense_fusion_rate(1.0, 1e4, "DT") / (n * E.reactivity("D", "T", 10.0) / 4)
    assert 0.8 < ratio < 1.2
    # Zero temperature: Salpeter & Van Horn's pycnonuclear rate.
    assert A.dense_fusion_rate(1e7, 0.0) == pytest.approx(A.pycnonuclear_rate_per_nucleus(1e7, "DT"), rel=0.1)
    # Liquid metallic hydrogen near freezing: negligible in every model.
    for model in A.DENSE_MODELS:
        assert A.dense_fusion_rate(40.0, 0.082, "DT", model) < 1e-50
    # Hot fuel ignites more cheaply than cold fuel at the same density.
    hot = A.ignition_ball(1e3, 16e3)["energy_J"]
    cold = A.ignition_ball(1e3, 1.0)["energy_J"]
    assert hot < 1e4 and cold > 1e20 * hot


def test_tamped_magnetised_ball():
    from nanofusion.improved import confinement as C
    from nanofusion.improved import energy as E

    # With no tamper, field or self-heating it reduces to the fixed-temperature burn of energy.array_burn.
    plasma = E.fuel_plasma("DT", 1.0, 20.0, 20.0)
    heat = 1.5 * (plasma.n_e + sum(plasma.ions.values())) * 20.0 * E.KEV
    old = plasma.fusion_power()[0] * 1e-3 / (4 * plasma.sound_speed()) / heat
    bare = C.burn(C.Ball(radius_m=1e-3, t_i_kev=20.0, alpha_heating=False))
    assert bare["gain"] == pytest.approx(old, rel=0.03)
    # Pieces: alpha gyro-radius 0.27 m in 1 T, Braginskii's unmagnetised and strongly magnetised limits,
    # Krokhin-Rozanov deposition joining smoothly, a gold wall receding at about a tenth of the sound speed.
    assert C.ALPHA_LARMOR_TM == pytest.approx(0.27, rel=0.02)
    assert C.kappa_perp_electron(0.0) == pytest.approx(3.16, rel=0.01)
    assert C.kappa_perp_ion(0.0) == pytest.approx(3.906, rel=0.01)
    assert C.kappa_perp_ion(30.0) == pytest.approx(2 / 30**2, rel=0.02)
    assert C.sphere_deposition(0.5 - 1e-9) == pytest.approx(C.sphere_deposition(0.5 + 1e-9), rel=1e-6)
    assert C.alpha_range_factor(0.0) == 1.0 and 0.6 < C.alpha_range_factor(0.3) < 0.7
    assert C.magnetic_thickness(1e-3, 0.0) == 0.0 < C.magnetic_thickness(1e-3, 3e3)
    pressure = 2 * C.SOLID_DT_IONS * 20.0 * C.KEV
    assert 0.05 < C.piston_speed(pressure) / plasma.sound_speed() < 0.12
    # Each trick helps, and tamper plus field bring a 0.3 mm ball to about breakeven.
    kw = dict(radius_m=3e-4, t_i_kev=20.0)
    none = C.burn(C.Ball(**kw))
    tamped = C.burn(C.Ball(tamper="gold", **kw))
    both = C.burn(C.Ball(tamper="gold", b_tesla=3e3, **kw))
    assert none["gain"] < tamped["gain"] < both["gain"]
    assert both["gain"] > 10 * none["gain"] and 0.8 < both["gain"] < 1.5
    # Energy bookkeeping: heat left + losses + work = heat at the start + alpha heating.
    spent = sum(both[k] for k in ("heat_left_J", "bremsstrahlung_J", "conduction_J", "work_on_wall_J", "carried_off_J"))
    assert spent == pytest.approx(both["heat_at_start_J"] + both["alpha_heating_J"], rel=1e-6)
    # A thin shell or an entrance hole holds less; ions heated first burn more for the same energy.
    assert C.burn(C.Ball(tamper="gold", b_tesla=3e3, tamper_mm=0.05, **kw))["gain"] < 0.8 * both["gain"]
    assert C.burn(C.Ball(tamper="gold", b_tesla=3e3, hole_fraction=0.05, **kw))["gain"] < 0.7 * both["gain"]
    hot_ions = C.burn(C.Ball(radius_m=3e-4, t_i_kev=39.0, t_e_kev=1.0, tamper="gold", b_tesla=3e3))
    assert hot_ions["energy_in_J"] == pytest.approx(both["energy_in_J"], rel=0.01)
    assert hot_ions["gain"] > 1.3 * both["gain"]
    # A column in a straight field also loses heat along the field to its ends; a ring field does not.
    column = dict(radius_m=3e-4, t_i_kev=20.0, tamper="gold", b_tesla=3e3, aspect=4.0)
    ring, straight = C.burn(C.Ball(**column)), C.burn(C.Ball(field="axial", **column))
    assert straight["gain"] < 0.8 * ring["gain"] and straight["conduction_J"] > 2 * ring["conduction_J"]
    assert C.Ball(**column).energy_in_J == pytest.approx(6 * both["energy_in_J"], rel=0.01)  # 2 a pi r^3 against 4/3 pi r^3
    with pytest.raises(ValueError):
        C.Ball(radius_m=3e-4, field="axial")
    # A large enough ball runs away (ignites); the cheapest breakeven is tens of kilojoules, not megajoules.
    big = C.burn(C.Ball(radius_m=4e-3, t_i_kev=7.0, tamper="gold", b_tesla=1e4))
    assert big["ignited"] and big["gain"] > 50 and big["stopped"] == "grew to the limit"
    assert both["stopped"] == "burn died out" and none["stopped"] == "flew apart"
    # The result does not depend on the integrator's tolerance, and slow burns in thin fuel are followed to the end.
    assert C.burn(C.Ball(radius_m=4e-3, t_i_kev=7.0, tamper="gold", b_tesla=1e4), rtol=1e-7)["gain"] == \
        pytest.approx(big["gain"], rel=2e-3)
    thin = C.burn(C.Ball(radius_m=3e-2, t_i_kev=20.0, fill=1e-4, tamper="gold", b_tesla=100.0))
    assert thin["gain"] > 8 and thin["stopped"] == "burn died out"
    # A gold shell of stated thickness around a column weighs what its geometry says: thicker shells hold longer.
    shells = [C.burn(C.Ball(radius_m=1e-3, t_i_kev=10.0, tamper="gold", b_tesla=3e3, aspect=8.0, tamper_mm=t))["gain"]
              for t in (0.1, 0.3, 1.0)]
    assert shells[0] < shells[1] < shells[2]
    best = C.smallest_ball(1.0, temperatures=(20.0,), tamper="gold", b_tesla=3e3)
    assert 2e4 < best["energy_in_J"] < 1e5
    assert C.smallest_ball(1.0, temperatures=(20.0,), r_hi=1e-3, tamper="none") is None


def test_ball_command_line(capsys):
    from nanofusion.__main__ import main

    assert main(["ball", "--radius", "0.3", "--temperature", "20", "--field", "3000"]) == 0
    assert "gain" in capsys.readouterr().out
    assert main(["ball", "--axial-field"]) == 2  # an axial field needs a column
    assert main(["ball", "--fill", "2"]) == 2


def test_chain_reaction_blanket_trap_and_projectiles():
    from nanofusion.improved import alternatives as A
    from nanofusion.improved import energy as E

    # Bosch-Hale's high-energy D-T fit joins the low-energy one and then falls.
    assert float(P.sigma_dt_mb(550.1, high_energy=True)) == pytest.approx(float(P.sigma_dt_mb(549.9)), rel=0.02)
    assert float(P.sigma_dt_mb(3000.0, high_energy=True)) < 0.4 * float(P.sigma_dt_mb(3000.0))
    assert float(P.sigma_dt_mb(64.0, high_energy=True)) == float(P.sigma_dt_mb(64.0))
    # No chain reaction in cold D-T. In hot D-T it nears 1 only in an infinite medium.
    assert E.dt_chain(neutrons=3000)["k"] < 0.01
    hot_plasma = E.fuel_plasma("DT", 1.0, 30.0, 30.0)
    hot = E.dt_chain(hot_plasma=hot_plasma, nuclear_boost=5.0, neutrons=3000)
    assert 0.5 < hot["k"] < 1.5 and hot["k_neutron"] > hot["k_alpha"]
    thin = E.dt_chain(hot_plasma=hot_plasma, nuclear_boost=5.0, first_collision_probability=0.01, neutrons=3000)
    assert thin["k_neutron"] < 0.01
    assert not E.DT_HIGH_ENERGY
    # Fission blanket: 22 MeV per fusion with lithium only. Depleted uranium with no chain gives about what a
    # natural-uranium pile measured per 14 MeV neutron, 290 MeV (McEnhill & Weale 1962).
    assert A.blanket_energy_per_fusion_mev(0.0, uranium=False) == pytest.approx(22.4, rel=0.01)
    assert A.blanket_energy_per_fusion_mev(0.0) - A.DT_FUSION_MEV + 14.07 == pytest.approx(290, rel=0.2)
    assert A.hybrid_k_needed(5e-4, 0.1) > 0.99  # exploding crystals: a critical reactor in all but name
    assert A.hybrid_k_needed(0.53, 0.4) == 0.0  # muon catalysis: plain uranium is enough for net electricity
    assert 0.6 < A.hybrid_k_needed(0.53, 0.4, plant_gain=4.0) < 0.8
    assert A.hybrid_k_needed(100.0, 0.1) is None  # a good enough fusion source needs no uranium
    k = A.hybrid_k_needed(0.15, 0.1)
    assert A.hybrid_plant_gain(0.15, 0.1, k) == pytest.approx(1.0)
    assert A.fission_share(k) > 0.9
    # Bare ions: the Brillouin limit makes a 10 T trap a milliwatt machine, and power goes as B^4.
    assert A.brillouin_density(10.0) == pytest.approx(1.06e17, rel=0.02)
    trap = A.ion_trap(10.0, 0.05, 1.0)
    assert trap["power_W"] < 1 and trap["voltage_V"] > 1e6
    assert A.ion_trap(20.0, 0.05, 1.0)["power_W"] == pytest.approx(16 * trap["power_W"])
    # Cluster impact: a 300 keV (D2O)100 cluster leaves 0.3 keV per deuteron, far too little to fuse. The yield
    # agrees with Carraro et al.'s estimate of 4e-39 per cluster (Phys. Rev. A 42, 1379 (1990)).
    cluster = A.cluster_impact_yield(100, 300.0)
    assert cluster["keV_per_deuteron"] == pytest.approx(0.3, rel=0.02)
    assert 4e-40 < cluster["fusions_per_cluster"] < 4e-38
    shot = A.crystal_projectile(15.0, 100.0)
    assert 4e7 < shot["voltage_V"] < 7e7 and shot["speed_m_s"] == pytest.approx(3.1e6, rel=0.02)
    # Nanowire arrays (Curtis et al., Phys. Rev. Research 3, 043181 (2021)): the simulated deuteron spectrum,
    # 3.8e12 above 0.2 MeV and 3e11 above 1 MeV, stopping in cold CD2 gives the measured 1.2e7 neutrons.
    slope = 800.0 / math.log(3.8e12 / 3e11)
    spectrum = A.spectrum_yield(slope, 13000.0, e_min_kev=200.0)
    assert 0.5 * 3.8e12 * spectrum["reactions_per_ion"] == pytest.approx(1.2e7, rel=0.5)
    assert 400 < spectrum["mean_keV"] < 600


def test_pauli_blocked_stopping_and_degenerate_routes():
    from nanofusion.improved import degenerate as D

    # The exact RPA low-speed coefficient matches Lindhard & Winther's Table 1, and the stopping number goes to
    # their friction law for slow ions and to Bethe's formula for fast ones.
    assert D.low_speed_coefficient(0.1) == pytest.approx(0.787, rel=0.005)
    assert D.low_speed_coefficient(0.01) == pytest.approx(1.823, rel=0.005)
    gas = D.fermi_gas(4.5e31)
    slow = float(D.stopping_number(0.1, gas["chi2"]))
    assert slow == pytest.approx(D.low_speed_coefficient(gas["chi2"]) * 1e-3, rel=0.02)
    bethe = math.log(2 * P.M_E * (20 * gas["v_F"]) ** 2 / (D.HBAR * gas["omega_p"]))
    assert float(D.stopping_number(20.0, gas["chi2"])) == pytest.approx(bethe, rel=0.01)
    # Slow protons in aluminium: linear response gives two-thirds of the measured stopping (NIST PSTAR, 10 keV).
    n_al = 3 * 2.699e3 / (26.9815 * P.AMU)
    al = float(D.electron_stopping(10.0, P.MASS_U["p"], 1, n_al)) * 1e-5 / 2.699
    assert 0.55 < al / 292.1 < 0.8
    # Squeezed cold D-T: a deuteron returns ~0.4 of its energy at 250 g/cm3, more than all of it above
    # ~1,000 g/cm3, and about 5 at most, once collisions with ions take over. Protons in squeezed boron stay below 1.
    assert 0.3 < D.best_in_flight("D", "DT", 250.0)["gain"] < 0.5
    assert D.best_in_flight("D", "DT", 2500.0)["gain"] > 1.5
    top = D.best_in_flight("D", "DT", 2.5e5)
    assert 4 < top["gain"] < 6 and top["electron_share"] < 0.1
    assert D.best_in_flight("p", "boron", 2.34e5, e_lo=100.0, e_hi=8000.0)["gain"] < 1
    D.ELECTRON_DRAG_FACTOR = 1.6
    try:
        assert D.best_in_flight("D", "DT", 1000.0)["gain"] < 0.8
    finally:
        D.ELECTRON_DRAG_FACTOR = 1.0
    # Son & Fisch: 100 keV ions in proton-boron keep themselves hot above ~1e5 g/cm3, with electrons near 30 keV.
    rho = D.self_heating_density("pB", 100.0)
    assert 1e5 < rho < 1e6
    balance = D.hot_ion_balance("pB", rho, 100.0)
    assert 20 < balance["T_e_keV"] < 45 and balance["T_e_over_T_F"] < 1
    assert balance["self_heating"] == pytest.approx(1.0, rel=0.05)
    # Muons: tuned to 35% reactivation at liquid density; squeezing to 100 g/cm3 frees over 90% of stuck muons,
    # mostly while the stuck ion is faster than the muon's orbit. A 50 um core at that density lasts ~0.1 ns.
    assert D.MU_HE_BINDING_KEV == pytest.approx(10.94, rel=0.002)
    assert float(D.hydrogen_ionisation_m2(100.0)) * 1e4 == pytest.approx(1.12e-16, rel=0.05)  # Rudd's fit
    assert D.muon_reactivation(2.5)["reactivation"] == pytest.approx(0.35, abs=0.05)
    squeezed = D.muon_reactivation(100.0)
    assert squeezed["reactivation"] > 0.9 and squeezed["fusions_limit"] > 1000
    assert D.muon_reactivation(100.0, fast_only=True)["reactivation"] > 0.85
    assert 1e-10 < D.squeeze_lifetime(100.0, 50e-6) < 2e-10


def test_xray_trapping_fibre_pinch_and_channelling():
    from nanofusion.improved import degenerate as D
    from nanofusion.improved import walls as W

    # Thermal pairs vanish in cool fuel and outnumber the fuel's own electrons in hot, thin fuel.
    assert W.pair_density(1e30, 10.0) < 1e-10 * 1e30
    assert W.pair_density(1e30, 300.0) > 1e30
    # Free-free opacity falls as T^-3.5; alphas heat ions in hot electrons and electrons in cold ones (Stix).
    ions, n_e = D.fuel_ions(D.pb_fuel(0.15), 1000.0)
    assert W.free_free_opacity(ions, n_e, 200.0) == pytest.approx(W.free_free_opacity(ions, n_e, 100.0) * 2**-3.5)
    assert W.alpha_share_classical(ions, n_e, 100.0) > 0.6 > W.alpha_share_classical(ions, n_e, 5.0)
    # The smallest proton-boron ball that repays its heat at 1e4 g/cm3 is transparent, holds ~MJ, rhoR ~13 g/cm2.
    ball = W.smallest_burning_ball(1e4)
    assert not ball["opaque"] and 1e5 < ball["energy_J"] < 1e7 and 8 < ball["rhoR_g_cm2"] < 25
    # At 1,000 g/cm3 an opaque ball is tens of kg; its X-rays and pairs hold most of its energy; it returns < 1.
    opaque = W.smallest_opaque_ball(1000.0)
    assert opaque["mass_kg"] > 1 and opaque["radiation_share"] > 0.9 and opaque["gain"] < 1
    # Z-pinch: Bennett's relation; alphas held above ~1.35 MA; at 26 MA a ~6-fold squeeze repays the heat in one
    # Alfven time, and ten Alfven times cut the squeeze tenfold.
    z = W.fibre_pinch(26e6)
    assert 2 * z["ions_per_m"] * 10.0 * P.J_PER_KEV == pytest.approx(4e-7 * math.pi * 26e6**2 / (8 * math.pi))
    assert W.alpha_confining_current() == pytest.approx(1.35e6, rel=0.01)
    assert 4 < z["squeeze"] < 8 and z["field_J_m"] > 5 * z["heat_J_m"]
    assert W.fibre_pinch(26e6, alfven_times=10.0)["squeeze"] == pytest.approx(z["squeeze"] / 10)
    assert W.fibre_pinch(1e6, fuel="6LiD")["squeeze"] > 1e4
    # Channelling: Lindhard's angle ~1.5 degrees for 300 keV deuterons along Ti rows, under 1% of an isotropic
    # source inside it, and at most ~6 times the random return.
    psi = W.critical_angle(300.0, 1, 22, 3.15e-10)
    assert 1.2 < math.degrees(psi) < 1.8 and W.axial_acceptance(psi, 26) < 0.01
    assert 0.05 < W.planar_acceptance(psi / 2.5, 13) < 0.2
    assert W.channelled_return(0.004, [(1.0, 3.0, 0.5, 1.0)]) == pytest.approx(0.024)
    assert W.channelled_return(0.004, [(0.01, 2.0, 0.9, 0.5), (0.1, 2.0, 0.9, 0.5)]) < 0.0045
    # Magnetised electrons need their cyclotron frequency above the plasma frequency: ~7e4 T in solid D-T.
    assert W.magnetised_electron_field(4.5e28, 1.0) == pytest.approx(6.8e4, rel=0.02)


def test_one_simulator_for_every_route(capsys):
    from nanofusion.__main__ import main
    from nanofusion.improved import simulator as S

    # Every route runs with its defaults and gives a finite energy account.
    for name in S.ROUTES:
        o = S.evaluate(name)
        assert math.isfinite(o.electric_gain) and o.electricity_in_J > 0 and o.stretch >= 1
    # Muon catalysis at the record returns ~0.1 of its electricity (Kelly et al. estimate ~0.14), with today's means.
    muon = S.evaluate("muon_catalysis")
    assert 0.08 < muon.electric_gain < 0.16 and muon.stretch == 1
    # Burn-up: fusion never exceeds what the fuel holds.
    assert S.deplete(1e9, 1.0) < 1.0 and S.deplete(1e-6, 1.0) == pytest.approx(1e-6, rel=1e-5)
    # Blankets: lithium adds 4.78 MeV per D-T fusion; a k_eff = 0.95 uranium blanket makes muons "net energy".
    lithium = S.evaluate("cold_catcher")
    bare = S.evaluate("cold_catcher", settings=S.Settings(blanket="none"))
    assert lithium.heat_J / bare.heat_J == pytest.approx(22.37 / 17.59, rel=1e-3)
    assert S.evaluate("muon_catalysis", settings=S.Settings(blanket="uranium", k_eff=0.95)).electric_gain > 10
    # Spin-polarised fuel raises fusion but is flagged as untested theory, and the default search leaves it out.
    spin = S.evaluate("fibre_pinch", {"spin": True})
    assert spin.fusion_J > S.evaluate("fibre_pinch").fusion_J and spin.untested and not S.evaluate("fibre_pinch").untested
    # A pinch too cool for its alphas to make up its X-rays must have them resupplied, which sinks it.
    cool = S.evaluate("fibre_pinch", {"current_ma": 100.0, "t_kev": 3.63, "alfven_times": 99.4, "squeeze": 140.0})
    assert cool.electric_gain < 0.2
    # A quick search: the pinch passes 1 only when it may stay stable far longer than any dense hot pinch has.
    found = S.search(["fibre_pinch"], samples=80, refine=40, seed=1)["fibre_pinch"]
    assert found[10.0].electric_gain < 1 < found[100.0].electric_gain
    assert found[100.0].limiting
    # Command line.
    assert main(["sim", "run", "tamped_ball", "radius_mm=0.3", "b_tesla=3000"]) == 0
    assert "electric gain" in capsys.readouterr().out
    assert main(["sim", "run", "tamped_ball", "nonsense=1"]) == 2
    assert main(["sim", "list"]) == 0


def test_simulator_new_routes_and_combinations():
    from nanofusion.improved import simulator as S

    # Ordinary inertial fusion, the reference: NIF-like inputs give ~0.15 mg of fuel, a ~600 Gbar isobaric hot spot,
    # and a target gain of ~2 (NIF's record: 4.1).
    nif = S.evaluate("compressed_icf", {"rho_g_cc": 950.0, "rho_r": 0.6, "adiabat": 3.0})
    assert 1.0 < nif.target_gain < 6.0
    reqs = {r.name: r.needed for r in nif.all_requirements}
    assert 400 < reqs["hot-spot pressure (Gbar)"] < 800
    assert 15e3 < nif.fuel_J < 50e3  # tens of kJ in the fuel, ~1-2% of NIF's laser energy
    # A crystal spark is worse than an ordinary hot spot. Pulsed-power X-rays put less of the wall-plug energy
    # into the fuel than a 10%-efficient laser: 0.15 x 0.35% against 0.10 x 1%.
    spark = S.evaluate("compressed_icf", {"ignitor": "crystal spark"})
    assert spark.electric_gain < nif.electric_gain
    xrays = S.evaluate("compressed_icf", {"rho_g_cc": 950.0, "rho_r": 0.6, "adiabat": 3.0,
                                          "driver": "pulsed-power X-rays"})
    assert xrays.electric_gain == pytest.approx(nif.electric_gain * 0.15 * 0.0035 / (0.10 * 0.01))
    # The other new routes run, with the expected scale.
    factory = S.evaluate("muon_factory")
    assert 20 < factory.drivers["laser"] / (1e3 * S.MEV) < 50  # ~34 GeV of light per muon (section 10)
    assert S.evaluate("hollow_shell").untested
    assert S.evaluate("lattice_confinement").electric_gain < 1e-12
    plain = S.evaluate("pB_catcher", {"avalanche": False})
    assert S.evaluate("pB_catcher").fusion_J == pytest.approx(plain.fusion_J / (1 - S.PB_AVALANCHE_K))
    assert S.evaluate("magnetised_liner").electric_gain < 0.01
    # Searches can hold an idea fixed, and start from known designs.
    found = S.search(["compressed_icf"], samples=40, refine=20, seed=1, fixed={"ignitor": "crystal spark"},
                     start=[{"rho_g_cc": 950.0, "rho_r": 0.6, "adiabat": 3.0}])["compressed_icf"]
    assert all(o is None or o.params["ignitor"] == "crystal spark" for o in found.values())


def test_simulator_untested_theory():
    from nanofusion.improved import simulator as S

    # Untested ideas are flagged, add no technology stretch, and the default search leaves them out.
    channelled = S.evaluate("opaque_pB", {"alpha_channelling": 0.9, "direct_conversion": True})
    assert len(channelled.untested) == 2 and not S.evaluate("opaque_pB").untested
    assert channelled.electric_gain > S.evaluate("opaque_pB").electric_gain
    plain = S.search(["opaque_pB"], samples=20, refine=10, seed=1)["opaque_pB"]
    assert all(o is None or not o.untested for o in plain.values())
    # Direct conversion: all of a proton-boron design's fusion is turned into electricity at 60% instead of 40%.
    d = S.evaluate("son_fisch_pB", {"direct_conversion": True})
    assert d.electricity_out_J == pytest.approx(d.fusion_J * S.DIRECT_EFFICIENCY)
    # Freeing stuck muons with fields doubles the fusions per muon, but the cooling it needs eats the gain.
    base = S.evaluate("muon_catalysis")
    tek = S.evaluate("muon_catalysis", {"rf_stripping": "Kulsrud/TEK"})
    assert tek.dt_fusions > 1.8 * base.dt_fusions and tek.electric_gain < 1.3 * base.electric_gain
    # Antimatter sparks: a NIF-like target needs ~1e15 antiprotons, over 100 times the most ever stored.
    spark = S.evaluate("antimatter_spark")
    need = {r.name: r.needed for r in spark.all_requirements}["antiprotons per shot"]
    assert 1e14 < need < 1e16 and spark.electric_gain < 1e-4 and spark.untested
    # Light-assisted tunnelling: the light costs far more than the extra fusion.
    assisted = S.evaluate("laser_assisted")
    assert assisted.electric_gain < 1e-3 and assisted.untested
    hot = S.evaluate("laser_assisted", {"rho_g_cc": 1000.0, "t_kev": 5.0, "radius_mm": 0.078})
    cold_model = S.evaluate("laser_assisted", {"rho_g_cc": 1000.0, "t_kev": 5.0, "radius_mm": 0.078,
                                               "model": "Thomson et al. 2024 (pessimistic)"})
    assert hot.fusion_J > cold_model.fusion_J


def test_icf_experiment_design_and_calibration():
    import importlib.util
    from pathlib import Path

    from nanofusion.improved import simulator as S

    # The sizing helpers invert the target model: density for a pressure, shell for a mass.
    rho = S.icf_shell_density(608.0, 2.45)
    rho_r = S.icf_shell_for_mass(608.0, 2.45, 0.22e-6)
    t = S.icf_target(rho, rho_r, 2.45)
    assert t["pressure_Pa"] / 1e14 == pytest.approx(608.0, rel=1e-6)
    assert t["mass_kg"] == pytest.approx(0.22e-6, rel=1e-6)
    # The best design keeps within its limits, and more fuel energy buys more fusion.
    d = S.icf_design(40e3, 700.0, 2.25)
    assert d["fuel_J"] <= 40e3 * (1 + 1e-9) and d["pressure_gbar"] <= 700.0 and d["adiabat"] >= 2.25
    assert S.icf_design(60e3, 700.0, 2.25)["fusion_J"] > d["fusion_J"]

    spec = importlib.util.spec_from_file_location("icf_experiment",
                                                  Path(__file__).parents[1] / "tools" / "icf_experiment.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cal = mod.calibrate()
    # NIF's first ignition (N221204): ~1.2 g/cm2 before burn (NIF's simulations: 1.35-1.41 without alpha heating),
    # a ceiling of ~11 MJ, so the shot reached ~30% of it; ~2% of the laser energy ends up in the fuel.
    n22 = cal["N221204"]
    assert 1.0 < n22["total_rho_r"] < 1.5 and 9e6 < n22["ceiling_J"] < 13e6
    assert 0.2 < n22["quality"] < 0.4 and 0.01 < n22["coupling"] < 0.03
    assert cal["N210808"]["quality"] < n22["quality"]
    rows = {r["step"][0]: r for r in mod.campaign(n22["coupling"], {"record": 0.79, "typical": 0.48})}
    # Repeating the record reproduces it; with today's fuel a lower adiabat does nothing, more pressure helps.
    assert rows["A"]["yield_record_J"] == pytest.approx(8.6e6, rel=0.01)
    assert rows["B"]["ceiling_J"] == pytest.approx(rows["A"]["ceiling_J"], rel=0.02)
    assert rows["C"]["ceiling_J"] > 1.1 * rows["A"]["ceiling_J"]
    # Re-sized fuel: the adiabat matters, and the upgrade's yields bracket published projections (20-30+ MJ).
    assert rows["E"]["ceiling_J"] > 1.05 * rows["D"]["ceiling_J"]
    assert 15e6 < rows["G"]["yield_record_J"] < rows["I"]["yield_record_J"] < 40e6
    assert all(r["electric_plant_laser"] < 1 for r in rows.values())
    assert mod.ice_um(0.22) == pytest.approx(65, abs=3)  # NIF's 0.22 mg is a 64.5-65 um ice layer


def test_burn1d_implosion_and_burn():
    import numpy as np

    from nanofusion.improved import burn1d as B
    from nanofusion.improved import energy as EN

    # Bosch-Hale D-T reactivity agrees with the project's own cross-section integral
    for t in (5.0, 10.0, 20.0):
        assert B.bosch_hale_dt(np.array([t]))[0] == pytest.approx(EN.reactivity("D", "T", t), rel=0.02)
    # the equation of state inverts, for D-T and for fully ionised carbon, hot and degenerate
    n_i = np.array([1e30, 1e31, 1e31])
    n_e = n_i * np.array([1.0, 1.0, 6.0])
    t = np.array([1.0, 0.05, 2.0]) * B.KEV
    assert B.temperature(B.pressure_of(t, n_i, n_e), n_i, n_e) == pytest.approx(t, rel=1e-6)
    # radiation transport conserves energy: emitted = absorbed + escaped
    rng = np.random.default_rng(0)
    t_kev = rng.uniform(0.05, 10.0, 30)
    n_e = rng.uniform(1e29, 1e32, 30)
    dr = rng.uniform(1e-7, 1e-5, 30)
    vol = rng.uniform(1e-15, 1e-12, 30)
    emitted, absorbed, escaped = B.radiation(t_kev, n_e, n_e, dr, vol)
    assert emitted.sum() == pytest.approx(absorbed.sum() + escaped, rel=1e-9)
    # N221204's clean implosion ignites; without alpha heating it doesn't (LLNL's 1-D: ~15 MJ for this design)
    clean = B.run(B.Capsule())
    cold = B.run(B.Capsule(alpha_heating=False))
    assert 10e6 < clean["yield_J"] < 40e6 and cold["yield_J"] < 0.2e6
    assert abs(clean["energy_error"]) < 0.01 and 0.15 < clean["burn_up"] < 0.4


def test_zpinch_dynamics_and_route():
    from nanofusion.improved import simulator as S
    from nanofusion.improved import zpinch as Z

    # Pease-Braginskii current near Haines' 0.433 sqrt(lnL) MA; Bennett's relation
    assert 1.0 < Z.pease_braginskii_ma(10) < 1.8
    assert Z.bennett_kev(1e6, 1e19) == pytest.approx(Z.MU0 * 1e12 / (16 * math.pi * 1e19) / Z.KEV)
    # a cold D-T fibre on Z can't take the current; preheated and fully stabilised it burns but doesn't ignite
    cold = Z.run(Z.Pinch(rise_stabilised=True))
    assert cold["peak_current_MA"] < 6 and cold["yield_J"] < 1.0
    sheared = Z.run(Z.Pinch(start_kev=0.1, shear=0.1))
    assert 0.3e6 < sheared["yield_J"] < 5e6 and sheared["burn_up"] < 0.05 and sheared["min_radius_um"] < 20
    # the 1980s NRL fibre: thermonuclear D-D neutrons far below the ~1e10 measured (mostly beam-target)
    nrl = Z.run(Z.Pinch(fuel="DD", radius_um=50, current_ma=0.64, rise_ns=130, rise_stabilised=True))
    assert nrl["neutrons"] < 1e8
    # the simulator route: untested stability is flagged and kept out of the default search
    plain = S.evaluate("dynamic_pinch")
    flow = S.evaluate("dynamic_pinch", {"stability": "sheared flow"})
    assert not plain.untested and flow.untested and flow.fusion_J > plain.fusion_J
    assert "stability" in S.THEORY_OFF
