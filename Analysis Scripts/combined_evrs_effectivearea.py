import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from crpropa import *


ENERGY_MIN = 1.0
ENERGY_MAX = 10.0
ENERGY_STEP = 0.5
TRIALS_PER_EVENT = 1000
MAX_ATTEMPTS = 5
TARGET_ENERGY = 5.0


def angdistance(l1, l2, b1, b2):
    # Clamp for numerical stability to avoid NaN from arccos domain drift.
    cos_d = np.sin(b1) * np.sin(b2) + np.cos(b1) * np.cos(b2) * np.cos(l2 - l1)
    cos_d = np.clip(cos_d, -1.0, 1.0)
    return np.rad2deg(np.arccos(cos_d))


def make_simulation(field_with_random):
    field = JF12Field(field_with_random)
    sim = ModuleList()
    sim.add(PropagationCK(field, 1e-8, 0.5 * parsec, 15 * parsec))
    obs = Observer()
    obs.add(ObserverSurface(Sphere(Vector3d(0), 20 * kpc)))
    sim.add(obs)
    return sim


def run_single_backtracking(sim, energy, lon, lat):
    pid = -nucleusId(1, 1)
    direction = Vector3d()
    direction.setRThetaPhi(1, lat, lon)
    position = Vector3d(-8.5, 0, 0) * kpc

    p = ParticleState(pid, energy, position, direction)
    c = Candidate(p)
    sim.run(c)
    d = c.current.getDirection()
    out_lon = d.getPhi()
    out_lat = np.pi / 2 - d.getTheta()
    return out_lon, out_lat


def run_trial_distances(energy, lon, lat, ref_lon, ref_lat, trials):
    distances = []
    for _ in range(trials):
        attempt = 0
        while attempt < MAX_ATTEMPTS:
            sim_rand = make_simulation(True)
            out_lon, out_lat = run_single_backtracking(sim_rand, energy, lon, lat)
            if np.isfinite(out_lon) and np.isfinite(out_lat):
                distances.append(angdistance(out_lon, ref_lon, out_lat, ref_lat))
                break
            attempt += 1
        if attempt == MAX_ATTEMPTS:
            distances.append(np.nan)
    return distances


def main():
    start = time.time()
    df = pd.read_csv("Cut50EeV.csv")
    l_vals = df["L"].values
    b_vals = df["B"].values

    energies = np.arange(ENERGY_MIN, ENERGY_MAX + 0.5 * ENERGY_STEP, ENERGY_STEP)
    n_values_all = list(range(1, 41))
    n_lines_all = {n: np.rad2deg(np.arccos(1 - 1 / (6 * n**2))) for n in n_values_all}
    # Visualization guide lines are kept as before.
    n_values_plot = [3, 6, 9, 12, 15, 18]
    mean_cut = []
    all_distances = {}

    with open("EvsR.csv", "w") as f_csv, open("EvsR.txt", "w") as f_txt:
        header = "E, mean, sigma, nearest_nside, nearest_nside_value\n"
        f_csv.write(header)
        f_txt.write(header)
        for energy in energies:
            energy_distances = []
            for lon, b in zip(l_vals, b_vals):
                lat = np.pi / 2 - b
                sim_ref = make_simulation(False)
                ref_lon, ref_lat = run_single_backtracking(sim_ref, energy, lon, lat)

                d_list = run_trial_distances(
                    energy=energy,
                    lon=lon,
                    lat=lat,
                    ref_lon=ref_lon,
                    ref_lat=ref_lat,
                    trials=TRIALS_PER_EVENT,
                )
                energy_distances.extend(d_list)

            energy_distances = np.asarray(energy_distances, dtype=float)
            mu = np.nanmean(energy_distances)
            sigma = np.nanstd(energy_distances)
            nearest_nside = min(n_values_all, key=lambda n: abs(mu - n_lines_all[n]))
            nearest_nside_value = n_lines_all[nearest_nside]
            mean_cut.append(mu)
            all_distances[energy] = energy_distances
            row = f"{energy:.6f}, {mu:.6f}, {sigma:.6f}, {nearest_nside}, {nearest_nside_value:.6f}\n"
            f_csv.write(row)
            f_txt.write(row)
            print(
                f"E={energy:.1f} EeV, mean={mu:.6f}, sigma={sigma:.6f}, "
                f"nearest_nside={nearest_nside}, nearest_nside_value={nearest_nside_value:.6f}"
            )

    target_distances = all_distances[TARGET_ENERGY]
    target_mu = np.nanmean(target_distances)
    target_sigma = np.nanstd(target_distances)

    plt.figure(figsize=(10, 6))
    plt.scatter(energies, mean_cut, c="r", s=30, label="mean distance")
    colors = ["orange", "yellow", "green", "indigo", "magenta", "purple"]
    for n, c in zip(n_values_plot, colors):
        plt.axhline(n_lines_all[n], color=c, linestyle="--", label=rf"$N_{{\mathrm{{side}}}}$={n}")
    plt.title("Average Distance vs Energy", fontsize=40)
    plt.xlabel("Energy (EeV)", fontsize=35)
    plt.ylabel("Average Distance (°)", fontsize=35)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(fontsize=25)

    hist_data_linear = target_distances[np.isfinite(target_distances)]
    hist_data_log = target_distances[np.isfinite(target_distances) & (target_distances > 0)]

    # Linear-style histogram (previous behavior)
    plt.figure(figsize=(10, 6))
    plt.hist(hist_data_linear, bins=30, edgecolor="blue", alpha=0.7, label="Distances")
    plt.axvline(target_mu, color="red", linestyle="--", linewidth=3, label=f"Mean : {target_mu:.2f}")
    plt.title(f"Effective Distance for {TARGET_ENERGY:g}EeV, 1000 Samples each location", fontsize=30)
    plt.xlabel("Angular Distance (°)", fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.yscale("log")

    # Log-scale histogram (log bins + log axes)
    if hist_data_log.size > 1 and hist_data_log.min() < hist_data_log.max():
        hist_bins_log = np.logspace(np.log10(hist_data_log.min()), np.log10(hist_data_log.max()), 31)
    else:
        hist_bins_log = 30

    plt.figure(figsize=(10, 6))
    plt.hist(hist_data_log, bins=hist_bins_log, edgecolor="blue", alpha=0.7, label="Distances")
    plt.axvline(target_mu, color="red", linestyle="--", linewidth=3, label=f"Mean : {target_mu:.2f}")
    plt.title(f"Effective Distance for {TARGET_ENERGY:g}EeV (Log-Log), 1000 samples each location", fontsize=30)
    plt.xlabel("Angular Distance (°)", fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.xscale("log")
    plt.yscale("log")

    print(f"Elapsed: {time.time() - start:.2f} s")
    plt.show()


if __name__ == "__main__":
    main()
