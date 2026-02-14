import argparse
import matplotlib.pyplot as plt
from nbody import carica_json, integr_rk4, radial_velocity

def main():
    arg = argparse.ArgumentParser()
    arg.add_argument("--config", required=True)
    arg.add_argument("--star_index", type=int, default=0)
    arg.add_argument("--axis", choices=["x","y","z"], default="x")
    args = arg.parse_args()

    mappa_assi = {"x": (1,0,0), "y": (0,1,0), "z": (0,0,1)}
    n_xyz = mappa_assi[args.axis]

    G, t, nomi, m, r0, v0 = carica_json(args.config)
    R, V = integr_rk4(r0, v0, m, t, G)

    vr = radial_velocity(V, args.star_index, n_xyz=n_xyz)

    plt.figure()
    plt.plot(t, vr)
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.title(f"Velocità radiale: {nomi[args.star_index]} lungo {args.axis}")
    plt.xlabel("t [yr]")
    plt.ylabel("v_r [AU/yr]")
    plt.show()

if __name__ == "__main__":
    main()
