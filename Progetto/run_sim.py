import argparse
import numpy as np
import matplotlib.pyplot as plt
from nbody import carica_json, integr_eulero, integr_scipy, integr_rk4, energia_totale, baricentro

def main():
    arg = argparse.ArgumentParser()
    arg.add_argument("--config", required=True)
    arg.add_argument("--method", choices=["euler","rk4","scipy"], default="rk4")
    arg.add_argument("--energy", action="store_true")
    arg.add_argument("--barycenter", action="store_true")
    args = arg.parse_args()

    G, t, nomi, m, r0, v0 = carica_json(args.config)

    if args.method == "euler":
        R, V = integr_eulero(r0, v0, m, t, G)
    elif args.method == "scipy":
        R, V = integr_scipy(r0, v0, m, t, G)
    else:
        R, V = integr_rk4(r0, v0, m, t, G)

    #orbite 2D
    plt.figure(figsize=(7,7))
    for i, name in enumerate(nomi):
        plt.plot(R[:, i, 0], R[:, i, 1], linewidth=1, label=name)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(fontsize=8)
    plt.title(f"Orbite ({args.method})")
    plt.show()

    if args.energy:
        E = energia_totale(R, V, m, G)
        Erel = (E - E[0]) / abs(E[0])

        plt.figure()
        plt.plot(t, Erel)
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.title("Errore relativo energia: (E(t)-E(0))/|E(0)|")
        plt.xlabel("t [yr]")
        plt.ylabel("Errore relativo")
        plt.show()

    if args.barycenter:
        Rcm = baricentro(R, m)

        # Traiettoria del baricentro
        dx = Rcm[:, 0] - Rcm[0, 0]
        dy = Rcm[:, 1] - Rcm[0, 1]

        larg = max(np.max(np.abs(dx)), np.max(np.abs(dy)))
        larg = max(larg, 1e-15)

        fig, ax = plt.subplots(figsize=(6,6))
        ax.plot(dx, dy)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-larg*1.3, larg*1.3)
        ax.set_ylim(-larg*1.3, larg*1.3)
        ax.grid(True, linestyle=":", alpha=0.6)
        ax.set_title("Traiettoria del baricentro")
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        plt.show()

        # Spostamento nel tempo
        spost = np.sqrt(np.sum((Rcm - Rcm[0])**2, axis=1))

        plt.figure()
        plt.plot(t, spost)
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.title("Spostamento del baricentro nel tempo")
        plt.xlabel("t [yr]")
        plt.ylabel("|R_CM(t) - R_CM(0)| [AU]")
        plt.show()

if __name__ == "__main__":
    main()
