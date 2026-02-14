import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from nbody import carica_json, integr_rk4

def main():
    arg = argparse.ArgumentParser()
    arg.add_argument("--config", required=True)
    arg.add_argument("--trail", type=int, default=400)
    arg.add_argument("--frames", type=int, default=400)
    arg.add_argument("--interval", type=int, default=30)
    args = arg.parse_args()

    G, t, nomi, m, r0, v0 = carica_json(args.config)
    R, V = integr_rk4(r0, v0, m, t, G)

    x = R[:, :, 0]
    y = R[:, :, 1]
    n_steps, N = x.shape

    #prendo solo x e y
    x = R[:, :, 0]
    y = R[:, :, 1]
    passi_totali, N = x.shape

    #quanti step salto per arrivare al numero di frame richiesti
    salto = max(1, passi_totali // args.frames)

    #limiti del grafico
    x_tutte = x.reshape(-1)
    y_tutte = y.reshape(-1)
    limite = 1.1 * max(np.percentile(np.abs(x_tutte), 99), np.percentile(np.abs(y_tutte), 99))

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-limite, limite)
    ax.set_ylim(-limite, limite)
    ax.grid(True, linestyle=":", alpha=0.6)

    #un punto + una scia per ogni corpo
    punti = [ax.plot([], [], marker="o", linestyle="")[0] for _ in range(N)]
    scie  = [ax.plot([], [], linewidth=1, alpha=0.6)[0] for _ in range(N)]

    ax.legend(punti, nomi, loc="lower right", fontsize=8)
    etichetta = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top")

    def init():
        for p, s in zip(punti, scie):
            p.set_data([], [])
            s.set_data([], [])
        etichetta.set_text("")
        return punti + scie + [etichetta]

    def update(frame):
        k = min(frame * salto, passi_totali - 1)

        etichetta.set_text(f"t = {t[k]:.3f} yr  |  step {k}/{passi_totali-1}")

        k0 = max(0, k - args.trail)
        for i in range(N):
            # IMPORTANTISSIMO: set_data vuole liste/array, non scalari
            punti[i].set_data([x[k, i]], [y[k, i]])
            scie[i].set_data(x[k0:k+1, i], y[k0:k+1, i])

        return punti + scie + [etichetta]

    ani = FuncAnimation(fig, update, frames=args.frames, init_func=init,
                        blit=True, interval=args.interval)
    plt.show()
    
if __name__ == "__main__":
    main()
