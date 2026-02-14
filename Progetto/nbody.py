import json
import numpy as np

# Unità di misura: AU (distanza), yr (tempo), Msun (massa)
# In queste unità: G = 4*pi^2
G_unita = 4.0 * np.pi**2


def carica_json(path: str):
    """
    Deve leggere un file JSON così:
    {
      "units": {"distance":"AU","time":"yr","mass":"Msun"},
      "duration_yr": float,
      "steps": int,
      "G": float,
      "bodies": [{"name":str,"m":float,"r":[x,y,z],"v":[vx,vy,vz]}]
    }
    """
    with open(path, "r", encoding="utf-8") as f:
        file_caricato = json.load(f)

    units = file_caricato.get("units", {"distance": "AU", "time": "yr", "mass": "Msun"})
    if units.get("distance") != "AU" or units.get("time") != "yr" or units.get("mass") != "Msun":
        raise ValueError("Questo progetto usa unità AU/yr/Msun")

    durata_anni = float(file_caricato["duration_yr"])
    passi = int(file_caricato["steps"])
    G = float(file_caricato.get("G", G_unita))

    nomi = [b["name"] for b in file_caricato["bodies"]]
    m = np.array([b["m"] for b in file_caricato["bodies"]], dtype=float)
    r0 = np.array([b["r"] for b in file_caricato["bodies"]], dtype=float)
    v0 = np.array([b["v"] for b in file_caricato["bodies"]], dtype=float)

    if r0.shape[1] != 3 or v0.shape[1] != 3:
        raise ValueError("r e v devono essere espresse in tre dimensioni: [x,y,z]")

    t = np.linspace(0.0, durata_anni, passi)
    return G, t, nomi, m, r0, v0

def accelerazioni(r, m, G, eps=0.0):
    n = len(r)
    a = np.zeros_like(r)
    
    for i in range(n):
        for j in range(n):
            if i != j:
                diff = r[j] - r[i]
                distanzaquad = float(np.dot(diff, diff)) + eps**2
                distanza = np.sqrt(distanzaquad)
                a[i] += G * m[j] * diff / (distanza**3)
    return a

def unisci(r, v):
    #metto posizioni e velocità in un unico vettore
    r_flat = r.reshape(-1)
    v_flat = v.reshape(-1)
    return np.concatenate((r_flat, v_flat))


def separa(y, n):
    #ricavo le posizioni
    r = y[:3*n].reshape(n, 3)
    #ricavo le velocità
    v = y[3*n:].reshape(n, 3)
    return r, v


def funzione(t, y, masse, G, eps=0.0):
    #numero dei corpi
    n = len(masse)
    #separo il vettore in r e v
    r, v = separa(y, n)
    # calcolo le accelerazioni
    a = accelerazioni(r, masse, G, eps)
    #ritorno derivata: dr/dt = v, dv/dt = a
    return unisci(v, a)

def integr_eulero(r0, v0, masse, tempi, G, epsilon=0.0):
    dt = tempi[1] - tempi[0]
    n = len(masse)

    # copio le condizioni iniziali
    r = r0.copy().astype(float)
    v = v0.copy().astype(float)

    #array dove salvo tutta l'evoluzione
    posizioni = np.zeros((len(tempi), n, 3))
    velocita = np.zeros((len(tempi), n, 3))

    posizioni[0] = r
    velocita[0] = v

    #ciclo nel tempo
    for i in range(len(tempi) - 1):

        # calcolo accelerazioni all'istante attuale
        a = accelerazioni(r, masse, G, epsilon)

        #aggiorno la velocità
        v = v + a * dt

        #aggiorno la posizione
        r = r + v * dt

        #nuovi valori
        posizioni[i + 1] = r
        velocita[i + 1] = v

    return posizioni, velocita


def integr_rk4(r0, v0, masse, tempi, G, epsilon=0.0):
    dt = tempi[1] - tempi[0]
    n = len(masse)

    #metto tutto nello stesso vettore perché RK4 lavora su un y unico
    y = unisci(r0, v0).astype(float)

    posizioni = np.zeros((len(tempi), n, 3))
    velocita = np.zeros((len(tempi), n, 3))
    posizioni[0] = r0
    velocita[0] = v0

    for i in range(len(tempi) - 1):
        t_i = float(tempi[i])

        s1 = funzione(t_i, y, masse, G, epsilon)
        s2 = funzione(t_i + dt/2, y + (dt/2)*s1, masse, G, epsilon)
        s3 = funzione(t_i + dt/2, y + (dt/2)*s2, masse, G, epsilon)
        s4 = funzione(t_i + dt,   y + dt*s3,      masse, G, epsilon)

        #media pesata
        y = y + (dt/6) * (s1 + 2*s2 + 2*s3 + s4)

        r, v = separa(y, n)
        posizioni[i + 1] = r
        velocita[i + 1] = v

    return posizioni, velocita

def integr_scipy(r0, v0, masse, tempi, G, epsilon=0.0, rtol=1e-9, atol=1e-12):
    #se scipy non è installato torno agli altri metodi
    try:
        from scipy.integrate import solve_ivp
    except ImportError:
        raise RuntimeError("SciPy non è installato: usa Eulero o RK4.")

    n = len(masse)
    y0 = unisci(r0, v0).astype(float)

    #funzione che l'integratore chiama: y' = f(t, y)
    def f(t, y):
        return funzione(t, y, masse, G, epsilon)

    risolvi = solve_ivp(f, (float(tempi[0]), float(tempi[-1])), y0, t_eval=tempi, rtol=rtol, atol=atol,method="DOP853")

    if not risolvi.success:
        raise RuntimeError(f"Errore in solve_ivp: {risolvi.message}")

    #risolvi.y è (6n, len(tempi)) quindi trasposto per avere una riga per ogni istante
    Y = risolvi.y.T

    posizioni = Y[:, :3*n].reshape(len(tempi), n, 3)
    velocita = Y[:, 3*n:].reshape(len(tempi), n, 3)

    return posizioni, velocita



def energia_totale(R, V, m, G):
    """
    Energia totale E(t) = K + U
    K = 1/2 sum m_i |v_i|^2
    U = - sum_{i<j} G m_i m_j / |r_i - r_j|
    """
    T = np.zeros(R.shape[0], dtype=float)
    U = np.zeros(R.shape[0], dtype=float)
    N = len(m)

    for k in range(R.shape[0]):
        v2 = np.sum(V[k]**2, axis=1)
        T[k] = 0.5 * np.sum(m * v2)

        u = 0.0
        for i in range(N):
            for j in range(i+1, N):
                diff = R[k, j] - R[k, i]
                dist = np.sqrt(np.dot(diff, diff))
                u += -G * m[i] * m[j] / dist
        U[k] = u

    return T + U


def baricentro(R, m):
    """
    Baricentro nel tempo: (len(t),3)
    """
    M = np.sum(m)
    return (R * m[None, :, None]).sum(axis=1) / M


def radial_velocity(V, indice_corpo, n_xyz=(1.0, 0.0, 0.0)):
    """
    v_r(t) = v_body(t) · n_xyz
    """
    n = np.array(n_xyz, dtype=float)
    n = n / np.linalg.norm(n)
    return (V[:, indice_corpo, :] * n[None, :]).sum(axis=1)
