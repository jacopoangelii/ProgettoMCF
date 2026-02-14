# Progetto Metodi Computazionali per la Fisica

Documentazione del progetto finale per il corso di Metodi Computazionali per la Fisica.

Questo progetto implementa un simulatore numerico di un sistema gravitazionale a N corpi basato sulla legge di gravitazione universale di Newton.

L’obiettivo è:

* Costruire un modulo Python in grado di integrare nel tempo le equazioni del moto di un sistema di corpi con masse e condizioni iniziali assegnate.
* Applicare il modulo allo studio del Sistema Solare.
* Simulare sistemi stellari alternativi e analizzare le velocità radiali della stella centrale, alla base del metodo Doppler per la scoperta di esopianeti.

Struttura del progetto
=========================================
```
nbody.py                            → modulo principale (dinamica N-corpi e integratori)
run_sim.py                          → simulazione e confronto tra metodi di integrazione
run_anim.py                         → animazione delle orbite
run_rv.py                           → analisi della velocità radiale
configs/                            → file JSON con configurazioni dei sistemi
   sistemasolare.json               → Sistema Solare
   sistemasolare_terra_veloce.json  → Sistema Solare con velocità della Terra modificata
   rv_giovecaldo.json               → Sistema Stella - Giove Caldo (velocità radiale)
   rv_eccentrica.json               → Sistema Stella - Pianeta con orbita eccentrica (velocità radiale)
   rv_duepianeti.json               → Sistema Stella - due pianeti (velocità radiale)
```
Il modulo `nbody.py` contiene:
* definizione delle equazioni del moto gravitazionali
*	integratori numerici:
    * Metodo di Eulero
    *	Runge–Kutta del quarto ordine (RK4)
    * solve_ivp di SciPy (opzionale)
* calcolo energia totale
* calcolo baricentro
* calcolo velocità radiale

Gli script `run_*.py` sono esempi di utilizzo del modulo per produrre i risultati richiesti.

Requisiti
=========================================
Python ≥ 3.10

Librerie necessarie:
```
pip install numpy scipy matplotlib
```

Unità utilizzate
=========================================
Le simulazioni sono effettuate in unità naturali:
* distanza: AU
* tempo: anni
* massa: masse solari

In tali unità:
$G = 4\pi^2$

Esecuzione delle Simulazioni
=========================================

Posizionarsi nella cartella principale del progetto:
```
cd .../Progetto
```

---
**1) Simulazione del Sistema Solare**

* Metodo Runge–Kutta (metodo principale):
```
python3.13 run_sim.py --config configs/sistemasolare.json --method rk4
```

* Metodo di Eulero:
```
python3.13 run_sim.py --config configs/sistemasolare.json --method euler
```

* Metodo Integrazione Scipy (opzionale):
```
python3.13 run_sim.py --config configs/sistemasolare.json --method scipy
```

* Con analisi energetica:
```
python3.13 run_sim.py --config configs/sistemasolare.json --method rk4 --energy
```

* Con analisi del baricentro:
```
python3.13 run_sim.py --config configs/sistemasolare.json --method rk4 --energy --barycenter
```
---
**2) Animazione delle Orbite**
```
python3.13 run_anim.py --config configs/sistemasolare.json
```

Con parametri opzionali:
```
python3.13 run_anim.py --config configs/sistemasolare.json --frames 500 --trail 500
```
---
**3) Studio della Sensibilità alle Condizioni Iniziali**

Simulazione con velocità terrestre modificata:
```
python3.13 run_sim.py --config configs/sistemasolare_terra_veloce.json --method rk4
```
---
**4) Velocità Radiale (Metodo Doppler)**

* Hot Jupiter:
  ```
  python3.13 run_rv.py --config configs/rv_giovecaldo.json --star_index 0 --axis x
  ```
  
* Orbita eccentrica:
  ```
  python3.13 run_rv.py --config configs/rv_eccentrica.json --star_index 0 --axis x
  ```
  
* Sistema a due pianeti:
  ```
  python3.13 run_rv.py --config configs/rv_duepianeti.json --star_index 0 --axis x
  ```
  
Metodologie Numeriche
=========================================
* Metodo di Eulero – utilizzato come riferimento per analisi di stabilità.
*	Runge–Kutta 4° ordine (RK4) – metodo principale adottato per l’integrazione.
* solve_ivp (SciPy) – metodo opzionale di confronto.

Output prodotti
=========================================
* Traiettorie orbitali nel piano XY
* Animazioni temporali delle orbite
* Errore relativo dell’energia totale
* Traiettoria del baricentro
* Spostamento nel tempo del baricentro
* Curve di velocità radiale
