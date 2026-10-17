#!/usr/bin/env python3
"""
Stima differenziale della velocita' di conduzione da una mappa LAT.

Gradiente tangenziale del LAT per faccia -> velocita' locale (1/|grad|) e
direzione di propagazione. Prototipo dello stadio 4 della pipeline.

IMPORTANTE — perche' NON si usa pyvista/VTK qui:
    leggere la mesh da .vtk e passarla per extract_surface() puo' riordinare
    o ricostruire le celle, rompendo l'associazione faccia <-> fibra. L'effetto
    dipende dalla versione di libreria: la STESSA identica simulazione (campo
    LAT uguale bit per bit) ha dato 1035.1 mm/s su una macchina e 998.7 su
    un'altra. I file testuali openCARP (.pts/.elem/.lon) sono inequivocabili:
    connettivita' esplicita, fibre nell'ordine degli elementi, precisione piena.

Unita': .pts in mm, LAT in ms  ->  1/|grad| e' in mm/ms = m/s.
        Moltiplicare per 1000 per avere mm/s.

Uso da riga di comando:
    python3 cv_from_lat.py MESH_BASENAME ACT_IGB [--validate CV_L CV_T]

Esempio (verifica della ground truth dello stadio 1):
    python3 cv_from_lat.py data/augmenta/LA_endo_with_fiber \
                           run_LSPV/act.igb --validate 1050 525

Uso come modulo:
    from cv_from_lat import read_carp_mesh, read_igb, cv_from_lat
"""

import argparse
import numpy as np


def read_carp_mesh(basename):
    """Legge una mesh openCARP dai file testuali.

    Ritorna (P, E, tags, f, t):
        P    (n_pts, 3)   coordinate dei vertici
        E    (n_el, 3)    indici dei vertici di ogni triangolo
        tags (n_el,)      tag di regione (ultima colonna del .elem)
        f    (n_el, 3)    direzione di fibra per elemento
        t    (n_el, 3)    direzione trasversale (None se .lon non ortotropo)
    """
    P = np.loadtxt(basename + '.pts', skiprows=1)

    # .elem: "Tr v0 v1 v2 tag" -> saltiamo la colonna di testo
    raw = np.loadtxt(basename + '.elem', skiprows=1,
                     usecols=(1, 2, 3, 4), dtype=np.int64)
    E, tags = raw[:, :3], raw[:, 3]

    # .lon: prima riga = numero di vettori per elemento (1 o 2)
    with open(basename + '.lon') as fh:
        n_vec = int(fh.readline().split()[0])
    lon = np.loadtxt(basename + '.lon', skiprows=1)
    f = lon[:, 0:3]
    t = lon[:, 3:6] if n_vec >= 2 else None

    return P, E, tags, f, t


def read_igb(path):
    """Legge un .igb di openCARP: header di 1024 byte, poi float32."""
    with open(path, 'rb') as fh:
        fh.seek(1024)
        return np.fromfile(fh, dtype=np.float32).astype(np.float64)


def face_gradient(P, E, u):
    """Gradiente tangenziale, per faccia, di un campo scalare per vertice.

    Formula per un campo lineare sul triangolo (p0,p1,p2) con normale n:
        grad u = [ (u1-u0) (n x (p0-p2)) + (u2-u0) (n x (p1-p0)) ] / 2A

    Ritorna (g, n, area) con g gradiente per faccia, n normali unitarie.
    """
    p0, p1, p2 = P[E[:, 0]], P[E[:, 1]], P[E[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    A2 = np.linalg.norm(n, axis=1)          # 2 * area
    n = n / A2[:, None]

    u0, u1, u2 = u[E[:, 0]], u[E[:, 1]], u[E[:, 2]]
    g = ((u1 - u0)[:, None] * np.cross(n, p0 - p2)
         + (u2 - u0)[:, None] * np.cross(n, p1 - p0)) / A2[:, None]

    return g, n, A2 / 2.0


def cv_from_lat(P, E, lat, eps=1e-9):
    """Velocita' e direzione di propagazione per faccia, dal campo LAT.

    Ritorna (speed, direction, valid, normals, area). Le facce con gradiente
    quasi nullo (fronte stazionario, degeneri) hanno valid=False e NaN.
    """
    g, n, area = face_gradient(P, E, lat)
    gn = np.linalg.norm(g, axis=1)
    valid = gn > eps

    speed = np.full(gn.shape, np.nan)
    speed[valid] = 1.0 / gn[valid]                      # mm/ms = m/s

    direction = np.full(g.shape, np.nan)
    direction[valid] = g[valid] / gn[valid][:, None]

    return speed, direction, valid, n, area


def validate_anisotropy(speed, direction, valid, fib, cv_l=None, cv_t=None):
    """Confronta la velocita' stimata con la direzione di fibra.

    Stampa la mediana sulle facce dove il fronte viaggia quasi parallelo e
    quasi perpendicolare alle fibre: devono avvicinarsi a CV_L e CV_T imposte.
    """
    v = speed[valid]
    d = direction[valid]
    f = fib[valid] / np.linalg.norm(fib[valid], axis=1)[:, None]
    cos = np.abs(np.sum(d * f, axis=1))

    out = {}
    for lo, hi, label in [(0.95, 1.01, 'parallelo'),
                          (0.00, 0.15, 'perpendicolare')]:
        sel = (cos >= lo) & (cos < hi)
        med = np.median(v[sel]) * 1000.0 if sel.sum() else np.nan
        out[label] = (int(sel.sum()), med)
        print(f'{label:16s} n={sel.sum():6d}  v mediana = {med:7.1f} mm/s')

    if cv_l and cv_t:
        el = abs(out['parallelo'][1] - cv_l) / cv_l * 100
        et = abs(out['perpendicolare'][1] - cv_t) / cv_t * 100
        ratio = out['parallelo'][1] / out['perpendicolare'][1]
        print(f'\natteso  CV_L={cv_l:.0f}  CV_T={cv_t:.0f}  rapporto={cv_l/cv_t:.2f}')
        print(f'errore  {el:.1f}%  {et:.1f}%   rapporto misurato={ratio:.2f}')

    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mesh', help='basename della mesh openCARP (senza estensione)')
    ap.add_argument('act', help='file act.igb con i tempi di attivazione')
    ap.add_argument('--validate', nargs=2, type=float, metavar=('CV_L', 'CV_T'),
                    help='confronta con le CV imposte, in mm/s')
    ap.add_argument('--save', metavar='NPZ',
                    help='salva speed/direction/valid in un .npz')
    args = ap.parse_args()

    P, E, tags, f, t = read_carp_mesh(args.mesh)
    lat = read_igb(args.act)
    print(f'pts {P.shape}  elem {E.shape}  lon {f.shape}  lat {lat.shape}')
    if lat.size != P.shape[0]:
        raise SystemExit(f'LAT ({lat.size}) e vertici ({P.shape[0]}) non corrispondono')
    if (lat < 0).any():
        print(f'ATTENZIONE: {int((lat < 0).sum())} nodi non attivati (-1): '
              f'aumentare dream.fim.max_iter')

    speed, direction, valid, n, area = cv_from_lat(P, E, lat)
    print(f'facce valide: {int(valid.sum())}/{valid.size}   '
          f'v mediana globale = {np.median(speed[valid]) * 1000:.1f} mm/s\n')

    cv_l, cv_t = args.validate if args.validate else (None, None)
    validate_anisotropy(speed, direction, valid, f, cv_l, cv_t)

    if args.save:
        np.savez(args.save, speed=speed, direction=direction,
                 valid=valid, normals=n, area=area, tags=tags)
        print(f'\nsalvato: {args.save}')


if __name__ == '__main__':
    main()
