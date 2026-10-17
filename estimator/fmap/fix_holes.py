#!/usr/bin/env python3
"""
fix_holes.py  -  Chiude i loop di bordo piccoli nelle mesh PLY atriali.
                 NON aggiunge vertici nuovi. Fan-triangolazione dei loop.

Uso: python3 fix_holes.py mesh1.ply mesh2.ply
     Sovrascrive i file in place (stesse property vertex invariate).
"""

import sys
import math
from collections import defaultdict


# ── PLY I/O ──────────────────────────────────────────────────────────────────

def read_ply(path):
    """Legge un PLY ASCII. Restituisce (header_lines, prop_names, verts, faces).
    verts: list of list[float]  (ogni riga = tutti i campi del vertice)
    faces: list of list[int]
    """
    with open(path, 'r') as fh:
        lines = fh.readlines()

    i = 0
    n_verts = n_faces = 0
    props = []
    elem = None
    header = []
    while lines[i].strip() != 'end_header':
        header.append(lines[i])
        tok = lines[i].strip().split()
        if tok[0] == 'element':
            elem = tok[1]
            if elem == 'vertex':
                n_verts = int(tok[2])
            elif elem == 'face':
                n_faces = int(tok[2])
        elif tok[0] == 'property' and elem == 'vertex':
            props.append(tok[-1])
        i += 1
    header.append(lines[i])  # end_header
    i += 1

    verts = []
    for _ in range(n_verts):
        verts.append([float(v) for v in lines[i].strip().split()])
        i += 1

    faces = []
    for _ in range(n_faces):
        tok = lines[i].strip().split()
        n = int(tok[0])
        faces.append([int(tok[j + 1]) for j in range(n)])
        i += 1

    return header, props, verts, faces


def write_ply(path, header, props, verts, faces):
    """Scrive PLY ASCII con conteggio aggiornato di vertici e facce."""
    with open(path, 'w') as fh:
        # Riscrivi header aggiornando i conteggi element
        for line in header:
            tok = line.strip().split()
            if tok[0] == 'element' and tok[1] == 'vertex':
                fh.write(f'element vertex {len(verts)}\n')
            elif tok[0] == 'element' and tok[1] == 'face':
                fh.write(f'element face {len(faces)}\n')
            else:
                fh.write(line)
        # Vertici
        for v in verts:
            fh.write(' '.join(f'{x}' for x in v) + '\n')
        # Facce
        for f in faces:
            fh.write(f'{len(f)} ' + ' '.join(str(vi) for vi in f) + '\n')


# ── Geometria di supporto ────────────────────────────────────────────────────

def tri_area(v0, v1, v2):
    ex, ey, ez = v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]
    fx, fy, fz = v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]
    cx = ey*fz - ez*fy
    cy = ez*fx - ex*fz
    cz = ex*fy - ey*fx
    return 0.5 * math.sqrt(cx*cx + cy*cy + cz*cz)


# ── Rilevamento boundary loops ───────────────────────────────────────────────

def find_boundary_loops(faces):
    """
    Restituisce:
      loops  : list of list[int] — vertici in ordine consistente con l'orientamento
               delle facce adiacenti (gli edge del loop sono gap-directed edges).
      n_bedge: int — numero totale di edge di bordo.

    Ordine loop: v0→v1→...→vN-1→v0 sono i 'gap directed edges', ovvero
    l'opposto dell'edge diretto usato dalla faccia adiacente.
    La fan-triangolazione (v0,v1,v2), (v0,v2,v3), ... è direttamente corretta.
    """
    edge_count = defaultdict(int)
    for f in faces:
        for j in range(len(f)):
            a, b = f[j], f[(j+1) % len(f)]
            edge_count[(min(a, b), max(a, b))] += 1

    # Per ogni edge di bordo, la faccia adiacente usa a→b → gap usa b→a
    nxt = {}  # gap_start → gap_end
    for f in faces:
        for j in range(len(f)):
            a, b = f[j], f[(j+1) % len(f)]
            canonical = (min(a, b), max(a, b))
            if edge_count[canonical] == 1:
                # faccia usa a→b; gap diretto: b→a
                nxt[b] = a

    n_bedge = len(nxt)
    visited = set()
    loops = []
    for start in nxt:
        if start in visited:
            continue
        loop = [start]
        visited.add(start)
        cur = start
        while True:
            nxt_v = nxt.get(cur)
            if nxt_v is None or nxt_v in visited:
                break
            loop.append(nxt_v)
            visited.add(nxt_v)
            cur = nxt_v
        loops.append(loop)

    return loops, n_bedge


# ── Diagnosi ─────────────────────────────────────────────────────────────────

def diagnose(verts, faces, label):
    """Stampa tabella dei loop. Restituisce (loops, total_area, bbox_diag)."""
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    zs = [v[2] for v in verts]
    bbox_diag = math.sqrt((max(xs)-min(xs))**2 +
                          (max(ys)-min(ys))**2 +
                          (max(zs)-min(zs))**2)
    total_area = sum(tri_area(verts[f[0]], verts[f[1]], verts[f[2]])
                     for f in faces if len(f) == 3)

    loops, n_bedge = find_boundary_loops(faces)

    print(f'\n  [{label}]  {len(verts)} vertici, {len(faces)} facce')
    print(f'  BBox diag={bbox_diag:.2f}, Area totale={total_area:.2f}')
    print(f'  Boundary edges: {n_bedge}, Loop di bordo: {len(loops)}')
    print(f'  {"Loop":>5}  {"Verts":>5}  {"Perimetro":>10}  {"Area proxy":>12}  {"% totale":>9}  Stato')
    for k, lp in enumerate(sorted(loops, key=len, reverse=True)):
        perim = sum(
            math.sqrt(sum((verts[lp[j]][ii] - verts[lp[(j+1) % len(lp)]][ii])**2
                          for ii in range(3)))
            for j in range(len(lp))
        )
        lx = [verts[v][0] for v in lp]
        ly = [verts[v][1] for v in lp]
        lz = [verts[v][2] for v in lp]
        dx, dy, dz = max(lx)-min(lx), max(ly)-min(ly), max(lz)-min(lz)
        area_proxy = 0.5*(dx*dy + dy*dz + dz*dx)
        pct = 100*area_proxy/total_area if total_area > 0 else 0
        stato = 'OK-piccolo' if (pct < 5.0 and len(lp) <= 20) else 'GRANDE-SEGNALARE'
        print(f'  {k:>5}  {len(lp):>5}  {perim:>10.4f}  {area_proxy:>12.6f}  {pct:>8.4f}%  {stato}')

    return loops, total_area, bbox_diag


# ── Weld near-coincident ─────────────────────────────────────────────────────

def weld_boundary_vertices(verts, faces, loops, bbox_diag, eps_rel=1e-3):
    """
    Controlla SOLO i vertici di bordo per quasi-coincidenza reciproca.
    Fonde solo coppie di vertici di bordo entro eps_rel * bbox_diag.
    I vertici interni non vengono mai toccati.
    Restituisce (new_verts, new_faces, n_fused, old_to_new).
    Se n_fused == 0 restituisce gli originali invariati.
    """
    threshold = eps_rel * bbox_diag
    boundary_verts = list({v for loop in loops for v in loop})

    old_to_new = list(range(len(verts)))
    merged = set()

    for ii, i in enumerate(boundary_verts):
        if i in merged:
            continue
        for j in boundary_verts[ii+1:]:
            if j in merged:
                continue
            d = math.sqrt(sum((verts[i][k] - verts[j][k])**2 for k in range(3)))
            if d < threshold:
                old_to_new[j] = old_to_new[i]
                merged.add(j)

    n_fused = len(merged)
    if n_fused == 0:
        return verts, faces, 0, old_to_new

    # Remap: solo i boundary vertex fusi cambiano indice
    counter = 0
    new_idx = {}
    new_verts = []
    for i in range(len(verts)):
        rep = old_to_new[i]
        if rep not in new_idx:
            new_idx[rep] = counter
            counter += 1
            new_verts.append(verts[rep])
    final_map = [new_idx[old_to_new[i]] for i in range(len(verts))]
    new_faces = []
    for f in faces:
        nf = [final_map[v] for v in f]
        if len(set(nf)) == len(nf):
            new_faces.append(nf)

    return new_verts, new_faces, n_fused, final_map


# ── Chiusura loop (fan triangolazione + gestione peninsula faces) ─────────────

LOOP_LARGE_VERT_THRESHOLD = 20
LOOP_LARGE_AREA_PCT = 5.0  # %


def build_edge_count(faces):
    ec = defaultdict(int)
    for f in faces:
        for j in range(len(f)):
            a, b = f[j], f[(j+1) % len(f)]
            ec[(min(a, b), max(a, b))] += 1
    return ec


def find_valid_fan_start(loop, edge_count):
    """
    Prova tutti i vertici del loop come punto di partenza della fan.
    Restituisce l'indice nel loop di un vertice valido (fan senza non-manifold),
    oppure -1 se nessuno funziona.
    """
    n = len(loop)
    for si in range(n):
        ok = True
        v0 = loop[si]
        for i in range(1, n - 1):
            vi  = loop[(si + i) % n]
            vi1 = loop[(si + i + 1) % n]
            for a, b in [(v0, vi), (vi, vi1), (vi1, v0)]:
                if edge_count.get((min(a, b), max(a, b)), 0) >= 2:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return si
    return -1


def fan_has_nonmanifold(loop, edge_count):
    """Vero se NESSUN punto di partenza nella fan evita edge non-manifold."""
    return find_valid_fan_start(loop, edge_count) < 0


def close_single_loop(loop, faces, edge_count, label, k):
    """
    Chiude un singolo loop.
    Se la fan standard crea non-manifold, cerca e rimuove la "peninsula face"
    (faccia che ha 2 edge consecutivi sul boundary), poi ritenta.
    Restituisce (new_faces, n_added, n_removed).
    """
    # Size check
    if len(loop) > LOOP_LARGE_VERT_THRESHOLD:
        print(f'\n  ERRORE [{label}] Loop {k}: {len(loop)} vertici > soglia.')
        sys.exit(2)

    si = find_valid_fan_start(loop, edge_count)
    if si >= 0:
        # Caso normale: fan da vertice valido
        n = len(loop)
        v0 = loop[si]
        ordered = [loop[(si + j) % n] for j in range(n)]
        new_tris = [[v0, ordered[i], ordered[i+1]] for i in range(1, n - 1)]
        new_faces = faces + new_tris
        for tri in new_tris:
            for j in range(3):
                a, b = tri[j], tri[(j+1) % 3]
                edge_count[(min(a, b), max(a, b))] += 1
        return new_faces, len(new_tris), 0

    # Caso degenere: cerca "peninsula face" (faccia con 2 edge consecutivi sul loop)
    boundary_edges = {(min(loop[j], loop[(j+1) % len(loop)]),
                       max(loop[j], loop[(j+1) % len(loop)]))
                      for j in range(len(loop))}

    peninsula_fi = None
    for fi, f in enumerate(faces):
        f_edges = [(min(f[j], f[(j+1) % 3]), max(f[j], f[(j+1) % 3]))
                   for j in range(3)]
        on_boundary = [e for e in f_edges if e in boundary_edges]
        if len(on_boundary) == 2:
            peninsula_fi = fi
            break

    if peninsula_fi is None:
        print(f'\n  ERRORE [{label}] Loop {k}: fan non-manifold ma nessuna peninsula face trovata.')
        print('  Caso topologicamente non gestibile senza nuovi vertici.')
        sys.exit(4)

    pf = faces[peninsula_fi]
    print(f'  Loop {k}: rimossa peninsula face fi={peninsula_fi}: {pf}')
    for j in range(3):
        a, b = pf[j], pf[(j+1) % 3]
        edge_count[(min(a, b), max(a, b))] -= 1

    new_faces_list = [f for i, f in enumerate(faces) if i != peninsula_fi]

    # Ri-rileva il loop (ora più piccolo)
    new_loops, _ = find_boundary_loops(new_faces_list)
    loop_vset = set(loop)
    new_loop = None
    for nl in new_loops:
        if any(v in loop_vset for v in nl):
            new_loop = nl
            break

    if new_loop is None:
        print(f'  ERRORE [{label}] Loop {k}: nessun loop residuo dopo rimozione peninsula.')
        sys.exit(4)

    if fan_has_nonmanifold(new_loop, edge_count):
        print(f'  ERRORE [{label}] Loop {k}: ancora non-manifold dopo rimozione peninsula.')
        sys.exit(4)

    v0 = new_loop[0]
    new_tris = [[v0, new_loop[i], new_loop[i+1]] for i in range(1, len(new_loop) - 1)]
    final_faces = new_faces_list + new_tris
    for tri in new_tris:
        for j in range(3):
            a, b = tri[j], tri[(j+1) % 3]
            edge_count[(min(a, b), max(a, b))] += 1
    return final_faces, len(new_tris), 1


def close_loops(verts, faces, loops, total_area, label):
    """
    Chiude tutti i loop.  NON aggiunge vertici nuovi.
    Gestisce il caso degenere di "peninsula face" che crea non-manifold.
    """
    edge_count = build_edge_count(faces)
    current_faces = list(faces)
    added = 0
    removed_peninsula = 0
    report = []

    for k, loop in enumerate(loops):
        # Soglia area
        lx = [verts[v][0] for v in loop]
        ly = [verts[v][1] for v in loop]
        lz = [verts[v][2] for v in loop]
        dx, dy, dz = max(lx)-min(lx), max(ly)-min(ly), max(lz)-min(lz)
        area_proxy = 0.5*(dx*dy + dy*dz + dz*dx)
        pct = 100*area_proxy/total_area if total_area > 0 else 0
        if pct > LOOP_LARGE_AREA_PCT:
            print(f'\n  ERRORE [{label}] Loop {k}: {pct:.4f}% > {LOOP_LARGE_AREA_PCT}% soglia area.')
            print('  Loop troppo grande: NON è solo connettività mancante.')
            sys.exit(2)

        current_faces, n_added, n_rem = close_single_loop(
            loop, current_faces, edge_count, label, k)
        added += n_added
        removed_peninsula += n_rem
        report.append((k, len(loop), n_added, n_rem))

    return current_faces, added, removed_peninsula, report


# ── Verifica integrità ───────────────────────────────────────────────────────

def verify(verts, faces, label):
    """Controlla: 0 loop, componenti connesse, non-manifold, facce degeneri."""
    loops, n_bedge = find_boundary_loops(faces)

    # Componenti connesse (union-find)
    parent = list(range(len(verts)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        parent[find(a)] = find(b)
    for f in faces:
        for j in range(1, len(f)):
            union(f[0], f[j])
    roots = len(set(find(i) for i in range(len(verts))))

    # Non-manifold: edge con >2 facce
    edge_count = defaultdict(int)
    for f in faces:
        for j in range(len(f)):
            a, b = f[j], f[(j+1) % len(f)]
            edge_count[(min(a,b), max(a,b))] += 1
    nm_edges = sum(1 for c in edge_count.values() if c > 2)

    # Facce degeneri
    deg = sum(1 for f in faces if len(set(f)) < len(f))

    ok = (len(loops) == 0 and roots == 1 and nm_edges == 0 and deg == 0)
    status = 'OK' if ok else 'PROBLEMI'
    print(f'\n  [{label}] Verifica integrità: {status}')
    print(f'    Loop di bordo residui : {len(loops)}  (atteso 0)')
    print(f'    Componenti connesse   : {roots}     (atteso 1)')
    print(f'    Edge non-manifold     : {nm_edges}     (atteso 0)')
    print(f'    Facce degeneri        : {deg}     (atteso 0)')
    return ok


# ── Main ─────────────────────────────────────────────────────────────────────

def repair(path, label):
    print(f'\n{"="*60}')
    print(f'Riparazione: {label}  ({path})')
    print('='*60)

    header, props, verts, faces = read_ply(path)
    n_verts_orig = len(verts)
    n_faces_orig = len(faces)

    # 1. Diagnosi
    print('\n--- 1. Diagnosi loop PRIMA della riparazione ---')
    loops, total_area, bbox_diag = diagnose(verts, faces, label)

    if not loops:
        print('  Nessun loop di bordo: mesh già chiusa.')
        return

    # 2. Weld (verifica near-coincident tra soli vertici di bordo)
    print('\n--- 2. Weld vertici di bordo quasi-coincidenti ---')
    verts_w, faces_w, n_fused, _ = weld_boundary_vertices(
        verts, faces, loops, bbox_diag)
    if n_fused == 0:
        print('  Nessun vertice di bordo quasi-coincidente (weld non necessario).')
    else:
        print(f'  Fusi {n_fused} vertici di bordo.')
        verts, faces = verts_w, faces_w
        loops, total_area, _ = diagnose(verts, faces, label + ' post-weld')

    # 3. Chiusura loop
    print('\n--- 3. Chiusura loop con fan triangolazione ---')
    new_faces, n_added, n_rem, rep = close_loops(verts, faces, loops, total_area, label)
    for k, n_vloop, n_tri, n_peninsula in rep:
        peninsula_note = f'  (rimossa {n_peninsula} peninsula face)' if n_peninsula else ''
        print(f'  Loop {k}: {n_vloop} vertici → {n_tri} triangoli aggiunti{peninsula_note}')
    print(f'  Totale triangoli aggiunti: {n_added}  (peninsula faces rimosse: {n_rem})')

    # 4. Verifica integrità
    print('\n--- 4. Verifica integrità POST-riparazione ---')
    print(f'  Vertici: {len(verts)} (prima: {n_verts_orig}, delta: {len(verts)-n_verts_orig})')
    print(f'  Facce  : {len(new_faces)} (prima: {n_faces_orig}, delta: {len(new_faces)-n_faces_orig})')
    ok = verify(verts, new_faces, label)

    if not ok:
        print('  ATTENZIONE: verifica fallita, file NON sovrascritto.')
        sys.exit(3)

    print('\n  Vertici originali e LAT (quality) invariati: SI')
    print(f'  Vertici aggiunti: {len(verts) - n_verts_orig} (atteso 0)')
    if n_rem > 0:
        print(f'  Peninsula faces rimosse: {n_rem} (facce con orientamento errato, sostituite)')

    # 5. Salvataggio
    write_ply(path, header, props, verts, new_faces)
    print(f'\n  Salvato: {path}')


def main():
    if len(sys.argv) < 3:
        print('Uso: python3 fix_holes.py source.ply target.ply')
        sys.exit(1)

    repair(sys.argv[1], 'SOURCE (ritmo normale)')
    repair(sys.argv[2], 'TARGET (FA)')

    print('\n\nRiepilogo: entrambe le mesh riparate e salvate.')


if __name__ == '__main__':
    main()
