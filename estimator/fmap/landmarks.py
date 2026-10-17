#!/usr/bin/env python3
"""
Picking interattivo di landmark (pyvista), in-process, pre-fmap.

I landmark corrispondono per POSIZIONE nella lista: il punto k della mesh
sorgente deve essere lo stesso punto anatomico del punto k della mesh
target. Per questo la sessione di picking non usa nomi anatomici fissi:
ogni punto è identificato da una lettera progressiva (A, B, C, ...),
assegnata nell'ordine di click, la stessa su entrambe le mesh. Il numero
di punti è libero, ma le due mesh devono finire con lo stesso numero di
punti nello stesso ordine.

Caching: gli indici vengono salvati accanto alla mesh in
`<mesh_basename>.landmarks.txt` (un indice per riga). Se il file esiste,
viene riusato senza aprire alcuna finestra; il picking interattivo avviene
solo se il file manca, se è richiesto esplicitamente un re-pick, oppure se
la cache risulta stale (vedi sotto).

La cache è una lista di *indici di vertice*, valida solo per la mesh su cui
è stata generata. Per rilevare una cache stale (rigenerata su una mesh
diversa) accanto agli indici viene salvato un riferimento alla mesh in un
file affiancato `<mesh_basename>.landmarks.meta` (numero di vertici + hash
delle coordinate). Se il riferimento non combacia con la mesh corrente, o
gli indici sono fuori range, la cache viene ignorata e si ripicca da capo
(niente assert/crash su un conteggio di vertici hardcoded).
"""

import hashlib
import json
import os
import sys

import numpy as np


# ─── Cache ──────────────────────────────────────────────────────────────────

def landmarks_cache_path(mesh_path):
    base, _ = os.path.splitext(mesh_path)
    return base + ".landmarks.txt"


def landmarks_meta_path(cache_path):
    base, _ = os.path.splitext(cache_path)
    return base + ".meta"


def load_landmark_cache(path):
    with open(path, "r") as f:
        return [int(line.strip()) for line in f if line.strip()]


def save_landmark_cache(path, indices, V=None):
    with open(path, "w") as f:
        for idx in indices:
            f.write(f"{int(idx)}\n")
    if V is not None:
        save_landmark_meta(path, V)


def compute_mesh_signature(V):
    """Riferimento compatto alla geometria di una mesh: numero di vertici +
    hash delle coordinate. Usato per capire se una landmark cache è stata
    generata su questa stessa mesh o su una diversa (stale)."""
    V = np.ascontiguousarray(np.asarray(V), dtype=np.float64)
    return {
        "n_vertices": int(len(V)),
        "vertex_hash": hashlib.sha256(V.tobytes()).hexdigest(),
    }


def save_landmark_meta(cache_path, V):
    with open(landmarks_meta_path(cache_path), "w") as f:
        json.dump(compute_mesh_signature(V), f)


def load_landmark_meta(cache_path):
    path = landmarks_meta_path(cache_path)
    if not os.path.exists(path):
        return None
    with open(path, "r") as f:
        return json.load(f)


def cache_matches_mesh(cache_path, V, indices):
    """True se `indices` è utilizzabile per la mesh corrente `V`.

    Falso (cache stale) se gli indici sono fuori range, oppure se esiste un
    riferimento mesh salvato (`.landmarks.meta`) e non combacia con `V`: in
    tal caso la cache va ignorata invece di far crashare la pipeline.
    """
    n = len(V)
    if any(idx < 0 or idx >= n for idx in indices):
        return False
    meta = load_landmark_meta(cache_path)
    if meta is None:
        # Cache legacy senza riferimento mesh: nessun modo di verificare,
        # ci si fida (comportamento identico a prima di questo controllo).
        return True
    sig = compute_mesh_signature(V)
    return (
        meta.get("n_vertices") == sig["n_vertices"]
        and meta.get("vertex_hash") == sig["vertex_hash"]
    )


# ─── Guardia di collinearità ────────────────────────────────────────────────

def check_landmark_spread(V, indices, mesh_path, ratio_threshold=0.1):
    """Avvisa se i landmark sono quasi collineari (poco spread fuori dall'asse
    dominante): un rapporto tra autovalori piccolo significa che i landmark
    non vincolano l'orientazione lungo quella direzione, e la fmap tende a
    flippare senza che l'orientation-preserving riesca a correggerla.
    """
    pts = np.asarray(V)[indices]
    n = len(pts)
    if n < 3:
        print(f"[landmarks] {mesh_path}: solo {n} landmark, spread 3D non valutabile (serve n>=3).")
        return None
    cov = np.cov((pts - pts.mean(axis=0)).T)
    eigvals = np.clip(np.sort(np.linalg.eigvalsh(cov))[::-1], 0, None)
    if eigvals[0] <= 1e-12:
        print(f"[landmarks] {mesh_path}: landmark degeneri (varianza ~0).")
        return eigvals
    r2, r3 = eigvals[1] / eigvals[0], eigvals[2] / eigvals[0]
    if r2 < ratio_threshold or r3 < ratio_threshold:
        print(
            f"[landmarks] AVVISO: landmark quasi collineari su {mesh_path} "
            f"(λ2/λ1={r2:.3f}, λ3/λ1={r3:.3f}, soglia={ratio_threshold}) — "
            f"non vincolano l'orientazione, atteso flip. Aggiungi punti fuori "
            f"dall'asse dominante (es. parete posteriore, tetto, setto)."
        )
    return eigvals


def _has_display():
    if sys.platform.startswith("linux"):
        # pyvista/VTK (pip wheels) render through X11/GLX even on a Wayland
        # session (via XWayland) — WAYLAND_DISPLAY alone is not sufficient,
        # vtkXOpenGLRenderWindow needs DISPLAY set to a reachable X server.
        return bool(os.environ.get("DISPLAY"))
    return True


def _has_trame():
    try:
        import trame  # noqa: F401
        import pyvista.trame  # noqa: F401
    except ImportError:
        return False
    return True


def _letter_for_index(i):
    """Lettera stile colonna-foglio-di-calcolo per l'indice i (0-based): A, B, ..., Z, AA, AB, ..."""
    letters = ""
    i += 1
    while i > 0:
        i, rem = divmod(i - 1, 26)
        letters = chr(65 + rem) + letters
    return letters


def _prompt_text(mesh_path, picked):
    done = ", ".join(_letter_for_index(i) for i in range(len(picked)))
    next_letter = _letter_for_index(len(picked))
    return (
        f"{os.path.basename(mesh_path)}  —  prossimo punto: {next_letter}\n"
        f"punti già messi: {done if picked else '(nessuno)'}"
    )


# ─── Picking interattivo (finestra nativa) ──────────────────────────────────

def pick_landmarks(mesh_path, V, F, cache_path=None):
    """Apre una finestra pyvista e raccoglie una lista ordinata di indici di vertice.

    Click sinistro = aggiungi punto più vicino al click ('find_closest_point').
    Tasto 'u'      = annulla ultimo punto.
    Tasto 's'      = salva e chiudi.
    """
    import pyvista as pv

    faces_pv = np.hstack(
        [np.full((len(F), 1), 3, dtype=np.int64), F.astype(np.int64)]
    ).ravel()
    mesh = pv.PolyData(V, faces_pv)

    picked = []
    state = {"finished": False}

    pl = pv.Plotter()
    pl.add_mesh(mesh, color="lightgray", show_edges=False)
    pl.add_axes()

    def _refresh():
        pl.add_text(
            _prompt_text(mesh_path, picked) +
            "\n[click sinistro = aggiungi   u = annulla   s = salva]",
            name="landmark_prompt", position="upper_left",
            font_size=12, color="white",
        )
        if picked:
            pts = V[picked]
            tags = [_letter_for_index(i) for i in range(len(picked))]
            pl.add_point_labels(
                pts, tags, name="landmark_markers", point_size=14, font_size=12,
                text_color="yellow", point_color="red", render_points_as_spheres=True,
                always_visible=True, shape_opacity=0.6,
            )
        else:
            pl.remove_actor("landmark_markers")

    def _on_pick(point):
        vidx = int(mesh.find_closest_point(point))
        picked.append(vidx)
        print(f"  [{_letter_for_index(len(picked) - 1)}] vertice {vidx}")
        _refresh()

    def _undo():
        if picked:
            removed = picked.pop()
            print(f"  annullato: vertice {removed}")
            _refresh()

    def _save():
        state["finished"] = True
        if cache_path is not None:
            save_landmark_cache(cache_path, picked, V)
        print(f"  salvati {len(picked)} punti" +
              (f" in {cache_path}" if cache_path is not None else ""))
        pl.close()

    # picker="cell" (vtkCellPicker, ray/surface intersection) — NOT "point"
    # (vtkPointPicker): on a curved/near-closed surface like an atrium,
    # vtkPointPicker snaps to the nearest rendered point within a screen-space
    # tolerance regardless of occlusion, so it can silently grab a vertex on
    # the far side of the mesh. vtkCellPicker respects depth, so the 3D point
    # it returns is genuinely on the clicked (front) surface.
    if hasattr(pl, "enable_surface_point_picking"):
        pl.enable_surface_point_picking(
            callback=_on_pick, show_message=True, left_clicking=True,
            show_point=False, picker="cell",
        )
    else:
        pl.enable_point_picking(
            callback=_on_pick, use_picker=False, show_message=True,
            left_clicking=True, show_point=False, picker="cell",
        )
    pl.add_key_event("u", _undo)
    pl.add_key_event("s", _save)
    _refresh()
    pl.show()

    if not state["finished"]:
        if not picked:
            raise RuntimeError(
                f"Picking annullato per {mesh_path}: nessun punto raccolto e finestra "
                f"chiusa senza premere 's'. Se la finestra non si è mai aperta "
                f"visibilmente, il server X non è raggiungibile (su sessioni Wayland "
                f"serve XWayland attivo con DISPLAY impostato, es. 'export DISPLAY=:0'; "
                f"verifica con 'echo $DISPLAY' in un terminale del desktop grafico). "
                f"Nessun file di cache scritto."
            )
        raise RuntimeError(
            f"Picking annullato per {mesh_path}: finestra chiusa senza premere 's' "
            f"(raccolti {len(picked)} punti, non salvati). Nessun file di cache scritto."
        )
    if len(picked) == 0:
        raise RuntimeError(f"Nessun landmark selezionato per {mesh_path}.")
    return picked


# ─── Picking interattivo (browser, via trame — non richiede X11/Wayland) ────

def pick_landmarks_browser(mesh_path, V, F, cache_path, port=None, _test_hooks=None):
    """Come ``pick_landmarks``, ma serve la scena in un server web locale.

    Il rendering è lato server (VTK reale, off-screen); il browser riceve
    solo le immagini e inoltra gli eventi mouse, quindi il picking funziona
    esattamente come nella finestra nativa — utile quando non c'è un server
    X raggiungibile (es. sessioni Wayland senza XWayland) ma è disponibile
    un browser qualsiasi (che essendo un'app Wayland nativa non ne ha bisogno).

    Il server gira in un thread separato e la funzione ritorna non appena
    l'utente preme 'Salva' (i punti sono già scritti su disco a quel
    momento): non si aspetta lo shutdown pulito del server web, che con
    alcune combinazioni browser/backend può restare bloccato — in tal caso
    resta comunque una scheda "connection closed" innocua in background.

    ``_test_hooks``, se passato un dict, viene popolato con i callback interni
    (``on_pick``, ``undo``, ``save``) per permettere di simulare un'interazione
    utente nei test automatici, senza un browser reale.
    """
    import asyncio
    import socket
    import threading

    import pyvista as pv
    from pyvista.trame import plotter_ui
    from trame.app import get_server
    from trame.ui.vuetify3 import SinglePageLayout
    from trame.widgets import html, vuetify3

    pv.OFF_SCREEN = True

    faces_pv = np.hstack(
        [np.full((len(F), 1), 3, dtype=np.int64), F.astype(np.int64)]
    ).ravel()
    mesh = pv.PolyData(V, faces_pv)

    picked = []
    done_event = threading.Event()

    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(mesh, color="lightgray", show_edges=False)
    pl.add_axes()

    server = get_server(client_type="vue3")
    state, ctrl = server.state, server.controller
    state.trame__title = f"Landmark picking — {os.path.basename(mesh_path)}"
    state.prompt = _prompt_text(mesh_path, picked)

    def _refresh():
        state.prompt = _prompt_text(mesh_path, picked)
        if picked:
            pts = V[picked]
            tags = [_letter_for_index(i) for i in range(len(picked))]
            pl.add_point_labels(
                pts, tags, name="landmark_markers", point_size=14, font_size=12,
                text_color="yellow", point_color="red", render_points_as_spheres=True,
                always_visible=True, shape_opacity=0.6,
            )
        else:
            pl.remove_actor("landmark_markers")
        ctrl.view_update()

    def _on_pick(point):
        vidx = int(mesh.find_closest_point(point))
        picked.append(vidx)
        print(f"  [{_letter_for_index(len(picked) - 1)}] vertice {vidx}")
        _refresh()

    def _undo():
        if picked:
            removed = picked.pop()
            print(f"  annullato: vertice {removed}")
            _refresh()

    def _save():
        save_landmark_cache(cache_path, picked, V)
        print(f"  salvati {len(picked)} punti in {cache_path}")
        state.prompt = (
            f"Salvato ({len(picked)} punti) in {cache_path}. "
            f"Puoi chiudere questa scheda: il terminale sta continuando."
        )
        ctrl.view_update()
        done_event.set()
        # Best-effort shutdown, scheduled as two independent tasks (not
        # awaited here) so neither blocks the function's return, and so a
        # hung server.stop() doesn't also block closing the VTK context.
        # pl is only ever touched from this event-loop thread — never close
        # it from the caller's (main) thread, VTK objects aren't thread-safe.
        asyncio.create_task(server.stop())

        async def _close_plotter():
            try:
                pl.close()
            except Exception:
                pass

        asyncio.create_task(_close_plotter())

    # picker="cell": see comment in pick_landmarks() — vtkPointPicker can
    # snap through to the far side of a curved/near-closed surface.
    if hasattr(pl, "enable_surface_point_picking"):
        pl.enable_surface_point_picking(
            callback=_on_pick, show_message=False, left_clicking=True,
            show_point=False, picker="cell",
        )
    else:
        pl.enable_point_picking(
            callback=_on_pick, use_picker=False, show_message=False,
            left_clicking=True, show_point=False, picker="cell",
        )

    if _test_hooks is not None:
        _test_hooks.update(on_pick=_on_pick, undo=_undo, save=_save,
                            picked=picked, done_event=done_event, mesh=mesh)

    with SinglePageLayout(server) as layout:
        layout.title.set_text(state.trame__title)
        with layout.toolbar:
            vuetify3.VSpacer()
            vuetify3.VBtn("Annulla", click=_undo, classes="mx-1")
            vuetify3.VBtn("Salva", click=_save, color="primary", classes="mx-1")
        with layout.content:
            html.Div(
                "{{ prompt }}",
                style="padding:8px 16px; font-family:monospace; white-space:pre-line;",
            )
            view = plotter_ui(pl)
            ctrl.view_update = view.update
        layout.footer.hide()

    if port is None:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
        s.close()

    url = f"http://127.0.0.1:{port}/"
    print(f"[landmarks] nessun server X raggiungibile: uso un server browser locale.")
    print(f"[landmarks] apri nel browser: {url}")
    print(f"  click sinistro sulla mesh = aggiungi punto; bottoni Annulla/Salva in alto.")
    print(f"  Il salvataggio scrive subito su disco; se dopo 'Salva' la scheda resta "
          f"su 'connection closed' puoi ignorarla, il terminale sta già proseguendo.")
    print(f"  Se invece chiudi la scheda senza premere 'Salva', interrompi con Ctrl+C "
          f"qui: nessuna cache verrà scritta.")

    server_thread = threading.Thread(
        target=server.start,
        kwargs=dict(port=port, thread=True, open_browser=False,
                    show_connection_info=False),
        daemon=True,
    )
    server_thread.start()

    try:
        done_event.wait()
    except KeyboardInterrupt:
        print(f"\n[landmarks] interrotto (Ctrl+C) per {mesh_path}.")
        raise

    if len(picked) == 0:
        raise RuntimeError(f"Nessun landmark selezionato per {mesh_path}.")
    return picked


# ─── Wrapper cache-aware ────────────────────────────────────────────────────

def get_landmarks_for_mesh(mesh_path, V, F, cache_path=None, repick=False):
    cache_path = cache_path or landmarks_cache_path(mesh_path)

    if not repick and os.path.exists(cache_path):
        idx = load_landmark_cache(cache_path)
        if cache_matches_mesh(cache_path, V, idx):
            print(f"[landmarks] uso cache: {cache_path} ({len(idx)} punti)")
            return idx
        print(
            f"[landmarks] cache stale per {mesh_path}: {cache_path} non combacia "
            f"con la mesh corrente (mesh diversa da quando è stata generata, o "
            f"indici fuori range). Ignoro la cache e ripicco."
        )

    if _has_display():
        print(f"[landmarks] avvio picking interattivo (finestra nativa) per {mesh_path} ...")
        idx = pick_landmarks(mesh_path, V, F, cache_path=cache_path)
    elif _has_trame():
        idx = pick_landmarks_browser(mesh_path, V, F, cache_path)
    else:
        raise RuntimeError(
            f"Landmark cache mancante ({cache_path}) e nessun display disponibile "
            f"per il picking interattivo (né un server X, né 'trame' installato "
            f"per il fallback via browser — 'pip install trame trame-vtk "
            f"trame-vuetify'). Esegui in una sessione con display, oppure fornisci "
            f"il file di cache, prima di rilanciare in batch."
        )

    # entrambi i picker scrivono già la cache al momento del salvataggio;
    # questa è solo una riscrittura di sicurezza idempotente.
    save_landmark_cache(cache_path, idx, V)
    print(f"[landmarks] {len(idx)} landmark pronti per {mesh_path}")
    return idx


def get_matched_landmarks(path_a, Va, Fa, path_b, Vb, Fb, repick=False):
    """Acquisisce (o carica dalla cache) i landmark per due mesh e verifica la corrispondenza.

    Ferma l'esecuzione con errore esplicito se le due liste non hanno la
    stessa lunghezza: un ancoraggio disallineato è peggio di nessun
    ancoraggio.
    """
    idx_a = get_landmarks_for_mesh(path_a, Va, Fa, repick=repick)
    idx_b = get_landmarks_for_mesh(path_b, Vb, Fb, repick=repick)

    check_landmark_spread(Va, idx_a, path_a)
    check_landmark_spread(Vb, idx_b, path_b)

    if len(idx_a) != len(idx_b):
        raise ValueError(
            f"Numero di landmark diverso tra le due mesh: "
            f"{path_a} ha {len(idx_a)} punti, {path_b} ha {len(idx_b)} punti. "
            f"Le liste devono corrispondere 1:1 per struttura anatomica "
            f"(stesso ordine). Ripeti il picking con --repick su entrambe le mesh."
        )
    return idx_a, idx_b
