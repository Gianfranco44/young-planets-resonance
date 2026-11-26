#!/usr/bin/env python3
"""
¿La EDAD de un sistema planetario CAUSA que sus planetas salgan de resonancia?

Análisis siguiendo la escalera de Pearl (el enfoque de KAIROS), en Python puro.
Uso:
    python3 -I causal_v1_catalog.py data/pscomppars_multi.csv

Antecedente: Dai et al. 2024 (AJ 168, 239) encuentran pares cerca de resonancia de primer orden
en ~70 % de sistemas < 100 Myr, ~24 % entre 0,1 y 1 Gyr y ~15 % en sistemas maduros.
Acá preguntamos si eso se sostiene como EFECTO CAUSAL una vez ajustados los confusores.

DECISIONES FIJADAS ANTES DE VER LOS RESULTADOS
  Unidad:       par de planetas vecinos en período (1 < r <= 5,5), sin planetas controvertidos.
  Tratamiento:  A = 1 si la estrella tiene edad < 1 Gyr ("joven"), 0 si >= 1 Gyr.
  Resultado:    Y = 1 si el par está a -0,5 % / +2,5 % de una resonancia de primer orden
                (2:1, 3:2, 4:3, 5:4, 6:5), la misma ventana del test de Fibonacci.

  DAG (supuestos causales; "no causes in, no causes out"):
      Facilidad de descubrimiento -> A   (TESS e imagen directa buscan estrellas jóvenes a propósito)
      Facilidad de descubrimiento -> Y   (las TTV detectan pares resonantes; imagen ve gigantes lejanos)
      Masa estelar                -> A   (la edad se mide distinto según la masa)
      Masa estelar                -> Y   (el disco y la migración dependen de la masa)
      Metalicidad                 -> A   (evolución química de la galaxia: lo viejo es más pobre en metales)
      Metalicidad                 -> Y   (más metales, más planetas grandes, otra migración)
      log P_interno               -> Y   (las mareas rompen resonancias en períodos cortos)
      log P_interno               <- Facilidad (sesgos de detección)
      A -> Y                             (el efecto que buscamos)
      A -> radio del planeta             (los planetas jóvenes están inflados)  => NO se ajusta
      A -> número de planetas            (inestabilidades con el tiempo)       => NO se ajusta (mediador)

  Conjunto de ajuste (puerta trasera): facilidad, masa estelar, metalicidad, log P_interno.
  Estimador:   g-computation con regresión logística; IC por bootstrap POR SISTEMA
               (los pares de un mismo sistema no son independientes).
  Control negativo: el mismo análisis con un "resultado placebo": estar cerca de valores
               corridos un 8 % de cada resonancia (no son resonancias). Ahí la edad no debería
               tener efecto. Si lo tiene, el ajuste está dejando sesgo.
  Sensibilidad: solo edades con error relativo < 50 %; ventana -1 % / +3 %; sin imagen directa;
               E-value (cuánta confusión no medida haría falta para anular el efecto).
"""
import csv
import math
import random
import sys
from collections import defaultdict, Counter
from pathlib import Path

R_MIN, R_MAX = 1.0, 5.5
PRIMER_ORDEN = (2 / 1, 3 / 2, 4 / 3, 5 / 4, 6 / 5)
VENTANA = (-0.005, 0.025)

# CONTROL NEGATIVO (corregido el 26-nov-2026). La primera versión corría cada resonancia un
# porcentaje fijo, y eso caía sobre OTRAS resonancias (4:3 x 1,12 = 1,493, casi 3:2). Ahora el
# placebo son ventanas del mismo ancho (3 %) puestas en los huecos entre las resonancias de 1er y
# 2do orden (con las ventanas de Dai et al.), en el mismo rango 1,2-2,1 que las resonancias de 1er orden.
F1 = [(j + 1) / j for j in range(1, 8)]          # 2:1 ... 8:7
F2 = [(j + 2) / j for j in (1, 3, 5, 7, 9)]      # 3:1, 5:3, 7:5, 9:7, 11:9


def ventanas_placebo(r_min=1.2, r_max=2.1, margen=0.005):
    prohibidas = [(c * 0.985, c * 1.03) for c in F1] + [(c * 0.985, c * 1.015) for c in F2]
    out, x = [], r_min
    while x * (1 + VENTANA[1]) <= r_max:
        lo, hi = x * (1 + VENTANA[0]), x * (1 + VENTANA[1])
        libre = all(hi * (1 + margen) < a or lo * (1 - margen) > b for a, b in prohibidas)
        if libre and (not out or lo > out[-1][1] * (1 + margen)):
            out.append((lo, hi))
        x *= 1.001
    return out


VENTANAS_PLACEBO = ventanas_placebo()


def en_placebo(r):
    return int(any(lo <= r <= hi for lo, hi in VENTANAS_PLACEBO))
EDAD_CORTE = 1.0
N_BOOT = 400
SEMILLA = 20261126


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def facilidad(f, metodo):
    if "Kepler" in f:
        return "Kepler"
    if "K2" in f:
        return "K2"
    if "TESS" in f or "Transiting Exoplanet Survey" in f:
        return "TESS"
    if metodo == "Imaging":
        return "Imagen"
    if metodo == "Radial Velocity":
        return "VR"
    return "Otra"


def cerca(r, centros, ventana):
    return int(any(c * (1 + ventana[0]) <= r <= c * (1 + ventana[1]) for c in centros))


def cargar(ruta):
    sist = defaultdict(list)
    with open(ruta, newline="") as f:
        for fila in csv.DictReader(f):
            if fila["pl_controv_flag"] == "1" or not fila["pl_orbper"]:
                continue
            sist[fila["hostname"]].append(fila)
    pares = []
    for host, pls in sist.items():
        pls.sort(key=lambda d: float(d["pl_orbper"]))
        h = pls[0]
        edad, e1, e2 = num(h["st_age"]), num(h["st_ageerr1"]), num(h["st_ageerr2"])
        err_rel = ((e1 - e2) / 2 / edad) if (edad and e1 is not None and e2 is not None) else None
        for a, b in zip(pls, pls[1:]):
            p1, p2 = float(a["pl_orbper"]), float(b["pl_orbper"])
            r = p2 / p1
            if not (R_MIN < r <= R_MAX):
                continue
            pares.append(dict(
                sistema=host, par=f"{a['pl_name']} / {b['pl_name']}", cociente=round(r, 5),
                edad_gyr=edad, edad_err_rel=None if err_rel is None else round(err_rel, 3),
                joven=None if edad is None else int(edad < EDAD_CORTE),
                masa_estelar=num(h["st_mass"]), metalicidad=num(h["st_met"]),
                log_p_interno=round(math.log10(p1), 4),
                facilidad=facilidad(h["disc_facility"], h["discoverymethod"]),
                n_planetas=int(h["sy_pnum"]),
                radio_interno=num(a["pl_rade"]), radio_externo=num(b["pl_rade"]),
                resonante=cerca(r, PRIMER_ORDEN, VENTANA),
                resonante_ancha=cerca(r, PRIMER_ORDEN, (-0.01, 0.03)),
                placebo_limpio=en_placebo(r),
                vecindario=int(any(0.90 <= r / c <= 1.10 for c in PRIMER_ORDEN)),
            ))
    return pares


# ------------------------------------------------------------ regresión logística (IRLS)
def resolver(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda i: abs(M[i][c]))
        M[c], M[piv] = M[piv], M[c]
        if abs(M[c][c]) < 1e-12:
            M[c][c] = 1e-12
        for i in range(n):
            if i != c:
                f = M[i][c] / M[c][c]
                if f:
                    for j in range(c, n + 1):
                        M[i][j] -= f * M[c][j]
    return [M[i][n] / M[i][i] for i in range(n)]


def sigm(z):
    return 1 / (1 + math.exp(-max(min(z, 30), -30)))


def logit_fit(X, y, ridge=1e-3, it=25):
    p = len(X[0])
    beta = [0.0] * p
    for _ in range(it):
        H = [[0.0] * p for _ in range(p)]
        g = [0.0] * p
        for xi, yi in zip(X, y):
            mu = sigm(sum(b * x for b, x in zip(beta, xi)))
            w = mu * (1 - mu)
            r = yi - mu
            for j in range(p):
                g[j] += xi[j] * r
                wj = w * xi[j]
                Hj = H[j]
                for k in range(j, p):
                    Hj[k] += wj * xi[k]
        for j in range(p):
            g[j] -= ridge * beta[j]
            H[j][j] += ridge
            for k in range(j):
                H[j][k] = H[k][j]
        paso = resolver(H, g)
        beta = [b + s for b, s in zip(beta, paso)]
        if max(abs(s) for s in paso) < 1e-7:
            break
    return beta


FACS = ("Kepler", "K2", "TESS", "VR", "Imagen", "Otra")


def diseño(d, a, ajustar=True):
    fila = [1.0, float(a)]
    if ajustar:
        fila += [float(d["facilidad"] == f) for f in FACS[1:]]
        fila += [d["masa_estelar"] - 1.0, d["metalicidad"], d["log_p_interno"] - 1.0]
    return fila


def g_comp(datos, var_y, ajustar=True):
    """Devuelve (P(Y|do(joven)), P(Y|do(maduro)), ATT en pp)."""
    X = [diseño(d, d["joven"], ajustar) for d in datos]
    y = [d[var_y] for d in datos]
    beta = logit_fit(X, y)
    pred = lambda d, a: sigm(sum(b * x for b, x in zip(beta, diseño(d, a, ajustar))))
    p1 = sum(pred(d, 1) for d in datos) / len(datos)
    p0 = sum(pred(d, 0) for d in datos) / len(datos)
    jovenes = [d for d in datos if d["joven"] == 1]
    att = (sum(d[var_y] for d in jovenes) - sum(pred(d, 0) for d in jovenes)) / len(jovenes)
    return p1, p0, att


def bootstrap(datos, var_y, rng, ajustar=True, n=N_BOOT):
    por_sist = defaultdict(list)
    for d in datos:
        por_sist[d["sistema"]].append(d)
    sistemas = list(por_sist)
    difs, rrs, atts = [], [], []
    for _ in range(n):
        muestra = [d for s in (rng.choice(sistemas) for _ in sistemas) for d in por_sist[s]]
        if sum(d["joven"] for d in muestra) < 5:
            continue
        p1, p0, att = g_comp(muestra, var_y, ajustar)
        difs.append(p1 - p0); rrs.append(p1 / p0 if p0 > 0 else float("inf")); atts.append(att)
    q = lambda v, a: sorted(v)[int(a * (len(v) - 1))]
    return dict(dif=(q(difs, .025), q(difs, .975)), rr=(q(rrs, .025), q(rrs, .975)),
                att=(q(atts, .025), q(atts, .975)), n=len(difs))


def e_value(rr):
    rr = rr if rr >= 1 else 1 / rr
    return rr + math.sqrt(rr * (rr - 1))


# ------------------------------------------------------------ informe
def main():
    ruta = Path(sys.argv[1])
    base = Path(__file__).resolve().parent
    salida = base / "results"
    salida.mkdir(exist_ok=True)
    rng = random.Random(SEMILLA)
    lin = []
    def say(s=""):
        print(s); lin.append(s)

    pares = cargar(ruta)
    with open(salida / "pairs_v1_catalog.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pares[0]))
        w.writeheader(); w.writerows(pares)

    completos = [d for d in pares if d["joven"] is not None and d["masa_estelar"] is not None
                 and d["metalicidad"] is not None]
    say(f"Pares en el rango: {len(pares)}  |  con edad, masa y metalicidad: {len(completos)} "
        f"(jóvenes: {sum(d['joven'] for d in completos)}, sistemas jóvenes: "
        f"{len({d['sistema'] for d in completos if d['joven']})})")

    say(); say("PELDAÑO 1 · ASOCIACIÓN (lo que se ve, sin ajustar)")
    for nombre, cond in (("< 0,1 Gyr", lambda e: e < 0.1), ("0,1-1 Gyr", lambda e: 0.1 <= e < 1),
                         (">= 1 Gyr", lambda e: e >= 1)):
        g = [d for d in pares if d["edad_gyr"] is not None and cond(d["edad_gyr"])]
        k = sum(d["resonante"] for d in g)
        say(f"   {nombre:>10}: {k:3d}/{len(g):4d} pares cerca de resonancia de 1er orden = {k/len(g)*100:5.1f} %")
    say("   (Dai et al. 2024: ~70 %, ~24 %, ~15 %, con su propia selección y ventana)")

    say(); say("   ¿Jóvenes y maduros son comparables? (positividad)")
    for grupo, sel in (("jóvenes", 1), ("maduros", 0)):
        g = [d for d in completos if d["joven"] == sel]
        c = Counter(d["facilidad"] for d in g)
        say(f"   {grupo:>8}: " + ", ".join(f"{f} {c.get(f, 0)}" for f in FACS)
            + f" | masa media {sum(d['masa_estelar'] for d in g)/len(g):.2f} Msol"
            + f" | [Fe/H] media {sum(d['metalicidad'] for d in g)/len(g):+.3f}"
            + f" | log P medio {sum(d['log_p_interno'] for d in g)/len(g):.2f}")

    say(); say("PELDAÑO 2 · INTERVENCIÓN: P(resonante | do(edad))  [g-computation, IC 95 % bootstrap por sistema]")

    def bloque(etiqueta, datos, var_y="resonante", ajustar=True, n=N_BOOT):
        p1, p0, att = g_comp(datos, var_y, ajustar)
        bs = bootstrap(datos, var_y, rng, ajustar, n)
        rr = p1 / p0
        say(f"   {etiqueta}")
        say(f"      n = {len(datos)} pares ({sum(d['joven'] for d in datos)} jóvenes)  |  "
            f"do(joven) {p1*100:.1f} %  vs  do(maduro) {p0*100:.1f} %")
        say(f"      diferencia {100*(p1-p0):+.1f} pp [{100*bs['dif'][0]:+.1f}, {100*bs['dif'][1]:+.1f}]  |  "
            f"razón {rr:.2f} [{bs['rr'][0]:.2f}, {bs['rr'][1]:.2f}]  |  E-value {e_value(rr):.2f}")
        return p1, p0, att, bs

    bloque("Sin ajustar (solo asociación, para comparar)", completos, ajustar=False)
    p1, p0, att, bs = bloque("AJUSTADO por facilidad, masa, metalicidad y período  <- estimación principal", completos)

    say(); say("PELDAÑO 3 · CONTRAFÁCTICO: efecto sobre los jóvenes (ATT)")
    jov = [d for d in completos if d["joven"]]
    k = sum(d["resonante"] for d in jov)
    say(f"   De {len(jov)} pares jóvenes, {k} están en resonancia ({k/len(jov)*100:.1f} %).")
    say(f"   Si esos mismos sistemas fueran maduros, estarían en resonancia {100*(k/len(jov)-att):.1f} % "
        f"→ la edad explica {100*att:+.1f} pp [{100*bs['att'][0]:+.1f}, {100*bs['att'][1]:+.1f}]")

    say(); say("CONTROL NEGATIVO (ventanas en los huecos entre resonancias: "
               + ", ".join(f"{lo:.3f}-{hi:.3f}" for lo, hi in VENTANAS_PLACEBO) + ")")
    bloque("Ajustado, placebo limpio", completos, "placebo_limpio", n=N_BOOT // 2)

    say(); say("SENSIBILIDAD")
    buenas = [d for d in completos if d["edad_err_rel"] is not None and d["edad_err_rel"] < 0.5]
    bloque("Solo edades con error relativo < 50 %", buenas, n=N_BOOT // 2)
    bloque("Ventana más ancha (-1 % / +3 %)", completos, "resonante_ancha", n=N_BOOT // 2)
    bloque("Sin sistemas de imagen directa", [d for d in completos if d["facilidad"] != "Imagen"], n=N_BOOT // 2)
    bloque("Solo Kepler + K2 + TESS (tránsitos)", [d for d in completos if d["facilidad"] in ("Kepler", "K2", "TESS")], n=N_BOOT // 2)

    say(); say("DISTRIBUCIÓN GENERAL Y DISEÑO LOCAL")
    say("   ¿Jóvenes y maduros reparten distinto los cocientes en general?")
    for grupo, sel in (("jóvenes", 1), ("maduros", 0)):
        g = sorted(d["cociente"] for d in completos if d["joven"] == sel)
        q = lambda a: g[int(a * (len(g) - 1))]
        say(f"   {grupo:>8}: mediana {q(.5):.2f} | cuartiles {q(.25):.2f}-{q(.75):.2f} | "
            f"% con r < 2,2: {100*sum(x < 2.2 for x in g)/len(g):.0f} %")
    local = [d for d in completos if d["vecindario"]]
    say("   Diseño local: solo pares a +-10 % de una resonancia de 1er orden; la pregunta pasa a ser")
    say("   'dado que el par está cerca, ¿cae justo en la resonancia?'. Así se neutraliza que jóvenes")
    say("   y maduros repartan distinto los cocientes a gran escala.")
    bloque("Ajustado, diseño local", local, n=N_BOOT // 2)
    bloque("Ajustado, diseño local, solo edades con error < 50 %",
           [d for d in local if d["edad_err_rel"] is not None and d["edad_err_rel"] < 0.5], n=N_BOOT // 2)

    (salida / "summary_v1_catalog.txt").write_text("\n".join(lin) + "\n")


if __name__ == "__main__":
    main()
