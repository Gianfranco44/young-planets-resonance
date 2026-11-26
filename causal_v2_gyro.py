#!/usr/bin/env python3
"""
Versión 2: ¿la edad causa la pérdida de resonancia? Ahora con edades por GIROCRONOLOGÍA.

Uso:
    python3 -I causal_v2_gyro.py data/kepler_koi_cumulative.csv data/bouma2024_koi_ages.csv

Datos:
  * Períodos: NASA Exoplanet Archive, tabla `cumulative` de Kepler (KOIs que el pipeline no
    descartó como falsos positivos; se excluye además lo que el archivo marca FALSE POSITIVE).
  * Edades: Bouma et al. 2024 ("Ages of Stars and Planets in the Kepler Field Younger Than Four
    Billion Years", arXiv:2410.06246), tabla completa del código fuente en arXiv. Son edades por
    rotación estelar calibradas con cúmulos abiertos de 0,04 a 2,5 Gyr.

DECISIONES FIJADAS ANTES DE VER LOS RESULTADOS DE ESTA VERSIÓN
  Unidad:       par de KOIs vecinos en período (1 < r <= 5,5).
  Tratamiento:  joven = edad por rotación (mediana) < 1000 Myr; maduro = mediana >= 1000 Myr
                o cota inferior "> 4000 Myr".
  Calidad:      criterio de los propios autores: ningún bit 0-10 de flag_gyro_quality levantado.
  Resultado:    -0,5 % / +2,5 % de 2:1, 3:2, 4:3, 5:4, 6:5 (igual que antes).
  DAG:          todo es Kepler, así que la facilidad de descubrimiento deja de ser confusor.
                Ajuste por puerta trasera: Teff (masa estelar), [Fe/H] y log P interno.
                No se ajusta por radio planetario ni multiplicidad (consecuencias de la edad).
  Control negativo: ventanas en los huecos entre resonancias de 1er y 2do orden (ver la v1; la
                primera versión corría las resonancias y caía sobre otras, error corregido).
  Diseño local: pares a +-10 % de una resonancia de 1er orden.
  Sensibilidad: sin filtro de calidad; "joven seguro" (mediana + 1 sigma < 1000 Myr).
"""
import csv
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from causal_v1_catalog import logit_fit, sigm, e_value, cerca, en_placebo, PRIMER_ORDEN, VENTANA  # noqa: E402

R_MIN, R_MAX = 1.0, 5.5
N_BOOT = 400
SEMILLA = 20261126


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def cargar(ruta_koi, ruta_edades):
    edades = {}
    for r in csv.DictReader(open(ruta_edades, newline="")):
        if r["kepid"] in edades:
            continue
        med = num(r["gyro_median"])
        cota = num(r["t_gyro_lowerlimit"].replace(">", "").strip()) if r["t_gyro_lowerlimit"] else None
        flag = int(r["flag_gyro_quality"]) if r["flag_gyro_quality"] else None
        edades[r["kepid"]] = dict(
            edad_myr=med if med is not None else cota,
            es_cota=med is None and cota is not None,
            mas1s=num(r["gyro_+1sigma"]),
            calidad_ok=flag is not None and (flag & 2047) == 0,
            teff=num(r["adopted_Teff"]),
        )
    sist = defaultdict(list)
    for r in csv.DictReader(open(ruta_koi, newline="")):
        if r["koi_disposition"] == "FALSE POSITIVE" or not r["koi_period"]:
            continue
        sist[r["kepid"]].append(r)
    pares = []
    for kepid, kois in sist.items():
        e = edades.get(kepid)
        if not e or e["edad_myr"] is None or len(kois) < 2:
            continue
        kois.sort(key=lambda d: float(d["koi_period"]))
        for a, b in zip(kois, kois[1:]):
            p1, p2 = float(a["koi_period"]), float(b["koi_period"])
            r = p2 / p1
            if not (R_MIN < r <= R_MAX):
                continue
            fila = dict(
                kepid=kepid, par=f"{a['kepoi_name']} / {b['kepoi_name']}",
                nombres=f"{a['kepler_name'] or '-'} / {b['kepler_name'] or '-'}", cociente=round(r, 5),
                edad_myr=e["edad_myr"], edad_es_cota=int(e["es_cota"]),
                joven=int((not e["es_cota"]) and e["edad_myr"] < 1000),
                joven_seguro=int((not e["es_cota"]) and e["mas1s"] is not None and e["edad_myr"] + e["mas1s"] < 1000),
                calidad_ok=int(e["calidad_ok"]), teff=e["teff"] or num(a["koi_steff"]),
                metalicidad=num(a["koi_smet"]), log_p_interno=round(math.log10(p1), 4),
                resonante=cerca(r, PRIMER_ORDEN, VENTANA),
                vecindario=int(any(0.90 <= r / c <= 1.10 for c in PRIMER_ORDEN)),
                placebo_limpio=en_placebo(r),
            )
            pares.append(fila)
    return pares


def diseño(d, a):
    return [1.0, float(a), (d["teff"] - 5500) / 500, d["metalicidad"], d["log_p_interno"] - 1.0]


def g_comp(datos, var_y, var_a="joven"):
    X = [diseño(d, d[var_a]) for d in datos]
    y = [d[var_y] for d in datos]
    beta = logit_fit(X, y)
    pred = lambda d, a: sigm(sum(b * x for b, x in zip(beta, diseño(d, a))))
    p1 = sum(pred(d, 1) for d in datos) / len(datos)
    p0 = sum(pred(d, 0) for d in datos) / len(datos)
    jov = [d for d in datos if d[var_a] == 1]
    att = (sum(d[var_y] for d in jov) - sum(pred(d, 0) for d in jov)) / len(jov)
    return p1, p0, att


def bootstrap(datos, var_y, rng, var_a="joven", n=N_BOOT):
    por = defaultdict(list)
    for d in datos:
        por[d["kepid"]].append(d)
    ks = list(por)
    difs, rrs, atts = [], [], []
    for _ in range(n):
        m = [d for k in (rng.choice(ks) for _ in ks) for d in por[k]]
        if sum(d[var_a] for d in m) < 3:
            continue
        p1, p0, att = g_comp(m, var_y, var_a)
        difs.append(p1 - p0); rrs.append(p1 / p0 if p0 > 0 else float("inf")); atts.append(att)
    q = lambda v, a: sorted(v)[int(a * (len(v) - 1))]
    return dict(dif=(q(difs, .025), q(difs, .975)), rr=(q(rrs, .025), q(rrs, .975)),
                att=(q(atts, .025), q(atts, .975)))


def main():
    base = Path(__file__).resolve().parent
    salida = base / "results"
    rng = random.Random(SEMILLA)
    lin = []
    def say(s=""):
        print(s); lin.append(s)

    pares = cargar(sys.argv[1], sys.argv[2])
    with open(salida / "pairs_v2_gyro.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pares[0])); w.writeheader(); w.writerows(pares)

    todos = [d for d in pares if d["teff"] is not None and d["metalicidad"] is not None]
    buenos = [d for d in todos if d["calidad_ok"]]
    say(f"Pares Kepler con edad por rotación: {len(todos)} (jóvenes {sum(d['joven'] for d in todos)}, "
        f"en {len({d['kepid'] for d in todos if d['joven']})} sistemas)")
    say(f"Con el criterio de calidad de Bouma et al.: {len(buenos)} (jóvenes {sum(d['joven'] for d in buenos)}, "
        f"en {len({d['kepid'] for d in buenos if d['joven']})} sistemas)  <- muestra principal")

    say(); say("PELDAÑO 1 · ASOCIACIÓN (muestra principal)")
    for nombre, cond in (("< 0,1 Gyr", lambda d: not d["edad_es_cota"] and d["edad_myr"] < 100),
                         ("0,1-1 Gyr", lambda d: not d["edad_es_cota"] and 100 <= d["edad_myr"] < 1000),
                         (">= 1 Gyr", lambda d: d["edad_es_cota"] or d["edad_myr"] >= 1000)):
        g = [d for d in buenos if cond(d)]
        if g:
            k = sum(d["resonante"] for d in g)
            say(f"   {nombre:>10}: {k:3d}/{len(g):4d} = {100*k/len(g):5.1f} %")
    for grupo, sel in (("jóvenes", 1), ("maduros", 0)):
        g = [d for d in buenos if d["joven"] == sel]
        say(f"   {grupo:>8}: Teff media {sum(d['teff'] for d in g)/len(g):.0f} K | [Fe/H] media "
            f"{sum(d['metalicidad'] for d in g)/len(g):+.3f} | log P medio {sum(d['log_p_interno'] for d in g)/len(g):.2f}")

    def bloque(etiqueta, datos, var_y="resonante", var_a="joven", n=N_BOOT):
        p1, p0, att = g_comp(datos, var_y, var_a)
        bs = bootstrap(datos, var_y, rng, var_a, n)
        say(f"   {etiqueta}")
        say(f"      n = {len(datos)} ({sum(d[var_a] for d in datos)} jóvenes) | do(joven) {100*p1:.1f} % vs "
            f"do(maduro) {100*p0:.1f} % | diferencia {100*(p1-p0):+.1f} pp [{100*bs['dif'][0]:+.1f}; "
            f"{100*bs['dif'][1]:+.1f}] | razón {p1/p0:.2f} [{bs['rr'][0]:.2f}; {bs['rr'][1]:.2f}] | "
            f"E-value {e_value(p1/p0):.2f}")
        return p1, p0, att, bs

    say(); say("PELDAÑO 2 · INTERVENCIÓN (ajuste por Teff, [Fe/H] y log P; bootstrap por sistema)")
    p1, p0, att, bs = bloque("Principal", buenos)

    say(); say("PELDAÑO 3 · CONTRAFÁCTICO (efecto sobre los jóvenes)")
    jov = [d for d in buenos if d["joven"]]
    k = sum(d["resonante"] for d in jov)
    say(f"   {k} de {len(jov)} pares jóvenes en resonancia ({100*k/len(jov):.1f} %); si fueran maduros: "
        f"{100*(k/len(jov)-att):.1f} % -> la edad explica {100*att:+.1f} pp [{100*bs['att'][0]:+.1f}; {100*bs['att'][1]:+.1f}]")

    say(); say("CONTROL NEGATIVO (ventanas en los huecos entre resonancias)")
    bloque("Placebo limpio", buenos, "placebo_limpio", n=N_BOOT // 2)

    say(); say("DISEÑO LOCAL Y SENSIBILIDAD")
    bloque("Local: solo pares a +-10 % de una resonancia", [d for d in buenos if d["vecindario"]], n=N_BOOT // 2)
    bloque("Sin filtro de calidad", todos, n=N_BOOT // 2)
    bloque("'Joven seguro' (mediana + 1 sigma < 1 Gyr)", buenos, var_a="joven_seguro", n=N_BOOT // 2)

    say(); say("Pares jóvenes en resonancia (muestra principal):")
    for d in sorted((d for d in buenos if d["joven"] and d["resonante"]), key=lambda d: d["edad_myr"]):
        say(f"   {d['nombres']:<32} {d['par']:<22} r = {d['cociente']:.4f}  edad {d['edad_myr']:.0f} Myr")

    (salida / "summary_v2_gyro.txt").write_text("\n".join(lin) + "\n")


if __name__ == "__main__":
    main()
