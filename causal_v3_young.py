#!/usr/bin/env python3
"""
Versión 3: los sistemas MUY jóvenes (< 100 Myr), donde la señal debería ser fuerte.

Uso:
    python3 -I causal_v3_young.py

Lee de data/: young_systems_archive.csv, pscomppars_multi.csv,
kepler_koi_cumulative.csv y bouma2024_koi_ages.csv.

MUESTRAS (fijadas antes de ver los resultados de esta versión)
  Muy jóvenes (< 100 Myr): los 7 sistemas multiplanetarios con planetas en tránsito cuya edad viene
    de pertenencia a cúmulo/asociación o de varios indicadores coincidentes (lista de Dai et al. 2024
    actualizada a nov-2026: entra HIP 67522, cuyo planeta c se publicó en 2024; sale TOI-1227, cuyo
    segundo planeta solo se infiere por TTV y no figura en el archivo).
  Adolescentes (0,1-1 Gyr): miembros de cúmulos citados por Dai et al. + Kepler con girocronología.
  Maduros A (>= 1 Gyr): Kepler con girocronología (Bouma et al. 2024, criterio de calidad de los autores).
  Maduros B (>= 1 Gyr): sistemas de TESS con edad del archivo. Edades menos precisas, pero MISMO
    instrumento que los jóvenes: si A y B dan lo mismo, la facilidad no explica el efecto.
  Soporte común: solo pares con P_externo <= 50 días (TESS y K2 no ven períodos largos).

DAG Y AJUSTE: igual que la v2 (Teff, [Fe/H], log P interno). Si falta [Fe/H] se imputa 0 (solar),
  y se cuenta cuántos casos.
RESULTADO: -0,5 % / +2,5 % de 2:1, 3:2, 4:3, 5:4, 6:5. Sensibilidad: ventana de Dai (-1,5 % / +3 %).
CONTROLES: placebo limpio (ventanas en los huecos entre resonancias, ver v1); diseño local
  (pares a +-10 % de una resonancia); solo planetas en tránsito (los descubiertos por TTV se encuentran
  justamente PORQUE están en resonancia: sesgo de selección); sin V1298 Tau; dejando afuera
  cada sistema joven de a uno; test de permutación a nivel sistema.
"""
import csv
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from causal_v1_catalog import logit_fit, sigm, e_value, cerca, en_placebo, PRIMER_ORDEN, VENTANA, VENTANAS_PLACEBO  # noqa: E402
import causal_v2_gyro as G  # noqa: E402

D = BASE / "data"
P_EXT_MAX = 50.0
VENTANA_DAI = (-0.015, 0.03)
N_BOOT = 1000
SEMILLA = 20261126

MUY_JOVENES = {  # edad en Myr y origen
    "HIP 67522": (17, "Sco-Cen (UCL)"), "AU Mic": (19, "grupo de beta Pic"),
    "V1298 Tau": (20, "Taurus-Auriga"), "HD 109833": (27, "LCC"),
    "TOI-942": (53, "varios indicadores"), "TIC 434398831": (61, "Theia 116"),
    "TOI-6109": (75, "Alfa Persei / Theia 133"),
}
ADOLESCENTES_CUMULO = {
    "TOI-451": (125, "Corriente Psc-Eri"), "TOI-2076": (210, "comóvil con TOI-1807"),
    "HD 63433": (414, "Osa Mayor"), "K2-136": (650, "Híades"), "K2-264": (650, "Pesebre"),
    "TOI-1136": (700, "Dai et al. 2023"),
}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def pares_de(planetas, base):
    """planetas: lista de (periodo, en_transito). Devuelve filas de pares vecinos."""
    planetas = sorted(p for p in planetas if p[0])
    out = []
    for (p1, t1), (p2, t2) in zip(planetas, planetas[1:]):
        r = p2 / p1
        if not (1 < r <= 5.5) or p2 > P_EXT_MAX:
            continue
        f = dict(base, cociente=round(r, 5), p_int=p1, p_ext=p2, ambos_transitan=int(t1 and t2),
                 log_p_interno=round(math.log10(p1), 4),
                 resonante=cerca(r, PRIMER_ORDEN, VENTANA), resonante_dai=cerca(r, PRIMER_ORDEN, VENTANA_DAI),
                 placebo_limpio=en_placebo(r), vecindario=int(any(0.90 <= r / c <= 1.10 for c in PRIMER_ORDEN)))
        out.append(f)
    return out


def cargar():
    filas, imputados = [], 0
    arch = defaultdict(list)
    for r in csv.DictReader(open(D / "young_systems_archive.csv", newline="")):
        if r["pl_controv_flag"] != "1":
            arch[r["hostname"]].append(r)
    multi = defaultdict(list)
    for r in csv.DictReader(open(D / "pscomppars_multi.csv", newline="")):
        if r["pl_controv_flag"] != "1":
            multi[r["hostname"]].append(r)

    def sistema(host, edad, origen, grupo, fuente):
        nonlocal imputados
        pls = arch.get(host) or multi.get(host)
        if not pls:
            print(f"   AVISO: {host} no está en los archivos descargados")
            return []
        h = pls[0]
        met = num(h["st_met"])
        if met is None:
            imputados += 1
            met = 0.0
        base = dict(sistema=host, grupo=grupo, edad_myr=edad, origen_edad=origen, fuente=fuente,
                    teff=num(h["st_teff"]), metalicidad=met, n_planetas=len(pls))
        plan = [(num(p["pl_orbper"]), p.get("tran_flag", "1") == "1") for p in pls]
        todos = pares_de(plan, base)
        # versión "solo tránsito": se rearman los pares SIN los planetas hallados por TTV
        solo_tr = pares_de([x for x in plan if x[1]], dict(base, grupo=base["grupo"] + "_tr"))
        return todos + solo_tr

    for host, (e, o) in MUY_JOVENES.items():
        filas += sistema(host, e, o, "muy_joven", "TESS/K2")
    for host, (e, o) in ADOLESCENTES_CUMULO.items():
        filas += sistema(host, e, o, "adolescente", "TESS/K2")

    # Kepler con girocronología (v2), recortado al mismo soporte de períodos
    for d in G.cargar(str(D / "kepler_koi_cumulative.csv"), str(D / "bouma2024_koi_ages.csv")):
        if not d["calidad_ok"] or d["teff"] is None or d["metalicidad"] is None:
            continue
        p1 = 10 ** d["log_p_interno"]
        p2 = p1 * d["cociente"]
        if p2 > P_EXT_MAX or d["kepid"] in ("11773022",):  # 11773022 = Kepler-51, ya contado
            continue
        grupo = "adolescente" if d["joven"] else "maduro_kepler"
        filas.append(dict(sistema=f"KIC {d['kepid']}", grupo=grupo, edad_myr=d["edad_myr"],
                          origen_edad="girocronología (Bouma+2024)", fuente="Kepler", teff=d["teff"],
                          metalicidad=d["metalicidad"], n_planetas=None, cociente=d["cociente"], p_int=p1,
                          p_ext=p2, ambos_transitan=1, log_p_interno=d["log_p_interno"],
                          resonante=d["resonante"], resonante_dai=cerca(d["cociente"], PRIMER_ORDEN, VENTANA_DAI),
                          placebo_limpio=d["placebo_limpio"], vecindario=d["vecindario"]))

    # TESS maduros (misma facilidad que los jóvenes)
    for host, pls in multi.items():
        h = pls[0]
        edad = num(h["st_age"])
        tess = "TESS" in h["disc_facility"] or "Transiting Exoplanet Survey" in h["disc_facility"]
        if not tess or edad is None or edad < 1.0 or host in MUY_JOVENES or host in ADOLESCENTES_CUMULO:
            continue
        met = num(h["st_met"])
        if met is None or num(h["st_teff"]) is None:
            continue
        base = dict(sistema=host, grupo="maduro_tess", edad_myr=edad * 1000, origen_edad="archivo NASA",
                    fuente="TESS", teff=num(h["st_teff"]), metalicidad=met, n_planetas=len(pls))
        filas += pares_de([(num(p["pl_orbper"]), True) for p in pls], base)
    return filas, imputados


def diseño(d, a):
    return [1.0, float(a), (d["teff"] - 5500) / 500, d["metalicidad"], d["log_p_interno"] - 1.0]


def g_comp(datos, y):
    X = [diseño(d, d["A"]) for d in datos]
    beta = logit_fit(X, [d[y] for d in datos])
    pred = lambda d, a: sigm(sum(b * x for b, x in zip(beta, diseño(d, a))))
    p1 = sum(pred(d, 1) for d in datos) / len(datos)
    p0 = sum(pred(d, 0) for d in datos) / len(datos)
    jov = [d for d in datos if d["A"]]
    att = (sum(d[y] for d in jov) - sum(pred(d, 0) for d in jov)) / len(jov)
    return p1, p0, att


def boot(datos, y, rng, n):
    por = defaultdict(list)
    for d in datos:
        por[d["sistema"]].append(d)
    tratados = [s for s in por if por[s][0]["A"]]
    control = [s for s in por if not por[s][0]["A"]]
    difs, atts = [], []
    for _ in range(n):  # bootstrap estratificado: se remuestrean jóvenes y maduros por separado
        m = [d for s in [rng.choice(tratados) for _ in tratados] + [rng.choice(control) for _ in control]
             for d in por[s]]
        p1, p0, att = g_comp(m, y)
        difs.append(p1 - p0); atts.append(att)
    q = lambda v, a: sorted(v)[int(a * (len(v) - 1))]
    return (q(difs, .025), q(difs, .975)), (q(atts, .025), q(atts, .975))


def permutacion(datos, y, rng, n=5000):
    """Test exacto aproximado: se permutan las etiquetas joven/maduro entre SISTEMAS."""
    por = defaultdict(list)
    for d in datos:
        por[d["sistema"]].append(d)
    sist = list(por)
    nt = sum(1 for s in sist if por[s][0]["A"])
    def dif(trat):
        t = [d[y] for s in sist if s in trat for d in por[s]]
        c = [d[y] for s in sist if s not in trat for d in por[s]]
        return sum(t) / len(t) - sum(c) / len(c)
    obs = dif({s for s in sist if por[s][0]["A"]})
    cnt = sum(dif(set(rng.sample(sist, nt))) >= obs for _ in range(n))
    return obs, (1 + cnt) / (1 + n)


def main():
    rng = random.Random(SEMILLA)
    lin = []
    def say(s=""):
        print(s); lin.append(s)

    filas, imputados = cargar()
    salida = BASE / "results"
    campos = list(dict.fromkeys(k for f in filas for k in f))
    with open(salida / "pairs_v3_young.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos); w.writeheader(); w.writerows(filas)

    say(f"Soporte común: P_externo <= {P_EXT_MAX:.0f} d  |  [Fe/H] imputada como solar en {imputados} sistemas")
    say(); say("PELDAÑO 1 · ASOCIACIÓN")
    for g in ("muy_joven", "muy_joven_tr", "adolescente", "maduro_kepler", "maduro_tess"):
        s = [d for d in filas if d["grupo"] == g]
        k, kd = sum(d["resonante"] for d in s), sum(d["resonante_dai"] for d in s)
        say(f"   {g:<14} {len({d['sistema'] for d in s}):3d} sistemas, {len(s):3d} pares | en resonancia "
            f"{k:3d} ({100*k/len(s):4.1f} %) | ventana de Dai {kd:3d} ({100*kd/len(s):4.1f} %)")
    say("   Pares muy jóvenes:")
    for d in sorted((d for d in filas if d["grupo"] == "muy_joven"), key=lambda d: (d["edad_myr"], d["p_int"])):
        marca = "RESONANTE" if d["resonante"] else ("(solo ventana Dai)" if d["resonante_dai"] else "")
        say(f"     {d['sistema']:<14} {d['edad_myr']:>3.0f} Myr  {d['p_int']:7.3f} -> {d['p_ext']:7.3f} d  "
            f"r = {d['cociente']:.4f}  {'' if d['ambos_transitan'] else '[incluye planeta por TTV]'} {marca}")

    def analisis(etiqueta, datos, y="resonante", n=N_BOOT, perm=False):
        p1, p0, att = g_comp(datos, y)
        ci, ci_att = boot(datos, y, rng, n)
        extra = ""
        if perm:
            obs, p = permutacion(datos, y, rng)
            extra = f" | permutación por sistema p = {p:.4f}"
        say(f"   {etiqueta}")
        say(f"      {sum(d['A'] for d in datos)} pares tratados / {len(datos)} | do(joven) {100*p1:.1f} % vs "
            f"do(maduro) {100*p0:.1f} % | diferencia {100*(p1-p0):+.1f} pp [{100*ci[0]:+.1f}; {100*ci[1]:+.1f}] | "
            f"razón {p1/p0:.2f} | E-value {e_value(p1/p0):.2f}{extra}")
        return p1, p0, att, ci_att

    def con_A(grupo_t, grupo_c, filtro=lambda d: True):
        return [dict(d, A=int(d["grupo"] == grupo_t)) for d in filas
                if d["grupo"] in (grupo_t, grupo_c) and filtro(d)]

    principal = con_A("muy_joven", "maduro_kepler")
    say(); say("PELDAÑO 2 · INTERVENCIÓN: muy jóvenes (< 100 Myr) contra maduros")
    p1, p0, att, ci_att = analisis("Principal: contra Kepler maduro (girocronología)", principal, perm=True)
    analisis("Triangulación: contra TESS maduro (mismo instrumento que los jóvenes)",
             con_A("muy_joven", "maduro_tess"), perm=True)

    say(); say("PELDAÑO 3 · CONTRAFÁCTICO (sobre los muy jóvenes, comparación principal)")
    jov = [d for d in principal if d["A"]]
    k = sum(d["resonante"] for d in jov)
    say(f"   {k} de {len(jov)} pares muy jóvenes en resonancia ({100*k/len(jov):.0f} %). Si esos sistemas fueran "
        f"maduros: {100*(k/len(jov)-att):.0f} %. La edad explica {100*att:+.0f} pp "
        f"[{100*ci_att[0]:+.0f}; {100*ci_att[1]:+.0f}]")

    say(); say("CONTROL NEGATIVO (ventanas en los huecos entre resonancias: "
               + ", ".join(f"{lo:.3f}-{hi:.3f}" for lo, hi in VENTANAS_PLACEBO) + ")")
    analisis("Placebo limpio", principal, "placebo_limpio", n=N_BOOT // 2, perm=True)
    analisis("Diseño local (solo pares a +-10 % de una resonancia)",
             con_A("muy_joven", "maduro_kepler", lambda d: d["vecindario"]), n=N_BOOT // 2, perm=True)

    say(); say("SENSIBILIDAD")
    analisis("Ventana de Dai (-1,5 % / +3 %)", principal, "resonante_dai", n=N_BOOT // 2)
    analisis("Solo planetas en tránsito (pares rearmados sin los hallados por TTV)",
             con_A("muy_joven_tr", "maduro_kepler"), n=N_BOOT // 2, perm=True)
    analisis("Sin V1298 Tau", con_A("muy_joven", "maduro_kepler", lambda d: d["sistema"] != "V1298 Tau"),
             n=N_BOOT // 2, perm=True)
    difs = []
    for s in MUY_JOVENES:
        sub = con_A("muy_joven", "maduro_kepler", lambda d, s=s: d["sistema"] != s)
        a, b, _ = g_comp(sub, "resonante")
        difs.append((100 * (a - b), s))
    say(f"   Dejando afuera un sistema joven por vez: diferencia entre {min(difs)[0]:+.1f} pp "
        f"(sin {min(difs)[1]}) y {max(difs)[0]:+.1f} pp (sin {max(difs)[1]})")

    say(); say("ADOLESCENTES (0,1-1 Gyr: cúmulos + girocronología) contra Kepler maduro")
    analisis("Adolescentes", con_A("adolescente", "maduro_kepler"), n=N_BOOT // 2, perm=True)

    (salida / "summary_v3_young.txt").write_text("\n".join(lin) + "\n")


if __name__ == "__main__":
    main()
