#!/usr/bin/env python3
"""
Respuesta a Hansen et al. (2025, arXiv:2510.01332), que recontaron la muestra joven de Dai et al.
sin los planetas hallados por TTV y obtuvieron 3/9 pares (33 +- 19 %), a menos de 1 sigma de Kepler.

Reproduzco su recuento y después agrego, de a uno, los datos publicados después:
  1. Período de V1298 Tau e resuelto por Livingston et al. (2026, Nature): 48,68 d -> b/e = 2,017 (2:1).
     Hansen usó la efeméride de Feinstein et al. (2022), con la que b/e queda lejos de 2:1.
  2. HIP 67522 b/c (planeta c publicado en 2024; Barber et al.).
Para cada variante: tasa cruda y efecto ajustado contra Kepler maduro (girocronología), con
permutación por sistema. Ventanas: la de Dai (-1,5 % / +3 %, solo 1er orden, como Hansen) y la mía.

Uso: python3 -I hansen_recount.py
"""
import csv
import math
import random
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
import causal_v3_young as V3  # noqa: E402
from causal_v1_catalog import cerca, PRIMER_ORDEN, VENTANA  # noqa: E402

VENTANA_DAI = (-0.015, 0.03)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    rng = random.Random(20261126)
    filas, _ = V3.cargar()
    maduros = [d for d in filas if d["grupo"] == "maduro_kepler"]
    base_tr = [d for d in filas if d["grupo"] == "muy_joven_tr"]  # pares solo con planetas en tránsito

    # AU Mic e (detectado por velocidad radial; el archivo lo marca controvertido, por eso no estaba)
    arch = [r for r in csv.DictReader(open(BASE / "data" / "young_systems_archive.csv", newline=""))
            if r["hostname"] == "AU Mic"]
    au = {r["pl_name"]: r for r in arch}
    print("AU Mic en el archivo:", {k: (round(float(v["pl_orbper"]), 3), "controvertido" if v["pl_controv_flag"] == "1" else "ok",
                                         v["discoverymethod"], "tránsito" if v["tran_flag"] == "1" else "no transita")
                                     for k, v in au.items()})
    au_c, au_e = num(au["AU Mic c"]["pl_orbper"]), num(au["AU Mic e"]["pl_orbper"])
    plantilla = next(d for d in base_tr if d["sistema"] == "AU Mic")
    r = au_e / au_c
    par_ce = dict(plantilla, cociente=round(r, 5), p_int=au_c, p_ext=au_e, log_p_interno=round(math.log10(au_c), 4),
                  resonante=cerca(r, PRIMER_ORDEN, VENTANA), resonante_dai=cerca(r, PRIMER_ORDEN, VENTANA_DAI))

    def es_v1298_be(d):
        return d["sistema"] == "V1298 Tau" and abs(d["p_int"] - 24.14) < 0.1

    def variante(con_livingston, con_hip, con_au_e=True):
        out = []
        for d in base_tr:
            if d["sistema"] == "HIP 67522" and not con_hip:
                continue
            d = dict(d)
            if es_v1298_be(d) and not con_livingston:
                d["resonante"] = d["resonante_dai"] = 0  # efeméride previa: b/e lejos de 2:1 (Hansen)
            out.append(d)
        if con_au_e:
            out.append(dict(par_ce))
        return out

    casos = [
        ("A. Recuento de Hansen et al. (efeméride vieja de V1298 Tau e, sin HIP 67522)", False, False),
        ("B. A + período de V1298 Tau e de Livingston et al. 2026", True, False),
        ("C. B + HIP 67522 (muestra actual completa, solo tránsito + AU Mic e)", True, True),
    ]
    lin = []
    def say(s=""):
        print(s); lin.append(s)
    say(); say("Kepler maduro (girocronología): ventana de Dai "
               f"{sum(d['resonante_dai'] for d in maduros)}/{len(maduros)} = "
               f"{100*sum(d['resonante_dai'] for d in maduros)/len(maduros):.1f} %, mi ventana "
               f"{sum(d['resonante'] for d in maduros)}/{len(maduros)} = {100*sum(d['resonante'] for d in maduros)/len(maduros):.1f} %")
    for etiqueta, liv, hip in casos:
        jov = variante(liv, hip)
        say(); say(etiqueta)
        for y, nombre in (("resonante_dai", "ventana de Dai (-1,5/+3 %)"), ("resonante", "mi ventana (-0,5/+2,5 %)")):
            k = sum(d[y] for d in jov)
            datos = [dict(d, A=1) for d in jov] + [dict(d, A=0) for d in maduros]
            p1, p0, _ = V3.g_comp(datos, y)
            ci, _ = V3.boot(datos, y, rng, 500)
            _, p = V3.permutacion(datos, y, rng, 5000)
            say(f"   {nombre:<28} jóvenes {k}/{len(jov)} = {100*k/len(jov):4.0f} % | ajustado {100*(p1-p0):+5.1f} pp "
                f"[{100*ci[0]:+.1f}; {100*ci[1]:+.1f}] | permutación p = {p:.4f}")
    (BASE / "results" / "summary_hansen_recount.txt").write_text("\n".join(lin) + "\n")


if __name__ == "__main__":
    main()
