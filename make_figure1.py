#!/usr/bin/env python3
"""Figure 1 of the Research Note (English, no title: the caption carries it).
Usage: python3 -I make_figure1.py   (reads results/pairs_v3_young.csv)"""
import csv
from math import sqrt
from pathlib import Path

BASE = Path(__file__).resolve().parent
AZUL, BERMELLON = "#0072B2", "#D55E00"
TINTA, TINTA2, TENUE, GRILLA = "#1a1a1a", "#555555", "#8a8a8a", "#e6e6e6"
FUENTE = "'Helvetica Neue', Helvetica, sans-serif"

filas = list(csv.DictReader(open(BASE / "results" / "pairs_v3_young.csv", newline="")))
GRUPOS = [("muy_joven", "&lt; 100 Myr", "7 cluster or association systems"),
          ("adolescente", "0.1–1 Gyr", "cluster members and gyrochronology"),
          ("maduro_kepler", "&gt; 1 Gyr, Kepler", "gyrochronology (Bouma et al. 2024)"),
          ("maduro_tess", "&gt; 1 Gyr, TESS", "same instrument as most young systems")]


def wilson(k, n, z=1.96):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return p, max(c - h, 0), min(c + h, 1)


W, H = 1000, 440
x0, x1, y0, y1 = 380, 940, 30, 350
X = lambda v: x0 + v * (x1 - x0)
s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
     f'font-family="{FUENTE}"><rect width="{W}" height="{H}" fill="#ffffff"/>']


def t(x, y, txt, size=12, col=TINTA, anchor="start", weight="normal"):
    s.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{col}" text-anchor="{anchor}" '
             f'font-weight="{weight}">{txt}</text>')


for v in (0, .2, .4, .6, .8, 1):
    s.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{y0}" y2="{y1}" stroke="{GRILLA}"/>')
    t(X(v), y1 + 18, f"{int(v * 100)}%", 12, TENUE, "middle")
t((x0 + x1) / 2, y1 + 38, "fraction of neighboring pairs", 12, TENUE, "middle")

fila_h = (y1 - y0) / len(GRUPOS)
for i, (g, nombre, sub) in enumerate(GRUPOS):
    sel = [f for f in filas if f["grupo"] == g]
    n = len(sel)
    yc = y0 + fila_h * (i + 0.5)
    t(40, yc - 4, nombre, 15, TINTA, weight="bold")
    t(40, yc + 14, f"{sub} · {n} pairs", 11.5, TENUE)
    for var, col, dy, relleno in (("resonante", AZUL, -9, True), ("placebo_limpio", BERMELLON, 9, False)):
        k = sum(int(f[var]) for f in sel)
        p, lo, hi = wilson(k, n)
        y = yc + dy
        s.append(f'<line x1="{X(lo):.1f}" x2="{X(hi):.1f}" y1="{y:.1f}" y2="{y:.1f}" stroke="{col}" '
                 f'stroke-width="2" stroke-linecap="round"/>')
        s.append(f'<circle cx="{X(p):.1f}" cy="{y:.1f}" r="5.5" fill="{col if relleno else "#ffffff"}" '
                 f'stroke="{col}" stroke-width="2"/>')
        if var == "resonante":
            t(X(hi) + 10, y + 4, f"{k}/{n} = {100 * p:.0f}%", 12, TINTA2)

ly = H - 22
s.append(f'<circle cx="46" cy="{ly - 4}" r="5.5" fill="{AZUL}" stroke="{AZUL}" stroke-width="2"/>')
t(58, ly, "near a first-order commensurability", 12, TINTA2)
s.append(f'<circle cx="326" cy="{ly - 4}" r="5.5" fill="#ffffff" stroke="{BERMELLON}" stroke-width="2"/>')
t(338, ly, "in the gaps between commensurabilities (negative control)", 12, TINTA2)
s.append("</svg>")
out = BASE / "results" / "figure1_youth_resonance.svg"
out.write_text("\n".join(s))
print(out)
