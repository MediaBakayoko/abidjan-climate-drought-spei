#!/usr/bin/env python3
"""
Redesign des visuels SPEI annuels : remplace les pseudo-cartes (10x9
pixels, resolution trop grossiere pour un rendu spatial credible -- le
lissage bilineaire precedent masquait ce manque de detail en produisant
un flou qui a ete juge illisible/peu professionnel) par des cartes-
valeur nettes : grand chiffre SPEI, statut, barre de positionnement sur
l'echelle -- honnete sur le niveau de detail reel disponible (donnee
zonale moyenne, pas une carte infra-urbaine).
"""
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

BASE_DIR = Path(__file__).resolve().parent.parent
MAPS_DIR = BASE_DIR / "maps"
DATA_DIR = BASE_DIR / "data"

FONT_BOLD = "C:/Windows/Fonts/arialbd.ttf"
FONT_REGULAR = "C:/Windows/Fonts/arial.ttf"

W, H = 1200, 900


def spei_color(v):
    if v <= -1.5:
        return (165, 0, 38)
    if v <= -1:
        return (215, 48, 39)
    if v <= -0.5:
        return (244, 109, 67)
    if v < 0.5:
        return (150, 150, 150)
    if v < 1:
        return (116, 173, 209)
    if v < 1.5:
        return (69, 117, 180)
    return (49, 54, 149)


def status_label(v):
    if v <= -1.5:
        return "SECHERESSE SEVERE"
    if v <= -1:
        return "SECHERESSE MODEREE"
    if v < 1:
        return "CONDITIONS NORMALES"
    if v < 1.5:
        return "HUMIDE"
    return "TRES HUMIDE"


def render_value_card(year, spei, precip, etp, bilan, out_path):
    color = spei_color(spei)
    img = Image.new("RGB", (W, H), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    f_year = ImageFont.truetype(FONT_BOLD, 54)
    f_title = ImageFont.truetype(FONT_BOLD, 30)
    f_value = ImageFont.truetype(FONT_BOLD, 160)
    f_status = ImageFont.truetype(FONT_BOLD, 38)
    f_detail = ImageFont.truetype(FONT_REGULAR, 26)
    f_scale = ImageFont.truetype(FONT_REGULAR, 22)

    # Bandeau haut colore
    draw.rectangle([0, 0, W, 170], fill=color)
    draw.text((50, 30), f"Abidjan {year}", fill=(255, 255, 255), font=f_year)
    draw.text((50, 100), "Indice SPEI (secheresse / humidite)", fill=(255, 255, 255), font=f_title)

    # Grande valeur SPEI
    value_text = f"{spei:+.2f}"
    tw = draw.textlength(value_text, font=f_value)
    draw.text(((W - tw) / 2, 220), value_text, fill=color, font=f_value)

    # Statut
    status_text = status_label(spei)
    tw2 = draw.textlength(status_text, font=f_status)
    draw.text(((W - tw2) / 2, 410), status_text, fill=(40, 40, 40), font=f_status)

    # Echelle horizontale avec curseur
    bar_x0, bar_x1, bar_y = 150, W - 150, 500
    bar_w = bar_x1 - bar_x0
    for x in range(bar_w):
        t = -3 + (x / bar_w) * 6
        c = spei_color(t)
        draw.line([(bar_x0 + x, bar_y), (bar_x0 + x, bar_y + 28)], fill=c)
    draw.rectangle([bar_x0, bar_y, bar_x1, bar_y + 28], outline=(60, 60, 60), width=2)
    # curseur
    cursor_x = bar_x0 + (max(min(spei, 3), -3) + 3) / 6 * bar_w
    draw.polygon([(cursor_x - 12, bar_y - 18), (cursor_x + 12, bar_y - 18), (cursor_x, bar_y + 2)], fill=(20, 20, 20))

    for val, label in [(-3, "-3"), (-1, "-1"), (0, "0"), (1, "+1"), (3, "+3")]:
        x_pos = bar_x0 + (val + 3) / 6 * bar_w
        tw3 = draw.textlength(label, font=f_scale)
        draw.text((x_pos - tw3 / 2, bar_y + 34), label, fill=(60, 60, 60), font=f_scale)

    draw.text((bar_x0, bar_y + 70), "Secheresse", fill=(165, 0, 38), font=f_scale)
    tw4 = draw.textlength("Tres humide", font=f_scale)
    draw.text((bar_x1 - tw4, bar_y + 70), "Tres humide", fill=(49, 54, 149), font=f_scale)

    # Details chiffres reels
    y = 650
    draw.line([(100, y), (W - 100, y)], fill=(220, 220, 220), width=2)
    y += 30
    details = [
        f"Precipitations totales (CHIRPS) : {precip:.0f} mm",
        f"Evapotranspiration totale (ERA5-Land) : {etp:.0f} mm",
        f"Bilan hydrique (P - ETP) : {bilan:+.0f} mm",
    ]
    for d in details:
        draw.text((100, y), d, fill=(50, 50, 50), font=f_detail)
        y += 40

    draw.text((100, H - 50), "District d'Abidjan, Cote d'Ivoire -- Google Earth Engine (CHIRPS + ERA5-Land)", fill=(130, 130, 130), font=ImageFont.truetype(FONT_REGULAR, 20))

    img.save(out_path)
    print(f"Carte-valeur ecrite : {out_path}")


def main():
    stats = json.loads((DATA_DIR / "stats.json").read_text(encoding="utf-8"))
    for r in stats["resultats_annuels"]:
        render_value_card(
            r["annee"], r["spei"],
            r["precipitation_totale_mm"], r["evapotranspiration_totale_mm"], r["bilan_hydrique_mm"],
            MAPS_DIR / f"spei_{r['annee']}.png",
        )


if __name__ == "__main__":
    main()
