#!/usr/bin/env python3
"""
Analyse climatique Abidjan 2015-2025 : secheresse et conditions
d'humidite par teledetection (CHIRPS + ERA5-Land), suivant la
methodologie SPEI (Standardized Precipitation Evapotranspiration Index,
Vicente-Serrano et al. 2010).

Zone d'etude : District d'Abidjan, Cote d'Ivoire -- meme emprise que la
serie 20-Jours-GEE-Cours.

Periode : annee complete (janvier-decembre), 2015-2025.
Climatologie de reference : 2010-2019 (meme fenetre que le modele
reproduit, adaptee a la disponibilite reelle des donnees).

Methode :
1. Calcul du bilan hydrique D = Precipitation - Evapotranspiration pour
   chaque annee (CHIRPS pour P, ERA5-Land pour ETP reelle -- faute de
   produit PET direct fiable sur GEE, utilise comme proxy, limite
   explicite documentee).
2. Ajustement d'une distribution log-logistique (fisk) sur la serie D
   de la periode de reference.
3. Standardisation de chaque annee par la fonction de repartition
   ajustee -> SPEI (z-score de la probabilite cumulee).
4. Cartes spatiales du SPEI par annee + serie temporelle moyenne sur
   la zone.

Produit :
- maps/spei_<annee>.tif / .png (une carte par annee, legende integree)
- maps/serie_temporelle_spei.png (graphique de la moyenne SPEI zone par annee)
- data/stats.json (statistiques reelles completes)
"""
import json
from pathlib import Path

import ee
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import requests
from PIL import Image, ImageDraw, ImageFont
from scipy import stats
import lmoments3.distr as lmom_distr

CREDENTIALS_PATH = "C:/Users/ASUS/Downloads/ornate-serenity-509219-d5-08760d05d66d.json"
BASE_DIR = Path(__file__).resolve().parent.parent
MAPS_DIR = BASE_DIR / "maps"
DATA_DIR = BASE_DIR / "data"
MAPS_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

STUDY_AREA_COORDS = [[-4.20, 5.20], [-3.80, 5.55]]
ANALYSIS_YEARS = list(range(2015, 2026))
BASELINE_YEARS = list(range(2010, 2020))  # 2010-2019, 10 ans

FONT_BOLD = "C:/Windows/Fonts/arialbd.ttf"
FONT_REGULAR = "C:/Windows/Fonts/arial.ttf"


def init_ee():
    with open(CREDENTIALS_PATH) as f:
        info = json.load(f)
    credentials = ee.ServiceAccountCredentials(info["client_email"], CREDENTIALS_PATH)
    ee.Initialize(credentials)
    print("GEE initialise avec le compte de service :", info["client_email"])


def build_study_area():
    return ee.Geometry.Rectangle(coords=STUDY_AREA_COORDS, geodesic=False)


def annual_water_balance(study_area, year):
    """Bilan hydrique annuel D = P - ETP pour une annee donnee.
    P : somme annuelle CHIRPS (mm). ETP : somme annuelle
    evapotranspiration reelle ERA5-Land (mm, converti depuis m).
    """
    start = f"{year}-01-01"
    end = f"{year}-12-31"

    chirps = (
        ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
        .filterBounds(study_area)
        .filterDate(start, end)
        .select("precipitation")
        .sum()
    )

    era5 = (
        ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
        .filterBounds(study_area)
        .filterDate(start, end)
        .select("total_evaporation_sum")
        .sum()
        .multiply(-1000)  # m (negatif, convention ERA5) -> mm positif
    )

    balance = chirps.subtract(era5).rename("D").clip(study_area)
    return balance, chirps.clip(study_area), era5.clip(study_area)


def reduce_mean(image, study_area, scale=5000):
    val = image.reduceRegion(reducer=ee.Reducer.mean(), geometry=study_area, scale=scale, maxPixels=1e9).getInfo()
    return val.get(list(val.keys())[0]) if val else None


def download_geotiff(image, band_name, filename, region, scale=5000):
    url = image.select(band_name).getDownloadURL({"region": region, "scale": scale, "format": "GEO_TIFF"})
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    path = MAPS_DIR / filename
    path.write_bytes(resp.content)
    return path


def fit_spei_distribution(baseline_values):
    """Ajuste une distribution logistique generalisee (GLO) sur les
    valeurs de bilan hydrique de la periode de reference, via la
    methode des L-moments -- methode STANDARD du SPEI original
    (Vicente-Serrano et al. 2010). Utilise le package 'lmoments3'
    (implementation testee des formules de Hosking 1990), pas une
    formule reimplementee a la main.

    Deux approches rejetees avant celle-ci, pour mémoire :
    - scipy.stats.fisk.fit() (MLE) : epingle le parametre de position
      sur le minimum exact de l'echantillon (n=10, petit et
      asymetrique) -- rend presque toutes les annees artificiellement
      "humides" (confirme : 11/11 annees positives, physiquement
      incoherent pour un indice standardise).
    - Formule PWM reimplementee a la main : parametres non-physiques
      (alpha, beta negatifs) -- probable erreur de signe dans la
      formule retranscrite.
    """
    params = lmom_distr.glo.lmom_fit(baseline_values)
    return params


def water_balance_to_spei(value, glo_params):
    """Convertit une valeur de bilan hydrique en indice SPEI (z-score
    de la probabilite cumulee selon la distribution logistique
    generalisee ajustee par L-moments)."""
    dist = lmom_distr.glo(**glo_params)
    cdf = float(dist.cdf(value))
    cdf = np.clip(cdf, 1e-6, 1 - 1e-6)
    return float(stats.norm.ppf(cdf))


def scaled_fonts(img_width):
    unit = max(img_width / 70, 22)
    return (
        ImageFont.truetype(FONT_BOLD, int(unit * 1.5)),
        ImageFont.truetype(FONT_REGULAR, int(unit * 1.05)),
        ImageFont.truetype(FONT_BOLD, int(unit * 1.15)),
    )


SPEI_COLOR_STOPS = [
    (-2.0, (103, 0, 13)), (-1.0, (239, 59, 44)), (0.0, (254, 240, 138)),
    (1.0, (65, 171, 93)), (2.0, (0, 68, 27)),
]


def interpolate_color(value, stops):
    positions = [s[0] for s in stops]
    r = np.interp(value, positions, [s[1][0] for s in stops])
    g = np.interp(value, positions, [s[1][1] for s in stops])
    b = np.interp(value, positions, [s[1][2] for s in stops])
    return r, g, b


def render_spei_png(tif_path, out_path, year, mean_spei):
    with rasterio.open(tif_path) as src:
        arr = src.read(1).astype(float)
        nodata_mask = src.read_masks(1) == 0

    norm_vals = np.clip(arr, -2, 2)
    r, g, b = interpolate_color(norm_vals, SPEI_COLOR_STOPS)
    rgb = np.stack([r, g, b], axis=-1).astype(np.uint8)
    rgb[nodata_mask] = (255, 255, 255)

    img = Image.fromarray(rgb, mode="RGB")
    # CHIRPS a une resolution native d'environ 5,5 km -- un
    # suréchantillonnage brut (NEAREST) sur une petite zone produit un
    # effet de mosaique en gros carres, pas lisible. Lissage bilineaire
    # pour un rendu visuel presentable -- n'ajoute aucune information
    # reelle, documente explicitement dans le README/article.
    target_width = 1000
    if img.width < target_width:
        scale = target_width / img.width
        img = img.resize((int(img.width * scale), int(img.height * scale)), Image.BILINEAR)

    font_title, font_subtitle, font_label = scaled_fonts(img.width)
    margin = int(img.width * 0.03)
    bar_h = int(img.width * 0.03)
    legend_h = int(img.width * 0.22)

    canvas = Image.new("RGB", (img.width, img.height + legend_h), (255, 255, 255))
    canvas.paste(img, (0, legend_h))
    draw = ImageDraw.Draw(canvas)

    y = margin
    draw.text((margin, y), f"SPEI {year} - District d'Abidjan", fill=(20, 20, 20), font=font_title)
    y += int(font_title.size * 1.3)
    draw.text((margin, y), "Bilan hydrique CHIRPS - ERA5-Land, standardise (methode log-logistique)", fill=(90, 90, 90), font=font_subtitle)
    y += int(font_subtitle.size * 1.4)
    status = "Secheresse" if mean_spei < -1 else ("Humidite excedentaire" if mean_spei > 1 else "Conditions normales")
    draw.text((margin, y), f"SPEI moyen zone : {mean_spei:+.2f} ({status})", fill=(20, 20, 20), font=font_label)
    y += int(font_label.size * 1.6)

    bar_w = img.width - 2 * margin
    bar_y = y
    for x in range(bar_w):
        t = -2 + (x / bar_w) * 4
        rr, gg, bb = interpolate_color(t, SPEI_COLOR_STOPS)
        draw.line([(margin + x, bar_y), (margin + x, bar_y + bar_h)], fill=(int(rr), int(gg), int(bb)))
    draw.rectangle([margin, bar_y, margin + bar_w, bar_y + bar_h], outline=(60, 60, 60), width=2)
    label_y = bar_y + bar_h + int(font_label.size * 0.3)
    for val, text in [(-2, "-2 (secheresse extreme)"), (0, "0 (normal)"), (2, "+2 (tres humide)")]:
        x_pos = margin + (val + 2) / 4 * bar_w
        tw = draw.textlength(text, font=font_label)
        draw.text((min(max(x_pos - tw / 2, margin), margin + bar_w - tw), label_y), text, fill=(20, 20, 20), font=font_label)

    canvas.save(out_path)
    print(f"PNG ecrit : {out_path}")


def main():
    init_ee()
    study_area = build_study_area()
    region = study_area.bounds()

    print("=== Calcul de la climatologie de reference (2010-2019) ===")
    baseline_balances = []
    for year in BASELINE_YEARS:
        balance_img, _, _ = annual_water_balance(study_area, year)
        mean_val = reduce_mean(balance_img, study_area)
        baseline_balances.append(mean_val)
        print(f"  {year}: bilan hydrique moyen zone = {mean_val:.1f} mm")

    fisk_params = fit_spei_distribution(baseline_balances)
    print(f"Parametres distribution logistique generalisee (GLO) ajustes par L-moments : {dict(fisk_params)}")

    print("\n=== Calcul SPEI annuel 2015-2025 ===")
    results = []
    for year in ANALYSIS_YEARS:
        balance_img, p_img, etp_img = annual_water_balance(study_area, year)
        mean_balance = reduce_mean(balance_img, study_area)
        mean_p = reduce_mean(p_img, study_area)
        mean_etp = reduce_mean(etp_img, study_area)
        mean_spei = float(water_balance_to_spei(mean_balance, fisk_params))

        # Carte spatiale SPEI : la distribution GLO ajustee sur la
        # moyenne zonale (seule serie temporelle assez longue dont on
        # dispose, n=10) est appliquee a chaque pixel du bilan
        # hydrique brut -- limite explicite : ce n'est pas un ajustement
        # pixel par pixel (qui demanderait une serie de reference par
        # pixel, hors de portee avec 10 ans de donnees), mais une
        # application spatiale de la meme fonction de standardisation.
        balance_tif = download_geotiff(balance_img, "D", f"_balance_{year}.tif", region, scale=5000)
        with rasterio.open(balance_tif) as src:
            balance_arr = src.read(1).astype(float)
            profile = src.profile
        dist = lmom_distr.glo(**fisk_params)
        cdf_arr = dist.cdf(balance_arr)
        cdf_arr = np.clip(cdf_arr, 1e-6, 1 - 1e-6)
        spei_arr_pixels = stats.norm.ppf(cdf_arr)

        spei_tif_path = MAPS_DIR / f"spei_{year}.tif"
        profile.update(dtype="float32", count=1)
        with rasterio.open(spei_tif_path, "w", **profile) as dst:
            dst.write(spei_arr_pixels.astype(np.float32), 1)
        balance_tif.unlink()  # fichier intermediaire, pas un livrable
        render_spei_png(spei_tif_path, MAPS_DIR / f"spei_{year}.png", year, mean_spei)

        results.append({
            "annee": year,
            "precipitation_totale_mm": round(mean_p, 1),
            "evapotranspiration_totale_mm": round(mean_etp, 1),
            "bilan_hydrique_mm": round(mean_balance, 1),
            "spei": round(mean_spei, 2),
            "categorie": (
                "secheresse_severe" if mean_spei < -1.5 else
                "secheresse_moderee" if mean_spei < -1 else
                "normal" if mean_spei < 1 else
                "humide" if mean_spei < 1.5 else
                "tres_humide"
            ),
        })
        print(f"  {year}: P={mean_p:.0f}mm ETP={mean_etp:.0f}mm bilan={mean_balance:+.0f}mm SPEI={mean_spei:+.2f}")

    # Serie temporelle
    years_arr = [r["annee"] for r in results]
    spei_arr = [r["spei"] for r in results]
    fig, ax = plt.subplots(figsize=(11, 5))
    colors = ["#a50026" if v < -1 else "#4575b4" if v > 1 else "#999999" for v in spei_arr]
    ax.bar(years_arr, spei_arr, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(-1, color="#a50026", linestyle="--", linewidth=0.8, label="Seuil secheresse moderee (SPEI < -1)")
    ax.axhline(1, color="#4575b4", linestyle="--", linewidth=0.8, label="Seuil humide (SPEI > 1)")
    ax.set_ylabel("SPEI (moyenne zone, annuel)")
    ax.set_title("SPEI annuel - District d'Abidjan (2015-2025)\nClimatologie de reference 2010-2019 (CHIRPS + ERA5-Land)")
    ax.legend(fontsize=9)
    ax.set_xticks(years_arr)
    fig.tight_layout()
    serie_path = MAPS_DIR / "serie_temporelle_spei.png"
    fig.savefig(serie_path, dpi=150)
    plt.close(fig)
    print(f"Serie temporelle ecrite : {serie_path}")

    stats_out = {
        "zone_etude": "District d'Abidjan, Cote d'Ivoire",
        "bbox": STUDY_AREA_COORDS,
        "periode_analyse": "2015-2025 (annee complete Jan-Dec)",
        "periode_reference_climatologie": "2010-2019",
        "sources": {
            "precipitation": "UCSB-CHG/CHIRPS/DAILY",
            "evapotranspiration": "ECMWF/ERA5_LAND/DAILY_AGGR (total_evaporation_sum, proxy ETP reelle)",
        },
        "methode_spei": "Distribution logistique generalisee (GLO) ajustee par L-moments (package lmoments3, formules de Hosking 1990) sur la serie de reference, standardisation par CDF -> z-score normal (Vicente-Serrano et al. 2010)",
        "parametres_distribution_reference": dict(fisk_params),
        "bilan_hydrique_reference_2010_2019_mm": baseline_balances,
        "resultats_annuels": results,
        "avertissement_methodologique": (
            "ERA5-Land total_evaporation_sum est l'evapotranspiration reelle (limitee par l'eau "
            "disponible), pas l'evapotranspiration potentielle (PET) utilisee dans le SPEI "
            "original -- utilisee ici comme proxy faute de produit PET direct fiable sur GEE pour "
            "cette zone. La distribution de reference (GLO, L-moments) est ajustee sur seulement "
            "10 annees (2010-2019) -- echantillon court pour calibrer fiablement les valeurs "
            "extremes (SPEI tres eleve/bas peu fiable en valeur absolue, la tendance relative "
            "entre annees reste interpretable). La carte spatiale applique la meme fonction de "
            "standardisation (ajustee sur la moyenne zonale) a chaque pixel -- ce n'est pas un "
            "ajustement pixel par pixel independant, limite assumee faute de serie de reference "
            "suffisamment longue par pixel."
        ),
    }
    stats_path = DATA_DIR / "stats.json"
    stats_path.write_text(json.dumps(stats_out, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(f"\nStatistiques ecrites : {stats_path}")


if __name__ == "__main__":
    main()
