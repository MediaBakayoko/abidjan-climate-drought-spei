# Analyse climatique du District d'Abidjan 2015-2025 : sécheresse et conditions d'humidité (SPEI)

Inspiré de l'analyse « Climate Change Analysis in Ethiopia Using Remote
Sensing (2015–2025) » — même méthodologie (SPEI, CHIRPS, climatologie de
référence), appliquée au District d'Abidjan, Côte d'Ivoire.

## Objectif

Étudier les variations interannuelles des conditions d'humidité et de
sécheresse à Abidjan sur 2015-2025, à partir du bilan hydrique
(précipitations − évapotranspiration) standardisé par rapport à une
climatologie de référence.

## Zone d'étude
District d'Abidjan, Côte d'Ivoire — rectangle `[-4.20, 5.20, -3.80, 5.55]`
(même emprise que la série `20-Jours-Google-Earth-Engine-Cours`).

## Données sources
| Source | Rôle | Résolution native |
|---|---|---|
| UCSB-CHG/CHIRPS/DAILY | Précipitations | ~5,5 km |
| ECMWF/ERA5_LAND/DAILY_AGGR (`total_evaporation_sum`) | Évapotranspiration (proxy) | ~9 km (rééchantillonné ERA5-Land à ~11 km natif) |

**Période analysée** : année complète (janvier-décembre), 2015-2025 (11 ans).
**Climatologie de référence** : 2010-2019 (10 ans).

## Méthode

1. **Bilan hydrique annuel** D = Précipitation − Évapotranspiration, moyenné
   sur la zone, pour chaque année de 2010 à 2025.
2. **Ajustement d'une distribution logistique généralisée (GLO)** sur les
   10 valeurs de la période de référence (2010-2019), via la méthode des
   **L-moments** (package `lmoments3`, formules de Hosking 1990) — méthode
   standard du SPEI original (Vicente-Serrano et al. 2010).
3. **Standardisation** : chaque valeur annuelle est convertie en SPEI par
   la fonction de répartition cumulée (CDF) de la distribution ajustée,
   transformée en z-score normal.
4. Application de la même transformation à chaque pixel du bilan hydrique
   pour produire une carte spatiale par année.

### Choix méthodologique documenté : L-moments plutôt que maximum de vraisemblance
Un premier essai avec `scipy.stats.fisk.fit()` (maximum de vraisemblance)
a produit un résultat non-physique : **11 années sur 11 ressortaient en
conditions humides**, ce qui est statistiquement incohérent pour un indice
standardisé. Diagnostic : le MLE « épingle » le paramètre de position sur
la valeur minimale exacte de l'échantillon de référence (n=10, petit et
asymétrique) — défaut connu de cette méthode sur petit échantillon. La
méthode des L-moments (utilisée par le package `lmoments3`, conforme à la
pratique standard du SPEI) corrige ce problème et produit des résultats
statistiquement cohérents (années de sécheresse réelles identifiées).

## Résultats réels

| Année | Précipitations (mm) | Évapotranspiration (mm) | Bilan (mm) | SPEI | Catégorie |
|---|---|---|---|---|---|
| 2015 | 1438 | 1359 | +93 | -0,76 | Normal (tendant sec) |
| 2016 | 1334 | 1339 | +14 | **-1,03** | **Sécheresse modérée** |
| 2017 | 1861 | 1374 | +504 | +1,12 | Humide |
| 2018 | 1446 | 1390 | +75 | -0,83 | Normal (tendant sec) |
| 2019 | 1672 | 1384 | +315 | +0,19 | Normal |
| 2020 | 1860 | 1355 | +534 | +1,27 | Humide |
| 2021 | 1731 | 1416 | +342 | +0,32 | Normal |
| 2022 | 1955 | 1344 | +628 | +1,74 | Très humide |
| 2023 | 2318 | 1336 | +1006 | **+3,64** | Très humide (extrême, voir limites) |
| 2024 | 1737 | 1317 | +435 | +0,78 | Normal |
| 2025 | 2088 | 1284 | +824 | +2,69 | Très humide |

**Tendance visible** : après une période 2015-2018 plutôt déficitaire
(dont 2016 en sécheresse modérée confirmée), les années 2020-2025 montrent
une humidité croissante et quasi continue, avec des précipitations totales
dépassant 2000 mm/an à partir de 2022.

## Cartes
- `maps/spei_<année>.png` / `.tif` — carte SPEI par année, échelle -2 à
  +2, palette rouge (sec) → vert (humide).
- `maps/serie_temporelle_spei.png` — graphique en barres de l'évolution
  du SPEI moyen zone, 2015-2025.

## Limites méthodologiques

- **ERA5-Land `total_evaporation_sum` est l'évapotranspiration réelle**
  (limitée par l'eau disponible au sol), pas l'évapotranspiration
  potentielle (PET) utilisée dans la formulation originale du SPEI —
  utilisée ici comme proxy faute de produit PET direct fiable sur Google
  Earth Engine pour cette zone.
- **Échantillon de référence court (n=10, 2010-2019)** : suffisant pour
  détecter la tendance relative entre années, mais insuffisant pour
  calibrer fiablement les valeurs extrêmes. La valeur SPEI de 2023
  (+3,64) dépasse l'échelle usuelle du SPEI (généralement bornée à ±3
  dans la littérature) — signe que l'extrapolation de la queue de
  distribution au-delà de l'échantillon de calibration n'est pas fiable
  en valeur absolue à ce niveau. **La tendance (2023 et 2025 nettement
  plus humides que la référence) est interprétable ; le chiffre exact
  +3,64 ne doit pas être cité comme une valeur SPEI calibrée précise.**
- **Résolution spatiale grossière** : CHIRPS (~5,5 km) et ERA5-Land
  (~9-11 km) produisent une grille de seulement 10×9 pixels sur cette
  zone. Les cartes sont lissées visuellement (interpolation bilinéaire)
  pour la lisibilité, mais ne montrent aucune variation fine
  infra-urbaine — à l'échelle de l'agglomération seulement.
- **Standardisation spatiale non pixel-par-pixel** : la même fonction de
  standardisation (ajustée sur la moyenne zonale) est appliquée à chaque
  pixel, faute de série de référence suffisamment longue par pixel sur
  cette durée.

## Script
Voir `scripts/run_spei_abidjan.py` (exécution réelle via l'API Python
Earth Engine + `lmoments3` + `scipy`).

## Article scientifique
Voir [`ARTICLE.md`](./ARTICLE.md).
