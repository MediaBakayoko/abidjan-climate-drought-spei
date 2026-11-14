# Variabilité interannuelle de l'humidité et de la sécheresse dans le District d'Abidjan (2015-2025) : approche par l'indice SPEI

**Auteur :** Média Marcel Bakayoko
**Zone d'étude :** District d'Abidjan, Côte d'Ivoire
**Données :** Google Earth Engine — CHIRPS Daily, ERA5-Land Daily Aggregated
**Date d'analyse :** Novembre 2026

## Résumé

Cette étude quantifie la variabilité interannuelle des conditions
d'humidité et de sécheresse dans le District d'Abidjan sur la période
2015-2025, à l'aide de l'indice SPEI (Standardized Precipitation
Evapotranspiration Index), calculé à partir du bilan hydrique
(précipitations CHIRPS moins évapotranspiration ERA5-Land) standardisé
par rapport à une climatologie de référence (2010-2019) via une
distribution logistique généralisée ajustée par la méthode des
L-moments. L'analyse révèle une période de déficit modéré en 2015-2018
(SPEI minimal de -1,03 en 2016), suivie d'une tendance à l'humidité
croissante et quasi continue de 2019 à 2025, culminant en 2023 (bilan
hydrique de +1006 mm, SPEI +3,64) et 2025 (+824 mm, SPEI +2,69).

## 1. Introduction

Le suivi des conditions climatiques par télédétection permet de
documenter, de façon reproductible et peu coûteuse, les variations
d'humidité et de sécheresse dans des zones où les stations
météorologiques au sol sont rares ou discontinues. L'indice SPEI
(Vicente-Serrano et al., 2010) est largement utilisé en climatologie
opérationnelle car il combine précipitations et demande évaporative,
contrairement à l'indice SPI plus ancien qui ne considère que les
précipitations — un atout particulièrement pertinent sous climat
équatorial où la demande évaporative reste élevée toute l'année.

## 2. Zone d'étude et données

La zone d'étude couvre le rectangle `[-4.20°, 5.20°, -3.80°, 5.55°]`,
englobant le District d'Abidjan et sa lagune Ébrié, cohérente avec les
analyses précédentes de la même agglomération (NDVI, NDBI, classification
d'occupation du sol).

| Source | Variable | Résolution native | Période utilisée |
|---|---|---|---|
| UCSB-CHG/CHIRPS/DAILY | Précipitation | ~5,5 km | 2010-2025 (sommes annuelles) |
| ECMWF/ERA5_LAND/DAILY_AGGR | Évapotranspiration réelle (`total_evaporation_sum`) | ~9-11 km | 2010-2025 (sommes annuelles) |

La climatologie de référence couvre 2010-2019 (10 années), fenêtre
dictée par la disponibilité conjointe des deux produits sur Google Earth
Engine et la volonté de conserver une décennie complète avant la période
d'analyse (2015-2025).

## 3. Méthode

### 3.1 Bilan hydrique annuel

Pour chaque année *i*, le bilan hydrique moyen sur la zone est calculé :

```
D_i = Σ(précipitation_journalière) - Σ(évapotranspiration_journalière)
```

### 3.2 Ajustement de la distribution de référence

Le SPEI standard ajuste une distribution log-logistique à 3 paramètres
sur la série du bilan hydrique de la période de référence. Cette étude
utilise la **distribution logistique généralisée (GLO)**, mathématiquement
équivalente, ajustée par la **méthode des L-moments** (moments pondérés
par probabilité, formules de Hosking 1990), implémentée par le package
Python `lmoments3` — approche désormais standard pour ce type d'ajustement
sur petit échantillon, préférée au maximum de vraisemblance.

**Note méthodologique sur un écueil évité** : un premier ajustement par
maximum de vraisemblance (`scipy.stats.fisk.fit()`) a produit un résultat
non-physique — la totalité des 11 années analysées ressortaient en
conditions humides, résultat statistiquement incohérent pour un indice
censé osciller autour de zéro par construction. L'inspection des
paramètres a révélé que le MLE avait fixé le paramètre de position
exactement sur la valeur minimale de l'échantillon de référence (n=10),
un artefact documenté du maximum de vraisemblance sur petits échantillons
asymétriques. Le passage aux L-moments a résolu ce problème et produit des
résultats cohérents avec la variabilité interannuelle réellement attendue.

### 3.3 Standardisation

Pour chaque année, la valeur du bilan hydrique est convertie en
probabilité cumulée (CDF) selon la distribution ajustée, puis transformée
en z-score de la loi normale standard — c'est la valeur SPEI.

## 4. Résultats

| Année | Précipitation (mm) | ETP (mm) | Bilan (mm) | SPEI |
|---|---|---|---|---|
| 2015 | 1438 | 1359 | +93 | -0,76 |
| 2016 | 1334 | 1339 | +14 | -1,03 |
| 2017 | 1861 | 1374 | +504 | +1,12 |
| 2018 | 1446 | 1390 | +75 | -0,83 |
| 2019 | 1672 | 1384 | +315 | +0,19 |
| 2020 | 1860 | 1355 | +534 | +1,27 |
| 2021 | 1731 | 1416 | +342 | +0,32 |
| 2022 | 1955 | 1344 | +628 | +1,74 |
| 2023 | 2318 | 1336 | +1006 | +3,64 |
| 2024 | 1737 | 1317 | +435 | +0,78 |
| 2025 | 2088 | 1284 | +824 | +2,69 |

La série temporelle (figure 1) révèle deux phases distinctes : une
période 2015-2018 globalement déficitaire ou proche de la normale, avec
un épisode de sécheresse modérée identifié en 2016 (SPEI -1,03), suivie
d'une transition marquée vers des conditions durablement humides à partir
de 2019-2020, qui s'intensifie nettement en fin de série (2022-2025).
L'évapotranspiration totale annuelle diminue légèrement sur la période
récente (de ~1420 mm en 2021 à ~1284 mm en 2025), pendant que les
précipitations augmentent fortement, accentuant l'écart du bilan
hydrique.

Les cartes spatiales par année (figures 2-12) confirment cette tendance
à l'échelle de l'agglomération, avec une relative homogénéité spatiale du
signal SPEI sur la zone — attendu compte tenu de la résolution native
grossière des produits sources (~5,5 à 11 km) relativement à la taille du
District d'Abidjan.

## 5. Discussion et limites méthodologiques

**Proxy d'évapotranspiration.** ERA5-Land `total_evaporation_sum` mesure
l'évapotranspiration réelle (contrainte par l'humidité disponible au
sol), et non l'évapotranspiration potentielle (PET, demande atmosphérique
non contrainte) utilisée dans la définition originale du SPEI. En
contexte équatorial humide comme Abidjan, où l'eau disponible est rarement
le facteur limitant, l'écart entre ETP réelle et PET est probablement
modéré, mais reste une approximation non validée pour cette zone
spécifique.

**Fiabilité des valeurs extrêmes.** L'échantillon de référence (10
années) permet une estimation robuste de la tendance centrale et des
écarts modérés, mais la valeur SPEI de 2023 (+3,64) dépasse largement
l'échelle usuelle de cet indice (généralement comprise entre -3 et +3
dans la littérature). Cela signale une extrapolation de la queue de
distribution au-delà du domaine couvert par l'échantillon de calibration
— la **tendance** (2023 nettement plus humide que la référence) est
robuste, mais le **chiffre précis** ne doit pas être interprété comme une
mesure SPEI calibrée fiable à ce niveau d'extrême.

**Résolution spatiale.** La grille native des deux produits sources
(~10×9 pixels sur cette zone) ne permet aucune analyse infra-urbaine ;
les cartes produites sont lissées visuellement pour la lisibilité sans
ajouter d'information réelle.

**Pas de validation terrain.** Aucune mesure pluviométrique au sol n'a
été utilisée pour valider les sommes CHIRPS sur cette zone spécifique.

## 6. Conclusion

Cette analyse documente, de façon reproductible, une transition climatique
nette dans le District d'Abidjan entre une période de déficit hydrique
modéré (2015-2018) et une période d'humidité croissante et marquée
(2020-2025), avec un pic d'humidité exceptionnel en 2023. Si la tendance
générale est robuste, la magnitude exacte des valeurs extrêmes
nécessiterait un allongement de la période de référence (données
antérieures à 2010, si disponibles sur un produit équivalent) pour être
calibrée avec davantage de confiance. Un prolongement naturel de cette
étude consisterait à croiser ces résultats avec des données agricoles
(rendements, calendriers culturaux) pour évaluer l'impact concret de
cette variabilité sur la sécurité alimentaire régionale — usage classique
du SPEI en appui à la planification agricole et à l'alerte précoce.

## Références

- Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A
  multiscalar drought index sensitive to global warming: the standardized
  precipitation evapotranspiration index. *Journal of Climate*, 23(7),
  1696-1718.
- Hosking, J. R. M. (1990). L-moments: analysis and estimation of
  distributions using linear combinations of order statistics. *Journal
  of the Royal Statistical Society: Series B*, 52(1), 105-124.
- Funk, C., et al. (2015). The climate hazards infrared precipitation with
  stations—a new environmental record for monitoring extremes. *Scientific
  Data*, 2, 150066.
