> 🌐 **Euskara** · [Castellano](README.es.md) · [English](README.en.md)

# Aralar GrAL — Larre-diagnostikoa teledetekzioz

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
[![Egoera](https://img.shields.io/badge/egoera-abian-yellow.svg)]()

Pipeline integratua Aralarko Natur Parkeko (ES2120011 KBE) larre-diagnostikorako, Sentinel-2 teledetekzioa, GAM bidezko klima-modelizazioa eta Larre Presioaren Indizearen (LPI) laukoitz triangulazioa erabiliz.

**Gradu Amaierako Lana**, Geografia eta Lurralde Antolakuntzako Gradua, UPV/EHU. LIFE Oreka Mendian proiektuaren baitan garatua (LIFE15 NAT/ES/000762).

> **Proiektuaren weba**: https://moxkix.github.io/ARALAR_TFG/
>
> **Dokumentazioa**: `docs/` karpeta automatikoki argitaratzen da weben. Hara igotzen den edozein fitxategi (Word, PDF, CSV, etab.) Dokumentuak atalean txartel gisa agertzen da HTMLa ukitu gabe.

---

## Lanaren laburpena

GrAL honek Aralarko mendi-azpiko larreetako larre-presioa diagnostikatzen du, PROSAIL-NNrekin alderantzikatutako Sentinel-2ren LAI serieak erabiliz, GAM bidez klimatikoki modelizatuak eta **lau estrategia independentek** baliozkotuak 2018 eta 2025 artean. Ikuspegirik originalena **D estrategia** da (pagadi-erreferentzia), ondoko baso-estalki klimazikoa klima-termometro gisa erabiltzen duena larre-seinalea klimatikotik bereizteko.

### Emaitza nagusiak

- **Ia bi urteko eredua**: larre-krisia 2020, 2022 eta 2025ean, 2023ko lasaitzearekin txandakatua (4/4 bateragarritasuna lau estrategien artean).
- **Bordetarako distantziarekin bat datorren gradiente espaziala**: larre hurbila > ertaina > urruna.
- **G6 klima-eredua** LOYO bidez (leave-one-year-out) baliozkotua 2018-2025 seriearen gainean.
- **Diagnostiko operatiboa** LIFE Oreka Mendianen jarraipen moldakorrerako.

---

## Biltegiaren egitura

```
ARALAR_TFG/
├── index.html               GitHub Pages orria (gunearen erroa, hiruleduna EUS/CAS/ING)
├── style.css                Webaren estiloak
├── README.md                Irakurri-nazazu nagusia (euskara, fitxategi hau)
├── README.es.md             Irakurri-nazazu gaztelaniaz
├── README.en.md             Irakurri-nazazu ingelesez
├── LICENSE                  CC BY-NC-ND 4.0
│
├── docs/                    Diagnostikoaren dokumentazioa eta CSVak (Word, PDF, CSV…)
│                            Weben automatikoki zerrendatua GitHubeko APIaren bidez.
│
├── data/                    Pipeline-aren irteera taularrak
│   ├── diagnostico_pastoral.csv
│   ├── validacion_loyo.csv
│   ├── comparacion_modelos.csv
│   ├── hayedo_serie_referencia.csv
│   └── coefficients/        Esportatutako PROSAIL LUTak eta NN koefizienteak
│                            ({PASTOS,HAYEDO}_gee.csv eta *_nn_coefficients.{js,json})
│                            Web-app-ak sortuak eta GEE scriptak kontsumituak.
│
├── figures/                 Lanaren sintesi-irudiak
│   ├── IPP_triangulacion_ABCD.png
│   ├── IPP_barras_interanual.png
│   ├── IPP_por_habitat.png
│   ├── evolucion_IPP_interanual.png
│   ├── gradiente_IPP_distancia_bordas.png
│   ├── diagnostico_integrado_global.png
│   ├── N1_series_LAI.png
│   ├── diag_A_pasto_{cercano,intermedio,remoto}.png
│   ├── diag_A_hayedo_{cercano,intermedio,remoto}.png
│   └── diag_C_anomalias.png
│
├── maps/                    Urteko larre-kartografia
│   ├── 2018/  …  2025/      Serie nagusia (PNG + GeoTIFF + resumen .txt urteko)
│   └── 2018_LIFE/  …  2025_LIFE/
│                            LIFE Oreka Mendianentzat etiketatutako aldaera
│
└── scripts/
    ├── python/
    │   ├── web-app/         FastAPI web-aplikazioa PROSAIL LUTak sortzeko
    │   │                    eta Sentinel-2 → aldagai biofisiko inbertsioa egiteko
    │   │   ├── main.py             FastAPI backend-a (REST + WebSocket)
    │   │   ├── prosail_pure.py     PROSPECT-D + 4SAIL eredua Python hutsean
    │   │   ├── requirements.txt    Mendekotasunak
    │   │   ├── README.md           Instalazio- eta erabilera-gida
    │   │   └── static/index.html   Frontend SPA
    │   ├── train_gee_coefficients.py   NN entrenamendua GEEn injektatzeko
    │   └── cartografia_pastoral_v5.py  Urteko mapen sorrera
    ├── javascripts/         Google Earth Engine pipeline-a
    │   └── GEE_Pipeline_Integrado_Aralar.js
    ├── R/                   Klima-modelizazioa eta laukoitz diagnostikoa
    │   └── diagnostico_pastoral_aralar_v9_4.R
    └── batch/               Windows abiarazleak
        └── ejecutar_cartografia5.bat
```

### PROSAIL web-app (Sentinel-2 inbertsioa)

`scripts/python/web-app/`-ko aplikazioa FastAPI zerbitzari bat da, frontend-a barne. LUTaren parametroak modu interaktiboan konfiguratzeko aukera ematen du (PROSAIL tarteak, geometria solarra, lagintze-banaketa, zarata), memorian sortzeko benetako PROSPECT-D + 4SAIL ereduarekin eta inbertsioa Sentinel-2 irudien gainean exekutatzeko (SNAPeko BEAM-DIMAP edo 10 bandako GeoTIFF) gutxieneko distantziaz cKDTree erabiliz.

```bash
cd scripts/python/web-app
python -m venv venv && source venv/bin/activate     # Windows: venv\Scripts\activate
pip install -r requirements.txt
python main.py
# Ireki http://localhost:8000
```

Azken produktuak (esportatutako LUTak eta GEEn injektatzeko prest dauden sare neuronalaren koefizienteak) `data/coefficients/`-en bertsionatzen dira. Runtime-aren katxea (`temp_prosail/`) biltegitik kanpo dago.

### Urteko produktu kartografikoak

`maps/<urtea>/` karpeta bakoitzak honako mapa hauek ditu PNGn (`aralar_mapa_*_<urtea>.png` / `aralar_pastos_mapa_*_<urtea>.png` izendegia), aldagai bakoitzaren GeoTIFFekin (`.tif`) batera:

| Familia | Mapak |
|---|---|
| LPI triangulazioa | `IPP_A_pastos`, `IPP_A_todos`, `IPP_B_pastos`, `IPP_C_pastos` |
| Sintesia | `IPP_consenso_pastos`, `IPP_concordancia_pastos`, `IPP_espacializado` |
| Kudeaketa | `mapa_biomasa`, `mapa_calidad`, `mapa_carga`, `mapa_balance`, `mapa_riesgo` |
| Kontrola | `mapa_pendiente` |

**D estrategia** diagnostikoan eta triangulazio-irudian integratzen da, baina ez da mapa espazial independente gisa esportatzen. Gainera, `aralar_resumen_<urtea>.txt` bat sortzen da urtearen metadatuekin.

---

## Nola erreproduzitu

### Eskakizunak

- **Python 3.10+** honekin: `prosail`, `scikit-learn`, `numpy`, `rasterio`, `matplotlib`, `pandas`
- **R 4.3+** honekin: `mgcv`, `SPEI`, `ggplot2`, `dplyr`, `tidyr`, `lubridate`
- **Google Earth Engine kontu** aktiboa

### Pipeline-aren urratsak

**1. LUTak eta sare neuronalen entrenamendua (Python lokala)**

```bash
# LUTaren sorrera web-app-arekin (web interfazea localhost:8000-en)
cd scripts/python/web-app
pip install -r requirements.txt
python main.py

# NN entrenamendua esportatutako LUTen gainean eta GEErako koefizienteen esportazioa
cd ../
python train_gee_coefficients.py
```

Esportatutako LUTak eta azken koefizienteak `data/coefficients/`-en geratzen dira.

**2. Eskala handiko inbertsioa GEEn (JavaScript)**

Igo Aralarko assetak GEE proiektura (muga, HIC, bordak, MDT5).
Itsatsi `scripts/javascripts/GEE_Pipeline_Integrado_Aralar.js` kode-editorean eta exekutatu. Driverako esportazioek LAI zonalaren CSVak, klima-CSVak eta urteko GeoTIFFak biltzen dituzte.

> ⚠️ **Dokumentatutako akats kritikoa**: murrizketa zonalean, `ee.Dictionary.combine()` deiek bigarren argumentu gisa `false` eraman behar dute (overwrite=false). Aurreneurri hori gabe, defektuzko balioek benetako balioak gainidazten dituzte eta esportatutako CSVa osorik -9999 ateratzen da.

**3. Klima-modelizazioa eta LPIaren lau estrategiak (R)**

```r
source("scripts/R/diagnostico_pastoral_aralar_v9_4.R")
```

Irteera: `data/diagnostico_pastoral.csv`, `data/validacion_loyo.csv`, `data/comparacion_modelos.csv` eta diagnostiko-irudiak `figures/`-en (grafikoak eta txostena euskaraz).

**4. Larre-kartografia (Python)**

```bash
cd scripts/python
python cartografia_pastoral_v5.py --year 2022
```

Edo exekutatu `scripts/batch/ejecutar_cartografia5.bat` 2018-2025 serie osoa prozesatzeko. Irteera `maps/<urtea>/`-n (mapak izenburu eta legendekin euskaraz).

---

## Proiektuaren weba

Weba (`index.html` + `style.css`) GitHub Pagesekin argitaratzen da `main` adarretik eta biltegiaren errotik. **Hiruleduna da (euskara / gaztelania / ingelesa)**, goiburuko hautagailu batekin; aukera nabigatzailean gogoratzen da. Dokumentuak eta Grafikoak atalen zerrenda orria kargatzean eraikitzen da GitHubeko APIa kontsultatuz (`/repos/Moxkix/ARALAR_TFG/contents/{docs,figures}`). Horrek esan nahi du **karpeta horietara fitxategi bat igotzea nahikoa dela ager dadin**, HTMLa ukitu gabe.

---

## Aipamena

Kode hau edo emaitzak lan eratorrietan erabiltzen badituzu, mesedez aipatu:

```
[Abizena, Izena] (2026). Aralarko Natur Parkeko larre-diagnostikoa
teledetekzioz eta klima-modelizazioz. Gradu Amaierako Lana, Geografia
eta Lurralde Antolakuntzako Gradua, Euskal Herriko Unibertsitatea
/ Universidad del País Vasco.
https://github.com/Moxkix/ARALAR_TFG
```

---

## Lizentzia

Lan hau [Creative Commons Aitortu-EzKomertziala-LanEratorririkGabe 4.0 Internazionala (CC BY-NC-ND 4.0)](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.eu) lizentziapean argitaratzen da.

Egin dezakezu:
- Materiala edozein euskarri edo formatutan partekatu
- Lan akademikoetan aipatu, beti aitortza zuzenarekin

Ezin duzu:
- Helburu komertzialetarako erabili egilearen baimen espresurik gabe
- Bertsio aldatu edo eratorriak banatu

Lizentzian aurreikusi gabeko erabileretarako, jarri harremanetan egilearekin.

---

## Esker onak

- **LIFE Oreka Mendian** (LIFE15 NAT/ES/000762), HAZI Fundazioak koordinatutako proiektua, lan honen testuinguru operatiboa eskaintzeagatik.
- **GrALaren zuzendaritza**, jarraipen iteratiboagatik eta LPIaren laukoitz triangulazioaren ikuspegia proposatzeagatik.
- **Natura Kontserbatzeko Zerbitzua, Gipuzkoako Foru Aldundia**: Aralarko Natur Parkeko HIC kartografia.
- **Nafarroako Gobernua**: Aralarko San Migel estazioaren serie meteorologikoa 1991-2025.
- **GeoEuskadi**: Aralarko mendigunearen MDT5.
- **ECMWF / Copernicus / ESA**: Sentinel-2 SR Harmonized eta ERA5-Land datuak.
- **Google Earth Engine**: hodeiko konputazio-plataforma.
