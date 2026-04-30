# TFG Aralar — Diagnóstico pastoral por teledetección

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
[![DOI](https://img.shields.io/badge/DOI-pendiente-blue.svg)](https://doi.org/)
[![Status](https://img.shields.io/badge/status-en%20curso-yellow.svg)]()

Pipeline integrado para el diagnóstico pastoral del Parque Natural de Aralar (ZEC ES2120011) mediante teledetección Sentinel-2, modelización climática con GAM y triangulación cuádruple del Índice de Presión Pastoral (IPP).

**Trabajo de Fin de Grado**, Grado en Geografía y Ordenación del Territorio, UPV/EHU. Desarrollado en el marco del proyecto LIFE Oreka Mendian (LIFE15 NAT/ES/000762).

> **Web del proyecto**: https://&lt;usuario&gt;.github.io/tfg-aralar
>
> **Documentación principal**: ver `/docs/Dossier_TFG_Aralar_v3.docx`
>
> **Material didáctico**: ver `/docs/Anexo_Didactico_TFG_Aralar_v2.docx`

---

## Resumen del trabajo

El TFG diagnostica la presión pastoral en los pastos subalpinos de Aralar usando series de LAI de Sentinel-2 invertido con PROSAIL-NN, modelizadas climáticamente con GAM y triangularmente validadas por **cuatro estrategias independientes** entre 2018 y 2025. La aproximación más original es la **Estrategia D** (referencia hayedo), que usa la cubierta forestal climácica adyacente como termómetro climático para aislar la señal pastoral de la climática.

### Resultados principales

- **Patrón cuasi-bianual**: crisis pastoral en 2020, 2022 y 2025 alternada con relajación en 2023 (concordancia 4/4 entre las cuatro estrategias).
- **Gradiente espacial coherente** con la distancia a las bordas: pasto cercano (-0.09, 3/3 Alta) > intermedio (-0.10, 3/3 Alta) > remoto (-0.03, 2/3 Baja).
- **Modelo climático G6** con AIC=-189.4, R²=0.78 y validación LOYO con r medio = 0.91 sobre 8 años.
- **Diagnóstico operativo** para el seguimiento adaptativo del LIFE Oreka Mendian.

---

## Estructura del repositorio

```
tfg-aralar/
├── docs/                    Documentación del TFG (Word, PDF, flujograma)
│   ├── Dossier_TFG_Aralar_v3.docx
│   ├── Anexo_Didactico_TFG_Aralar_v2.docx
│   ├── Explicacion_Estrategia_D.docx
│   └── Flujograma_Metodologico_v2.pdf
├── scripts/
│   ├── python/              PROSAIL local + entrenamiento NN + cartografía
│   │   ├── main.py
│   │   ├── main_aralar_hayedo.py
│   │   ├── prosail_pure.py
│   │   ├── train_gee_coefficients.py
│   │   └── cartografia_pastoral_v5.py
│   ├── javascript/          Pipeline GEE
│   │   └── Pipeline_Integrado_Aralar.js
│   ├── R/                   Modelización climática y diagnóstico
│   │   ├── diagnostico_pastoral_aralar_v9_3.R
│   │   ├── Aralar_SPEI_SerieLarga_v2.R
│   │   └── comparativa_era5_vs_estacion_serie_completa.R
│   └── batch/               Lanzadores Windows
│       └── ejecutar_cartografia5.bat
├── data/                    CSVs de salida del pipeline
│   ├── diagnostico_pastoral.csv
│   ├── validacion_loyo.csv
│   ├── comparacion_modelos.csv
│   └── hayedo_serie_referencia.csv
├── maps/                    Cartografía pastoral por año
│   ├── 2018/
│   ├── 2019/
│   ├── ...
│   └── 2025/
├── figures/                 Figuras síntesis del trabajo
│   ├── heatmap_IPP_Aralar.png
│   ├── esquema_4_estrategias.png
│   ├── esquema_pipeline.png
│   └── ...
├── web/                     Web GitHub Pages
│   ├── index.html
│   ├── style.css
│   └── mapas/
├── README.md                Este archivo
└── LICENSE                  CC BY-NC-ND 4.0
```

---

## Cómo reproducir

### Requisitos

- **Python 3.10+** con paquetes: `prosail`, `scikit-learn`, `numpy`, `rasterio`, `matplotlib`, `pandas`
- **R 4.3+** con paquetes: `mgcv`, `SPEI`, `ggplot2`, `dplyr`, `tidyr`, `lubridate`
- **Cuenta Google Earth Engine** activa
- **Node.js 18+** (solo para regenerar la documentación; no necesario para el pipeline)

### Pasos del pipeline

**1. Generación de las LUTs y entrenamiento de las redes neuronales (Python local)**

```bash
cd scripts/python
python main.py            # genera LUT del pasto + entrena NN_PASTO
python main_aralar_hayedo.py  # genera LUT del hayedo + entrena NN_HAYEDO
```

Salida: `LUT_PASTO.pkl`, `LUT_HAYEDO.pkl` y los pesos de las NN listos para inyectar en GEE.

**2. Inversión a escala en GEE (JavaScript)**

Subir los assets de Aralar al proyecto GEE (límite, HIC, bordas, MDT5).
Pegar `Pipeline_Integrado_Aralar.js` en el editor de código de GEE y ejecutar.
Las exportaciones a Drive incluyen el CSV de LAI zonal, el CSV climático y los GeoTIFFs anuales.

> ⚠️ **Bug crítico documentado**: en el bloque de reducción zonal, las llamadas a `ee.Dictionary.combine()` deben llevar el segundo argumento `false` (overwrite=false). Sin esa precaución, los defaults sobrescriben los valores reales y el CSV exportado sale entero a -9999. Ver línea 487 del script para detalles.

**3. Modelización climática y cálculo de las cuatro estrategias del IPP (R)**

```r
source("scripts/R/Aralar_SPEI_SerieLarga_v2.R")    # SPEI a partir de ERA5-Land
source("scripts/R/diagnostico_pastoral_aralar_v9_3.R")  # Modelo G6 + 4 estrategias
```

Salida: `diagnostico_pastoral.csv`, `validacion_loyo.csv`, `comparacion_modelos.csv` y figuras de diagnóstico.

**4. Cartografía pastoral (Python)**

```bash
cd scripts/python
python cartografia_pastoral_v5.py --year 2022
```

O bien ejecutar `scripts/batch/ejecutar_cartografia5.bat` para procesar la serie 2018-2025 completa.

Salida en `/maps/<año>/`: 12 mapas pastorales por año (biomasa, calidad, carga, balance, riesgo, IPP_A, IPP_B/D zonal, IPP_C píxel, consenso integrado, concordancia píxel a píxel) con su correspondiente .txt con metadatos.

---

## Citación

Si usas este código o los resultados en trabajos derivados, por favor cita:

```
[Apellido, Nombre del alumno] (2026). Diagnóstico pastoral del Parque Natural de
Aralar mediante teledetección y modelización climática. Trabajo de Fin de Grado,
Grado en Geografía y Ordenación del Territorio, Universidad del País Vasco / Euskal
Herriko Unibertsitatea. https://github.com/<usuario>/tfg-aralar
```

Para la metodología específica de la Estrategia D (referencia hayedo), citar también el dossier técnico v3 (`/docs/Dossier_TFG_Aralar_v3.docx`).

---

## Licencia

Este trabajo se publica bajo licencia [Creative Commons Atribución-NoComercial-SinDerivadas 4.0 Internacional (CC BY-NC-ND 4.0)](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.es).

Esto significa que puedes:

- ✅ **Compartir** el material en cualquier medio o formato
- ✅ **Citarlo** en trabajos académicos siempre con atribución correcta

Pero **no puedes**:

- ❌ Usarlo con fines comerciales sin autorización expresa del autor
- ❌ Distribuir versiones modificadas o derivadas

Para usos no contemplados en la licencia, contactar con el autor.

---

## Limitaciones del trabajo (transparencia)

Este TFG ha sido revisado iterativamente por el director y se han identificado y corregido tres imprecisiones que es importante documentar:

1. **Modelo G6**: está entrenado con la serie agregada del AOI completo (`LAI_aoi`), NO con las zonas de pasto remoto. La pseudo-referencia interna (zonas remotas) corresponde a la Estrategia B, no a la A.
2. **Resiliencia diferencial pasto-hayedo**: el hayedo es ecológicamente más resiliente al estrés climático que el pasto (sistema radicular, reservas, microclima). Una parte del IPP_D negativo puede reflejar esta resiliencia diferencial fisiológica, no solo pastoreo. El método se mantiene defendible mediante tres argumentos: variabilidad interanual incompatible con fisiología constante, gradiente espacial imposible de explicar fisiológicamente, y convergencia con A, B, C que no usan el hayedo.
3. **R² de las redes neuronales**: el R² = 0.778 obtenido sobre la LUT de pasto es coherente con la literatura para arquitecturas comparables (MLP de una capa oculta sobre LUT PROSAIL, R² publicados entre 0.55 y 0.73 según Croci et al. 2022 y SNAP ATBD). Los R² superiores a 0.85 que aparecen en literatura corresponden a Gaussian Process Regression, método distinto.

Otras limitaciones reconocidas del trabajo:

- Serie temporal de 8 años todavía corta para variabilidad climática plurianual (literatura recomienda 20-30 años).
- Ausencia de validación in situ con datos ganaderos reales (UGM, fechas).
- Año 2024 con cobertura observacional reducida (n=9 escenas) afecta el LOYO de ese año (r=0.77).
- LUT específica para Fagus sylvatica del País Vasco sin precedente bibliográfico publicado: rangos parametrizados a partir de literatura genérica de caducifolios europeos.

---

## Agradecimientos

- **LIFE Oreka Mendian** (LIFE15 NAT/ES/000762), proyecto coordinado por HAZI Fundazioa, por proporcionar el contexto operativo de este trabajo.
- **Dirección del TFG**: por el seguimiento iterativo, la detección de errores conceptuales y la propuesta del enfoque de triangulación cuádruple del IPP.
- **Servicio de Conservación de la Naturaleza, Diputación Foral de Guipúzcoa**: cartografía HIC del Parque Natural de Aralar.
- **Gobierno de Navarra**: serie meteorológica de la estación de San Miguel de Aralar 1991-2025.
- **GeoEuskadi**: MDT5 del macizo de Aralar.
- **ECMWF / Copernicus / ESA**: datos Sentinel-2 SR Harmonized y ERA5-Land.
- **Google Earth Engine**: plataforma de cómputo en la nube.

---

## Contacto

[Nombre del alumno] — [email institucional UPV/EHU]
[Nombre del director] — [email institucional UPV/EHU]
