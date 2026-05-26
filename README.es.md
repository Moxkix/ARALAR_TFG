> 🌐 [Euskara](README.md) · **Castellano** · [English](README.en.md)

# TFG Aralar — Diagnóstico pastoral por teledetección

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
[![Status](https://img.shields.io/badge/status-en%20curso-yellow.svg)]()

Pipeline integrado para el diagnóstico pastoral del Parque Natural de Aralar (ZEC ES2120011) mediante teledetección Sentinel-2, modelización climática con GAM y triangulación cuádruple del Índice de Presión Pastoral (IPP).

**Trabajo de Fin de Grado**, Grado en Geografía y Ordenación del Territorio, UPV/EHU. Desarrollado en el marco del proyecto LIFE Oreka Mendian (LIFE15 NAT/ES/000762).

> **Web del proyecto**: https://moxkix.github.io/ARALAR_TFG/
>
> **Documentación**: la carpeta `docs/` se publica en la web automáticamente. Cualquier archivo que se suba ahí (Word, PDF, CSV, etc.) aparece como tarjeta en la sección Documentos sin tocar el HTML.

---

## Resumen del trabajo

El TFG diagnostica la presión pastoral en los pastos subalpinos de Aralar usando series de LAI de Sentinel-2 invertido con PROSAIL-NN, modelizadas climáticamente con GAM y validadas por **cuatro estrategias independientes** entre 2018 y 2025. La aproximación más original es la **Estrategia D** (referencia hayedo), que usa la cubierta forestal climácica adyacente como termómetro climático para aislar la señal pastoral de la climática.

### Resultados principales

- **Patrón cuasi-bianual**: crisis pastoral en 2020, 2022 y 2025 alternada con relajación en 2023 (concordancia 4/4 entre las cuatro estrategias).
- **Gradiente espacial coherente** con la distancia a las bordas: pasto cercano > intermedio > remoto.
- **Modelo climático G6** validado mediante LOYO (leave-one-year-out) sobre la serie 2018-2025.
- **Diagnóstico operativo** para el seguimiento adaptativo del LIFE Oreka Mendian.

---

## Estructura del repositorio

```
ARALAR_TFG/
├── index.html               Página de GitHub Pages (raíz del sitio, trilingüe EUS/CAS/ING)
├── style.css                Estilos de la web
├── README.md                Léeme principal (euskera)
├── README.es.md             Léeme en castellano (este archivo)
├── README.en.md             Léeme en inglés
├── LICENSE                  CC BY-NC-ND 4.0
│
├── docs/                    Documentación y CSVs del diagnóstico (Word, PDF, CSV…)
│                            Listado automático en la web vía API de GitHub.
│
├── data/                    Salidas tabulares del pipeline
│   ├── diagnostico_pastoral.csv
│   ├── validacion_loyo.csv
│   ├── comparacion_modelos.csv
│   ├── hayedo_serie_referencia.csv
│   └── coefficients/        LUTs PROSAIL exportadas y coeficientes NN
│                            ({PASTOS,HAYEDO}_gee.csv y *_nn_coefficients.{js,json})
│                            Generados por la web-app y consumidos por el script GEE.
│
├── figures/                 Figuras de síntesis del trabajo
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
├── maps/                    Cartografía pastoral por año
│   ├── 2018/  …  2025/      Serie principal (PNG + GeoTIFF + resumen .txt por año)
│   └── 2018_LIFE/  …  2025_LIFE/
│                            Variante etiquetada para LIFE Oreka Mendian
│
└── scripts/
    ├── python/
    │   ├── web-app/         Aplicación web FastAPI para generar LUTs PROSAIL
    │   │                    y hacer la inversión Sentinel-2 → variables biofísicas
    │   │   ├── main.py             Backend FastAPI (REST + WebSocket)
    │   │   ├── prosail_pure.py     Modelo PROSPECT-D + 4SAIL puro Python
    │   │   ├── requirements.txt    Dependencias
    │   │   ├── README.md           Guía de instalación y uso
    │   │   └── static/index.html   Frontend SPA
    │   ├── train_gee_coefficients.py   Entrenamiento de la NN para inyectar en GEE
    │   └── cartografia_pastoral_v5.py  Generación de los mapas anuales
    ├── javascripts/         Pipeline en Google Earth Engine
    │   └── GEE_Pipeline_Integrado_Aralar.js
    ├── R/                   Modelización climática y diagnóstico cuádruple
    │   └── diagnostico_pastoral_aralar_v9_4.R
    └── batch/               Lanzadores Windows
        └── ejecutar_cartografia5.bat
```

### Web-app PROSAIL (inversión Sentinel-2)

La aplicación de `scripts/python/web-app/` es un servidor FastAPI con frontend incluido. Permite configurar de forma interactiva los parámetros de la LUT (rangos PROSAIL, geometría solar, distribución de muestreo, ruido), generarla en memoria con el modelo PROSPECT-D + 4SAIL real y ejecutar la inversión sobre imágenes Sentinel-2 (BEAM-DIMAP de SNAP o GeoTIFF de 10 bandas) por mínima distancia con cKDTree.

```bash
cd scripts/python/web-app
python -m venv venv && source venv/bin/activate     # Windows: venv\Scripts\activate
pip install -r requirements.txt
python main.py
# Abrir http://localhost:8000
```

Los productos finales (LUTs exportadas y coeficientes de la red neuronal listos para inyectar en GEE) se versionan en `data/coefficients/`. La caché de runtime (`temp_prosail/`) está excluida del repositorio.

### Productos cartográficos por año

Cada carpeta `maps/<año>/` contiene los siguientes mapas en PNG (nomenclatura `aralar_mapa_*_<año>.png` / `aralar_pastos_mapa_*_<año>.png`), junto con los GeoTIFF (`.tif`) de cada variable:

| Familia | Mapas |
|---|---|
| Triangulación IPP | `IPP_A_pastos`, `IPP_A_todos`, `IPP_B_pastos`, `IPP_C_pastos` |
| Síntesis | `IPP_consenso_pastos`, `IPP_concordancia_pastos`, `IPP_espacializado` |
| Gestión | `mapa_biomasa`, `mapa_calidad`, `mapa_carga`, `mapa_balance`, `mapa_riesgo` |
| Control | `mapa_pendiente` |

La **Estrategia D** se integra en el diagnóstico y en la figura de triangulación, pero no se exporta como mapa espacial independiente. Además se genera un `aralar_resumen_<año>.txt` con los metadatos del año.

---

## Cómo reproducir

### Requisitos

- **Python 3.10+** con: `prosail`, `scikit-learn`, `numpy`, `rasterio`, `matplotlib`, `pandas`
- **R 4.3+** con: `mgcv`, `SPEI`, `ggplot2`, `dplyr`, `tidyr`, `lubridate`
- **Cuenta Google Earth Engine** activa

### Pasos del pipeline

**1. LUTs y entrenamiento de redes neuronales (Python local)**

```bash
# Generación de la LUT con la web-app (interfaz web en localhost:8000)
cd scripts/python/web-app
pip install -r requirements.txt
python main.py

# Entrenamiento de la NN sobre las LUTs exportadas y exportación de coeficientes para GEE
cd ../
python train_gee_coefficients.py
```

Las LUTs exportadas y los coeficientes finales quedan en `data/coefficients/`.

**2. Inversión a escala en GEE (JavaScript)**

Subir los assets de Aralar al proyecto GEE (límite, HIC, bordas, MDT5).
Pegar `scripts/javascripts/GEE_Pipeline_Integrado_Aralar.js` en el editor de código y ejecutar. Las exportaciones a Drive incluyen los CSV de LAI zonal, los CSV climáticos y los GeoTIFF anuales.

> ⚠️ **Bug crítico documentado**: en la reducción zonal, las llamadas a `ee.Dictionary.combine()` deben llevar el segundo argumento `false` (overwrite=false). Sin esa precaución, los defaults sobrescriben los valores reales y el CSV exportado sale entero a -9999.

**3. Modelización climática y cuatro estrategias del IPP (R)**

```r
source("scripts/R/diagnostico_pastoral_aralar_v9_4.R")
```

Salida: `data/diagnostico_pastoral.csv`, `data/validacion_loyo.csv`, `data/comparacion_modelos.csv` y las figuras de diagnóstico en `figures/` (gráficos e informe en euskera).

**4. Cartografía pastoral (Python)**

```bash
cd scripts/python
python cartografia_pastoral_v5.py --year 2022
```

O ejecutar `scripts/batch/ejecutar_cartografia5.bat` para procesar la serie 2018-2025 completa. Salida en `maps/<año>/` (mapas con títulos y leyendas en euskera).

---

## Web del proyecto

La web (`index.html` + `style.css`) se publica con GitHub Pages desde la rama `main` y la raíz del repositorio. Es **trilingüe (euskera / castellano / inglés)** con un selector en la cabecera; la elección se recuerda en el navegador. El listado de las secciones Documentos y Gráficos se construye al cargar la página consultando la API de GitHub (`/repos/Moxkix/ARALAR_TFG/contents/{docs,figures}`). Esto significa que **basta con subir un archivo a esas carpetas para que aparezca**, sin tocar el HTML.

---

## Citación

Si usas este código o los resultados en trabajos derivados, por favor cita:

```
[Apellido, Nombre del alumno] (2026). Diagnóstico pastoral del Parque Natural
de Aralar mediante teledetección y modelización climática. Trabajo de Fin de
Grado, Grado en Geografía y Ordenación del Territorio, Universidad del País
Vasco / Euskal Herriko Unibertsitatea.
https://github.com/Moxkix/ARALAR_TFG
```

---

## Licencia

Este trabajo se publica bajo licencia [Creative Commons Atribución-NoComercial-SinDerivadas 4.0 Internacional (CC BY-NC-ND 4.0)](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.es).

Puedes:
- Compartir el material en cualquier medio o formato
- Citarlo en trabajos académicos siempre con atribución correcta

No puedes:
- Usarlo con fines comerciales sin autorización expresa del autor
- Distribuir versiones modificadas o derivadas

Para usos no contemplados en la licencia, contactar con el autor.

---

## Agradecimientos

- **LIFE Oreka Mendian** (LIFE15 NAT/ES/000762), proyecto coordinado por HAZI Fundazioa, por proporcionar el contexto operativo de este trabajo.
- **Dirección del TFG**, por el seguimiento iterativo y la propuesta del enfoque de triangulación cuádruple del IPP.
- **Servicio de Conservación de la Naturaleza, Diputación Foral de Guipúzcoa**: cartografía HIC del Parque Natural de Aralar.
- **Gobierno de Navarra**: serie meteorológica de la estación de San Miguel de Aralar 1991-2025.
- **GeoEuskadi**: MDT5 del macizo de Aralar.
- **ECMWF / Copernicus / ESA**: datos Sentinel-2 SR Harmonized y ERA5-Land.
- **Google Earth Engine**: plataforma de cómputo en la nube.
