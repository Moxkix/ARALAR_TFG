#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CARTOGRAFÍA PASTORAL v5 — Sierra de Aralar (ZEC ES2120011)
Mapas píxel a píxel restringidos a hábitats pastorales (HIC TIPO=1).

NOVEDADES v5 (Nivel 2 — abril 2026):
  - Soporte para Estrategia D (referencia hayedo) además de A y C
  - Cálculo automático de IPP_consenso y concordancia_n cuando se aportan
    los rasters de las cuatro estrategias
  - Etiqueta del modelo climático actualizada a G6: s(GDA)+s(P60)+s(SM_root)
  - Detección de la banda modelo_inv del pipeline GEE dual (PROSAIL pasto +
    PROSAIL hayedo) y mapa diagnóstico opcional

Productos (8-12 mapas + 9-13 GeoTIFFs por año, según estrategias disponibles):
  1. Biomasa forrajera (kg MS/ha) — solo pastos
  2. Calidad pastoral (5 clases) — solo pastos
  3. Capacidad de acogida ganadera (UGM/ha) — solo pastos
  4. Balance forrajero (oferta − demanda) — solo pastos
  5. Diagnóstico de riesgo/desequilibrio — solo pastos, con pendiente 5m
  6. Pendiente — solo pastos
  7. IPP espacializado Estrategia C (todos los hábitats / solo pastos)
  8. IPP Estrategia A si se aporta --ipp_a (todos / solo pastos)  [NUEVO v5: G6]
  9. IPP Estrategia D si se aporta --ipp_d (todos / solo pastos)  [NUEVO v5]
 10. IPP consenso y concordancia_n si están las 4 estrategias    [NUEVO v5]
 11. Mapa de modelo invertido (pasto vs hayedo) si la banda existe [NUEVO v5]

Uso:
  python cartografia_pastoral_v5.py ^
    --prosail PROSAIL_NN_Aralar_2025.tif ^
    --hic HIC_tipo_aralar.tif ^
    --bordas dist_borda_aralar.tif ^
    --dem MDT05_aralar.tif ^
    --output_dir ./mapas/2025/ ^
    --carga 0.5 --temporada 180 ^
    --ipp_a IPP_estrategia_A_2025.tif ^
    --ipp_d IPP_estrategia_D_2025.tif
"""

import argparse, os, sys, csv, datetime
import numpy as np


class TeeLogger:
    """Duplica stdout a un fichero de texto. Permite seguir usando print()
    en el resto del script sin cambios."""
    def __init__(self, path, mode='w'):
        self.terminal = sys.stdout
        self.log = open(path, mode, encoding='utf-8')
        self.path = path
    def write(self, msg):
        self.terminal.write(msg)
        self.log.write(msg)
        self.log.flush()
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    def close(self):
        try: self.log.close()
        except Exception: pass

try:
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import reproject
except ImportError:
    sys.exit("ERROR: pip install rasterio --break-system-packages")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import Patch
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    sys.exit("ERROR: pip install matplotlib --break-system-packages")


# ══════════════════════════════════════════════════════
# CONFIGURACIÓN
# ══════════════════════════════════════════════════════

class Config:
    # ── LAI → Biomasa (pastos atlánticos de montaña) ──
    BIO_A = 520.0          # kg MS/ha por unidad LAI
    BIO_B = 1.10           # exponente alométrico
    BIO_MAX = 3500.0       # techo biomasa en pie (kg MS/ha)

    # ── Renovación del pasto (turnover) ──
    # En pastos atlánticos de montaña, la biomasa se renueva
    # 2.5-3.5 veces durante la temporada de pastoreo.
    # Aldezábal & Moragues (2020), Canals & Sebastià (2000)
    TURNOVER = 3.0

    # ── Capacidad de acogida ──
    UGM_KG_DIA = 12.0     # kg MS/día/UGM
    APROVECH = 0.55        # fracción palatable en pasto de montaña
    TEMPORADA = 180        # días mayo–octubre (sistema mixto)

    # ── Balance ──
    CARGA_REF = 0.5        # UGM/ha carga media real estimada

    # ── Riesgo ──
    SLOPE_EROSION = 25.0   # grados — umbral erosión

    # ── Anillos de distancia a bordas ──
    DIST_CERCANO = 500
    DIST_INTERMEDIO = 1500


# ══════════════════════════════════════════════════════
# NUEVO v5 — CARGA Y ESPACIALIZACIÓN DEL CSV ZONAL
# ══════════════════════════════════════════════════════

# Mapeo nombre de zona (script R) → (HIC TIPO, anillo)
ZONA_NAME_MAP = {
    'pasto_cercano':    (1, 1),  'pasto_intermedio': (1, 2),  'pasto_remoto': (1, 3),
    'hayedo_cercano':   (3, 1),  'hayedo_intermedio': (3, 2), 'hayedo_remoto': (3, 3),
    'brezal_cercano':   (2, 1),  'brezal_intermedio': (2, 2), 'brezal_remoto': (2, 3),
    'encinar_cercano':  (4, 1),  'encinar_intermedio': (4, 2),'encinar_remoto': (4, 3),
}
# También por código numérico (11, 12, 13, 31, 32, 33...)
def _zona_code_to_tuple(z):
    s = str(z).strip()
    if s in ZONA_NAME_MAP: return ZONA_NAME_MAP[s]
    if s.isdigit() and len(s) == 2:
        return (int(s[0]), int(s[1]))
    return None


def load_ipp_zonal_csv(csv_path, year):
    """Lee el CSV de diagnóstico zonal del script R v9 y devuelve un dict
    { (tipo, anillo): {'IPP_A':float, 'IPP_B':..., 'IPP_C':..., 'IPP_D':...,
                       'IPP_consenso':..., 'concord_k':int, 'concord_n':int} }
    Solo para las filas correspondientes al año pedido. Tolera diferentes
    nombres de columna (year/anyo/anio/AÑO) y de zona.
    """
    if not csv_path or not os.path.exists(csv_path):
        return None
    out = {}
    try:
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            # Detectar separador
            sample = f.read(4096); f.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=',;\t')
            except csv.Error:
                dialect = csv.excel
            reader = csv.DictReader(f, dialect=dialect)
            # Normalizar nombres de columna
            cols_lower = {c.lower().strip(): c for c in reader.fieldnames}
            year_col = next((cols_lower[k] for k in ('year','anyo','anio','año','ano') if k in cols_lower), None)
            zona_col = next((cols_lower[k] for k in ('zona','zone','code_zona','codigo','code') if k in cols_lower), None)
            if zona_col is None:
                print(f"  AVISO: CSV {csv_path} no tiene columna 'zona' reconocible")
                return None
            def get(row, *keys):
                for k in keys:
                    if k in cols_lower and row.get(cols_lower[k], '') not in ('', 'NA', 'NaN', None):
                        try: return float(row[cols_lower[k]])
                        except ValueError: return None
                return None
            for row in reader:
                if year_col and str(row.get(year_col,'')).strip() not in (str(year), year):
                    continue
                ztup = _zona_code_to_tuple(row[zona_col])
                if ztup is None: continue
                # Parsear concordancia "k/n"
                ck = cn = None
                conc_raw = row.get(cols_lower.get('concordancia_n',''),
                              row.get(cols_lower.get('concordancia',''), ''))
                if isinstance(conc_raw, str) and '/' in conc_raw:
                    try:
                        ck, cn = [int(x) for x in conc_raw.split('/')]
                    except ValueError:
                        ck = cn = None
                out[ztup] = {
                    'IPP_A': get(row, 'ipp_a'),
                    'IPP_B': get(row, 'ipp_b'),
                    'IPP_C': get(row, 'ipp_c'),
                    'IPP_D': get(row, 'ipp_d'),
                    'IPP_consenso': get(row, 'ipp_consenso','consenso'),
                    'concord_k': ck, 'concord_n': cn,
                }
    except Exception as e:
        print(f"  AVISO: error leyendo CSV zonal: {e}")
        return None
    if not out:
        print(f"  AVISO: CSV cargado pero sin filas para el año {year}")
        return None
    print(f"  CSV zonal cargado: {len(out)} zonas para el año {year}")
    return out


def spatialize_from_csv(zonal_dict, key, hic_tipo, anillo, valid):
    """Construye un raster aplicando el valor zonal de `key` (ej 'IPP_D')
    a todos los píxeles de cada zona (HIC TIPO × anillo). Devuelve None
    si ninguna zona tiene valor para esa key."""
    raster = np.full(hic_tipo.shape, np.nan, dtype=np.float32)
    used = 0
    for (tipo, ani), vals in zonal_dict.items():
        v = vals.get(key)
        if v is None or (isinstance(v, float) and np.isnan(v)): continue
        m = (hic_tipo == tipo) & (anillo == ani) & valid
        if np.any(m):
            raster[m] = v
            used += 1
    if used == 0: return None
    return raster



# ══════════════════════════════════════════════════════

def resample_to_match(src_path, ref_meta, ref_transform, ref_shape, method=Resampling.nearest):
    """Remuestrea un raster para que coincida con la rejilla de referencia."""
    with rasterio.open(src_path) as src:
        data = np.empty((1, ref_shape[0], ref_shape[1]), dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=data[0],
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_transform,
            dst_crs=ref_meta['crs'],
            resampling=method
        )
    return data[0]


# ══════════════════════════════════════════════════════
# CÁLCULOS PASTORALES
# ══════════════════════════════════════════════════════

def lai_to_biomass(lai):
    """LAI → biomasa forrajera en pie (kg MS/ha)."""
    b = Config.BIO_A * np.power(np.maximum(lai, 0), Config.BIO_B)
    return np.minimum(b, Config.BIO_MAX)


def compute_quality(lai, cab):
    """Calidad pastoral 1-5 basada en LAI + clorofila."""
    q = np.ones_like(lai, dtype=np.int8)
    q[(lai >= 0.3) & (cab >= 8)] = 2   # pobre
    q[(lai >= 0.6) & (cab >= 18)] = 3  # moderado
    q[(lai >= 1.0) & (lai <= 4.5) & (cab >= 30)] = 4  # bueno
    q[(lai >= 1.5) & (lai <= 3.5) & (cab >= 45)] = 5  # excelente
    return q


def biomass_to_ugm(biomass):
    """Biomasa en pie → capacidad de acogida (UGM/ha) con turnover."""
    produccion_total = biomass * Config.TURNOVER
    util = produccion_total * Config.APROVECH
    demanda = Config.UGM_KG_DIA * Config.TEMPORADA
    return util / demanda


def compute_risk(balance, slope):
    """Diagnóstico de riesgo en 7 clases."""
    r = np.full_like(balance, 4, dtype=np.int8)  # 4 = equilibrio
    r[balance > 0.12] = 3   # subpastoreo leve
    r[balance > 0.35] = 2   # matorralización
    r[balance < -0.08] = 5  # presión elevada
    r[balance < -0.20] = 6  # degradación
    if slope is not None:
        r[(balance < -0.20) & (slope > Config.SLOPE_EROSION)] = 7
    return r


def compute_spatial_ipp(lai, hic_tipo, dist_borda):
    """IPP espacializado: anomalía LAI por zona HIC × anillo."""
    ipp = np.full_like(lai, np.nan, dtype=np.float32)

    anillo = np.zeros_like(dist_borda, dtype=np.int8)
    anillo[dist_borda < Config.DIST_CERCANO] = 1
    anillo[(dist_borda >= Config.DIST_CERCANO) & (dist_borda < Config.DIST_INTERMEDIO)] = 2
    anillo[dist_borda >= Config.DIST_INTERMEDIO] = 3

    for tipo in [1, 2, 3, 4]:
        for ani in [1, 2, 3]:
            mask = (hic_tipo == tipo) & (anillo == ani) & ~np.isnan(lai) & (lai > 0)
            if np.sum(mask) < 10:
                continue
            med = np.median(lai[mask])
            if med > 0:
                ipp[mask] = (lai[mask] - med) / med

    return ipp, anillo


# ══════════════════════════════════════════════════════
# VISUALIZACIÓN
# ══════════════════════════════════════════════════════

def setup_plot():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.size': 11,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.dpi': 200, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15,
    })


def get_extent(transform, shape):
    return [transform[2], transform[2] + shape[1] * transform[0],
            transform[5] + shape[0] * transform[4], transform[5]]


def plot_continuous(data, transform, title, label, cmap, vmin, vmax, outpath,
                    hic_tipo=None, pasto_only=False, stats=True):
    setup_plot()
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    extent = get_extent(transform, data.shape)

    if pasto_only and hic_tipo is not None:
        bg = np.where((hic_tipo > 0) & (hic_tipo != 1), 0.85, np.nan)
        ax.imshow(bg, extent=extent, cmap='gray', vmin=0, vmax=1,
                  interpolation='nearest', aspect='equal', alpha=0.4)

    masked = np.ma.masked_where(np.isnan(data), data)
    im = ax.imshow(masked, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation='nearest', aspect='equal')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    plt.colorbar(im, cax=cax, label=label)

    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel('Ekialdea (m)')
    ax.set_ylabel('Iparraldea (m)')
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(6, 6))

    if stats:
        valid = data[~np.isnan(data)]
        if len(valid) > 0:
            txt = (f"Batezbestekoa: {np.mean(valid):.2f}\nMediana: {np.median(valid):.2f}\n"
                   f"Desbid.: {np.std(valid):.2f}\nMin: {np.min(valid):.2f}\n"
                   f"Max: {np.max(valid):.2f}")
            ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=9,
                    va='bottom', bbox=dict(boxstyle='round', fc='white', alpha=0.85))

    plt.savefig(outpath)
    plt.close()
    print(f"  -> {outpath}")


def plot_classified(data, transform, title, classes, colors, outpath,
                    hic_tipo=None, pasto_only=False):
    setup_plot()
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    extent = get_extent(transform, data.shape)

    if pasto_only and hic_tipo is not None:
        bg = np.where((hic_tipo > 0) & (hic_tipo != 1), 0.85, np.nan)
        ax.imshow(bg, extent=extent, cmap='gray', vmin=0, vmax=1,
                  interpolation='nearest', aspect='equal', alpha=0.4)

    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(-0.5, len(colors) + 0.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    masked = np.ma.masked_where((data == 0) | np.isnan(data.astype(float)), data - 1)
    ax.imshow(masked, extent=extent, cmap=cmap, norm=norm,
              interpolation='nearest', aspect='equal')

    patches = [Patch(facecolor=c, edgecolor='gray', linewidth=0.5, label=l)
               for c, l in zip(colors, classes)]
    ax.legend(handles=patches, loc='upper right', fontsize=9, framealpha=0.9,
              edgecolor='gray', bbox_to_anchor=(1.42, 1.0))

    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel('Ekialdea (m)')
    ax.set_ylabel('Iparraldea (m)')
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(6, 6))

    plt.savefig(outpath)
    plt.close()
    print(f"  -> {outpath}")


def plot_ipp_spatial(ipp, transform, title, outpath, hic_tipo=None):
    setup_plot()
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    extent = get_extent(transform, ipp.shape)

    if hic_tipo is not None:
        bg = np.where(hic_tipo > 0, 0.92, np.nan)
        ax.imshow(bg, extent=extent, cmap='gray', vmin=0, vmax=1,
                  interpolation='nearest', aspect='equal', alpha=0.3)

    masked = np.ma.masked_where(np.isnan(ipp), ipp)
    im = ax.imshow(masked, extent=extent, cmap='RdYlGn', vmin=-0.4, vmax=0.4,
                   interpolation='nearest', aspect='equal')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('LPI (gorria = gainlarratzea | berdea = azpierabilera)')

    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel('Ekialdea (m)')
    ax.set_ylabel('Iparraldea (m)')
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(6, 6))

    if hic_tipo is not None:
        nombres = {1: 'Larrea', 2: 'Txilardia', 3: 'Pagadia', 4: 'Artadia'}
        lines = []
        for tipo, nombre in nombres.items():
            m = (hic_tipo == tipo) & ~np.isnan(ipp)
            if np.sum(m) > 0:
                lines.append(f"{nombre}: LPI={np.mean(ipp[m]):+.3f} (n={np.sum(m):,})")
        if lines:
            ax.text(0.02, 0.02, '\n'.join(lines), transform=ax.transAxes, fontsize=9,
                    va='bottom', bbox=dict(boxstyle='round', fc='white', alpha=0.85))

    plt.savefig(outpath)
    plt.close()
    print(f"  -> {outpath}")


# ══════════════════════════════════════════════════════
# PROCESAMIENTO PRINCIPAL
# ══════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Cartografia pastoral v3 — Aralar")
    parser.add_argument('--prosail', required=True, help='GeoTIFF PROSAIL-NN')
    parser.add_argument('--hic', required=True, help='Raster HIC TIPO')
    parser.add_argument('--bordas', required=True, help='Raster distancia a bordas (m)')
    parser.add_argument('--dem', required=True, help='DEM 5m GeoEuskadi')
    parser.add_argument('--ipp_a', default=None, help='GeoTIFF IPP Estrategia A (GAM clim. G6) de GEE/R [opcional]')
    parser.add_argument('--ipp_d', default=None, help='GeoTIFF IPP Estrategia D (referencia hayedo) de R [opcional, NUEVO v5]')
    parser.add_argument('--ipp_b', default=None, help='GeoTIFF IPP Estrategia B (pseudo-referencia) [opcional]')
    parser.add_argument('--ipp_csv', default=None,
        help='CSV de diagnostico zonal del script R v9 (cols: zona, year/anyo, IPP_A/B/C/D, '
             'IPP_consenso, concordancia_n). Si se aporta, espacializa B, D y consenso por zona '
             'cuando no haya GeoTIFFs propios. [NUEVO v5]')
    parser.add_argument('--output_dir', default='.', help='Directorio salida')
    parser.add_argument('--carga', type=float, default=Config.CARGA_REF, help='UGM/ha referencia')
    parser.add_argument('--temporada', type=int, default=Config.TEMPORADA, help='Dias de pastoreo')
    args = parser.parse_args()

    Config.CARGA_REF = args.carga
    Config.TEMPORADA = args.temporada
    os.makedirs(args.output_dir, exist_ok=True)

    # Extraer año del nombre de fichero
    year = ''.join(c for c in os.path.basename(args.prosail) if c.isdigit())[-4:]

    # NUEVO v5: logging dual a fichero TXT
    log_path = os.path.join(args.output_dir, f"aralar_resumen_{year}.txt")
    sys.stdout = TeeLogger(log_path, mode='w')
    print(f"# Cartografia pastoral v5 — Aralar {year}")
    print(f"# Generado: {datetime.datetime.now().isoformat(timespec='seconds')}")
    print(f"# Comando: python {' '.join(sys.argv)}")
    print(f"# Log: {log_path}")

    print(f"\n{'='*65}")
    print(f"  CARTOGRAFIA PASTORAL v5 — Aralar {year}")
    print(f"  Carga ref: {Config.CARGA_REF} UGM/ha | Temporada: {Config.TEMPORADA} dias")
    print(f"  Turnover: {Config.TURNOVER} | Aprovechamiento: {Config.APROVECH}")
    print(f"  Estrategias IPP disponibles: C (siempre)" +
          (" + A" if args.ipp_a else "") +
          (" + B" if args.ipp_b else "") +
          (" + D" if args.ipp_d else ""))
    print(f"{'='*65}")

    # ── 1. Cargar PROSAIL ──
    print("\n1. Cargando PROSAIL GeoTIFF...")
    with rasterio.open(args.prosail) as src:
        ref_meta = src.meta.copy()
        ref_transform = src.transform
        ref_shape = (src.height, src.width)

        band_names = ['LAI', 'Cab', 'fCover', 'CCC']
        bands = {}
        for i, name in enumerate(band_names):
            if i < src.count:
                d = src.read(i + 1).astype(np.float32)
                d[(d == 0) | (d < -900)] = np.nan
                bands[name] = d

        # NUEVO v5: banda modelo_inv del pipeline GEE dual (1=pasto, 3=hayedo)
        modelo_inv = None
        if src.count >= 5:
            modelo_inv = src.read(5).astype(np.int8)
            n_pasto_inv = int(np.sum(modelo_inv == 1))
            n_hayedo_inv = int(np.sum(modelo_inv == 3))
            print(f"   Banda modelo_inv detectada (pipeline dual):")
            print(f"     {n_pasto_inv:,} px invertidos con NN_PASTO")
            print(f"     {n_hayedo_inv:,} px invertidos con NN_HAYEDO")

    lai = bands['LAI']
    cab = bands.get('Cab', np.full_like(lai, 30.0))
    valid = ~np.isnan(lai) & (lai > 0)
    print(f"   Pixeles validos: {np.sum(valid):,}")
    print(f"   LAI rango: [{np.nanmin(lai[valid]):.2f}, {np.nanmax(lai[valid]):.2f}]")

    # ── 2. Cargar capas auxiliares ──
    print("\n2. Cargando capas auxiliares...")

    print(f"   HIC: {args.hic}")
    hic_tipo = resample_to_match(args.hic, ref_meta, ref_transform, ref_shape,
                                  Resampling.nearest).astype(np.int8)
    for t, n in [(1, 'Pasto'), (2, 'Brezal'), (3, 'Hayedo'), (4, 'Encinar')]:
        cnt = np.sum(hic_tipo == t)
        if cnt > 0:
            print(f"     TIPO {t} ({n}): {cnt:,} px ({cnt * 400 / 10000:.0f} ha)")

    print(f"   Bordas: {args.bordas}")
    dist_borda = resample_to_match(args.bordas, ref_meta, ref_transform, ref_shape,
                                    Resampling.bilinear)

    print(f"   DEM 5m: {args.dem}")
    dem = resample_to_match(args.dem, ref_meta, ref_transform, ref_shape,
                             Resampling.bilinear)
    # Pendiente con protección overflow
    res = abs(ref_transform[0])
    dem64 = dem.astype(np.float64)
    dem64[dem64 < -100] = np.nan
    dy, dx = np.gradient(dem64, res)
    slope = np.degrees(np.arctan(np.clip(np.sqrt(dx**2 + dy**2), 0, 100))).astype(np.float32)
    print(f"   Pendiente rango: [{np.nanmin(slope):.1f}, {np.nanmax(slope):.1f}] grados")

    # ── 3. Máscara pastoral ──
    pasto_mask = (hic_tipo == 1) & valid
    n_pasto = np.sum(pasto_mask)
    print(f"\n   Pixeles PASTO validos: {n_pasto:,} ({n_pasto * 400 / 10000:.0f} ha)")

    if n_pasto == 0:
        sys.exit("ERROR: No hay pixeles de pasto (TIPO=1).")

    # ══════════════════════════════════════════════════
    # PRODUCTOS CARTOGRÁFICOS
    # ══════════════════════════════════════════════════

    # ── 4. BIOMASA FORRAJERA ──
    print(f"\n3. Biomasa forrajera {year} (solo pastos)...")
    biomass = np.full_like(lai, np.nan)
    biomass[pasto_mask] = lai_to_biomass(lai[pasto_mask])

    plot_continuous(biomass, ref_transform,
        f"Belar-biomasa estimatua — Aralarko larreak {year}",
        "kg ML/ha", 'YlGnBu', 0, Config.BIO_MAX,
        os.path.join(args.output_dir, f"aralar_pastos_mapa_biomasa_{year}.png"),
        hic_tipo=hic_tipo, pasto_only=True)

    # ── 5. CALIDAD PASTORAL ──
    print(f"4. Calidad pastoral {year} (solo pastos)...")
    quality = np.zeros_like(lai, dtype=np.int8)
    quality[pasto_mask] = compute_quality(lai[pasto_mask], cab[pasto_mask])

    plot_classified(quality, ref_transform,
        f"Larre-kalitatea — Aralarko larreak {year}",
        ['Oso urria', 'Urria', 'Ertaina', 'Ona', 'Bikaina'],
        ['#d32f2f', '#ff9800', '#fdd835', '#4caf50', '#1b5e20'],
        os.path.join(args.output_dir, f"aralar_pastos_mapa_calidad_{year}.png"),
        hic_tipo=hic_tipo, pasto_only=True)

    # ── 6. CAPACIDAD DE ACOGIDA ──
    print(f"5. Capacidad de acogida {year} (solo pastos, turnover={Config.TURNOVER})...")
    capacity = np.full_like(lai, np.nan)
    capacity[pasto_mask] = biomass_to_ugm(biomass[pasto_mask])

    plot_continuous(capacity, ref_transform,
        f"Abere-hartzeko gaitasuna — Aralarko larreak {year}\n"
        f"Birsorkuntza={Config.TURNOVER} | Aprob.={Config.APROVECH} | {Config.TEMPORADA} egun",
        "ALU/ha", 'RdYlGn', 0, 2.0,
        os.path.join(args.output_dir, f"aralar_pastos_mapa_carga_{year}.png"),
        hic_tipo=hic_tipo, pasto_only=True)

    # ── 7. BALANCE FORRAJERO ──
    print(f"6. Balance forrajero {year} (demanda={Config.CARGA_REF} UGM/ha)...")
    balance = np.full_like(lai, np.nan)
    balance[pasto_mask] = capacity[pasto_mask] - Config.CARGA_REF

    plot_continuous(balance, ref_transform,
        f"Belar-balantzea — Aralarko larreak {year}\n"
        f"Eskaria={Config.CARGA_REF} ALU/ha | Birsorkuntza={Config.TURNOVER} | {Config.TEMPORADA} egun",
        "Balantzea (ALU/ha)", 'RdYlGn', -0.5, 0.5,
        os.path.join(args.output_dir, f"aralar_pastos_mapa_balance_{year}.png"),
        hic_tipo=hic_tipo, pasto_only=True)

    # ── 8. DIAGNÓSTICO DE RIESGO ──
    print(f"7. Diagnostico de riesgo {year} (pendiente DEM 5m)...")
    risk = np.zeros_like(lai, dtype=np.int8)
    risk[pasto_mask] = compute_risk(balance[pasto_mask], slope[pasto_mask])

    plot_classified(risk, ref_transform,
        f"Desorekaren diagnostikoa — Aralarko larreak {year}\n"
        f"DEM 5m GeoEuskadi | Erref. karga {Config.CARGA_REF} ALU/ha",
        ['Sastrakatzea (azpilarratze larria)',
         'Azpilarratze arina', 'Larre-oreka',
         'Presio handia', 'Degradazioa (gainlarratze larria)',
         'Higadura (gainlarratzea + malda)'],
        ['#0d47a1', '#64b5f6', '#4caf50', '#ff9800', '#d32f2f', '#3e2723'],
        os.path.join(args.output_dir, f"aralar_pastos_mapa_riesgo_{year}.png"),
        hic_tipo=hic_tipo, pasto_only=True)

    # ── 9. PENDIENTE ──
    print(f"8. Mapa de pendientes {year} (DEM 5m)...")
    slope_pasto = np.full_like(lai, np.nan)
    slope_pasto[pasto_mask] = slope[pasto_mask]

    plot_continuous(slope_pasto, ref_transform,
        f"Malda — Aralarko larreak {year} (DEM 5m GeoEuskadi)",
        "Malda (graduak)", 'YlOrRd', 0, 45,
        os.path.join(args.output_dir, f"aralar_pastos_mapa_pendiente_{year}.png"),
        hic_tipo=hic_tipo, pasto_only=True)

    # ── 10. IPP ESPACIALIZADO ──
    print(f"\n9. IPP espacializado {year} (todos los habitats)...")
    ipp, anillo = compute_spatial_ipp(lai, hic_tipo, dist_borda)

    plot_ipp_spatial(ipp, ref_transform,
        f"LPI espazializatua — Aralar {year}\n"
        f"LAI anomalia habitataren eta bordetarainoko distantzia-eraztunaren arabera",
        os.path.join(args.output_dir, f"aralar_mapa_IPP_espacializado_{year}.png"),
        hic_tipo=hic_tipo)

    # IPP solo pastos
    ipp_pasto = ipp.copy()
    ipp_pasto[hic_tipo != 1] = np.nan

    plot_ipp_spatial(ipp_pasto, ref_transform,
        f"LPI espazializatua (C estr.) — Aralarko larreak soilik {year}\n"
        f"LAI anomalia bordetarainoko distantzia-eraztunaren arabera",
        os.path.join(args.output_dir, f"aralar_mapa_IPP_C_pastos_{year}.png"),
        hic_tipo=hic_tipo)

    # ── NUEVO v5: cargar CSV zonal del script R si se ha aportado ──
    zonal_dict = None
    if args.ipp_csv:
        print(f"\n  Cargando CSV zonal: {args.ipp_csv}")
        zonal_dict = load_ipp_zonal_csv(args.ipp_csv, year)

    # ── 10b. IPP ESTRATEGIA A (con corrección climática GAM G6) ──
    ipp_a = None
    ipp_a_pasto = None
    
    if args.ipp_a and os.path.exists(args.ipp_a):
        print(f"\n9b. IPP Estrategia A {year} (correccion climatica GAM G6)...")
        ipp_a = resample_to_match(args.ipp_a, ref_meta, ref_transform, ref_shape,
                                   Resampling.bilinear)
        ipp_a[~valid | (hic_tipo == 0)] = np.nan

        plot_ipp_spatial(ipp_a, ref_transform,
            f"LPI A estrategia (G6 zuzenketa klimatikoa) — Aralar {year}\n"
            f"LAI_norm hondarra GAM ereduarekiko: s(GDA) + s(P60) + s(SM_root)",
            os.path.join(args.output_dir, f"aralar_mapa_IPP_A_todos_{year}.png"),
            hic_tipo=hic_tipo)

        ipp_a_pasto = ipp_a.copy()
        ipp_a_pasto[hic_tipo != 1] = np.nan

        plot_ipp_spatial(ipp_a_pasto, ref_transform,
            f"LPI A estrategia (G6 zuzenketa klimatikoa) — Aralarko larreak soilik {year}\n"
            f"LAI_norm hondarra GAM ereduarekiko: s(GDA) + s(P60) + s(SM_root)",
            os.path.join(args.output_dir, f"aralar_mapa_IPP_A_pastos_{year}.png"),
            hic_tipo=hic_tipo)

        # NOTA v5: la comparativa entre estrategias se hace de forma unificada
        # en el bloque 10e, donde se calcula consenso + concordancia con todas
        # las estrategias disponibles (A, B, C, D).
    else:
        if args.ipp_a:
            print(f"\n  AVISO: No se encontro {args.ipp_a}")
        print(f"  IPP Estrategia A no disponible.")

    # ── 10c. IPP ESTRATEGIA B (pseudo-referencia interna) — NUEVO v5 ──
    ipp_b = None
    ipp_b_pasto = None
    if args.ipp_b and os.path.exists(args.ipp_b):
        print(f"\n9c. IPP Estrategia B {year} (pseudo-referencia interna, desde GeoTIFF)...")
        ipp_b = resample_to_match(args.ipp_b, ref_meta, ref_transform, ref_shape,
                                   Resampling.bilinear)
        ipp_b[~valid | (hic_tipo == 0)] = np.nan
    elif zonal_dict is not None:
        ipp_b = spatialize_from_csv(zonal_dict, 'IPP_B', hic_tipo, anillo, valid)
        if ipp_b is not None:
            print(f"\n9c. IPP Estrategia B {year} (espacializado por zona desde CSV R)...")

    if ipp_b is not None:
        ipp_b_pasto = ipp_b.copy()
        ipp_b_pasto[hic_tipo != 1] = np.nan
        plot_ipp_spatial(ipp_b_pasto, ref_transform,
            f"LPI B estrategia (pseudo-erreferentzia) — Aralarko larreak soilik {year}\n"
            f"Larre urruna + malda handia oinarri gisa erabilita entrenatutako eredua",
            os.path.join(args.output_dir, f"aralar_mapa_IPP_B_pastos_{year}.png"),
            hic_tipo=hic_tipo)

    # ── 10d. IPP ESTRATEGIA D (referencia hayedo) — NUEVO v5 ──
    ipp_d = None
    ipp_d_pasto = None
    src_d = None  # 'tif' o 'csv'
    if args.ipp_d and os.path.exists(args.ipp_d):
        print(f"\n9d. IPP Estrategia D {year} (referencia hayedo, desde GeoTIFF)...")
        ipp_d = resample_to_match(args.ipp_d, ref_meta, ref_transform, ref_shape,
                                   Resampling.bilinear)
        ipp_d[~valid | (hic_tipo == 0)] = np.nan
        src_d = 'tif'
    elif zonal_dict is not None:
        ipp_d = spatialize_from_csv(zonal_dict, 'IPP_D', hic_tipo, anillo, valid)
        if ipp_d is not None:
            print(f"\n9d. IPP Estrategia D {year} (espacializado por zona desde CSV R)...")
            src_d = 'csv'

    if ipp_d is not None:
        nota_src = "(GeoTIFFetik)" if src_d == 'tif' else "(R CSVren balio zonala pixelez pixel)"
        plot_ipp_spatial(ipp_d, ref_transform,
            f"LPI D estrategia (pagadi-erreferentzia) — Aralar {year}\n"
            f"Larre-anomalia - pagadi-anomalia {nota_src}",
            os.path.join(args.output_dir, f"aralar_mapa_IPP_D_todos_{year}.png"),
            hic_tipo=hic_tipo)

        ipp_d_pasto = ipp_d.copy()
        ipp_d_pasto[hic_tipo != 1] = np.nan

        plot_ipp_spatial(ipp_d_pasto, ref_transform,
            f"LPI D estrategia (pagadi-erreferentzia) — Aralarko larreak soilik {year}\n"
            f"Larre-anomalia - pagadi-anomalia {nota_src}",
            os.path.join(args.output_dir, f"aralar_mapa_IPP_D_pastos_{year}.png"),
            hic_tipo=hic_tipo)
    elif args.ipp_d:
        print(f"\n  AVISO: No se encontro {args.ipp_d} ni se ha podido obtener desde CSV.")

    # ── 10e. CONSENSO + CONCORDANCIA si hay >= 2 estrategias adicionales — NUEVO v5 ──
    ipp_consenso = None
    ipp_consenso_pasto = None
    concordancia_n = None
    estrategias_disp = []
    estrategias_arr = []
    if ipp_a is not None: estrategias_disp.append('A'); estrategias_arr.append(ipp_a)
    if ipp_b is not None: estrategias_disp.append('B'); estrategias_arr.append(ipp_b)
    estrategias_disp.append('C'); estrategias_arr.append(ipp)  # C siempre
    if ipp_d is not None: estrategias_disp.append('D'); estrategias_arr.append(ipp_d)

    n_estr = len(estrategias_arr)
    if n_estr >= 3:
        print(f"\n9e. IPP consenso + concordancia ({n_estr} estrategias: {'+'.join(estrategias_disp)})...")

        stack = np.stack(estrategias_arr, axis=0)  # shape (n_estr, H, W)
        # Consenso = media ignorando NaN
        ipp_consenso = np.nanmean(stack, axis=0).astype(np.float32)
        ipp_consenso[~valid | (hic_tipo == 0)] = np.nan

        # Concordancia: nº de estrategias coincidentes en signo (sobre las no-NaN)
        # 1) signo de cada estrategia (-1, 0, +1); NaN -> 0 contado aparte
        signos = np.sign(stack)
        signos[np.isnan(stack)] = 0  # las NaN no votan
        signo_consenso = np.sign(np.nansum(stack, axis=0))
        # Para cada píxel, contar cuántas estrategias coinciden en signo con el consenso
        coincidentes = np.sum(signos == signo_consenso[None, :, :], axis=0).astype(np.int8)
        # Si una estrategia es NaN no cuenta como coincidente; restar las NaN
        n_nan = np.sum(np.isnan(stack), axis=0).astype(np.int8)
        # n_validas para ese píxel = n_estr - n_nan; concordancia_n = coincidentes - n_nan
        # (porque las NaN entran en coincidentes con signo 0 que puede coincidir con signo_consenso=0)
        n_validas = (n_estr - n_nan).astype(np.int8)
        concordancia_n = (coincidentes - n_nan).astype(np.int8)
        # Para señalizar el píxel inválido
        concordancia_n[~valid | (hic_tipo == 0)] = 0
        n_validas[~valid | (hic_tipo == 0)] = 0

        # Mapa consenso (solo pastos)
        ipp_consenso_pasto = ipp_consenso.copy()
        ipp_consenso_pasto[hic_tipo != 1] = np.nan
        plot_ipp_spatial(ipp_consenso_pasto, ref_transform,
            f"LPI ADOSTASUNA ({n_estr} estrategia) — Aralarko larreak soilik {year}\n"
            f"{'+'.join(estrategias_disp)} batezbestekoa (NaN-ak baztertuta)",
            os.path.join(args.output_dir, f"aralar_mapa_IPP_consenso_pastos_{year}.png"),
            hic_tipo=hic_tipo)

        # Mapa concordancia (en pastos)
        concord_pasto = concordancia_n.copy().astype(np.float32)
        concord_pasto[hic_tipo != 1] = np.nan
        # Plot clasificado: 1..n_estr
        labels = [f"{k}/{n_estr} estrategia" for k in range(1, n_estr+1)]
        # Paleta gradiente del rojo (poca concordancia) al verde (4/4)
        if n_estr == 4:
            colors = ['#d32f2f', '#ff9800', '#fdd835', '#1b5e20']
        elif n_estr == 3:
            colors = ['#d32f2f', '#fdd835', '#1b5e20']
        else:
            colors = ['#d32f2f'] * n_estr
        # Reemplazar NaN por 0 antes del cast para evitar RuntimeWarning
        concord_pasto_int = np.where(np.isnan(concord_pasto), 0, concord_pasto).astype(np.int8)
        plot_classified(concord_pasto_int, ref_transform,
            f"Zeinu-bateragarritasuna (n={n_estr} estrategia) — Aralarko larreak {year}\n"
            f"LPIaren zeinuan bat datozen estrategien kopurua",
            labels, colors,
            os.path.join(args.output_dir, f"aralar_mapa_IPP_concordancia_pastos_{year}.png"),
            hic_tipo=hic_tipo, pasto_only=True)

        # Resumen comparativo en consola (solo pastos por anillo)
        print(f"\n  Comparativa entre estrategias por anillo (solo pastos):")
        hdr = f"  {'Anillo':25s}" + "".join([f"  IPP_{e:>1s}" for e in estrategias_disp]) + "  Consen."
        print(hdr)
        for ani, nombre in [(1, 'Cercano <500m'), (2, 'Intermedio 500-1500m'), (3, 'Remoto >1500m')]:
            m = (hic_tipo == 1) & (anillo == ani)
            row = f"  {nombre:25s}"
            for arr in estrategias_arr:
                mv = m & ~np.isnan(arr)
                row += f"  {np.mean(arr[mv]):+6.3f}" if np.sum(mv) > 0 else f"  {'   NA':>6s}"
            mc = m & ~np.isnan(ipp_consenso)
            row += f"  {np.mean(ipp_consenso[mc]):+6.3f}" if np.sum(mc) > 0 else f"  {'   NA':>6s}"
            print(row)

        # Resumen concordancia en pastos
        print(f"\n  Distribucion de concordancia (% pixeles pasto, n={n_estr} estrategias):")
        for k in range(1, n_estr+1):
            pct = 100 * np.sum((concordancia_n == k) & pasto_mask) / max(n_pasto, 1)
            print(f"    {k}/{n_estr} ({'robusto' if k==n_estr else 'parcial' if k>=n_estr-1 else 'sensible'}): {pct:5.1f}%")

    # ── 10f. MAPA MODELO INVERTIDO (pipeline dual) — NUEVO v5 ──
    if modelo_inv is not None:
        print(f"\n9f. Mapa de modelo PROSAIL invertido por pixel...")
        # Reclasificar 1=pasto, 3=hayedo a 1,2 para visualización
        modelo_vis = np.zeros_like(modelo_inv, dtype=np.int8)
        modelo_vis[modelo_inv == 1] = 1
        modelo_vis[modelo_inv == 3] = 2
        plot_classified(modelo_vis, ref_transform,
            f"Pixelez pixel aplikatutako PROSAIL eredua — Aralar {year}\n"
            f"Kanalizazio bikoitza: NN_PASTO HIC TIPO=1 gainean, NN_HAYEDO HIC TIPO=3 gainean",
            ['NN_PASTO (HIC 1)', 'NN_HAYEDO (HIC 3)'],
            ['#4caf50', '#7b1fa2'],
            os.path.join(args.output_dir, f"aralar_mapa_modelo_invertido_{year}.png"),
            hic_tipo=hic_tipo, pasto_only=False)

    # ── 11. EXPORTAR GeoTIFFs ──
    print(f"\n10. Exportando GeoTIFFs {year}...")

    out_f32 = ref_meta.copy()
    out_f32.update(dtype='float32', count=1, nodata=-9999)
    out_u8 = ref_meta.copy()
    out_u8.update(dtype='uint8', count=1, nodata=0)

    exports_f32 = {
        f'aralar_biomasa_pastos_{year}': biomass,
        f'aralar_capacidad_pastos_{year}': capacity,
        f'aralar_balance_pastos_{year}': balance,
        f'aralar_IPP_C_espacializado_{year}': ipp,
        f'aralar_IPP_C_pastos_{year}': ipp_pasto,
        f'aralar_pendiente_{year}': slope,
    }
    if ipp_a is not None:
        exports_f32[f'aralar_IPP_A_todos_{year}'] = ipp_a
        exports_f32[f'aralar_IPP_A_pastos_{year}'] = ipp_a_pasto
    if ipp_b is not None:
        exports_f32[f'aralar_IPP_B_pastos_{year}'] = ipp_b_pasto
    if ipp_d is not None:
        exports_f32[f'aralar_IPP_D_todos_{year}'] = ipp_d
        exports_f32[f'aralar_IPP_D_pastos_{year}'] = ipp_d_pasto
    if ipp_consenso is not None:
        exports_f32[f'aralar_IPP_consenso_pastos_{year}'] = ipp_consenso_pasto
    for name, data in exports_f32.items():
        d = data.copy()
        d[np.isnan(d)] = -9999
        path = os.path.join(args.output_dir, f"{name}.tif")
        with rasterio.open(path, 'w', **out_f32) as dst:
            dst.write(d, 1)
        print(f"  -> {path}")

    exports_u8 = {
        f'aralar_calidad_pastos_{year}': quality,
        f'aralar_riesgo_pastos_{year}': risk,
        f'aralar_anillo_dist_{year}': anillo.astype(np.uint8),
    }
    if concordancia_n is not None:
        exports_u8[f'aralar_IPP_concordancia_{year}'] = concordancia_n.astype(np.uint8)
    if modelo_inv is not None:
        exports_u8[f'aralar_modelo_invertido_{year}'] = modelo_inv.astype(np.uint8)
    for name, data in exports_u8.items():
        path = os.path.join(args.output_dir, f"{name}.tif")
        with rasterio.open(path, 'w', **out_u8) as dst:
            dst.write(data.astype(np.uint8), 1)
        print(f"  -> {path}")

    # ══════════════════════════════════════════════════
    # RESUMEN
    # ══════════════════════════════════════════════════
    print(f"\n{'='*65}")
    print(f"  RESUMEN — Pastos de Aralar {year}")
    print(f"  Carga: {Config.CARGA_REF} UGM/ha | Temporada: {Config.TEMPORADA} d")
    print(f"  Turnover: {Config.TURNOVER} | Aprovech.: {Config.APROVECH}")
    print(f"{'='*65}")

    bp = biomass[pasto_mask]
    cp = capacity[pasto_mask]
    blp = balance[pasto_mask & ~np.isnan(balance)]
    sp = slope[pasto_mask]

    print(f"  Superficie pastoral:     {n_pasto * 400 / 10000:.0f} ha")
    print(f"  Biomasa en pie media:    {np.mean(bp):.0f} kg MS/ha")
    print(f"  Produccion total media:  {np.mean(bp) * Config.TURNOVER:.0f} kg MS/ha")
    print(f"  Capacidad media:         {np.mean(cp):.2f} UGM/ha")
    print(f"  Balance medio:           {np.mean(blp):+.3f} UGM/ha")
    print(f"  Pendiente media:         {np.mean(sp):.1f} grados")
    print(f"\n  % pixeles deficit:       {100*np.sum(blp<0)/len(blp):.1f}%")
    print(f"  % pixeles equilibrio:    {100*np.sum(np.abs(blp)<=0.10)/len(blp):.1f}%")
    print(f"  % pixeles excedente:     {100*np.sum(blp>0.10)/len(blp):.1f}%")

    # IPP por anillo (solo pastos)
    print(f"\n  IPP medio por anillo (solo pastos):")
    for ani, nombre in [(1, 'Cercano <500m'), (2, 'Intermedio 500-1500m'), (3, 'Remoto >1500m')]:
        m = (hic_tipo == 1) & (anillo == ani) & ~np.isnan(ipp)
        if np.sum(m) > 0:
            print(f"    {nombre:25s}: IPP = {np.mean(ipp[m]):+.3f}  (n={np.sum(m):,})")

    # IPP por tipo de hábitat
    print(f"\n  IPP medio por tipo de habitat:")
    for tipo, nombre in [(1, 'Pasto'), (2, 'Brezal'), (3, 'Hayedo'), (4, 'Encinar')]:
        m = (hic_tipo == tipo) & ~np.isnan(ipp)
        if np.sum(m) > 0:
            print(f"    {nombre:10s}: IPP = {np.mean(ipp[m]):+.3f}  (n={np.sum(m):,})")

    # Calidad
    qp = quality[pasto_mask]
    print(f"\n  Calidad pastoral:")
    for qi, ql in enumerate(['Muy pobre', 'Pobre', 'Moderado', 'Bueno', 'Excelente'], 1):
        pct = 100 * np.sum(qp == qi) / len(qp) if len(qp) > 0 else 0
        print(f"    {ql:12s}: {pct:5.1f}%")

    # Riesgo
    rp = risk[pasto_mask]
    print(f"\n  Diagnostico de riesgo:")
    risk_labels = {2: 'Matorralizacion', 3: 'Subpastoreo leve', 4: 'Equilibrio',
                   5: 'Presion elevada', 6: 'Degradacion', 7: 'Erosion'}
    for ri, rl in risk_labels.items():
        pct = 100 * np.sum(rp == ri) / len(rp) if len(rp) > 0 else 0
        if pct > 0:
            print(f"    {rl:20s}: {pct:5.1f}%")

    print(f"\n  Ficheros en: {args.output_dir}/")
    print(f"  Resumen completo en: {log_path}")
    print(f"{'='*65}\n")

    # Cerrar logger
    if isinstance(sys.stdout, TeeLogger):
        sys.stdout.close()
        sys.stdout = sys.stdout.terminal if hasattr(sys.stdout, 'terminal') else sys.__stdout__


if __name__ == '__main__':
    main()
