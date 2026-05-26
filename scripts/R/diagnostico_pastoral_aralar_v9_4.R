# ============================================================================
# DIAGNÓSTICO PASTORAL INTEGRADO v9.4 — Sierra de Aralar
# Pipeline PROSAIL-NN + ERA5-Land + SPEI → IPP cuádruple (A/B/C/D)
# ============================================================================
#
# ENTRADA: panel CSV exportado desde GEE (pipeline_integrado_aralar.js)
#          + SPEI calculado en R (Aralar_SPEI_SerieLarga.R)
#
# NOVEDADES v9.4 (sobre v9.3):
#   - G6m NUEVO: candidato GAM con interaccion tensorial te(GDA, P60).
#       Captura la dependencia del efecto del calor acumulado sobre el LAI
#       segun la disponibilidad hidrica previa. Justificacion empirica en
#       la comparativa de modelos estacionales (abril 2026): G6m mejora
#       sobre G6 en dAIC ~24 unidades sin absorber senal pastoral.
#   - Chequeo de hull convexo en exportacion a GEE: si hay fechas con
#       (GDA, P60, SM_root) fuera del rango de entrenamiento del modelo
#       seleccionado, se marcan como NA y se filtran del CSV. Relevante
#       para candidatos con termino tensorial te().
#   - Chequeo anti-absorcion (salvaguarda): si un futuro candidato con
#       termino estacional s(doy) es seleccionado, se verifica que el
#       diferencial cercano-remoto en agosto se preserva. Si se colapsa
#       (<0.08 LAI_norm), se emite aviso sobre degradacion silenciosa
#       de la capacidad diagnostica del IPP.
#
# NOVEDADES v9 (sobre v8):
#   - Sección 3.5: Hayedo como termómetro climático sin pastoreo.
#       Genera LAI_hayedo_anom (anomalía respecto a fenología media multianual).
#   - Fix regla ΔAIC<2: parsimonia ahora se respeta (en v8 se sobreescribía).
#   - Sección 7b: Estrategia D = LAI_pasto_norm − LAI_hayedo_anom.
#       Cuarta estrategia de diagnóstico independiente de los modelos A/B/C.
#   - IPP_consenso = media(A, B, C, D) ignorando NAs.
#   - Métrica concordancia_n (formato "k/n" de coincidencias de signo).
#   - Gráfico 9h: IPP_triangulacion_ABCD.png.
#
# SALIDAS:
#   - comparacion_modelos.csv         (candidatos GAM/LM + AIC/R²adj)
#   - diagnostico_pastoral.csv        (IPP por zona/año/estrategia + consenso)
#   - hayedo_serie_referencia.csv     (serie temporal hayedo + anomalías)
#   - gam_predicciones_fecha.csv      (para IPP espacializado en GEE)
#   - validacion_loyo.csv             (leave-one-year-out)
#   - *.png                           (gráficos de diagnóstico, 9 figuras)
#
# ============================================================================

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 0. CONFIGURACIÓN
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# ── Rutas de entrada ──
LAI_CSV   <- "LAI_zonal_HIC_aralar.csv"        # ← Export 1 de GEE
CLIMA_CSV <- "clima_diario_aralar.csv"        # ← Export 2 de GEE
SPEI_CSV  <- "Aralar_SPEI_all_1981_present.csv"  # ← SPEI calculado en R
ZONAS_CSV <- "zonas_caracterizacion.csv"      # ← Caracterización de zonas (alt, pend, dist)

# ── Zona pseudo-referencia (Estrategia B) ──
# pasto_remoto: >1500m de bordas, mínima presión ganadera
ZONA_REFERENCIA <- "pasto_remoto"

# ── Zonas diana para diagnóstico pastoral ──
# Solo los pastos. Brezales y hayedos se analizan por separado.
ZONAS_PASTO  <- c("pasto_cercano", "pasto_intermedio", "pasto_remoto")
ZONAS_BREZAL <- c("brezal_cercano", "brezal_intermedio", "brezal_remoto")
ZONAS_HAYEDO <- c("hayedo_cercano", "hayedo_intermedio", "hayedo_remoto")

# ── Umbrales IPP ──
IPP_THRESHOLDS <- c(
  severo_neg = -0.20,
  moderado_neg = -0.10,
  leve_neg = -0.05,
  leve_pos = 0.05,
  moderado_pos = 0.10,
  severo_pos = 0.20
)

# ── Tolerancia para gráficos (banda de equilibrio) ──
TOLERANCIA <- 0.15  # ±15% respecto a LAI esperado

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PAQUETES
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

pkgs <- c("dplyr", "tidyr", "ggplot2", "lubridate", "patchwork", "zoo", "mgcv")
for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. LECTURA Y PREPARACIÓN
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("=" , rep("=", 60), "\n")
cat("  LARRE-DIAGNOSTIKOA — Aralarko mendilerroa\n")
cat("  PROSAIL-NN + ERA5-Land + SPEI\n")
cat("=", rep("=", 60), "\n\n")

# ── LAI zonal (Export 1 de GEE) ──
cat("▸ LAI zonala irakurtzen:", LAI_CSV, "\n")
lai_wide <- read.csv(LAI_CSV, stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

cat("  Eszenak:", nrow(lai_wide), "\n")

# Nombres de zona (z11..z43 = TIPO*10 + anillo)
zona_map <- c(
  "LAI_z11" = "pasto_cercano",    "LAI_z12" = "pasto_intermedio",   "LAI_z13" = "pasto_remoto",
  "LAI_z21" = "brezal_cercano",   "LAI_z22" = "brezal_intermedio",  "LAI_z23" = "brezal_remoto",
  "LAI_z31" = "hayedo_cercano",   "LAI_z32" = "hayedo_intermedio",  "LAI_z33" = "hayedo_remoto",
  "LAI_z41" = "encinar_cercano",  "LAI_z42" = "encinar_intermedio", "LAI_z43" = "encinar_remoto"
)

# Pivotar a formato largo (zona × fecha)
# Primero LAI
lai_long <- lai_wide %>%
  pivot_longer(cols = starts_with("LAI_z"),
               names_to = "zona_col", values_to = "LAI") %>%
  mutate(zona = zona_map[zona_col],
         zona_id = gsub("LAI_z", "", zona_col)) %>%
  filter(!is.na(LAI), LAI > 0, LAI != -9999)

# Pivotar conteos (N_z11, N_z12, ...)
n_cols <- grep("^N_z", names(lai_wide), value = TRUE)
if (length(n_cols) > 0) {
  n_long <- lai_wide %>%
    select(date, all_of(n_cols)) %>%
    pivot_longer(cols = starts_with("N_z"),
                 names_to = "n_col", values_to = "n_pixels") %>%
    mutate(zona_id = gsub("N_z", "", n_col)) %>%
    select(date, zona_id, n_pixels)
  
  lai_long <- lai_long %>%
    left_join(n_long, by = c("date", "zona_id"))
} else {
  lai_long$n_pixels <- NA
}

# ── Leer caracterización de zonas para obtener N total por zona ──
n_total_zona <- NULL
if (file.exists(ZONAS_CSV)) {
  zonas_char <- read.csv(ZONAS_CSV, stringsAsFactors = FALSE)
  n_total_zona <- zonas_char %>% select(zona_id, zona_name, n_pixels_total = n_pixels)
  cat("\n  Zonen karakterizazioa:\n")
  print(as.data.frame(zonas_char))
}

# ── FILTRO DE COBERTURA NUBOSA ──
# Descartar observaciones donde <50% de los píxeles de la zona son válidos
# Esto elimina escenas con nubes parciales sobre la zona
MIN_COVERAGE_PCT <- 50  # % mínimo de píxeles válidos

if (!is.null(n_total_zona) && !all(is.na(lai_long$n_pixels))) {
  lai_long <- lai_long %>%
    left_join(n_total_zona %>% mutate(zona_id = as.character(zona_id)),
              by = "zona_id") %>%
    mutate(
      coverage_pct = ifelse(n_pixels_total > 0,
                            100 * n_pixels / n_pixels_total, NA)
    )
  
  n_before <- nrow(lai_long)
  lai_long <- lai_long %>%
    filter(is.na(coverage_pct) | coverage_pct >= MIN_COVERAGE_PCT)
  n_after <- nrow(lai_long)
  
  cat("\n  Estaldura iragazkia (>=", MIN_COVERAGE_PCT, "%):",
      n_before, "->", n_after, "erregistro\n")
  if (n_before > n_after) {
    cat("  Ezabatuta:", n_before - n_after, "erregistro hodei partzialekin\n")
  }
} else {
  cat("\n  Pixel-zenbaketarik gabe — estaldura iragazkia ez da aplikatu.\n")
  cat("  (GEEko v3 pipeline-a exekutatu N_z* CSVan eskuratzeko)\n")
}

lai_long <- lai_long %>%
  select(date, year, doy, zona, LAI) %>%
  arrange(zona, date)

cat("  LAI erregistroak (formatu luzea):", nrow(lai_long), "\n")
cat("  Zonak:", paste(sort(unique(lai_long$zona)), collapse = ", "), "\n")

# ── Clima diario (Export 2 de GEE) ──
cat("\n▸ Klima eguneroa irakurtzen:", CLIMA_CSV, "\n")
clima <- read.csv(CLIMA_CSV, stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

cat("  Egunak:", nrow(clima), "\n")

# Calcular variables acumuladas
has_rad_col <- "Rad_MJ" %in% names(clima)
if (!has_rad_col) {
  cat("  OHARRA: Rad_MJ ez da aurkitu CSV-an — GEE pipeline v3 berriz exekutatu\n")
  clima$Rad_MJ <- NA
}

clima <- clima %>%
  group_by(year) %>%
  mutate(
    P30 = zoo::rollsumr(P_mm, k = 30, fill = NA),
    P60 = zoo::rollsumr(P_mm, k = 60, fill = NA),
    GD_dia = pmax(T2m_C - 5, 0),
    GDA = cumsum(ifelse(doy >= 91, GD_dia, 0)),
    Rad7 = if (has_rad_col) zoo::rollmeanr(Rad_MJ, k = 7, fill = NA) else NA_real_,
    Rad30 = if (has_rad_col) zoo::rollsumr(Rad_MJ, k = 30, fill = NA) else NA_real_
  ) %>%
  ungroup() %>%
  select(date, P30, P60, T2m_C, GDA, SM_sup, SM_root, Rad7, Rad30)

# ── Join LAI + Clima por fecha ──
cat("\n▸ LAI + klima dataren arabera elkartzen...\n")
panel <- lai_long %>%
  left_join(clima, by = "date")

cat("  Panel erregistroak:", nrow(panel), "\n")
cat("  P60rekin:", sum(!is.na(panel$P60)), "\n")

# Verificar
cat("  Datak zonaka:\n")
panel %>% count(zona, year) %>%
  pivot_wider(names_from = year, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print()

# ── SPEI (join por mes) ──
spei_ok <- FALSE
if (file.exists(SPEI_CSV)) {
  cat("\n▸ SPEI irakurtzen:", SPEI_CSV, "\n")
  spei <- read.csv(SPEI_CSV, stringsAsFactors = FALSE) %>%
    mutate(date = as.Date(date),
           ym = floor_date(date, "month"))

  # Seleccionar columnas SPEI disponibles
  spei_cols <- grep("^SPEI_", names(spei), value = TRUE)
  cat("  Eskala eskuragarriak:", paste(spei_cols, collapse = ", "), "\n")

  if (length(spei_cols) > 0) {
    spei_join <- spei %>% select(ym, all_of(spei_cols))

    panel <- panel %>%
      mutate(ym = floor_date(date, "month")) %>%
      left_join(spei_join, by = "ym") %>%
      select(-ym)

    spei_ok <- TRUE
    cat("  SPEI panelarekin hilabeteka elkartuta.\n")
    cat("  SPEI_3 NA ez diren balioak:", sum(!is.na(panel$SPEI_3)), "/", nrow(panel), "\n")
  }
} else {
  cat("\n⚠ SPEI ez da aurkitu (", SPEI_CSV, "). SPEIdun ereduak baztertuko dira.\n")
}

# ── Colores por zona ──
zonas <- sort(unique(panel$zona))
n_zonas <- length(zonas)
pal_pasto  <- c("pasto_cercano" = "#d73027", "pasto_intermedio" = "#fee08b", "pasto_remoto" = "#1a9850")
pal_brezal <- c("brezal_cercano" = "#e31a1c", "brezal_intermedio" = "#fb9a99", "brezal_remoto" = "#fdbf6f")
pal_hayedo <- c("hayedo_cercano" = "#1f78b4", "hayedo_intermedio" = "#a6cee3", "hayedo_remoto" = "#b2df8a")
pal_encinar <- c("encinar_cercano" = "#ff7f00", "encinar_intermedio" = "#ffc966", "encinar_remoto" = "#cab2d6")
colores <- c(pal_pasto, pal_brezal, pal_hayedo, pal_encinar)

# ── Etiquetas de zonas en euskera (solo display; la columna `zona` se mantiene
#    intacta para preservar compatibilidad con CSV/index.html/GEE) ──
nombres_zona <- c(
  "pasto_cercano"     = "Larrea (hurbila)",   "pasto_intermedio"   = "Larrea (ertaina)",   "pasto_remoto"     = "Larrea (urruna)",
  "brezal_cercano"    = "Txilardia (hurbila)", "brezal_intermedio" = "Txilardia (ertaina)", "brezal_remoto"    = "Txilardia (urruna)",
  "hayedo_cercano"    = "Pagadia (hurbila)",   "hayedo_intermedio"= "Pagadia (ertaina)",    "hayedo_remoto"   = "Pagadia (urruna)",
  "encinar_cercano"   = "Artadia (hurbila)",   "encinar_intermedio"="Artadia (ertaina)",    "encinar_remoto"  = "Artadia (urruna)"
)

# ── Zonas pastoreadas (para IPP) ──
zonas_past <- intersect(zonas, ZONAS_PASTO)
zonas_past <- setdiff(zonas_past, ZONA_REFERENCIA)
cat("\n  Erreferentzia-zona (B estr.):", ZONA_REFERENCIA, "\n")
cat("  Larre-zonak LPI kalkulatzeko:", paste(zonas_past, collapse = ", "), "\n")
cat("  Txilardi-zonak:", paste(intersect(zonas, ZONAS_BREZAL), collapse = ", "), "\n")
cat("  Pagadi-zonak:", paste(intersect(zonas, ZONAS_HAYEDO), collapse = ", "), "\n")


# zonas_char ya leído en sección 2 (filtro cobertura)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2b. NORMALIZACIÓN DEL LAI POR ZONA-AÑO
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Elimina diferencias estructurales entre tipos de cubierta (hayedo LAI~4
# vs pasto LAI~1.5). Cada zona se escala a [0,1] por su máximo anual.
# El modelo climático trabaja sobre LAI_norm → residuos comparables.

cat("\n▸ LAI zona-urteka normalizatzen...\n")

panel <- panel %>%
  group_by(zona, year) %>%
  mutate(
    LAI_max_zy = max(LAI, na.rm = TRUE),
    LAI_norm = ifelse(LAI_max_zy > 0, LAI / LAI_max_zy, NA)
  ) %>%
  ungroup()

cat("  LAI_norm tartea: [", round(min(panel$LAI_norm, na.rm = TRUE), 3), ",",
    round(max(panel$LAI_norm, na.rm = TRUE), 3), "]\n")

# Resumen de normalización
panel %>%
  group_by(zona) %>%
  summarise(
    LAI_abs_mean = round(mean(LAI, na.rm = TRUE), 2),
    LAI_max_mean = round(mean(LAI_max_zy, na.rm = TRUE), 2),
    LAI_norm_mean = round(mean(LAI_norm, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>% as.data.frame() %>% print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. NIVEL 1 — ANÁLISIS EXPLORATORIO
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  1. MAILA — Esplorazioa\n")
cat(rep("=", 60), "\n")

# Estadísticas por zona y año
resumen <- panel %>%
  group_by(zona, year) %>%
  summarise(
    n = n(),
    LAI_mean = round(mean(LAI, na.rm = TRUE), 2),
    LAI_sd   = round(sd(LAI, na.rm = TRUE), 2),
    LAI_max  = round(max(LAI, na.rm = TRUE), 2),
    P60_mean = round(mean(P60, na.rm = TRUE), 1),
    GDA_max  = round(max(GDA, na.rm = TRUE), 0),
    .groups = "drop"
  )
print(as.data.frame(resumen))

# ── Gráfico: Series temporales de LAI por zona, facetado por año ──
p_series <- ggplot(panel, aes(x = doy, y = LAI, color = zona)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = colores, labels = nombres_zona) +
  facet_wrap(~year, scales = "free_x") +
  labs(title = "LAI PROSAIL-NN zonaka — Aralar",
       subtitle = "Urteka faketatua | Iturria: Sentinel-2 NN-rekin alderantzikatua",
       x = "Urteko eguna", y = "LAI (m\u00b2/m\u00b2)", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 13))

ggsave("N1_series_LAI.png", p_series, width = 14, height = 6, dpi = 150)
print(p_series)

# ── Correlaciones LAI ~ clima ──
clim_vars <- c("P30", "P60", "T2m_C", "GDA", "SM_sup", "SM_root", "Rad7", "Rad30")
if (spei_ok) clim_vars <- c(clim_vars, "SPEI_3")
clim_present <- clim_vars[clim_vars %in% names(panel)]

cat("\n━━━ LAI ~ klima korrelazioak ━━━\n")
cor_mat <- panel %>%
  select(LAI, all_of(clim_present)) %>%
  cor(use = "pairwise.complete.obs")
print(round(cor_mat, 3))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3.5. HAYEDO COMO TERMÓMETRO CLIMÁTICO (referencia sin pastoreo)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Nivel 2 de explotación del modelo dual PROSAIL: el hayedo (HIC 9120/9150)
# recibe exactamente el mismo clima que los pastos de Aralar pero NO está
# sometido a presión ganadera. Por tanto, su LAI normalizado es una sonda
# del forzamiento climático "limpio" que podemos usar para construir un IPP
# diferencial independiente de los modelos estadísticos (A, B, C).
#
# Estrategia:
#   1. Serie temporal LAI_hayedo(t) = media ponderada por N de las tres
#      zonas de hayedo (cercano + intermedio + remoto) en cada fecha.
#   2. LAI_hayedo_norm(t) = escalado por LAI máximo anual del hayedo.
#   3. LAI_hayedo_anom(t,y) = LAI_hayedo_norm(t) menos la fenología media
#      multi-anual del hayedo (curva suavizada por DOY). Queda una serie
#      centrada en 0 que representa la anomalía climática de cada fecha.
#   4. Estrategia D: IPP_D = (LAI_pasto_norm_res respecto al hayedo) / LAI
#      donde LAI_pasto_norm_res = LAI_pasto_norm − LAI_hayedo_anom.
#      La lógica: si el hayedo baja 0.1 por sequía, el pasto "debería" bajar
#      lo mismo; la parte que NO se explica por la anomalía del hayedo es
#      la atribuible a presión ganadera.

cat("\n\n", rep("=", 60), "\n")
cat("  3.5. PAGADIA larretu gabeko klima-erreferentzia gisa\n")
cat(rep("=", 60), "\n")

hayedo_zonas_pres <- intersect(zonas, ZONAS_HAYEDO)
cat("\nPagadi-zona eskuragarriak:", paste(hayedo_zonas_pres, collapse = ", "), "\n")

if (length(hayedo_zonas_pres) >= 1) {
  # 1. Serie temporal del hayedo — media simple por fecha de las 3 zonas
  #    (El panel no lleva n_pixels; las tres zonas de hayedo tienen
  #     áreas comparables, así que el sesgo por no ponderar es marginal.)
  hayedo_serie <- panel %>%
    filter(zona %in% hayedo_zonas_pres) %>%
    group_by(date) %>%
    summarise(
      LAI_hayedo = mean(LAI, na.rm = TRUE),
      n_zonas_hayedo = sum(!is.na(LAI)),
      year = first(year),
      doy  = first(doy),
      .groups = "drop"
    ) %>%
    filter(!is.na(LAI_hayedo), is.finite(LAI_hayedo))

  cat("Pagadi seriea:", nrow(hayedo_serie), "data\n")
  cat("LAI pagadi tartea:", round(min(hayedo_serie$LAI_hayedo), 2),
      "—", round(max(hayedo_serie$LAI_hayedo), 2), "\n")

  # 2. Normalización: escalar por el máximo anual del hayedo
  hayedo_serie <- hayedo_serie %>%
    group_by(year) %>%
    mutate(LAI_hayedo_max  = max(LAI_hayedo, na.rm = TRUE),
           LAI_hayedo_norm = LAI_hayedo / LAI_hayedo_max) %>%
    ungroup()

  # 3. Fenología media multi-anual del hayedo vs DOY (curva de referencia)
  #    Se ajusta un GAM suave para tener una curva estable de la fenología
  #    típica, y la anomalía es la desviación de cada fecha respecto a ella.
  hayedo_pheno_ok <- FALSE
  if (nrow(hayedo_serie) >= 15) {
    tryCatch({
      pheno_mod <- gam(LAI_hayedo_norm ~ s(doy, k = 6), data = hayedo_serie)
      hayedo_serie <- hayedo_serie %>%
        mutate(
          LAI_hayedo_norm_pheno = predict(pheno_mod, newdata = .),
          LAI_hayedo_anom = LAI_hayedo_norm - LAI_hayedo_norm_pheno
        )
      hayedo_pheno_ok <- TRUE
      cat("Pagadi fenologia ertaina doitua: R2 =",
          round(summary(pheno_mod)$r.sq, 3),
          " | dev.expl =", round(summary(pheno_mod)$dev.expl * 100, 1), "%\n")
      cat("Pagadi anomalia tartea:",
          round(min(hayedo_serie$LAI_hayedo_anom, na.rm = TRUE), 3),
          "—",
          round(max(hayedo_serie$LAI_hayedo_anom, na.rm = TRUE), 3), "\n")
    }, error = function(e) {
      cat("  ! Ezin izan da pagadi fenologia doitu:", e$message, "\n")
    })
  }

  if (!hayedo_pheno_ok) {
    cat("  (Fenologia ertaina ez dago eskuragarri; D estrategia baztertzen da)\n")
    hayedo_serie$LAI_hayedo_anom <- NA_real_
    hayedo_serie$LAI_hayedo_norm_pheno <- NA_real_
  }

  # 4. Incorporar al panel la anomalía y la serie normalizada del hayedo
  #    como variables a nivel de fecha (no dependen de zona).
  panel <- panel %>%
    left_join(
      hayedo_serie %>%
        select(date, LAI_hayedo, LAI_hayedo_norm, LAI_hayedo_anom),
      by = "date"
    )

  # Export de la serie hayedo para inspección
  write.csv(hayedo_serie, "hayedo_serie_referencia.csv", row.names = FALSE)
  cat("\n✓ Esportatua: hayedo_serie_referencia.csv\n")

  hayedo_disponible <- TRUE

} else {
  cat("⚠ Ez dago pagadi-zonarik eskuragarri. D estrategia baztertuko da.\n")
  panel$LAI_hayedo       <- NA_real_
  panel$LAI_hayedo_norm  <- NA_real_
  panel$LAI_hayedo_anom  <- NA_real_
  hayedo_disponible <- FALSE
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. NIVEL 2 — MODELOS CLIMÁTICOS: LM + GAM
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  2. MAILA — Klima-ereduak (LM + GAM)\n")
cat(rep("=", 60), "\n")

lai_aoi <- panel %>%
  group_by(date, year, doy) %>%
  summarise(
    LAI_aoi = mean(LAI_norm, na.rm = TRUE),
    LAI_aoi_abs = mean(LAI, na.rm = TRUE),
    across(all_of(clim_present), ~mean(., na.rm = TRUE)),
    .groups = "drop"
  )

lai_aoi <- lai_aoi %>% mutate(GDA2 = GDA^2)
panel <- panel %>% mutate(GDA2 = GDA^2)

cat("AOI-batezbesteko behaketak:", nrow(lai_aoi), "\n")

# ── Función auxiliar para extraer métricas de lm o gam ──
get_metrics <- function(m) {
  if (inherits(m, "gam")) {
    s <- summary(m)
    return(list(
      R2     = round(s$r.sq, 4),
      R2_adj = round(s$r.sq, 4),  # GAM reporta R² ajustado como r.sq
      dev_expl = round(s$dev.expl, 4),
      AIC    = round(AIC(m), 1),
      n      = nrow(m$model),
      tipo   = "GAM"
    ))
  } else {
    s <- summary(m)
    return(list(
      R2     = round(s$r.squared, 4),
      R2_adj = round(s$adj.r.squared, 4),
      dev_expl = round(s$r.squared, 4),
      AIC    = round(AIC(m), 1),
      n      = nrow(m$model),
      tipo   = "LM"
    ))
  }
}

# ── Modelos candidatos: LM lineales + GAM con splines ──
models <- list()
labels <- list()

# Datos para modelos con SPEI
has_spei <- spei_ok && "SPEI_3" %in% names(lai_aoi) && sum(!is.na(lai_aoi$SPEI_3)) >= 10
has_rad  <- "Rad7" %in% names(lai_aoi) && sum(!is.na(lai_aoi$Rad7)) >= 10

d_full <- lai_aoi %>% filter(!is.na(P60), !is.na(GDA))
if (has_spei) d_spei <- d_full %>% filter(!is.na(SPEI_3)) else d_spei <- d_full
if (has_rad)  d_rad  <- d_full %>% filter(!is.na(Rad7))  else d_rad  <- d_full

# --- LINEALES ---
cat("\n  Eredu linealak doitzen...\n")

tryCatch({
  models[["L1"]] <- lm(LAI_aoi ~ P60 + GDA, data = d_full)
  labels[["L1"]] <- "LM: P60 + GDA"
}, error = function(e) cat("  L1 huts egin du\n"))

tryCatch({
  models[["L2"]] <- lm(LAI_aoi ~ P60 + GDA + GDA2, data = d_full)
  labels[["L2"]] <- "LM: P60 + GDA + GDA2"
}, error = function(e) cat("  L2 huts egin du\n"))

if (has_spei) tryCatch({
  models[["L3"]] <- lm(LAI_aoi ~ P60 + GDA + SPEI_3, data = d_spei)
  labels[["L3"]] <- "LM: P60 + GDA + SPEI_3"
}, error = function(e) cat("  L3 huts egin du\n"))

if (has_rad) tryCatch({
  models[["L4"]] <- lm(LAI_aoi ~ P60 + GDA + Rad7, data = d_rad)
  labels[["L4"]] <- "LM: P60 + GDA + Rad7"
}, error = function(e) cat("  L4 huts egin du\n"))

if (has_spei && has_rad) tryCatch({
  models[["L5"]] <- lm(LAI_aoi ~ P60 + GDA + SPEI_3 + Rad7, data = d_spei %>% filter(!is.na(Rad7)))
  labels[["L5"]] <- "LM: P60 + GDA + SPEI_3 + Rad7"
}, error = function(e) cat("  L5 huts egin du\n"))

# --- GAM con splines suavizados ---
cat("  GAM ereduak doitzen...\n")

# G1: s(GDA) + s(P60) — fenología + precipitación no lineales
tryCatch({
  models[["G1"]] <- gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5), data = d_full)
  labels[["G1"]] <- "GAM: s(GDA) + s(P60)"
}, error = function(e) cat("  G1 huts egin du:", e$message, "\n"))

# G2: s(GDA) + s(P60) + SPEI_3
if (has_spei) tryCatch({
  models[["G2"]] <- gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5) + SPEI_3, data = d_spei)
  labels[["G2"]] <- "GAM: s(GDA) + s(P60) + SPEI_3"
}, error = function(e) cat("  G2 huts egin du\n"))

# G3: s(GDA) + s(P60) + s(Rad7)
if (has_rad) tryCatch({
  models[["G3"]] <- gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5) + s(Rad7, k = 5),
                        data = d_rad)
  labels[["G3"]] <- "GAM: s(GDA) + s(P60) + s(Rad7)"
}, error = function(e) cat("  G3 huts egin du\n"))

# G4: s(GDA) + s(P60) + SPEI_3 + s(Rad7) — modelo completo
if (has_spei && has_rad) tryCatch({
  d_all <- d_spei %>% filter(!is.na(Rad7))
  models[["G4"]] <- gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5) + SPEI_3 + s(Rad7, k = 5),
                        data = d_all)
  labels[["G4"]] <- "GAM: s(GDA) + s(P60) + SPEI_3 + s(Rad7)"
}, error = function(e) cat("  G4 huts egin du\n"))

# G5: s(GDA) + P60 + SPEI_3 + Rad7 — solo GDA no lineal
if (has_spei && has_rad) tryCatch({
  d_all <- d_spei %>% filter(!is.na(Rad7))
  models[["G5"]] <- gam(LAI_aoi ~ s(GDA, k = 5) + P60 + SPEI_3 + Rad7,
                        data = d_all)
  labels[["G5"]] <- "GAM: s(GDA) + P60 + SPEI_3 + Rad7"
}, error = function(e) cat("  G5 huts egin du\n"))

# G6: s(GDA) + s(P60) + s(SM_root)
if ("SM_root" %in% names(d_full) && sum(!is.na(d_full$SM_root)) >= 10) tryCatch({
  models[["G6"]] <- gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5) + s(SM_root, k = 5),
                        data = d_full %>% filter(!is.na(SM_root)))
  labels[["G6"]] <- "GAM: s(GDA) + s(P60) + s(SM_root)"
}, error = function(e) cat("  G6 huts egin du\n"))

# G6m: s(GDA) + s(P60) + s(SM_root) + te(GDA, P60)
# Interaccion tensorial entre calor acumulado y precipitacion reciente.
# Captura la dependencia del efecto del estres termico sobre el LAI segun
# la disponibilidad hidrica previa. Mejora sobre G6 documentada en la
# comparativa de modelos estacionales (abril 2026).
# k = c(4, 4) produce 16 bases de tensor; si el modelo no converge con el
# panel actual, se reintenta con k = c(3, 3).
if ("SM_root" %in% names(d_full) && sum(!is.na(d_full$SM_root)) >= 10) {
  d_sm <- d_full %>% filter(!is.na(SM_root))
  mod_g6m <- tryCatch(
    gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5) + s(SM_root, k = 5) +
                  te(GDA, P60, k = c(4, 4)),
        data = d_sm, method = "REML"),
    error = function(e) {
      cat("  G6m k=c(4,4) huts egin du, k=c(3,3)rekin saiatzen\n")
      tryCatch(
        gam(LAI_aoi ~ s(GDA, k = 5) + s(P60, k = 5) + s(SM_root, k = 5) +
                      te(GDA, P60, k = c(3, 3)),
            data = d_sm, method = "REML"),
        error = function(e2) { cat("  G6m behin betiko huts egin du\n"); NULL }
      )
    }
  )
  if (!is.null(mod_g6m)) {
    models[["G6m"]] <- mod_g6m
    labels[["G6m"]] <- "GAM: s(GDA) + s(P60) + s(SM_root) + te(GDA, P60)"
  }
}

# ── Comparar todos los modelos ──
cat("\n  Doitutako ereduak:", length(models), "\n")

if (length(models) > 0) {
  metrics_list <- lapply(models, get_metrics)
  
  # ── Función para extraer p-valor global del modelo ──
  get_pvalue <- function(m) {
    if (inherits(m, "gam")) {
      # Para GAM: p-valor de los smooth terms (el mayor = menos significativo)
      s <- summary(m)
      p_smooth <- if (nrow(s$s.table) > 0) max(s$s.table[, "p-value"]) else 1
      p_param  <- if (nrow(s$p.table) > 1) max(s$p.table[-1, "Pr(>|t|)"]) else 1
      return(min(p_smooth, p_param))
    } else {
      # Para LM: p-valor del F-test global
      f <- summary(m)$fstatistic
      if (is.null(f)) return(NA)
      return(pf(f[1], f[2], f[3], lower.tail = FALSE))
    }
  }
  
  model_comp <- data.frame(
    id       = names(models),
    formula  = unlist(labels[names(models)]),
    tipo     = sapply(metrics_list, function(m) m$tipo),
    R2       = sapply(metrics_list, function(m) m$R2),
    R2_adj   = sapply(metrics_list, function(m) m$R2_adj),
    dev_expl = sapply(metrics_list, function(m) m$dev_expl),
    AIC      = sapply(metrics_list, function(m) m$AIC),
    p_global = sapply(models, function(m) signif(get_pvalue(m), 4)),
    n        = sapply(metrics_list, function(m) m$n),
    row.names = NULL
  ) %>% arrange(AIC)
  
  cat("\n", rep("=", 85), "\n")
  cat("  EREDUEN KONPARAZIOA (LM + GAM)\n")
  cat(rep("=", 85), "\n\n")
  print(as.data.frame(model_comp))
  write.csv(model_comp, "comparacion_modelos.csv", row.names = FALSE)

  # ── Significancia detallada de TODOS los modelos ──
  cat("\n", rep("-", 85), "\n")
  cat("  EREDUEN ESANGURATASUN XEHATUA\n")
  cat(rep("-", 85), "\n")
  
  signif_rows <- list()
  
  for (mid in model_comp$id) {
    m <- models[[mid]]
    cat("\n  ── ", labels[[mid]], " ──\n")
    
    if (inherits(m, "gam")) {
      s <- summary(m)
      
      # Términos paramétricos (incluyendo intercepto)
      if (nrow(s$p.table) > 0) {
        cat("  Termino parametrikoak:\n")
        for (j in 1:nrow(s$p.table)) {
          pval <- s$p.table[j, "Pr(>|t|)"]
          sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
                 ifelse(pval < 0.05, "*", ifelse(pval < 0.1, ".", "ns"))))
          cat(sprintf("    %-15s  est=%+.4f  se=%.4f  t=%.3f  p=%.2e  %s\n",
              rownames(s$p.table)[j], s$p.table[j, "Estimate"],
              s$p.table[j, "Std. Error"], s$p.table[j, "t value"], pval, sig))
          
          signif_rows[[length(signif_rows) + 1]] <- data.frame(
            modelo = mid, termino = rownames(s$p.table)[j], tipo_termino = "parametrico",
            edf = NA, estimacion = round(s$p.table[j, "Estimate"], 5),
            F_chi = round(s$p.table[j, "t value"], 3),
            p_valor = signif(pval, 4), significancia = sig,
            stringsAsFactors = FALSE)
        }
      }
      
      # Smooth terms
      if (nrow(s$s.table) > 0) {
        cat("  Termino leunduak (splines):\n")
        for (j in 1:nrow(s$s.table)) {
          pval <- s$s.table[j, "p-value"]
          sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
                 ifelse(pval < 0.05, "*", ifelse(pval < 0.1, ".", "ns"))))
          cat(sprintf("    %-15s  edf=%.2f  ref.df=%.2f  F=%.3f  p=%.2e  %s\n",
              rownames(s$s.table)[j], s$s.table[j, "edf"],
              s$s.table[j, "Ref.df"], s$s.table[j, "F"], pval, sig))
          
          signif_rows[[length(signif_rows) + 1]] <- data.frame(
            modelo = mid, termino = rownames(s$s.table)[j], tipo_termino = "smooth",
            edf = round(s$s.table[j, "edf"], 2), estimacion = NA,
            F_chi = round(s$s.table[j, "F"], 3),
            p_valor = signif(pval, 4), significancia = sig,
            stringsAsFactors = FALSE)
        }
      }
      
      cat(sprintf("  R2(adj)=%.4f  Dev.expl=%.1f%%  n=%d\n",
          s$r.sq, s$dev.expl * 100, nrow(m$model)))
      
    } else {
      # LM
      s <- summary(m)
      ct <- s$coefficients
      
      for (j in 1:nrow(ct)) {
        pval <- ct[j, "Pr(>|t|)"]
        sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
               ifelse(pval < 0.05, "*", ifelse(pval < 0.1, ".", "ns"))))
        cat(sprintf("    %-15s  est=%+.5f  se=%.5f  t=%.3f  p=%.2e  %s\n",
            rownames(ct)[j], ct[j, "Estimate"], ct[j, "Std. Error"],
            ct[j, "t value"], pval, sig))
        
        signif_rows[[length(signif_rows) + 1]] <- data.frame(
          modelo = mid, termino = rownames(ct)[j], tipo_termino = "parametrico",
          edf = NA, estimacion = round(ct[j, "Estimate"], 5),
          F_chi = round(ct[j, "t value"], 3),
          p_valor = signif(pval, 4), significancia = sig,
          stringsAsFactors = FALSE)
      }
      
      f <- s$fstatistic
      p_global <- pf(f[1], f[2], f[3], lower.tail = FALSE)
      cat(sprintf("  R2adj=%.4f  F(%d,%d)=%.2f  p_global=%.2e  n=%d\n",
          s$adj.r.squared, f[2], f[3], f[1], p_global, nrow(m$model)))
    }
  }
  
  # Exportar tabla de significancia completa
  if (length(signif_rows) > 0) {
    signif_table <- do.call(rbind, signif_rows)
    write.csv(signif_table, "significancia_modelos.csv", row.names = FALSE)
    cat("\n  -> Taula esportatua: significancia_modelos.csv\n")
  }

  cat("\n  Kodeak: *** p<0.001  ** p<0.01  * p<0.05  . p<0.1  ns ez esanguratsua\n")
  
  # Seleccionar el mejor
  # Regla de parsimonia: si hay modelos dentro de ΔAIC < 2 del mejor,
  # escoger el de menor número de términos (edf efectivo)
  aic_min <- min(model_comp$AIC)
  candidatos <- model_comp %>% filter(AIC < aic_min + 2)
  if (nrow(candidatos) > 1) {
    # contar términos en la fórmula como proxy de complejidad
    candidatos$n_terminos <- stringr::str_count(candidatos$formula, "\\+") + 1
    candidatos <- candidatos %>% arrange(n_terminos, AIC)
    cat(sprintf("\n  ΔAIC<2: %d eredu berdinduta; parsimoniatsuena hautatu da (%s)\n",
                nrow(candidatos), candidatos$id[1]))
  }
  best_id <- candidatos$id[1]
  best_mod <- models[[best_id]]
  best_label <- labels[[best_id]]
  # Recuperar tipo desde model_comp usando el id seleccionado
  best_tipo <- model_comp$tipo[match(best_id, model_comp$id)]
  uses_spei <- grepl("SPEI", best_label)
  uses_rad  <- grepl("Rad", best_label)
  
  cat("\n", rep("=", 85), "\n")
  cat("  HAUTATUTAKO EREDUA:", best_label, "\n")
  cat("  Mota:", best_tipo, "| R2:", model_comp$R2[1],
      "| Dev.expl:", model_comp$dev_expl[1],
      "| AIC:", model_comp$AIC[1],
      "| p:", model_comp$p_global[1], "\n")
  cat(rep("=", 85), "\n")
  
  # ── Exportar predicciones GAM por fecha (para IPP espacializado en GEE) ──
  # Si el modelo seleccionado contiene termino tensorial te(), la extrapolacion
  # fuera del hull convexo (GDA, P60) puede producir predicciones inestables.
  # Se chequea el hull del entrenamiento y se marcan como NA las fechas fuera.
  tiene_tensor <- grepl("te\\(", best_label)

  if (tiene_tensor) {
    cat("\n  Termino tentsorialdun eredua detektatuta: oskol konbexua egiaztatzen...\n")
    rng_gda_train <- range(best_mod$model$GDA, na.rm = TRUE)
    rng_p60_train <- range(best_mod$model$P60, na.rm = TRUE)

    gam_export <- lai_aoi %>%
      mutate(
        fuera_hull    = GDA < rng_gda_train[1] | GDA > rng_gda_train[2] |
                        P60 < rng_p60_train[1] | P60 > rng_p60_train[2],
        LAI_norm_pred = ifelse(fuera_hull, NA_real_, predict(best_mod, newdata = .)),
        residuo_medio = LAI_aoi - LAI_norm_pred
      )

    n_fuera <- sum(gam_export$fuera_hull, na.rm = TRUE)
    n_total <- nrow(gam_export)
    if (n_fuera > 0) {
      cat(sprintf("  OHARRA: %d data %d-tik (%.1f%%) oskoletik kanpo (GDA, P60).\n",
                  n_fuera, n_total, 100 * n_fuera / n_total))
      cat("         NA gisa markatuta eta GEErako CSVtik kanpo utziak.\n")
    } else {
      cat("  Data guztiak entrenamendu-oskolaren barruan.\n")
    }

    gam_export <- gam_export %>%
      filter(!is.na(LAI_norm_pred)) %>%
      select(date, year, doy, LAI_norm_obs = LAI_aoi, LAI_norm_pred, residuo_medio)

  } else {
    # Modelo sin tensor: comportamiento estandar del v9.3
    gam_export <- lai_aoi %>%
      mutate(
        LAI_norm_pred = predict(best_mod, newdata = .),
        residuo_medio = LAI_aoi - LAI_norm_pred
      ) %>%
      filter(!is.na(LAI_norm_pred)) %>%
      select(date, year, doy, LAI_norm_obs = LAI_aoi, LAI_norm_pred, residuo_medio)
  }

  write.csv(gam_export, "gam_predicciones_fecha.csv", row.names = FALSE)
  cat("\n  -> Esportatua: gam_predicciones_fecha.csv (", nrow(gam_export), " data)\n")
  cat("     GEEn aktibora igo A estrategiaren LPI espazializatua kalkulatzeko\n")

  # ── Chequeo anti-absorcion: si el modelo tiene termino estacional, verificar
  # que no esta absorbiendo la senal de pastoreo diferencial entre anillos ──
  tiene_estacional <- grepl("doy|month|season", best_label, ignore.case = TRUE)
  if (tiene_estacional && exists("panel_zonal")) {
    cat("\n  [!] Sasoiko terminodun eredua: hurbil-urrun diferentziala egiaztatzen...\n")
    panel_ago <- panel_zonal %>%
      filter(grepl("pasto", zona), month(date) == 8, !is.na(LAI_norm)) %>%
      mutate(anillo = case_when(
        grepl("cercano", zona)    ~ "cer",
        grepl("intermedio", zona) ~ "int",
        grepl("remoto", zona)     ~ "rem",
        TRUE ~ NA_character_))
    if (nrow(panel_ago) > 0) {
      panel_ago$pred <- predict(best_mod, newdata = panel_ago)
      panel_ago$res  <- panel_ago$LAI_norm - panel_ago$pred
      medias <- aggregate(res ~ anillo, data = panel_ago, FUN = mean, na.rm = TRUE)
      if (all(c("cer", "rem") %in% medias$anillo)) {
        dif <- abs(medias$res[medias$anillo == "cer"] -
                   medias$res[medias$anillo == "rem"])
        cat(sprintf("      |hurbila - urruna| diferentziala abuztuan: %.3f\n", dif))
        if (dif < 0.08) {
          cat("      [!!] DIFERENTZIALA KOLAPSATUTA. Ereduak larratze-diferentziala xurgatzen du.\n")
          cat("           LPIren diagnostiko-ahalmena kaltetuta egon daiteke.\n")
        } else {
          cat("      [ok] Diferentziala mantenduta.\n")
        }
      }
    }
  }
  
} else {
  stop("Ezin izan da inolako eredurik doitu.")
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. ESTRATEGIA A — Modelo climático directo (sobre LAI medio AOI)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  A ESTRATEGIA — Klima-eredua LAI normalizatuaren gainean\n")
cat(rep("=", 60), "\n")

panel_A <- panel
if (uses_spei) panel_A <- panel_A %>% filter(!is.na(SPEI_3))
if (uses_rad)  panel_A <- panel_A %>% filter(!is.na(Rad7))

panel_A <- panel_A %>%
  mutate(
    LAI_norm_clim_A = predict(best_mod, newdata = .),
    LAI_norm_res_A  = LAI_norm - LAI_norm_clim_A,
    LAI_clim_A = LAI_norm_clim_A * LAI_max_zy,
    LAI_res_A  = LAI - LAI_clim_A
  )

write.csv(panel_A, "panel_A_con_residuos.csv", row.names = FALSE)

cat("A hondarrarekin proiekzioak:", nrow(panel_A), "erregistro\n")
cat("Batezbesteko hondar normalizatua:", round(mean(panel_A$LAI_norm_res_A, na.rm = TRUE), 4), "\n")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. ESTRATEGIA B — Pseudo-referencia por accesibilidad
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  B ESTRATEGIA — Pseudo-erreferentzia (", ZONA_REFERENCIA, ")\n")
cat(rep("=", 60), "\n")

ref_data <- panel %>% filter(zona == ZONA_REFERENCIA, !is.na(P60))
cat("Erreferentzia datuak:", nrow(ref_data), "beh\n")

if (nrow(ref_data) >= 10) {
  # Modelo base LM
  mod_B <- lm(LAI_norm ~ P60, data = ref_data)

  # LM mas complejo
  tryCatch({
    mod_B2 <- lm(LAI_norm ~ P60 + GDA, data = ref_data)
    if (AIC(mod_B2) < AIC(mod_B) - 2) mod_B <- mod_B2
  }, error = function(e) {})

  tryCatch({
    mod_B3 <- lm(LAI_norm ~ P60 + GDA + GDA2, data = ref_data)
    if (AIC(mod_B3) < AIC(mod_B) - 2) mod_B <- mod_B3
  }, error = function(e) {})

  # GAM sobre referencia
  tryCatch({
    mod_B_gam <- gam(LAI_norm ~ s(GDA, k = 5) + s(P60, k = 5), data = ref_data)
    if (AIC(mod_B_gam) < AIC(mod_B) - 2) {
      mod_B <- mod_B_gam
      cat("  B estrategia: GAM hautatuta (AIC hobea)\n")
    }
  }, error = function(e) {})

  cat("B eredua:", ifelse(inherits(mod_B, "gam"), "GAM", "LM"), "\n")
  if (inherits(mod_B, "gam")) {
    cat("  R2 =", round(summary(mod_B)$r.sq, 3), "\n")
  } else {
    cat("  R2 =", round(summary(mod_B)$r.squared, 3), "\n")
  }

  panel_A <- panel_A %>%
    mutate(
      LAI_norm_clim_B = predict(mod_B, newdata = .),
      LAI_norm_res_B  = LAI_norm - LAI_norm_clim_B,
      LAI_clim_B = LAI_norm_clim_B * LAI_max_zy,
      LAI_res_B  = LAI - LAI_clim_B
    )
} else {
  cat("Datu nahikorik ez B estrategiarako. A eta C bakarrik erabiliko dira.\n")
  panel_A$LAI_norm_clim_B <- NA
  panel_A$LAI_norm_res_B  <- NA
  panel_A$LAI_clim_B <- NA
  panel_A$LAI_res_B  <- NA
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. ESTRATEGIA C — Anomalía espacial
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  C ESTRATEGIA — Anomalia espaziala\n")
cat(rep("=", 60), "\n")

panel_A <- panel_A %>%
  group_by(date) %>%
  mutate(
    # Anomalía espacial sobre LAI normalizado
    LAI_norm_mediana = median(LAI_norm, na.rm = TRUE),
    LAI_norm_anomalia = LAI_norm - LAI_norm_mediana,
    LAI_anomalia_rel = ifelse(LAI_norm_mediana > 0,
                              (LAI_norm - LAI_norm_mediana) / LAI_norm_mediana, 0),
    # Conservar anomalía absoluta para gráficos
    LAI_mediana = median(LAI, na.rm = TRUE),
    LAI_anomalia = LAI - LAI_mediana
  ) %>%
  ungroup()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7b. ESTRATEGIA D — Referencia hayedo (nivel 2 del modelo dual)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Usa la anomalía del hayedo como termómetro climático. Si el hayedo baja
# por sequía, el pasto debería bajar lo mismo; lo que no se explica por
# la anomalía del hayedo es atribuible a presión ganadera.
#
# CRÍTICO: hay que comparar anomalías con anomalías, no niveles con anomalías.
# Para eso calculamos PRIMERO la anomalía del pasto (respecto a su propia
# fenología media por zona) y DESPUÉS se la restamos la anomalía del hayedo:
#   LAI_pasto_anom = LAI_norm_pasto − fenologia_media_pasto_por_zona(doy)
#   LAI_norm_res_D = LAI_pasto_anom − LAI_hayedo_anom
#   IPP_D          = mean(LAI_norm_res_D) / mean(LAI_norm_pasto)
# Complementa A/B/C sin depender de ajustes estadísticos climáticos.

cat("\n\n", rep("=", 60), "\n")
cat("  D ESTRATEGIA — Pagadi-erreferentzia (larretu gabea)\n")
cat(rep("=", 60), "\n")

if (hayedo_disponible && "LAI_hayedo_anom" %in% names(panel_A)) {

  # 1. Fenología media multianual del pasto POR ZONA
  #    Ajustamos un GAM suave por zona para tener la curva típica de cada una.
  cat("\n  Larrearen anomalietarako batezbesteko fenologia zonaka doitzen...\n")

  panel_A <- panel_A %>%
    group_by(zona) %>%
    mutate(
      LAI_norm_pheno = {
        obs <- sum(!is.na(LAI_norm))
        if (obs >= 10) {
          mod_ph <- tryCatch(
            gam(LAI_norm ~ s(doy, k = 6), data = data.frame(LAI_norm, doy)),
            error = function(e) NULL)
          if (!is.null(mod_ph)) {
            predict(mod_ph, newdata = data.frame(doy = doy))
          } else rep(NA_real_, length(LAI_norm))
        } else rep(NA_real_, length(LAI_norm))
      },
      LAI_norm_anom = LAI_norm - LAI_norm_pheno
    ) %>%
    ungroup()

  # 2. Estrategia D: diferencia de anomalías (pasto respecto a su fenología
  #    menos la anomalía climática capturada por el hayedo)
  panel_A <- panel_A %>%
    mutate(
      LAI_norm_res_D = LAI_norm_anom - LAI_hayedo_anom,
      # Versión en unidades absolutas (retro-escalada por LAI_max_zy)
      LAI_res_D      = LAI_norm_res_D * LAI_max_zy
    )

  n_validos_D <- sum(!is.na(panel_A$LAI_norm_res_D))
  cat("LPI_D kalkulagarriak diren errenkadak:", n_validos_D, "/", nrow(panel_A), "\n")
  if (n_validos_D > 0) {
    cat("LAI_norm_res_D tartea (larre-anomalia - pagadi-anomalia):",
        round(min(panel_A$LAI_norm_res_D, na.rm = TRUE), 3),
        "—",
        round(max(panel_A$LAI_norm_res_D, na.rm = TRUE), 3), "\n")
    cat("LAI_norm_anom tartea (larrea vs bere batezbesteko fenologia):",
        round(min(panel_A$LAI_norm_anom, na.rm = TRUE), 3),
        "—",
        round(max(panel_A$LAI_norm_anom, na.rm = TRUE), 3), "\n")
  }
} else {
  cat("Pagadia ez dago eskuragarri: D estrategia baztertzen da.\n")
  panel_A$LAI_norm_anom  <- NA_real_
  panel_A$LAI_norm_pheno <- NA_real_
  panel_A$LAI_norm_res_D <- NA_real_
  panel_A$LAI_res_D      <- NA_real_
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 8. DIAGNÓSTICO INTEGRADO — IPP
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  3. MAILA — DIAGNOSTIKO INTEGRATUA\n")
cat(rep("=", 60), "\n")

# Función clasificación IPP
clasificar_ipp <- function(ipp) {
  case_when(
    ipp < IPP_THRESHOLDS["severo_neg"]   ~ "Gainlarratze larria",
    ipp < IPP_THRESHOLDS["moderado_neg"] ~ "Gainlarratze ertaina",
    ipp < IPP_THRESHOLDS["leve_neg"]     ~ "Gainlarratze arina",
    ipp <= IPP_THRESHOLDS["leve_pos"]    ~ "Larre-oreka",
    ipp <= IPP_THRESHOLDS["moderado_pos"]~ "Azpierabilera arina",
    ipp <= IPP_THRESHOLDS["severo_pos"]  ~ "Azpierabilera ertaina",
    TRUE                                  ~ "Sastrakatzea / utzikeria"
  )
}

# ── IPP por zona y año (sobre LAI normalizado) ──
compute_ipp <- function(df, periodo_label) {
  df %>%
    group_by(zona) %>%
    summarise(
      periodo = periodo_label,
      n = n(),
      LAI_medio = round(mean(LAI, na.rm = TRUE), 3),
      LAI_norm_medio = round(mean(LAI_norm, na.rm = TRUE), 3),
      LAI_max   = round(max(LAI, na.rm = TRUE), 3),
      # Estrategia A (sobre LAI normalizado)
      IPP_A = ifelse(mean(LAI_norm, na.rm = TRUE) > 0,
                     round(mean(LAI_norm_res_A, na.rm = TRUE) / mean(LAI_norm, na.rm = TRUE), 4), NA),
      # Estrategia B (sobre LAI normalizado)
      IPP_B = ifelse(!all(is.na(LAI_norm_res_B)) & mean(LAI_norm, na.rm = TRUE) > 0,
                     round(mean(LAI_norm_res_B, na.rm = TRUE) / mean(LAI_norm, na.rm = TRUE), 4), NA),
      # Estrategia C (anomalía espacial, ya relativa)
      IPP_C = round(mean(LAI_anomalia_rel, na.rm = TRUE), 4),
      # Estrategia D (referencia hayedo): LAI_norm_res_D ya es una diferencia
      # de anomalías (anom_pasto - anom_hayedo), adimensional y centrada en 0.
      # No se divide por mean(LAI_norm) porque el numerador no es un nivel.
      #
      # PARCHE v9.3: para el periodo "Global" la Estrategia D es, por construcción
      # matemática, ≈0 (la media de los residuos de un GAM multianual sobre su
      # serie de ajuste es cero). Eso contamina el IPP_consenso global y la
      # métrica de concordancia. Se devuelve NA en el Global para que rowMeans
      # la excluya y la concordancia pase a contarse sobre A, B y C (k/3).
      # En los periodos anuales IPP_D sí se calcula porque ahí es informativa.
      IPP_D = ifelse(periodo_label == "Global" | all(is.na(LAI_norm_res_D)),
                     NA_real_,
                     round(mean(LAI_norm_res_D, na.rm = TRUE), 4)),
      # Proporción fechas con residuo negativo
      prop_neg_A = round(mean(LAI_norm_res_A < 0, na.rm = TRUE), 2),
      .groups = "drop"
    ) %>%
    mutate(
      # Consenso = media de las estrategias disponibles (ignora NAs)
      IPP_consenso = round(rowMeans(cbind(IPP_A, IPP_B, IPP_C, IPP_D), na.rm = TRUE), 4),
      # Concordancia: cuántas de las estrategias disponibles coinciden en signo
      concordancia_n = mapply(function(a, b, c, d) {
        vals <- c(a, b, c, d)
        vals <- vals[!is.na(vals)]
        if (length(vals) < 2) return(NA_character_)
        neg <- sum(vals < 0); pos <- sum(vals > 0); zero <- sum(vals == 0)
        max_concord <- max(neg, pos, zero)
        sprintf("%d/%d", max_concord, length(vals))
      }, IPP_A, IPP_B, IPP_C, IPP_D),
      diagnostico = clasificar_ipp(IPP_consenso),
      concordancia = ifelse(
        !is.na(IPP_A) & !is.na(IPP_C),
        ifelse(sign(IPP_A) == sign(IPP_C) &
               (is.na(IPP_B) | sign(IPP_A) == sign(IPP_B)) &
               (is.na(IPP_D) | sign(IPP_A) == sign(IPP_D)),
               "Handia", "Txikia"),
        "N/A"
      )
    )
}

# Global (todos los años)
diag_global <- compute_ipp(panel_A, "Global")

# Por año
diag_anual <- bind_rows(
  lapply(sort(unique(panel_A$year)), function(yr) {
    compute_ipp(panel_A %>% filter(year == yr), as.character(yr))
  })
)

# Unir
diag_all <- bind_rows(diag_global, diag_anual) %>%
  arrange(periodo, zona)

cat("\n━━━ DIAGNOSTIKOA ZONA ETA EPEKA ━━━\n")
print(as.data.frame(diag_all %>%
  select(periodo, zona, n, LAI_medio, LAI_norm_medio,
         IPP_A, IPP_B, IPP_C, IPP_D, IPP_consenso,
         diagnostico, concordancia, concordancia_n)))

write.csv(diag_all, "diagnostico_pastoral.csv", row.names = FALSE)
cat("\n✓ Esportatua: diagnostico_pastoral.csv\n")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 9. GRÁFICOS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  GRAFIKOAK\n")
cat(rep("=", 60), "\n")

# ── 9a. LAI observado vs. esperado (Estrategia A) por zona, facetado ──
for (z in zonas) {
  sub <- panel_A %>% filter(zona == z, !is.na(LAI_clim_A))
  if (nrow(sub) < 2) next

  p_obs <- ggplot(sub, aes(x = date)) +
    geom_ribbon(aes(ymin = LAI_clim_A * (1 - TOLERANCIA),
                    ymax = LAI_clim_A * (1 + TOLERANCIA)),
                fill = "grey80", alpha = 0.4) +
    geom_line(aes(y = LAI_clim_A), color = "#984ea3", linewidth = 1, linetype = "dashed") +
    geom_line(aes(y = LAI), color = colores[z], linewidth = 1.1) +
    geom_point(aes(y = LAI), color = colores[z], size = 3) +
    facet_wrap(~year, scales = "free_x") +
    labs(title = paste0(nombres_zona[z], ": LAI behatua vs. itxarondakoa (A estr.)"),
         subtitle = paste0("Eredua: ", best_label, " | Banda grisa: \u00b1", TOLERANCIA*100, "%"),
         y = "LAI (m\u00b2/m\u00b2)", x = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  p_res <- ggplot(sub, aes(x = date, y = LAI_res_A, fill = LAI_res_A > 0)) +
    geom_col(width = 5, show.legend = FALSE) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = c("TRUE" = "#4daf4a", "FALSE" = "#e41a1c")) +
    facet_wrap(~year, scales = "free_x") +
    labs(title = paste0("LAI hondarra: ", nombres_zona[z]),
         y = "LAI hondarra", x = "Data") +
    theme_minimal(base_size = 12)

  combined <- p_obs / p_res
  ggsave(paste0("diag_A_", z, ".png"), combined, width = 14, height = 9, dpi = 150)
  print(combined)
}

# ── 9b. Estrategia C: anomalías espaciales ──
p_C <- ggplot(panel_A, aes(x = date, y = LAI_anomalia, fill = zona)) +
  geom_col(position = position_dodge(width = 6), width = 5) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colores, labels = nombres_zona) +
  facet_wrap(~year, scales = "free_x") +
  labs(title = "C estrategia: LAI anomalia AOI medianarekiko",
       subtitle = "Medianatik behera = presio handiagoa | Medianatik gora = presio txikiagoa",
       y = "LAI - mediana(LAI)", x = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("diag_C_anomalias.png", p_C, width = 14, height = 7, dpi = 150)
print(p_C)

# ── 9c. Panel de diagnóstico integrado (IPP consenso) ──
p_diag <- ggplot(diag_global, aes(x = reorder(zona, IPP_consenso),
                                   y = IPP_consenso,
                                   fill = IPP_consenso < 0)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = c(-0.20, 0.20), linetype = "dotted", color = "grey40") +
  geom_text(aes(label = diagnostico),
            hjust = ifelse(diag_global$IPP_consenso < 0, 1.05, -0.05),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#c44e52", "FALSE" = "#2b8c2b")) +
  scale_x_discrete(labels = nombres_zona) +
  coord_flip() +
  labs(title = "Larre-diagnostiko integratua — Aralar (Globala)",
       subtitle = "LPI adostasuna = batezbestekoa(A, B, C, D estr.)",
       x = NULL, y = "LPI (Larre-Presio Indizea)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("diagnostico_integrado_global.png", p_diag, width = 11, height = 7, dpi = 150)
print(p_diag)

# ── 9d. Evolución interanual del IPP ──
diag_years <- diag_anual %>%
  filter(periodo != "Global") %>%
  mutate(year = as.integer(periodo))

if (nrow(diag_years) > 0 && length(unique(diag_years$year)) > 1) {
  p_evol <- ggplot(diag_years, aes(x = year, y = IPP_consenso, color = zona, group = zona)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -0.05, ymax = 0.05,
             fill = "#4daf4a", alpha = 0.1) +
    geom_hline(yintercept = c(-0.20, 0.20), linetype = "dotted", color = "grey50") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    geom_text(aes(label = round(IPP_consenso, 2)), vjust = -1, size = 3) +
    scale_color_manual(values = colores, labels = nombres_zona) +
    scale_x_continuous(breaks = unique(diag_years$year)) +
    labs(title = "LPIaren urtearteko bilakaera — Aralar",
         subtitle = "Banda berdea = oreka (\u00b10.05) | Marra puntukatuak = larria (\u00b10.20)",
         x = "Urtea", y = "LPI adostasuna", color = NULL) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  ggsave("evolucion_IPP_interanual.png", p_evol, width = 11, height = 7, dpi = 150)
  print(p_evol)
}

# ── 9e. IPP interanual barras agrupadas ──
if (nrow(diag_years) > 0) {
  p_bars <- ggplot(diag_years, aes(x = zona, y = IPP_consenso, fill = factor(year))) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = "grey60") +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -0.05, ymax = 0.05,
             fill = "#4daf4a", alpha = 0.08) +
    scale_fill_brewer(palette = "Set2") +
    scale_x_discrete(labels = nombres_zona) +
    coord_flip() +
    labs(title = "LPIa zonaka eta urteka — Aralar",
         x = NULL, y = "LPI adostasuna", fill = "Urtea") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave("IPP_barras_interanual.png", p_bars, width = 11, height = 7, dpi = 150)
  print(p_bars)
}

# ── 9f. Gradiente IPP vs. distancia a bordas ──
# Las zonas cercanas a bordas deberían tener IPP más negativo
diag_pasto <- diag_global %>%
  filter(grepl("pasto", zona)) %>%
  mutate(
    anillo = case_when(
      grepl("cercano", zona)     ~ "Hurbila (<500m)",
      grepl("intermedio", zona)  ~ "Ertaina (500-1500m)",
      grepl("remoto", zona)      ~ "Urruna (>1500m)"
    ),
    anillo = factor(anillo, levels = c("Hurbila (<500m)",
                                        "Ertaina (500-1500m)",
                                        "Urruna (>1500m)"))
  )

if (nrow(diag_pasto) > 0) {
  p_grad <- ggplot(diag_pasto, aes(x = anillo, y = IPP_consenso, fill = anillo)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = "grey60") +
    geom_text(aes(label = round(IPP_consenso, 3)),
              vjust = ifelse(diag_pasto$IPP_consenso < 0, 1.5, -0.5),
              size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("Hurbila (<500m)" = "#d73027",
                                  "Ertaina (500-1500m)" = "#fee08b",
                                  "Urruna (>1500m)" = "#1a9850")) +
    labs(title = "Larre-LPIa eta bordetarainoko distantzia — Aralar",
         subtitle = "LPIa negatiboagoa bordetatik hurbil = larre-seinale koherentea",
         x = "Bordetarainoko distantzia-eraztuna", y = "LPI adostasuna") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("gradiente_IPP_distancia_bordas.png", p_grad, width = 10, height = 7, dpi = 150)
  print(p_grad)
}

# ── 9f-bis. IPP por tipo de habitat (pasto vs brezal vs hayedo) ──
diag_tipo <- diag_global %>%
  mutate(
    tipo_hab = case_when(
      grepl("pasto", zona)   ~ "Larrea",
      grepl("brezal", zona)  ~ "Txilardia",
      grepl("hayedo", zona)  ~ "Pagadia",
      grepl("encinar", zona) ~ "Artadia",
      TRUE ~ "Bestelakoa"
    )
  )

p_tipo <- ggplot(diag_tipo, aes(x = reorder(zona, IPP_consenso),
                                 y = IPP_consenso,
                                 fill = tipo_hab)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = "grey60") +
  scale_fill_manual(values = c("Larrea" = "#33a02c", "Txilardia" = "#e31a1c",
                                "Pagadia" = "#1f78b4", "Artadia" = "#ff7f00")) +
  scale_x_discrete(labels = nombres_zona) +
  coord_flip() +
  labs(title = "LPIa habitat eta eraztunka — Aralar (Globala)",
       x = NULL, y = "LPI adostasuna", fill = "Habitat mota") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("IPP_por_habitat.png", p_tipo, width = 11, height = 8, dpi = 150)
print(p_tipo)


# ── 9h. Triangulación de estrategias A/B/C/D por zona ──
# Muestra los cuatro IPP lado a lado para cada zona pastoreada.
# Si los cuatro convergen en signo → diagnóstico robusto.
# Si divergen → advertir al lector que el resultado es sensible al método.
diag_long <- diag_global %>%
  filter(grepl("pasto", zona)) %>%
  select(zona, IPP_A, IPP_B, IPP_C, IPP_D) %>%
  pivot_longer(cols = c(IPP_A, IPP_B, IPP_C, IPP_D),
               names_to = "estrategia", values_to = "IPP") %>%
  mutate(
    estrategia = factor(estrategia,
      levels = c("IPP_A", "IPP_B", "IPP_C", "IPP_D"),
      labels = c("A: GAM klimatikoa",
                 "B: Barneko pseudo-erreferentzia",
                 "C: Anomalia espaziala",
                 "D: Pagadi-erreferentzia"))
  ) %>%
  filter(!is.na(IPP))

if (nrow(diag_long) > 0) {
  p_tri <- ggplot(diag_long, aes(x = zona, y = IPP, fill = estrategia)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.6) +
    geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed",
               color = "grey60") +
    scale_fill_manual(values = c(
      "A: GAM klimatikoa"             = "#1b9e77",
      "B: Barneko pseudo-erreferentzia" = "#d95f02",
      "C: Anomalia espaziala"         = "#7570b3",
      "D: Pagadi-erreferentzia"       = "#e7298a"
    )) +
    geom_text(aes(label = round(IPP, 3)),
              position = position_dodge(width = 0.8),
              vjust = ifelse(diag_long$IPP < 0, 1.3, -0.4),
              size = 2.8) +
    scale_x_discrete(labels = nombres_zona) +
    labs(title = "LPIaren triangulazio metodologikoa — Aralarko larreak",
         subtitle = "4 estrategia independente; zeinuan bat etortzea = diagnostiko sendoa",
         x = NULL, y = "LPI", fill = "Estrategia") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top",
          axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave("IPP_triangulacion_ABCD.png", p_tri, width = 12, height = 7, dpi = 150)
  print(p_tri)
  cat("\n✓ 9h grafikoa: IPP_triangulacion_ABCD.png\n")
}


# ── 9g. Validación cruzada leave-one-year-out ──
cat("\n━━━ Leave-one-year-out balidazioa ━━━\n")
years_avail <- sort(unique(lai_aoi$year))

if (length(years_avail) >= 3) {
  loyo_results <- data.frame()

  for (yr_out in years_avail) {
    train <- lai_aoi %>% filter(year != yr_out)
    test  <- lai_aoi %>% filter(year == yr_out)

    if (nrow(train) < 5 || nrow(test) < 3) next

    tryCatch({
      # Usar la misma clase de modelo (lm o gam) que el seleccionado
      if (inherits(best_mod, "gam")) {
        mod_cv <- gam(formula(best_mod), data = train)
      } else {
        mod_cv <- lm(formula(best_mod), data = train)
      }
      pred <- predict(mod_cv, newdata = test)
      obs  <- test$LAI_aoi
      valid <- !is.na(pred) & !is.na(obs)
      rmse <- sqrt(mean((pred[valid] - obs[valid])^2))
      r <- cor(pred[valid], obs[valid])

      loyo_results <- bind_rows(loyo_results, data.frame(
        year_out = yr_out,
        n_train = nrow(train),
        n_test = sum(valid),
        RMSE = round(rmse, 3),
        r = round(r, 3)
      ))
    }, error = function(e) {
      cat("  LOYO", yr_out, "huts egin du:", e$message, "\n")
    })
  }

  if (nrow(loyo_results) > 0) {
    cat("\n")
    print(as.data.frame(loyo_results))
    cat("LOYO batezbesteko RMSE:", round(mean(loyo_results$RMSE), 3), "\n")
    cat("LOYO batezbesteko r:", round(mean(loyo_results$r), 3), "\n")
    write.csv(loyo_results, "validacion_loyo.csv", row.names = FALSE)
  }
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 10. RESUMEN FINAL
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n\n", rep("=", 60), "\n")
cat("  LABURPENA\n")
cat(rep("=", 60), "\n\n")

cat("Hautatutako klima-eredua:", best_label, "\n")
cat("R2adj:", round(summary(best_mod)$adj.r.squared, 3), "\n\n")

cat("━━━ Diagnostiko globala ━━━\n")
diag_global %>%
  select(zona, IPP_consenso, diagnostico, concordancia) %>%
  as.data.frame() %>% print()

# Test de coherencia ecológica: gradiente de distancia a bordas
ipp_cercano <- diag_global$IPP_consenso[diag_global$zona == "pasto_cercano"]
ipp_interm  <- diag_global$IPP_consenso[diag_global$zona == "pasto_intermedio"]
ipp_remoto  <- diag_global$IPP_consenso[diag_global$zona == "pasto_remoto"]

cat("\nLPI larre hurbila (<500m):",   round(ipp_cercano, 3), "\n")
cat("LPI larre ertaina (500-1500m):", round(ipp_interm, 3), "\n")
cat("LPI larre urruna (>1500m):",    round(ipp_remoto, 3), "\n")

if (length(ipp_cercano) > 0 && length(ipp_remoto) > 0) {
  if (ipp_cercano < ipp_remoto) {
    cat("-> KOHERENTEA: LPI negatiboagoa bordetatik hurbil = larre-seinalea berretsia.\n")
    cat("   Gradientea: ", round(ipp_remoto - ipp_cercano, 3),
        " LPI unitate eraztuneko.\n")
  } else {
    cat("-> KONTUZ: itxarondakoaren aurkako patroia. Distantzia-atalaseak berrikusi.\n")
  }
}

# Coherencia por tipo de habitat
ipp_pastos  <- mean(diag_global$IPP_consenso[grepl("pasto", diag_global$zona)], na.rm = TRUE)
ipp_brezal  <- mean(diag_global$IPP_consenso[grepl("brezal", diag_global$zona)], na.rm = TRUE)
ipp_hayedo  <- mean(diag_global$IPP_consenso[grepl("hayedo", diag_global$zona)], na.rm = TRUE)

cat("\nLPI batezbestekoa habitat motaren arabera:\n")
cat("  Larreak:", round(ipp_pastos, 3), "\n")
if (!is.na(ipp_brezal)) cat("  Txilardiak:", round(ipp_brezal, 3), "\n")
if (!is.na(ipp_hayedo)) cat("  Pagadiak:", round(ipp_hayedo, 3), "\n")

cat("\n━━━ Esportatutako fitxategiak ━━━\n")
cat("  comparacion_modelos.csv\n")
cat("  diagnostico_pastoral.csv\n")
for (f in list.files(pattern = "\\.png$")) cat(" ", f, "\n")

cat("\n", rep("=", 60), "\n")
cat("  Diagnostikoa burututa.\n")
cat(rep("=", 60), "\n")
