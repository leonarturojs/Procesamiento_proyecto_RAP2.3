"""
=============================================================================
ANÁLISIS DE COEXPRESIÓN GÉNICA — RAP2.3 / AT3G16770
Dataset: PEP1 en raíz Arabidopsis thaliana (Mock y PEP1 × 3H, 6H, 12H, 24H)
=============================================================================

QUÉ ANALIZA ESTE SCRIPT:

ETAPA 1 — Carga y preprocesamiento
  • Lee 16 archivos .xlsx (TPM por muestra)
  • Filtra genes con >20% NAs; imputa con mediana del gen
  • log2(TPM+1); selecciona top 5000 genes por varianza
  → Fig01: Dendrograma de muestras (QC)

ETAPA 2 — Red de coexpresión (WGCNA-equivalente)
  • Correlación Pearson → Adjacency firmada → TOM
  • Clustering jerárquico sobre 1-TOM → módulos
  → Fig02: Selección del soft-thresholding power (R² vs β)
  → Fig03: Dendrograma de genes con barra de colores de módulos

ETAPA 3 — Módulo de RAP2.3, eigengenes, correlación con traits
  → Fig04: Heatmap correlación módulo–PEP1/tiempo

ETAPA 4 — Puntuaciones y figuras de análisis centradas en GOI
  • GS = correlación gen–PEP1
  • MM = correlación gen–eigengene del módulo
  • kWithin = conectividad intra-modular
  • Pearson r de todos los genes vs AT3G16770 + p-valor + FDR
  → Fig05: GS vs MM (posición de RAP2.3 resaltada)
  → Fig06: Perfil temporal del eigengene (Mock vs PEP1)
  → Fig07: Heatmap top 50 genes del módulo (Z-score)
  → Fig08: Volcano r vs -log10(FDR)
  → Fig09: Barplot top 20 co-expresados + 20 anti-correlacionados
  → Fig10: Grilla de scatter GOI vs top vecinos
  → Fig11: kWithin vs MM (hub genes)
  → Tabla_Correlaciones_RAP2-3.csv
  → Tabla_Modulo_RAP2-3.csv

ETAPA 5 — Redes de coexpresión
  → Fig12: Red vecindad directa RAP2.3 (TOM ≥ 0.10)
  → Fig13: Red módulo completo top 150 genes (color = GS)
  → Fig14: Red de módulos (eigengenes como nodos)
  → Fig15: Red vecindad coloreada por Fold Change PEP1 vs Mock 24H
  + Cytoscape (.graphml + _edges.csv + _nodes.csv) para Red 1, 2 y 4

CÓMO CORRER:
  1. Abre cmd o PowerShell en la carpeta del proyecto
  2. python coexpresion_WGCNA.py
  Si se interrumpe: vuelve a correr — retoma desde el último checkpoint.
=============================================================================
"""

# ─── INSTALACIÓN AUTOMÁTICA DE DEPENDENCIAS ──────────────────────────────────
import subprocess, sys, os

# ── Recursión: dendrogramas de 5000 genes requieren profundidad mayor ─────────
sys.setrecursionlimit(100_000)

_PKGS = {
    "numpy": "numpy", "pandas": "pandas", "scipy": "scipy",
    "sklearn": "scikit-learn", "matplotlib": "matplotlib",
    "seaborn": "seaborn", "openpyxl": "openpyxl",
    "networkx": "networkx", "statsmodels": "statsmodels",
    "threadpoolctl": "threadpoolctl",
}

def _install(pkg_dict):
    import importlib
    for imp, pip in pkg_dict.items():
        try:
            importlib.import_module(imp)
        except ImportError:
            print(f"  Instalando {pip}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", pip, "-q"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"  ✓ {pip}")

print("=" * 65)
print("  Verificando dependencias...")
_install(_PKGS)
print("  ✓ Todo listo\n")

# ─── IMPORTS ─────────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
import seaborn as sns
import networkx as nx
import pickle, warnings, gc
from pathlib import Path
from scipy import stats
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ─── OPTIMIZACIÓN DE CPU ──────────────────────────────────────────────────────
# Usa todos los núcleos físicos para numpy/scipy (BLAS/LAPACK multihilo)
_N_CPUS = os.cpu_count() or 1
_N_WORKERS = max(1, _N_CPUS // 2)          # mitad de núcleos (conservador)
os.environ.setdefault("OMP_NUM_THREADS",    str(_N_WORKERS))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(_N_WORKERS))
os.environ.setdefault("MKL_NUM_THREADS",    str(_N_WORKERS))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(_N_WORKERS))
# threadpoolctl aplica el límite a BLAS en tiempo de ejecución
try:
    from threadpoolctl import threadpool_limits
    _POOL_CTX = threadpool_limits(limits=_N_WORKERS)
except Exception:
    _POOL_CTX = None

# ─── ESTILO GLOBAL ───────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.framealpha": 0.85, "legend.fontsize": 9,
})

PAL = {
    "GOI": "#D32F2F", "pos": "#1565C0", "neg": "#2E7D32",
    "hub": "#6A1B9A", "alta_gs": "#E65100", "otros": "#78909C",
    "Mock": "#1976D2", "PEP1": "#C62828",
}

# ─── PARÁMETROS ──────────────────────────────────────────────────────────────
DATA_DIR = Path(r"C:\Users\LeonA\OneDrive\Documents\Proyecto\Transcriptómica\PEP1 en raíz AT\WGCNA")
OUT_DIR  = DATA_DIR

GOI          = "AT3G16770"
GOI_NAME     = "RAP2.3"
N_TOP_GENES  = 5000
MIN_MOD_SIZE = 30
TOM_THRESH   = 0.10
TOM_STRICT   = 0.15
MAX_NODES    = 150
TOP_BAR      = 20
TIME_ORDER   = ["3H", "6H", "12H", "24H"]
TP_COLORS    = {"3H": "#FFF176", "6H": "#FFB300", "12H": "#E65100", "24H": "#B71C1C"}
TRT_COLORS   = {"Mock": PAL["Mock"], "PEP1": PAL["PEP1"]}

print(f"  Carpeta : {DATA_DIR}")
print(f"  GOI     : {GOI} ({GOI_NAME})")
print(f"  CPUs    : {_N_WORKERS} de {_N_CPUS} (núcleos BLAS)\n")

# ─── UTILIDADES ──────────────────────────────────────────────────────────────
def cp_save(obj, name):
    with open(OUT_DIR / f"{name}.pkl", "wb") as f:
        pickle.dump(obj, f, protocol=4)
    print(f"    ✓ Checkpoint → {name}.pkl")

def cp_load(name):
    p = OUT_DIR / f"{name}.pkl"
    if p.exists():
        print(f"  ✓ Checkpoint '{name}.pkl' encontrado — etapa saltada")
        with open(p, "rb") as f:
            return pickle.load(f)
    return None

def savefig(fig, name):
    fig.savefig(OUT_DIR / f"{name}.pdf")
    fig.savefig(OUT_DIR / f"{name}.png", dpi=150)
    plt.close(fig)
    print(f"    → {name}.pdf / .png")

def axlabels(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, labelpad=8)
    ax.set_ylabel(ylabel, labelpad=8)
    if title:
        ax.set_title(title, pad=12)

# =============================================================================
# ETAPA 1 — CARGA Y PREPROCESAMIENTO
# =============================================================================
cp1 = cp_load("cp1_datExpr")

if cp1 is None:
    print("\n" + "=" * 65)
    print("  ETAPA 1: Carga y preprocesamiento")
    print("=" * 65)

    FILE_META = [
        ("GSM8001046_Mock1_24H.xlsx",   "Mock_24H",   "Mock", "24H"),
        ("GSM8001047_PEP1_24H_1.xlsx",  "PEP1_24H_1", "PEP1", "24H"),
        ("GSM8001048_PEP1_24H_2.xlsx",  "PEP1_24H_2", "PEP1", "24H"),
        ("GSM8001049_PEP1_24H_3.xlsx",  "PEP1_24H_3", "PEP1", "24H"),
        ("GSM8001050_Mock2_12H.xlsx",   "Mock_12H",   "Mock", "12H"),
        ("GSM8001051_PEP1_12H_1.xlsx",  "PEP1_12H_1", "PEP1", "12H"),
        ("GSM8001052_PEP1_12H_2.xlsx",  "PEP1_12H_2", "PEP1", "12H"),
        ("GSM8001053_PEP1_12H_3.xlsx",  "PEP1_12H_3", "PEP1", "12H"),
        ("GSM8001054_Mock3_3H.xlsx",    "Mock_3H",    "Mock",  "3H"),
        ("GSM8001055_PEP1_3H_1.xlsx",   "PEP1_3H_1",  "PEP1",  "3H"),
        ("GSM8001056_PEP1_3H_2.xlsx",   "PEP1_3H_2",  "PEP1",  "3H"),
        ("GSM8001057_PEP1_3H_3.xlsx",   "PEP1_3H_3",  "PEP1",  "3H"),
        ("GSM8001058_Mock4_6H.xlsx",    "Mock_6H",    "Mock",  "6H"),
        ("GSM8001059_PEP1_6H_1.xlsx",   "PEP1_6H_1",  "PEP1",  "6H"),
        ("GSM8001060_PEP1_6H_2.xlsx",   "PEP1_6H_2",  "PEP1",  "6H"),
        ("GSM8001061_PEP1_6H_3.xlsx",   "PEP1_6H_3",  "PEP1",  "6H"),
    ]
    sample_meta = pd.DataFrame(FILE_META,
                               columns=["file","sample","treatment","timepoint"])

    print("  Leyendo Excel (paralelo)...")
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _read_one(row):
        df = (pd.read_excel(DATA_DIR / row["file"], usecols=["Gene ID","TPM"])
                .rename(columns={"Gene ID":"gene_id","TPM":row["sample"]})
                .dropna(subset=["gene_id"]).set_index("gene_id"))
        return row["sample"], row["file"], df

    futures_map = {}
    with ThreadPoolExecutor(max_workers=min(_N_WORKERS, len(sample_meta))) as ex:
        for _, r in sample_meta.iterrows():
            fut = ex.submit(_read_one, r)
            futures_map[fut] = r["file"]

    frames_dict = {}
    for fut in as_completed(futures_map):
        sample_name, fname, df = fut.result()
        frames_dict[sample_name] = df
        print(f"    ✓ {fname}")
    # Reordenar según sample_meta
    frames = [frames_dict[s] for s in sample_meta["sample"]]

    expr = pd.concat(frames, axis=1); del frames; gc.collect()
    print(f"  Genes: {len(expr):,} | Muestras: {len(expr.columns)}")

    expr = expr[expr.isna().mean(axis=1) <= 0.20]
    expr = expr.apply(lambda r: r.fillna(r.median()), axis=1)
    expr_log = np.log2(expr + 1); del expr; gc.collect()

    variances = expr_log.var(axis=1)
    top_genes = variances.nlargest(N_TOP_GENES).index.tolist()
    if GOI not in top_genes:
        top_genes = top_genes[:N_TOP_GENES-1] + [GOI]
        print(f"  {GOI} añadido manualmente.")

    datExpr = expr_log.loc[top_genes].T
    del expr_log; gc.collect()
    print(f"  datExpr: {datExpr.shape[0]} muestras × {datExpr.shape[1]} genes")

    # Fig01: Dendrograma de muestras
    Z_s = linkage(datExpr.values, method="average")
    fig, ax = plt.subplots(figsize=(14, 5))
    dendrogram(Z_s, labels=datExpr.index.tolist(), ax=ax,
               leaf_rotation=40, leaf_font_size=9,
               color_threshold=0.7 * max(Z_s[:, 2]))
    ax.axhline(15, color="red", ls="--", lw=1, label="Umbral outliers (h=15)")
    axlabels(ax,
        xlabel="Muestra",
        ylabel="Distancia euclidiana (sobre log₂(TPM+1), adim.)",
        title="Fig. 01 — Control de calidad: clustering de muestras\n"
              "Ramas cortas = muestras similares. Línea roja = umbral sugerido de outliers.")
    ax.legend(loc="upper right")
    fig.tight_layout(); savefig(fig, "Fig01_SampleClustering")

    cp_save({"datExpr": datExpr, "sample_meta": sample_meta}, "cp1_datExpr")
    print("  ✓ ETAPA 1 completada\n")
else:
    datExpr     = cp1["datExpr"]
    sample_meta = cp1["sample_meta"]
    del cp1; gc.collect()


# =============================================================================
# ETAPA 2 — RED DE COEXPRESIÓN (TOM + MÓDULOS)
# =============================================================================
cp2 = cp_load("cp2_network")

if cp2 is None:
    print("\n" + "=" * 65)
    print("  ETAPA 2: Red de coexpresión (paso más lento)")
    print("=" * 65)

    X = datExpr.values.astype(np.float32)
    genes = datExpr.columns.tolist()

    # Soft-thresholding power
    print("  Evaluando soft-thresholding powers...")
    cor_mat = np.corrcoef(X.T).astype(np.float32)

    def scale_free_r2(C, pw):
        A = ((1 + C) / 2) ** pw; np.fill_diagonal(A, 0)
        k = A.sum(axis=1); k = k[k > 0]
        lk = np.log10(k)
        h, ed = np.histogram(lk, bins=20)
        bc = (ed[:-1] + ed[1:]) / 2; mask = h > 0
        if mask.sum() < 3: return 0.0, A.mean()
        r, *_ = stats.linregress(bc[mask], np.log10(h[mask]+1))
        return r**2, A.mean()

    rows = []
    for pw in range(1, 21):
        r2, mk = scale_free_r2(cor_mat, pw)
        rows.append({"power": pw, "R2": r2, "MeanK": mk})
        print(f"    β={pw:2d}  R²={r2:.3f}  MeanK={mk:.1f}")
    fit_df = pd.DataFrame(rows)

    cands = fit_df[fit_df["R2"] >= 0.85]
    soft_power = int(cands["power"].iloc[0]) if len(cands) else \
                 int(fit_df.loc[fit_df["R2"].idxmax(), "power"])
    print(f"  → Power elegido: {soft_power}")

    # Fig02: Soft-thresholding
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.plot(fit_df["power"], fit_df["R2"], "o-", color=PAL["pos"], lw=1.8, ms=6)
    ax1.axhline(0.85, color=PAL["PEP1"], ls="--", lw=1, label="Umbral R²=0.85")
    ax1.axvline(soft_power, color=PAL["neg"], ls="--", lw=1.2,
                label=f"β elegido = {soft_power}")
    ax1.set_ylim(0, 1.05)
    axlabels(ax1,
        xlabel="Soft-thresholding power β (adim.)",
        ylabel="R² de ajuste a red libre de escala (adim.)\n(más alto = mejor topología libre de escala)",
        title="Ajuste a red libre de escala\nSe busca R² ≥ 0.85 con β mínimo")
    ax1.legend()

    ax2.plot(fit_df["power"], fit_df["MeanK"], "o-", color=PAL["alta_gs"], lw=1.8, ms=6)
    ax2.axvline(soft_power, color=PAL["neg"], ls="--", lw=1.2, label=f"β = {soft_power}")
    axlabels(ax2,
        xlabel="Soft-thresholding power β (adim.)",
        ylabel="Conectividad media (suma de pesos de adyacencia, adim.)\n"
               "Cae al aumentar β (red más dispersa)",
        title="Conectividad media de la red\nCaída suave = buen umbralizado")
    ax2.legend()
    fig.suptitle("Fig. 02 — Selección del soft-thresholding power (WGCNA)\n"
                 "β óptimo: la red tiene topología libre de escala (R²≥0.85) "
                 "sin perder demasiada conectividad",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout(); savefig(fig, "Fig02_SoftThreshold")

    # TOM
    print("  Calculando TOM...")
    adj = ((1 + cor_mat) / 2).astype(np.float32) ** soft_power
    np.fill_diagonal(adj, 0); del cor_mat; gc.collect()
    k = adj.sum(axis=1)
    # Multiplicación en bloques para reducir picos de memoria en matrices grandes
    _CHUNK = 1000
    n_genes_tom = adj.shape[0]
    adj_sq = np.zeros_like(adj)
    for _i in range(0, n_genes_tom, _CHUNK):
        adj_sq[_i:_i+_CHUNK] = adj[_i:_i+_CHUNK] @ adj
    num = adj_sq + adj; del adj_sq; gc.collect()
    den = np.minimum.outer(k, k) + 1 - adj
    den = np.where(den == 0, 1e-10, den)
    TOM_full = (num / den).astype(np.float32)
    np.fill_diagonal(TOM_full, 1.0)
    del adj, num, den; gc.collect()
    np.save(OUT_DIR / "TOM_full.npy", TOM_full)

    # Clustering
    print("  Clustering jerárquico...")
    dissim = squareform((1 - TOM_full).astype(np.float64), checks=False)
    Z_gene = linkage(dissim, method="average"); del dissim; gc.collect()

    cut = 0.25
    for t in np.arange(0.1, 1.0, 0.05):
        lbl = fcluster(Z_gene, t=t, criterion="distance")
        sz = pd.Series(lbl).value_counts()
        if sz.max() < len(genes)*0.5 and (sz >= MIN_MOD_SIZE).sum() >= 3:
            cut = t; break

    raw = fcluster(Z_gene, t=cut, criterion="distance")
    ls  = pd.Series(raw, index=genes)
    small = ls.value_counts()[ls.value_counts() < MIN_MOD_SIZE].index
    ls[ls.isin(small)] = 0
    ls = ls.map({old: new for new, old in enumerate(sorted(ls.unique()))})

    CLIST = ["grey","turquoise","blue","brown","yellow","green","red","black",
             "pink","magenta","purple","greenyellow","tan","salmon","cyan",
             "midnightblue","lightcyan","grey60","lightgreen","royalblue",
             "darkred","darkgreen","darkturquoise","orange","darkorange",
             "steelblue","paleturquoise","violet","darkolivegreen","coral"]
    module_colors = ls.map(lambda l: "grey" if l == 0 else CLIST[l % (len(CLIST)-1)])
    print(f"  Módulos: {(module_colors != 'grey').nunique()}")
    print(module_colors.value_counts().head(10).to_string())

    # Fig03: Dendrograma + módulos
    # truncate_mode="lastp" colapsa el árbol en los p nodos finales,
    # evitando la recursión de ~5000 niveles que causa RecursionError.
    _DEND_P = min(300, len(genes))   # nodos visibles; ajusta si quieres más detalle
    # Orden completo de hojas sin recursión (usa leaves_list de scipy)
    from scipy.cluster.hierarchy import leaves_list as _leaves_list
    _full_order = _leaves_list(Z_gene).tolist()
    cols_ord = module_colors.iloc[_full_order].values
    fig, axes = plt.subplots(2, 1, figsize=(16, 7),
                             gridspec_kw={"height_ratios": [5, 0.55], "hspace": 0.04})
    dendrogram(Z_gene, no_labels=True, ax=axes[0],
               link_color_func=lambda _: "#546E7A",
               truncate_mode="lastp", p=_DEND_P)
    axes[0].set_ylabel("Distancia (1 − TOM, adim.)\n0 = idéntico | 1 = sin solapamiento",
                       fontsize=10)
    axes[0].set_title(
        "Fig. 03 — Dendrograma de genes y módulos de coexpresión\n"
        "Cada hoja = un gen. Ramas cortas = genes con patrón de expresión similar. "
        "Colores = módulo asignado.",
        fontweight="bold", fontsize=12)
    axes[0].set_xticks([])
    for i, c in enumerate(cols_ord):
        axes[1].add_patch(plt.Rectangle((i, 0), 1, 1, color=c, lw=0))
    axes[1].set_xlim(0, len(genes)); axes[1].set_ylim(0, 1)
    axes[1].set_yticks([0.5]); axes[1].set_yticklabels(["Módulo"], fontsize=9)
    axes[1].set_xticks([])
    fig.tight_layout(); savefig(fig, "Fig03_ModuleDendrogram")

    cp_save({"soft_power": soft_power, "module_colors": module_colors,
             "Z_gene": Z_gene, "gene_names": genes}, "cp2_network")
    del Z_gene; gc.collect()
    print("  ✓ ETAPA 2 completada\n")
else:
    soft_power    = cp2["soft_power"]
    module_colors = cp2["module_colors"]
    gene_names    = cp2["gene_names"]
    del cp2; gc.collect()


# =============================================================================
# ETAPA 3 — MÓDULO GOI, EIGENGENES, CORRELACIÓN CON TRAITS
# =============================================================================
cp3 = cp_load("cp3_tom_module")

if cp3 is None:
    print("\n" + "=" * 65)
    print("  ETAPA 3: Módulo de RAP2.3, eigengenes y TOM")
    print("=" * 65)

    goi_mod   = module_colors[GOI]
    goi_genes = module_colors[module_colors == goi_mod].index.tolist()
    print(f"  {GOI_NAME} → módulo '{goi_mod}'  ({len(goi_genes)} genes)")

    # Eigengenes
    unique_mods = [m for m in module_colors.unique() if m != "grey"]
    ME_dict = {}
    for mod in unique_mods:
        gm = module_colors[module_colors == mod].index
        Xm = datExpr[gm].values; Xm -= Xm.mean(axis=0)
        pc = PCA(n_components=1).fit_transform(Xm)[:, 0]
        if np.corrcoef(pc, Xm.mean(axis=1))[0, 1] < 0: pc = -pc
        ME_dict[f"ME_{mod}"] = pc
    ME_df = pd.DataFrame(ME_dict, index=datExpr.index)

    # Traits
    smi = sample_meta.set_index("sample").reindex(datExpr.index).copy()
    smi["tp_num"] = smi["timepoint"].str.replace("H","").astype(float)
    smi["is_PEP1"] = (smi["treatment"] == "PEP1").astype(float)
    pep1_s = pd.Series(smi["is_PEP1"].values, index=ME_df.index)

    # Correlación módulo–trait
    me_r, me_p = {}, {}
    for col in ME_df.columns:
        for tr in ["is_PEP1","tp_num"]:
            r, p = pearsonr(ME_df[col], smi[tr])
            me_r.setdefault(col, {})[tr] = r
            me_p.setdefault(col, {})[tr] = p
    mod_r = pd.DataFrame(me_r).T.rename(columns={"is_PEP1":"PEP1","tp_num":"Tiempo"})
    mod_p = pd.DataFrame(me_p).T.rename(columns={"is_PEP1":"PEP1","tp_num":"Tiempo"})
    mod_r.index = mod_r.index.str.replace("ME_","")
    mod_p.index = mod_p.index.str.replace("ME_","")
    order_r = mod_r["PEP1"].abs().sort_values(ascending=False).index
    mod_r = mod_r.loc[order_r]; mod_p = mod_p.loc[order_r]

    annot = mod_r.round(2).astype(str) + "\n(p=" + \
            mod_p.applymap(lambda v: f"{v:.1e}").astype(str) + ")"

    # Fig04
    fig, ax = plt.subplots(figsize=(7, max(8, len(mod_r)*0.38)))
    sns.heatmap(mod_r.astype(float), annot=annot, fmt="",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                linewidths=0.4, ax=ax, cbar_kws={"label": "r de Pearson (adim.)"},
                annot_kws={"size": 8},
                xticklabels=["PEP1\n(1=PEP1, 0=Mock)", "Tiempo\n(horas)"],
                yticklabels=mod_r.index)
    if goi_mod in list(mod_r.index):
        pos = list(mod_r.index).index(goi_mod)
        ax.add_patch(plt.Rectangle((0, pos), 2, 1, fill=False,
                                   edgecolor="gold", lw=2.5))
        ax.text(2.05, pos+0.5, f"← {GOI_NAME}", va="center",
                fontsize=8, color="goldenrod", fontweight="bold")
    axlabels(ax,
        xlabel="Trait experimental",
        ylabel="Módulo de coexpresión (colores WGCNA)",
        title=f"Fig. 04 — Correlación módulo–tratamiento/tiempo\n"
              f"r de Pearson entre el eigengene y el trait experimental\n"
              f"Recuadro dorado = módulo que contiene {GOI_NAME} ({GOI})")
    fig.tight_layout(); savefig(fig, "Fig04_ModuleTrait_Heatmap")

    # TOM del módulo
    tom_path = OUT_DIR / "TOM_full.npy"
    if tom_path.exists():
        print("  Cargando TOM del disco...")
        TOM_full = np.load(tom_path).astype(np.float32)
        all_g = gene_names
        idx = [all_g.index(g) for g in goi_genes if g in all_g]
        sub_g = [all_g[i] for i in idx]
        TOM_module = pd.DataFrame(TOM_full[np.ix_(idx, idx)],
                                  index=sub_g, columns=sub_g)
        del TOM_full; gc.collect()
    else:
        print(f"  Recalculando TOM para {len(goi_genes)} genes...")
        Xm  = datExpr[goi_genes].values.T.astype(np.float32)
        cm  = np.corrcoef(Xm).astype(np.float32)
        am  = ((1+cm)/2)**soft_power; np.fill_diagonal(am, 0)
        km  = am.sum(axis=1)
        nm  = am@am+am; dm = np.minimum.outer(km,km)+1-am
        dm  = np.where(dm==0,1e-10,dm)
        tm  = (nm/dm).astype(np.float32); np.fill_diagonal(tm,1.0)
        TOM_module = pd.DataFrame(tm, index=goi_genes, columns=goi_genes)
        del cm,am,km,nm,dm,tm; gc.collect()

    cp_save({"goi_mod": goi_mod, "goi_genes": goi_genes,
             "ME_df": ME_df, "smi": smi, "pep1_s": pep1_s,
             "TOM_module": TOM_module}, "cp3_tom_module")
    print("  ✓ ETAPA 3 completada\n")
else:
    goi_mod    = cp3["goi_mod"]
    goi_genes  = cp3["goi_genes"]
    ME_df      = cp3["ME_df"]
    smi        = cp3["smi"]
    pep1_s     = cp3["pep1_s"]
    TOM_module = cp3["TOM_module"]
    del cp3; gc.collect()


# =============================================================================
# ETAPA 4 — GS, MM, CORRELACIONES Y FIGURAS DE ANÁLISIS
# =============================================================================
cp4 = cp_load("cp4_scores")

if cp4 is None:
    print("\n" + "=" * 65)
    print("  ETAPA 4: GS, MM, conectividad, correlaciones y figuras")
    print("=" * 65)

    ME_goi = ME_df[f"ME_{goi_mod}"].values
    pep1_v = smi["is_PEP1"].values
    X   = datExpr.values.astype(np.float64)
    Xc  = X - X.mean(axis=0)
    pc  = pep1_v - pep1_v.mean()
    mc  = ME_goi  - ME_goi.mean()
    Xn  = np.sqrt((Xc**2).sum(axis=0))

    GS  = pd.Series((Xc*pc[:,None]).sum(0)/(Xn*np.sqrt((pc**2).sum())),
                    index=datExpr.columns)
    MM  = pd.Series((Xc*mc[:,None]).sum(0)/(Xn*np.sqrt((mc**2).sum())),
                    index=datExpr.columns)

    # Correlación directa de todos los genes con GOI
    print("  Calculando correlaciones con GOI...")
    gv  = datExpr[GOI].values.astype(np.float64)
    gc_ = gv - gv.mean()
    gn  = np.sqrt((gc_**2).sum())
    r_all = pd.Series((Xc*gc_[:,None]).sum(0)/(Xn*gn), index=datExpr.columns)
    n  = len(datExpr)
    t  = r_all * np.sqrt(n-2) / np.sqrt(1 - r_all**2 + 1e-15)
    pv = pd.Series(2*stats.t.sf(np.abs(t), df=n-2), index=datExpr.columns)
    _, fdr_v, _, _ = multipletests(pv.values, method="fdr_bh")
    fdr = pd.Series(fdr_v, index=datExpr.columns)

    cor_df = (pd.DataFrame({
        "gene_id": r_all.index, "pearson_r": r_all.values,
        "pvalue": pv.values, "FDR": fdr.values,
        "module": module_colors.reindex(r_all.index).values,
    })[lambda d: d["gene_id"] != GOI]
    .sort_values("pearson_r", ascending=False).reset_index(drop=True))
    cor_df.to_csv(OUT_DIR / "Tabla_Correlaciones_RAP2-3.csv", index=False)
    print("    ✓ Tabla_Correlaciones_RAP2-3.csv")

    # Conectividad intra-modular
    print("  Conectividad intra-modular...")
    Xm  = datExpr[goi_genes].values.T.astype(np.float32)
    am  = ((1+np.corrcoef(Xm))/2).astype(np.float32)**soft_power
    np.fill_diagonal(am,0)
    kW  = am.sum(axis=1); del am,Xm; gc.collect()

    kIM = pd.DataFrame({
        "gene_id": goi_genes, "kWithin": kW,
        "MM": MM[goi_genes].values, "GS": GS[goi_genes].values,
        "pearson_r_GOI": r_all[goi_genes].values,
        "is_GOI": [g==GOI for g in goi_genes],
    }).sort_values("kWithin", ascending=False).reset_index(drop=True)
    kIM["hub"]    = kIM["kWithin"] >= kIM["kWithin"].quantile(0.9)
    kIM["module"] = goi_mod
    kIM.to_csv(OUT_DIR / "Tabla_Modulo_RAP2-3.csv", index=False)
    print("    ✓ Tabla_Modulo_RAP2-3.csv")

    gsm = kIM[["gene_id","GS","MM","is_GOI"]].copy()
    del Xc, Xn, X; gc.collect()

    # ── Fig05: GS vs MM ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 7))
    oth = gsm[~gsm["is_GOI"]]
    gpt = gsm[gsm["is_GOI"]]
    sc5 = ax.scatter(oth["MM"], oth["GS"], c=oth["GS"].abs(), cmap="plasma",
                     alpha=0.55, s=18, vmin=0, vmax=1,
                     label=f"Otros genes del módulo (n={len(oth)})")
    plt.colorbar(sc5, ax=ax,
                 label="|GS| = |correlación gen–PEP1| (adim., 0–1)")
    if len(gpt):
        ax.scatter(gpt["MM"], gpt["GS"], s=220, color=PAL["GOI"],
                   zorder=6, marker="*", label=f"★ {GOI_NAME} ({GOI})")
        ax.annotate(f"  ★ {GOI_NAME}\n  ({GOI})",
                    (float(gpt["MM"]), float(gpt["GS"])),
                    fontsize=10, fontweight="bold", color=PAL["GOI"])
    ax.axhline(0, color="grey", lw=0.7, ls="--")
    ax.axvline(0, color="grey", lw=0.7, ls="--")
    ax.axvline(0.8, color="grey", lw=0.8, ls=":", label="MM=0.8 (umbral hub)")
    axlabels(ax,
        xlabel="Module Membership (MM)\n"
               "Correlación del gen con el eigengene del módulo (adim., -1 a 1)\n"
               "Mayor MM = gen más representativo del módulo",
        ylabel="Gene Significance (GS)\n"
               "Correlación del gen con el tratamiento PEP1 (adim., -1 a 1)\n"
               "Mayor |GS| = gen más responde a PEP1",
        title=f"Fig. 05 — GS vs MM en el módulo '{goi_mod}'\n"
              f"Genes en la esquina sup. derecha = co-expresados Y responden a PEP1")
    ax.legend(loc="upper left")
    fig.tight_layout(); savefig(fig, "Fig05_GS_vs_MM")

    # ── Fig06: Perfil temporal del eigengene ─────────────────────────────────
    me_plot = pd.DataFrame({
        "eigengene": ME_goi,
        "treatment": smi["treatment"].values,
        "timepoint": pd.Categorical(smi["timepoint"].values,
                                    categories=TIME_ORDER, ordered=True)
    })
    summ = (me_plot.groupby(["timepoint","treatment"], observed=True)["eigengene"]
            .agg(["mean","sem"]).reset_index())

    fig, ax = plt.subplots(figsize=(8, 5))
    for trt, col in TRT_COLORS.items():
        s = summ[summ["treatment"]==trt].sort_values("timepoint")
        ax.plot(s["timepoint"].astype(str), s["mean"], "o-",
                color=col, lw=2, ms=8, label=trt)
        ax.errorbar(s["timepoint"].astype(str), s["mean"],
                    yerr=s["sem"], fmt="none", color=col, capsize=5, lw=1.5)
    ax.axhline(0, color="grey", lw=0.7, ls="--")
    axlabels(ax,
        xlabel="Tiempo post-tratamiento (horas)",
        ylabel="Valor del eigengene modular (adim.)\n"
               "= 1ª componente principal del módulo\n"
               "(representa el patrón dominante de expresión)",
        title=f"Fig. 06 — Perfil temporal del eigengene — módulo '{goi_mod}'\n"
              f"Módulo que contiene {GOI_NAME}. Puntos = media ± SEM de réplicas.")
    ax.legend(title="Tratamiento")
    fig.tight_layout(); savefig(fig, "Fig06_Eigengene_Profile")

    # ── Fig07: Heatmap top 50 genes del módulo ───────────────────────────────
    top50 = kIM.head(50)["gene_id"].tolist()
    if GOI not in top50: top50 = [GOI] + top50[:49]

    col_ord = smi.sort_values(["treatment","timepoint"]).index.tolist()
    hm_raw  = datExpr[top50].T
    hm_z    = pd.DataFrame(StandardScaler().fit_transform(hm_raw.T).T,
                           index=hm_raw.index, columns=hm_raw.columns)[col_ord]
    col_c   = pd.DataFrame({
        "Tratamiento": smi.loc[col_ord,"treatment"].map(TRT_COLORS),
        "Tiempo":      smi.loc[col_ord,"timepoint"].map(TP_COLORS),
    }, index=col_ord)
    row_lbl = [f"★ {g} ({GOI_NAME})" if g==GOI else g for g in hm_z.index]

    cg = sns.clustermap(
        hm_z, col_colors=col_c, row_cluster=True, col_cluster=False,
        cmap="RdBu_r", center=0, vmin=-3, vmax=3,
        yticklabels=row_lbl, xticklabels=True,
        figsize=(15, 17), dendrogram_ratio=0.12, colors_ratio=0.025,
        cbar_pos=(0.01, 0.78, 0.025, 0.15),
        method="average", metric="correlation")
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_yticklabels(),
                                   fontsize=6.5, fontstyle="italic")
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xticklabels(),
                                   fontsize=7.5, rotation=40, ha="right")
    cg.ax_heatmap.set_xlabel("Muestra (ordenadas: Mock→PEP1, 3H→24H)", labelpad=8)
    cg.ax_heatmap.set_ylabel("Gen (★ = RAP2.3, top 50 por kWithin)", labelpad=8)
    cg.ax_cbar.set_ylabel("Z-score log₂(TPM+1)\n(adim.; azul=baja, rojo=alta expresión)",
                          fontsize=8)
    cg.figure.suptitle(
        f"Fig. 07 — Heatmap del módulo '{goi_mod}': top 50 genes por kWithin\n"
        f"Color = Z-score de log₂(TPM+1) | ★ = {GOI_NAME} ({GOI})\n"
        f"Barra superior: azul=Mock, rojo=PEP1 | amarillo=3H → rojo oscuro=24H",
        y=1.01, fontsize=11, fontweight="bold")
    cg.figure.savefig(OUT_DIR / "Fig07_Heatmap_Modulo.pdf", bbox_inches="tight")
    cg.figure.savefig(OUT_DIR / "Fig07_Heatmap_Modulo.png", bbox_inches="tight", dpi=150)
    plt.close(); print("    → Fig07_Heatmap_Modulo.pdf / .png")

    # ── Fig08: Volcano ───────────────────────────────────────────────────────
    vc = cor_df.copy()
    vc["neg_log10_FDR"] = -np.log10(vc["FDR"].clip(lower=1e-300))
    vc["cat"] = "No significativo"
    vc.loc[(vc["pearson_r"] >  0.5) & (vc["FDR"] < 0.05), "cat"] = "Co-expresado (r>0.5)"
    vc.loc[(vc["pearson_r"] < -0.5) & (vc["FDR"] < 0.05), "cat"] = "Anti-correlado (r<-0.5)"
    cat_col = {"Co-expresado (r>0.5)": PAL["pos"],
               "Anti-correlado (r<-0.5)": PAL["neg"],
               "No significativo": PAL["otros"]}

    top_lab = pd.concat([vc[vc["cat"]=="Co-expresado (r>0.5)"].head(12),
                         vc[vc["cat"]=="Anti-correlado (r<-0.5)"].tail(12)])

    fig, ax = plt.subplots(figsize=(10, 7))
    for cat, col in cat_col.items():
        s = vc[vc["cat"]==cat]
        ax.scatter(s["pearson_r"], s["neg_log10_FDR"],
                   color=col, s=8 if "sig" in cat.lower() else 18,
                   alpha=0.4 if "sig" in cat.lower() else 0.75,
                   label=f"{cat} (n={len(s):,})")
    for _, row in top_lab.iterrows():
        ax.annotate(row["gene_id"], (row["pearson_r"], row["neg_log10_FDR"]),
                    fontsize=6.5, xytext=(3,3), textcoords="offset points")
    ax.axhline(-np.log10(0.05), color=PAL["PEP1"], ls="--", lw=0.9,
               label="FDR = 0.05")
    ax.axvline(0.5,  color=PAL["pos"], ls=":", lw=0.8)
    ax.axvline(-0.5, color=PAL["neg"], ls=":", lw=0.8)
    ax.text(0.51,  ax.get_ylim()[1]*0.97, "r=+0.5", fontsize=8, color=PAL["pos"])
    ax.text(-0.51, ax.get_ylim()[1]*0.97, "r=−0.5", fontsize=8,
            color=PAL["neg"], ha="right")
    axlabels(ax,
        xlabel=f"Coeficiente de correlación de Pearson r con {GOI_NAME} ({GOI})\n"
               "(adim.; -1 = anticorrelación perfecta, 0 = sin correlación, +1 = correlación perfecta)",
        ylabel="-log₁₀(FDR)\nFDR = tasa de falso descubrimiento (Benjamini-Hochberg)\n"
               "Valores mayores = más significativo estadísticamente",
        title=f"Fig. 08 — Volcano de correlaciones de todos los genes con {GOI_NAME} ({GOI})\n"
              f"n = {len(vc):,} genes | Línea horizontal = FDR 0.05 | Líneas verticales = |r| 0.5")
    ax.legend(loc="upper left", markerscale=2)
    fig.tight_layout(); savefig(fig, "Fig08_Volcano_Correlaciones")

    # ── Fig09: Barplot top correlaciones ─────────────────────────────────────
    top_p = cor_df.head(TOP_BAR).copy(); top_p["tipo"] = "Co-expresado"
    top_n = cor_df.tail(TOP_BAR).copy(); top_n["tipo"] = "Anti-correlado"
    top_bar = pd.concat([top_n, top_p])

    fig, ax = plt.subplots(figsize=(10, 9))
    bc = [PAL["neg"] if t=="Anti-correlado" else PAL["pos"] for t in top_bar["tipo"]]
    ax.barh(range(len(top_bar)), top_bar["pearson_r"],
            color=bc, height=0.72, edgecolor="white", lw=0.3)
    ax.set_yticks(range(len(top_bar)))
    ax.set_yticklabels(top_bar["gene_id"], fontsize=8, fontstyle="italic")
    for i, (_, row) in enumerate(top_bar.iterrows()):
        ha = "left" if row["pearson_r"] < 0 else "right"
        xo = -0.005 if row["pearson_r"] < 0 else 0.005
        ax.text(row["pearson_r"]+xo, i, f"{row['pearson_r']:.3f}",
                va="center", ha=ha, fontsize=7)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlim(-1.05, 1.05)
    pp = mpatches.Patch(color=PAL["pos"], label=f"Co-expresados (top {TOP_BAR})")
    np_ = mpatches.Patch(color=PAL["neg"], label=f"Anti-correlados (top {TOP_BAR})")
    ax.legend(handles=[pp, np_], loc="lower right")
    axlabels(ax,
        xlabel=f"Coeficiente de Pearson r con {GOI_NAME} ({GOI})\n"
               "(adim.; calculado sobre 16 muestras, log₂(TPM+1))\n"
               "-1=anticorrelación perfecta | 0=sin correlación | +1=correlación perfecta",
        ylabel="Gen (ID TAIR de Arabidopsis thaliana)",
        title=f"Fig. 09 — Top {TOP_BAR} genes co-expresados y anti-correlacionados "
              f"con {GOI_NAME} ({GOI})\n"
              f"Verde = anti-correlación | Azul = co-expresión | "
              f"Valor sobre cada barra = r exacto")
    fig.tight_layout(); savefig(fig, "Fig09_Barplot_Correlaciones")

    # ── Fig10: Scatter GOI vs top vecinos (grilla 3×4) ───────────────────────
    tp6 = cor_df.head(6)["gene_id"].tolist()
    tn6 = cor_df.tail(6)["gene_id"].tolist()
    panel = tp6 + tn6
    goi_v = datExpr[GOI].values

    fig, axes = plt.subplots(3, 4, figsize=(16, 12))
    for ai, (gene, ax) in enumerate(zip(panel, axes.flatten())):
        gv2    = datExpr[gene].values
        r_val  = float(cor_df.loc[cor_df["gene_id"]==gene,"pearson_r"].values[0])
        fdr_val= float(cor_df.loc[cor_df["gene_id"]==gene,"FDR"].values[0])
        col    = PAL["pos"] if r_val >= 0 else PAL["neg"]
        for trt, mkr in [("Mock","o"),("PEP1","^")]:
            msk = smi["treatment"]==trt
            ax.scatter(goi_v[msk], gv2[msk], color=TRT_COLORS[trt],
                       marker=mkr, s=35, alpha=0.8, label=trt, zorder=3)
        xf = np.linspace(goi_v.min(), goi_v.max(), 100)
        sl, ic, *_ = stats.linregress(goi_v, gv2)
        ax.plot(xf, sl*xf+ic, color=col, lw=1.5, alpha=0.7)
        ax.set_title(f"{gene}\nr={r_val:.3f}  FDR={fdr_val:.1e}",
                     fontsize=8.5, fontweight="bold", color=col)
        ax.set_xlabel(f"{GOI_NAME} — log₂(TPM+1)", fontsize=7.5)
        ax.set_ylabel(f"{gene} — log₂(TPM+1)", fontsize=7.5)
        ax.tick_params(labelsize=7)
        if ai == 0: ax.legend(title="Tratamiento", fontsize=7, title_fontsize=7)
    for ax_x in axes.flatten()[len(panel):]: ax_x.set_visible(False)
    fig.suptitle(
        f"Fig. 10 — Scatter: {GOI_NAME} ({GOI}) vs sus top vecinos\n"
        f"Filas 1–2: top 6 co-expresados | Fila 3: top 6 anti-correlacionados\n"
        f"Eje X e Y: log₂(TPM+1), adim. | ○=Mock | △=PEP1 | Línea=regresión lineal\n"
        f"Color título: azul=co-expresado, verde=anti-correlado",
        fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout(); savefig(fig, "Fig10_Scatter_Top_Vecinos")

    # ── Fig11: kWithin vs MM ─────────────────────────────────────────────────
    gm  = kIM["is_GOI"]
    fig, ax = plt.subplots(figsize=(9, 7))
    sc11 = ax.scatter(kIM.loc[~gm,"kWithin"], kIM.loc[~gm,"MM"],
                      c=kIM.loc[~gm,"GS"].abs(), cmap="plasma",
                      s=15, alpha=0.6, vmin=0, vmax=1)
    ax.scatter(kIM.loc[gm,"kWithin"], kIM.loc[gm,"MM"],
               color=PAL["GOI"], s=220, zorder=6, marker="*",
               label=f"★ {GOI_NAME} ({GOI})")
    cb11 = plt.colorbar(sc11, ax=ax)
    cb11.set_label("|GS| = |correlación gen–PEP1| (adim., 0–1)", fontsize=10)
    for _, row in kIM.head(20).iterrows():
        ax.annotate(row["gene_id"], (row["kWithin"], row["MM"]),
                    fontsize=6.5, alpha=0.85, xytext=(2,2),
                    textcoords="offset points")
    ax.axvline(kIM["kWithin"].quantile(0.9), ls=":", lw=0.9, color="grey",
               label="Percentil 90 kWithin (hub genes)")
    ax.axhline(0.8, ls=":", lw=0.9, color="#546E7A", label="MM=0.8")
    ax.legend(loc="lower right")
    axlabels(ax,
        xlabel="kWithin — conectividad intra-modular (adim.)\n"
               "Suma de los pesos de adyacencia (soft) dentro del módulo\n"
               "Mayor kWithin = hub gene (gen más central del módulo)",
        ylabel="Module Membership (MM) (adim., -1 a 1)\n"
               "Correlación del gen con el eigengene del módulo\n"
               "Mayor MM = gen más representativo del módulo",
        title=f"Fig. 11 — Conectividad intra-modular vs MM — módulo '{goi_mod}'\n"
              f"Color = |GS con PEP1| | ★ = {GOI_NAME} | Líneas punteadas = umbrales hub")
    fig.tight_layout(); savefig(fig, "Fig11_kWithin_vs_MM")

    cp_save({"GS": GS, "MM": MM, "ME_goi": ME_goi,
             "kIM": kIM, "gsm": gsm, "cor_df": cor_df}, "cp4_scores")
    print("  ✓ ETAPA 4 completada\n")
else:
    GS     = cp4["GS"]
    MM     = cp4["MM"]
    ME_goi = cp4["ME_goi"]
    kIM    = cp4["kIM"]
    gsm    = cp4["gsm"]
    cor_df = cp4["cor_df"]
    del cp4; gc.collect()

# Recargar TOM si fue liberada
if not isinstance(TOM_module, pd.DataFrame):
    cp3r = cp_load("cp3_tom_module")
    TOM_module = cp3r["TOM_module"]
    del cp3r; gc.collect()


# =============================================================================
# ETAPA 5 — REDES DE COEXPRESIÓN
# =============================================================================
print("\n" + "=" * 65)
print("  ETAPA 5: Redes de coexpresión")
print("=" * 65)


def draw_net(G, title, filename, node_attr="cat",
             fill_values=None, cbar_label=None, figsize=(14, 11)):
    """Dibuja red de coexpresión con layout FR. Etiqueta solo GOI + hubs."""
    if G.number_of_nodes() == 0:
        print(f"    ⚠ {filename}: red vacía."); return
    k_sp = max(0.5, 3.0 / np.sqrt(G.number_of_nodes()))
    pos  = nx.spring_layout(G, seed=42, k=k_sp, iterations=60)
    nds  = list(G.nodes())
    ew   = np.array([G[u][v]["weight"] for u,v in G.edges()])
    e_lw = np.clip(ew*5, 0.3, 4.0)
    e_al = np.clip(ew*6, 0.1, 0.85)
    ns   = np.array([350 if n==GOI else
                     max(25, 15+abs(float(MM.get(n,0)))*180) for n in nds])
    kw90 = kIM["kWithin"].quantile(0.85)
    lbl_set = {n for n in nds if n==GOI or
               (n in kIM["gene_id"].values and
                float(kIM.loc[kIM["gene_id"]==n,"kWithin"].values[:1]) >= kw90)}

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor("#F8F9FA")

    # Aristas
    for (u,v), lw, al in zip(G.edges(), e_lw, e_al):
        ax.plot([pos[u][0],pos[v][0]], [pos[u][1],pos[v][1]],
                color="#90A4AE", lw=lw, alpha=al, zorder=1, solid_capstyle="round")

    # Nodos
    xs = [pos[n][0] for n in nds]; ys = [pos[n][1] for n in nds]
    if fill_values is not None:
        vals = np.array([float(fill_values.get(n,0)) for n in nds])
        sc = ax.scatter(xs, ys, c=vals, cmap="RdBu_r", vmin=-1, vmax=1,
                        s=ns, zorder=3,
                        edgecolors=["gold" if n==GOI else "white" for n in nds],
                        linewidths=[2.5 if n==GOI else 0.5 for n in nds])
        plt.colorbar(sc, ax=ax, label=cbar_label or "", shrink=0.65)
    else:
        nc = [PAL["GOI"] if n==GOI
              else PAL["alta_gs"] if abs(float(GS.get(n,0))) >= 0.6
              else PAL["hub"]     if float(MM.get(n,0))      >= 0.8
              else PAL["otros"] for n in nds]
        ax.scatter(xs, ys, c=nc, s=ns, zorder=3,
                   edgecolors=["gold" if n==GOI else "white" for n in nds],
                   linewidths=[2.5 if n==GOI else 0.5 for n in nds])
        patches = [
            mpatches.Patch(color=PAL["GOI"],      label=f"★ {GOI_NAME} (GOI)"),
            mpatches.Patch(color=PAL["alta_gs"],  label="|GS PEP1| ≥ 0.6"),
            mpatches.Patch(color=PAL["hub"],      label="Hub (MM ≥ 0.8)"),
            mpatches.Patch(color=PAL["otros"],    label="Co-expresado"),
        ]
        ax.legend(handles=patches, loc="lower left", fontsize=9,
                  title="Tipo de nodo", title_fontsize=9)

    # Etiquetas (solo nodos clave)
    for nd in lbl_set:
        lbl = f"★{GOI_NAME}" if nd==GOI else nd
        ax.annotate(lbl, pos[nd], ha="center", va="bottom",
                    fontsize=7 if nd!=GOI else 9,
                    fontweight="bold" if nd==GOI else "normal",
                    fontstyle="italic", xytext=(0,6), textcoords="offset points",
                    color=PAL["GOI"] if nd==GOI else "#212121",
                    bbox=dict(boxstyle="round,pad=0.15",fc="white",alpha=0.75,lw=0)
                         if nd==GOI else None)

    ax.annotate("Grosor de arista ∝ TOM (similitud topológica)\n"
                "Tamaño de nodo ∝ |Module Membership|",
                xy=(0.01,0.01), xycoords="axes fraction",
                fontsize=8, color="#546E7A", style="italic")
    ax.set_title(title, fontweight="bold", fontsize=11, pad=12)
    ax.axis("off")
    fig.tight_layout(); savefig(fig, filename)


def export_cy(G, name, extra_cols=None):
    nx.write_graphml(G, str(OUT_DIR / f"{name}.graphml"))
    rows_e = [(u,v,d["weight"]) for u,v,d in G.edges(data=True)]
    pd.DataFrame(rows_e, columns=["source","target","TOM_weight"]).to_csv(
        OUT_DIR / f"{name}_edges.csv", index=False)
    rows_n = []
    for nd in G.nodes():
        r = {"id": nd, "is_GOI": nd==GOI,
             "GS_PEP1": float(GS.get(nd, np.nan)),
             "MM": float(MM.get(nd, np.nan)),
             "pearson_r_to_GOI": float(
                 cor_df.loc[cor_df["gene_id"]==nd,"pearson_r"].values[0])
                 if nd!=GOI and nd in cor_df["gene_id"].values else 1.0,
             "module": str(module_colors.get(nd,"grey")),
             "kWithin": float(kIM.loc[kIM["gene_id"]==nd,"kWithin"].values[0])
                        if nd in kIM["gene_id"].values else np.nan,
             "node_type": ("GOI" if nd==GOI
                           else "Alta_GS" if abs(float(GS.get(nd,0))) >= 0.6
                           else "Hub"     if float(MM.get(nd,0))       >= 0.8
                           else "Co-expresado")}
        if extra_cols: r.update({k: float(v.get(nd, np.nan))
                                  for k,v in extra_cols.items()})
        rows_n.append(r)
    pd.DataFrame(rows_n).to_csv(OUT_DIR / f"{name}_nodes.csv", index=False)
    print(f"    Cytoscape: {name}.graphml + _edges.csv + _nodes.csv")


# RED 1: Vecindad directa
print("  Red 1: vecindad directa de RAP2.3...")
G1, edges1, nodes1 = nx.Graph(), [], []
if GOI in TOM_module.index:
    nbr = TOM_module.loc[GOI][TOM_module.loc[GOI] >= TOM_THRESH].index.tolist()
    net1_g = list(dict.fromkeys([GOI]+nbr))
    TN1    = TOM_module.loc[net1_g, net1_g]
    G1.add_nodes_from(net1_g)
    for i,ga in enumerate(net1_g):
        for gb in net1_g[i+1:]:
            w = float(TN1.loc[ga,gb])
            if w >= TOM_THRESH: G1.add_edge(ga,gb,weight=w)
    draw_net(G1,
        title=(f"Fig. 12 — Red de coexpresión: vecindad directa de {GOI_NAME} ({GOI})\n"
               f"Aristas: TOM ≥ {TOM_THRESH} (similitud topológica de solapamiento)\n"
               f"n={G1.number_of_nodes()} genes | {G1.number_of_edges()} conexiones\n"
               f"★=GOI | Naranja=|GS|≥0.6 | Violeta=Hub(MM≥0.8) | Gris=co-expresado"),
        filename="Fig12_Red1_Vecindad_GOI", node_attr="cat")
    edges1 = [(u,v,d["weight"]) for u,v,d in G1.edges(data=True)]
    nodes1 = net1_g
    export_cy(G1, "Cytoscape_Red1_Vecindad_GOI")
    del TN1
else:
    print(f"    ⚠ {GOI} no está en TOM_module")


# RED 2: Módulo completo
print("  Red 2: módulo completo...")
top_m = kIM.head(MAX_NODES)["gene_id"].tolist()
if GOI not in top_m: top_m = [GOI]+top_m[:-1]
avail = [g for g in top_m if g in TOM_module.index]
TN2   = TOM_module.loc[avail, avail]
G2 = nx.Graph(); G2.add_nodes_from(avail)
for i,ga in enumerate(avail):
    for gb in avail[i+1:]:
        w = float(TN2.loc[ga,gb])
        if w >= TOM_STRICT: G2.add_edge(ga,gb,weight=w)
draw_net(G2,
    title=(f"Fig. 13 — Red del módulo '{goi_mod}' (top {MAX_NODES} genes por kWithin)\n"
           f"Aristas: TOM ≥ {TOM_STRICT} | n={G2.number_of_nodes()} genes | "
           f"{G2.number_of_edges()} conexiones\n"
           f"Color nodo = GS con PEP1 (azul=reprimido, rojo=inducido)\n"
           f"Tamaño ∝ |MM| | ★ = {GOI_NAME}"),
    filename="Fig13_Red2_Modulo_Completo",
    fill_values=GS.to_dict(),
    cbar_label="GS = r Pearson gen–PEP1 (adim., -1 a 1)\n"
               "Azul=anticorrelado con PEP1 | Rojo=co-expresado con PEP1",
    figsize=(16,13))
export_cy(G2, "Cytoscape_Red2_Modulo_Completo")
del TN2, G2; gc.collect()


# RED 3: Red de módulos
print("  Red 3: red de módulos...")
ME_THRESH = 0.5
me_cor = ME_df.corr(); me_cor.columns = me_cor.columns.str.replace("ME_","")
me_cor.index = me_cor.index.str.replace("ME_","")
mod_sz = module_colors.value_counts()
me_pep = {c.replace("ME_",""): float(ME_df[c].corr(pep1_s)) for c in ME_df.columns}

G3 = nx.Graph()
for mod in me_cor.index:
    G3.add_node(mod, n_genes=int(mod_sz.get(mod,0)),
                is_goi=(mod==goi_mod), me_pep1=me_pep.get(mod,0.0))
for i,m1 in enumerate(me_cor.index):
    for m2 in me_cor.index[i+1:]:
        r = me_cor.loc[m1,m2]
        if abs(r) >= ME_THRESH:
            G3.add_edge(m1,m2,weight=float(abs(r)),sign="pos" if r>0 else "neg")

if G3.number_of_edges() > 0:
    pos3 = nx.spring_layout(G3, seed=42, k=2.5/np.sqrt(G3.number_of_nodes()))
    n3   = list(G3.nodes())
    pv3  = np.array([G3.nodes[n]["me_pep1"] for n in n3])
    ng3  = np.array([G3.nodes[n]["n_genes"] for n in n3])
    ns3  = np.clip(ng3/max(ng3)*2800, 150, 3000)

    fig, ax = plt.subplots(figsize=(13,10))
    ax.set_facecolor("#F8F9FA")
    for (u,v,d) in G3.edges(data=True):
        col_e = PAL["PEP1"] if d["sign"]=="pos" else PAL["Mock"]
        ax.plot([pos3[u][0],pos3[v][0]],[pos3[u][1],pos3[v][1]],
                color=col_e, lw=d["weight"]*4, alpha=0.55, zorder=1)
    sc3 = ax.scatter([pos3[n][0] for n in n3], [pos3[n][1] for n in n3],
                     c=pv3, cmap="RdBu_r", vmin=-1, vmax=1, s=ns3, zorder=3,
                     edgecolors=["gold" if G3.nodes[n]["is_goi"] else "white" for n in n3],
                     linewidths=[3.0 if G3.nodes[n]["is_goi"] else 0.7 for n in n3])
    plt.colorbar(sc3, ax=ax, shrink=0.65,
                 label="Correlación eigengene–PEP1\n(r de Pearson, adim., -1 a 1)\n"
                       "Azul=anticorrelado | Rojo=co-expresado con PEP1")
    for nd in n3:
        ng = G3.nodes[nd]["n_genes"]
        ax.annotate(f"{nd}\n(n={ng})", pos3[nd], ha="center", va="center",
                    fontsize=7, fontweight="bold" if G3.nodes[nd]["is_goi"] else "normal",
                    bbox=dict(boxstyle="round,pad=0.2",fc="white",alpha=0.75,lw=0))
    ep = mpatches.Patch(color=PAL["PEP1"],label="Correlación positiva entre módulos")
    en = mpatches.Patch(color=PAL["Mock"], label="Correlación negativa entre módulos")
    eg = mpatches.Patch(fc="white",ec="gold",lw=2.5,label=f"Módulo de {GOI_NAME}")
    ax.legend(handles=[ep,en,eg],loc="lower left",fontsize=9)
    ax.set_title(
        f"Fig. 14 — Red de módulos de coexpresión\n"
        f"Cada nodo = un módulo WGCNA | Aristas: |r eigengenes| ≥ {ME_THRESH}\n"
        f"Tamaño ∝ n genes del módulo | Color = correlación eigengene–PEP1\n"
        f"Borde dorado = módulo que contiene {GOI_NAME} ({GOI})",
        fontsize=11, fontweight="bold", pad=12)
    ax.axis("off"); fig.tight_layout(); savefig(fig,"Fig14_Red3_Red_Modulos")
else:
    print(f"    Sin conexiones con ME_THRESH={ME_THRESH}")
del G3; gc.collect()


# RED 4: Vecindad coloreada por FC
print("  Red 4: vecindad con Fold Change PEP1 vs Mock 24H...")
def calc_fc(tp):
    mc = [s for s in datExpr.index if "Mock" in s and tp in s]
    pc = [s for s in datExpr.index if "PEP1" in s and tp in s]
    if not mc or not pc: return pd.Series(np.nan, index=datExpr.columns)
    return datExpr.loc[pc].mean() - datExpr.loc[mc].mean()

fc24 = calc_fc("24H")

if edges1 and GOI in TOM_module.index:
    G4 = nx.Graph(); G4.add_nodes_from(nodes1)
    for u,v,w in edges1: G4.add_edge(u,v,weight=w)
    draw_net(G4,
        title=(f"Fig. 15 — Vecindad de {GOI_NAME} ({GOI}) coloreada por "
               f"Fold Change PEP1 vs Mock a 24H\n"
               f"Color nodo = log₂FC (azul=reprimido, rojo=inducido por PEP1)\n"
               f"Tamaño ∝ |MM| | Grosor arista ∝ TOM | Borde dorado = {GOI_NAME}\n"
               f"n={G4.number_of_nodes()} genes | {G4.number_of_edges()} conexiones"),
        filename="Fig15_Red4_FC_24H",
        fill_values=fc24.to_dict(),
        cbar_label="log₂ FC PEP1 vs Mock a 24H (adim.)\n"
                   "Azul=reprimido por PEP1 | Rojo=inducido por PEP1")
    fc_extra = {"FC_"+tp: calc_fc(tp).to_dict() for tp in TIME_ORDER}
    export_cy(G4, "Cytoscape_Red4_FC_24H", extra_cols=fc_extra)
    del G4
del TOM_module; gc.collect()

print("  ✓ ETAPA 5 completada")


# =============================================================================
# RESUMEN FINAL
# =============================================================================
goi_gs   = float(GS.get(GOI, np.nan))
goi_mm   = float(MM.get(GOI, np.nan))
goi_kw_v = kIM.loc[kIM["gene_id"]==GOI,"kWithin"].values
goi_kw   = float(goi_kw_v[0]) if len(goi_kw_v) else np.nan
n_co  = len(cor_df[(cor_df["pearson_r"]> 0.5)&(cor_df["FDR"]<0.05)])
n_anti= len(cor_df[(cor_df["pearson_r"]<-0.5)&(cor_df["FDR"]<0.05)])

print("\n" + "=" * 65)
print(f"   ANÁLISIS COMPLETO — {GOI_NAME} ({GOI})")
print("=" * 65)
print(f"  Genes analizados     : {datExpr.shape[1]:,}")
print(f"  Muestras             : {datExpr.shape[0]}")
print(f"  Soft-threshold power : {soft_power}")
print(f"  Módulos detectados   : {(module_colors!='grey').nunique()}")
print(f"  Módulo de {GOI_NAME}        : {goi_mod}")
print(f"  Genes en el módulo   : {len(goi_genes)}")
print(f"  MM de {GOI_NAME}            : {goi_mm:.4f}")
print(f"  GS de {GOI_NAME} (PEP1)    : {goi_gs:.4f}")
print(f"  kWithin de {GOI_NAME}      : {goi_kw:.2f}")
print(f"  Co-expresados sig.   : {n_co}  (r>0.5, FDR<0.05)")
print(f"  Anti-correlados sig. : {n_anti}  (r<-0.5, FDR<0.05)")
print()
print("  FIGURAS (PDF + PNG):")
FIGS = [
    ("Fig01_SampleClustering",      "Dendrograma QC de muestras"),
    ("Fig02_SoftThreshold",         "Selección del soft-thresholding power (R² y conectividad)"),
    ("Fig03_ModuleDendrogram",      "Dendrograma de genes + barra de colores de módulos"),
    ("Fig04_ModuleTrait_Heatmap",   "Correlación módulo–PEP1/tiempo (heatmap)"),
    ("Fig05_GS_vs_MM",              "Gene Significance vs Module Membership (scatter)"),
    ("Fig06_Eigengene_Profile",     "Perfil temporal del eigengene Mock vs PEP1"),
    ("Fig07_Heatmap_Modulo",        "Heatmap top 50 genes del módulo (Z-score)"),
    ("Fig08_Volcano_Correlaciones", "Volcano: r vs -log10(FDR) todos los genes vs GOI"),
    ("Fig09_Barplot_Correlaciones", f"Barplot top {TOP_BAR} co-expresados + anti-correlados"),
    ("Fig10_Scatter_Top_Vecinos",   "Grilla de scatter: GOI vs sus top vecinos"),
    ("Fig11_kWithin_vs_MM",         "Conectividad intra-modular (kWithin) vs MM"),
    ("Fig12_Red1_Vecindad_GOI",     f"Red: vecindad directa de {GOI_NAME} (TOM≥{TOM_THRESH})"),
    ("Fig13_Red2_Modulo_Completo",  f"Red: módulo completo top {MAX_NODES} genes (color=GS)"),
    ("Fig14_Red3_Red_Modulos",      "Red: módulos como nodos (eigengenes como aristas)"),
    ("Fig15_Red4_FC_24H",           "Red: vecindad coloreada por Fold Change 24H"),
]
for fname, desc in FIGS:
    print(f"    {fname}.pdf/.png  —  {desc}")
print()
print("  TABLAS CSV:")
for t in ["Tabla_Correlaciones_RAP2-3.csv","Tabla_Modulo_RAP2-3.csv"]:
    print(f"    {t}")
print()
print("  ARCHIVOS CYTOSCAPE (.graphml + _edges.csv + _nodes.csv):")
for cy in ["Cytoscape_Red1_Vecindad_GOI",
           "Cytoscape_Red2_Modulo_Completo",
           "Cytoscape_Red4_FC_24H"]:
    print(f"    {cy}")
print()
print("  CHECKPOINTS (borrar para reiniciar):")
for cp in ["cp1_datExpr.pkl","cp2_network.pkl",
           "cp3_tom_module.pkl","cp4_scores.pkl","TOM_full.npy"]:
    print(f"    {cp}")
print("=" * 65)
