# -*- coding: utf-8 -*-
"""
生成峰形演化分析的所有科学可视化图表 (整合版 v5, 面向期刊投稿)
图号按论文最终阅读顺序连续编号 (Fig.1-8)：
  Fig.1 谱带偏度-载荷曲线                   [3.1]
  Fig.2 主峰FWHM-载荷迟滞回线               [3.1]
  Fig.3 D2/main-载荷曲线                    [3.1]
  Fig.4 峰位 vs 偏度灵敏度对比 (TO加载段)    [3.1]
  Fig.5 代表性拉曼光谱分解 (原始+基线+分峰)  [3.2]
  Fig.6 全载荷序列拟合质量瀑布图 (TO)        [3.2]
  Fig.7 加载前/完全卸载后峰形参数对比        [3.4]
  Fig.8 PCA得分图                          [3.5]

本版本相对上一版的修改：
  1) 所有图内不再绘制图题 (标题/子图标题一律移除，图注留待正文/图片说明处添加)。
  2) 原先含有并列子图的图 (Fig.1, Fig.2, Fig.3, Fig.6, Fig.7) 全部合并为单幅图，
     用不同颜色 (以及必要时的线型/图案) 区分原来分属不同子图的曲线/柱状。
"""
from pathlib import Path
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

BASE_DIR = Path(__file__).resolve().parent.parent
FIGDIR = BASE_DIR / 'figures'
OUTDIR = BASE_DIR / 'outputs'
DATA_ROOT = BASE_DIR / 'mds2-2281'
FIGDIR.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)

import sys
sys.path.insert(0, str(BASE_DIR / 'src'))
from peak_analysis import (load_raman_sequence, subtract_local_baseline,
                            normalize_area, fit_three_peaks, _als_baseline)

csv_path = OUTDIR / 'peak_shape_results.csv'
if not csv_path.exists():
    raise FileNotFoundError(f"找不到基础数据文件 {csv_path}，请先运行 run_analysis.py。")
df = pd.read_csv(csv_path)

# ---------------------------------------------------------------------------
# 颜色 / 标记 / 线型约定
#   - OBJ_COLORS: 用颜色区分物镜配置 (TO / CO)  —— 原先分属两个子图，现合并后用颜色区分
#   - SEG_STYLES: 用线型区分加载/卸载段
#   - MARKERS:    用标记形状进一步区分物镜配置 (辅助颜色，双重编码更易读)
# ---------------------------------------------------------------------------
OBJ_COLORS = {'TO': '#c0392b', 'CO': '#2980b9'}
SEG_STYLES = {'loading': '-', 'unloading': '--'}
MARKERS = {'TO': 'o', 'CO': 's'}
# 仍保留原 loading/unloading 配色，供未合并子图的图 (如 Fig.4) 使用
COLORS = {'loading': '#c0392b', 'unloading': '#2980b9'}


def style_ax(ax, xlabel, ylabel, title=None):
    """设置坐标轴标签与网格样式；不在图内绘制标题 (title 参数保留仅为兼容，不使用)。"""
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ===========================================================================
# Fig.1  Envelope skewness vs load  [Sec. 3.1]
#   颜色区分物镜配置，线型区分加载/卸载
# ===========================================================================
fig, ax = plt.subplots(figsize=(7.5, 5.2))
for obj in ['TO', 'CO']:
    for seg in ['loading', 'unloading']:
        sub = df[(df.objective == obj) & (df.segment == seg)].sort_values('load_mN')
        ax.plot(sub.load_mN, sub.env_skewness, marker=MARKERS[obj], color=OBJ_COLORS[obj],
                linestyle=SEG_STYLES[seg], label=f'{obj} \u2013 {seg.capitalize()}',
                linewidth=2, markersize=6, alpha=0.9)
ax.axhline(0, color='gray', lw=1.0, ls=':')
style_ax(ax, 'Indentation load (mN)', 'Envelope skewness (dimensionless)')
ax.legend(frameon=False, ncol=1, fontsize=10, handlelength=4)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig1_skewness_vs_load.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.2  Main-band FWHM hysteresis loop vs load  [Sec. 3.1]
#   颜色区分物镜配置，线型区分加载/卸载
# ===========================================================================
fig, ax = plt.subplots(figsize=(7.5, 5.2))
for obj in ['TO', 'CO']:
    for seg in ['loading', 'unloading']:
        sub = df[(df.objective == obj) & (df.segment == seg)].sort_values('load_mN')
        ax.plot(sub.load_mN, sub.main_fwhm, marker=MARKERS[obj], color=OBJ_COLORS[obj],
                linestyle=SEG_STYLES[seg], label=f'{obj} \u2013 {seg.capitalize()}',
                linewidth=2, markersize=6, alpha=0.9)
style_ax(ax, 'Indentation load (mN)', 'Main-band FWHM (cm$^{-1}$)')
ax.legend(frameon=False, ncol=1, fontsize=10, handlelength=4)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig2_fwhm_hysteresis.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.3  D2/main area ratio vs load  [Sec. 3.1]
#   原为 TO / CO 两个并列子图，现合并为单幅：颜色区分物镜配置，线型区分加载/卸载
# ===========================================================================
fig, ax = plt.subplots(figsize=(7.5, 6.6))
for obj in ['TO', 'CO']:
    for seg in ['loading', 'unloading']:
        sub = df[(df.objective == obj) & (df.segment == seg)].sort_values('load_mN')
        ax.plot(sub.load_mN, sub.ratio_D2_main, marker=MARKERS[obj], color=OBJ_COLORS[obj],
                linestyle=SEG_STYLES[seg], label=f'{obj} \u2013 {seg.capitalize()}',
                linewidth=2, markersize=6, alpha=0.9)
style_ax(ax, 'Indentation load (mN)', 'D2 / Main-band area ratio')
# ax.legend(frameon=False, ncol=2, fontsize=12, loc='center right', bbox_to_anchor=(1.0, 0.5))
ax.legend(frameon=False, ncol=1, fontsize=10, loc='upper right', handlelength=4)

fig.tight_layout()
fig.savefig(FIGDIR / 'fig3_D2main_ratio.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.4  Peak position vs skewness sensitivity comparison  [Sec. 3.1]
#   (单一双 y 轴图，无并列子图，此处仅去除图题)
# ===========================================================================
fig, ax1 = plt.subplots(figsize=(7.5, 5.2))
sub = df[(df.objective == 'TO') & (df.segment == 'loading')].sort_values('load_mN')
l1, = ax1.plot(sub.load_mN, sub.env_centroid, 'o-', color='#2980b9', linewidth=2,
               markersize=6, label='Envelope centroid position')
ax1.set_xlabel('Indentation load (mN)', fontsize=12)
ax1.set_ylabel('Envelope centroid position (cm$^{-1}$)', color='#2980b9', fontsize=12)
ax1.tick_params(axis='y', labelcolor='#2980b9')
ax1.grid(alpha=0.3, linestyle='--')
ax1.spines['top'].set_visible(False)
ax2 = ax1.twinx()
l2, = ax2.plot(sub.load_mN, sub.env_skewness, 's--', color="#c0392b", linewidth=2,
               markersize=6, label='envelope skewness')
ax2.set_ylabel('Envelope skewness (dimensionless)', color='#c0392b', fontsize=12)
ax2.tick_params(axis='y', labelcolor='#c0392b')
ax2.spines['top'].set_visible(False)
fig.legend([l1, l2], ['envelope centroid position', 'envelope skewness'],
           loc='center right', bbox_to_anchor=(0.90, 0.55), frameon=True, framealpha=0.8, handlelength=4)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig4_position_vs_skewness.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.5  Representative Raman spectrum decomposition  [Sec. 3.2]
#   (a) raw spectrum + ALS baseline
#   (b) baseline-subtracted spectrum + 3-peak Pseudo-Voigt deconvolution
#   代表性光谱: TO物镜, 加载段, 100 mN (信噪比适中, 三分量分离清晰, redchi低)
#   这两个子图内容不同 (原始谱 vs 分峰谱)，不属于"同类曲线拆分在不同子图"的情形，
#   故保留为两个子图，仅去除图题文字，改用简洁的 (a)/(b) 面板标签。
# ===========================================================================
REP_RAMAN = DATA_ROOT / 'Dataset C/TO_load_curve_raw_data.csv'
REP_LOAD = 100.0

wn_full, forces_full, I_full = load_raman_sequence(REP_RAMAN)
idx_rep = forces_full.index(REP_LOAD)

ctx_mask = (wn_full >= 150) & (wn_full <= 900)
wn_ctx = wn_full[ctx_mask]
y_ctx_raw = I_full[ctx_mask, idx_rep]
baseline_ctx = _als_baseline(y_ctx_raw, lam=1e7, p=0.001, niter=15)

wn_win, y_win = subtract_local_baseline(wn_full, I_full[:, idx_rep])
y_norm = normalize_area(wn_win, y_win)
fit_rep = fit_three_peaks(wn_win, y_norm, keep_result_obj=True)
result_rep = fit_rep['_result']

comp_colors = {'main': '#2980b9', 'D1': '#27ae60', 'D2': '#e74c3c'}
comp_labels = {'main': 'Main-band (Si\u2013O\u2013Si bending)', 'D1': 'D1 (4-membered ring)',
               'D2': 'D2 (3-membered ring)'}

fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.5, 8.0))

axA.plot(wn_ctx, y_ctx_raw, color='#2c3e50', lw=1.3, label='Experimental spectrum (raw)')
axA.plot(wn_ctx, baseline_ctx, color='#f39c12', lw=1.8, ls='--', label='ALS baseline')
axA.axvspan(350, 700, color='gray', alpha=0.08, label='Fitting window (350\u2013700 cm$^{-1}$)')
style_ax(axA, 'Raman shift (cm$^{-1}$)', 'Intensity (counts)')
axA.text(0.02, 0.95, '(a)', transform=axA.transAxes, fontsize=12, fontweight='bold', va='top')
axA.legend(frameon=False, fontsize=10, loc='upper right')

axB.plot(wn_win, y_norm, 'o', color='#2c3e50', ms=3.5, alpha=0.6,
         label='Baseline-subtracted spectrum')
comps = result_rep.eval_components(x=wn_win)
total_fit = result_rep.eval(x=wn_win)
for name in ['main', 'D1', 'D2']:
    y_comp = comps[f'{name}_']
    axB.fill_between(wn_win, 0, y_comp, color=comp_colors[name], alpha=0.35,
                      label=comp_labels[name])
    axB.plot(wn_win, y_comp, color=comp_colors[name], lw=1.2)
axB.plot(wn_win, total_fit, color='black', lw=2.0, label='Fitted envelope (sum)')
style_ax(axB, 'Raman shift (cm$^{-1}$)', 'Normalized intensity')
axB.text(0.02, 0.95, '(b)', transform=axB.transAxes, fontsize=12, fontweight='bold', va='top')
axB.legend(frameon=False, fontsize=10, loc='upper right', bbox_to_anchor=(1.0, 1.0))

fig.tight_layout()
fig.savefig(FIGDIR / 'fig5_spectrum_decomposition.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"[Fig.5] Representative spectrum: TO configuration, loading, {REP_LOAD:.0f} mN, "
      f"redchi={fit_rep['redchi']:.2e}")



# ===========================================================================
# Fig.6  Pseudo-Voigt fitting results for the TO loading-unloading sequence
#        (waterfall)  [Sec. 3.2, \label{fig:waterfall}]
#
#   Caption (paper):
#     "Pseudo-Voigt fitting results for the TO loading-unloading sequence.
#      Experimental spectra and fitted envelopes are shown for the seven
#      loading steps and seven unloading steps."
#
#   Body text (paper):
#     "The fitting behavior was also examined over the complete TO
#      loading-unloading sequence. Figure~\ref{fig:waterfall} shows the
#      measured spectra and fitted envelopes at each load level. The fitted
#      curves reproduce the experimental spectra throughout the sequence,
#      including the broadened and asymmetric profiles observed at the
#      highest load."
#
#   Implementation notes:
#     - Each spectrum in the TO load/unload sequence is baseline-subtracted,
#       area-normalized, and fit with the same 3-peak Pseudo-Voigt model
#       used elsewhere (with sequential warm-starting for stability).
#     - Left panel = loading segment, right panel = unloading segment, each
#       stacked vertically ("waterfall") in load order with experimental
#       points (markers) overlaid by the fitted envelope (line).
#     - The stacking offset is derived from the data itself (not a fixed
#       constant) so the panels remain legible regardless of absolute
#       normalized-intensity scale.
#     - A per-spectrum fit-quality table (redchi + RMSE between the
#       experimental and fitted envelope) is exported alongside the figure
#       so the "fitted curves reproduce the experimental spectra throughout
#       the sequence" claim is backed by numbers, not just visual inspection.
# ===========================================================================
EXPECTED_N_STEPS = 7  # per the figure caption: seven loading + seven unloading steps


def _fit_series(raman_path, expected_n=EXPECTED_N_STEPS):
    """Baseline-subtract, normalize, and 3-peak-fit every spectrum in a
    load or unload sequence file, returning one record per step in load
    order together with a simple experimental-vs-fit residual metric."""
    wn_s, forces_s, I_s = load_raman_sequence(raman_path)

    if expected_n is not None and len(forces_s) != expected_n:
        warnings.warn(
            f"{Path(raman_path).name}: found {len(forces_s)} load steps, "
            f"expected {expected_n} (as stated in the Fig.6 caption). "
            f"Proceeding with the steps actually present in the file."
        )

    order = np.argsort(forces_s)  # ensure strictly increasing load order for the waterfall stack
    out = []
    warm = None
    for i in order:
        f = forces_s[i]
        wn_win_s, y_win_s = subtract_local_baseline(wn_s, I_s[:, i])
        y_norm_s = normalize_area(wn_win_s, y_win_s)
        fit_s = fit_three_peaks(wn_win_s, y_norm_s, warm_start=warm, keep_result_obj=True)
        warm = fit_s if (fit_s['redchi'] < 0.01) else None
        total_s = fit_s['_result'].eval(x=wn_win_s)
        rmse_s = float(np.sqrt(np.mean((y_norm_s - total_s) ** 2)))
        out.append(dict(load=f, wn=wn_win_s, y=y_norm_s, fit=total_s,
                         redchi=fit_s['redchi'], rmse=rmse_s))
    return out


def _stack_offset(series):
    """Offset between successive traces in the waterfall, sized so stacked
    spectra don't overlap: a fraction of the typical peak height."""
    typical_peak = np.median([np.max(rec['y']) for rec in series])
    return 0.55 * typical_peak


def plot_fig6_waterfall(series_load, series_unload, out_path):
    """Render the two-panel (loading / unloading) waterfall of experimental
    spectra + fitted Pseudo-Voigt envelopes for the TO sequence."""
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.5, 6.2))

    for ax, series, panel_tag, seg_label in [(axL, series_load, '(a)', 'Loading'),
                                              (axR, series_unload, '(b)', 'Unloading')]:
        offset = _stack_offset(series)
        for k, rec in enumerate(series):
            off = k * offset
            ax.plot(rec['wn'], rec['y'] + off, 'o', color='#7f8c8d', ms=2.2, alpha=0.55,
                    markeredgewidth=0)
            ax.plot(rec['wn'], rec['fit'] + off, color='#c0392b', lw=1.4)
            ax.text(705, rec['fit'][-1] + off, f"{rec['load']:.0f} mN",
                    fontsize=8.5, va='center', color='#2c3e50', clip_on=False)

        style_ax(ax, 'Raman shift (cm$^{-1}$)', 'Normalized intensity (stacked, a.u.)')
        ax.text(0.02, 0.99, panel_tag, transform=ax.transAxes, fontsize=12,
                fontweight='bold', va='top')
        ax.set_xlim(345, 775)
        ax.set_yticks([])

    legend_elems = [Line2D([0], [0], marker='o', color='none', markerfacecolor='#7f8c8d',
                            markersize=6, label='Experimental spectrum'),
                    Line2D([0], [0], color='#c0392b', lw=1.8, label='3-peak Pseudo-Voigt fit (sum)')]
    axL.legend(handles=legend_elems, frameon=False, fontsize=9, loc='upper center',
               bbox_to_anchor=(0.6, 1.10), ncol=1)
    axR.legend(handles=legend_elems, frameon=False, fontsize=9, loc='upper center',
               bbox_to_anchor=(0.6, 1.10), ncol=1)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


series_load = _fit_series(DATA_ROOT / 'Dataset C/TO_load_curve_raw_data.csv')
series_unload = _fit_series(DATA_ROOT / 'Dataset C/TO_unload_curve_raw_data.csv')

plot_fig6_waterfall(series_load, series_unload, FIGDIR / 'fig6_fit_quality_waterfall.png')

# Per-step fit-quality table (backs up the "fitted curves reproduce the
# experimental spectra throughout the sequence" claim quantitatively).
fitq_rows = []
for seg_label, series in [('loading', series_load), ('unloading', series_unload)]:
    for rec in series:
        fitq_rows.append(dict(segment=seg_label, load_mN=rec['load'],
                               redchi=rec['redchi'], rmse=rec['rmse']))
fitq_df = pd.DataFrame(fitq_rows)
fitq_df.to_csv(OUTDIR / 'TO_sequence_fit_quality.csv', index=False)

mean_redchi_load = np.mean([r['redchi'] for r in series_load])
mean_redchi_unload = np.mean([r['redchi'] for r in series_unload])
highest_load_rmse = series_load[-1]['rmse'] if series_load else np.nan
print(f"[Fig.6] steps -- loading: {len(series_load)}, unloading: {len(series_unload)} "
      f"(caption states seven each)")
print(f"[Fig.6] Mean redchi -- loading: {mean_redchi_load:.2e}, unloading: {mean_redchi_unload:.2e}")
print(f"[Fig.6] Highest-load step RMSE (experimental vs. fitted envelope): {highest_load_rmse:.2e}")

# ===========================================================================
# Fig.7  Pristine vs residual peak-shape comparison  [Sec. 3.4]
#   原为三个并列子图 (main-FWHM / skewness / D2-main ratio)，现合并为单幅分组柱状图：
#   每个指标各自做 min-max 归一化以便共享同一纵轴，用颜色区分指标 (原来的三个子图)，
#   用填充图案 (无填充=pristine，斜线=residual) 区分加载前/完全卸载后。
# ===========================================================================
metrics = [('main_fwhm', 'Main-FWHM (cm$^{-1}$)'),
           ('env_skewness', 'Envelope skewness'),
           ('ratio_D2_main', 'D2/main area ratio')]
metric_colors = {'main_fwhm': '#2980b9', 'env_skewness': '#27ae60', 'ratio_D2_main': '#e74c3c'}


def get_zero_load_val(obj_name, seg_name, col_name):
    sub_df = df[(df.objective == obj_name) & (df.segment == seg_name)]
    if sub_df.empty:
        return np.nan
    idx = np.abs(sub_df.load_mN).idxmin()
    return sub_df.loc[idx, col_name]


# 收集原始数值，并按指标做 min-max 归一化到 [0, 1]，使不同量纲的指标可以共享一个纵轴
raw_vals = {}   # (metric, obj, cond) -> raw value
for col, _ in metrics:
    vals = []
    for obj_cfg in ['TO', 'CO']:
        pristine = get_zero_load_val(obj_cfg, 'loading', col)
        residual = get_zero_load_val(obj_cfg, 'unloading', col)
        raw_vals[(col, obj_cfg, 'pristine')] = pristine
        raw_vals[(col, obj_cfg, 'residual')] = residual
        vals += [pristine, residual]
    vmin, vmax = np.nanmin(vals), np.nanmax(vals)
    span = (vmax - vmin) if (vmax - vmin) > 0 else 1.0
    for obj_cfg in ['TO', 'CO']:
        for cond in ['pristine', 'residual']:
            raw_vals[(col, obj_cfg, cond, 'norm')] = (raw_vals[(col, obj_cfg, cond)] - vmin) / span

fig, ax = plt.subplots(figsize=(7.5, 6.2))
x = np.arange(2)  # TO, CO
bar_w = 0.12
cond_hatch = {'pristine': '', 'residual': '//'}
cond_alpha = {'pristine': 0.55, 'residual': 0.95}

for m_idx, (col, metric_label) in enumerate(metrics):
    for c_idx, cond in enumerate(['pristine', 'residual']):
        offset = (m_idx * 2 + c_idx - 2.5) * bar_w
        for i, obj_cfg in enumerate(['TO', 'CO']):
            norm_val = raw_vals[(col, obj_cfg, cond, 'norm')]
            raw_val = raw_vals[(col, obj_cfg, cond)]
            bar = ax.bar(x[i] + offset, norm_val, bar_w, color=metric_colors[col],
                         alpha=cond_alpha[cond], hatch=cond_hatch[cond],
                         edgecolor='black', linewidth=0.5,
                         label=(metric_label if (i == 0 and c_idx == 0) else None))
            ax.text(x[i] + offset, norm_val + 0.02, f'{raw_val:.2g}', ha='center', va='bottom',
                    fontsize=12, rotation=90)

ax.set_xticks(x)
ax.set_xticklabels(['TO', 'CO'])
style_ax(ax, 'Objective configuration', 'Min-max normalized value (per metric)')
ax.set_ylim(0, 1.25)

metric_legend = [Line2D([0], [0], marker='s', color='none', markerfacecolor=metric_colors[c],
                         markersize=10, label=lbl) for c, lbl in metrics]
cond_legend = [
    plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', hatch='', label='Pristine (before loading)'),
    plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', hatch='//', label='Residual (after full unload)'),
]
leg1 = ax.legend(handles=metric_legend, frameon=False, fontsize=10, loc='upper left', bbox_to_anchor=(0.02, 1.08))
ax.add_artist(leg1)
ax.legend(handles=cond_legend, frameon=False, fontsize=10, loc='upper right', bbox_to_anchor=(1.05, 1.08))

fig.tight_layout()
fig.savefig(FIGDIR / 'fig7_pristine_vs_residual.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.8  PCA of band-shape feature matrix  [Sec. 3.5]
#   (单一散点图，无并列子图，此处仅去除图题)
# ===========================================================================
feature_cols = ['env_sigma', 'env_skewness', 'env_kurtosis', 'env_asymmetry',
                 'main_fwhm', 'main_area', 'D1_fwhm', 'D1_area',
                 'D2_fwhm', 'D2_area', 'ratio_D1_main', 'ratio_D2_main']
pca_df = df.dropna(subset=feature_cols).copy()
X = pca_df[feature_cols].values
Xs = StandardScaler().fit_transform(X)
pca = PCA(n_components=2)
scores = pca.fit_transform(Xs)
pca_df['PC1'] = scores[:, 0]
pca_df['PC2'] = scores[:, 1]

fig, ax = plt.subplots(figsize=(7.5, 5.2))
for obj in ['TO', 'CO']:
    for seg in ['loading', 'unloading']:
        sub2 = pca_df[(pca_df.objective == obj) & (pca_df.segment == seg)]
        sc = ax.scatter(sub2.PC1, sub2.PC2, c=sub2.load_mN, cmap='viridis',
                         marker=MARKERS[obj], s=80, alpha=0.9,
                         edgecolors=('none' if seg == 'loading' else '#e74c3c'),
                         linewidths=1.5, label=f'{obj} ({seg})')
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label('Indentation load (mN)', fontsize=12)
style_ax(ax, f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',
         f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
fig.tight_layout()
fig.savefig(FIGDIR / 'fig8_pca_scores.png', dpi=300, bbox_inches='tight')
plt.close(fig)

loadings = pd.DataFrame(pca.components_.T, index=feature_cols,
                         columns=['PC1_loading', 'PC2_loading'])
loadings.to_csv(OUTDIR / 'pca_loadings.csv')
pca_df.to_csv(OUTDIR / 'pca_scores.csv', index=False)

# ===========================================================================
# 统计检验: 各配置内 Load 与关键峰形参数的 Pearson 相关 (r, p)
# ===========================================================================
stat_rows = []
for obj in ['TO', 'CO']:
    sub3 = df[(df.objective == obj) & (df.segment == 'loading')].sort_values('load_mN')
    # print(sub3[['load_mN', 'env_skewness', 'env_centroid', 'main_fwhm', 'ratio_D2_main', 'env_kurtosis', 'env_asymmetry']].to_string(index=False))
    # pass
    for col in ['env_skewness', 'env_centroid', 'main_fwhm', 'ratio_D2_main', 'env_kurtosis', 'env_asymmetry']:
        r, p = stats.pearsonr(sub3['load_mN'], sub3[col])
        stat_rows.append(dict(objective=obj, metric=col, n=len(sub3), pearson_r=r, p_value=p))
stat_df = pd.DataFrame(stat_rows)
stat_df.to_csv(OUTDIR / 'correlation_stats.csv', index=False)
print(stat_df.to_string(index=False))

print('\nPCA explained variance ratio:', pca.explained_variance_ratio_)
print(loadings)
print(f'\n所有图表 (Fig.1-8, 300 dpi) 已保存至: {FIGDIR}')