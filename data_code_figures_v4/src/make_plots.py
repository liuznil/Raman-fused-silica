# -*- coding: utf-8 -*-
"""
生成峰形演化分析的所有科学可视化图表 (整合版 v3, 面向期刊投稿)
图号按论文最终阅读顺序连续编号 (Fig.1-8)：
  Fig.1 主峰FWHM-载荷迟滞回线               [3.1]
  Fig.2 谱带偏度-载荷曲线                   [3.1]
  Fig.3 D2/主峰面积比-载荷曲线              [3.1]
  Fig.4 峰位 vs 偏度灵敏度对比 (TO加载段)    [3.1]
  Fig.5 代表性拉曼光谱分解 (原始+基线+分峰)   [3.2]
  Fig.6 全载荷序列拟合质量瀑布图 (TO)        [3.2]
  Fig.7 加载前/完全卸载后峰形参数对比        [3.4]
  Fig.8 PCA得分图                          [3.5]
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

COLORS = {'loading': '#c0392b', 'unloading': '#2980b9'}
MARKERS = {'TO': 'o', 'CO': 's'}


def style_ax(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.grid(alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ===========================================================================
# Fig.1  Main band FWHM hysteresis loop vs load  [Sec. 3.1]
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
for ax, obj in zip(axes, ['TO', 'CO']):
    for seg in ['loading', 'unloading']:
        sub = df[(df.objective == obj) & (df.segment == seg)].sort_values('load_mN')
        ax.plot(sub.load_mN, sub.main_fwhm, marker=MARKERS[obj], color=COLORS[seg],
                label=seg.capitalize(), linewidth=2, markersize=6, alpha=0.9)
    style_ax(ax, 'Indentation load (mN)', 'Main band FWHM (cm$^{-1}$)', f'{obj} objective')
    ax.legend(frameon=False)
fig.suptitle('Fig. 1  Hysteresis loop of the Main band FWHM versus load',
             fontsize=13, fontweight='bold', y=1.03)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig1_fwhm_hysteresis.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.2  Envelope skewness vs load  [Sec. 3.1]
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
for ax, obj in zip(axes, ['TO', 'CO']):
    for seg in ['loading', 'unloading']:
        sub = df[(df.objective == obj) & (df.segment == seg)].sort_values('load_mN')
        ax.plot(sub.load_mN, sub.env_skewness, marker=MARKERS[obj], color=COLORS[seg],
                label=seg.capitalize(), linewidth=2, markersize=6, alpha=0.9)
    style_ax(ax, 'Indentation load (mN)', 'Envelope skewness (dimensionless)', f'{obj} objective')
    ax.axhline(0, color='gray', lw=1.0, ls='--')
    ax.legend(frameon=False)
fig.suptitle('Fig. 2  Overall Raman band skewness versus load', fontsize=13, fontweight='bold', y=1.03)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig2_skewness_vs_load.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.3  D2/Main area ratio vs load  [Sec. 3.1]
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
for ax, obj in zip(axes, ['TO', 'CO']):
    for seg in ['loading', 'unloading']:
        sub = df[(df.objective == obj) & (df.segment == seg)].sort_values('load_mN')
        ax.plot(sub.load_mN, sub.ratio_D2_main, marker=MARKERS[obj], color=COLORS[seg],
                label=seg.capitalize(), linewidth=2, markersize=6, alpha=0.9)
    style_ax(ax, 'Indentation load (mN)', 'D2 / Main band area ratio', f'{obj} objective')
    ax.legend(frameon=False)
fig.suptitle('Fig. 3  D2/Main band area ratio versus load \u2014 densification/ring structure indicator',
             fontsize=13, fontweight='bold', y=1.03)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig3_D2main_ratio.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.4  Peak position vs skewness sensitivity comparison  [Sec. 3.1]
# ===========================================================================
fig, ax1 = plt.subplots(figsize=(7.5, 4.8))
sub = df[(df.objective == 'TO') & (df.segment == 'loading')].sort_values('load_mN')
l1, = ax1.plot(sub.load_mN, sub.env_centroid, 'o-', color='#8e44ad', linewidth=2,
               markersize=6, label='Peak centroid position')
ax1.set_xlabel('Indentation load (mN)', fontsize=11)
ax1.set_ylabel('Peak centroid position (cm$^{-1}$)', color='#8e44ad', fontsize=11)
ax1.tick_params(axis='y', labelcolor='#8e44ad')
ax1.grid(alpha=0.3, linestyle='--')
ax2 = ax1.twinx()
l2, = ax2.plot(sub.load_mN, sub.env_skewness, 's--', color='#16a085', linewidth=2,
               markersize=6, label='Envelope skewness')
ax2.set_ylabel('Envelope skewness (dimensionless)', color='#16a085', fontsize=11)
ax2.tick_params(axis='y', labelcolor='#16a085')
ax2.spines['top'].set_visible(False)
ax1.set_title('Fig. 4  TO loading: band position vs. skewness \u2014 which metric responds first?',
              fontsize=12, fontweight='bold', pad=12)
fig.legend([l1, l2], ['band centroid position', 'Envelope skewness'],
           loc='lower right', bbox_to_anchor=(0.88, 0.15), frameon=True, framealpha=0.8)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig4_position_vs_skewness.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.5  Representative Raman spectrum decomposition  [Sec. 3.2, NEW]
#   (a) raw spectrum + ALS baseline
#   (b) baseline-subtracted spectrum + 3-peak Pseudo-Voigt deconvolution
#   代表性光谱: TO物镜, 加载段, 100 mN (信噪比适中, 三分量分离清晰, redchi低)
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
comp_labels = {'main': 'Main band (Si\u2013O\u2013Si bending)', 'D1': 'D1 (4-membered ring)',
               'D2': 'D2 (3-membered ring)'}

fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.5, 9.0))

axA.plot(wn_ctx, y_ctx_raw, color='#2c3e50', lw=1.3, label='Experimental spectrum (raw)')
axA.plot(wn_ctx, baseline_ctx, color='#f39c12', lw=1.8, ls='--', label='ALS baseline')
axA.axvspan(350, 700, color='gray', alpha=0.08, label='Fitting window (350\u2013700 cm$^{-1}$)')
style_ax(axA, 'Raman shift (cm$^{-1}$)', 'Intensity (counts)',
         '(a) Raw spectrum and ALS baseline (TO objective, 100 mN, loading)')
axA.legend(frameon=False, fontsize=9, loc='upper right')

axB.plot(wn_win, y_norm, 'o', color='#2c3e50', ms=3.5, alpha=0.6,
         label='Experimental spectrum (baseline-subtracted)')
comps = result_rep.eval_components(x=wn_win)
total_fit = result_rep.eval(x=wn_win)
for name in ['main', 'D1', 'D2']:
    y_comp = comps[f'{name}_']
    axB.fill_between(wn_win, 0, y_comp, color=comp_colors[name], alpha=0.35,
                      label=comp_labels[name])
    axB.plot(wn_win, y_comp, color=comp_colors[name], lw=1.2)
axB.plot(wn_win, total_fit, color='black', lw=2.0, label='Fitted envelope (sum)')
style_ax(axB, 'Raman shift (cm$^{-1}$)', 'Normalized intensity',
         '(b) Baseline-subtracted spectrum and 3-peak Pseudo-Voigt decomposition')
axB.legend(frameon=False, fontsize=8.5, loc='upper right')

fig.suptitle('Fig. 5  Representative Raman spectrum decomposition', fontsize=13,
             fontweight='bold', y=1.0)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig5_spectrum_decomposition.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"[Fig.5] Representative spectrum: TO objective, loading, {REP_LOAD:.0f} mN, "
      f"redchi={fit_rep['redchi']:.2e}")

# ===========================================================================
# Fig.6  Fitting quality across the full load-unload sequence (waterfall)  [Sec. 3.2]
#   展示: TO配置加载段(7条)与卸载段(6条)全部原始光谱(点) + 三峰拟合总包络(线)
#   叠加显示，验证拟合质量在整个载荷序列上保持稳定，而非仅代表性单谱良好。
# ===========================================================================
def _fit_series(raman_path, segment_forces=None):
    wn_s, forces_s, I_s = load_raman_sequence(raman_path)
    out = []
    warm = None
    for i, f in enumerate(forces_s):
        wn_win_s, y_win_s = subtract_local_baseline(wn_s, I_s[:, i])
        y_norm_s = normalize_area(wn_win_s, y_win_s)
        fit_s = fit_three_peaks(wn_win_s, y_norm_s, warm_start=warm, keep_result_obj=True)
        warm = fit_s if (fit_s['redchi'] < 0.01) else None
        total_s = fit_s['_result'].eval(x=wn_win_s)
        out.append(dict(load=f, wn=wn_win_s, y=y_norm_s, fit=total_s, redchi=fit_s['redchi']))
    return out


series_load = _fit_series(DATA_ROOT / 'Dataset C/TO_load_curve_raw_data.csv')
series_unload = _fit_series(DATA_ROOT / 'Dataset C/TO_unload_curve_raw_data.csv')

fig, (axL, axU) = plt.subplots(1, 2, figsize=(12, 8), sharey=False)
OFFSET = 0.008

for ax, series, seg_label in [(axL, series_load, 'Loading'), (axU, series_unload, 'Unloading')]:
    for k, rec in enumerate(series):
        off = k * OFFSET
        ax.plot(rec['wn'], rec['y'] + off, 'o', color='#7f8c8d', ms=2.2, alpha=0.55,
                markeredgewidth=0)
        ax.plot(rec['wn'], rec['fit'] + off, color='#c0392b', lw=1.4)
        ax.text(705, rec['fit'][-1] + off, f"{rec['load']:.0f} mN",
                fontsize=8.5, va='center', color='#2c3e50', clip_on=False)
    style_ax(ax, 'Raman shift (cm$^{-1}$)', 'Normalized intensity (stacked, a.u.)',
              f'{seg_label} segment (TO objective)')
    ax.set_xlim(345, 775)
    ax.set_yticks([])

# 图例 (仅需一份，放在左图)
from matplotlib.lines import Line2D
legend_elems = [Line2D([0], [0], marker='o', color='none', markerfacecolor='#7f8c8d',
                        markersize=6, label='Experimental data'),
                Line2D([0], [0], color='#c0392b', lw=1.8, label='3-peak Pseudo-Voigt fit (sum)')]
axL.legend(handles=legend_elems, frameon=False, fontsize=9, loc='upper left')

fig.suptitle('Fig. 6  Fitting quality across the full load-unload sequence (TO objective)',
             fontsize=13, fontweight='bold', y=1.0)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig6_fit_quality_waterfall.png', dpi=300, bbox_inches='tight')
plt.close(fig)

mean_redchi_load = np.mean([r['redchi'] for r in series_load])
mean_redchi_unload = np.mean([r['redchi'] for r in series_unload])
print(f"[Fig.6] Mean redchi -- loading: {mean_redchi_load:.2e}, unloading: {mean_redchi_unload:.2e}")

# ===========================================================================
# Fig.7  Pristine vs residual peak-shape comparison  [Sec. 3.4]
# ===========================================================================
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
metrics = [('main_fwhm', 'Main FWHM (cm$^{-1}$)'),
           ('env_skewness', 'Envelope skewness'),
           ('ratio_D2_main', 'D2/Main area ratio')]


def get_zero_load_val(obj_name, seg_name, col_name):
    sub_df = df[(df.objective == obj_name) & (df.segment == seg_name)]
    if sub_df.empty:
        return np.nan
    idx = np.abs(sub_df.load_mN).idxmin()
    return sub_df.loc[idx, col_name]


for ax, (col, metric_label) in zip(axes, metrics):
    x = np.arange(2)
    width = 0.35
    for i, obj_cfg in enumerate(['TO', 'CO']):
        pristine = get_zero_load_val(obj_cfg, 'loading', col)
        residual = get_zero_load_val(obj_cfg, 'unloading', col)
        ax.bar(x[i] - width / 2, pristine, width, color='#95a5a6',
               label='Pristine (before loading)' if i == 0 else None)
        ax.bar(x[i] + width / 2, residual, width, color='#e67e22',
               label='Residual (after full unload)' if i == 0 else None)
    ax.set_xticks(x)
    ax.set_xticklabels(['TO', 'CO'])
    style_ax(ax, 'Objective configuration', metric_label, metric_label)
    if ax is axes[0]:
        ax.legend(frameon=False, fontsize=9, loc='best')
fig.suptitle('Fig. 7  band shape parameters before loading vs. after full unloading (at 0 mN)',
             fontsize=13, fontweight='bold', y=1.03)
fig.tight_layout()
fig.savefig(FIGDIR / 'fig7_pristine_vs_residual.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# ===========================================================================
# Fig.8  PCA of band-shape feature matrix  [Sec. 3.5]
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

fig, ax = plt.subplots(figsize=(6.8, 5.4))
for obj in ['TO', 'CO']:
    for seg in ['loading', 'unloading']:
        sub2 = pca_df[(pca_df.objective == obj) & (pca_df.segment == seg)]
        sc = ax.scatter(sub2.PC1, sub2.PC2, c=sub2.load_mN, cmap='viridis',
                         marker=MARKERS[obj], s=80, alpha=0.9,
                         edgecolors=('none' if seg == 'loading' else '#e74c3c'),
                         linewidths=1.5, label=f'{obj} ({seg})')
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label('Indentation load (mN)', fontsize=11)
style_ax(ax, f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',
         f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)',
         'Fig. 8  PCA of band shape features\n(circles: TO, squares: CO; red edges: unloading)')
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
    for col in ['env_skewness', 'env_centroid', 'main_fwhm', 'ratio_D2_main', 'env_kurtosis']:
        r, p = stats.pearsonr(sub3['load_mN'], sub3[col])
        stat_rows.append(dict(objective=obj, metric=col, n=len(sub3), pearson_r=r, p_value=p))
stat_df = pd.DataFrame(stat_rows)
stat_df.to_csv(OUTDIR / 'correlation_stats.csv', index=False)
print(stat_df.to_string(index=False))

print('\nPCA explained variance ratio:', pca.explained_variance_ratio_)
print(loadings)
print(f'\n所有图表 (Fig.1-8, 300 dpi) 已保存至: {FIGDIR}')
