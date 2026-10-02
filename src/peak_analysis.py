"""
NIST 熔融石英压痕原位拉曼数据 - 峰形演化分析核心模块 (整合优化版)
======================================================

数据来源: NIST mds2-2281 "In-situ Raman spectra of indented fused silica"
          (Gerbig & Michaels, 2020) https://doi.org/10.18434/M32281
"""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
from lmfit.models import PseudoVoigtModel

# 兼容 NumPy 1.x 和 2.0+ 的梯形积分
try:
    from scipy.integrate import trapezoid
except ImportError:
    trapezoid = getattr(np, 'trapezoid', getattr(np, 'trapz'))

# ----------------------------------------------------------------------
# 1. 数据读取与映射
# ----------------------------------------------------------------------

def load_raman_sequence(csv_path):
    """
    读取 Dataset C / F 格式的 load 或 unload 序列光谱文件。
    某些文件 (如 CO_unload) 含非力值的参考列 (如 'Pre-test')，自动跳过。
    返回: wavenumber (n,), forces (list[float] mN), intensities (n, m) ndarray
    """
    df = pd.read_csv(csv_path)
    wavenumber = df.iloc[:, 0].values.astype(float)

    forces, cols = [], []
    for c in df.columns[1:]:
        try:
            val = float(str(c).replace(' mN', '').strip())
            forces.append(val)
            cols.append(c)
        except ValueError:
            continue  # 跳过参考列如 'Pre-test'

    intensities = df[cols].values.astype(float)
    return wavenumber, forces, intensities


def load_indentation_curve(csv_path):
    """读取 Dataset D / G 格式的 Force-Displacement-Time 曲线"""
    df = pd.read_csv(csv_path)
    df.columns = ['time_s', 'force_mN', 'displacement_nm']
    return df


def force_to_depth_map(indent_df, target_forces, segment='loading', force_tol=5.0):
    """
    把目标力值 (mN) 映射为对应的位移 (nm)。
    force_tol: 匹配容差 (mN)，若最近点与目标力值偏差超过该值则发出警告
               (用于提示压痕曲线采样率是否足以覆盖所有目标载荷点)。
    """
    f = indent_df['force_mN'].values
    d = indent_df['displacement_nm'].values
    imax = np.argmax(f)

    if segment == 'loading':
        seg_f, seg_d = f[:imax + 1], d[:imax + 1]
    else:
        seg_f, seg_d = f[imax:], d[imax:]

    depths = []
    for tf in target_forces:
        if tf == 0 and segment == 'unloading':
            depths.append(seg_d[-1])  # 残余压痕深度
            continue
        idx = np.argmin(np.abs(seg_f - tf))
        if np.abs(seg_f[idx] - tf) > force_tol:
            import warnings
            warnings.warn(f"目标力值 {tf} mN 与最近匹配点 {seg_f[idx]:.1f} mN 偏差较大！")
        depths.append(seg_d[idx])

    return np.array(depths)


# ----------------------------------------------------------------------
# 2. 预处理: ALS 基线校正 + 面积归一化
# ----------------------------------------------------------------------

def _als_baseline(y, lam=1e6, p=0.01, niter=10):
    """
    非对称最小二乘基线 (ALS, Eilers & Boelens 2005)。
    lam: 平滑度惩罚；p: 非对称权重。CSC 格式加速稀疏求解。
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2), dtype=float)
    D = (lam * D.dot(D.transpose())).tocsc()

    w = np.ones(L)
    z = y.copy()
    for _ in range(niter):
        W = sparse.diags(w, 0, shape=(L, L), format='csc')
        Z = W + D
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def subtract_local_baseline(wavenumber, intensity, window=(350, 700),
                             fit_range=(180, 1450), lam=1e7, p=0.001, niter=15,
                             return_qc=False):
    """
    在原始采集范围内 (默认180-1450 cm-1) 用 ALS 估计基线，
    再截取分析窗口 (默认350-700 cm-1，覆盖主峰/D1/D2) 并扣除基线。
    默认参数经过手动调参验证：0 mN 光谱校正后极大值落在 ~443 cm-1，
    与文献报道的熔融石英主峰位置吻合 (若基线过度侵入峰复合体，极大值
    会错误漂移至 D1 区域 ~490 cm-1，见正文方法节说明)。
    """
    lo, hi = window
    flo, fhi = fit_range
    fit_mask = (wavenumber >= flo) & (wavenumber <= fhi)
    wn_fit = wavenumber[fit_mask]
    y_fit = intensity[fit_mask]

    baseline_fit = _als_baseline(y_fit, lam=lam, p=p, niter=niter)

    win_mask_local = (wn_fit >= lo) & (wn_fit <= hi)
    wn_win = wn_fit[win_mask_local]
    corr_win = y_fit[win_mask_local] - baseline_fit[win_mask_local]
    negative_fraction = float(np.mean(corr_win < 0))
    corrected = np.clip(corr_win, 0, None)
    if return_qc:
        return wn_win, corrected, {'negative_fraction': negative_fraction,
                                   'fit_range': (float(wn_fit[0]), float(wn_fit[-1]))}
    return wn_win, corrected


def normalize_area(wavenumber, intensity):
    """面积归一化，消除激光功率/积分时间漂移对总强度的影响"""
    area = trapezoid(intensity, wavenumber)
    if area <= 0 or not np.isfinite(area):
        return intensity
    return intensity / area


# ----------------------------------------------------------------------
# 3. 谱带整体矩分析 (envelope-level moments，模型无关)
# ----------------------------------------------------------------------

def envelope_moments(wn, y):
    """
    计算强度加权矩：中心位置(一阶矩)、方差、偏度和超值峰度。
    """
    y = np.clip(y, 0, None)
    area = trapezoid(y, wn)

    empty_res = dict(centroid=np.nan, sigma=np.nan, skewness=np.nan,
                     kurtosis=np.nan, peak_max_wn=np.nan, area=area)
    if area <= 0 or not np.isfinite(area):
        return empty_res

    centroid = trapezoid(wn * y, wn) / area
    var = trapezoid(((wn - centroid) ** 2) * y, wn) / area
    sigma = np.sqrt(var) if var > 0 else 0.0

    if sigma > 0:
        m3 = trapezoid(((wn - centroid) ** 3) * y, wn) / area
        m4 = trapezoid(((wn - centroid) ** 4) * y, wn) / area
        skewness = m3 / (sigma ** 3)
        kurtosis = (m4 / (sigma ** 4)) - 3  # excess kurtosis
    else:
        skewness = kurtosis = np.nan

    return dict(centroid=centroid, sigma=sigma, skewness=skewness,
                kurtosis=kurtosis, peak_max_wn=wn[np.argmax(y)], area=area)

# ----------------------------------------------------------------------
# 4. 三峰 Pseudo-Voigt 联合拟合 (缓存模型 + 安全热启动)
# ----------------------------------------------------------------------

PEAK_INIT = {
    # name: 初始center, center允许下限/上限, 初始sigma, sigma上限
    # 边界比常压文献值适当放宽：高载荷下主峰/D2峰会因局域应力和致密化发生
    # 显著的展宽与偏移，这正是本研究要刻画的现象，边界过窄会人为压制该信号。
    'main': dict(center=440, cmin=380, cmax=480, sigma=35, smax=140),
    'D1':   dict(center=495, cmin=470, cmax=530, sigma=15, smax=60),
    'D2':   dict(center=606, cmin=580, cmax=700, sigma=15, smax=110),
}

_CACHED_MODEL = None
_CACHED_PARAMS = None


def _get_three_peak_model():
    """获取/初始化单例复合模型，避免批量循环中重复解析 lmfit 表达式"""
    global _CACHED_MODEL, _CACHED_PARAMS
    if _CACHED_MODEL is None:
        model = None
        for name in PEAK_INIT:
            pv = PseudoVoigtModel(prefix=f'{name}_')
            model = pv if model is None else (model + pv)
        _CACHED_MODEL = model
        _CACHED_PARAMS = model.make_params()
    return _CACHED_MODEL, _CACHED_PARAMS.copy()


def fit_three_peaks(wn, y, warm_start=None, keep_result_obj=False):
    """
    对预处理后的窗口光谱 (350-700 cm-1) 做三峰 Pseudo-Voigt 联合拟合。

    warm_start: 上一载荷步的拟合结果 dict，其峰位/宽度 (经边界裁剪保护)
                作为本次初值，提升载荷序列拟合的连续性和收敛稳定性。
    keep_result_obj: 是否在返回字典中附带 lmfit ModelResult 对象 ('_result')；
                批量处理长序列时建议 False 以节省内存。
    """
    y = np.clip(y, 1e-6, None)
    model, params = _get_three_peak_model()

    max_y = np.max(y)
    for name, init in PEAK_INIT.items():
        pref = f'{name}_'
        c0, s0 = init['center'], init['sigma']

        if warm_start is not None and name in warm_start:
            c0 = np.clip(warm_start[name]['center'], init['cmin'], init['cmax'])
            s0 = np.clip(warm_start[name]['fwhm'] / 2.355, 3, init['smax'])

        amp_guess = max_y * s0 * 2.0
        params[f'{pref}center'].set(value=c0, min=init['cmin'], max=init['cmax'])
        params[f'{pref}sigma'].set(value=s0, min=3, max=init['smax'])
        params[f'{pref}amplitude'].set(value=amp_guess, min=0)
        params[f'{pref}fraction'].set(value=0.5, min=0, max=1)

    result = model.fit(y, params, x=wn)

    out = {}
    for name in PEAK_INIT:
        pref = f'{name}_'
        out[name] = dict(
            center=result.params[f'{pref}center'].value,
            fwhm=result.params[f'{pref}fwhm'].value,
            height=result.params[f'{pref}height'].value,
            area=result.params[f'{pref}amplitude'].value,
            fraction=result.params[f'{pref}fraction'].value
        )
    out['redchi'] = result.redchi
    if keep_result_obj:
        out['_result'] = result
    return out


def fit_d1_d2_and_residual(wn, y, keep_result_obj=False):
    """Fit D1 and D2 locally and estimate the residual R-band area.

    This is a complementary, less model-dependent analysis. D1 and D2 are
    fitted only in their local spectral regions; the remaining area of the
    measured envelope is reported as ``R_area = total_area - D1_area -
    D2_area`` rather than being assigned to a third fitted peak.
    """
    y = np.clip(y, 1e-6, None)
    total_area = float(trapezoid(y, wn))
    regions = {'D1': (460, 550, 495, 470, 530, 15, 60),
               'D2': (550, 700, 606, 580, 700, 15, 110)}
    out = {}
    for name, (lo, hi, center, cmin, cmax, sigma, smax) in regions.items():
        mask = (wn >= lo) & (wn <= hi)
        model = PseudoVoigtModel(prefix=f'{name}_')
        params = model.make_params()
        params[f'{name}_center'].set(value=center, min=cmin, max=cmax)
        params[f'{name}_sigma'].set(value=sigma, min=3, max=smax)
        params[f'{name}_amplitude'].set(value=max(np.max(y[mask]) * sigma * 2, 1e-6), min=0)
        params[f'{name}_fraction'].set(value=0.5, min=0, max=1)
        result = model.fit(y[mask], params, x=wn[mask])
        fitted_local = result.eval(x=wn[mask])
        local_area = float(trapezoid(fitted_local, wn[mask]))
        out[name] = dict(center=result.params[f'{name}_center'].value,
                         fwhm=result.params[f'{name}_fwhm'].value,
                         height=result.params[f'{name}_height'].value,
                 area=local_area,
                         fraction=result.params[f'{name}_fraction'].value)
        if keep_result_obj:
            out[f'_{name}_result'] = result

    defect_area = out['D1']['area'] + out['D2']['area']
    out['total_area'] = total_area
    out['R_area'] = total_area - defect_area
    out['R_area_nonnegative'] = max(out['R_area'], 0.0)
    out['D1_fraction_total'] = out['D1']['area'] / total_area if total_area > 0 else np.nan
    out['D2_fraction_total'] = out['D2']['area'] / total_area if total_area > 0 else np.nan
    out['R_fraction_total'] = out['R_area'] / total_area if total_area > 0 else np.nan
    return out
