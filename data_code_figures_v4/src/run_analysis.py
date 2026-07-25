"""
主运行脚本：对 NIST 熔融石英压痕原位拉曼数据集执行完整的峰形演化分析流程。
生成:
  outputs/peak_shape_results.csv   -- 逐光谱峰形参数总表
  outputs/summary_stats.csv        -- 各配置摘要统计
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# 锚定脚本所在目录，使相对路径不受调用时工作目录影响
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR / 'src'))

try:
    from peak_analysis import (load_raman_sequence, load_indentation_curve,
                                force_to_depth_map, subtract_local_baseline,
                                normalize_area, envelope_moments, fit_three_peaks)
except ImportError as e:
    raise RuntimeError(f"无法导入核心分析模块，请确保 '{BASE_DIR/'src'/'peak_analysis.py'}' 存在。") from e

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

DATA_ROOT = BASE_DIR / 'mds2-2281'
OUTPUT_DIR = BASE_DIR / 'outputs'

CONFIGS = [
    dict(objective='TO', specimen='thick_FS', segment='loading',
         raman=DATA_ROOT / 'Dataset C/TO_load_curve_raw_data.csv',
         indent=DATA_ROOT / 'Dataset D/TO_thick_FS_loadsequence_indentation.csv'),
    dict(objective='TO', specimen='thick_FS', segment='unloading',
         raman=DATA_ROOT / 'Dataset C/TO_unload_curve_raw_data.csv',
         indent=DATA_ROOT / 'Dataset D/TO_thick_FS_loadsequence_indentation.csv'),
    dict(objective='CO', specimen='thin_FS', segment='loading',
         raman=DATA_ROOT / 'Dataset F/CO_load_curve_raw_data.csv',
         indent=DATA_ROOT / 'Dataset G/CO_thin_FS_loadsequence_indentation.csv'),
    dict(objective='CO', specimen='thin_FS', segment='unloading',
         raman=DATA_ROOT / 'Dataset F/CO_unload_curve_raw_data.csv',
         indent=DATA_ROOT / 'Dataset G/CO_thin_FS_loadsequence_indentation.csv'),
]


def process_config(cfg):
    """处理单个配置序列的数据读取、基线扣除、特征提取与分峰拟合"""
    if not cfg['raman'].exists() or not cfg['indent'].exists():
        raise FileNotFoundError(f"找不到数据文件，请检查 DATA_ROOT 路径:\n"
                                 f"  Raman: {cfg['raman']}\n  Indent: {cfg['indent']}")

    print(f"正在处理: [{cfg['objective']}-{cfg['specimen']} | {cfg['segment']}] ...", end='', flush=True)

    wn, forces, I = load_raman_sequence(cfg['raman'])
    indent_df = load_indentation_curve(cfg['indent'])
    depths = force_to_depth_map(indent_df, forces, segment=cfg['segment'])

    rows = []
    warm = None
    for i, (f, depth) in enumerate(zip(forces, depths)):
        wn_win, y_win = subtract_local_baseline(wn, I[:, i])
        y_norm = normalize_area(wn_win, y_win)

        env = envelope_moments(wn_win, y_norm)
        fit = fit_three_peaks(wn_win, y_norm, warm_start=warm, keep_result_obj=False)

        # 安全热启动传播：仅在拟合合理时更新初值，阻断错误拟合参数的污染
        redchi = fit.get('redchi', np.inf)
        if np.isfinite(redchi) and redchi < 0.01 and np.isfinite(fit['main']['center']):
            warm = fit
        else:
            warm = None

        row = dict(objective=cfg['objective'], specimen=cfg['specimen'],
                   segment=cfg['segment'], load_mN=f, depth_nm=depth)
        for k in ['centroid', 'sigma', 'skewness', 'kurtosis',
                  'fwhm_left', 'fwhm_right', 'asymmetry', 'peak_max_wn']:
            row[f'env_{k}'] = env.get(k, np.nan)
        for name in ['main', 'D1', 'D2']:
            for k in ['center', 'fwhm', 'height', 'area', 'fraction']:
                row[f'{name}_{k}'] = fit[name].get(k, np.nan)
        row['fit_redchi'] = redchi
        rows.append(row)

    df = pd.DataFrame(rows)
    main_area = df['main_area'].replace(0, np.nan)
    df['ratio_D1_main'] = df['D1_area'] / main_area
    df['ratio_D2_main'] = df['D2_area'] / main_area

    print(f" 完成 ({len(df)} 谱图)")
    return df


def generate_summary_statistics(full_df):
    """按配置聚合核心参数的统计摘要 (mean/std/min/max)"""
    target_cols = [
        'ratio_D1_main', 'ratio_D2_main',
        'env_centroid', 'env_sigma', 'env_asymmetry',
        'main_center', 'main_fwhm', 'D2_center', 'fit_redchi'
    ]
    valid_cols = [c for c in target_cols if c in full_df.columns]
    summary = (full_df
               .groupby(['objective', 'specimen', 'segment'])[valid_cols]
               .agg(['mean', 'std', 'min', 'max'])
               .round(4))
    summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
    return summary.reset_index()


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("NIST 熔融石英压痕拉曼原位分析 - 批量自动化执行")
    print("=" * 60)

    all_dfs = []
    for cfg in CONFIGS:
        try:
            df_cfg = process_config(cfg)
            all_dfs.append(df_cfg)
        except Exception as e:
            print(f"\n[错误] 配置 {cfg['objective']}-{cfg['segment']} 处理失败: {e}")

    if not all_dfs:
        print("[致命错误] 所有配置处理均失败，未生成有效数据。")
        return None

    full = pd.concat(all_dfs, ignore_index=True)
    res_path = OUTPUT_DIR / 'peak_shape_results.csv'
    full.to_csv(res_path, index=False)
    print(f"\n[输出] 逐光谱参数总表已保存至: {res_path} (共 {len(full)} 行)")

    summary_df = generate_summary_statistics(full)
    sum_path = OUTPUT_DIR / 'summary_stats.csv'
    summary_df.to_csv(sum_path, index=False)
    print(f"[输出] 各配置摘要统计表已保存至: {sum_path}")

    print("\n--- 实验配置数据记录统计 ---")
    summary_counts = full.groupby(['objective', 'specimen', 'segment'])['load_mN'].count().reset_index()
    summary_counts.rename(columns={'load_mN': 'spectra_count'}, inplace=True)
    print(summary_counts.to_string(index=False))
    print("=" * 60)

    return full


if __name__ == '__main__':
    main()
