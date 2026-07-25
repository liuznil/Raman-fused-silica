交付文件说明

```
outputs/
├── 研究报告_拉曼峰形演化分析_SCI版.md   本文件
├── peak_shape_results.csv               27条光谱记录 × 全部峰形/载荷/深度参数
├── summary_stats.csv                     各配置(TO/CO × load/unload)摘要统计
├── correlation_stats.csv                 载荷-峰形参数 Pearson相关系数与p值
├── pca_loadings.csv                      PCA载荷矩阵 (12特征在PC1/PC2上的权重)
├── pca_scores.csv                        PCA得分 + 原始特征
├── figures/ (300 dpi, 英文标注, 适配期刊投稿；图号按论文阅读顺序编号)
│   ├── fig1_fwhm_hysteresis.png            Fig.1 [3.1]
│   ├── fig2_skewness_vs_load.png           Fig.2 [3.1]
│   ├── fig3_D2main_ratio.png               Fig.3 [3.1]
│   ├── fig4_position_vs_skewness.png       Fig.4 [3.1]
│   ├── fig5_spectrum_decomposition.png     Fig.5 [3.2, 新增] 代表性光谱分解
│   ├── fig6_fit_quality_waterfall.png      Fig.6 [3.2, 新增] 全序列拟合质量瀑布图
│   ├── fig7_pristine_vs_residual.png       Fig.7 [3.4]
│   └── fig8_pca_scores.png                 Fig.8 [3.5]
└── src/
    ├── peak_analysis.py   核心模块: 数据读取/ALS基线/矩分析/三峰拟合
    ├── run_analysis.py    批量处理TO/CO的load/unload序列，生成结果表+摘要统计
    └── make_plots.py      生成全部图表 + PCA + 相关性统计检验
```

复现方式：将 NIST 数据集解压后的 `mds2-2281` 文件夹放在与 `src/` 同级目录下，
依次运行 `python src/run_analysis.py` 与 `python src/make_plots.py` 即可。
