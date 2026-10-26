# 検証レポート (Validation Report)

Created: 2026-10-25

---

## AQ-1: SubGraph.tsv と 128状態の対応 — 確認済みバグ

### 結論

**`build_patient_scores.R` の行36 `subgraph[, pattern_id := .I]` に確認済みのバグがあり、
47/51疾患でbasin_idの割り当てが誤っている。** 複数basinを持つ30疾患では、
観察の12〜61%が誤ったbasinに分類されている。

### 根拠

#### 1. SubGraph.tsvのフォーマット

Landscaperは `StatusNetwork.tsv`（128×128 有向隣接行列）の構造に基づき、
各状態の最急降下近傍を記録する。StatusNetworkに出力辺（1）がない状態は
「孤立局所最小（isolated local minimum）」であり、主要basinに到達しない。

- **47疾患**: SubGraph.tsv < 128行。孤立状態がスキップされ、非孤立状態のみが
  state-ID順に格納される（**SKIP形式**）。
- **4疾患** (dementia, delirium, diabetes, sepsis): SubGraph.tsv = 128行。
  孤立状態も含めて全128状態にbasinラベルが付与されている（**FULL形式**）。

SKIP形式の確認は甲状腺（thyroid_disorder）で実施:
- 孤立状態: 55, 115, 119（StatusNetwork行がすべて0）
- SubGraph: 125行
- SubGraph行55のラベル=4 → state 56（非孤立の55番目）に対応
- patient_scores.tsv: pattern_id=55 → basin_id='4'（state 55は孤立だが basin 4 が付与）
- patient_scores.tsv: pattern_id=56 → basin_id='1'（state 56はbasin 4のはずが basin 1）
- patient_scores.tsv: pattern_id=126-128 → basin_id=''（非孤立だがSubGraph範囲外）

#### 2. build_patient_scores.Rのロジック

```r
# L34-36: SubGraphを読み込み、行番号をpattern_idに割り当て
subgraph <- fread(file.path(ls_dir, "SubGraph.tsv"), header = FALSE)
setnames(subgraph, "basin_id")
subgraph[, pattern_id := .I]  # 1, 2, ..., N (N < 128)

# L54-66: 患者のpattern_idはAllstates.tsv（128行）から計算 → 1..128
# L69: mergeでpattern_idを結合
scores <- merge(scores, subgraph, by = "pattern_id", all.x = TRUE)
```

SKIP形式のSubGraphでは、行k（pattern_id=k）は「k番目の非孤立状態」であり
「state k」ではない。最初の孤立状態の後、すべての状態でラベルがずれる。

#### 3. 影響の定量化（全51疾患）

| 影響レベル | 疾患数 | 説明 |
|-----------|--------|------|
| 影響なし (full_128) | 4 | dementia, delirium, diabetes, sepsis |
| 単一basin (ラベルずれのみ) | 17 | 誤basin割り当てなし、孤立状態の分類のみ影響 |
| **複数basin (誤basin割り当て)** | **30** | **観察の1.5〜61%が誤ったbasinに分類** |

#### 影響上位の複数basin疾患

| 疾患 | Basins | 孤立 | 誤割当状態 | 誤割当観察% | 誤basin数 |
|------|--------|------|----------|-----------|----------|
| acute_renal_failure | 2 | 9 | 78 | 61.0% | 60 |
| malnutrition | 2 | 15 | 69 | 53.9% | 39 |
| heart_failure | 2 | 6 | 65 | 52.0% | 54 |
| ckd | 2 | 12 | 61 | 46.1% | 37 |
| hypertension | 2 | 4 | 57 | 45.8% | 49 |
| urinary_tract_infection | 3 | 10 | 52 | 40.5% | 33 |
| asthma | 2 | 3 | 56 | 40.1% | 50 |
| rheumatoid_arthritis | 2 | 8 | 50 | 40.3% | 35 |

完全な一覧: `subgraph_reconciliation.tsv`

### 影響を受ける下流解析

1. **Per-basin Cox (43/51 PASS率)**: basin内のエネルギーzスコア、basin内HR計算が
   誤ったbasin構成に基づく → PASSカウントが変動する可能性
2. **Figure 5 (縦断軌跡)**: basin_id=0（Non-basin）の判定が不正確 → 安定率・遷移率に影響
3. **Figure 6-9 (疾患別詳細)**: basin間遷移率、basin別死亡率・入院率が影響
4. **Tagcloud/enrichment**: basin別のコード頻度集計は正しい（SubGraphを使わず、
   Allstatesとbasinの対応をLandscaper出力から直接構築しているため、要確認）
5. **Table 1 (分類)**: 遷移率・非対称度に影響 → パターン分類が変動する可能性

### 影響を受けない解析

1. **全体Cox (atlas_summary direction_ok)**: disease全体のHR。basin_idを使わない。
2. **PCA (Scores, Eigen_vectors)**: Landscaper上流の処理。
3. **Landscaper出力そのもの (Energy, Basin, Allstates, StatusNetwork)**: 正しい。
4. **4疾患 (dementia, delirium, diabetes, sepsis)**: SubGraph=128行で影響なし。

### 実施した修正

#### 1. `src/build_state_basin_map.R`（新規作成）

StatusNetwork.tsv から入次数・出次数をともに計算し、真の孤立頂点（in=0 AND out=0）と
出次数ゼロの局所最小（in>0, out=0）を区別。SubGraph.tsv の行を非孤立状態IDに正しく対応
させる。出力は `state_basin_map.tsv`（128行、state_id / basin_id / energy / isolated）。

全51疾患で生成済み。合成8状態テスト → thyroid_disorder（SKIP, 3孤立） →
dementia（FULL, 0孤立） → 全51疾患で検証。

#### 2. 修正済みスクリプト（6ファイル）

| スクリプト | 修正内容 |
|-----------|---------|
| `src/build_patient_scores.R` | SubGraph直接読み → state_basin_map.tsv |
| `src/cox_sanity_check_simple.R` | pattern_id := .I → state_basin_map.tsv |
| `src/tagcloud_per_basin.R` | SubGraph直接読み → state_basin_map.tsv |
| `src/tagcloud_per_state.R` | SubGraph直接読み → state_basin_map.tsv |
| `src/trajectory_to_death.R` | SubGraph直接読み → state_basin_map.tsv |
| `src/trajectory_to_hospitalization.R` | SubGraph直接読み → state_basin_map.tsv |

変更前のスクリプトは `paper/claude/backup_20261025_pre_subgraph_fix/` に退避済み。

#### 3. 出次数ゼロの再検証

初期の検証レポートでは「StatusNetworkの行和ゼロ＝孤立」としていたが、
StatusNetworkは有向グラフ（各行の1は最急降下近傍）であり、行和=0は出次数=0を
意味する。入次数（列和）も0のとき初めて真の孤立頂点となる。

- **47疾患**: 出次数0の状態はすべて入次数も0（真の孤立）→ SubGraphからスキップが正当
- **dementia**: 状態45, 59が出次数0・入次数>0（局所最小）→ SubGraph 128行に含む
- **diabetes**: 状態25, 33, 72, 96が出次数0・入次数>0 → SubGraph 128行に含む
- **delirium, sepsis**: SubGraph 128行（全状態が非孤立）

`build_state_basin_map.R` はこの区別を正しく実装している。

### 分類

| 項目 | 分類 |
|------|------|
| 深刻度 | **Critical** — 47/51疾患で平均30/128状態のbasin_idが誤り |
| 確認方法 | コード検査 + 合成テスト + 全51疾患突き合わせ + StatusNetwork解析 |
| 修正範囲 | 6スクリプト修正 + state_basin_map.tsv生成 + 全疾患の下流再計算 |
| 原稿への影響 | 全数値の再検証が必要（43/51 PASSなどを仮定しない） |
| 修正状態 | **修正完了・下流再計算完了** — patient_scores, Cox(全体+per-basin), summary すべて再生成済み。結果: per-basin mortality 43→45 PASS, hosp 40→42 PASS。Mortality FAIL変動: breast_cancer, heart_failure, UTI → PASS; lung_cancer → FAIL。Hosp FAIL変動: hypertension, obesity → PASS |

---

## AQ-5: Figure 7 入院KM — 再調査

### 結論

**`km_per_basin_hosp.R` はtime-to-first-admissionを正しく使用しているが、
Figure 7 (`Fig07a_dementia_paired_km.png`) はこのスクリプトで生成されていない。**

#### 1. `km_per_basin_hosp.R` の確認

コード L52-55:
```r
pt[, adm_y := parse_ymd(first_admission_ymd)]
pt[, has_adm := !is.na(adm_y) & adm_y >= start_y & adm_y <= end_y]
pt[, time_years := ifelse(has_adm, adm_y - start_y, end_y - start_y)]
pt[, event_hosp := as.integer(has_adm)]
```
→ time-to-first-admissionを使用。出力は `km_per_basin_hosp.png`（energy tertile
within each basin）と `km_basin_energy_combined_hosp.png`。Basin間比較KMは生成しない。

#### 2. Figure 7の生成元

リポジトリ内にFigure 7を生成するスクリプトは存在しない。git履歴によると：
- commit 337a026 で `plot/Figures/Fig07/Fig06_dementia_paired_km.png` として新規追加（0→420,639 bytes）
- commit 331e9be で `Fig07a_dementia_paired_km.png` にリネーム
- ファイルハッシュはどの既存KM出力とも一致しない → ad hoc作成

#### 3. number-at-risk の不整合

Figure 7(a) Mortality: Basin 1 = 496,150 / Basin 2 = 433,801
Figure 7(b) Hospitalization: Basin 1 = 496,150 / Basin 2 = 433,801 ← **死亡分析と同数**

`km_per_basin_hosp.R` の出力 (plot/dementia/km_per_basin_hosp.png):
Basin 1 = 496,072 / Basin 2 = 433,782 ← 入院データ欠損で78+19=97名少ない

Figure 7(b) が死亡分析と同じnumber-at-riskを示していることは、
(a) 入院パネルが死亡パネルと同じコホートフィルタで作られた、または
(b) number-at-risk表示が死亡パネルからコピーされた
可能性を示唆する。いずれの場合もパネル間のイベント定義の独立性に疑問が残る。

#### 4. 対応方針

- Figure 7は再生成が必要（Basin間比較KM、mortality + hospitalization）
- 入院パネルは `km_per_basin_hosp.R` のフィルタに従い正しいnumber-at-riskを表示すべき
- 再生成スクリプトを作成するか、生成元を特定して修正が必要

### 分類: **解決済み**。`km_inter_basin_paired.R` を survminer 非依存（ggplot2 + survival）に書き直し、
dementia で再生成済み。入院パネルの number-at-risk が正しく表示される（496,072 / 433,782、
死亡パネルの 496,150 / 433,801 と異なる）。

---

## AQ-3: Table 1 分類の整合性

### 現状

AQ-1のSubGraphバグにより、遷移率・非対称度の計算が影響を受ける疾患がある。
AQ-1修正前にTable 1の分類を確定することはできない。

### 既知の問題（AQ-1と独立）

- liver_cancer, UTI: Predictiveに分類されるがper-basin mortality FAIL
- 境界値付近: peptic_ulcer (6.6%), delirium (4.6%), AF (asym 0.300), AKI (10.3%)
- これらの問題はAQ-1修正後に再評価が必要

### 分類: AQ-1修正に依存、保留

---

## AQ-4: PCA variance 12.5% vs 13.3% — 解決済み

### 結論

差異は定義の違いによるもので、両値ともに正しい。

- **13.3%** (eigenvalue-based): `output/*/all_comorbid/ooc_cov/eigenvalues.csv` から
  Σ(top 7 eigenvalues) / Σ(all eigenvalues) → 全51疾患で計算、中央値13.3%
- **12.5%** (Figure 2c): 患者代表スコア（`which.max(old_row_id)` による重複排除後）の
  分散を total variance で割った値

重複排除により観察窓ベースの分散構造が変化するため差が生じる。

Methods §4.3とLimitation 3に両値と定義差を明記完了。

### 分類: 解決済み。
