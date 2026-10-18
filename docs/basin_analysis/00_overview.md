# Basin Analysis Overview — 51 Disease Families

Date: 2026-10-18
Data source: `summary/basin_exploration/per_disease_basin_report.tsv` (142 rows, 51 diseases)

## Cross-Cutting Findings

### 1. BasinNA（非Basin群）の二面性 ★★★ 最重要新規知見

BasinNA（どのBasinにも属さない患者群）の死亡率は**疾患によって正反対の方向**を示す。

**BasinNA = 最高死亡率（22疾患）** — 複雑・非典型パターン型:
anxiety (55.1%), atrial_fibrillation (56.9%), ckd (52.1%), hypertension (53.7%),
stroke (54.9%), hematologic_malignancy (51.7%), hip_fracture (51.5%),
lung_cancer (52.0%), pancreatic_cancer (54.6%), colorectal_cancer (48.3%),
interstitial_lung_disease (45.9%), ischemic_heart_disease (44.8%),
depression (38.5%), dyslipidemia (46.8%), malignancy (45.1%),
asthma (40.6%), autoimmune (38.3%), cirrhosis (41.5%),
liver_disease (44.5%), peptic_ulcer (42.2%)

特徴: 高齢、高CCI、多コード → 複雑すぎて典型パターンに当てはまらない重症群

**BasinNA = 最低死亡率（18疾患）** — 若年・スパース型:
breast_cancer (28.2%), copd (34.6%), schizophrenia (30.5%),
epilepsy (32.1%), osteoporosis (38.3%), rheumatoid_arthritis (34.7%),
sleep_disorders (36.4%), anemia (36.3%), aspiration_pneumonia (42.4%),
heart_failure (37.6%), inflammatory_bowel_disease (27.2%),
liver_cancer (44.1%), metastatic_cancer (43.1%), obesity (29.4%),
osteoarthritis (41.8%), stomach_cancer (38.9%), tuberculosis (36.3%),
acute_renal_failure (40.4%), gout (38.1%)

特徴: 若年、低CCI、少コード → 併存症が少なすぎてパターン未確立

**BasinNA なし（4疾患）**: delirium, dementia, diabetes, sepsis

**解釈**: ELAのBasinは「典型的な併存症パターン」を表す。BasinNA患者は2群に分かれる:
- Type A（高リスク）: 複雑すぎて分類不能 → 死亡率最高
- Type B（低リスク）: シンプルすぎて分類対象外 → 死亡率最低

### 2. Basin間死亡率ギャップ Top 10

| Disease | Max death% | Min death% | Gap | N basins |
|---|---|---|---|---|
| thyroid_disorder | 50.4 (Basin2) | 30.0 (Basin4) | 20.4% | 5 |
| hematologic_malignancy | 51.7 (BasinNA) | 31.7 (Basin2) | 20.0% | 3 |
| diabetes | 54.5 (Basin4) | 38.4 (Basin2) | 16.1% | 4 |
| anxiety | 55.1 (BasinNA) | 37.9 (Basin2) | 17.2% | 3 |
| stroke | 54.9 (BasinNA) | 43.4 (Basin1) | 11.5% | 2 |
| sleep_disorders | 46.2 (Basin1) | 36.4 (BasinNA) | 9.8% | 3 |
| osteoporosis | 50.7 (Basin2) | 38.3 (BasinNA) | 12.4% | 3 |
| copd | 45.0 (Basin1) | 34.6 (BasinNA) | 10.4% | 2 |
| asthma | 40.6 (BasinNA) | 31.3 (Basin3) | 9.3% | 3 |
| breast_cancer | 39.2 (Basin1) | 28.2 (BasinNA) | 11.0% | 2 |

### 3. Cross-Disease Basin Similarity

| Pair | Jaccard | Shared codes | Interpretation |
|---|---|---|---|
| heart_failure_B1 ~ ckd_B1 | 0.33 | Z955, I252 | 心腎症候群の自動検出 |
| parkinson_B1 ~ dementia_B1 | 0.33 | Z931, J690, L899 等10件 | 神経変性疾患の終末期共通パターン |
| aspiration_pneumonia_B1 ~ anemia_B1 | 0.25 | 8件 | 衰弱・フレイルパターン |
| atrial_fibrillation_B1 ~ anemia_B1 | 0.25 | 5件 | 心不全＋貧血（機能性鉄欠乏） |
| sepsis_B1 ~ atrial_fibrillation_B1 | 0.25 | 5件 | 重症心疾患パターン |
| epilepsy_B1 ~ anemia_B1 | 0.23 | 7件 | 慢性疾患全般パターン |

### 4. 疾患横断的「フレイル/歯科シグネチャ」★★ (2026-10-18追記)

糖尿病Basin4、甲状腺Basin2、栄養失調Basin1の濃縮コードを詳細に調べた結果、
**いずれも歯科・眼科コードが支配的で、疾患特異的コードはむしろ枯渇**していた。

| Disease Basin | Death% | Top enriched codes |
|---|---|---|
| diabetes Basin4 | 54.5% | K081(歯牙喪失)2.85x, K029(う歯)2.52x, K053(歯周炎)2.16x, H041(涙腺)1.97x |
| thyroid Basin2 | 50.4% | K081 2.70x, T888(手術合併症)2.64x, K052(急性歯周炎)2.37x, K029 2.28x |
| malnutrition Basin1 | 43.4% | H045(結膜炎), H020(眼瞼炎), H524(屈折障害), H353(黄斑変性) |

**解釈**: これらの小規模高死亡率Basinは「疾患の重症サブタイプ」ではなく**フレイル高齢者の
一般的データ指紋**。歯科受診パターン＝多くの医療機関を受診するフレイルマーカー。
口腔フレイルと全身予後の関連（Tanaka 2018 J Dent Res等）と整合的。

**論文の意義**:
- ELAが疾患特異的シグナルだけでなく、**フレイルシグナルも自動検出**することを示す
- 糖尿病Basin2とBasin4は同一コードプロファイルで死亡率16ポイント差 → 年齢・性別が追加の分離軸
- 疾患ごとに見ると「意味不明なBasin」が、横断的に見ると「フレイル共通パターン」として一貫

### 5. 特に注目すべき疾患（専門医相談推奨）

- **糖尿病Basin 4**: 3.1%の患者で死亡率54.5%, ncode中央値48（最多）→ 糖尿病内科
- **甲状腺Basin 2 vs Basin 4**: 50.4% vs 30.0%の20%ギャップ → 内分泌科
- **統合失調症BasinNA**: 最年少(49.6歳)・最低死亡率(30.5%) → 精神科
- **てんかんBasinNA**: 53.4歳で死亡率32.1%（Basin1の40%より低い）→ 神経内科
- **認知症Basin 2**: Basin 1より若い(81.2 vs 83.7)のに死亡率が高い(53.9% vs 49.3%) → 老年科
- **心不全BasinNA**: 2.7%のみ、死亡率最低(37.6%) → 循環器内科
- **敗血症Basin 3**: 入院率93.9%（全Basin中最高）、年齢82.2歳 → 救急/集中治療

## ファイル構成

- `01_cardiovascular.md` — 循環器系5疾患
- `02_cancer.md` — がん10疾患
- `03_neuropsychiatric.md` — 神経精神8疾患
- `04_respiratory.md` — 呼吸器5疾患
- `05_renal_urological.md` — 腎泌尿器4疾患
- `06_metabolic_endocrine.md` — 代謝内分泌5疾患
- `07_gi_hepatic.md` — 消化器肝臓5疾患
- `08_musculoskeletal.md` — 運動器4疾患
- `09_other.md` — その他5疾患
