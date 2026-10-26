# 著者確認事項 (Author Queries)

Updated: 2026-10-24（改稿版）

---

## I. 結論・解析結果に影響する事項

### AQ-1: SubGraph.tsvと128状態の対応（確認済みバグ — 最優先修正）

**ステータス: 確認済みバグ**

#### 発見事項

`build_patient_scores.R` L36 (`subgraph[, pattern_id := .I]`) が、SubGraph.tsvの行番号をそのままstate IDとして使用している。しかし47/51疾患では、SubGraph.tsvは孤立局所最小（StatusNetworkに出力辺がない状態）を**スキップ**して非孤立状態のみをstate-ID順に格納している。最初の孤立状態以降、全状態のbasinラベルがずれる。

#### 検証方法

1. 結核のStatusNetwork.tsv（128×128有向隣接行列）を解析 → 10状態が孤立（行の合計=0）
2. SubGraph=118行 = 128 - 10孤立 → **SKIP形式**確認
3. 甲状腺（4 basin, 3孤立: states 55,115,119）で突き合わせ:
   - SubGraph行55のラベル=4 → これはstate 56（55番目の非孤立状態）のラベル
   - patient_scores.tsv: pattern_id=55 → basin_id='4'（state 55は孤立、正しくはnon-basin）
   - patient_scores.tsv: pattern_id=56 → basin_id='1'（state 56はbasin 4、正しくは4）
4. 全51疾患で自動検証 → 30の複数basin疾患で誤割当確認（詳細: `validation/subgraph_reconciliation.tsv`）

#### 影響範囲

- **影響あり**: 47疾患（うち30は複数basinで誤basin割当）
- **影響なし**: 4疾患（dementia, delirium, diabetes, sepsis — SubGraph=128行）
- **影響規模**: 複数basin疾患で観察の1.5%〜61%が誤分類（急性腎不全61%, 心不全52%, CKD 46%）
- **影響する解析**: per-basin Cox (43/51 PASS率), 縦断軌跡, basin別enrichment, Table 1分類

#### 必要な修正

1. `build_patient_scores.R`: StatusNetwork.tsvから孤立状態を特定し、正しいstate-ID対応を構築
2. 47疾患のpatient_scores.tsv再生成
3. 下流解析（Cox, trajectory, tagcloud, enrichment, Table 1）の再実行
4. 再実行後に43/51 PASS率、6パターン分類、basin別統計を再確認

修正コード案: `validation/validation_report.md` に記載

### AQ-2: Table 1のliver_cancerとUTIの分類根拠（AQ-1修正後に再評価）

**ステータス: AQ-1に依存 — 保留**

Table 1ではliver_cancer（Convergent = Predictive）、UTI（Progressive = Predictive）に分類されているが、per-basin mortality majority基準ではいずれも1/3 PASS = FAIL。

**注意**: liver_cancer, UTI共にSubGraphバグの影響を受ける（それぞれ5孤立, 10孤立。multiple basins）。
per-basin HRの計算自体が誤ったbasin構成に基づいている可能性があるため、AQ-1修正後に再評価が必要。

**必要な対応**: AQ-1修正 → patient_scores再生成 → per-basin Cox再実行 → 分類再確認。

### AQ-3: Table 1の境界値付近の疾患

以下の疾患が名目上の閾値と不整合：
- Peptic ulcer: Progressive（<5%/mo）だが遷移率6.636%
- Delirium: Convergent（5–11%/mo）だが遷移率4.622%
- Atrial fibrillation: Convergent（asymmetry 0.3–0.6）だがasymmetry 0.2999
- Acute renal failure: Symmetric（5–10%/mo）だが遷移率10.320%

**必要な対応**: 閾値を連続的目安として扱うか、厳密な境界として記載するか決定。現行原稿では「一部の疾患は境界付近に位置する」と注記。

### AQ-4: Figure 2cの寄与率12.5% vs 固有値計算13.3%（解決済み）

**ステータス: 解決済み — 定義差として本文に明記**

- Figure 2c表示値（中央値12.5%）: 患者代表スコア（dedup後）の分散
- 固有値ファイル（`ooc_cov/eigenvalues.csv`）: 中央値13.3% = Σ(top7 eigenvalues) / Σ(all eigenvalues)
- 差異の原因: `which.max(old_row_id)` による患者代表化（dedup）により、観察窓ベースの分散構造が変化

Methods §4.3とLimitation 3に両値と定義差を明記。Figure 2c captionはそのまま（表示値を正確に参照）。

### AQ-5: 入院KMのイベント定義（Figure 7）（検証済み — 当初の懸念は非該当）

**ステータス: 解決済み**

`km_per_basin_hosp.R` L52-55を検証:
- `first_admission_ymd` から time-to-first-admission を計算
- `event_hosp := as.integer(has_adm)` = 初回入院イベント
- `time_years := ifelse(has_adm, adm_y - start_y, end_y - start_y)` = 適切な打ち切り

Codex版が指摘した「any-admission flag + total observation duration」は実際のコードと一致しない。
`km_per_basin_hosp.R`と`figure_km_comparison_hosp.R`は共にtime-to-first-admission形式。
なお、認知症はfull_128形式（SubGraphバグの影響なし）。

原稿のFigure 7 legendは修正不要（既に適切な記述）。

### AQ-6: 死亡前「月数」の定義

`trajectory_to_death.R` の `periods_before_death` は暦上の月数ではなく、保持された観察レコードの逆順位。欠測窓があっても1と数えられる。Figure 6d、Figure 9aの横軸と、本文中の "months before death" という表現の整合を確認してください。

**必要な対応**: 図の横軸ラベルが "observation positions before last retained record" であるべきか確認。

### AQ-7: Landscaperのバージョン固定

container tag `ghcr.io/chiba-ai-med/landscaper:main` は不変版でない。再現性のためimmutable digestの記録が望ましい。

**必要な対応**: 実行時のimage digestを確認し、Methodsに追記。対応する公開引用があれば提供。

---

## II. 著者情報・投稿準備の不足

### AQ-8: Lead contact

Patternsの投稿規定では、責任著者とは別にLead contactを1名指定する必要がある場合がある。Tsuyuzaki先生とKawakami先生のどちらかを指定してください。

### AQ-9: 倫理承認

倫理委員会名、承認番号、日付、同意免除の有無を提供してください。

### AQ-10: 研究費

Grant番号と受領者を含む研究費情報を提供してください。

### AQ-11: 著者貢献（CRediT）

CRediT形式の著者貢献を提供してください。

### AQ-12: 利益相反

全著者の利益相反の有無を確認してください。

### AQ-13: DeSCデータベースの正式名称・引用

DeSCの正式名称（Descriptive Studies of Claimsの略か）、正式な引用文献を提供してください。

### AQ-14: リポジトリURL

解析コード・ワークフローの公開リポジトリURLを提供してください。

### AQ-15: データアクセスポリシー

DeSCデータの第三者アクセス方法、利用条件を記載してください。

### AQ-16: OnlinePCA.jl引用形式

OnlinePCA.jlのJOSS 2026論文が存在する場合は正式な引用情報を提供してください。GitHubリポジトリのみの場合はその旨確認。

### AQ-17: AI利用声明

Claude以外のAIツール使用（ChatGPTによるFigure 1ワークフロー図ドラフト作成など）の追記が必要か確認してください。著者全員による最終確認・承認後に声明を確定。

### AQ-18: Patternsの投稿規定

公式投稿規定（2026-10-24時点でHTTP 403により取得不可）の確認事項：
- Article typeの選択（Article vs Resource）
- Highlights, eTOC/In Brief, Bigger Picture等の追加必須項目
- 図表数上限（現行: Main 9図 + Table 1 + Supplementary 8図）
- 引用形式
- Key Resources Table, STAR Methodsの要否

### AQ-19: Schizophreniaの除外理由

Figure 4のJaccard解析からschizophreniaを除外した理由を本文に記載すべきか判断してください。

### AQ-20: ユニーク患者総数

約1,100万人の正確なユニーク患者数、データrelease版、保険者種別・年齢範囲を提供してください。
